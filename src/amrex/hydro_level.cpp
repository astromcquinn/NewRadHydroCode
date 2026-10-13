// The AMReX backend of the unsplit hydro step (hydro_backend = amrex;
// docs/stage_e_hydro.md): hydro_step_unsplit of src/hydro_unsplit.cpp on
// one uniform AMReX level, MPI + OpenMP, built on the same kernel
// (include/radray/kernels/hydro_mh.hpp). Stage E1 of docs/ROADMAP.md.
//
//   * state   MultiFab U (NC = 10 doubles per cell: the 5 conserved
//             variables, then the NS scalars: the internal energy density
//             and the 4 species densities; NG = 2 ghost layers,
//             FillBoundary at the start of every substep);
//   * Es      MultiFab (1 component, NG ghosts): the dual-energy selected
//             thermal energy density, FillBoundary after the selection;
//   * En      MultiFab (1 component, NG ghosts): the total energy after the
//             update and before finalize_cell, whose neighbour maximum the
//             eta1/eta2 rules need (FillBoundary);
//   * Q       FabArray<BaseFab<QCell>> (the 20 predictor doubles per cell,
//             1 ghost layer): computed on the valid box grown by one from
//             the ghost-filled U and Es, so it is never communicated;
//   * faces   per-box FArrayBox temporaries of the FaceData of each face of
//             the box (the HLLC flux, the scalar fluxes, the interface
//             velocity), computed by the box that owns the cells on both
//             sides from the same Q values, so they are never communicated.
//
// Every value is the same function of the same stencil as in the serial
// driver (rrh::face_neighbour_max, select_thermal, signal_sum,
// predict_cell, face_all, update_from_faces, finalize_cell), the only
// reductions are a max (the CFL speed), a min and integer counters, and
// the sub-cycling arithmetic is copied, so the result is BITWISE IDENTICAL
// to hydro_backend = serial for any rank count, box size and thread count
// (tests: ctest -L amrex_hydro).
//
// Between the two half steps the host's GasGrid (u, w, e, on every rank)
// holds the state: a step copies the rank's boxes in from it (no
// communication) and gathers the result back onto every rank
// (each rank broadcasts its boxes in global order, as floats). The redundant
// host state goes away when the PM, the halo finder and the I/O are
// distributed too (Stage E).
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include <AMReX.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_FabArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Reduce.H>
#ifdef AMREX_USE_OMP
#include <omp.h>
#endif

#include "gas_state.hpp"
#include "radray/hydro_backend.hpp"
#include "radray/kernels/hydro_mh.hpp"

namespace radray {

namespace {

namespace ax = amrex;
using ax::Long;
using rrh::NG;
using rrh::NQ;
using rrh::NS;
using rrh::NV;
constexpr int NC = NV + NS;  // components of U: conserved, then the scalars
constexpr int NF = NV + NS + 1;  // components of a face: F, Fs, v

static_assert(std::is_same_v<ax::Real, double>, "build AMReX with double precision (AMReX_PRECISION=DOUBLE)");

// The predictor states of one cell (AoS: the kernels take contiguous Q[NQ])
struct QCell {
  double v[NQ];
};
using QFab = ax::FabArray<ax::BaseFab<QCell>>;

double now() {
  return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

struct Timers {
  double copyin = 0, fill = 0, select = 0, cfl = 0, predict = 0, flux = 0, update = 0, finalize = 0,
         reduce = 0, copyout = 0, allgather = 0;
  Long substeps = 0, calls = 0;
};

// FaceData in a face-centred FArrayBox: NF components
AMREX_FORCE_INLINE void put_face(ax::Array4<double> const& f, int i, int j, int k, const rrh::FaceData& d) {
  for (int m = 0; m < NV; ++m) f(i, j, k, m) = d.F[m];
  for (int s = 0; s < NS; ++s) f(i, j, k, NV + s) = d.Fs[s];
  f(i, j, k, NV + NS) = d.v;
}
AMREX_FORCE_INLINE void get_face(ax::Array4<double const> const& f, int i, int j, int k, rrh::FaceData& d) {
  for (int m = 0; m < NV; ++m) d.F[m] = f(i, j, k, m);
  for (int s = 0; s < NS; ++s) d.Fs[s] = f(i, j, k, NV + s);
  d.v = f(i, j, k, NV + NS);
}

class HydroLevel {
 public:
  // state: amrex_state = 1, the gas state on the level (gas_state.hpp): the
  // level takes its decomposition and copies its boxes in and out locally
  HydroLevel(int N_, int max_grid_size, const gs::GasState* state_)
      : N(N_), domain(ax::IntVect(0), ax::IntVect(N_ - 1)), state(state_) {
    const ax::RealBox rb({0., 0., 0.}, {1., 1., 1.});
    const ax::Array<int, 3> per{1, 1, 1};
    geom.define(domain, rb, ax::CoordSys::cartesian, per);
    period = geom.periodicity();
    if (state) {
      if (state->N != N) throw std::runtime_error("hydro_backend = amrex: the gas state has another N");
      ba = state->ba;
      dm = state->dm;
    } else {
      ba.define(domain);
      ba.maxSize(max_grid_size);
      dm.define(ba);
    }
    U.define(ba, dm, NC, NG);
    Es.define(ba, dm, 1, NG);
    En.define(ba, dm, 1, NG);
    Q.define(ba, dm, 1, 1);
    const int nloc = static_cast<int>(U.local_size());
    for (int d = 0; d < 3; ++d) faces[d].resize(nloc);
    if (state) return;  // no host gather
    // the gather plan: every box in global order, owner and size
    const int nranks = ax::ParallelDescriptor::NProcs();
    counts.assign(nranks, 0);
    displs.assign(nranks, 0);
    box_offset.assign(ba.size(), 0);
    std::vector<Long> per_rank(nranks, 0);
    for (int b = 0; b < static_cast<int>(ba.size()); ++b) per_rank[dm[b]] += ba[b].numPts() * NC;
    Long off = 0;
    for (int r = 0; r < nranks; ++r) {
      counts[r] = static_cast<int>(per_rank[r]);
      displs[r] = static_cast<int>(off);
      if (per_rank[r] > INT32_MAX || off > INT32_MAX)
        throw std::runtime_error("hydro_backend = amrex: the host gather needs > 2^31 floats per rank; use more ranks");
      off += per_rank[r];
    }
    // the offset of each box inside its owner's block, in the global box order
    std::vector<Long> next(nranks);
    for (int r = 0; r < nranks; ++r) next[r] = displs[r];
    for (int b = 0; b < static_cast<int>(ba.size()); ++b) {
      box_offset[b] = next[dm[b]];
      next[dm[b]] += ba[b].numPts() * NC;
    }
    total = off;
  }

  int N;
  ax::Box domain;
  ax::Geometry geom;
  ax::Periodicity period;
  ax::BoxArray ba;
  ax::DistributionMapping dm;
  ax::MultiFab U, Es, En;
  QFab Q;
  ax::Vector<ax::FArrayBox> faces[3];
  std::vector<int> counts, displs;
  std::vector<Long> box_offset;
  Long total = 0;
  std::vector<float> sendbuf, recvbuf;
  const gs::GasState* state = nullptr;
  Timers T;

  // ---- the gas state on the level (amrex_state = 1): local copies -----------

  void copy_in(const gs::GasState& st) {
    const double t0 = now();
    for (ax::MFIter mfi(U); mfi.isValid(); ++mfi) {
      auto u = U.array(mfi);
      auto g = st.G.const_array(mfi);
      ax::ParallelForOMP(mfi.validbox(), NC, [=] AMREX_GPU_DEVICE(int i, int j, int k, int m) {
        u(i, j, k, m) = g(i, j, k, m);  // the same component order: u(5), e, w(4)
      });
    }
    T.copyin += now() - t0;
  }
  void copy_out(gs::GasState& st) {
    const double t0 = now();
    for (ax::MFIter mfi(U); mfi.isValid(); ++mfi) {
      auto u = U.const_array(mfi);
      auto g = st.G.array(mfi);
      ax::ParallelForOMP(mfi.validbox(), NC, [=] AMREX_GPU_DEVICE(int i, int j, int k, int m) {
        g(i, j, k, m) = static_cast<float>(u(i, j, k, m));  // rounded once, as the serial step
      });
    }
    T.copyout += now() - t0;
  }

  // ---- host <-> level ------------------------------------------------------

  // The rank's boxes from the host arrays (float -> double)
  void copy_in(const GasGrid& gas) {
    const double t0 = now();
    const float *gu = gas.u.data(), *gw = gas.w.data(), *ge = gas.e.data();
    const int NN = N;
    for (ax::MFIter mfi(U); mfi.isValid(); ++mfi) {
      auto u = U.array(mfi);
      ax::ParallelForOMP(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const std::int64_t cg = i + static_cast<std::int64_t>(NN) * (j + static_cast<std::int64_t>(NN) * k);
        for (int m = 0; m < NV; ++m) u(i, j, k, m) = gu[cg * 5 + m];
        u(i, j, k, NV) = ge[cg];
        for (int s = 0; s < 4; ++s) u(i, j, k, NV + 1 + s) = gw[cg * 4 + s];
      });
    }
    T.copyin += now() - t0;
  }

  // The whole level onto every rank's host arrays (double -> float once,
  // as the serial step rounds once): each rank packs its boxes in global
  // box order and broadcasts the block; every rank unpacks each block as it
  // arrives, so the extra memory is one block, not the level.
  void copy_out(GasGrid& gas) {
    const double t0 = now();
    const int nranks = ax::ParallelDescriptor::NProcs();
    const int myproc = ax::ParallelDescriptor::MyProc();
    sendbuf.resize(static_cast<std::size_t>(counts[myproc]));
    for (ax::MFIter mfi(U); mfi.isValid(); ++mfi) {
      auto u = U.const_array(mfi);
      const ax::Box bx = mfi.validbox();
      const auto lo = ax::lbound(bx), len = ax::length(bx);
      float* dst = sendbuf.data() + (box_offset[mfi.index()] - displs[myproc]);
      ax::ParallelForOMP(bx, NC, [=] AMREX_GPU_DEVICE(int i, int j, int k, int m) {
        const Long c = (i - lo.x) + static_cast<Long>(len.x) * ((j - lo.y) + static_cast<Long>(len.y) * (k - lo.z));
        dst[c * NC + m] = static_cast<float>(u(i, j, k, m));
      });
    }
    float *gu = gas.u.data(), *gw = gas.w.data(), *ge = gas.e.data();
    const int NN = N;
    const int nbox = static_cast<int>(ba.size());
    // the boxes of each rank, in global order (what its block holds)
    std::vector<std::vector<int>> boxes_of(nranks);
    for (int b = 0; b < nbox; ++b) boxes_of[dm[b]].push_back(b);
    for (int r = 0; r < nranks; ++r) {
      const float* blk = sendbuf.data();
      if (nranks > 1) {
        const double t1 = now();
        if (r != myproc) {
          recvbuf.resize(static_cast<std::size_t>(counts[r]));
          blk = recvbuf.data();
        }
#ifdef AMREX_USE_MPI
        MPI_Bcast(const_cast<float*>(blk), counts[r], MPI_FLOAT, r, ax::ParallelDescriptor::Communicator());
#endif
        T.allgather += now() - t1;
      }
      const std::vector<int>& bl = boxes_of[r];
      const int nb = static_cast<int>(bl.size());
#ifdef AMREX_USE_OMP
#pragma omp parallel for schedule(dynamic)
#endif
      for (int ib = 0; ib < nb; ++ib) {
        const int b = bl[ib];
        const ax::Box bx = ba[b];
        const auto lo = ax::lbound(bx), hi = ax::ubound(bx), len = ax::length(bx);
        const float* src = blk + (box_offset[b] - displs[r]);
        for (int k = lo.z; k <= hi.z; ++k)
          for (int j = lo.y; j <= hi.y; ++j)
            for (int i = lo.x; i <= hi.x; ++i) {
              const Long c = (i - lo.x) + static_cast<Long>(len.x) * ((j - lo.y) + static_cast<Long>(len.y) * (k - lo.z));
              const float* v = src + c * NC;
              const std::int64_t cg = i + static_cast<std::int64_t>(NN) * (j + static_cast<std::int64_t>(NN) * k);
              for (int m = 0; m < NV; ++m) gu[cg * 5 + m] = v[m];
              ge[cg] = v[NV];
              for (int s = 0; s < 4; ++s) gw[cg * 4 + s] = v[NV + 1 + s];
            }
      }
    }
    T.copyout += now() - t0;
  }

  // ---- the step -------------------------------------------------------------

  // The selected thermal energy density of every valid cell, then its ghosts
  void select(const rrh::DualEnergy& de) {
    const double t0 = now();
    for (ax::MFIter mfi(U); mfi.isValid(); ++mfi) {
      auto u = U.const_array(mfi);
      auto es = Es.array(mfi);
      const double eta1 = de.eta1;
      const int radius = de.nb_radius;
      const bool local = de.eta1_local;
      ax::ParallelForOMP(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        double Uc[NV];
        for (int m = 0; m < NV; ++m) Uc[m] = u(i, j, k, m);
        double Eref = Uc[4], e1;
        if (!local) {
          auto Ea = [&](int ii, int jj, int kk) { return u(ii, jj, kk, 4); };
          Eref = rrh::face_neighbour_max(Ea, i, j, k, radius, e1);
        }
        es(i, j, k) = rrh::select_thermal(Uc, u(i, j, k, NV), Eref, eta1);
      });
    }
    T.select += now() - t0;
    const double t1 = now();
    Es.FillBoundary(period);
    T.fill += now() - t1;
  }

  // Largest signal speed over the level (max: exact, order-independent)
  double max_signal() {
    const double t0 = now();
    double smax = 0;
    for (ax::MFIter mfi(U); mfi.isValid(); ++mfi) {
      auto u = U.const_array(mfi);
      auto es = Es.const_array(mfi);
      const ax::Box bx = mfi.validbox();
      double bmax = 0;
#ifdef AMREX_USE_OMP
#pragma omp parallel reduction(max : bmax)
#endif
      {
        const auto lo = ax::lbound(bx), hi = ax::ubound(bx);
#ifdef AMREX_USE_OMP
#pragma omp for schedule(static) collapse(2)
#endif
        for (int k = lo.z; k <= hi.z; ++k)
          for (int j = lo.y; j <= hi.y; ++j)
            for (int i = lo.x; i <= hi.x; ++i) {
              double Uc[NV];
              for (int m = 0; m < NV; ++m) Uc[m] = u(i, j, k, m);
              bmax = rrf::max(bmax, rrh::signal_sum(Uc, es(i, j, k)));
            }
      }
      smax = rrf::max(smax, bmax);
    }
    ax::ParallelDescriptor::ReduceRealMax(smax);
    T.cfl += now() - t0;
    return smax;
  }

  struct Counts {
    Long v[5] = {0, 0, 0, 0, 0};  // from_e, synced, floor_T, floor_rho, reset
  };

  // One substep of length dt on the ghost-filled U (Es selected and filled)
  Counts substep(double dt, const rrh::DualEnergy& de, double rho_min, double Tmin, bool p_old) {
    const double dtdx = dt;  // dx = 1 cell
    // 1. predictor on the valid box plus one ghost layer
    double t0 = now();
    for (ax::MFIter mfi(U); mfi.isValid(); ++mfi) {
      auto u = U.const_array(mfi);
      auto es = Es.const_array(mfi);
      auto q = Q.array(mfi);
      ax::ParallelForOMP(ax::grow(mfi.validbox(), 1), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        auto Ua = [&](int ii, int jj, int kk, int m) { return u(ii, jj, kk, m); };
        auto Ea = [&](int ii, int jj, int kk) { return es(ii, jj, kk); };
        rrh::predict_cell(Ua, Ea, i, j, k, dtdx, q(i, j, k).v);
      });
    }
    double t1 = now();
    T.predict += t1 - t0;
    // 2. the faces of every valid cell (box-local) and the update in place
    for (ax::MFIter mfi(U); mfi.isValid(); ++mfi) {
      const ax::Box bx = mfi.validbox();
      const int li = mfi.LocalIndex();
      auto q = Q.const_array(mfi);
      auto u = U.array(mfi);
      auto uc = U.const_array(mfi);
      for (int d = 0; d < 3; ++d) {
        const ax::Box fbx = ax::surroundingNodes(bx, d);
        faces[d][li].resize(fbx, NF);
        auto f = faces[d][li].array();
        const int ei = d == 0, ej = d == 1, ek = d == 2;
        ax::ParallelForOMP(fbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
          auto Qa = [&](int ii, int jj, int kk) -> const double* { return q(ii, jj, kk).v; };
          auto Ua = [&](int ii, int jj, int kk, int m) { return uc(ii, jj, kk, m); };
          auto Sa = [&](int ii, int jj, int kk, int s) { return uc(ii, jj, kk, NV + s); };
          rrh::FaceData fd;
          rrh::face_all(Qa, Ua, Sa, i - ei, j - ej, k - ek, d, dtdx, fd);
          put_face(f, i, j, k, fd);
        });
      }
      double t2 = now();
      T.flux += t2 - t1;
      t1 = t2;
      auto fx = faces[0][li].const_array(), fy = faces[1][li].const_array(), fz = faces[2][li].const_array();
      auto es = Es.const_array(mfi);
      auto en = En.array(mfi);
      ax::ParallelForOMP(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        rrh::FaceData x0, x1, y0, y1, z0, z1;
        get_face(fx, i, j, k, x0);
        get_face(fx, i + 1, j, k, x1);
        get_face(fy, i, j, k, y0);
        get_face(fy, i, j + 1, k, y1);
        get_face(fz, i, j, k, z0);
        get_face(fz, i, j, k + 1, z1);
        double Uc[NV], S[NS], Un[NV], Sn[NS];
        for (int m = 0; m < NV; ++m) Uc[m] = u(i, j, k, m);
        for (int s = 0; s < NS; ++s) S[s] = u(i, j, k, NV + s);
        // the cell pressure of the work term: time-centred (Hancock) or at
        // the old time (Bryan et al. 2014 eq. 43)
        const double p_half = p_old ? (rrh::kGamma - 1) * rrf::max(es(i, j, k), 0.0) : q(i, j, k).v[4];
        rrh::update_from_faces(Uc, S, x0, x1, y0, y1, z0, z1, p_half, dtdx, Un, Sn);
        for (int m = 0; m < NV; ++m) u(i, j, k, m) = Un[m];
        for (int s = 0; s < NS; ++s) u(i, j, k, NV + s) = Sn[s];
        en(i, j, k) = Un[4];
      });
      t2 = now();
      T.update += t2 - t1;
      t1 = t2;
    }
    // 3. floors, dual-energy selection and synchronisation (the neighbour
    // maximum of the updated E through En's ghosts)
    double t3 = now();
    En.FillBoundary(period);
    double t4 = now();
    T.fill += t4 - t3;
    Counts c;
    for (ax::MFIter mfi(U); mfi.isValid(); ++mfi) {
      auto u = U.array(mfi);
      auto en = En.const_array(mfi);
      Long* cv = c.v;
      const int radius = de.nb_radius;
      ax::ParallelForOMP(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        auto Ena = [&](int ii, int jj, int kk) { return en(ii, jj, kk); };
        double Enb1;
        const double Enb = rrh::face_neighbour_max(Ena, i, j, k, radius, Enb1);
        double Uc[NV];
        for (int m = 0; m < NV; ++m) Uc[m] = u(i, j, k, m);
        double eD = u(i, j, k, NV);
        const rrh::CellFlags f = rrh::finalize_cell(Uc, eD, Enb, Enb1, de, rho_min, Tmin);
        for (int m = 0; m < NV; ++m) u(i, j, k, m) = Uc[m];
        u(i, j, k, NV) = eD;
        if (f.from_e) ax::HostDevice::Atomic::Add(cv + 0, Long(1));
        if (f.synced) ax::HostDevice::Atomic::Add(cv + 1, Long(1));
        if (f.floor_T) ax::HostDevice::Atomic::Add(cv + 2, Long(1));
        if (f.floor_rho) ax::HostDevice::Atomic::Add(cv + 3, Long(1));
        if (f.reset) ax::HostDevice::Atomic::Add(cv + 4, Long(1));
      });
    }
    T.finalize += now() - t4;
    return c;
  }

  void step(const Config& cfg, const Units& un, double dt, GasGrid& gas, HydroState& hs) {
    if (state)
      copy_in(*state);
    else
      copy_in(gas);
    const rrh::DualEnergy de{cfg.dual_energy ? cfg.dual_eta1 : 0.0, cfg.dual_eta2, cfg.dual_eta1_local,
                             cfg.dual_nb_radius, cfg.dual_reset_rule, cfg.dual_eta2_local};
    const double rho_min = un.Dgas_min;
    const double Tmin = 1 / un.temp1_unit;
    UnsplitCounters& C = hs.ucounters;
    Counts total_c;
    double remaining = dt;
    Long nsub_call = 0;
    while (remaining > 0) {
      const double t0 = now();
      U.FillBoundary(period);
      T.fill += now() - t0;
      select(de);
      const double smax = max_signal();
      if (!(smax > 0) || !std::isfinite(smax))
        throw std::runtime_error("hydro_step_unsplit (amrex): non-finite or zero signal speed");
      // the sub-cycling arithmetic of hydro_step_unsplit, unchanged
      const double dt_cfl = cfg.unsplit_cfl / smax;
      C.dt_cfl_min = std::min(C.dt_cfl_min, dt_cfl);
      int nsub = static_cast<int>(std::ceil(remaining / dt_cfl - 1e-9));
      if (nsub < 1) nsub = 1;
      const double dt_sub = nsub == 1 ? remaining : remaining / nsub;
      const Counts sc = substep(dt_sub, de, rho_min, Tmin, cfg.dual_p_old);
      for (int m = 0; m < 5; ++m) total_c.v[m] += sc.v[m];
      ++nsub_call;
      remaining = nsub == 1 ? 0.0 : remaining - dt_sub;
    }
    const double t1 = now();
    ax::ParallelDescriptor::ReduceLongSum(total_c.v, 5);
    T.reduce += now() - t1;
    C.cells_from_e += total_c.v[0];
    C.cells_synced += total_c.v[1];
    C.floor_T += total_c.v[2];
    C.floor_rho += total_c.v[3];
    C.cells_reset += total_c.v[4];
    C.calls += 1;
    C.substeps += nsub_call;
    C.max_substeps = std::max(C.max_substeps, static_cast<std::int64_t>(nsub_call));
    T.substeps += nsub_call;
    T.calls += 1;
    if (state)
      copy_out(*const_cast<gs::GasState*>(state));
    else
      copy_out(gas);
  }

  std::string timers() const {
    char buf[512];
    std::snprintf(buf, sizeof buf,
                  "AMReX hydro timers [s]: copyin %.3f fill %.3f select %.3f cfl %.3f predict %.3f flux %.3f "
                  "update %.3f finalize %.3f reduce %.3f copyout %.3f (bcast %.3f); %lld substeps in %lld "
                  "calls, %d boxes",
                  T.copyin, T.fill, T.select, T.cfl, T.predict, T.flux, T.update, T.finalize, T.reduce, T.copyout,
                  T.allgather, static_cast<long long>(T.substeps), static_cast<long long>(T.calls),
                  static_cast<int>(ba.size()));
    return buf;
  }
};

std::unique_ptr<HydroLevel> g_level;

HydroLevel& level_of(const Config& cfg, int N) {
  if (!amrex::Initialized())
    throw std::runtime_error("hydro_backend = amrex: backend_init was not called (AMReX is not initialised)");
  if (!g_level || g_level->N != N) {
    check_hydro_backend_param(cfg);
    const int mgs = cfg.amrex_hydro_max_grid_size ? cfg.amrex_hydro_max_grid_size : cfg.amrex_max_grid_size;
    if (N % 2 != 0 && mgs < N)
      throw std::runtime_error("hydro_backend = amrex: N must be even or amrex_hydro_max_grid_size >= N");
    const gs::GasState* st = cfg.amrex_state ? gs::active() : nullptr;
    if (cfg.amrex_state && !st) throw std::runtime_error("amrex_state = 1: amrex_state_init was not called");
    g_level = std::make_unique<HydroLevel>(N, mgs, st);
  }
  return *g_level;
}

}  // namespace

void amrex_hydro_step(const Config& cfg, const Units& un, double dt, GasGrid& gas, HydroState& hs) {
  if (!gas.has_e()) throw std::runtime_error("amrex_hydro_step: call init_unsplit_state first");
  level_of(cfg, gas.N).step(cfg, un, dt, gas, hs);
}

std::string amrex_hydro_timers() { return g_level ? g_level->timers() : ""; }
void amrex_hydro_release() { g_level.reset(); }

}  // namespace radray
