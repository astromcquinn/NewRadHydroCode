// The AMReX backend of the unsplit hydro step (hydro_backend = amrex;
// docs/stage_e_hydro.md, docs/stage_f.md): hydro_step_unsplit of
// src/hydro_unsplit.cpp on an AMReX level hierarchy, MPI + OpenMP, built on
// the same kernel (include/radray/kernels/hydro_mh.hpp). Stage E1 (one
// level) and Stage F1 (amrex_max_level > 0) of docs/ROADMAP.md.
//
// Per level:
//   * state   MultiFab U (NC = 10 doubles per cell: the 5 conserved
//             variables, then the NS scalars: the internal energy density
//             and the 4 species densities; NG = 2 ghost layers);
//   * Uold    the state before the level's substep (the time interpolation
//             of the finer level's ghosts; levels below the finest only);
//   * Es      MultiFab (1 component, NG ghosts): the dual-energy selected
//             thermal energy density;
//   * En      MultiFab (1 component, NG ghosts): the total energy after the
//             update and before finalize_cell, whose neighbour maximum the
//             eta1/eta2 rules need;
//   * Q       FabArray<BaseFab<QCell>> (the 20 predictor doubles per cell,
//             1 ghost layer): computed on the valid box grown by one from
//             the ghost-filled U and Es, so it is never communicated;
//   * flux    three face-centred MultiFabs of the FaceData of each face of
//             the level's boxes (the HLLC flux, the scalar fluxes, the
//             interface velocity), computed by the box that owns the cells
//             on both sides from the same Q values, so they are never
//             communicated; levels > 0 add them to the flux register with
//             their parent;
//   * fr      FluxRegister with the parent (levels > 0).
//
// Every value is the same function of the same stencil as in the serial
// driver (rrh::face_neighbour_max, select_thermal, signal_sum,
// predict_cell, face_all, update_from_faces, finalize_cell), the only
// reductions are a max (the CFL speed), a min and integer counters, and
// the sub-cycling arithmetic is copied, so with amrex_max_level = 0 the
// result is BITWISE IDENTICAL to hydro_backend = serial for any rank
// count, box size and thread count (tests: ctest -L amrex_hydro).
//
// With amrex_max_level > 0 (docs/stage_f.md section 1): the refined
// levels are tagged on the gas overdensity (amrex_refine_mass * 8^l) and
// regridded by AmrMesh every amrex_regrid_every calls; a substep of length
// dt on level l is followed by two substeps of dt/2 on level l+1 whose
// ghosts come from FillPatchTwoLevels (cell-conservative linear
// interpolation in space, linear in time between Uold and U of the
// parent), then the parent is refluxed and the finer level averaged down
// onto it. Because dx and dt halve together the kernel's dtdx is the same
// on every level. The RT and the PM act on level 0 only until Stage F3
// and Stage G: their changes to the level-0 gas between hydro calls are
// pushed down to the fine cells by ratio (push_down below).
//
// Between the two half steps the gas state (amrex_state = 1) or the
// host's GasGrid (u, w, e, on every rank) holds level 0: a step copies
// the rank's boxes in and the result out (the host case gathers the whole
// level onto every rank). The refined levels persist in the hierarchy
// (and in the checkpoint through amrex_hydro_write_levels).
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include <AMReX.H>
#include <AMReX_AmrMesh.H>
#include <AMReX_BCRec.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_FabArray.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_FluxRegister.H>
#include <AMReX_Geometry.H>
#include <AMReX_Interpolater.H>
#include <AMReX_MFInterpolater.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_Reduce.H>
#include <AMReX_TagBox.H>
#include <AMReX_Utility.H>
#include <AMReX_VisMF.H>
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
constexpr int kRef = 2;          // refinement ratio per direction

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
         reduce = 0, copyout = 0, allgather = 0, regrid = 0, reflux = 0, pushdown = 0, sink = 0;
  Long substeps = 0, calls = 0, regrids = 0, sink_events = 0;
  double sink_mass = 0;  // this rank's, in mean-matter x level-0 cell volumes
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

struct Counts {
  Long v[5] = {0, 0, 0, 0, 0};  // from_e, synced, floor_T, floor_rho, reset
};

// One level of the hierarchy
struct LevelData {
  ax::BoxArray ba;
  ax::DistributionMapping dm;
  ax::MultiFab U, Uold, Es, En;
  QFab Q;
  ax::MultiFab flux[3];
  std::unique_ptr<ax::FluxRegister> fr;  // with the parent (lev > 0)
  double t_old = 0, t_new = 0;           // the times of Uold and U (children's interpolation)
  bool defined = false;

  void define(const ax::BoxArray& ba_, const ax::DistributionMapping& dm_, int lev) {
    ba = ba_;
    dm = dm_;
    U.define(ba, dm, NC, NG);
    Uold.define(ba, dm, NC, 0);
    Es.define(ba, dm, 1, NG);
    En.define(ba, dm, 1, NG);
    Q.define(ba, dm, 1, 1);
    for (int d = 0; d < 3; ++d) flux[d].define(ax::convert(ba, ax::IntVect::TheDimensionVector(d)), dm, NF, 0);
    fr.reset();
    if (lev > 0) fr = std::make_unique<ax::FluxRegister>(ba, dm, ax::IntVect(kRef), lev, NC);
    defined = true;
  }
  void clear() {
    U.clear();
    Uold.clear();
    Es.clear();
    En.clear();
    Q.clear();
    for (auto& f : flux) f.clear();
    fr.reset();
    ba = ax::BoxArray();
    dm = ax::DistributionMapping();
    defined = false;
  }
};

class HydroLevel;

// The gridding (AmrMesh): tags from the owner's densities
struct Mesh : ax::AmrMesh {
  HydroLevel* owner = nullptr;
  Mesh(const ax::Geometry& g0, const ax::AmrInfo& info) : ax::AmrMesh(g0, info) {}
  void ErrorEst(int lev, ax::TagBoxArray& tags, ax::Real, int) override;
};

class HydroLevel {
 public:
  // state: amrex_state = 1, the gas state on the level (gas_state.hpp): the
  // level takes its decomposition and copies its boxes in and out locally
  HydroLevel(const Config& cfg, int N_, int max_grid_size, gs::GasState* state_)
      : N(N_), domain(ax::IntVect(0), ax::IntVect(N_ - 1)), state(state_), max_level(cfg.amrex_max_level),
        refine_mass(cfg.amrex_refine_mass), regrid_every(cfg.amrex_regrid_every) {
    for (int m = 0; m < 24; ++m) refine_box[m] = cfg.amrex_refine_box[m];
    refine_nbox = cfg.amrex_refine_nbox;
    const ax::RealBox rb({0., 0., 0.}, {double(N_), double(N_), double(N_)});  // dx = 1 (cell units)
    const ax::Array<int, 3> per{1, 1, 1};
    geom.define(domain, rb, ax::CoordSys::cartesian, per);
    period = geom.periodicity();
    levels.resize(static_cast<std::size_t>(max_level) + 1);
    LevelData& L0 = levels[0];
    ax::BoxArray ba;
    ax::DistributionMapping dm;
    if (state) {
      if (state->N != N) throw std::runtime_error("hydro_backend = amrex: the gas state has another N");
      ba = state->ba;
      dm = state->dm;
    } else {
      ba.define(domain);
      ba.maxSize(max_grid_size);
      dm.define(ba);
    }
    L0.define(ba, dm, 0);
    if (max_level > 0) {
      ax::AmrInfo info;
      info.max_level = max_level;
      info.ref_ratio.assign(1, ax::IntVect(kRef));
      info.blocking_factor.assign(1, ax::IntVect(cfg.amrex_blocking_factor));
      // level 0: the level's own box size; finer levels: amrex_fine_max_grid_size
      info.max_grid_size = {ax::IntVect(max_grid_size), ax::IntVect(cfg.amrex_fine_max_grid_size)};
      info.n_error_buf.assign(1, ax::IntVect(cfg.amrex_n_error_buf));
      info.grid_eff = cfg.amrex_grid_eff;
      if (N % cfg.amrex_blocking_factor != 0)
        throw std::runtime_error("amrex_max_level > 0: N must be a multiple of amrex_blocking_factor");
      mesh = std::make_unique<Mesh>(geom, info);
      mesh->owner = this;
      mesh->SetBoxArray(0, ba);
      mesh->SetDistributionMap(0, dm);
      mesh->SetFinestLevel(0);
      prev0.define(ba, dm, NC, 0);
      bcs.assign(NC, ax::BCRec(ax::BCType::int_dir, ax::BCType::int_dir, ax::BCType::int_dir, ax::BCType::int_dir,
                               ax::BCType::int_dir, ax::BCType::int_dir));
    }
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
  ax::Geometry geom;  // level 0
  ax::Periodicity period;
  std::vector<LevelData> levels;
  std::vector<int> counts, displs;
  std::vector<Long> box_offset;
  Long total = 0;
  std::vector<float> sendbuf, recvbuf;
  gs::GasState* state = nullptr;
  Timers T;
  // the hierarchy (amrex_max_level > 0)
  int max_level = 0;
  int finest = 0;
  int finest_ever = 0;  // the deepest level used so far (the timers line)
  double refine_mass = 8;
  int regrid_every = 1;
  int refine_box[24] = {0, 0, 0, -1, -1, -1};
  int refine_nbox = 0;
  std::unique_ptr<Mesh> mesh;
  ax::MultiFab prev0;  // level 0 as last written out (the push-down reference)
  bool have_prev0 = false;
  ax::Vector<ax::BCRec> bcs;
  double rho_avg = 0;  // the mean gas density (tagging)
  std::vector<SinkEvent> events;  // the sink's events of the last call (this rank)
  double mass0 = 0;  // the gas mass of level 0 after the last call (double sum; diagnostic)

  const ax::Geometry& geom_of(int lev) const { return max_level > 0 ? mesh->Geom(lev) : geom; }
  // The coarse-to-fine interpolation: cell-conservative linear with the
  // slope vector limited so that no fine value leaves the range of the
  // coarse neighbours (positivity in strongly contrasted fields; the
  // per-direction limiters of cell_cons_interp / lincc_interp produced
  // negative densities in random3d). FillPatch dispatches the MF
  // interpolaters by the exact pointer type.
  static ax::MFInterpolater* interpolater() { return &ax::mf_linear_slope_minmax_interp; }

  // RADRAY_AMR_DEBUG=1: scan a level's U (valid cells, or with ghosts) for
  // non-finite or non-positive densities after a phase and report the first
  bool debug = std::getenv("RADRAY_AMR_DEBUG") != nullptr;
  void check(int lev, const char* where, bool ghosts) {
    if (!debug) return;
    const ax::MultiFab& Ul = U(lev);
    for (ax::MFIter mfi(Ul); mfi.isValid(); ++mfi) {
      auto u = Ul.const_array(mfi);
      const ax::Box bx = ghosts ? ax::grow(mfi.validbox(), NG) : mfi.validbox();
      const auto lo = ax::lbound(bx), hi = ax::ubound(bx);
      for (int k = lo.z; k <= hi.z; ++k)
        for (int j = lo.y; j <= hi.y; ++j)
          for (int i = lo.x; i <= hi.x; ++i) {
            bool bad = !(u(i, j, k, 0) > 0);
            for (int m = 0; m < NC && !bad; ++m) bad = !std::isfinite(u(i, j, k, m));
            if (bad) {
              std::fprintf(stderr, "AMR debug: level %d %s: cell (%d,%d,%d)%s rho %g p (%g,%g,%g) E %g e %g w (%g,%g,%g,%g)\n",
                           lev, where, i, j, k, mfi.validbox().contains(ax::IntVect(i, j, k)) ? "" : " [ghost]",
                           u(i, j, k, 0), u(i, j, k, 1), u(i, j, k, 2), u(i, j, k, 3), u(i, j, k, 4), u(i, j, k, 5),
                           u(i, j, k, 6), u(i, j, k, 7), u(i, j, k, 8), u(i, j, k, 9));
              if (lev > 0) {  // the periodic image, its fine coverage and its parent cell (rank 0's data)
                const int NL = N << lev;
                const ax::IntVect im((i % NL + NL) % NL, (j % NL + NL) % NL, (k % NL + NL) % NL);
                const ax::IntVect cc(im[0] >> 1, im[1] >> 1, im[2] >> 1);
                const bool covered = levels[static_cast<std::size_t>(lev)].ba.contains(im);
                std::fprintf(stderr, "  image (%d,%d,%d) covered by level %d: %d; parent (%d,%d,%d):", im[0], im[1],
                             im[2], lev, covered ? 1 : 0, cc[0], cc[1], cc[2]);
                const ax::MultiFab& Uc = U(lev - 1);
                ax::MFIter::allowMultipleMFIters(true);
                for (ax::MFIter cm(Uc); cm.isValid(); ++cm)
                  if (ax::grow(cm.validbox(), NG).contains(cc)) {
                    auto c = Uc.const_array(cm);
                    const ax::Box gb = ax::grow(cm.validbox(), NG);
                    for (int dk = -2; dk <= 2; ++dk) {
                      const ax::IntVect q(cc[0], cc[1], cc[2] + dk);
                      if (gb.contains(q))
                        std::fprintf(stderr, " [k=%d rho %g E %g e %g]", q[2], c(q[0], q[1], q[2], 0), c(q[0], q[1], q[2], 4),
                                     c(q[0], q[1], q[2], 5));
                    }
                    break;
                  }
                std::fprintf(stderr, "\n  fine column (%d,%d,k): ", i, j);
                for (int kk = -2; kk <= 1; ++kk) std::fprintf(stderr, "[k=%d rho %g] ", kk, u(i, j, kk, 0));
                std::fprintf(stderr, "\n");
              }
              return;
            }
          }
    }
  }
  ax::MultiFab& U(int lev) { return levels[static_cast<std::size_t>(lev)].U; }

  // ---- the gas state on the level (amrex_state = 1): local copies -----------

  void copy_in(const gs::GasState& st) {
    const double t0 = now();
    ax::MultiFab& U0 = U(0);
    for (ax::MFIter mfi(U0); mfi.isValid(); ++mfi) {
      auto u = U0.array(mfi);
      auto g = st.G.const_array(mfi);
      ax::ParallelForOMP(mfi.validbox(), NC, [=] AMREX_GPU_DEVICE(int i, int j, int k, int m) {
        u(i, j, k, m) = g(i, j, k, m);  // the same component order: u(5), e, w(4)
      });
    }
    T.copyin += now() - t0;
  }
  void copy_out(gs::GasState& st) {
    const double t0 = now();
    const ax::MultiFab& U0 = U(0);
    for (ax::MFIter mfi(U0); mfi.isValid(); ++mfi) {
      auto u = U0.const_array(mfi);
      auto g = st.G.array(mfi);
      ax::ParallelForOMP(mfi.validbox(), NC, [=] AMREX_GPU_DEVICE(int i, int j, int k, int m) {
        g(i, j, k, m) = static_cast<float>(u(i, j, k, m));  // rounded once, as the serial step
      });
    }
    T.copyout += now() - t0;
  }

  // The refined levels from and to the state (amrex_state = 1): the state's
  // float copy is the gas between steps; the hydro's double U is refilled
  // at every call (F3a, docs/stage_f.md section 3)
  void copy_in_fine() {
    if (!state) return;
    const double t0 = now();
    for (int lev = 1; lev <= finest; ++lev) {
      ax::MultiFab& Ul = U(lev);
      const gs::FFab& G = state->fine[static_cast<std::size_t>(lev) - 1].G;
      for (ax::MFIter mfi(Ul); mfi.isValid(); ++mfi) {
        auto u = Ul.array(mfi);
        auto g = G.const_array(mfi);
        ax::ParallelForOMP(mfi.validbox(), NC, [=] AMREX_GPU_DEVICE(int i, int j, int k, int m) {
          u(i, j, k, m) = g(i, j, k, m);
        });
      }
    }
    T.copyin += now() - t0;
  }
  void copy_out_fine() {
    if (!state) return;
    const double t0 = now();
    for (int lev = 1; lev <= finest; ++lev) {
      const ax::MultiFab& Ul = U(lev);
      gs::FFab& G = state->fine[static_cast<std::size_t>(lev) - 1].G;
      for (ax::MFIter mfi(Ul); mfi.isValid(); ++mfi) {
        auto u = Ul.const_array(mfi);
        auto g = G.array(mfi);
        ax::ParallelForOMP(mfi.validbox(), NC, [=] AMREX_GPU_DEVICE(int i, int j, int k, int m) {
          g(i, j, k, m) = static_cast<float>(u(i, j, k, m));
        });
      }
    }
    T.copyout += now() - t0;
  }
  // The state's fine level lev (re)defined from the hydro's (ba, dm) and
  // filled from U(lev), rounded
  void state_define_fine(int lev) {
    if (!state) return;
    if (static_cast<int>(state->fine.size()) < lev) state->fine.resize(static_cast<std::size_t>(lev));
    gs::GasState::FineLevel& F = state->fine[static_cast<std::size_t>(lev) - 1];
    const LevelData& L = levels[static_cast<std::size_t>(lev)];
    F.ba = L.ba;
    F.dm = L.dm;
    F.G.define(L.ba, L.dm, NC, 0);
    const ax::MultiFab& Ul = L.U;
    for (ax::MFIter mfi(Ul); mfi.isValid(); ++mfi) {
      auto u = Ul.const_array(mfi);
      auto g = F.G.array(mfi);
      ax::ParallelForOMP(mfi.validbox(), NC, [=] AMREX_GPU_DEVICE(int i, int j, int k, int m) {
        g(i, j, k, m) = static_cast<float>(u(i, j, k, m));
      });
    }
  }

  // ---- host <-> level ------------------------------------------------------

  // The rank's boxes from the host arrays (float -> double)
  void copy_in(const GasGrid& gas) {
    const double t0 = now();
    const float *gu = gas.u.data(), *gw = gas.w.data(), *ge = gas.e.data();
    const int NN = N;
    ax::MultiFab& U0 = U(0);
    for (ax::MFIter mfi(U0); mfi.isValid(); ++mfi) {
      auto u = U0.array(mfi);
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
    const ax::MultiFab& U0 = U(0);
    const ax::BoxArray& ba = levels[0].ba;
    const ax::DistributionMapping& dm = levels[0].dm;
    sendbuf.resize(static_cast<std::size_t>(counts[myproc]));
    for (ax::MFIter mfi(U0); mfi.isValid(); ++mfi) {
      auto u = U0.const_array(mfi);
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

  // ---- the hierarchy (amrex_max_level > 0) -----------------------------------

  // Tag the level's cells whose gas mass exceeds refine_mass base-cell
  // masses: rho/rho_avg > refine_mass * 8^lev
  void tag(int lev, ax::TagBoxArray& tags) {
    const double thr = rho_avg * refine_mass * std::ldexp(1.0, 3 * lev);
    const ax::MultiFab& Ul = U(lev);
    // amrex_refine_box (level-0 cells, inclusive) at this level's resolution:
    // box lev, or the last one for the deeper levels
    const int* rb = refine_box + 6 * std::min(lev, std::max(refine_nbox - 1, 0));
    const bool have_box = refine_nbox > 0 && rb[3] >= rb[0] && rb[4] >= rb[1] && rb[5] >= rb[2];
    const ax::Box fixed(ax::IntVect(rb[0] << lev, rb[1] << lev, rb[2] << lev),
                        ax::IntVect(((rb[3] + 1) << lev) - 1, ((rb[4] + 1) << lev) - 1, ((rb[5] + 1) << lev) - 1));
    for (ax::MFIter mfi(Ul); mfi.isValid(); ++mfi) {
      auto u = Ul.const_array(mfi);
      auto t = tags.array(mfi);
      ax::ParallelForOMP(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        if (u(i, j, k, 0) > thr || (have_box && fixed.contains(ax::IntVect(i, j, k)))) t(i, j, k) = ax::TagBox::SET;
      });
    }
  }

  // FillPatch of a level's MultiFab (ncomp components from scomp of the
  // sources) from its own data and the parent's (one or two times)
  void fill_patch(int lev, ax::MultiFab& mf, int nghost, double t, const ax::Vector<ax::MultiFab*>& cmf,
                  const ax::Vector<double>& ct, const ax::Vector<ax::MultiFab*>& fmf, int scomp, int ncomp) {
    ax::PhysBCFunctNoOp cbc, fbc;
    ax::Vector<double> ft{t};
    ax::FillPatchTwoLevels(mf, ax::IntVect(nghost), t, cmf, ct, fmf, ft, scomp, 0, ncomp, geom_of(lev - 1),
                           geom_of(lev), cbc, 0, fbc, 0, ax::IntVect(kRef), interpolater(), bcs, 0);
  }

  // Regrid levels 1..max_level from the tags (the data where a level
  // existed, the parent interpolated elsewhere)
  void regrid() {
    const double t0 = now();
    int new_finest = 0;
    ax::Vector<ax::BoxArray> new_grids(static_cast<std::size_t>(max_level) + 1);
    mesh->MakeNewGrids(0, 0.0, new_finest, new_grids);
    for (int lev = 1; lev <= new_finest; ++lev) {
      LevelData& L = levels[static_cast<std::size_t>(lev)];
      const ax::BoxArray& nba = new_grids[static_cast<std::size_t>(lev)];
      if (L.defined && nba == L.ba) continue;
      ax::DistributionMapping ndm(nba);
      ax::MultiFab nU(nba, ndm, NC, NG);
      ax::MultiFab& P = U(lev - 1);
      ax::Vector<ax::MultiFab*> cmf{&P};
      ax::Vector<double> ct{0.0};
      if (L.defined) {
        ax::Vector<ax::MultiFab*> fmf{&L.U};
        fill_patch(lev, nU, 0, 0.0, cmf, ct, fmf, 0, NC);
      } else {
        ax::PhysBCFunctNoOp cbc, fbc;
        ax::InterpFromCoarseLevel(nU, ax::IntVect(0), 0.0, P, 0, 0, NC, geom_of(lev - 1), geom_of(lev), cbc, 0,
                                  fbc, 0, ax::IntVect(kRef), interpolater(), bcs, 0);
      }
      L.define(nba, ndm, lev);
      ax::MultiFab::Copy(L.U, nU, 0, 0, NC, 0);
      mesh->SetBoxArray(lev, nba);
      mesh->SetDistributionMap(lev, ndm);
      state_define_fine(lev);
    }
    for (int lev = new_finest + 1; lev <= max_level; ++lev) {
      levels[static_cast<std::size_t>(lev)].clear();
      mesh->ClearBoxArray(lev);
      mesh->ClearDistributionMap(lev);
    }
    if (state) {
      state->fine.resize(static_cast<std::size_t>(new_finest));
      state->finest = new_finest;
    }
    finest = new_finest;
    finest_ever = std::max(finest_ever, finest);
    mesh->SetFinestLevel(finest);
    for (int lev = 1; lev <= finest; ++lev) check(lev, "after regrid", false);
    if (debug && ax::ParallelDescriptor::IOProcessor())
      std::fprintf(stderr, "AMR debug: regrid: finest %d, level boxes:%s\n", finest,
                   [&] {
                     std::string s;
                     for (int lev = 0; lev <= finest; ++lev) {
                       const ax::BoxArray& ba = levels[static_cast<std::size_t>(lev)].ba;
                       const ax::Box mb = ba.minimalBox();
                       s += " " + std::to_string(ba.size()) + "(" + std::to_string(ba.numPts()) + " cells, " +
                            std::to_string(mb.smallEnd(0)) + "," + std::to_string(mb.smallEnd(1)) + "," +
                            std::to_string(mb.smallEnd(2)) + ".." + std::to_string(mb.bigEnd(0)) + "," +
                            std::to_string(mb.bigEnd(1)) + "," + std::to_string(mb.bigEnd(2)) + ")";
                     }
                     return s;
                   }()
                       .c_str());
    T.regrid += now() - t0;
    T.regrids += 1;
  }

  // The changes others made to level 0 since the last copy-out (the RT:
  // w, u(5), e; the PM kick: the momenta and E; the metals) pushed down to
  // the fine cells of every refined level (docs/stage_f.md section 1):
  // density, e and the species by the coarse ratio (additive where the old
  // value is 0), the momenta by the coarse velocity change, E by the fine
  // cell's kinetic and thermal changes. An unchanged coarse cell leaves its
  // fine cells untouched exactly.
  void push_down() {
    if (finest == 0 || !have_prev0) return;
    const double t0 = now();
    const ax::MultiFab& U0 = U(0);
    for (int lev = 1; lev <= finest; ++lev) {
      LevelData& L = levels[static_cast<std::size_t>(lev)];
      const int r = 1 << lev;
      // the parent cells of the level's boxes, old and new, on the level's mapping
      const ax::BoxArray cba = ax::coarsen(L.ba, r);
      ax::MultiFab cnew(cba, L.dm, NC, 0), cold(cba, L.dm, NC, 0);
      cnew.ParallelCopy(U0, 0, 0, NC, 0, 0, period);
      cold.ParallelCopy(prev0, 0, 0, NC, 0, 0, period);
      for (ax::MFIter mfi(L.U); mfi.isValid(); ++mfi) {
        auto u = L.U.array(mfi);
        auto cn = cnew.const_array(mfi);
        auto co = cold.const_array(mfi);
        ax::ParallelForOMP(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
          const int ci = i >> lev, cj = j >> lev, ck = k >> lev;
          bool changed = false;
          for (int m = 0; m < NC; ++m)
            if (cn(ci, cj, ck, m) != co(ci, cj, ck, m)) changed = true;
          if (!changed) return;
          auto ratio = [&](int m, double x) {
            const double o = co(ci, cj, ck, m), n = cn(ci, cj, ck, m);
            if (n == o) return x;
            // the ratio only when the old value carries it (a species that
            // was ~0 and grew: additive, or n/o overflows)
            const double r = o != 0 ? n / o : 0.0;
            return (o != 0 && std::fabs(r) < 1e6) ? x * r : x + (n - o);
          };
          const double rho_old = u(i, j, k, 0);
          const double rho_new = ratio(0, rho_old);
          const double co0 = co(ci, cj, ck, 0), cn0 = cn(ci, cj, ck, 0);
          double p_old[3], p_new[3];
          for (int m = 0; m < 3; ++m) {
            p_old[m] = u(i, j, k, 1 + m);
            const double dv = (co0 > 0 && cn0 > 0) ? cn(ci, cj, ck, 1 + m) / cn0 - co(ci, cj, ck, 1 + m) / co0 : 0.0;
            p_new[m] = p_old[m] + rho_new * dv;
          }
          const double ke_old = rho_old > 0 ? 0.5 * (p_old[0] * p_old[0] + p_old[1] * p_old[1] + p_old[2] * p_old[2]) / rho_old : 0.0;
          const double ke_new = rho_new > 0 ? 0.5 * (p_new[0] * p_new[0] + p_new[1] * p_new[1] + p_new[2] * p_new[2]) / rho_new : 0.0;
          const double e_old = u(i, j, k, NV);
          const double e_new = ratio(NV, e_old);
          u(i, j, k, 0) = rho_new;
          for (int m = 0; m < 3; ++m) u(i, j, k, 1 + m) = p_new[m];
          u(i, j, k, 4) = u(i, j, k, 4) + (ke_new - ke_old) + (e_new - e_old);
          u(i, j, k, NV) = e_new;
          for (int s = 1; s < NS; ++s) u(i, j, k, NV + s) = ratio(NV + s, u(i, j, k, NV + s));
        });
      }
    }
    T.pushdown += now() - t0;
  }

  // prev0 from the state's level 0 (the RT hierarchy averaged the refined
  // levels down itself)
  void level0_is_current() {
    if (max_level == 0 || !state) return;
    for (ax::MFIter mfi(prev0); mfi.isValid(); ++mfi) {
      auto p = prev0.array(mfi);
      auto g = state->G.const_array(mfi);
      ax::ParallelForOMP(mfi.validbox(), NC, [=] AMREX_GPU_DEVICE(int i, int j, int k, int m) {
        p(i, j, k, m) = static_cast<double>(g(i, j, k, m));
      });
    }
    have_prev0 = true;
  }
  void remember_level0() {
    if (max_level == 0) return;
    const ax::MultiFab& U0 = U(0);
    for (ax::MFIter mfi(U0); mfi.isValid(); ++mfi) {
      auto u = U0.const_array(mfi);
      auto p = prev0.array(mfi);
      ax::ParallelForOMP(mfi.validbox(), NC, [=] AMREX_GPU_DEVICE(int i, int j, int k, int m) {
        p(i, j, k, m) = static_cast<double>(static_cast<float>(u(i, j, k, m)));  // what was written out
      });
    }
    have_prev0 = true;
  }

  // ---- the step -------------------------------------------------------------

  // The ghosts of U (and of the children's time interpolation)
  void fill_U(int lev, double t) {
    const double t0 = now();
    LevelData& L = levels[static_cast<std::size_t>(lev)];
    if (lev == 0) {
      L.U.FillBoundary(period);
    } else {
      LevelData& P = levels[static_cast<std::size_t>(lev - 1)];
      ax::MultiFab tmp(L.ba, L.dm, NC, NG);
      ax::Vector<ax::MultiFab*> cmf{&P.Uold, &P.U};
      ax::Vector<double> ct{P.t_old, P.t_new};
      if (!(P.t_new > P.t_old)) {  // no interval yet (before the parent's first substep)
        cmf = {&P.U};
        ct = {t};
      }
      ax::Vector<ax::MultiFab*> fmf{&L.U};
      if (debug) tmp.setVal(-12345.0);  // mark untouched cells
      fill_patch(lev, tmp, NG, t, cmf, ct, fmf, 0, NC);
      if (debug) {
        const int NL = N << lev;
        const ax::Box dom(ax::IntVect(0), ax::IntVect(NL - 1));
        Long n_in = 0, n_out = 0, n_in_bad = 0, n_out_bad = 0, n_untouched_in = 0, n_untouched_out = 0;
        for (ax::MFIter mfi(tmp); mfi.isValid(); ++mfi) {
          auto u = tmp.const_array(mfi);
          const ax::Box gb = ax::grow(mfi.validbox(), NG);
          const auto lo = ax::lbound(gb), hi = ax::ubound(gb);
          for (int k = lo.z; k <= hi.z; ++k)
            for (int j = lo.y; j <= hi.y; ++j)
              for (int i = lo.x; i <= hi.x; ++i) {
                if (mfi.validbox().contains(ax::IntVect(i, j, k))) continue;
                const bool in = dom.contains(ax::IntVect(i, j, k));
                const double r = u(i, j, k, 0);
                (in ? n_in : n_out) += 1;
                if (r == -12345.0) (in ? n_untouched_in : n_untouched_out) += 1;
                else if (!(r > 0)) (in ? n_in_bad : n_out_bad) += 1;
              }
        }
        std::fprintf(stderr, "AMR debug: fill level %d t %g: ghosts inside the domain %lld (bad %lld, untouched %lld), outside %lld (bad %lld, untouched %lld)\n",
                     lev, t, (long long)n_in, (long long)n_in_bad, (long long)n_untouched_in, (long long)n_out, (long long)n_out_bad, (long long)n_untouched_out);
      }
      std::swap(L.U, tmp);
    }
    T.fill += now() - t0;
  }

  // The selected thermal energy density of every valid cell, then its ghosts
  void select(int lev, double t, const rrh::DualEnergy& de) {
    const double t0 = now();
    LevelData& L = levels[static_cast<std::size_t>(lev)];
    for (ax::MFIter mfi(L.U); mfi.isValid(); ++mfi) {
      auto u = L.U.const_array(mfi);
      auto es = L.Es.array(mfi);
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
    if (lev == 0) {
      L.Es.FillBoundary(period);
    } else {
      // the parent's current selection in the coarse/fine ghosts
      LevelData& P = levels[static_cast<std::size_t>(lev - 1)];
      ax::MultiFab tmp(L.ba, L.dm, 1, NG);
      ax::Vector<ax::MultiFab*> cmf{&P.Es};
      ax::Vector<double> ct{t};
      ax::Vector<ax::MultiFab*> fmf{&L.Es};
      fill_patch(lev, tmp, NG, t, cmf, ct, fmf, 0, 1);
      std::swap(L.Es, tmp);
    }
    T.fill += now() - t1;
  }

  // Largest signal speed over the level (max: exact, order-independent)
  double max_signal(int lev) {
    LevelData& L = levels[static_cast<std::size_t>(lev)];
    double smax = 0;
    for (ax::MFIter mfi(L.U); mfi.isValid(); ++mfi) {
      auto u = L.U.const_array(mfi);
      auto es = L.Es.const_array(mfi);
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
    return smax;
  }

  // One substep on the ghost-filled U (Es selected and filled). dtdx is the
  // level-0 substep length: the kernel works in units of its own cell, and
  // a level-l cell (dx = 2^-l) takes dt_0 2^-l, so dt_l / dx_l = dt_0 on
  // every level (velocities and fluxes stay in level-0 cell units).
  Counts substep(int lev, double dtdx, const rrh::DualEnergy& de, double rho_min, double Tmin, bool p_old) {
    LevelData& L = levels[static_cast<std::size_t>(lev)];
    // 1. predictor on the valid box plus one ghost layer
    double t0 = now();
    for (ax::MFIter mfi(L.U); mfi.isValid(); ++mfi) {
      auto u = L.U.const_array(mfi);
      auto es = L.Es.const_array(mfi);
      auto q = L.Q.array(mfi);
      ax::ParallelForOMP(ax::grow(mfi.validbox(), 1), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        auto Ua = [&](int ii, int jj, int kk, int m) { return u(ii, jj, kk, m); };
        auto Ea = [&](int ii, int jj, int kk) { return es(ii, jj, kk); };
        rrh::predict_cell(Ua, Ea, i, j, k, dtdx, q(i, j, k).v);
      });
    }
    double t1 = now();
    T.predict += t1 - t0;
    // 2. the faces of every valid cell (box-local) and the update in place
    for (ax::MFIter mfi(L.U); mfi.isValid(); ++mfi) {
      const ax::Box bx = mfi.validbox();
      auto q = L.Q.const_array(mfi);
      auto u = L.U.array(mfi);
      auto uc = L.U.const_array(mfi);
      auto esf = L.Es.const_array(mfi);
      const bool e_self = de.e_self;
      for (int d = 0; d < 3; ++d) {
        const ax::Box fbx = ax::surroundingNodes(bx, d);
        auto f = L.flux[d].array(mfi);
        const int ei = d == 0, ej = d == 1, ek = d == 2;
        ax::ParallelForOMP(fbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
          auto Qa = [&](int ii, int jj, int kk) -> const double* { return q(ii, jj, kk).v; };
          auto Ua = [&](int ii, int jj, int kk, int m) { return uc(ii, jj, kk, m); };
          auto Sa = [&](int ii, int jj, int kk, int s) { return uc(ii, jj, kk, NV + s); };
          rrh::FaceData fd;
          // e_self: the cells whose internal energy is evolved from its own value
          const bool own_a = e_self && rrh::e_own(uc(i - ei, j - ej, k - ek, NV), esf(i - ei, j - ej, k - ek));
          const bool own_b = e_self && rrh::e_own(uc(i, j, k, NV), esf(i, j, k));
          rrh::face_all(Qa, Ua, Sa, i - ei, j - ej, k - ek, d, dtdx, fd, own_a, own_b);
          put_face(f, i, j, k, fd);
        });
      }
      double t2 = now();
      T.flux += t2 - t1;
      t1 = t2;
      auto fx = L.flux[0].const_array(mfi), fy = L.flux[1].const_array(mfi), fz = L.flux[2].const_array(mfi);
      auto es = L.Es.const_array(mfi);
      auto en = L.En.array(mfi);
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
        // the old time (Bryan et al. 2014 eq. 43); with e_self that of the
        // cell's own internal energy
        const double p_half = (p_old ? (rrh::kGamma - 1) * rrf::max(es(i, j, k), 0.0) : q(i, j, k).v[4]) *
                              (e_self ? rrh::e_ratio(S[0], es(i, j, k)) : 1.0);
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
    if (lev == 0) {
      L.En.FillBoundary(period);
    } else {
      // the parent's updated E in the coarse/fine ghosts
      LevelData& P = levels[static_cast<std::size_t>(lev - 1)];
      ax::MultiFab PE(P.U, ax::make_alias, 4, 1);
      ax::MultiFab tmp(L.ba, L.dm, 1, NG);
      ax::Vector<ax::MultiFab*> cmf{&PE};
      ax::Vector<double> ct{0.0};
      ax::Vector<ax::MultiFab*> fmf{&L.En};
      fill_patch(lev, tmp, NG, 0.0, cmf, ct, fmf, 0, 1);
      std::swap(L.En, tmp);
    }
    double t4 = now();
    T.fill += t4 - t3;
    Counts c;
    for (ax::MFIter mfi(L.U); mfi.isValid(); ++mfi) {
      auto u = L.U.array(mfi);
      auto en = L.En.const_array(mfi);
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

  // Level lev from t to t + dt (its own time step: dt_0 2^-lev), and its
  // children twice by dt/2 (recursive). dtdx = dt 2^lev is the level-0
  // substep, the same on every level (substep above); the flux register
  // holds dt A F / V_c in level-0 cell units: -dtdx from the coarse side
  // and dtdx/8 per fine face and fine substep (dt_f A_f / V_c with
  // dt_f = dt_c/2, A_f = A_c/4, V_c = A_c dx_c).
  void advance(int lev, double t, double dt, const rrh::DualEnergy& de, double rho_min, double Tmin, bool p_old,
               Counts& total_c, bool prepared) {
    LevelData& L = levels[static_cast<std::size_t>(lev)];
    const double dtdx = dt * (1 << lev);
    if (!prepared) {  // the CFL pass filled and selected every level at this time
      fill_U(lev, t);
      select(lev, t, de);
    }
    check(lev, "before substep", true);
    if (lev < finest) {
      ax::MultiFab::Copy(L.Uold, L.U, 0, 0, NC, 0);
      L.t_old = t;
      L.t_new = t + dt;
    }
    const Counts sc = substep(lev, dtdx, de, rho_min, Tmin, p_old);
    check(lev, "after substep", false);
    for (int m = 0; m < 5; ++m) total_c.v[m] += sc.v[m];
    if (lev < finest) {
      LevelData& F = levels[static_cast<std::size_t>(lev) + 1];
      const double t0 = now();
      F.fr->setVal(0.0);
      for (int d = 0; d < 3; ++d) F.fr->CrseInit(L.flux[d], d, 0, 0, NC, -dtdx);
      T.reflux += now() - t0;
      for (int n = 0; n < kRef; ++n)
        advance(lev + 1, t + n * (dt / kRef), dt / kRef, de, rho_min, Tmin, p_old, total_c, prepared && n == 0);
      const double t1 = now();
      // unit volumes: the register holds dt_c A_c F_c and dt_f A_f F_f in
      // level-0 cell units, already divided by V_c (the factors above)
      ax::MultiFab vol(L.ba, L.dm, 1, 0);
      vol.setVal(1.0);
      F.fr->Reflux(L.U, vol, 1.0, 0, 0, NC, geom_of(lev));
      check(lev, "after reflux", false);
      ax::average_down(F.U, L.U, geom_of(lev + 1), geom_of(lev), 0, NC, kRef);
      check(lev, "after average_down", false);
      T.reflux += now() - t1;
    }
    if (lev > 0) {
      // this level's fluxes into the register with its parent: dtdx / 8 per
      // fine face and fine substep (above)
      const double t0 = now();
      for (int d = 0; d < 3; ++d) L.fr->FineAdd(L.flux[d], d, 0, 0, NC, dtdx / 8.0);
      T.reflux += now() - t0;
    }
  }

  void step(const Config& cfg, const Units& un, double dt, double scalefactor, GasGrid& gas, HydroState& hs) {
    events.clear();
    if (state)
      copy_in(*state);
    else
      copy_in(gas);
    copy_in_fine();
    rho_avg = un.Dgas_avg;
    if (max_level > 0) {
      push_down();
      if (T.calls % regrid_every == 0 || finest == 0) regrid();
    }
    const rrh::DualEnergy de{cfg.dual_energy ? cfg.dual_eta1 : 0.0, cfg.dual_eta2, cfg.dual_eta1_local,
                             cfg.dual_nb_radius, cfg.dual_reset_rule, cfg.dual_eta2_local,
                             cfg.dual_e_self && cfg.dual_energy};
    const double rho_min = un.Dgas_min;
    const double Tmin = un.Thydro_min;  // thermal_floor (1 K by default)
    UnsplitCounters& C = hs.ucounters;
    Counts total_c;
    double remaining = dt;
    double t = 0;  // the time within this call: the levels' interpolation intervals start empty
    for (LevelData& L : levels) L.t_old = L.t_new = 0;
    Long nsub_call = 0;
    while (remaining > 0) {
      // the CFL speed over every level (a level-l cell with signal s bounds
      // the level-0 substep by cfl/s too, since dt and dx halve together)
      const double tc = now();
      double smax = 0;
      for (int lev = 0; lev <= finest; ++lev) {
        fill_U(lev, t);
        select(lev, t, de);
        smax = rrf::max(smax, max_signal(lev));
      }
      ax::ParallelDescriptor::ReduceRealMax(smax);
      T.cfl += now() - tc;
      if (!(smax > 0) || !std::isfinite(smax))
        throw std::runtime_error("hydro_step_unsplit (amrex): non-finite or zero signal speed");
      // the sub-cycling arithmetic of hydro_step_unsplit, unchanged
      const double dt_cfl = cfg.unsplit_cfl / smax;
      C.dt_cfl_min = std::min(C.dt_cfl_min, dt_cfl);
      int nsub = static_cast<int>(std::ceil(remaining / dt_cfl - 1e-9));
      if (nsub < 1) nsub = 1;
      const double dt_sub = nsub == 1 ? remaining : remaining / nsub;
      advance(0, t, dt_sub, de, rho_min, Tmin, cfg.dual_p_old, total_c, true);
      ++nsub_call;
      t += dt_sub;
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
    sink(cfg, dt, scalefactor);
    mass0 = U(0).sum(0);  // (MultiFab::sum: valid cells, reduced; order-dependent at 1e-16)
    remember_level0();
    if (state)
      copy_out(*state);
    else
      copy_out(gas);
    copy_out_fine();
  }

  // The star-formation sink (docs/stage_f.md section 2): on every level's
  // cells not covered by a finer level, a cell with rho/rho_avg >
  // sf_overdensity loses f = min(1, sf_efficiency dt/t_ff) of its mass at
  // constant specific quantities (rho, the momenta, e and the species
  // scale by 1 - f; E loses f (KE + e)), t_ff = (3 pi/4)/sqrt(a x_tot) in
  // code time with x_tot = rho/Dgas_avg (the gas traces the total matter;
  // 4 pi G rho_mean = 2a/3 in the PM's units, pm_level.cpp gravity). One
  // event per cell; then the finer levels are averaged down.
  void sink(const Config& cfg, double dt, double a) {
    if (!(cfg.sf_overdensity > 0)) return;
    const double t0 = now();
    const double thr = rho_avg * cfg.sf_overdensity;
    const double eff = cfg.sf_efficiency;
    const double inv_avg = 1.0 / rho_avg;
    for (int lev = finest; lev >= 0; --lev) {
      LevelData& L = levels[static_cast<std::size_t>(lev)];
      std::unique_ptr<ax::iMultiFab> mask;
      if (lev < finest)
        mask = std::make_unique<ax::iMultiFab>(
            ax::makeFineMask(L.ba, L.dm, levels[static_cast<std::size_t>(lev) + 1].ba, ax::IntVect(kRef), 0, 1));
      const double vol = std::ldexp(1.0, -3 * lev);  // level-0 cell volumes
      const double h = std::ldexp(1.0, -lev);        // the cell size
      for (ax::MFIter mfi(L.U); mfi.isValid(); ++mfi) {
        auto u = L.U.array(mfi);
        const ax::Box bx = mfi.validbox();
        const auto lo = ax::lbound(bx), hi = ax::ubound(bx);
        ax::Array4<int const> mk;
        if (mask) mk = mask->const_array(mfi);
        for (int k = lo.z; k <= hi.z; ++k)
          for (int j = lo.y; j <= hi.y; ++j)
            for (int i = lo.x; i <= hi.x; ++i) {
              if (mask && mk(i, j, k)) continue;  // covered by a finer level
              const double rho = u(i, j, k, 0);
              if (!(rho > thr)) continue;
              const double xtot = rho * inv_avg;
              const double tff = (3 * M_PI / 4) / std::sqrt(a * xtot);
              const double f = std::min(1.0, eff * dt / tff);
              if (!(f > 0)) continue;
              const double keep = 1 - f;
              const double ke = 0.5 * (u(i, j, k, 1) * u(i, j, k, 1) + u(i, j, k, 2) * u(i, j, k, 2) +
                                       u(i, j, k, 3) * u(i, j, k, 3)) / rho;
              SinkEvent ev;
              for (int m = 0; m < 3; ++m) {
                ev.x[m] = static_cast<float>(((m == 0 ? i : m == 1 ? j : k) + 0.5) * h);
                ev.v[m] = static_cast<float>(u(i, j, k, 1 + m) / rho);
              }
              ev.mass = static_cast<float>(f * rho * vol);
              u(i, j, k, 4) = u(i, j, k, 4) - f * (ke + u(i, j, k, NV));
              for (int m = 0; m < NC; ++m)
                if (m != 4) u(i, j, k, m) = u(i, j, k, m) * keep;
              events.push_back(ev);
              T.sink_events += 1;
              T.sink_mass += f * rho * vol;
            }
      }
    }
    for (int lev = finest; lev > 0; --lev)
      ax::average_down(U(lev), U(lev - 1), geom_of(lev), geom_of(lev - 1), 0, NC, kRef);
    T.sink += now() - t0;
  }

  // Level 0 onto a new mapping (amrex_balance = 1): its data is copied in
  // from the gas state at every step; prev0 (the push-down reference) moves
  void remap(const ax::DistributionMapping& dm2) {
    if (state && !(dm2 == state->dm)) throw std::runtime_error("amrex_hydro_remap: the gas state has another mapping");
    LevelData& L0 = levels[0];
    const ax::BoxArray ba = L0.ba;
    ax::MultiFab p2;
    if (max_level > 0) {
      p2.define(ba, dm2, NC, 0);
      p2.ParallelCopy(prev0, 0, 0, NC, 0, 0, period);
    }
    L0.clear();
    L0.define(ba, dm2, 0);
    if (max_level > 0) {
      prev0 = std::move(p2);
      mesh->SetDistributionMap(0, dm2);
    }
  }

  // ---- the refined levels on disk (the checkpoint) ---------------------------

  void write_levels(const std::string& dir) const {
    if (max_level == 0) return;
    if (ax::ParallelDescriptor::IOProcessor()) {
      ax::UtilCreateDirectory(dir, 0755);
      std::ofstream h(dir + "/Header");
      h << "radray amr levels 1\n" << finest << "\n";
      for (int lev = 1; lev <= finest; ++lev) {
        const LevelData& L = levels[static_cast<std::size_t>(lev)];
        h << lev << "\n";
        L.ba.writeOn(h);
        h << "\n";
      }
    }
    ax::ParallelDescriptor::Barrier();
    for (int lev = 1; lev <= finest; ++lev) {
      const LevelData& L = levels[static_cast<std::size_t>(lev)];
      if (state) {
        // the state's float copy is the gas between steps: write it (as
        // doubles, VisMF's format), so a restart is bit-identical
        ax::MultiFab tmp(L.ba, L.dm, NC, 0);
        const gs::FFab& G = state->fine[static_cast<std::size_t>(lev) - 1].G;
        for (ax::MFIter mfi(tmp); mfi.isValid(); ++mfi) {
          auto t = tmp.array(mfi);
          auto g = G.const_array(mfi);
          ax::ParallelForOMP(mfi.validbox(), NC, [=] AMREX_GPU_DEVICE(int i, int j, int k, int m) {
            t(i, j, k, m) = g(i, j, k, m);
          });
        }
        ax::VisMF::Write(tmp, dir + "/Level_" + std::to_string(lev) + "_U");
      } else {
        ax::VisMF::Write(L.U, dir + "/Level_" + std::to_string(lev) + "_U");
      }
    }
    // The push-down reference (level 0 as the fine levels last saw it): the
    // changes made to level 0 after it was taken (the summaries' clamps,
    // the halo metals of a step without RT) are pushed down to the fine
    // cells at the next hydro call. A restart that dropped it skipped that
    // push-down, so a restart through such a step was not bit-identical
    // (docs/stage_g.md section 3). Its values are floats (as written out),
    // stored as VisMF doubles: exact.
    if (finest > 0 && have_prev0) ax::VisMF::Write(prev0, dir + "/Level_0_prev");
  }

  void read_levels(const std::string& dir) {
    if (max_level == 0) return;
    std::ifstream h(dir + "/Header");
    if (!h) throw std::runtime_error("amrex_hydro_read_levels: cannot read " + dir + "/Header");
    std::string tag;
    std::getline(h, tag);
    int nf = 0;
    h >> nf;
    if (nf > max_level) throw std::runtime_error("amrex_hydro_read_levels: the checkpoint has more levels than amrex_max_level");
    for (int lev = 1; lev <= nf; ++lev) {
      int l = 0;
      h >> l;
      ax::BoxArray ba;
      ba.readFrom(h);
      LevelData& L = levels[static_cast<std::size_t>(lev)];
      L.define(ba, ax::DistributionMapping(ba), lev);
      ax::VisMF::Read(L.U, dir + "/Level_" + std::to_string(lev) + "_U");
      mesh->SetBoxArray(lev, ba);
      mesh->SetDistributionMap(lev, L.dm);
    }
    for (int lev = nf + 1; lev <= max_level; ++lev) levels[static_cast<std::size_t>(lev)].clear();
    finest = nf;
    if (state) {
      state->fine.clear();
      for (int lev = 1; lev <= nf; ++lev) state_define_fine(lev);  // the floats written above, exactly
      state->finest = nf;
    }
    mesh->SetFinestLevel(finest);
    // the push-down reference of the checkpoint (write_levels); checkpoints
    // without it (before 2026-10-29) restart as before, without the
    // pending push-down
    have_prev0 = false;
    if (nf > 0 && ax::VisMF::Exist(dir + "/Level_0_prev")) {
      ax::MultiFab p;
      ax::VisMF::Read(p, dir + "/Level_0_prev");
      if (p.nComp() != NC || !(p.boxArray().minimalBox() == prev0.boxArray().minimalBox()))
        throw std::runtime_error("amrex_hydro_read_levels: the checkpoint's level-0 reference does not match");
      prev0.ParallelCopy(p, 0, 0, NC);
      have_prev0 = true;
    }
  }

  std::string timers() const {
    char buf[768];
    std::snprintf(buf, sizeof buf,
                  "AMReX hydro timers [s]: copyin %.3f fill %.3f select %.3f cfl %.3f predict %.3f flux %.3f "
                  "update %.3f finalize %.3f reduce %.3f copyout %.3f (bcast %.3f); %lld substeps in %lld "
                  "calls, %d boxes%s",
                  T.copyin, T.fill, T.select, T.cfl, T.predict, T.flux, T.update, T.finalize, T.reduce, T.copyout,
                  T.allgather, static_cast<long long>(T.substeps), static_cast<long long>(T.calls),
                  static_cast<int>(levels[0].ba.size()),
                  max_level > 0 ? (", levels " + std::to_string(finest_ever + 1) + " (now " + std::to_string(finest + 1) +
                                   ") regrid " + std::to_string(T.regrid) +
                                   " reflux " + std::to_string(T.reflux) + " pushdown " + std::to_string(T.pushdown))
                                      .c_str()
                                : "");
    char bm[64];
    std::snprintf(bm, sizeof bm, "; mass0 %.17g", mass0);
    std::string out = std::string(buf) + bm;
    if (T.sink_events > 0) {
      char b2[160];
      std::snprintf(b2, sizeof b2, "; sink %lld events, mass %.6g (this rank, %.3f s)",
                    static_cast<long long>(T.sink_events), T.sink_mass, T.sink);
      out += b2;
    }
    return out;
  }
};

void Mesh::ErrorEst(int lev, ax::TagBoxArray& tags, ax::Real, int) { owner->tag(lev, tags); }

std::unique_ptr<HydroLevel> g_level;

HydroLevel& level_of(const Config& cfg, int N) {
  if (!amrex::Initialized())
    throw std::runtime_error("hydro_backend = amrex: backend_init was not called (AMReX is not initialised)");
  if (!g_level || g_level->N != N) {
    check_hydro_backend_param(cfg);
    const int mgs = cfg.amrex_hydro_max_grid_size ? cfg.amrex_hydro_max_grid_size : cfg.amrex_max_grid_size;
    if (N % 2 != 0 && mgs < N)
      throw std::runtime_error("hydro_backend = amrex: N must be even or amrex_hydro_max_grid_size >= N");
    gs::GasState* st = cfg.amrex_state ? gs::active() : nullptr;
    if (cfg.amrex_state && !st) throw std::runtime_error("amrex_state = 1: amrex_state_init was not called");
    g_level = std::make_unique<HydroLevel>(cfg, N, mgs, st);
  }
  return *g_level;
}

}  // namespace

void amrex_hydro_step(const Config& cfg, const Units& un, double dt, double scalefactor, GasGrid& gas,
                      HydroState& hs) {
  if (!gas.has_e()) throw std::runtime_error("amrex_hydro_step: call init_unsplit_state first");
  level_of(cfg, gas.N).step(cfg, un, dt, scalefactor, gas, hs);
}
void amrex_hydro_take_sink_events(std::vector<SinkEvent>& out) {
  out.clear();
  if (g_level) std::swap(out, g_level->events);
}

std::string amrex_hydro_timers() { return g_level ? g_level->timers() : ""; }
void amrex_hydro_release() { g_level.reset(); }
void amrex_hydro_remap(const amrex::DistributionMapping& dm) {
  if (g_level) g_level->remap(dm);
}
void amrex_hydro_level0_is_current() {
  if (g_level) g_level->level0_is_current();
}
int amrex_hydro_finest_level() { return g_level ? g_level->finest : 0; }
void amrex_hydro_write_levels(const std::string& dir) {
  if (g_level) g_level->write_levels(dir);
}
void amrex_hydro_read_levels(const Config& cfg, int N, const std::string& dir) { level_of(cfg, N).read_levels(dir); }

}  // namespace radray
