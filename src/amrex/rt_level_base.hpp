// The cells of an RT level on the AMReX backend (rt_backend = amrex), shared
// by the ray level (rt_level.cpp, docs/stage_d.md) and the moment level
// (moment_level.cpp, rt_method = m1 | otvet, docs/stage_e_moments.md):
// the decomposition, the cell FabArray (rrk::CellData), the gas and
// sources into the cells (prep), and the transfers between the cells and
// the host's grids (sync_fields; gather_fields / scatter_fields of
// amrex_state = 1). Internal to src/amrex.
#pragma once
#include <algorithm>
#include <array>
#include <chrono>
#include <cstddef>
#include <cstring>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <AMReX.H>
#include <AMReX_FabArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_MFIter.H>
#include <AMReX_ParallelDescriptor.H>

#include "gas_state.hpp"
#include "radray/config.hpp"
#include "radray/kernels/ray_kernels.hpp"
#include "radray/kernels/rt_cell.hpp"
#include "radray/simulation.hpp"

namespace radray {
namespace rtl {

namespace ax = amrex;
using RTCell = rrk::CellData;  // not CellData: amrex::CellData exists too
using rrk::RayP;
using ax::Long;

constexpr int NC = static_cast<int>(sizeof(RTCell) / sizeof(float));
static_assert(sizeof(RTCell) == 46 * sizeof(float), "RTCell must be 46 packed floats");
constexpr int kLuv = offsetof(RTCell, Luv) / sizeof(float);
constexpr int kSfr = offsetof(RTCell, sfr) / sizeof(float);
constexpr int kTau = offsetof(RTCell, tau) / sizeof(float);

using FFab = ax::FabArray<ax::BaseFab<float>>;
using LFab = ax::FabArray<ax::BaseFab<Long>>;
using CFab = ax::FabArray<ax::BaseFab<RTCell>>;

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE float cf(ax::Array4<RTCell const> const& a, int i, int j, int k, int m) {
  float v;
  std::memcpy(&v, reinterpret_cast<const char*>(&a(i, j, k)) + m * sizeof(float), sizeof v);
  return v;
}

inline double now() {
  return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

// scales (bits) of the fixed-point statistics: resolution 2^-bits, range
// 2^(127-bits); the addends are the per-cell values the serial port adds
constexpr int kBitsRate = 100;    // nion, ncol, nrec, nrab, ngam (per H atom, /Ncell): <= 1e8
constexpr int kBitsGamma = 120;   // gamma [1/s] /Ncell: <= 1e2
constexpr int kBitsLambda = 160;  // lambda [erg/s] /Ncell: <= 1e-10
constexpr int kBitsTau = 110;     // dtau: <= 1e5
constexpr int kBitsNgam = 40;     // photons per ray, summed over rays: <= 1e26
// kSPlaneAbs, kSPlaneEsc: plane_parallel = 1, the photons the cells absorbed
// and the rays deleted at their domain's far face carried (the rays'
// units, 2^kBitsNgam fixed point = PlaneState::kFxBits, plane_parallel.hpp)
enum { kSNgam = 0, kSNion = 3, kSNcol = 6, kSNrec = 9, kSNrab = 12, kSGamma = 15, kSLambda = 18, kSDtau = 21,
       kSPlaneAbs = 22, kSPlaneEsc = 23, kSNstat = 24 };
inline int stat_bits(int a) {
  if (a >= kSPlaneAbs) return kBitsNgam;
  if (a >= kSGamma && a < kSLambda) return kBitsGamma;
  if (a >= kSLambda && a < kSDtau) return kBitsLambda;
  if (a == kSDtau) return kBitsTau;
  return kBitsRate;
}

// Per-cell statistics of radiative_physics written by the physics kernel
// and folded into the sums (fixed point; or in cell order, the moment
// level's amrex_moment_stats = exact)
struct CellPhys {
  rrk::PhotoStats st;
  double ngam_src[3][5];  // photo_source addends of a processed source cell
  double plane_abs;       // plane_parallel = 1: deposited - transmitted photons
  int flags;              // 1 processed, 2 photons (ngfill), 4 source (nsfill)
  int err;
};

struct Timers {
  double cell = 0, hist = 0, step1 = 0, cast1 = 0, phys = 0, cast3 = 0, sort = 0, migrate = 0,
         post = 0, prep = 0, sync = 0;
  double fillbnd = 0, sumbnd = 0, zero = 0;
  double bins = 0, gather = 0, writeback = 0, depred = 0;  // parts of sort / cast1
  // the moment level: calc_rt_1, injection, transport (and its halo
  // exchange), deposit, step 3, the exact sums (chains, gathers)
  double rt1 = 0, inject = 0, trans = 0, halo = 0, deposit = 0, absorb = 0, sums = 0, tensor = 0, tflux = 0, tglf = 0, tcopy = 0;
  double c_thin = 0, c_near = 0, c_mip = 0, c_march = 0, c_budget = 0, c_prolong = 0, c_blend = 0;
};

// The level as rt_level.cpp drives it
struct LevelBase {
  virtual ~LevelBase() = default;
  virtual void prep(Simulation& s) = 0;
  virtual void substep(Simulation& s) = 0;
  virtual void post(Simulation& s) = 0;
  virtual void sync_fields(Simulation& s) = 0;
  // amrex_state = 1: the cells into rank 0's mirror (all fields, or the
  // five of write_stats), and rank 0's RT grid onto the cells (a restart)
  virtual void gather_fields(Simulation& s, bool stats_only) = 0;
  // a movie frame: rhog, fHI, fHeI, fHeII, temp and gamma(1) only (the
  // moment level; the ray levels gather every field)
  virtual void gather_movie(Simulation& s) { gather_fields(s, false); }
  virtual void scatter_fields(Simulation& s) = 0;
  // amrex_state_stream_rt (amrex_state.hpp): the float components comps
  // (in that order) of every cell onto rank 0, slab by slab. Collective.
  virtual void stream(const std::vector<int>& comps, const StreamFn& fn) = 0;
  // amrex_state_fill_rt: rank 0's slabs onto the components comps of every
  // cell. Collective.
  virtual void fill(const std::vector<int>& comps, const FillFn& fn) = 0;
  // amrex_balance = 1 (gas_state.hpp): the live rays per box, and the cells
  // and rays moved onto a new mapping of the same BoxArray
  virtual void box_ray_counts(ax::Vector<ax::Real>& count) const = 0;
  virtual void remap(const ax::DistributionMapping& dm2) = 0;
  virtual std::int64_t nrays() const = 0;
  virtual int bytes_per_ray() const = 0;
  virtual void gather_rays(std::vector<RayP>& out) const = 0;
  virtual void gather_rays_chunks(const std::function<void(const std::vector<RayP>&)>& sink) const = 0;
  virtual void scatter_rays(const std::vector<RayP>& in, bool first) = 0;
  virtual std::string timers() const = 0;
  // fields_levels = 1: amrex_rt_level_slice (amrex_state.hpp). Collective.
  virtual void level_slice(int axis, int base_index, std::vector<float>& out, int& nf, int& finest,
                           int& plane) const = 0;
  // the moment level's checkpoint blocks (rt_backend.hpp)
  virtual void ckpt_write(const std::function<void(const void*, std::size_t)>&) const {
    throw std::runtime_error("ckpt_write: not the moment level");
  }
  virtual void ckpt_read(const std::function<void(void*, std::size_t)>&, int) {
    throw std::runtime_error("ckpt_read: not the moment level");
  }
  // gammaHI_avg (Simulation::gamma_avg): the substep-averaged Gamma_HI of
  // every cell onto rank 0 slab by slab (amrex_gamma_avg_stream,
  // amrex_state.hpp). Collective.
  virtual void gamma_avg_stream(const Simulation&, const StreamFn&) const {
    throw std::runtime_error("gammaHI_avg on rt_backend = amrex: implemented for the moment RT only");
  }
};

// Zero slabs of one component in the order of gs::stream_slabs, on rank 0
// only (no communication)
inline void stream_zeros(int N, const StreamFn& fn) {
  if (!ax::ParallelDescriptor::IOProcessor()) return;
  std::vector<float> buf;
  for (int k0 = 0; k0 < N; k0 += gs::kSlabPlanes) {
    StreamSlab sl;
    sl.k0 = k0;
    sl.nk = std::min(gs::kSlabPlanes, N - k0);
    sl.N = N;
    sl.ncomp = 1;
    buf.assign(static_cast<std::size_t>(sl.ncell()), 0.0f);
    sl.data = buf.data();
    fn(sl);
  }
}

// amrex_state_stream_rt for a level's cells C: the components comps
// copied into a float FabArray of the same decomposition, then streamed to
// rank 0 slab by slab (gs::stream_slabs)
inline void stream_cells(const CFab& C, const ax::BoxArray& ba, const ax::DistributionMapping& dm, int N,
                         const std::vector<int>& comps, const StreamFn& fn) {
  const int nc = static_cast<int>(comps.size());
  for (int m : comps)
    if (m < 0 || m >= NC) throw std::runtime_error("amrex_state_stream_rt: component out of range");
  if (nc < 1) throw std::runtime_error("amrex_state_stream_rt: no components");
  FFab Cf(ba, dm, nc, 0);
  std::array<int, NC> cm{};
  for (int m = 0; m < nc && m < NC; ++m) cm[static_cast<std::size_t>(m)] = comps[static_cast<std::size_t>(m)];
  if (nc > NC) throw std::runtime_error("amrex_state_stream_rt: too many components");
  for (ax::MFIter mfi(C); mfi.isValid(); ++mfi) {
    auto c = C.const_array(mfi);
    auto f = Cf.array(mfi);
    ax::ParallelForOMP(mfi.validbox(), nc, [=] AMREX_GPU_DEVICE(int i, int j, int k, int m) {
      f(i, j, k, m) = cf(c, i, j, k, cm[m]);
    });
  }
  std::vector<float> buf;
  gs::stream_slabs(Cf, N, 0, nc, buf, fn);
}

// amrex_state_fill_rt for a level's cells C: rank 0's slabs into a float
// FabArray of the same decomposition (gs::fill_slabs), then into the
// components comps of the cells
inline void fill_cells(CFab& C, const ax::BoxArray& ba, const ax::DistributionMapping& dm, int N,
                       const std::vector<int>& comps, const FillFn& fn) {
  const int nc = static_cast<int>(comps.size());
  if (nc < 1 || nc > NC) throw std::runtime_error("amrex_state_fill_rt: bad component count");
  std::array<int, NC> cm{};
  for (int m = 0; m < nc; ++m) {
    if (comps[static_cast<std::size_t>(m)] < 0 || comps[static_cast<std::size_t>(m)] >= NC)
      throw std::runtime_error("amrex_state_fill_rt: component out of range");
    cm[static_cast<std::size_t>(m)] = comps[static_cast<std::size_t>(m)];
  }
  FFab Cf(ba, dm, nc, 0);
  std::vector<float> buf;
  gs::fill_slabs(Cf, N, 0, nc, buf, fn);
  for (ax::MFIter mfi(C); mfi.isValid(); ++mfi) {
    auto c = C.array(mfi);
    auto f = Cf.const_array(mfi);
    ax::ParallelForOMP(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      float* d = reinterpret_cast<float*>(&c(i, j, k));
      for (int m = 0; m < nc; ++m) d[cm[m]] = f(i, j, k, m);
    });
  }
}

// The cells of a level and their transfers
struct CellLevel : LevelBase {
  int N = 0;
  bool has_e = false;
  gs::GasState* state = nullptr;  // amrex_state = 1 (gas_state.hpp)
  ax::Box domain;
  ax::Geometry geom;
  ax::Periodicity period;
  ax::BoxArray ba;
  ax::DistributionMapping dm;
  CFab C;
  Timers T;

  CellLevel(const Config& cfg, bool with_e) {
    N = cfg.Nrt;
    has_e = with_e;
    domain = ax::Box(ax::IntVect(0), ax::IntVect(N - 1));
    const ax::RealBox rb({0., 0., 0.}, {double(N), double(N), double(N)});  // dx = 1
    const ax::Array<int, 3> per{1, 1, 1};
    geom = ax::Geometry(domain, rb, ax::CoordSys::cartesian, per);
    period = geom.periodicity();
    gs::make_decomposition(cfg, N, ba, dm);
    if (cfg.amrex_state) {
      // the gas state's decomposition (gas_state.hpp): the gas is read and
      // written locally
      state = gs::active();
      if (!state) throw std::runtime_error("amrex_state = 1: amrex_state_init was not called");
      if (state->N != N) throw std::runtime_error("rt_backend = amrex: the gas state has another N");
      ba = state->ba;
      dm = state->dm;
    }
    C.define(ba, dm, 1, 0);
    if (state) {
      // the cells persist between steps (the host arrays start from zero)
      for (ax::MFIter mfi(C); mfi.isValid(); ++mfi) {
        auto c = C.array(mfi);
        ax::ParallelForOMP(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) { c(i, j, k) = RTCell{}; });
      }
    }
  }

  // amrex_state = 1: the gas the RT changed (w, u(5), e) back into the
  // level's state, locally (the host mirror is refreshed on demand,
  // gather_fields)
  void write_gas_back(Simulation& s) {
    const double t1 = now();
    const bool he2 = has_e;
    for (ax::MFIter mfi(C); mfi.isValid(); ++mfi) {
      auto c = C.const_array(mfi);
      auto ga = state->G.array(mfi);
      ax::ParallelForOMP(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const RTCell& d = c(i, j, k);
        for (int m = 0; m < 5; ++m) ga(i, j, k, gs::kU + m) = d.u[m];
        for (int m = 0; m < 4; ++m) ga(i, j, k, gs::kW + m) = d.w[m];
        if (he2) ga(i, j, k, gs::kE) = d.e;
      });
    }
    T.sync += now() - t1;
    s.timers.sync += now() - t1;
  }

  // ---- host state <-> level ------------------------------------------------

  using Box_t = ax::Box;

  void prep(Simulation& s) override {
    const double t0 = now();
    if (has_e != s.gas.has_e()) throw std::runtime_error("rt_backend = amrex: the dual-energy field appeared after the level was built");
    const RTGrid& G = s.rt;
    const GasGrid& g = s.gas;
    const Units u = s.units;
    const bool he = has_e;
    const int NN = N;
    if (state) {
      // the gas from the level's state, the radiation fields of the
      // previous step from the cells themselves, the sources from calc_halo's
      // sparse list (G.src_*, on every rank; ascending cell index)
      const std::int64_t N2 = static_cast<std::int64_t>(NN) * NN;
      for (ax::MFIter mfi(C); mfi.isValid(); ++mfi) {
        auto c = C.array(mfi);
        auto ga = state->G.const_array(mfi);
        const Box_t bx = mfi.validbox();
        ax::ParallelForOMP(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
          RTCell d = c(i, j, k);
          for (int m = 0; m < 5; ++m) d.u[m] = ga(i, j, k, gs::kU + m);
          for (int m = 0; m < 4; ++m) d.w[m] = ga(i, j, k, gs::kW + m);
          d.e = he ? ga(i, j, k, gs::kE) : 0.0f;
          d.Luv = 0.0f;
          d.sfr = 0.0f;
          rrk::prep_cell(d, u, he);
          c(i, j, k) = d;
        });
        // the source cells in the box's k-planes
        const auto b0 = std::lower_bound(G.src_cell.begin(), G.src_cell.end(), bx.smallEnd(2) * N2);
        const auto b1 = std::lower_bound(b0, G.src_cell.end(), (bx.bigEnd(2) + 1) * N2);
        for (auto it = b0; it != b1; ++it) {
          const std::int64_t gc = *it;
          const int i = static_cast<int>(gc % NN), j = static_cast<int>((gc / NN) % NN),
                    k = static_cast<int>(gc / N2);
          if (!bx.contains(ax::IntVect(i, j, k))) continue;
          const std::size_t n = static_cast<std::size_t>(it - G.src_cell.begin());
          c(i, j, k).Luv = G.src_Luv[n];
          c(i, j, k).sfr = G.src_sfr[n];
        }
      }
      T.prep += now() - t0;
      s.timers.prep += now() - t0;
      return;
    }
    // (host pointers: the gas grid is a host array)
    const float *gu = g.u.data(), *gw = g.w.data(), *ge = he ? g.e.data() : nullptr;
    const float *gLuv = G.Luv.data(), *gsfr = G.sfr.data(), *gdx = G.dx.data(), *ggam = G.gamma.data(),
                *glam = G.lambda.data(), *gnu = G.nu.data(), *gng = G.ngam.data(), *gfg = G.fgam.data(),
                *gtau = G.tau.data();
    for (ax::MFIter mfi(C); mfi.isValid(); ++mfi) {
      auto c = C.array(mfi);
      const Box_t bx = mfi.validbox();
      ax::ParallelForOMP(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const std::int64_t gc = i + static_cast<std::int64_t>(NN) * (j + static_cast<std::int64_t>(NN) * k);
        RTCell d{};
        for (int m = 0; m < 5; ++m) d.u[m] = gu[gc * 5 + m];
        for (int m = 0; m < 4; ++m) d.w[m] = gw[gc * 4 + m];
        d.e = he ? ge[gc] : 0.0f;
        d.Luv = gLuv[gc];
        d.sfr = gsfr[gc];
        // the radiation fields of the previous hydro step (calc_rt_1 resets
        // them; kept so sync_fields before the first substep is consistent)
        for (int m = 0; m < 3; ++m) {
          d.dx[m] = gdx[gc * 3 + m];
          d.gamma[m] = ggam[gc * 3 + m];
          d.lambda[m] = glam[gc * 3 + m];
        }
        for (int m = 0; m < 5; ++m) {
          d.nu[m] = gnu[gc * 5 + m];
          d.ngam[m] = gng[gc * 5 + m];
          d.fgam[m] = gfg[gc * 5 + m];
          d.tau[m] = gtau[gc * 5 + m];
        }
        rrk::prep_cell(d, u, he);
        c(i, j, k) = d;
      });
    }
    // the host's RT grid gas state of prep_rt (rt_driver's outputs and the
    // halo/source code read it between the substeps only through sync)
    T.prep += now() - t0;
    s.timers.prep += now() - t0;
  }

  // The level's cells into s.rt and s.gas on every rank: gathered on rank 0
  // and broadcast slab by slab (kSyncPlanes k-planes at a time: a one-box
  // FabArray of the slab on rank 0 by ParallelCopy, flattened, Bcast), so no
  // rank holds a second copy of the whole grid (46 x N^3 floats, 25 GB at
  // 512^3, which with the redundant serial state of every rank overflowed
  // a 1.5 TB node at 24 ranks).
  static constexpr int kSyncPlanes = 8;
  void sync_fields(Simulation& s) override {
    if (state) {  // amrex_state = 1: rank 0's mirror only
      gather_fields(s, false);
      return;
    }
    const double t0 = now();
    FFab Cf(ba, dm, NC, 0);
    for (ax::MFIter mfi(C); mfi.isValid(); ++mfi) {
      auto c = C.const_array(mfi);
      auto f = Cf.array(mfi);
      ax::ParallelForOMP(mfi.validbox(), NC, [=] AMREX_GPU_DEVICE(int i, int j, int k, int m) {
        f(i, j, k, m) = cf(c, i, j, k, m);
      });
    }
    RTGrid& G = s.rt;
    GasGrid& g = s.gas;
    const ax::DistributionMapping dm0(ax::Vector<int>{0});
    const std::size_t nplane = static_cast<std::size_t>(N) * N;
    std::vector<float> flat(nplane * kSyncPlanes * NC);
    for (int k0 = 0; k0 < N; k0 += kSyncPlanes) {
      const int nk = std::min(kSyncPlanes, N - k0);
      const ax::Box slab(ax::IntVect(0, 0, k0), ax::IntVect(N - 1, N - 1, k0 + nk - 1));
      FFab C0(ax::BoxArray(slab), dm0, NC, 0);
      C0.ParallelCopy(Cf, 0, 0, NC);
      if (ax::ParallelDescriptor::IOProcessor()) {
        const auto a = C0.const_array(0);
        for (int k = 0; k < nk; ++k)
          for (int j = 0; j < N; ++j)
            for (int i = 0; i < N; ++i) {
              float* b = &flat[(i + static_cast<std::size_t>(N) * (j + static_cast<std::size_t>(N) * k)) * NC];
              for (int m = 0; m < NC; ++m) b[m] = a(i, j, k0 + k, m);
            }
      }
      if (ax::ParallelDescriptor::NProcs() > 1)
        ax::ParallelDescriptor::Bcast(flat.data(), nplane * nk * NC, ax::ParallelDescriptor::IOProcessorNumber());
      const std::int64_t c0 = static_cast<std::int64_t>(k0) * nplane;
      const std::int64_t c1 = c0 + static_cast<std::int64_t>(nk) * nplane;
#pragma omp parallel for schedule(static)
      for (std::int64_t c = c0; c < c1; ++c) {
        RTCell d;
        std::memcpy(&d, &flat[static_cast<std::size_t>(c - c0) * NC], sizeof d);
        G.rhog[c] = d.rhog;
        G.fHI[c] = d.fHI;
        G.fHeI[c] = d.fHeI;
        G.fHeII[c] = d.fHeII;
        G.temp[c] = d.temp;
        for (int m = 0; m < 3; ++m) {
          G.dx[c * 3 + m] = d.dx[m];
          G.gamma[c * 3 + m] = d.gamma[m];
          G.lambda[c * 3 + m] = d.lambda[m];
        }
        for (int m = 0; m < 5; ++m) {
          G.nu[c * 5 + m] = d.nu[m];
          G.ngam[c * 5 + m] = d.ngam[m];
          G.fgam[c * 5 + m] = d.fgam[m];
          G.tau[c * 5 + m] = d.tau[m];
        }
        for (int m = 0; m < 5; ++m) g.u[c * 5 + m] = d.u[m];
        for (int m = 0; m < 4; ++m) g.w[c * 4 + m] = d.w[m];
        if (has_e) g.e[c] = d.e;
      }
    }
    T.sync += now() - t0;
  }

  // ---- amrex_state = 1: the cells and rank 0's mirror -------------------------

  // The first nc components of every cell (the rt_grid fields, then u, w,
  // e) into rank 0's s.rt and s.gas, slab by slab (gas_state.hpp)
  // A movie frame's fields (rhog, fHI, fHeI, fHeII, temp, gamma(1): the
  // slabs and the whole-box statistics of write_frame) into rank 0's s.rt
  // (the other components of gamma are not refreshed)
  void gather_movie(Simulation& s) override {
    const double t0 = now();
    constexpr int kGamma = offsetof(RTCell, gamma) / sizeof(float);
    FFab Cf(ba, dm, 6, 0);
    for (ax::MFIter mfi(C); mfi.isValid(); ++mfi) {
      auto c = C.const_array(mfi);
      auto f = Cf.array(mfi);
      ax::ParallelForOMP(mfi.validbox(), 6, [=] AMREX_GPU_DEVICE(int i, int j, int k, int m) {
        f(i, j, k, m) = cf(c, i, j, k, m < 5 ? m : kGamma);
      });
    }
    RTGrid& G = s.rt;
    const int NN = N;
    gs::gather_slabs(Cf, N, 6, [&](int k0, int nk, ax::Array4<float const> const& a) {
      if (static_cast<std::int64_t>(G.rhog.size()) != G.ncell)
        throw std::runtime_error("amrex_state = 1: rank 0's RT grid is not allocated");
#pragma omp parallel for schedule(static) collapse(2)
      for (int k = k0; k < k0 + nk; ++k)
        for (int j = 0; j < NN; ++j)
          for (int i = 0; i < NN; ++i) {
            const std::int64_t c = i + static_cast<std::int64_t>(NN) * (j + static_cast<std::int64_t>(NN) * k);
            G.rhog[c] = a(i, j, k, 0);
            G.fHI[c] = a(i, j, k, 1);
            G.fHeI[c] = a(i, j, k, 2);
            G.fHeII[c] = a(i, j, k, 3);
            G.temp[c] = a(i, j, k, 4);
            G.gamma[c * 3] = a(i, j, k, 5);
          }
    });
    T.sync += now() - t0;
    s.timers.sync += now() - t0;
  }

  void stream(const std::vector<int>& comps, const StreamFn& fn) override {
    const double t0 = now();
    stream_cells(C, ba, dm, N, comps, fn);
    T.sync += now() - t0;
  }
  void fill(const std::vector<int>& comps, const FillFn& fn) override {
    const double t0 = now();
    fill_cells(C, ba, dm, N, comps, fn);
    T.sync += now() - t0;
  }

  void gather_fields(Simulation& s, bool stats_only) override {
    const double t0 = now();
    const int nc = stats_only ? 5 : NC;
    FFab Cf(ba, dm, nc, 0);
    for (ax::MFIter mfi(C); mfi.isValid(); ++mfi) {
      auto c = C.const_array(mfi);
      auto f = Cf.array(mfi);
      ax::ParallelForOMP(mfi.validbox(), nc, [=] AMREX_GPU_DEVICE(int i, int j, int k, int m) {
        f(i, j, k, m) = cf(c, i, j, k, m);
      });
    }
    RTGrid& G = s.rt;
    GasGrid& g = s.gas;
    const bool he = has_e;
    const int NN = N;
    gs::gather_slabs(Cf, N, nc, [&](int k0, int nk, ax::Array4<float const> const& a) {
      if (static_cast<std::int64_t>(G.rhog.size()) != G.ncell)
        throw std::runtime_error("amrex_state = 1: rank 0's RT grid is not allocated");
#pragma omp parallel for schedule(static) collapse(2)
      for (int k = k0; k < k0 + nk; ++k)
        for (int j = 0; j < NN; ++j)
          for (int i = 0; i < NN; ++i) {
            const std::int64_t c = i + static_cast<std::int64_t>(NN) * (j + static_cast<std::int64_t>(NN) * k);
            G.rhog[c] = a(i, j, k, 0);
            G.fHI[c] = a(i, j, k, 1);
            G.fHeI[c] = a(i, j, k, 2);
            G.fHeII[c] = a(i, j, k, 3);
            G.temp[c] = a(i, j, k, 4);
            if (stats_only) continue;
            RTCell d;
            float* dp = reinterpret_cast<float*>(&d);
            for (int m = 0; m < NC; ++m) dp[m] = a(i, j, k, m);
            for (int m = 0; m < 3; ++m) {
              G.dx[c * 3 + m] = d.dx[m];
              G.gamma[c * 3 + m] = d.gamma[m];
              G.lambda[c * 3 + m] = d.lambda[m];
            }
            for (int m = 0; m < 5; ++m) {
              G.nu[c * 5 + m] = d.nu[m];
              G.ngam[c * 5 + m] = d.ngam[m];
              G.fgam[c * 5 + m] = d.fgam[m];
              G.tau[c * 5 + m] = d.tau[m];
            }
            for (int m = 0; m < 5; ++m) g.u[c * 5 + m] = d.u[m];
            for (int m = 0; m < 4; ++m) g.w[c * 4 + m] = d.w[m];
            if (he) g.e[c] = d.e;
          }
    });
    T.sync += now() - t0;
    s.timers.sync += now() - t0;
  }

  // Rank 0's rt_grid fields (the first 36 floats of a cell) onto the cells
  // (a restart); u, w, e come from the gas state at the next prep
  void scatter_fields(Simulation& s) override {
    constexpr int nf = kTau + 5;  // rhog .. tau: the rt_grid fields
    static_assert(nf == 36, "rt_grid fields");
    FFab Cf(ba, dm, nf, 0);
    const RTGrid& G = s.rt;
    const int NN = N;
    gs::scatter_slabs(Cf, N, nf, [&](int k0, int nk, ax::Array4<float> const& a) {
      if (static_cast<std::int64_t>(G.rhog.size()) != G.ncell)
        throw std::runtime_error("amrex_state = 1: rank 0's RT grid is not allocated");
#pragma omp parallel for schedule(static) collapse(2)
      for (int k = k0; k < k0 + nk; ++k)
        for (int j = 0; j < NN; ++j)
          for (int i = 0; i < NN; ++i) {
            const std::int64_t c = i + static_cast<std::int64_t>(NN) * (j + static_cast<std::int64_t>(NN) * k);
            RTCell d{};
            d.rhog = G.rhog[c];
            d.fHI = G.fHI[c];
            d.fHeI = G.fHeI[c];
            d.fHeII = G.fHeII[c];
            d.temp = G.temp[c];
            d.Luv = G.Luv[c];
            d.sfr = G.sfr[c];
            for (int m = 0; m < 3; ++m) {
              d.dx[m] = G.dx[c * 3 + m];
              d.gamma[m] = G.gamma[c * 3 + m];
              d.lambda[m] = G.lambda[c * 3 + m];
            }
            for (int m = 0; m < 5; ++m) {
              d.nu[m] = G.nu[c * 5 + m];
              d.ngam[m] = G.ngam[c * 5 + m];
              d.fgam[m] = G.fgam[c * 5 + m];
              d.tau[m] = G.tau[c * 5 + m];
            }
            const float* dp = reinterpret_cast<const float*>(&d);
            for (int m = 0; m < nf; ++m) a(i, j, k, m) = dp[m];
          }
    });
    for (ax::MFIter mfi(C); mfi.isValid(); ++mfi) {
      auto c = C.array(mfi);
      auto f = Cf.const_array(mfi);
      ax::ParallelForOMP(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        float* dp = reinterpret_cast<float*>(&c(i, j, k));
        for (int m = 0; m < nf; ++m) dp[m] = f(i, j, k, m);
      });
    }
  }
};

// The moment level (moment_level.cpp): rt_method = m1 | otvet on the slabs
std::unique_ptr<LevelBase> make_moment_level(const Config& cfg, bool with_e);

}  // namespace rtl
}  // namespace radray
