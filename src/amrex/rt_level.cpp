// The AMReX backend of the ray tracing (rt_backend = amrex; docs/stage_d.md):
// the sorted-ray-array RT substep (docs/design_parallel_and_encoding.md
// section 2) on one uniform AMReX level, MPI + OpenMP, built on the same
// kernels as the serial port (include/radray/kernels/). Provenance:
// spike/amrex/rt_amrex.cpp of the framework spike, turned into a library
// component with the full source model, the pool guard, the RT statistics,
// redshifting and post_rt inside the hydro loop, field synchronisation for
// movie frames, checkpoints and the R-struct ray encoding.
//
//   * cells     FabArray<BaseFab<rrk::CellData>> C (46 floats per cell:
//               rt_grid fields, u_gas, w_gas, e; no ghosts) and
//               FabArray<BaseFab<float>> TF (tau(5), fgam(5); amrex_ghosts
//               ghost layers, FillBoundary before cast_1 / cast_3);
//   * deposits  FabArray<BaseFab<Long>> A (13 int64 fixed-point
//               accumulators, the same ghosts): atomic adds into the box
//               including its ghosts, then SumBoundary (periodic);
//   * rays      ParticleContainerPureSoA (position = the centre of the home
//               cell; the payload in float components, or packed into int
//               components with ray_encoding = rstruct), gathered per box
//               into a cell-sorted rrk::RayP scratch array (DenseBins),
//               migrated by Redistribute() after cast_3;
//   * sums      integer histograms, counters and 128-bit fixed-point sums
//               of the double statistics (exact, associative), reduced
//               with ReduceLongSum (as 32-bit limbs for the 128-bit sums).
//
// Every cross-cell sum is an integer, every merge is in the canonical order
// (cone, payload hash, payload bytes), so the result is bitwise identical
// for any rank count, box size and thread count, and identical to
// spike/serial/rt_serial (the same arithmetic). Against the serial port
// (rt_driver --set legacy_rt=0) the differences are round-off from the
// summation order of deposits and merges (docs/framework_decision.md 2).
//
// Between hydro steps the host's Simulation (rt and gas grids, on every
// rank) holds the cells: amrex_rt_prep copies them in, amrex_rt_post gathers
// them on rank 0 and broadcasts them (the hydro and PM of Stage D run
// redundantly on every rank; Stage E distributes them). The rays live only
// on the level.
#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include <AMReX.H>
#include <AMReX_DenseBins.H>
#include <AMReX_FabArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuAtomic.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Particles.H>
#include <AMReX_Reduce.H>
#include <AMReX_Scan.H>
#ifdef AMREX_USE_OMP
#include <omp.h>
#endif

#include "gas_state.hpp"
#include "radray/amrex_state.hpp"
#include "radray/hydro_backend.hpp"
#include "radray/kernels/ray_kernels.hpp"
#include "radray/kernels/rt_cell.hpp"
#include "radray/ray_tracing.hpp"
#include "radray/rt.hpp"
#include "radray/rt_backend.hpp"
#include "radray/simulation.hpp"
#include "ray_encoding.hpp"

namespace radray {

namespace {

namespace ax = amrex;
using RTCell = rrk::CellData;  // not CellData: amrex::CellData exists too
using rrk::RayP;
using ax::Long;

constexpr int NC = static_cast<int>(sizeof(RTCell) / sizeof(float));
static_assert(sizeof(RTCell) == 46 * sizeof(float), "RTCell must be 46 packed floats");
constexpr int kLuv = offsetof(RTCell, Luv) / sizeof(float);
constexpr int kSfr = offsetof(RTCell, sfr) / sizeof(float);
constexpr int kTau = offsetof(RTCell, tau) / sizeof(float);
static_assert(std::is_same_v<ax::ParticleReal, float>, "build AMReX with AMReX_PARTICLES_PRECISION=SINGLE");
static_assert(sizeof(Long) == sizeof(long long), "int64 accumulators");

using FFab = ax::FabArray<ax::BaseFab<float>>;
using LFab = ax::FabArray<ax::BaseFab<Long>>;
using CFab = ax::FabArray<ax::BaseFab<RTCell>>;

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE float cf(ax::Array4<RTCell const> const& a, int i, int j, int k, int m) {
  float v;
  std::memcpy(&v, reinterpret_cast<const char*>(&a(i, j, k)) + m * sizeof(float), sizeof v);
  return v;
}
// Wrapped global 0-based index g -> index of its periodic image nearest to
// the home index h (h is a valid cell of the box).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE int near_image(int g, int h, int N) {
  int d = g - h;
  if (2 * d > N) d -= N;
  else if (2 * d < -N) d += N;
  return h + d;
}

double now() {
  return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

// ---------------------------------------------------------------------------
// 128-bit fixed point: exact, associative sums of the double statistics
// ---------------------------------------------------------------------------
using i128 = __int128;
inline i128 to_fx128(double x, int bits) {
  return static_cast<i128>(std::floor(std::ldexp(x, bits) + 0.5));
}
inline double from_fx128(i128 v, int bits) { return std::ldexp(static_cast<double>(v), -bits); }
// MPI sum of 128-bit integers as four 32-bit limbs in int64 (each limb sum
// fits for < 2^31 ranks); two's complement makes the limb sum of signed
// values the signed sum modulo 2^128.
void allreduce_i128(i128* v, int n) {
  if (ax::ParallelDescriptor::NProcs() == 1) return;
  std::vector<Long> limbs(static_cast<std::size_t>(4 * n));
  for (int i = 0; i < n; ++i) {
    const unsigned __int128 u = static_cast<unsigned __int128>(v[i]);
    for (int k = 0; k < 4; ++k) limbs[4 * i + k] = static_cast<Long>((u >> (32 * k)) & 0xffffffffULL);
  }
  ax::ParallelDescriptor::ReduceLongSum(limbs.data(), 4 * n);
  for (int i = 0; i < n; ++i) {
    unsigned __int128 r = 0;
    for (int k = 0; k < 4; ++k)
      r += static_cast<unsigned __int128>(static_cast<unsigned long long>(limbs[4 * i + k])) << (32 * k);
    v[i] = static_cast<i128>(r);
  }
}
// scales (bits) of the fixed-point statistics: resolution 2^-bits, range
// 2^(127-bits); the addends are the per-cell values the serial port adds
constexpr int kBitsRate = 100;    // nion, ncol, nrec, nrab, ngam (per H atom, /Ncell): <= 1e8
constexpr int kBitsGamma = 120;   // gamma [1/s] /Ncell: <= 1e2
constexpr int kBitsLambda = 160;  // lambda [erg/s] /Ncell: <= 1e-10
constexpr int kBitsTau = 110;     // dtau: <= 1e5
constexpr int kBitsNgam = 40;     // photons per ray, summed over rays: <= 1e26
enum { kSNgam = 0, kSNion = 3, kSNcol = 6, kSNrec = 9, kSNrab = 12, kSGamma = 15, kSLambda = 18, kSDtau = 21, kSNstat = 22 };
int stat_bits(int a) {
  if (a >= kSGamma && a < kSLambda) return kBitsGamma;
  if (a >= kSLambda && a < kSDtau) return kBitsLambda;
  if (a == kSDtau) return kBitsTau;
  return kBitsRate;
}

// Per-cell statistics of radiative_physics written by the physics kernel
// and folded serially per box into the fixed-point sums
struct CellPhys {
  rrk::PhotoStats st;
  double ngam_src[3][5];  // photo_source addends of a processed source cell
  int flags;              // 1 processed, 2 photons (ngfill), 4 source (nsfill)
  int err;
};

struct Timers {
  double cell = 0, hist = 0, step1 = 0, cast1 = 0, phys = 0, cast3 = 0, sort = 0, migrate = 0,
         post = 0, prep = 0, sync = 0;
  double fillbnd = 0, sumbnd = 0, zero = 0;
  double bins = 0, gather = 0, writeback = 0, depred = 0;  // parts of sort / cast1
};

// Per-box scratch, kept across substeps
struct BoxScratch {
  ax::DenseBins<int> bins;
  ax::Gpu::DeviceVector<int> hcell;              // local home cell of each particle
  ax::Gpu::DeviceVector<RayP> rays, out;         // cell-sorted input / output rays
  ax::Gpu::DeviceVector<int> idx;                // step1_write scratch
  ax::Gpu::DeviceVector<rrk::MergeKey> key;      // step1_write scratch
  ax::Gpu::DeviceVector<Long> cnt, ostart;       // per-cell output counts / offsets
  ax::Gpu::DeviceVector<int> ocell;              // emitting (home) cell of each output ray
  ax::Gpu::DeviceVector<int> keep, wpos;         // cast_3 survivors / compaction offsets
  ax::Gpu::DeviceVector<Long> cmax;              // per-cell cast/split/merge counts (maxima)
  ax::Gpu::DeviceVector<CellPhys> phys;          // per-cell statistics
  ax::Gpu::DeviceVector<Long> depbuf;            // per-thread deposit buffers (CPU)
  Long nout = 0;
};

// A parallel loop over n items with dynamic scheduling on the CPU (the
// per-cell work of step 1 and of the physics is very uneven: source cells
// emit 192 rays, most cells hold none); ParallelFor on the GPU.
template <class F>
void par_for_dyn(int n, F&& f) {
#if defined(AMREX_USE_OMP) && !defined(AMREX_USE_GPU)
#pragma omp parallel for schedule(dynamic, 64)
  for (int i = 0; i < n; ++i) f(i);
#else
  ax::ParallelFor(n, std::forward<F>(f));
#endif
}
#if defined(AMREX_USE_OMP) && !defined(AMREX_USE_GPU)
constexpr bool kThreadBuffers = true;  // per-thread deposit buffers instead of atomics
inline int nthreads_max() { return omp_get_max_threads(); }
inline int thread_id() { return omp_get_thread_num(); }
#else
constexpr bool kThreadBuffers = false;
inline int nthreads_max() { return 1; }
inline int thread_id() { return 0; }
#endif

void zero(ax::Gpu::DeviceVector<Long>& v) {
  Long* p = v.data();
  ax::ParallelFor(static_cast<int>(v.size()), [=] AMREX_GPU_DEVICE(int n) { p[n] = 0; });
}

// The level with one ray encoding (ray_encoding.hpp: FloatEnc, RStructEnc)
struct LevelBase {
  virtual ~LevelBase() = default;
  virtual void prep(Simulation& s) = 0;
  virtual void substep(Simulation& s) = 0;
  virtual void post(Simulation& s) = 0;
  virtual void sync_fields(Simulation& s) = 0;
  // amrex_state = 1: the cells into rank 0's mirror (all fields, or the
  // five of write_stats), and rank 0's RT grid onto the cells (a restart)
  virtual void gather_fields(Simulation& s, bool stats_only) = 0;
  virtual void scatter_fields(Simulation& s) = 0;
  virtual std::int64_t nrays() const = 0;
  virtual int bytes_per_ray() const = 0;
  virtual void gather_rays(std::vector<RayP>& out) const = 0;
  virtual void scatter_rays(const std::vector<RayP>& in) = 0;
  virtual std::string timers() const = 0;
};

template <class Enc>
struct Level : LevelBase {
  using RayPC = ax::ParticleContainerPureSoA<Enc::kNR, Enc::kNI>;
  using PTD = typename RayPC::ParticleTileType::ParticleTileDataType;

  int N = 0, ng = 2, redist_local = 3;
  bool has_e = false;
  gs::GasState* state = nullptr;  // amrex_state = 1 (gas_state.hpp)
  std::int64_t substep_count = 0;  // global RT substep counter (dithering keys)
  ax::Box domain;
  ax::Geometry geom;
  ax::Periodicity period;
  ax::BoxArray ba;
  ax::DistributionMapping dm;
  CFab C;
  FFab TF;
  LFab A;
  std::unique_ptr<RayPC> pc;
  std::vector<BoxScratch> scratch;
  ax::Gpu::DeviceVector<Long> hist;   // ngam_pdf[kNpdfMax] + spl_pdf[kNpdfMax] + npos
  ax::Gpu::DeviceVector<Long> cstat;  // ncst, nspl, nmrg, ndel
  ax::Gpu::DeviceVector<int> err;
  // the atomic-rate tables in device memory (host memory on CPU builds),
  // refreshed when the host tables change (legacy_physics)
  ax::Gpu::DeviceVector<double> d_cool, d_alpha, d_gcol;
  const double* rates_src = nullptr;
  double sig_init[3][5];
  Timers T;
  Long ray_substeps = 0;

  rrk::RatesView rates_view(const Simulation& s) {
    const AtomicRates& R = s.rates;
    if (rates_src != R.coef_cool_.data() || d_cool.size() != R.coef_cool_.size()) {
      d_cool.resize(R.coef_cool_.size());
      d_alpha.resize(R.alpha_rec_.size());
      d_gcol.resize(R.gamma_col_.size());
      ax::Gpu::copy(ax::Gpu::hostToDevice, R.coef_cool_.begin(), R.coef_cool_.end(), d_cool.begin());
      ax::Gpu::copy(ax::Gpu::hostToDevice, R.alpha_rec_.begin(), R.alpha_rec_.end(), d_alpha.begin());
      ax::Gpu::copy(ax::Gpu::hostToDevice, R.gamma_col_.begin(), R.gamma_col_.end(), d_gcol.begin());
      rates_src = R.coef_cool_.data();
    }
    rrk::RatesView v;
    v.cool = d_cool.data();
    v.alpha = d_alpha.data();
    v.gcol = d_gcol.data();
    v.legacy_physics = s.cfg.legacy_physics;
    return v;
  }

  Level(const Config& cfg, bool with_e) {
    N = cfg.Nrt;
    // The reach of one substep is k+1 cells for rt_cells_per_substep = k (a
    // split moves a ray by up to one cell, then it is cast k cells; a ray
    // moves k cells, plus one at a cell face): the ghost layers and the
    // local Redistribute grow with k (docs/rt_multicell.md). k = 1: the
    // parameters as given.
    const int reach = cfg.rt_cells_per_substep + 1;
    ng = std::max(cfg.amrex_ghosts, reach);
    redist_local = cfg.amrex_redist_local > 0 ? std::max(cfg.amrex_redist_local, reach) : 0;
    has_e = with_e;
    domain = ax::Box(ax::IntVect(0), ax::IntVect(N - 1));
    const ax::RealBox rb({0., 0., 0.}, {double(N), double(N), double(N)});  // dx = 1
    const ax::Array<int, 3> per{1, 1, 1};
    geom = ax::Geometry(domain, rb, ax::CoordSys::cartesian, per);
    period = geom.periodicity();
    ba = ax::BoxArray(domain);
    ba.maxSize(cfg.amrex_max_grid_size);
    dm = ax::DistributionMapping(ba);
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
    TF.define(ba, dm, 10, ng);
    A.define(ba, dm, rrk::kNdep, ng);
    pc = std::make_unique<RayPC>(geom, dm, ba);
    scratch.resize(C.local_size());
    hist.resize(2 * rrk::kNpdfMax + 1);
    cstat.resize(4);
    err.resize(1, 0);
    double nu[5];
    for (int n = 0; n < 5; ++n) nu[n] = cp::freq_init[n];
    rrk::photoion_sigma(nu, sig_init);
    if (rrk::pdf_nbins() + 2 > rrk::kNpdfMax) throw std::runtime_error("rrk::kNpdfMax too small");
    for (ax::MFIter mfi(C); mfi.isValid(); ++mfi)
      pc->DefineAndReturnParticleTile(0, mfi.index(), mfi.LocalTileIndex());
  }

  void check_err(const char* where) {
    int e = 0;
    ax::Gpu::copy(ax::Gpu::deviceToHost, err.begin(), err.end(), &e);
    ax::ParallelDescriptor::ReduceIntMax(e);
    if (e)
      throw std::runtime_error(std::string("rt_backend = amrex: a ray reached outside the ghost region "
                                           "or a count mismatch in ") +
                               where + " (raise amrex_ghosts?)");
  }

  // Full store of a ray: payload, valid id, and the AMReX position = centre
  // of the new home cell cell_of(x), taken as the periodic image nearest to
  // the emitting cell h (it may lie just outside the domain): Redistribute
  // applies the periodic shift itself; a pre-wrapped position would look
  // like a jump of ~N cells and break the local Redistribute.
  AMREX_GPU_HOST_DEVICE static void store_ray(PTD const& p, int n, const RayP& r, int N, int cpu,
                                              const int h[3], const typename Enc::Ctx& ctx) {
    Enc::store(p, n, r, ctx);
    for (int d = 0; d < 3; ++d)
      p.m_rdata[d][n] = static_cast<float>(near_image(rrk::cell_of(r.x[d], N), h[d], N)) + 0.5f;
    p.m_idcpu[n] = ax::SetParticleIDandCPU(1, cpu);
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
      // the gas from the level's state, the sources from the host (computed
      // on every rank), the radiation fields of the previous step from the
      // cells themselves
      const float *gLuv = G.Luv.data(), *gsfr = G.sfr.data();
      for (ax::MFIter mfi(C); mfi.isValid(); ++mfi) {
        auto c = C.array(mfi);
        auto ga = state->G.const_array(mfi);
        ax::ParallelForOMP(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
          const std::int64_t gc = i + static_cast<std::int64_t>(NN) * (j + static_cast<std::int64_t>(NN) * k);
          RTCell d = c(i, j, k);
          for (int m = 0; m < 5; ++m) d.u[m] = ga(i, j, k, gs::kU + m);
          for (int m = 0; m < 4; ++m) d.w[m] = ga(i, j, k, gs::kW + m);
          d.e = he ? ga(i, j, k, gs::kE) : 0.0f;
          d.Luv = gLuv[gc];
          d.sfr = gsfr[gc];
          rrk::prep_cell(d, u, he);
          c(i, j, k) = d;
        });
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

  // ---- one RT substep: rt_step_1, rt_step_2, rt_step_3 ---------------------

  void substep(Simulation& s) override {
    s.timers.substeps += 1;
    substep_count += 1;
    const int myproc = ax::ParallelDescriptor::MyProc();
    Long* const hist_p = hist.data();
    Long* const cstat_p = cstat.data();
    int* const err_p = err.data();
    const int NN = N;
    const int ngh = ng;
    typename Enc::Ctx ectx = Enc::make_ctx(N, substep_count);

    // rt_stat(1), ray_stat(1) of this substep
    {
      RTStat& t = s.rt_stat;
      t.dtau = 0;
      t.ngam = {};
      t.nion = {};
      t.ncol = {};
      t.nrec = {};
      t.nrab = {};
      t.gamma = {};
      t.lambda = {};
      t.nsfill = 0;
      t.ngfill = 0;
      RayStat& r = s.rays.stat;
      r.ncst1 = r.ncst2 = r.nmrg1 = r.nmrg2 = r.nspl1 = r.nspl2 = r.ndel1 = 0;
    }
    const double sfr_min = s.units.sfr_min;

    double t0 = now();
    // ---- 1a. calc_rt_1: reset and tau (legacy_rt = 0: B3, B4); tau -> TF ghosts ----
    {
      const Units u = s.units;
      double sig[3][5];
      std::memcpy(sig, sig_init, sizeof sig);
      for (ax::MFIter mfi(C); mfi.isValid(); ++mfi) {
        auto c = C.array(mfi);
        auto tf = TF.array(mfi);
        ax::ParallelForOMP(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
          RTCell& d = c(i, j, k);
          rrk::calc_rt_1_cell(rrk::cell_ref(d), u, sig, false);
          for (int q = 0; q < 5; ++q) tf(i, j, k, q) = d.tau[q];
        });
      }
      const double tf0 = now();
      TF.FillBoundary(0, 5, period);
      T.fillbnd += now() - tf0;
    }
    double t1 = now();
    T.cell += t1 - t0;
    s.timers.calc_rt1 += t1 - t0;

    // ---- gather rays per box into cell-sorted scratch (DenseBins): one
    // decode per ray per substep; the histograms below read the scratch ----
    for (ax::MFIter mfi(C); mfi.isValid(); ++mfi) {
      BoxScratch& S = scratch[mfi.LocalIndex()];
      auto& pt = pc->DefineAndReturnParticleTile(0, mfi.index(), mfi.LocalTileIndex());
      const int np = static_cast<int>(pt.numParticles());
      const PTD p = pt.getParticleTileData();
      const Box_t bx = mfi.validbox();
      const auto lo = ax::lbound(bx), len = ax::length(bx);
      S.hcell.resize(np);
      int* hc = S.hcell.data();
      const double tg0 = now();
      ax::ParallelForOMP(np, [=] AMREX_GPU_DEVICE(int n) {
        const int i = static_cast<int>(p.m_rdata[0][n]) - lo.x;
        const int j = static_cast<int>(p.m_rdata[1][n]) - lo.y;
        const int k = static_cast<int>(p.m_rdata[2][n]) - lo.z;
        if (i < 0 || i >= len.x || j < 0 || j >= len.y || k < 0 || k >= len.z) {
          ax::HostDevice::Atomic::Add(err_p, 1);
          hc[n] = 0;
        } else {
          hc[n] = i + len.x * (j + len.y * k);
        }
      });
#if defined(AMREX_USE_GPU) || !defined(AMREX_USE_OMP)
      S.bins.build(np, S.hcell.data(), static_cast<int>(bx.numPts()), [=] AMREX_GPU_DEVICE(int h) { return h; });
#else
      S.bins.build(ax::BinPolicy::OpenMP, np, S.hcell.data(), static_cast<int>(bx.numPts()),
                   [=](int h) { return h; });
#endif
      const double tg1 = now();
      T.bins += tg1 - tg0;
      S.rays.resize(np);
      RayP* rays = S.rays.data();
      const int* perm = S.bins.permutationPtr();
      ax::ParallelForOMP(np, [=] AMREX_GPU_DEVICE(int n) { rays[n] = Enc::load(p, perm[n], ectx); });
      // the tile is rewritten from the survivors at the end of the substep:
      // release it, so at most two copies of the rays exist at any time
      // (tile + scratch here, scratch + output in step 1, output + tile at
      // the write-back); the arena reuses the freed blocks
      pt.resize(0);
      pt.shrink_to_fit();
      T.gather += now() - tg1;
    }
    check_err("gather");
    double t1b = now();
    T.sort += t1b - t1;
    // ---- 1b. calc_ray_1: histograms of the live rays -> thresholds, pool guard ----
    i128 ngam_sum = 0;
    {
      const int Np = rrk::pdf_nbins();
      zero(hist);
      for (ax::MFIter mfi(C); mfi.isValid(); ++mfi) {
        BoxScratch& S = scratch[mfi.LocalIndex()];
        const int np = static_cast<int>(S.rays.size());
        const RayP* rays = S.rays.data();
        if (kThreadBuffers) {
          // CPU: per-thread histograms (a shared histogram with atomics
          // serialises the threads), folded serially; the exact photon sum
          // per thread in 128 bits
          const int nt = nthreads_max();
          std::vector<Long> h(static_cast<std::size_t>(nt) * (2 * rrk::kNpdfMax + 1), 0);
          std::vector<i128> gs(nt, 0);
#pragma omp parallel
          {
            Long* ht = h.data() + static_cast<std::size_t>(thread_id()) * (2 * rrk::kNpdfMax + 1);
            i128 g = 0;
#pragma omp for schedule(static)
            for (int n = 0; n < np; ++n) {
              const RayP& r = rays[n];
              const double ntot = rrk::sum_ngam(r);
              const int b = rrk::pdf_bin(ntot, Np);
              ht[b] += 1;
              if (r.dr < 0 && r.nside < rayp::nside_max) ht[rrk::kNpdfMax + b] += 1;
              if (ntot > 0) ht[2 * rrk::kNpdfMax] += 1;
              g += to_fx128(ntot, kBitsNgam);
            }
            gs[thread_id()] = g;
          }
          for (int t = 0; t < nt; ++t) {
            for (int b = 0; b < 2 * rrk::kNpdfMax + 1; ++b) hist_p[b] += h[static_cast<std::size_t>(t) * (2 * rrk::kNpdfMax + 1) + b];
            ngam_sum += gs[t];
          }
        } else {
          ax::ParallelForOMP(np, [=] AMREX_GPU_DEVICE(int n) {
            const RayP& r = rays[n];
            const double ntot = rrk::sum_ngam(r);
            const int b = rrk::pdf_bin(ntot, Np);
            ax::HostDevice::Atomic::Add(hist_p + b, Long(1));
            if (r.dr < 0 && r.nside < rayp::nside_max) ax::HostDevice::Atomic::Add(hist_p + rrk::kNpdfMax + b, Long(1));
            if (ntot > 0) ax::HostDevice::Atomic::Add(hist_p + 2 * rrk::kNpdfMax, Long(1));
          });
          // exact sum of the photons (ngam_avg of ray.log), serial per box
          for (int n = 0; n < np; ++n) ngam_sum += to_fx128(rrk::sum_ngam(rays[n]), kBitsNgam);
        }
      }
      std::vector<Long> h(2 * rrk::kNpdfMax + 1);
      ax::Gpu::copy(ax::Gpu::deviceToHost, hist.begin(), hist.end(), h.begin());
      ax::ParallelDescriptor::ReduceLongSum(h.data(), 2 * rrk::kNpdfMax + 1);
      allreduce_i128(&ngam_sum, 1);
      std::vector<std::int64_t> ngam_pdf(h.begin(), h.begin() + rrk::kNpdfMax);
      std::vector<std::int64_t> spl_pdf(h.begin() + rrk::kNpdfMax, h.begin() + 2 * rrk::kNpdfMax);
      std::int64_t L = 0;
      for (int n = 1; n <= Np; ++n) L += ngam_pdf[n];
      const double ngam_avg = h[2 * rrk::kNpdfMax] > 0 ? from_fx128(ngam_sum, kBitsNgam) : 0.0;
      calc_ngam_ray(ngam_pdf.data(), Np, ngam_avg, L > 0 ? L : 1, s.rays.Nray_max, false,
                    s.rays.ngam_ray, s.raylog);
      // F3 pool guard: source cells (stellar) of the level
      std::int64_t nsrc = 0;
      for (ax::MFIter mfi(C); mfi.isValid(); ++mfi) {
        auto c = C.const_array(mfi);
        const Box_t bx = mfi.validbox();
        const auto lo = ax::lbound(bx), hi = ax::ubound(bx);
        for (int k = lo.z; k <= hi.z; ++k)
          for (int j = lo.y; j <= hi.y; ++j)
            for (int i = lo.x; i <= hi.x; ++i)
              if (static_cast<double>(cf(c, i, j, k, kSfr)) > sfr_min) ++nsrc;
      }
      ax::ParallelDescriptor::ReduceLongSum(nsrc);
      pool_guard(spl_pdf.data(), Np, L, nsrc, 0, s.rays.Nray_max, s.rays.ngam_ray, s.rays.npool_limited,
                 s.rays.npool_forced, s.raylog,
                 rrk::nray_from_nside(rrk::emit_nside(s.cfg.rt_cells_per_substep)));
    }
    const rrk::RTP P = rt_params(s);
    double t2 = now();
    T.hist += t2 - t1b;
    s.timers.calc_ray1 += t2 - t1b;

    double t2b = t2;

    // ---- 1d. per-cell emit/classify/split/merge: count, scan, write ----
    zero(cstat);
    const double tz = now();
    A.setVal(0);
    T.zero += now() - tz;
    Long cmax[3] = {0, 0, 0};
    for (ax::MFIter mfi(C); mfi.isValid(); ++mfi) {
      BoxScratch& S = scratch[mfi.LocalIndex()];
      const Box_t bx = mfi.validbox();
      const auto lo = ax::lbound(bx), len = ax::length(bx);
      const int ncb = static_cast<int>(bx.numPts());
      const int np = static_cast<int>(S.rays.size());
      auto c = C.const_array(mfi);
      auto acc = A.array(mfi);
      const int* off = S.bins.offsetsPtr();
      const RayP* rays = S.rays.data();
      S.cnt.resize(ncb + 1);
      S.ostart.resize(ncb + 1);
      S.cmax.resize(3 * ncb);
      Long* cnt = S.cnt.data();
      Long* cm = S.cmax.data();
      par_for_dyn(ncb + 1, [=] AMREX_GPU_DEVICE(int q) {
        if (q == ncb) {
          cnt[q] = 0;
          return;
        }
        const int i = lo.x + q % len.x, j = lo.y + (q / len.x) % len.y, k = lo.z + q / (len.x * len.y);
        float tau[5];
        for (int m = 0; m < 5; ++m) tau[m] = cf(c, i, j, k, kTau + m);
        const bool src = static_cast<double>(cf(c, i, j, k, kSfr)) > sfr_min;
        const bool thick = src && rrk::source_is_thick(tau);
        rrk::SourceSpec spec;
        if (src && !thick) spec = rrk::source_spec(i, j, k, cf(c, i, j, k, kLuv), P);
        if (src && !thick && P.ray_rotate) {
          // B16 fix ray_rotate: the serial code's frame rotation (1-based cell)
          double rot[3][3];
          rrk::frame_rotation(P, i + 1, j + 1, k + 1, 0, rot);
          rrk::set_rotation(spec, rot);
        }
        rrk::CellStats cs;
        cnt[q] = rrk::step1_count(rays + off[q], off[q + 1] - off[q], (src && !thick) ? &spec : nullptr, P, cs);
        cm[3 * q] = cs.ncst;
        cm[3 * q + 1] = cs.nspl;
        cm[3 * q + 2] = cs.nmrg;
        if (cs.ncst) ax::HostDevice::Atomic::Add(cstat_p + 0, Long(cs.ncst));
        if (cs.nspl) ax::HostDevice::Atomic::Add(cstat_p + 1, Long(cs.nspl));
        if (cs.nmrg) ax::HostDevice::Atomic::Add(cstat_p + 2, Long(cs.nmrg));
        if (cs.ndel) ax::HostDevice::Atomic::Add(cstat_p + 3, Long(cs.ndel));
      });
      for (int q = 0; q < ncb; ++q)
        for (int m = 0; m < 3; ++m) cmax[m] = std::max(cmax[m], cm[3 * q + m]);
      S.nout = ax::Scan::ExclusiveSum(ncb + 1, S.cnt.data(), S.ostart.data(), ax::Scan::retSum);
      S.out.resize(S.nout);
      S.ocell.resize(S.nout);
      S.idx.resize(np);
      S.key.resize(np);
      const Long* ost = S.ostart.data();
      RayP* out = S.out.data();
      int* ocell = S.ocell.data();
      int* idx = S.idx.data();
      rrk::MergeKey* key = S.key.data();
      par_for_dyn(ncb, [=] AMREX_GPU_DEVICE(int q) {
        const int i = lo.x + q % len.x, j = lo.y + (q / len.x) % len.y, k = lo.z + q / (len.x * len.y);
        float tau[5];
        for (int m = 0; m < 5; ++m) tau[m] = cf(c, i, j, k, kTau + m);
        const float Luv = cf(c, i, j, k, kLuv);
        const bool src = static_cast<double>(cf(c, i, j, k, kSfr)) > sfr_min;
        const bool thick = src && rrk::source_is_thick(tau);
        if (thick) {
          // optically thick source cell: the photons go straight into the cell
          double v[rrk::kNdep];
          rrk::source_thick_deposit(Luv, P, v);
          for (int m = 0; m < rrk::kNdep; ++m)
            ax::HostDevice::Atomic::Add(&acc(i, j, k, m), Long(rrk::to_fx(v[m], rrk::dep_scale(m))));
        }
        const Long n0 = ost[q], nc = ost[q + 1] - ost[q];
        if (nc == 0) return;
        rrk::SourceSpec spec;
        if (src && !thick) spec = rrk::source_spec(i, j, k, Luv, P);
        if (src && !thick && P.ray_rotate) {
          // B16 fix ray_rotate: the serial code's frame rotation (1-based cell)
          double rot[3][3];
          rrk::frame_rotation(P, i + 1, j + 1, k + 1, 0, rot);
          rrk::set_rotation(spec, rot);
        }
        const int w = rrk::step1_write(rays + off[q], off[q + 1] - off[q], (src && !thick) ? &spec : nullptr, i,
                                       j, k, P, idx + off[q], key + off[q], out + n0);
        if (w != nc) ax::HostDevice::Atomic::Add(err_p, 1);
        for (Long r = n0; r < n0 + nc; ++r) ocell[r] = q;
      });
    }
    check_err("step1");
    // the cell-sorted input and its sort scratch are no longer needed
    for (ax::MFIter mfi(C); mfi.isValid(); ++mfi) {
      BoxScratch& S = scratch[mfi.LocalIndex()];
      for (auto* v : {&S.rays}) {
        v->clear();
        v->shrink_to_fit();
      }
      S.idx.clear();
      S.idx.shrink_to_fit();
      S.key.clear();
      S.key.shrink_to_fit();
      S.hcell.clear();
      S.hcell.shrink_to_fit();
    }
    double t3 = now();
    T.step1 += t3 - t2b;

    // ---- 1d'. cast every output ray once (deposits incl. ghosts) ----
    for (ax::MFIter mfi(C); mfi.isValid(); ++mfi) {
      BoxScratch& S = scratch[mfi.LocalIndex()];
      const Box_t bx = mfi.validbox();
      const Box_t gbx = ax::grow(bx, ngh);
      const auto lo = ax::lbound(bx), len = ax::length(bx), glo = ax::lbound(gbx), ghi = ax::ubound(gbx);
      auto tfa = TF.const_array(mfi);
      auto acc = A.array(mfi);
      const RayP* out = S.out.data();
      const int* ocell = S.ocell.data();
      // CPU: each thread deposits into its own buffer of the grown box
      // (no atomics; adjacent rays target the same cells), the buffers are
      // then summed into A cell by cell (integers: any order)
      const int nt = kThreadBuffers ? nthreads_max() : 1;
      const Long ngc = gbx.numPts();
      const auto glen = ax::length(gbx);
      if (kThreadBuffers) {
        S.depbuf.resize(static_cast<std::size_t>(nt) * ngc * rrk::kNdep);
        std::fill(S.depbuf.begin(), S.depbuf.end(), Long(0));
      }
      Long* const dbuf = S.depbuf.data();
      ax::ParallelForOMP(S.nout, [=] AMREX_GPU_DEVICE(Long r) {
        const int q = ocell[r];
        const int hi = lo.x + q % len.x, hj = lo.y + (q / len.x) % len.y, hk = lo.z + q / (len.x * len.y);
        bool bad = false;
        auto map = [&](int i, int j, int k, int& ii, int& jj, int& kk) {
          ii = near_image(i, hi, NN);
          jj = near_image(j, hj, NN);
          kk = near_image(k, hk, NN);
          if (ii < glo.x || ii > ghi.x || jj < glo.y || jj > ghi.y || kk < glo.z || kk > ghi.z) {
            bad = true;
            ii = hi;
            jj = hj;
            kk = hk;
          }
        };
        float tb[5];
        auto tau = [&](int i, int j, int k) -> const float* {
          int ii, jj, kk;
          map(i, j, k, ii, jj, kk);
          for (int m = 0; m < 5; ++m) tb[m] = tfa(ii, jj, kk, m);
          return tb;
        };
        Long* mine = kThreadBuffers ? dbuf + static_cast<std::size_t>(thread_id()) * ngc * rrk::kNdep : nullptr;
        auto dep = [&](int i, int j, int k, const double* v) {
          int ii, jj, kk;
          map(i, j, k, ii, jj, kk);
          if (kThreadBuffers) {
            Long* d = mine + (static_cast<Long>(ii - glo.x) + glen.x * (static_cast<Long>(jj - glo.y) + glen.y * static_cast<Long>(kk - glo.z))) * rrk::kNdep;
            for (int m = 0; m < rrk::kNdep; ++m) d[m] += rrk::to_fx(v[m], rrk::dep_scale(m));
          } else {
            for (int m = 0; m < rrk::kNdep; ++m)
              ax::HostDevice::Atomic::Add(&acc(ii, jj, kk, m), Long(rrk::to_fx(v[m], rrk::dep_scale(m))));
          }
        };
        rrk::cast_1(out[r], P, tau, dep);
        if (bad) ax::HostDevice::Atomic::Add(err_p, 1);
      });
      if (kThreadBuffers) {
        const double tr0 = now();
        ax::ParallelForOMP(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
          const Long cidx = (static_cast<Long>(i - glo.x) + glen.x * (static_cast<Long>(j - glo.y) + glen.y * static_cast<Long>(k - glo.z))) * rrk::kNdep;
          for (int m = 0; m < rrk::kNdep; ++m) {
            Long sum = 0;
            for (int t = 0; t < nt; ++t) sum += dbuf[static_cast<std::size_t>(t) * ngc * rrk::kNdep + cidx + m];
            acc(i, j, k, m) += sum;
          }
        });
        T.depred += now() - tr0;
      }
      ray_substeps += S.nout;
    }
    const double tsb = now();
    A.SumBoundary(period);  // fold ghost deposits into their owners (int64: exact)
    T.sumbnd += now() - tsb;
    check_err("cast_1");
    double t4 = now();
    T.cast1 += t4 - t3;
    s.timers.trace1 += t4 - t2;

    // ---- 2. radiative_physics; statistics; fgam -> TF ghosts ----
    i128 stats[kSNstat];
    for (int a = 0; a < kSNstat; ++a) stats[a] = 0;
    Long nsfill = 0, ngfill = 0;
    {
      rrk::PhysCtx X;
      X.u = s.units;
      X.rates = rates_view(s);
      X.dtsec_rt = s.clk.dtsec_rt;
      X.cells = s.cfg.rt_cells_per_substep;
      X.photon_conserve = s.cfg.photon_conserve;
      X.ia = s.cfg.recomb_case_A ? 0 : 3;
      X.want_nrab = false;
      const bool he = has_e;
      const double Ncell = s.cfg.ncellr;
      const double nrtcr = s.cfg.nrtcr;
      const Units u = s.units;
      const double dtsec_rt = s.clk.dtsec_rt;
      for (ax::MFIter mfi(C); mfi.isValid(); ++mfi) {
        BoxScratch& S = scratch[mfi.LocalIndex()];
        const Box_t bx = mfi.validbox();
        const int ncb = static_cast<int>(bx.numPts());
        const auto lo = ax::lbound(bx), len = ax::length(bx);
        auto c = C.array(mfi);
        auto tf = TF.array(mfi);
        auto acc = A.const_array(mfi);
        S.phys.resize(ncb);
        CellPhys* ph = S.phys.data();
        par_for_dyn(ncb, [=] AMREX_GPU_DEVICE(int q) {
          const int i = lo.x + q % len.x, j = lo.y + (q / len.x) % len.y, k = lo.z + q / (len.x * len.y);
          RTCell& d = c(i, j, k);
          for (int m = 0; m < 5; ++m) d.ngam[m] = static_cast<float>(rrk::from_fx(acc(i, j, k, m), rrk::kScaleNgam));
          for (int m = 0; m < 5; ++m) d.nu[m] = static_cast<float>(rrk::from_fx(acc(i, j, k, 5 + m), rrk::kScaleNu));
          for (int m = 0; m < 3; ++m) d.dx[m] = static_cast<float>(rrk::from_fx(acc(i, j, k, 10 + m), rrk::kScaleNgam));
          CellPhys& cp_ = ph[q];
          cp_ = CellPhys{};
          bool processed = false;
          double ngam_tot = 0;
          cp_.err = rrk::physics_cell(d, X, &cp_.st, he, &processed, &ngam_tot);
          if (processed) {
            cp_.flags |= 1;
            if (ngam_tot > u.ngam_min) cp_.flags |= 2;
            if (static_cast<double>(d.sfr) > u.sfr_min) {
              cp_.flags |= 4;
              rrk::source_stats(d.Luv, dtsec_rt, u, cp_.ngam_src);
            }
          }
          for (int m = 0; m < 5; ++m) tf(i, j, k, 5 + m) = d.fgam[m];
        });
        // fold the cell statistics into the fixed-point sums (serial per
        // box: order independent, so any order will do)
        for (int q = 0; q < ncb; ++q) {
          const CellPhys& cp_ = ph[q];
          if (cp_.err) throw std::runtime_error("radiative_physics: temperature index out of table");
          if (!(cp_.flags & 1)) continue;
          if (cp_.flags & 2) ngfill += 1;
          if (cp_.flags & 4) {
            nsfill += 1;
            for (int inu = 0; inu < 5; ++inu)
              for (int sp = 0; sp < 3; ++sp)
                if (cp_.ngam_src[sp][inu] != 0) stats[kSNgam + sp] += to_fx128(cp_.ngam_src[sp][inu] / nrtcr, kBitsRate);
          }
          const double navg[3] = {u.nH_avg, u.nHe_avg, u.nHe_avg};
          const rrk::PhotoStats& st = cp_.st;
          for (int sp = 0; sp < 3; ++sp) {
            stats[kSNion + sp] += to_fx128(st.nion[sp] / navg[sp] / Ncell, kBitsRate);
            stats[kSNcol + sp] += to_fx128(st.ncol[sp] / navg[sp] / Ncell, kBitsRate);
            stats[kSNrec + sp] += to_fx128(st.nrec[sp] / navg[sp] / Ncell, kBitsRate);
            stats[kSNrab + sp] += to_fx128(st.nrab[sp] / navg[sp] / Ncell, kBitsRate);
            stats[kSGamma + sp] += to_fx128(st.gam[sp] / Ncell, kBitsGamma);
            stats[kSLambda + sp] += to_fx128(st.lam[sp] / Ncell, kBitsLambda);
          }
          stats[kSDtau] += to_fx128(st.dtau / Ncell, kBitsTau);
        }
      }
      const double tf0 = now();
      TF.FillBoundary(5, 5, period);
      T.fillbnd += now() - tf0;
    }
    double t5 = now();
    T.phys += t5 - t4;
    s.timers.physics += t5 - t4;

    // ---- 3. advance (cast_3) ----
    Long nout_loc = 0, nkeep_loc = 0;
    for (ax::MFIter mfi(C); mfi.isValid(); ++mfi) {
      BoxScratch& S = scratch[mfi.LocalIndex()];
      const Box_t bx = mfi.validbox();
      const Box_t gbx = ax::grow(bx, ngh);
      const auto lo = ax::lbound(bx), len = ax::length(bx), glo = ax::lbound(gbx), ghi = ax::ubound(gbx);
      auto tfa = TF.const_array(mfi);
      RayP* out = S.out.data();
      const int* ocell = S.ocell.data();
      S.keep.resize(S.nout + 1);
      S.wpos.resize(S.nout + 1);
      int* keep = S.keep.data();
      const Long nout = S.nout;
      ax::ParallelForOMP(nout + 1, [=] AMREX_GPU_DEVICE(Long r) {
        if (r == nout) {
          keep[r] = 0;
          return;
        }
        const int q = ocell[r];
        const int hi = lo.x + q % len.x, hj = lo.y + (q / len.x) % len.y, hk = lo.z + q / (len.x * len.y);
        bool bad = false;
        float tb[5], fb[5];
        auto tf = [&](int i, int j, int k, const float*& tau, const float*& fgam) {
          int ii = near_image(i, hi, NN), jj = near_image(j, hj, NN), kk = near_image(k, hk, NN);
          if (ii < glo.x || ii > ghi.x || jj < glo.y || jj > ghi.y || kk < glo.z || kk > ghi.z) {
            bad = true;
            ii = hi;
            jj = hj;
            kk = hk;
          }
          for (int m = 0; m < 5; ++m) {
            tb[m] = tfa(ii, jj, kk, m);
            fb[m] = tfa(ii, jj, kk, 5 + m);
          }
          tau = tb;
          fgam = fb;
        };
        keep[r] = rrk::cast_3(out[r], P, tf) ? 1 : 0;
        if (bad) ax::HostDevice::Atomic::Add(err_p, 1);
      });
      nout_loc += nout;
    }
    check_err("cast_3");
    double t6 = now();
    T.cast3 += t6 - t5;
    // ---- write survivors back into the particle tiles ----
    for (ax::MFIter mfi(C); mfi.isValid(); ++mfi) {
      BoxScratch& S = scratch[mfi.LocalIndex()];
      const Long nout = S.nout;
      const Long nkeep = ax::Scan::ExclusiveSum(nout + 1, S.keep.data(), S.wpos.data(), ax::Scan::retSum);
      nkeep_loc += nkeep;
      auto& pt = pc->DefineAndReturnParticleTile(0, mfi.index(), mfi.LocalTileIndex());
      pt.resize(nkeep);
      const PTD p = pt.getParticleTileData();
      const int* keep = S.keep.data();
      const int* wpos = S.wpos.data();
      const RayP* out = S.out.data();
      const int* ocell = S.ocell.data();
      const auto lo = ax::lbound(mfi.validbox()), len = ax::length(mfi.validbox());
      ax::ParallelForOMP(nout, [=] AMREX_GPU_DEVICE(Long r) {
        if (!keep[r]) return;
        const int q = ocell[r];
        const int h[3] = {lo.x + q % len.x, lo.y + (q / len.x) % len.y, lo.z + q / (len.x * len.y)};
        store_ray(p, wpos[r], out[r], NN, myproc, h, ectx);
      });
    }
    // the output rays are in the tiles now
    for (ax::MFIter mfi(C); mfi.isValid(); ++mfi) {
      BoxScratch& S = scratch[mfi.LocalIndex()];
      S.out.clear();
      S.out.shrink_to_fit();
      S.ocell.clear();
      S.ocell.shrink_to_fit();
      S.keep.clear();
      S.keep.shrink_to_fit();
      S.wpos.clear();
      S.wpos.shrink_to_fit();
    }
    double t7 = now();
    T.sort += t7 - t6;
    T.writeback += t7 - t6;
    // ---- migrate rays to the owners of their new home cells ----
    if (redist_local > 0)
      pc->Redistribute(0, 0, 0, redist_local);
    else
      pc->Redistribute();
    double t8 = now();
    T.migrate += t8 - t7;
    s.timers.trace3 += t8 - t5;

    // ---- the substep's statistics into s (every rank the same) ----
    {
      Long sv[9];
      ax::Gpu::copy(ax::Gpu::deviceToHost, cstat.begin(), cstat.end(), sv);
      sv[4] = nout_loc - nkeep_loc;  // deleted in cast_3
      sv[5] = nkeep_loc;             // live rays
      sv[6] = nsfill;
      sv[7] = ngfill;
      sv[8] = 0;
      ax::ParallelDescriptor::ReduceLongSum(sv, 8);
      ax::ParallelDescriptor::ReduceLongMax(cmax, 3);
      allreduce_i128(stats, kSNstat);
      RayStat& r = s.rays.stat;
      r.ncst1 = sv[0];
      r.nspl1 = sv[1];
      r.nmrg1 = sv[2];
      r.ndel1 = sv[3] + sv[4];
      r.ncst2 = cmax[0];
      r.nspl2 = cmax[1];
      r.nmrg2 = cmax[2];
      r.nfree = 0;
      s.rays.Nray = sv[5];
      s.rays.save.nfree = 0;
      RTStat& t = s.rt_stat;
      t.nsfill = sv[6];
      t.ngfill = sv[7];
      for (int sp = 0; sp < 3; ++sp) {
        t.ngam[sp] = from_fx128(stats[kSNgam + sp], kBitsRate);
        t.nion[sp] = from_fx128(stats[kSNion + sp], kBitsRate);
        t.ncol[sp] = from_fx128(stats[kSNcol + sp], kBitsRate);
        t.nrec[sp] = from_fx128(stats[kSNrec + sp], kBitsRate);
        t.nrab[sp] = from_fx128(stats[kSNrab + sp], kBitsRate);
        t.gamma[sp] = from_fx128(stats[kSGamma + sp], kBitsGamma);
        t.lambda[sp] = from_fx128(stats[kSLambda + sp], kBitsLambda);
      }
      t.dtau = from_fx128(stats[kSDtau], kBitsTau);
      accumulate_rt_stats(s);
    }
    s.timers.sync += now() - t8;
  }

  // ---- post_rt: redshift, calc_rt_f, photon total, fields back to s ------

  void post(Simulation& s) override {
    const double t0 = now();
    const RedshiftTable RT = make_redshift_table();
    // (a - da)/a, or per distance with redshift_per_distance (cosmology.hpp)
    const double fscale = redshift_fscale(s.cfg, s.clk);
    const double ngam_min = s.units.ngam_min;
    const Units u = s.units;
    const bool he = has_e;
    double sig[3][5];
    std::memcpy(sig, sig_init, sizeof sig);
    typename Enc::Ctx ectx = Enc::make_ctx(N, substep_count);
    i128 tot = 0;
    Long nr = 0;
    for (ax::MFIter mfi(C); mfi.isValid(); ++mfi) {
      auto& pt = pc->DefineAndReturnParticleTile(0, mfi.index(), mfi.LocalTileIndex());
      const int np = static_cast<int>(pt.numParticles());
      const PTD p = pt.getParticleTileData();
      ax::ParallelForOMP(np, [=] AMREX_GPU_DEVICE(int n) {
        RayP r = Enc::load(p, n, ectx);
        rrk::redshift_ray(r, RT, fscale, ngam_min);
        Enc::store(p, n, r, ectx);
      });
      for (int n = 0; n < np; ++n) {
        const RayP r = Enc::load(p, n, ectx);
        // each group in fixed point (as spike/serial rt_serial: exact sum)
        for (int m = 0; m < 5; ++m) tot += to_fx128(static_cast<double>(r.ngam[m]), kBitsNgam);
      }
      nr += np;
      auto c = C.array(mfi);
      ax::ParallelForOMP(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        rrk::post_cell(c(i, j, k), u, sig, he);
      });
    }
    allreduce_i128(&tot, 1);
    ax::ParallelDescriptor::ReduceLongSum(nr);
    s.ngam_total_backend = from_fx128(tot, kBitsNgam);
    s.rays.Nray = nr;
    s.rays.save.nfree = 0;
    T.post += now() - t0;
    s.timers.post += now() - t0;
    if (state) {
      // the gas the RT changed (w, u(5), e) back into the level's state,
      // locally; the host mirror is refreshed on demand (gather_fields)
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
      return;
    }
    sync_fields(s);
  }

  // ---- amrex_state = 1: the cells and rank 0's mirror -------------------------

  // The first nc components of every cell (the rt_grid fields, then u, w,
  // e) into rank 0's s.rt and s.gas, slab by slab (gas_state.hpp)
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

  std::int64_t nrays() const override {
    Long nr = 0;
    for (ax::MFIter mfi(C); mfi.isValid(); ++mfi)
      nr += pc->DefineAndReturnParticleTile(0, mfi.index(), mfi.LocalTileIndex()).numParticles();
    ax::ParallelDescriptor::ReduceLongSum(nr);
    return nr;
  }
  int bytes_per_ray() const override { return Enc::kBytes; }

  // All rays in the float layout, gathered on rank 0: every rank sends its
  // list to rank 0 in chunks of kGatherChunk rays (an MPI_Gatherv would be
  // limited to 2 GB, 25 M rays; a 512^3 production run holds 1.1 G)
  static constexpr std::size_t kGatherChunk = std::size_t(1) << 23;  // 8 M rays, 640 MB
  void gather_rays(std::vector<RayP>& out) const override {
    typename Enc::Ctx ectx = Enc::make_ctx(N, substep_count);
    std::vector<RayP> mine;
    for (ax::MFIter mfi(C); mfi.isValid(); ++mfi) {
      auto& pt = pc->DefineAndReturnParticleTile(0, mfi.index(), mfi.LocalTileIndex());
      const int np = static_cast<int>(pt.numParticles());
      const PTD p = pt.getParticleTileData();
      for (int n = 0; n < np; ++n) mine.push_back(Enc::load(p, n, ectx));
    }
    const int nranks = ax::ParallelDescriptor::NProcs();
    out.clear();
    if (nranks == 1) {
      out.swap(mine);
      return;
    }
#ifdef AMREX_USE_MPI
    const MPI_Comm comm = ax::ParallelDescriptor::Communicator();
    const int root = ax::ParallelDescriptor::IOProcessorNumber();
    const int myproc = ax::ParallelDescriptor::MyProc();
    const bool io = ax::ParallelDescriptor::IOProcessor();
    long long count = static_cast<long long>(mine.size());
    std::vector<long long> counts(nranks);
    MPI_Gather(&count, 1, MPI_LONG_LONG, counts.data(), 1, MPI_LONG_LONG, root, comm);
    if (io) {
      std::size_t total = 0;
      for (long long c : counts) total += static_cast<std::size_t>(c);
      out.reserve(total);
      out.insert(out.end(), mine.begin(), mine.end());  // rank 0's own first
    }
    for (int r = 0; r < nranks; ++r) {
      if (r == root) continue;
      if (myproc == r) {
        for (std::size_t off = 0; off < mine.size(); off += kGatherChunk) {
          const int n = static_cast<int>(std::min(kGatherChunk, mine.size() - off));
          MPI_Send(mine.data() + off, n * static_cast<int>(sizeof(RayP)), MPI_BYTE, root, 11, comm);
        }
      } else if (io) {
        const std::size_t total = static_cast<std::size_t>(counts[r]);
        for (std::size_t off = 0; off < total; off += kGatherChunk) {
          const int n = static_cast<int>(std::min(kGatherChunk, total - off));
          const std::size_t n0 = out.size();
          out.resize(n0 + static_cast<std::size_t>(n));
          MPI_Recv(out.data() + n0, n * static_cast<int>(sizeof(RayP)), MPI_BYTE, r, 11, comm, MPI_STATUS_IGNORE);
        }
      }
    }
#endif
  }

  // Rank 0 adds the rays (home-cell centres as positions), a global
  // Redistribute moves them to their owners
  // Rank 0 holds the whole list (a checkpoint): the rays go to their owners
  // in chunks of kScatterChunk, each written into rank 0's first tile and
  // redistributed, so rank 0's peak holds one chunk of particles and its
  // migration buffers, not the whole pool (1.1 G rays = 110 GB of
  // particles plus their send buffers at 512^3). Every rank calls
  // Redistribute the same number of times (collective).
  static constexpr std::size_t kScatterChunk = std::size_t(1) << 25;
  void scatter_rays(const std::vector<RayP>& in) override {
    typename Enc::Ctx ectx = Enc::make_ctx(N, substep_count);
    for (ax::MFIter mfi(C); mfi.isValid(); ++mfi)
      pc->DefineAndReturnParticleTile(0, mfi.index(), mfi.LocalTileIndex()).resize(0);
    const bool root = ax::ParallelDescriptor::IOProcessor();
    ax::Long nchunk = root ? static_cast<ax::Long>((in.size() + kScatterChunk - 1) / kScatterChunk) : 0;
    ax::ParallelDescriptor::ReduceLongMax(nchunk);
    if (nchunk == 0) nchunk = 1;
    for (ax::Long ic = 0; ic < nchunk; ++ic) {
      if (root && !in.empty()) {
        ax::MFIter mfi(C);
        if (!mfi.isValid()) throw std::runtime_error("checkpoint: rank 0 owns no box");
        auto& pt = pc->DefineAndReturnParticleTile(0, mfi.index(), mfi.LocalTileIndex());
        const std::size_t n0 = static_cast<std::size_t>(ic) * kScatterChunk;
        const std::size_t n1 = std::min(in.size(), n0 + kScatterChunk);
        const int old = pt.numParticles();
        pt.resize(old + (n1 - n0));
        const PTD p = pt.getParticleTileData();
        const int myproc = ax::ParallelDescriptor::MyProc();
        const int NN = N;
        for (std::size_t n = n0; n < n1; ++n) {
          const int h[3] = {rrk::cell_of(in[n].x[0], NN), rrk::cell_of(in[n].x[1], NN),
                            rrk::cell_of(in[n].x[2], NN)};
          store_ray(p, old + static_cast<int>(n - n0), in[n], NN, myproc, h, ectx);
        }
      }
      pc->Redistribute();
    }
  }

  std::string timers() const override {
    char buf[512];
    std::snprintf(buf, sizeof buf,
                  "AMReX RT timers [s]: prep %.3f cell %.3f hist %.3f step1 %.3f cast1 %.3f phys %.3f "
                  "cast3 %.3f sort %.3f migrate %.3f post %.3f sync %.3f (FillBoundary %.3f SumBoundary %.3f "
                  "setVal %.3f; sort: bins %.3f gather %.3f writeback %.3f; cast1: deposit reduce %.3f); "
                  "ray-substeps %lld",
                  T.prep, T.cell, T.hist, T.step1, T.cast1, T.phys, T.cast3, T.sort, T.migrate, T.post, T.sync,
                  T.fillbnd, T.sumbnd, T.zero, T.bins, T.gather, T.writeback, T.depred,
                  static_cast<long long>(ray_substeps));
    return buf;
  }
};

}  // namespace

struct AmrexRT {
  std::unique_ptr<LevelBase> level;
};

namespace {
bool g_initialized = false;

LevelBase& level_of(Simulation& s) {
  if (!s.amrex) {
    check_rt_backend_param(s.cfg);
    s.amrex = std::make_shared<AmrexRT>();
    const bool with_e = s.gas.has_e();
    if (s.cfg.ray_encoding == RayEncoding::RStruct)
      s.amrex->level = std::make_unique<Level<amrex_enc::RStructEnc>>(s.cfg, with_e);
    else
      s.amrex->level = std::make_unique<Level<amrex_enc::FloatEnc>>(s.cfg, with_e);
  }
  return *s.amrex->level;
}
}  // namespace

bool amrex_available() { return true; }

void backend_init(int& argc, char**& argv) {
  if (g_initialized) return;
  amrex::Initialize(argc, argv, false, MPI_COMM_WORLD, []() {
    amrex::ParmParse pp("amrex");
    pp.add("verbose", 0);
    pp.add("the_arena_is_managed", 0);
  });
  g_initialized = true;
}
void backend_finalize() {
  if (!g_initialized) return;
  amrex_hydro_release();  // the hydro level's FabArrays before the arenas go
  amrex_state_release();
  amrex::Finalize();
  g_initialized = false;
}
bool io_rank() { return !g_initialized || amrex::ParallelDescriptor::IOProcessor(); }
int backend_rank() { return g_initialized ? amrex::ParallelDescriptor::MyProc() : 0; }
int backend_nranks() { return g_initialized ? amrex::ParallelDescriptor::NProcs() : 1; }
void backend_barrier() {
  if (g_initialized) amrex::ParallelDescriptor::Barrier();
}
void backend_reduce_minmax(double& lo, double& hi) {
  if (!g_initialized) return;
  amrex::ParallelDescriptor::ReduceRealMin(lo);
  amrex::ParallelDescriptor::ReduceRealMax(hi);
}
void backend_bcast(double* v, int n) {
  if (!g_initialized || amrex::ParallelDescriptor::NProcs() == 1) return;
  amrex::ParallelDescriptor::Bcast(v, n, amrex::ParallelDescriptor::IOProcessorNumber());
}
void backend_bcast_bytes(void* p, std::size_t n) {
  if (!g_initialized || amrex::ParallelDescriptor::NProcs() == 1) return;
  const std::size_t chunk = std::size_t(1) << 30;
  for (std::size_t off = 0; off < n; off += chunk) {
    const int cnt = static_cast<int>(std::min(chunk, n - off));
    MPI_Bcast(static_cast<char*>(p) + off, cnt, MPI_BYTE, amrex::ParallelDescriptor::IOProcessorNumber(),
              amrex::ParallelDescriptor::Communicator());
  }
}
void amrex_state_gather_rt(Simulation& s, bool stats_only) { level_of(s).gather_fields(s, stats_only); }
void amrex_state_scatter_rt(Simulation& s) { level_of(s).scatter_fields(s); }

void amrex_rt_prep(Simulation& s) { level_of(s).prep(s); }
void amrex_rt_substep(Simulation& s) { level_of(s).substep(s); }
void amrex_rt_post(Simulation& s) { level_of(s).post(s); }
void amrex_rt_sync_fields(Simulation& s) { level_of(s).sync_fields(s); }
std::int64_t amrex_rt_nrays(const Simulation& s) { return s.amrex ? s.amrex->level->nrays() : 0; }
int amrex_rt_bytes_per_ray(const Simulation& s) { return s.amrex ? s.amrex->level->bytes_per_ray() : 0; }
void amrex_rt_gather_rays(const Simulation& s, std::vector<rrk::RayP>& out) {
  out.clear();
  if (s.amrex) s.amrex->level->gather_rays(out);
}
void amrex_rt_scatter_rays(Simulation& s, const std::vector<rrk::RayP>& in) { level_of(s).scatter_rays(in); }
std::string amrex_rt_timers(const Simulation& s) { return s.amrex ? s.amrex->level->timers() : ""; }

}  // namespace radray
