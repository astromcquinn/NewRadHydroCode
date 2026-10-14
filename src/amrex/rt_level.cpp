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
#include <functional>
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

#include <csignal>
#include <execinfo.h>
#include <fcntl.h>
#include <unistd.h>
#endif

#include "exact_sum.hpp"
#include "gas_state.hpp"
#include "rt_level_base.hpp"
#include "radray/amrex_state.hpp"
#include "radray/hydro_backend.hpp"
#include "radray/kernels/ray_kernels.hpp"
#include "radray/kernels/rt_cell.hpp"
#include "radray/ray_tracing.hpp"
#include "radray/rt.hpp"
#include "radray/rt_backend.hpp"
#include "radray/rt_moments.hpp"
#include "radray/simulation.hpp"
#include "ray_encoding.hpp"

namespace radray {

namespace {

using namespace rtl;
static_assert(std::is_same_v<ax::ParticleReal, float>, "build AMReX with AMReX_PARTICLES_PRECISION=SINGLE");
static_assert(sizeof(Long) == sizeof(long long), "int64 accumulators");


// Wrapped global 0-based index g -> index of its periodic image nearest to
// the home index h (h is a valid cell of the box).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE int near_image(int g, int h, int N) {
  int d = g - h;
  if (2 * d > N) d -= N;
  else if (2 * d < -N) d += N;
  return h + d;
}


// ---------------------------------------------------------------------------
// 128-bit fixed point: exact, associative sums of the double statistics
// ---------------------------------------------------------------------------
// (exact_sum.hpp)
using xs::allreduce_i128;
using xs::from_fx128;
using xs::i128;
using xs::to_fx128;



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
template <class Enc>
struct Level : CellLevel {
  using RayPC = ax::ParticleContainerPureSoA<Enc::kNR, Enc::kNI>;
  using PTD = typename RayPC::ParticleTileType::ParticleTileDataType;

  int ng = 2, redist_local = 3;
  std::int64_t substep_count = 0;  // global RT substep counter (dithering keys)
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

  Level(const Config& cfg, bool with_e) : CellLevel(cfg, with_e) {
    // The reach of one substep is k+1 cells for rt_cells_per_substep = k (a
    // split moves a ray by up to one cell, then it is cast k cells; a ray
    // moves k cells, plus one at a cell face): the ghost layers and the
    // local Redistribute grow with k (docs/rt_multicell.md). k = 1: the
    // parameters as given.
    const int reach = cfg.rt_cells_per_substep + 1;
    ng = std::max(cfg.amrex_ghosts, reach);
    redist_local = cfg.amrex_redist_local > 0 ? std::max(cfg.amrex_redist_local, reach) : 0;
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
      write_gas_back(s);
      return;
    }
    sync_fields(s);
  }

  // ---- amrex_balance = 1: costs and the move to a new mapping ----------------

  void box_ray_counts(ax::Vector<ax::Real>& count) const override {
    if (static_cast<int>(count.size()) != static_cast<int>(ba.size())) count.resize(ba.size(), 0.0);
    std::vector<Long> n(static_cast<std::size_t>(ba.size()), 0);
    for (ax::MFIter mfi(C); mfi.isValid(); ++mfi)
      n[static_cast<std::size_t>(mfi.index())] =
          pc->DefineAndReturnParticleTile(0, mfi.index(), mfi.LocalTileIndex()).numParticles();
    ax::ParallelDescriptor::ReduceLongSum(n.data(), static_cast<int>(n.size()));
    for (std::size_t b = 0; b < n.size(); ++b) count[static_cast<int>(b)] += static_cast<ax::Real>(n[b]);
  }

  // The cells by ParallelCopy, the rays by a Redistribute under the new
  // mapping (their positions are unchanged); the work arrays and the
  // per-box scratch rebuilt
  void remap(const ax::DistributionMapping& dm2) override {
    const double t0 = now();
    // the cells through float FabArrays (ParallelCopy needs an arithmetic
    // FAB type): pack, move, unpack
    FFab Cf(ba, dm, NC, 0);
    for (ax::MFIter mfi(C); mfi.isValid(); ++mfi) {
      auto c = C.const_array(mfi);
      auto f = Cf.array(mfi);
      ax::ParallelForOMP(mfi.validbox(), NC, [=] AMREX_GPU_DEVICE(int i, int j, int k, int m) {
        f(i, j, k, m) = cf(c, i, j, k, m);
      });
    }
    FFab Cf2(ba, dm2, NC, 0);
    Cf2.ParallelCopy(Cf, 0, 0, NC);
    CFab C2(ba, dm2, 1, 0);
    for (ax::MFIter mfi(C2); mfi.isValid(); ++mfi) {
      auto c = C2.array(mfi);
      auto f = Cf2.const_array(mfi);
      ax::ParallelForOMP(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        float* dp = reinterpret_cast<float*>(&c(i, j, k));
        for (int m = 0; m < NC; ++m) dp[m] = f(i, j, k, m);
      });
    }
    C = std::move(C2);
    dm = dm2;
    TF.define(ba, dm, 10, ng);
    A.define(ba, dm, rrk::kNdep, ng);
    pc->SetParticleDistributionMap(0, dm);
    pc->Redistribute();
    for (ax::MFIter mfi(C); mfi.isValid(); ++mfi)
      pc->DefineAndReturnParticleTile(0, mfi.index(), mfi.LocalTileIndex());
    scratch.clear();
    scratch.resize(C.local_size());
    T.sync += now() - t0;
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
    out.clear();
    gather_rays_chunks([&](const std::vector<RayP>& chunk) { out.insert(out.end(), chunk.begin(), chunk.end()); });
  }

  // The streaming gather, in a canonical order: box by box in the global
  // box order, each box's rays sorted by home cell and then by their bytes
  // (a total order; identical rays are interchangeable), so the checkpoint
  // does not depend on the mapping or the rank count (the tiles' order
  // after a step depends on the arrival order of the migration). The
  // owner of a box sends it to rank 0 in chunks of kGatherChunk; rank 0
  // hands every chunk to sink.
  void gather_rays_chunks(const std::function<void(const std::vector<RayP>&)>& sink) const override {
    typename Enc::Ctx ectx = Enc::make_ctx(N, substep_count);
    const int nranks = ax::ParallelDescriptor::NProcs();
    const int myproc = ax::ParallelDescriptor::MyProc();
    const int root = ax::ParallelDescriptor::IOProcessorNumber();
    const bool io = ax::ParallelDescriptor::IOProcessor();
    const int NN = N;
    auto key_less = [NN](const RayP& a, const RayP& b) {
      const std::int64_t ca = rrk::cell_of(a.x[0], NN) +
                              static_cast<std::int64_t>(NN) * (rrk::cell_of(a.x[1], NN) +
                                                               static_cast<std::int64_t>(NN) * rrk::cell_of(a.x[2], NN));
      const std::int64_t cb = rrk::cell_of(b.x[0], NN) +
                              static_cast<std::int64_t>(NN) * (rrk::cell_of(b.x[1], NN) +
                                                               static_cast<std::int64_t>(NN) * rrk::cell_of(b.x[2], NN));
      if (ca != cb) return ca < cb;
      return std::memcmp(&a, &b, sizeof(RayP)) < 0;
    };
    std::vector<RayP> box_rays, chunk;
#ifdef AMREX_USE_MPI
    const MPI_Comm comm = ax::ParallelDescriptor::Communicator();
#endif
    for (int b = 0; b < static_cast<int>(ba.size()); ++b) {
      const int owner = dm[b];
      if (myproc == owner) {
        auto& pt = pc->DefineAndReturnParticleTile(0, b, 0);
        const int np = static_cast<int>(pt.numParticles());
        const PTD p = pt.getParticleTileData();
        box_rays.resize(static_cast<std::size_t>(np));
        for (int n = 0; n < np; ++n) box_rays[static_cast<std::size_t>(n)] = Enc::load(p, n, ectx);
        std::sort(box_rays.begin(), box_rays.end(), key_less);
      }
      if (owner == root) {
        if (io)
          for (std::size_t off = 0; off < box_rays.size(); off += kGatherChunk) {
            const std::size_t n = std::min(kGatherChunk, box_rays.size() - off);
            chunk.assign(box_rays.begin() + off, box_rays.begin() + off + n);
            sink(chunk);
          }
        continue;
      }
      if (nranks == 1) continue;
#ifdef AMREX_USE_MPI
      if (myproc == owner) {
        long long count = static_cast<long long>(box_rays.size());
        MPI_Send(&count, 1, MPI_LONG_LONG, root, 12, comm);
        for (std::size_t off = 0; off < box_rays.size(); off += kGatherChunk) {
          const int n = static_cast<int>(std::min(kGatherChunk, box_rays.size() - off));
          MPI_Send(box_rays.data() + off, n * static_cast<int>(sizeof(RayP)), MPI_BYTE, root, 11, comm);
        }
      } else if (io) {
        long long count = 0;
        MPI_Recv(&count, 1, MPI_LONG_LONG, owner, 12, comm, MPI_STATUS_IGNORE);
        const std::size_t total = static_cast<std::size_t>(count);
        for (std::size_t off = 0; off < total; off += kGatherChunk) {
          const int n = static_cast<int>(std::min(kGatherChunk, total - off));
          chunk.resize(static_cast<std::size_t>(n));
          MPI_Recv(chunk.data(), n * static_cast<int>(sizeof(RayP)), MPI_BYTE, owner, 11, comm, MPI_STATUS_IGNORE);
          sink(chunk);
        }
      }
#endif
    }
  }

  // Rank 0 holds the whole list (a checkpoint): the rays go to their owners
  // in chunks of kScatterChunk, each written into rank 0's first tile and
  // redistributed, so rank 0's peak holds one chunk of particles and its
  // migration buffers, not the whole pool (1.1 G rays = 110 GB of
  // particles plus their send buffers at 512^3). Every rank calls
  // Redistribute the same number of times (collective).
  static constexpr std::size_t kScatterChunk = std::size_t(1) << 25;
  void scatter_rays(const std::vector<RayP>& in, bool first) override {
    typename Enc::Ctx ectx = Enc::make_ctx(N, substep_count);
    if (first)
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

// SIGTERM (scancel, the time limit, srun --kill-on-bad-exit): every rank
// writes its call stack to Backtrace.<rank> in the working directory,
// unbuffered, and dies by the signal. Where a hung run was
// (docs/stage_e_state.md section 5; tools/hyak/backtrace_lines.py resolves
// the addresses). AMReX's own handler goes through stdio and MPI_Abort,
// and the first rank done killed the others before their files were
// written.
int g_sig_rank = 0;
void sigterm_backtrace(int sig) {
  void* frames[64];
  const int n = backtrace(frames, 64);
  char name[32];
  std::snprintf(name, sizeof name, "Backtrace.%d", g_sig_rank);
  const int fd = open(name, O_WRONLY | O_CREAT | O_TRUNC, 0644);
  if (fd >= 0) {
    backtrace_symbols_fd(frames, n, fd);
    close(fd);
  }
  std::signal(sig, SIG_DFL);
  raise(sig);
}

LevelBase& level_of(Simulation& s) {
  if (!s.amrex) {
    check_rt_backend_param(s.cfg);
    s.amrex = std::make_shared<AmrexRT>();
    const bool with_e = s.gas.has_e();
    if (use_moments(s.cfg))
      s.amrex->level = make_moment_level(s.cfg, with_e);  // rt_method = m1 | otvet (moment_level.cpp)
    else if (s.cfg.ray_encoding == RayEncoding::RStruct)
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
  g_sig_rank = amrex::ParallelDescriptor::MyProc();
  std::signal(SIGTERM, sigterm_backtrace);  // Backtrace.<rank> on scancel or the time limit
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
void rt_box_ray_counts(Simulation& s, amrex::Vector<amrex::Real>& count) { level_of(s).box_ray_counts(count); }
void rt_remap(Simulation& s, const amrex::DistributionMapping& dm) { level_of(s).remap(dm); }
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
void amrex_rt_scatter_rays(Simulation& s, const std::vector<rrk::RayP>& in) { level_of(s).scatter_rays(in, true); }
void amrex_rt_scatter_rays(Simulation& s, const std::vector<rrk::RayP>& in, bool first) {
  level_of(s).scatter_rays(in, first);
}
void amrex_rt_gather_rays_chunks(const Simulation& s,
                                 const std::function<void(const std::vector<rrk::RayP>&)>& sink) {
  if (s.amrex) s.amrex->level->gather_rays_chunks(sink);
}
std::string amrex_rt_timers(const Simulation& s) { return s.amrex ? s.amrex->level->timers() : ""; }

}  // namespace radray
