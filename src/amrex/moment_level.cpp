// The moment RT on the AMReX backend (rt_method = m1 | otvet with
// rt_backend = amrex; docs/stage_e_moments.md): the serial solver of
// src/rt_moments.cpp on the z-slabs of amrex_slab = 1, MPI + OpenMP,
// bitwise identical to it for any rank count, slab thickness and thread
// count.
//
//   * cells     the RT level's CellData FabArray (rt_level_base.hpp), with
//               the gas, the sources, the radiation fields of rt_grid;
//   * moments   M: 5 groups x (N, F_x, F_y, F_z, E), one ghost layer
//               (FillBoundary of one group before every GLF subcycle); the
//               per-cell arithmetic is moment_cell.hpp, shared with the
//               serial solver;
//   * sums      every sum in the serial order: the per-plane partials of
//               step 3 and of the redshift are local (whole planes on a
//               slab) and summed in plane order; the source photons are
//               gathered on rank 0 in cell order; the photons in flight and
//               the photochemistry statistics (amrex_moment_stats = exact)
//               run as a chain over the ranks in z order (slab s lives on
//               rank s * nranks / nslab); amrex_moment_stats = fixed sums
//               the statistics in 128-bit fixed point instead (exact and
//               decomposition-independent, but not the serial rounding).
//
//   * tensor    OTVET (thin, split, screen): the direct sum per local cell
//               over the global source list (every rank holds it), or the
//               FFT convolution of the serial code one component per task
//               (OtvetConvolver: the same plans and arithmetic) on rank
//               m * nranks / 6, copied onto the slabs (ParallelCopy) and
//               normalised there; h and ndir with one ghost layer.
//
//               The cone: the blend's thin tensor and the near field as 12
//               FFT tasks; the mip-maps built per slab down to the finest
//               level a march samples (lmin), that level gathered on every
//               rank and the coarser ones rebuilt there; the march of each
//               slab's receivers; the budget sums by plane and, over the
//               receivers, as a chain; the prolongation with one ghost
//               receiver layer; the blend's Gaussian along x and y on the
//               slabs and along z on y-slabs (a transpose).
//
//   * recomb    recomb_rad = 1: every cell in launch mode (case A) counts
//               its ground-state recombinations (recomb_count_cell) into a
//               buffer on the slabs and emits them at the next substep, the
//               emitted total by plane.
//   * checkpoint the serial format (checkpoint.cpp, version 7): the moment and
//               recombination blocks streamed through rank 0 slab by slab.
#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <fftw3.h>
#include <omp.h>

#include <AMReX_BoxList.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_MultiFab.H>

#include "exact_sum.hpp"
#include "radray/constants.hpp"
#include "radray/cosmology.hpp"
#include "radray/kernels/moment_cell.hpp"
#include "radray/otvet_closure.hpp"
#include "radray/recomb_rad.hpp"
#include "radray/ray_tracing.hpp"
#include "radray/rt.hpp"
#include "radray/rt_moments.hpp"
#include "radray/rt_rates.hpp"
#include "rt_level_base.hpp"

namespace radray {
namespace rtl {

namespace {

using xs::allreduce_i128;
using xs::from_fx128;
using xs::i128;
using xs::to_fx128;

// A receiver's 6 components in a FAB (otvet_cone_prolong reads q[m])
struct PrRef {
  ax::Array4<double const> a;
  int i, j, k;
  double operator[](int m) const { return a(i, j, k, m); }
};

constexpr int kMC = 5 * kNfreq;  // the moment components: 5 g + q (q: 0 N, 1-3 F, 4 E)
constexpr int kChainTag = 7101;

MPI_Comm comm() { return ax::ParallelDescriptor::Communicator(); }

// The distance between two components of a FAB (its points, ghosts included)
template <class A>
std::int64_t comp_stride(const A& a) {
  const auto lo = ax::lbound(a);
  return static_cast<std::int64_t>(&a(lo.x, lo.y, lo.z, 1) - &a(lo.x, lo.y, lo.z, 0));
}

// A sequential sum over the ranks in z order: rank r receives the running
// values v[0..n) of rank r-1, adds its own part (local) and sends them on;
// the last rank's values are then broadcast. Every rank calls it with the
// starting values.
template <class F>
void chain(double* v, int n, F&& local) {
  const int me = ax::ParallelDescriptor::MyProc(), np = ax::ParallelDescriptor::NProcs();
  if (me > 0) MPI_Recv(v, n, MPI_DOUBLE, me - 1, kChainTag, comm(), MPI_STATUS_IGNORE);
  local(v);
  if (me < np - 1) MPI_Send(v, n, MPI_DOUBLE, me + 1, kChainTag, comm());
  if (np > 1) MPI_Bcast(v, n, MPI_DOUBLE, np - 1, comm());
}

struct MomentLevel : CellLevel {
  int stats_mode = 0;  // amrex_moment_stats: 0 exact, 1 fixed
  ax::FabArray<ax::BaseFab<float>> M;   // the moments, 1 ghost
  ax::FabArray<ax::BaseFab<float>> Wn;  // the new values of one group
  std::vector<int> zorder;              // the local boxes (local indices) in z order
  // per local box and rt_stat component (kSNgam + sp: kNfreq per cell, in
  // the serial inu order; the others one per cell): the addends of
  // radiative_physics in cell order, 0 where the serial code adds nothing
  // (a sum starting at +0.0 never becomes -0.0, so adding +0.0 is exact)
  std::vector<std::vector<std::vector<double>>> adds;
  // the local source cells: local box, cell (i, j, k), Luv, sfr (ascending
  // cell index)
  struct Src {
    int li, i, j, k;
    float Luv, sfr;
  };
  std::vector<Src> src;
  // OTVET: the tensor h (6) and the split's ndir (1), one ghost layer
  ax::FabArray<ax::BaseFab<float>> H, Nd;
  bool split = false;      // otvet_split = 1 (Nd defined)
  bool split_now = false;  // the transport uses N_dir this hydro step (not with the cone)
  std::unique_ptr<OtvetConvolver> conv;  // the FFT tasks of this rank
  int fft_threads = 1;                   // of conv's plans
  bool fft_checked = false;              // the guard has run
  int fft_threads_cfg = 0;               // amrex_moment_fft_threads
  // recomb_rad = 1: the ground-state recombination photons of the last
  // substep per cell, ngam (5) and egam (5), emitted at the next substep
  // (RecombState::ngam / egam on the slabs; every cell in launch mode)
  FFab RB;
  bool recomb = false;
  double launched = 0;
  std::int64_t nsrc = 0;
  // the photon budget (MomentState's, rank 0's values on every rank)
  double injected = 0, absorbed = 0, lost = 0;
  double sig_init[3][5];

  MomentLevel(const Config& cfg, bool with_e) : CellLevel(cfg, with_e) {
    if (!cfg.amrex_slab)
      throw std::runtime_error("rt_method = m1 | otvet with rt_backend = amrex needs amrex_slab = 1 (z-slabs)");
    stats_mode = cfg.amrex_moment_stats;
    fft_threads_cfg = cfg.amrex_moment_fft_threads;
    // whole planes per box, and the ranks in z order (the chains)
    for (int b = 0; b < static_cast<int>(ba.size()); ++b) {
      const ax::Box& bx = ba[b];
      if (bx.smallEnd(0) != 0 || bx.bigEnd(0) != N - 1 || bx.smallEnd(1) != 0 || bx.bigEnd(1) != N - 1)
        throw std::runtime_error("the moment level needs z-slabs (whole i, j planes per box)");
    }
    for (int b = 0; b < static_cast<int>(ba.size()); ++b)
      for (int b2 = 0; b2 < static_cast<int>(ba.size()); ++b2)
        if (ba[b2].smallEnd(2) > ba[b].smallEnd(2) && dm[b2] < dm[b])
          throw std::runtime_error("the moment level needs the slabs mapped to the ranks in z order");
    M.define(ba, dm, kMC, 1);
    M.setVal(0.0f);
    Wn.define(ba, dm, 5, 0);
    const int nloc = static_cast<int>(C.local_size());
    zorder.resize(nloc);
    for (int l = 0; l < nloc; ++l) zorder[l] = l;
    std::sort(zorder.begin(), zorder.end(), [&](int a, int b) {
      return ba[C.IndexArray()[a]].smallEnd(2) < ba[C.IndexArray()[b]].smallEnd(2);
    });
    adds.assign(nloc, std::vector<std::vector<double>>(kSNstat));
    if (cfg.recomb_rad) {
      recomb = true;
      RB.define(ba, dm, 2 * kNfreq, 0);
      RB.setVal(0.0f);
    }
    if (cfg.rt_method == RTMethod::OTVET) {
      H.define(ba, dm, 6, 1);
      split = otvet_closure_params(cfg).split != 0;
      if (split) Nd.define(ba, dm, 1, 1);
    }
    double nu[5];
    for (int n = 0; n < 5; ++n) nu[n] = cp::freq_init[n];
    rrk::photoion_sigma(nu, sig_init);
  }

  ax::Box local_box(int li) const { return ba[C.IndexArray()[li]]; }

  // the plane partials (part[k] of the planes this rank owns, 0 elsewhere)
  // over the ranks, summed in plane order as plane_sum
  double plane_total(std::vector<double>& part) {
    if (ax::ParallelDescriptor::NProcs() > 1) ax::ParallelDescriptor::ReduceRealSum(part.data(), static_cast<int>(part.size()));
    double t = 0;
    for (double p : part) t = t + p;
    return t;
  }

  // ---- prep: the cells, then the local source list -----------------------

  void prep(Simulation& s) override {
    if (s.tdiag) throw std::runtime_error("rt_backend = amrex: tensor_dump is not implemented for the moment RT");
    CellLevel::prep(s);
    const double t0 = now();
    // the source cells of the slabs in ascending cell order, from the cells
    // (prep set Luv and sfr from the sources)
    src.clear();
    const double sfr_min = s.units.sfr_min;
    for (int li : zorder) {
      const ax::Box bx = local_box(li);
      auto c = C.atLocalIdx(li).const_array();
      const RTGrid& G = s.rt;
      const std::int64_t N2 = static_cast<std::int64_t>(N) * N;
      const auto b0 = std::lower_bound(G.src_cell.begin(), G.src_cell.end(), bx.smallEnd(2) * N2);
      const auto b1 = std::lower_bound(b0, G.src_cell.end(), (bx.bigEnd(2) + 1) * N2);
      for (auto it = b0; it != b1; ++it) {
        const std::int64_t gc = *it;
        const int i = static_cast<int>(gc % N), j = static_cast<int>((gc / N) % N), k = static_cast<int>(gc / N2);
        const RTCell& d = c(i, j, k);
        if (!(static_cast<double>(d.sfr) > sfr_min)) continue;
        src.push_back({li, i, j, k, d.Luv, d.sfr});
      }
    }
    T.prep += now() - t0;
    s.timers.prep += now() - t0;
    if (s.cfg.rt_method == RTMethod::OTVET) tensor(s);
  }

  // ---- the OTVET tensor (once per hydro step) -------------------------------

  // Photons injected into group g per unit Luv per substep, in the units of N
  static double inject_per_L(const Simulation& s, int g) {
    return cp::fResc * cp::pop2_init[g] * s.clk.dtsec_rt / s.units.vol_rt / s.units.nH_avg;
  }

  // The convolver of this rank's FFT tasks: plans with amrex_moment_fft_threads
  // threads (threaded FFTW), else the serial plans
  void make_convolver() {
    int nt = fft_threads_cfg > 0 ? fft_threads_cfg : omp_get_max_threads();
#ifdef RADRAY_FFTW_OMP
    if (nt > 1) {
      static const bool init = fftwf_init_threads() != 0;
      if (!init) throw std::runtime_error("fftwf_init_threads failed");
      fftwf_plan_with_nthreads(nt);
    }
    conv = std::make_unique<OtvetConvolver>(N);
    if (nt > 1) fftwf_plan_with_nthreads(1);
#else
    nt = 1;
    conv = std::make_unique<OtvetConvolver>(N);
#endif
    fft_threads = nt;
  }

  // The FFT convolution P[m] = (K_m * rho) / N^3 of the serial code onto the
  // slabs (comp m of P6), one component per task (task base + m of ntask on
  // rank (base + m) * nranks / ntask): rho = L / Lmax at the sources.
  // Collective.
  void fft_components(const std::vector<OtvetSource>& gs, double Lmax, OtvetKernelFill fill, const void* ctx,
                      const char* tag, double param, FFab& P6, int base, int ntask) {
    const int np = ax::ParallelDescriptor::NProcs(), me = ax::ParallelDescriptor::MyProc();
    const std::int64_t nc = static_cast<std::int64_t>(N) * N * N;
    auto owner = [&](int m) { return static_cast<int>(static_cast<long long>(base + m) * np / ntask); };
    bool mine = false;
    for (int m = 0; m < 6; ++m) mine = mine || owner(m) == me;
    if (mine) {
      if (!conv) make_convolver();
      float* r = conv->r();
#pragma omp parallel for schedule(static)
      for (std::int64_t c = 0; c < nc; ++c) r[c] = 0.0f;
      for (const OtvetSource& q : gs) {
        double rho = 0.0;
        rho += q.L / Lmax;  // (the serial rho[c] += L / Lmax on a zero grid)
        r[q.i + static_cast<std::int64_t>(N) * (q.j + static_cast<std::int64_t>(N) * q.k)] = static_cast<float>(rho);
      }
      conv->forward();
    }
    for (int m = 0; m < 6; ++m) {
      const ax::DistributionMapping dm1(ax::Vector<int>{owner(m)});
      FFab one(ax::BoxArray(domain), dm1, 1, 0);
      if (owner(m) == me) {
        const void* ks = conv->kernel(m, fill, ctx, tag, param, conv->out());
        if (fft_threads > 1 && !fft_checked) {
          // the guard: threaded plans bitwise the serial ones (docs/stage_e_moments.md 4.1)
          if (!conv->convolve_checked(ks, conv->out()))
            throw std::runtime_error("rt_backend = amrex: the threaded OTVET FFT plans differ from the serial "
                                     "ones on this machine; set amrex_moment_fft_threads = 1");
          fft_checked = true;
          std::printf(" OTVET FFT task (rank %d): %d-thread plans bitwise the serial ones\n", me, fft_threads);
          std::fflush(stdout);
        } else {
          conv->convolve(ks, conv->out());
        }
        const float* o = conv->out();
        auto a = one.array(0);
#pragma omp parallel for schedule(static)
        for (int k = 0; k < N; ++k)
          for (int j = 0; j < N; ++j)
            for (int i = 0; i < N; ++i) a(i, j, k) = o[i + static_cast<std::int64_t>(N) * (j + static_cast<std::int64_t>(N) * k)];
      }
      P6.ParallelCopy(one, 0, m, 1);
    }
  }

  // The thin tensor h (6 comps of Hd) and N_dir = tr P / (4 pi) (NDd, may be
  // null) on the valid cells (compute_otvet_tensor / thin_tensor: the
  // direct sum, method 1, or the FFT, method 2); gs not empty
  void thin_slabs(const std::vector<OtvetSource>& gs, double Lmax, int method, FFab& Hd, FFab* NDd, int base,
                  int ntask) {
    const std::int64_t nc = static_cast<std::int64_t>(N) * N * N;
    if (method == 1) {
      for (ax::MFIter mfi(Hd); mfi.isValid(); ++mfi) {
        auto h = Hd.array(mfi);
        ax::Array4<float> nd = NDd ? NDd->array(mfi) : ax::Array4<float>{};
        const bool sp = NDd != nullptr;
        const int NN = N;
        const OtvetSource* sv = gs.data();
        const std::size_t ns = gs.size();
        ax::ParallelForOMP(mfi.validbox(), [=](int i, int j, int k) {
          auto mimg = [NN](int d) {
            d %= NN;
            if (d < 0) d += NN;
            if (2 * d > NN) d -= NN;
            return d;
          };
          double P[6] = {0, 0, 0, 0, 0, 0}, K[6];
          for (std::size_t n = 0; n < ns; ++n) {
            const OtvetSource& q2 = sv[n];
            otvet_kernel(NN, mimg(i - q2.i), mimg(j - q2.j), mimg(k - q2.k), K);
            for (int m = 0; m < 6; ++m) P[m] += q2.L * K[m];
          }
          float hc[6];
          otvet_normalise(P, hc);
          for (int m = 0; m < 6; ++m) h(i, j, k, m) = hc[m];
          if (sp) nd(i, j, k) = static_cast<float>(((P[0] + P[1]) + P[2]) / (4 * M_PI));
        });
      }
      return;
    }
    FFab P6(ba, dm, 6, 0);
    fft_components(gs, Lmax, otvet_thin_fill, nullptr, "thin", 0.0, P6, base, ntask);
    const double ncd = static_cast<double>(nc);
    for (ax::MFIter mfi(Hd); mfi.isValid(); ++mfi) {
      auto h = Hd.array(mfi);
      auto pp = P6.const_array(mfi);
      ax::Array4<float> nd = NDd ? NDd->array(mfi) : ax::Array4<float>{};
      const bool sp = NDd != nullptr;
      ax::ParallelForOMP(mfi.validbox(), [=](int i, int j, int k) {
        double Pc[6];
        for (int m = 0; m < 6; ++m) Pc[m] = static_cast<double>(pp(i, j, k, m)) / ncd;
        float hc[6];
        otvet_normalise(Pc, hc);
        for (int m = 0; m < 6; ++m) h(i, j, k, m) = hc[m];
        if (sp) nd(i, j, k) = static_cast<float>(std::max(0.0, (Pc[0] + Pc[1]) + Pc[2]) * Lmax / (4 * M_PI));
      });
    }
  }

  // The method of the thin tensor (otvet_tensor; 0: direct when Nsrc Ncell
  // is small)
  int thin_method(const Simulation& s, std::size_t ns) const {
    const double nc = static_cast<double>(N) * N * N;
    int method = s.cfg.otvet_tensor;
    if (method == 0) method = (!otvet_fft_available() || static_cast<double>(ns) * nc <= 2.7e8) ? 1 : 2;
    if (method == 2 && !otvet_fft_available()) throw std::runtime_error("otvet_tensor = fft: built without FFTW3f");
    return method;
  }

  // trmin = 1e-6 sum rho, by plane in cell order (the empty cells add 0)
  double trmin_of(const std::vector<OtvetSource>& gs, double Lmax) const {
    std::vector<double> part(static_cast<std::size_t>(N), 0.0);
    for (const OtvetSource& q2 : gs) {
      double rho = 0.0;
      rho += q2.L / Lmax;
      part[static_cast<std::size_t>(q2.k)] += rho;
    }
    double t = 0;
    for (double pk : part) t = t + pk;
    return kOtvetNoiseFloor * t;
  }

  // rt_moments otvet_update_tensor / compute_otvet_tensor /
  // otvet_closure_tensor on the slabs
  void tensor(Simulation& s) {
    const double t0 = now();
    const OtvetClosureParams p = otvet_closure_params(s.cfg);
    const Units& u = s.units;
    const RTGrid& G = s.rt;
    const std::int64_t nc = static_cast<std::int64_t>(N) * N * N;
    // the sources in cell order (source_cells: sfr > sfr_min and Luv > 0)
    std::vector<OtvetSource> gs;
    for (std::size_t n = 0; n < G.src_cell.size(); ++n) {
      if (!(static_cast<double>(G.src_sfr[n]) > u.sfr_min && G.src_Luv[n] > 0)) continue;
      const std::int64_t c = G.src_cell[n];
      gs.push_back({static_cast<int>(c % N), static_cast<int>((c / N) % N),
                    static_cast<int>(c / (static_cast<std::int64_t>(N) * N)), static_cast<double>(G.src_Luv[n])});
    }
    nsrc = static_cast<std::int64_t>(gs.size());
    const SigmaTable sig = photoion_sigma(cp::freq_init);
    std::array<double, kNfreq> q{}, sigma18{};
    for (int g = 0; g < kNfreq; ++g) {
      q[g] = inject_per_L(s, g);
      sigma18[g] = sig[0][g] * 1e18;
    }
    double qw = 0;
    for (int g = 0; g < kNfreq; ++g) qw += q[g] * sigma18[g];
    // Gamma-weighted photons of the solution per cell (nw)
    auto nw_at = [&](const ax::Array4<float const>& m, int i, int j, int k) {
      double w = 0;
      for (int g = 0; g < kNfreq; ++g) w += sigma18[g] * static_cast<double>(m(i, j, k, 5 * g));
      return w;
    };
    // the budget: the Gamma-weighted photons of the solution, by plane
    double target = 0;
    if (p.active()) {
      std::vector<double> part(static_cast<std::size_t>(N), 0.0);
      for (ax::MFIter mfi(M); mfi.isValid(); ++mfi) {
        auto m = M.const_array(mfi);
        const ax::Box bx = mfi.validbox();
#pragma omp parallel for schedule(static)
        for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
          double t = 0;
          for (int j = 0; j < N; ++j)
            for (int i = 0; i < N; ++i) t += nw_at(m, i, j, k);
          part[static_cast<std::size_t>(k)] = t;
        }
      }
      target = plane_total(part);
    }
    double Lsum = 0, Lmax = 0;
    for (const OtvetSource& q2 : gs) Lsum += q2.L;
    for (const OtvetSource& q2 : gs) Lmax = std::max(Lmax, q2.L);
    double lam = 0, lam_x = 0, ratio = 0;

    auto iso = [&]() {
      H.setVal(0.0f);
      H.setVal(static_cast<float>(1.0 / 3), 0, 3, H.nGrow());
      if (split) Nd.setVal(0.0f);
    };
    split_now = split;
    if (p.closure == 0) {
      // thin, or thin + split (compute_otvet_tensor / thin_tensor)
      if (gs.empty()) {
        iso();
      } else {
        thin_slabs(gs, Lmax, thin_method(s, gs.size()), H, split ? &Nd : nullptr, 0, 6);
      }
      if (split && target > 0 && !gs.empty()) {
        std::vector<double> part(static_cast<std::size_t>(N), 0.0);
        for (ax::MFIter mfi(Nd); mfi.isValid(); ++mfi) {
          auto nd = Nd.const_array(mfi);
          const ax::Box bx = mfi.validbox();
          for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
            double t = 0;
            for (int j = 0; j < N; ++j)
              for (int i = 0; i < N; ++i) t += static_cast<double>(nd(i, j, k));
            part[static_cast<std::size_t>(k)] = t;
          }
        }
        ratio = qw * plane_total(part) / target;
      }
    } else if (!otvet_fft_available()) {
      throw std::runtime_error("otvet_closure = screen | cone need FFTW3f");
    } else if (gs.empty() || !(Lsum > 0)) {
      iso();
    } else if (p.closure == 1) {
      // screen (otvet_closure_tensor closure 1)
      lam = p.screen_lambda;
      if (!(lam > 0)) lam = (target > 0 && qw > 0) ? target / (qw * Lsum) : 1.0;
      lam = std::max(lam, 1.0);
      const OtvetScreenCtx ctx{lam};
      FFab P6(ba, dm, 6, 0);
      fft_components(gs, Lmax, otvet_screen_fill, &ctx, p.screen_lambda > 0 ? "screen" : nullptr, lam, P6, 0, 6);
      const double trmin = trmin_of(gs, Lmax);
      const double ncd = static_cast<double>(nc);
      std::vector<double> part(static_cast<std::size_t>(N), 0.0);
      for (ax::MFIter mfi(H); mfi.isValid(); ++mfi) {
        auto h = H.array(mfi);
        auto pp = P6.const_array(mfi);
        ax::Array4<float> nd = split ? Nd.array(mfi) : ax::Array4<float>{};
        const bool sp = split;
        const ax::Box bx = mfi.validbox();
#pragma omp parallel for schedule(static)
        for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
          double t = 0;
          for (int j = 0; j < N; ++j)
            for (int i = 0; i < N; ++i) {
              double Pc[6];
              for (int m = 0; m < 6; ++m) Pc[m] = static_cast<double>(pp(i, j, k, m)) / ncd;
              float hc[6];
              otvet_normalise_safe(Pc, trmin, hc);
              for (int m = 0; m < 6; ++m) h(i, j, k, m) = hc[m];
              if (sp) nd(i, j, k) = static_cast<float>(std::max(0.0, (Pc[0] + Pc[1]) + Pc[2]) * Lmax / (4 * M_PI));
              t += (Pc[0] + Pc[1]) + Pc[2];
            }
          part[static_cast<std::size_t>(k)] = t;
        }
      }
      if (target > 0) ratio = qw * plane_total(part) * Lmax / (4 * M_PI) / target;
    } else {
      // the cone writes no N_dir (the transport is not split)
      split_now = false;
      cone_tensor(s, p, gs, Lmax, q, sigma18, qw, nw_at, lam_x, ratio);
    }
    H.FillBoundary(period);
    if (split) Nd.FillBoundary(period);
    const double secs = now() - t0;
    T.tensor += secs;
    if (p.active() && s.raylog)
      *s.raylog << " OTVET      : " << otvet_closure_summary(s.cfg) << "; lambda " << static_cast<float>(lam)
                << " lambda_x " << static_cast<float>(lam_x) << " N_dir/N " << static_cast<float>(ratio) << " sources "
                << nsrc << " " << static_cast<float>(secs) << " s\n";
  }

  // ---- the cone (otvet_closure_tensor closure 2) on the slabs ---------------

  template <class NwAt>
  void cone_tensor(const Simulation& s, const OtvetClosureParams& p, const std::vector<OtvetSource>& gs, double Lmax,
                   const std::array<double, kNfreq>& q, const std::array<double, kNfreq>& sigma18, double qw,
                   NwAt&& nw_at, double& lam_x, double& ratio) {
    const std::int64_t nc = static_cast<std::int64_t>(N) * N * N;
    const double ncd = static_cast<double>(nc);
    const OtvetCone cone(p, q, sigma18, N);
    const int co = cone.co, ncr = cone.ncr, nb = cone.nb;
    const int lstep = 1 << cone.lmin;
    for (int b = 0; b < static_cast<int>(ba.size()); ++b) {
      const ax::Box bx = ba[b];
      if (bx.smallEnd(2) % co != 0 || bx.length(2) % co != 0 || bx.smallEnd(2) % lstep != 0 ||
          bx.length(2) % lstep != 0)
        throw std::runtime_error("otvet_closure = cone: amrex_slab_planes must be a multiple of otvet_cone_coarse (" +
                                 std::to_string(co) + ") and of 2^" + std::to_string(cone.lmin));
    }
    // the Gamma-weighted photons and the ionised cells (x_HII = 1 - fHI > 0.5)
    auto ion_at = [](const ax::Array4<RTCell const>& c, int i, int j, int k) {
      return 1.0f - c(i, j, k).fHI > 0.5f;
    };

    // h_s: the split thin tensor (blend), from the Gamma-weighted N
    FFab HS;
    if (p.cone_blend > 0) {
      HS.define(ba, dm, 6, 0);
      FFab ND(ba, dm, 1, 0);
      thin_slabs(gs, Lmax, thin_method(s, gs.size()), HS, &ND, 0, 12);
      for (ax::MFIter mfi(HS); mfi.isValid(); ++mfi) {
        auto h = HS.array(mfi);
        auto nd = ND.const_array(mfi);
        auto m = M.const_array(mfi);
        ax::ParallelForOMP(mfi.validbox(), [=, &nw_at](int i, int j, int k) {
          float hc[6];
          for (int c = 0; c < 6; ++c) hc[c] = h(i, j, k, c);
          otvet_blend_hs(hc, static_cast<double>(nd(i, j, k)), nw_at(m, i, j, k), qw);
          for (int c = 0; c < 6; ++c) h(i, j, k, c) = hc[c];
        });
      }
    }

    // the near field r < 2 r_c (FFT), in L/Lmax units
    FFab PN(ba, dm, 6, 0);
    {
      const OtvetNearCtx ctx{p.cone_rc};
      fft_components(gs, Lmax, otvet_near_fill, &ctx, "near", p.cone_rc, PN, 6, 12);
    }

    // the mip-maps of the emissivity and of the capped HI column: each slab
    // builds its planes down to level lmin, which is gathered on every rank,
    // then the coarser levels (the same arithmetic as pyramid())
    std::vector<OtvetMip> R(cone.nlev), Cm(cone.nlev);
    {
      const double cap = sigma18[0] > 0 ? p.cone_taucap / sigma18[0] : 0.0;
      const Units u = s.units;
      // the local fine planes, local boxes in z order
      std::vector<float> fr, fc;
      for (int li : zorder) {
        const ax::Box bx = local_box(li);
        auto c = C.atLocalIdx(li).const_array();
        const std::size_t off = fr.size();
        const std::size_t npl = static_cast<std::size_t>(N) * N * bx.length(2);
        fr.resize(off + npl, 0.0f);
        fc.resize(off + npl, 0.0f);
        float* pr = fr.data() + off;
        float* pc = fc.data() + off;
        const int z0 = bx.smallEnd(2);
        const int NN = N;
        ax::ParallelForOMP(bx, [=](int i, int j, int k) {
          const RTCell& d = c(i, j, k);
          // HI column per cell (1e-18 cm^-2), as calc_rt_1's tau
          const double nH = static_cast<double>(d.rhog) * u.rhob_avg * cp::XHy / pc::mH_cgs;
          const float col = static_cast<float>(nH * static_cast<double>(d.fHI) * u.len_rt * 1e-18);
          pc[i + static_cast<std::size_t>(NN) * (j + static_cast<std::size_t>(NN) * (k - z0))] =
              static_cast<float>(std::min(static_cast<double>(col), cap));
        });
        for (const OtvetSource& q2 : gs) {
          if (q2.k < bx.smallEnd(2) || q2.k > bx.bigEnd(2)) continue;
          double rho = 0.0;
          rho += q2.L / Lmax;
          pr[q2.i + static_cast<std::size_t>(N) * (q2.j + static_cast<std::size_t>(N) * (q2.k - z0))] =
              static_cast<float>(rho);
        }
      }
      // down to lmin (the local boxes are contiguous in z)
      int n = N;
      for (int l = 0; l < cone.lmin; ++l) {
        const int npl = static_cast<int>(fr.size() / (static_cast<std::size_t>(n) * n));
        const int nn = n / 2;
        std::vector<float> cr(static_cast<std::size_t>(nn) * nn * (npl / 2)), cc(cr.size());
        otvet_mip_coarsen(fr.data(), n, 0, npl / 2, cr.data());
        otvet_mip_coarsen(fc.data(), n, 0, npl / 2, cc.data());
        fr.swap(cr);
        fc.swap(cc);
        n = nn;
      }
      // the level lmin on every rank (ranks in z order)
      const int np = ax::ParallelDescriptor::NProcs();
      for (int which = 0; which < 2; ++which) {
        std::vector<float>& loc = which == 0 ? fr : fc;
        OtvetMip& L = which == 0 ? R[cone.lmin] : Cm[cone.lmin];
        L.n = n;
        L.a.assign(static_cast<std::size_t>(n) * n * n, 0.0f);
        if (np == 1) {
          L.a = loc;
        } else {
          int cnt = static_cast<int>(loc.size());
          std::vector<int> cnts(np), displ(np);
          MPI_Allgather(&cnt, 1, MPI_INT, cnts.data(), 1, MPI_INT, comm());
          long long tot = 0;
          for (int r = 0; r < np; ++r) {
            displ[r] = static_cast<int>(tot);
            tot += cnts[r];
          }
          if (tot != static_cast<long long>(L.a.size())) throw std::runtime_error("otvet cone: mip level gather");
          MPI_Allgatherv(loc.data(), cnt, MPI_FLOAT, L.a.data(), cnts.data(), displ.data(), MPI_FLOAT, comm());
        }
      }
      for (int l = cone.lmin + 1; l < cone.nlev; ++l)
        for (std::vector<OtvetMip>* P : {&R, &Cm}) {
          const OtvetMip& b = (*P)[l - 1];
          OtvetMip& c = (*P)[l];
          c.n = b.n / 2;
          c.a.assign(static_cast<std::size_t>(c.n) * c.n * c.n, 0.0f);
          otvet_mip_coarsen(b.a.data(), b.n, 0, c.n, c.a.data());
        }
    }
    const double trmin = trmin_of(gs, Lmax);

    // the march of the local receivers (whole receivers per slab)
    std::vector<std::vector<float>> Pb(zorder.size());
    std::vector<std::array<int, 2>> rz(zorder.size());  // the receiver planes of each local box
    for (std::size_t a = 0; a < zorder.size(); ++a) {
      const ax::Box bx = local_box(zorder[a]);
      const int iz0 = bx.smallEnd(2) / co, iz1 = (bx.bigEnd(2) + 1) / co;
      rz[a] = {iz0, iz1};
      const std::int64_t nr = static_cast<std::int64_t>(ncr) * ncr * (iz1 - iz0);
      Pb[a].assign(static_cast<std::size_t>(nr) * nb * 6, 0.0f);
      float* pb = Pb[a].data();
#pragma omp parallel for schedule(dynamic, 64)
      for (std::int64_t r = 0; r < nr; ++r) {
        const int ix = static_cast<int>(r % ncr), iy = static_cast<int>((r / ncr) % ncr),
                  iz = iz0 + static_cast<int>(r / (static_cast<std::int64_t>(ncr) * ncr));
        cone.march(R.data(), Cm.data(), ix, iy, iz, pb + static_cast<std::size_t>(r) * nb * 6);
      }
    }
    std::vector<OtvetMip>().swap(R);
    std::vector<OtvetMip>().swap(Cm);

    // the budget (target_ion, near_tr by plane; tb in receiver order)
    std::vector<double> fac(nb, 1.0);
    double target_ion = 0;
    {
      std::vector<double> part(static_cast<std::size_t>(N), 0.0);
      for (ax::MFIter mfi(M); mfi.isValid(); ++mfi) {
        auto m = M.const_array(mfi);
        auto c = C.const_array(mfi);
        const ax::Box bx = mfi.validbox();
#pragma omp parallel for schedule(static)
        for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
          double t = 0;
          for (int j = 0; j < N; ++j)
            for (int i = 0; i < N; ++i)
              if (ion_at(c, i, j, k)) t += nw_at(m, i, j, k);
          part[static_cast<std::size_t>(k)] = t;
        }
      }
      target_ion = plane_total(part);
    }
    if (target_ion > 0 && qw > 0) {
      const double unit = qw * Lmax / (4 * M_PI);  // photons per (L/Lmax) unit of tr P
      std::vector<double> part(static_cast<std::size_t>(N), 0.0);
      for (ax::MFIter mfi(PN); mfi.isValid(); ++mfi) {
        auto pn = PN.const_array(mfi);
        auto c = C.const_array(mfi);
        const ax::Box bx = mfi.validbox();
#pragma omp parallel for schedule(static)
        for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
          double t = 0;
          for (int j = 0; j < N; ++j)
            for (int i = 0; i < N; ++i)
              if (ion_at(c, i, j, k))
                t += (static_cast<double>(pn(i, j, k, 0)) / ncd + static_cast<double>(pn(i, j, k, 1)) / ncd) +
                     static_cast<double>(pn(i, j, k, 2)) / ncd;
          part[static_cast<std::size_t>(k)] = t;
        }
      }
      const double near_tr = plane_total(part);
      // tb[b] += wr[r] tr(Pb[r][b]) over the receivers in order (a chain)
      std::vector<double> tb(nb, 0.0);
      const double cvol = static_cast<double>(co) * co * co;
      chain(tb.data(), nb, [&](double* t) {
        for (std::size_t a = 0; a < zorder.size(); ++a) {
          auto c = C.atLocalIdx(zorder[a]).const_array();
          const std::int64_t nr = static_cast<std::int64_t>(ncr) * ncr * (rz[a][1] - rz[a][0]);
          for (std::int64_t r = 0; r < nr; ++r) {
            const int ix = static_cast<int>(r % ncr), iy = static_cast<int>((r / ncr) % ncr),
                      iz = rz[a][0] + static_cast<int>(r / (static_cast<std::int64_t>(ncr) * ncr));
            int n = 0;
            for (int dz = 0; dz < co; ++dz)
              for (int dy = 0; dy < co; ++dy)
                for (int dx = 0; dx < co; ++dx) n += ion_at(c, ix * co + dx, iy * co + dy, iz * co + dz);
            const double wr = static_cast<double>(n) / cvol;
            for (int b = 0; b < nb; ++b) {
              const float* o = &Pb[a][(static_cast<std::size_t>(r) * nb + b) * 6];
              t[b] += wr * ((static_cast<double>(o[0]) + o[1]) + o[2]);
            }
          }
        }
      });
      lam_x = cone.budget(p.cone_budget, tb, near_tr, target_ion, unit, fac, ratio);
    }

    // the receivers' P on the coarse slabs, one ghost receiver layer
    ax::BoxArray bac = ba;
    bac.coarsen(co);
    ax::MultiFab Pr(bac, dm, 6, 1);
    Pr.setVal(0.0);
    for (std::size_t a = 0; a < zorder.size(); ++a) {
      auto pr = Pr.atLocalIdx(zorder[a]).array();
      const std::int64_t nr = static_cast<std::int64_t>(ncr) * ncr * (rz[a][1] - rz[a][0]);
      const float* pb = Pb[a].data();
      const int iz0 = rz[a][0];
      const std::vector<double>& fc = fac;
#pragma omp parallel for schedule(static)
      for (std::int64_t r = 0; r < nr; ++r) {
        const int ix = static_cast<int>(r % ncr), iy = static_cast<int>((r / ncr) % ncr),
                  iz = iz0 + static_cast<int>(r / (static_cast<std::int64_t>(ncr) * ncr));
        double v[6] = {0, 0, 0, 0, 0, 0};
        for (int b = 0; b < nb; ++b) {
          const float* o = pb + (static_cast<std::size_t>(r) * nb + b) * 6;
          for (int m = 0; m < 6; ++m) v[m] += fc[b] * static_cast<double>(o[m]);
        }
        for (int m = 0; m < 6; ++m) pr(ix, iy, iz, m) = v[m];
      }
    }
    std::vector<std::vector<float>>().swap(Pb);
    Pr.FillBoundary(ax::Geometry(ax::Box(ax::IntVect(0), ax::IntVect(ncr - 1)),
                                 ax::RealBox({0., 0., 0.}, {double(ncr), double(ncr), double(ncr)}), ax::CoordSys::cartesian,
                                 ax::Array<int, 3>{1, 1, 1})
                        .periodicity());

    // h_cone = normalise(P_near + prolonged P_r) on the fine cells, then the
    // blend: D = h_cone - h_s, blurred (x, y local; z on y-slabs)
    FFab& D = H;  // h_cone, then D, then h
    for (ax::MFIter mfi(H); mfi.isValid(); ++mfi) {
      auto h = H.array(mfi);
      auto pn = PN.const_array(mfi);
      auto pr = Pr.const_array(mfi);
      const ax::Box cb = ax::grow(bac[mfi.index()], 1);
      const int czlo = cb.smallEnd(2), czhi = cb.bigEnd(2);
      const int NN = N, NR = ncr, CO = co;
      ax::ParallelForOMP(mfi.validbox(), [=](int i, int j, int k) {
        double P[6];
        for (int m = 0; m < 6; ++m) P[m] = static_cast<double>(pn(i, j, k, m)) / ncd;
        otvet_cone_prolong(CO, NR, i, j, k,
                           [&](int a, int b2, int c2) {
                             // the periodic index c2 as the copy within the grown coarse box
                             int cz = c2;
                             if (cz > czhi) cz -= NR;
                             if (cz < czlo) cz += NR;
                             return PrRef{pr, a, b2, cz};
                           },
                           P);
        float hc[6];
        otvet_normalise_safe(P, trmin, hc);
        for (int m = 0; m < 6; ++m) h(i, j, k, m) = hc[m];
        (void)NN;
      });
    }
    if (!(p.cone_blend > 0)) return;
    for (ax::MFIter mfi(D); mfi.isValid(); ++mfi) {
      auto d = D.array(mfi);
      auto hs = HS.const_array(mfi);
      ax::ParallelForOMP(mfi.validbox(), 6, [=](int i, int j, int k, int m) { d(i, j, k, m) = d(i, j, k, m) - hs(i, j, k, m); });
    }
    const std::vector<double> w = otvet_gauss_weights(p.cone_blend);
    // axes x and y on the slabs (whole planes)
    for (int axis = 0; axis < 2; ++axis)
      for (ax::MFIter mfi(D); mfi.isValid(); ++mfi) {
        auto d = D.array(mfi);
        const ax::Box bx = mfi.validbox();
        const int NN = N;
#pragma omp parallel for schedule(static) collapse(2)
        for (int m = 0; m < 6; ++m)
          for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
            std::vector<double> line(NN), outl(NN);
            for (int a2 = 0; a2 < NN; ++a2) {
              for (int t = 0; t < NN; ++t) line[t] = axis == 0 ? d(t, a2, k, m) : d(a2, t, k, m);
              for (int t = 0; t < NN; ++t) outl[t] = otvet_gauss_point(w, t, NN, [&](int v) { return line[v]; });
              for (int t = 0; t < NN; ++t) (axis == 0 ? d(t, a2, k, m) : d(a2, t, k, m)) = static_cast<float>(outl[t]);
            }
          }
      }
    // axis z on y-slabs (the transpose of the slabs: whole z lines)
    {
      ax::BoxList bly;
      for (int b = 0; b < static_cast<int>(ba.size()); ++b) {
        const ax::Box& bz = ba[b];
        bly.push_back(ax::Box(ax::IntVect(0, bz.smallEnd(2), 0), ax::IntVect(N - 1, bz.bigEnd(2), N - 1)));
      }
      const ax::BoxArray bay(bly);
      FFab DY(bay, dm, 6, 0);
      DY.ParallelCopy(D, 0, 0, 6);
      for (ax::MFIter mfi(DY); mfi.isValid(); ++mfi) {
        auto d = DY.array(mfi);
        const ax::Box bx = mfi.validbox();
        const int NN = N;
#pragma omp parallel for schedule(static) collapse(2)
        for (int m = 0; m < 6; ++m)
          for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
            std::vector<double> line(NN), outl(NN);
            for (int i = 0; i < NN; ++i) {
              for (int t = 0; t < NN; ++t) line[t] = d(i, j, t, m);
              for (int t = 0; t < NN; ++t) outl[t] = otvet_gauss_point(w, t, NN, [&](int v) { return line[v]; });
              for (int t = 0; t < NN; ++t) d(i, j, t, m) = static_cast<float>(outl[t]);
            }
          }
      }
      D.ParallelCopy(DY, 0, 0, 6);
    }
    for (ax::MFIter mfi(H); mfi.isValid(); ++mfi) {
      auto h = H.array(mfi);
      auto hs = HS.const_array(mfi);
      ax::ParallelForOMP(mfi.validbox(), [=](int i, int j, int k) {
        float a[6], dd[6], o[6];
        for (int m = 0; m < 6; ++m) {
          a[m] = hs(i, j, k, m);
          dd[m] = h(i, j, k, m);
        }
        otvet_blend_finish(a, dd, o);
        for (int m = 0; m < 6; ++m) h(i, j, k, m) = o[m];
      });
    }
  }

  // ---- one RT substep: moment_step_1, rt_step_2, rt_step_3 ----------------

  // The sources (moment_step_1's inject): the photons into the source cells;
  // the injected total in cell order (the addends gathered on rank 0)
  void inject(Simulation& s) {
    const double t0 = now();
    std::vector<double> rec(src.size() * kNfreq);
    double dummy = 0;
    for (std::size_t n = 0; n < src.size(); ++n) {
      const Src& q = src[n];
      auto m = M.atLocalIdx(q.li).array();
      const std::int64_t gs_ = 5 * comp_stride(m);  // between groups
      mk::inject_cell(static_cast<double>(q.Luv), s.clk.dtsec_rt, s.units, &m(q.i, q.j, q.k, 0),
                      &m(q.i, q.j, q.k, 4), gs_, dummy, &rec[n * kNfreq]);
    }
    // rank 0: the addends of every rank in rank (= z = cell) order
    const int np = ax::ParallelDescriptor::NProcs(), me = ax::ParallelDescriptor::MyProc();
    double inj = 0;
    if (np == 1) {
      for (double a : rec) inj += a;
    } else {
      int cnt = static_cast<int>(rec.size());
      std::vector<int> cnts(me == 0 ? np : 0), displ(me == 0 ? np : 0);
      MPI_Gather(&cnt, 1, MPI_INT, cnts.data(), 1, MPI_INT, 0, comm());
      std::vector<double> all;
      if (me == 0) {
        long long tot = 0;
        for (int r = 0; r < np; ++r) {
          displ[r] = static_cast<int>(tot);
          tot += cnts[r];
        }
        all.resize(static_cast<std::size_t>(tot));
      }
      MPI_Gatherv(rec.data(), cnt, MPI_DOUBLE, all.data(), cnts.data(), displ.data(), MPI_DOUBLE, 0, comm());
      for (double a : all) inj += a;
      MPI_Bcast(&inj, 1, MPI_DOUBLE, 0, comm());
    }
    injected += inj / s.cfg.nrtcr;
    // recomb_rad = 1: the ground-state recombination photons photo_cell
    // counted in the last substep, emitted isotropically where they were
    // made; the photons added by plane (moment_step_1's plane_sum)
    if (recomb) {
      std::vector<double> part(static_cast<std::size_t>(N), 0.0);
      for (ax::MFIter mfi(M); mfi.isValid(); ++mfi) {
        auto m = M.array(mfi);
        auto rb = RB.array(mfi);
        const ax::Box bx = mfi.validbox();
#pragma omp parallel for schedule(static)
        for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
          double t = 0;
          for (int j = 0; j < N; ++j)
            for (int i = 0; i < N; ++i)
              for (int g = 0; g < kNfreq; ++g) {
                float& rn = rb(i, j, k, g);
                float& re = rb(i, j, k, kNfreq + g);
                if (!(rn > 0)) continue;
                float& n = m(i, j, k, 5 * g);
                float& e = m(i, j, k, 5 * g + 4);
                const float n0 = n;
                n = static_cast<float>(static_cast<double>(n0) + static_cast<double>(rn));
                e = static_cast<float>(static_cast<double>(e) + static_cast<double>(re));
                t += static_cast<double>(n) - static_cast<double>(n0);
                rn = 0;
                re = 0;
              }
          part[static_cast<std::size_t>(k)] = t;
        }
      }
      const double rec = plane_total(part);
      injected += rec / s.cfg.nrtcr;
      launched += rec / s.cfg.nrtcr;
    }
    T.inject += now() - t0;
  }

  // One GLF subcycle of group g (rt_moments transport_group): the face
  // ghosts of the group (FillBoundary, cross: the 7-point stencil), the
  // fluxes of the grown box into a per-cell contiguous scratch (as the
  // serial P), the update of the valid cells into Wn, copied back. Plain
  // OpenMP loops on raw indices; the arithmetic is moment_cell.hpp's.
  std::vector<std::vector<float>> pscratch;  // per local box: 9 fluxes per grown cell
  void transport_group(int g, mk::Closure cl, double C, double qg) {
    double t0 = now();
    M.FillBoundary(5 * g, 5, period, true);
    T.halo += now() - t0;
    if (pscratch.size() != M.local_size()) pscratch.resize(M.local_size());
    for (ax::MFIter mfi(M); mfi.isValid(); ++mfi) {
      double ta = now();
      auto ma = M.array(mfi, 5 * g);
      const auto glo = ax::lbound(ma), glen = ax::length(ma);  // the grown box
      const std::int64_t js = glen.x, ks = static_cast<std::int64_t>(glen.x) * glen.y, cs = comp_stride(ma);
      const std::int64_t npg = ks * glen.z;
      float* mb = &ma(glo.x, glo.y, glo.z, 0);
      const float* hb = nullptr;
      std::int64_t hcs = 0;
      if (cl != mk::Closure::M1) {
        auto ha = H.const_array(mfi);
        hb = &ha(glo.x, glo.y, glo.z, 0);
        hcs = comp_stride(ha);
      }
      const float* nb = nullptr;
      if (cl == mk::Closure::Split) {
        auto na = Nd.const_array(mfi);
        nb = &na(glo.x, glo.y, glo.z, 0);
      }
      std::vector<float>& Pv = pscratch[mfi.LocalIndex()];
      if (static_cast<std::int64_t>(Pv.size()) != 9 * npg) Pv.assign(static_cast<std::size_t>(9 * npg), 0.0f);
      float* P = Pv.data();
#pragma omp parallel for collapse(2) schedule(static)
      for (int k = 0; k < glen.z; ++k)
        for (int j = 0; j < glen.y; ++j) {
          const std::int64_t row = j * js + k * ks;
          for (int i = 0; i < glen.x; ++i) {
            const std::int64_t c = row + i;
            float hc[6];
            if (hb)
              for (int q = 0; q < 6; ++q) hc[q] = hb[c + q * hcs];
            const double ndc = nb ? qg * static_cast<double>(nb[c]) : 0.0;
            mk::flux_cell(cl, mb[c], mb[c + cs], mb[c + 2 * cs], mb[c + 3 * cs], mb[c + 4 * cs], hc, ndc, &P[9 * c]);
          }
        }
      double tb = now();
      T.tflux += tb - ta;
      // the update of the valid cells (offsets 1 into the grown box)
      auto wa = Wn.array(mfi);
      const auto vlen = ax::length(wa);
      const auto vlo = ax::lbound(wa);
      const std::int64_t wjs = vlen.x, wks = static_cast<std::int64_t>(vlen.x) * vlen.y, wcs = comp_stride(wa);
      float* wb = &wa(vlo.x, vlo.y, vlo.z, 0);
      const std::int64_t off[7] = {0, -1, 1, -js, js, -ks, ks};
#pragma omp parallel for collapse(2) schedule(static)
      for (int k = 0; k < vlen.z; ++k)
        for (int j = 0; j < vlen.y; ++j) {
          const std::int64_t row = (j + 1) * js + (k + 1) * ks + 1;
          const std::int64_t wrow = j * wjs + k * wks;
          for (int i = 0; i < vlen.x; ++i) {
            const std::int64_t c = row + i;
            float out[5];
            mk::glf_cell([&](int q, int p) { return mb[c + off[p] + q * cs]; },
                         [&](int p, int q) { return P[9 * (c + off[p]) + q]; }, C, out);
            const std::int64_t w = wrow + i;
            for (int q = 0; q < 5; ++q) wb[w + q * wcs] = out[q];
          }
        }
      double tc = now();
      T.tglf += tc - tb;
#pragma omp parallel for collapse(2) schedule(static)
      for (int k = 0; k < vlen.z; ++k)
        for (int j = 0; j < vlen.y; ++j) {
          const std::int64_t row = (j + 1) * js + (k + 1) * ks + 1;
          const std::int64_t wrow = j * wjs + k * wks;
          for (int q = 0; q < 5; ++q)
            for (int i = 0; i < vlen.x; ++i) mb[row + i + q * cs] = wb[wrow + i + q * wcs];
        }
      T.tcopy += now() - tc;
    }
  }

  void substep(Simulation& s) override {
    s.timers.substeps += 1;
    const Units u = s.units;
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
    double tw = now();
    auto lap = [&](double& t, double& t2) {
      const double tn = now();
      t += tn - tw;
      t2 += tn - tw;
      tw = tn;
    };

    // calc_rt_1
    double sig[3][kNfreq];
    std::memcpy(sig, sig_init, sizeof sig);
    for (ax::MFIter mfi(C); mfi.isValid(); ++mfi) {
      auto c = C.array(mfi);
      ax::ParallelForOMP(mfi.validbox(), [=, &sig](int i, int j, int k) {
        mk::calc_rt_1_cell(rrk::cell_ref(c(i, j, k)), u, sig);
      });
    }
    lap(s.timers.calc_rt1, T.rt1);

    // the sources (moment_inject = 0: before the transport), the transport
    // (moment_nsub GLF subcycles of Courant number 1/moment_nsub), the
    // sources after it (moment_inject = 1)
    if (s.cfg.moment_inject == 0) inject(s);
    {
      const double t0 = now();
      const int nsub = s.cfg.moment_nsub;
      const double Cn = 1.0 / nsub;
      const bool otvet = s.cfg.rt_method == RTMethod::OTVET;
      const mk::Closure cl = !otvet ? mk::Closure::M1 : split_now ? mk::Closure::Split : mk::Closure::Otvet;
      // otvet_split = 1: N_dir,g = q_g ndir
      std::array<double, kNfreq> qg{};
      for (int g = 0; g < kNfreq; ++g) qg[g] = cl == mk::Closure::Split ? inject_per_L(s, g) : 0.0;
      for (int g = 0; g < kNfreq; ++g) {
        bool any = false;
        for (ax::MFIter mfi(M); mfi.isValid() && !any; ++mfi) {
          auto m = M.const_array(mfi, 5 * g);
          const ax::Box bx = mfi.validbox();
          const auto lo = ax::lbound(bx), hi = ax::ubound(bx);
          for (int k = lo.z; k <= hi.z && !any; ++k)
            for (int j = lo.y; j <= hi.y && !any; ++j)
              for (int i = lo.x; i <= hi.x; ++i)
                if (m(i, j, k, 0) > 0) {
                  any = true;
                  break;
                }
        }
        ax::ParallelDescriptor::ReduceBoolOr(any);
        if (!any) continue;
        for (int it = 0; it < nsub; ++it) transport_group(g, cl, Cn, qg[g]);
      }
      T.trans += now() - t0;
    }
    if (s.cfg.moment_inject == 1) inject(s);
    {
      double dummy = 0;
      lap(s.timers.trace1, dummy);
    }

    // the deposit for radiative_physics: ngam = N, nu = E, dx = sum F
    for (ax::MFIter mfi(C); mfi.isValid(); ++mfi) {
      auto c = C.array(mfi);
      auto m = M.const_array(mfi);
      const std::int64_t ns = comp_stride(m), gs_ = 5 * ns;
      ax::ParallelForOMP(mfi.validbox(), [=](int i, int j, int k) {
        const float* p = &m(i, j, k, 0);
        mk::deposit_cell(p, p + 4 * ns, p + ns, p + 2 * ns, p + 3 * ns, gs_,
                         rrk::cell_ref(c(i, j, k)));
      });
    }
    lap(s.timers.sync, T.deposit);

    physics(s);
    tw = now();

    // rt_step_3: N, F, E *= fgam; the absorbed photons by plane
    {
      std::vector<double> part(static_cast<std::size_t>(N), 0.0);
      for (ax::MFIter mfi(C); mfi.isValid(); ++mfi) {
        auto c = C.const_array(mfi);
        auto m = M.array(mfi);
        const ax::Box bx = mfi.validbox();
#pragma omp parallel for schedule(static)
        for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
          double t = 0;
          for (int j = 0; j < N; ++j)
            for (int i = 0; i < N; ++i)
              for (int g = 0; g < kNfreq; ++g) {
                if (m(i, j, k, 5 * g) == 0) continue;
                t += mk::absorb_group(m(i, j, k, 5 * g), m(i, j, k, 5 * g + 4), m(i, j, k, 5 * g + 1),
                                      m(i, j, k, 5 * g + 2), m(i, j, k, 5 * g + 3), c(i, j, k).fgam[g]);
              }
          part[static_cast<std::size_t>(k)] = t;
        }
      }
      absorbed += plane_total(part) / s.cfg.nrtcr;
    }
    lap(s.timers.trace3, T.absorb);
    accumulate_rt_stats(s);
  }

  // rt_step_2: radiative_physics on every cell (physics_cell), and the
  // rt_stat sums in the serial order (exact) or in fixed point
  void physics(Simulation& s) {
    const double t0 = now();
    rrk::PhysCtx X;
    X.u = s.units;
    X.rates.cool = s.rates.coef_cool_.data();
    X.rates.alpha = s.rates.alpha_rec_.data();
    X.rates.gcol = s.rates.gamma_col_.data();
    X.rates.legacy_physics = s.cfg.legacy_physics;
    X.dtsec_rt = s.clk.dtsec_rt;
    X.cells = s.cfg.rt_cells_per_substep;
    X.photon_conserve = s.cfg.photon_conserve;
    X.ia = s.cfg.recomb_case_A ? 0 : 3;
    X.want_nrab = recomb;
    const RecombState* Rc = recomb ? s.recomb.get() : nullptr;
    if (recomb && !Rc) throw std::runtime_error("recomb_rad = 1: recomb_init was not called");
    const bool lines = s.cfg.recomb_rad_lines;
    const bool he = has_e;
    const Units u = s.units;
    const double dtsec_rt = s.clk.dtsec_rt;
    const double Ncell = s.cfg.ncellr;
    const double nrtcr = s.cfg.nrtcr;
    const double navg[3] = {u.nH_avg, u.nHe_avg, u.nHe_avg};
    int err = 0;
    Long nsf = 0, ngf = 0;  // nsfill, ngfill
    for (ax::MFIter mfi(C); mfi.isValid(); ++mfi) {
      const ax::Box bx = mfi.validbox();
      const int ncb = static_cast<int>(bx.numPts());
      const auto lo = ax::lbound(bx), len = ax::length(bx);
      auto c = C.array(mfi);
      ax::Array4<float> rb = recomb ? RB.array(mfi) : ax::Array4<float>{};
      const std::int64_t rstride = recomb ? comp_stride(rb) : 0;
      std::vector<std::vector<double>>& ad = adds[mfi.LocalIndex()];
      double* A[kSNstat];
      for (int a2 = 0; a2 < kSNstat; ++a2) {
        ad[a2].resize(static_cast<std::size_t>(a2 < kSNion ? kNfreq : 1) * ncb);
        A[a2] = ad[a2].data();
      }
#pragma omp parallel for schedule(dynamic, 64) reduction(max : err) reduction(+ : nsf, ngf)
      for (int q = 0; q < ncb; ++q) {
        const int i = lo.x + q % len.x, j = lo.y + (q / len.x) % len.y, k = lo.z + q / (len.x * len.y);
        RTCell& d = c(i, j, k);
        rrk::PhotoStats st;
        double ngsrc[3][5] = {};
        bool processed = false, source = false;
        double ngam_tot = 0;
        const int e = rrk::physics_cell(d, X, &st, he, &processed, &ngam_tot);
        err = std::max(err, e);
        if (processed && !e) {
          if (Rc) {
            float ots = 0;
            recomb_count_cell(*Rc, u, false, st.nrab, lines ? st.nrB : nullptr, st.T_in, u.vol_gas, &rb(i, j, k, 0),
                              &rb(i, j, k, kNfreq), rstride, &ots);
          }
          if (ngam_tot > u.ngam_min) ++ngf;
          if (static_cast<double>(d.sfr) > u.sfr_min) {
            ++nsf;
            source = true;
            rrk::source_stats(d.Luv, dtsec_rt, u, ngsrc);
          }
        }
        const bool on = processed && !e;
        for (int sp = 0; sp < 3; ++sp) {
          for (int inu = 0; inu < kNfreq; ++inu)
            A[kSNgam + sp][static_cast<std::size_t>(q) * kNfreq + inu] =
                (source && ngsrc[sp][inu] != 0) ? ngsrc[sp][inu] / nrtcr : 0.0;
          A[kSNion + sp][q] = on ? st.nion[sp] / navg[sp] / Ncell : 0.0;
          A[kSNcol + sp][q] = on ? st.ncol[sp] / navg[sp] / Ncell : 0.0;
          A[kSNrec + sp][q] = on ? st.nrec[sp] / navg[sp] / Ncell : 0.0;
          A[kSNrab + sp][q] = on ? st.nrab[sp] / navg[sp] / Ncell : 0.0;
          A[kSGamma + sp][q] = on ? st.gam[sp] / Ncell : 0.0;
          A[kSLambda + sp][q] = on ? st.lam[sp] / Ncell : 0.0;
        }
        A[kSDtau][q] = on ? st.dtau / Ncell : 0.0;
      }
    }
    ax::ParallelDescriptor::ReduceIntMax(err);
    if (err) throw std::runtime_error("radiative_physics: temperature index out of table");
    const double t1 = now();
    T.phys += t1 - t0;
    s.timers.physics += t1 - t0;

    // the sums
    Long cnt[2] = {nsf, ngf};
    ax::ParallelDescriptor::ReduceLongSum(cnt, 2);
    double v[kSNstat] = {};
    if (stats_mode == 0) {
      // exact: each component's addends in cell order, a sequential sum as
      // RPSink's (the components in parallel)
      chain(v, kSNstat, [&](double* a) {
#pragma omp parallel for schedule(dynamic, 1)
        for (int c2 = 0; c2 < kSNstat; ++c2) {
          double x = a[c2];
          for (int li : zorder)
            for (double y : adds[li][c2]) x = x + y;
          a[c2] = x;
        }
      });
    } else {
      // fixed: 128-bit fixed point (associative), as the ray level
      i128 fx[kSNstat] = {};
#pragma omp parallel for schedule(dynamic, 1)
      for (int c2 = 0; c2 < kSNstat; ++c2) {
        i128 x = 0;
        for (const auto& ad : adds)
          for (double y : ad[c2]) x += to_fx128(y, stat_bits(c2));
        fx[c2] = x;
      }
      allreduce_i128(fx, kSNstat);
      for (int a2 = 0; a2 < kSNstat; ++a2) v[a2] = from_fx128(fx[a2], stat_bits(a2));
    }
    RTStat& t = s.rt_stat;
    t.nsfill = cnt[0];
    t.ngfill = cnt[1];
    for (int sp = 0; sp < 3; ++sp) {
      t.ngam[sp] = v[kSNgam + sp];
      t.nion[sp] = v[kSNion + sp];
      t.ncol[sp] = v[kSNcol + sp];
      t.nrec[sp] = v[kSNrec + sp];
      t.nrab[sp] = v[kSNrab + sp];
      t.gamma[sp] = v[kSGamma + sp];
      t.lambda[sp] = v[kSLambda + sp];
    }
    t.dtau = v[kSDtau];
    T.sums += now() - t1;
  }

  // ---- post_rt: the redshift, calc_rt_f, the photons in flight -----------

  void post(Simulation& s) override {
    const double t0 = now();
    const double fscale = redshift_fscale(s.cfg, s.clk);
    const RedshiftTable RT = make_redshift_table();
    const double ngam_min = s.units.ngam_min;
    std::vector<double> part(static_cast<std::size_t>(N), 0.0);
    for (ax::MFIter mfi(M); mfi.isValid(); ++mfi) {
      auto m = M.array(mfi);
      const ax::Box bx = mfi.validbox();
      const std::int64_t gs_ = 5 * comp_stride(m);
      const std::int64_t ns = comp_stride(m);
#pragma omp parallel for schedule(static)
      for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
        double t = 0;
        for (int j = 0; j < N; ++j)
          for (int i = 0; i < N; ++i) {
            float* p = &m(i, j, k, 0);
            double l = 0;
            if (mk::redshift_cell(RT, fscale, ngam_min, p, p + 4 * ns, p + ns, p + 2 * ns, p + 3 * ns, gs_, l)) t += l;
          }
        part[static_cast<std::size_t>(k)] = t;
      }
    }
    lost += plane_total(part) / s.cfg.nrtcr;
    // the photons in flight, in cell order
    const double ts = now();
    double tot = 0;
    chain(&tot, 1, [&](double* a) {
      for (int li : zorder) {
        auto m = M.atLocalIdx(li).const_array();
        const ax::Box bx = local_box(li);
        const std::int64_t gs_ = 5 * comp_stride(m);
        for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k)
          for (int j = 0; j < N; ++j)
            for (int i = 0; i < N; ++i) a[0] = a[0] + mk::photons_cell(&m(i, j, k, 0), gs_);
      }
    });
    T.sums += now() - ts;
    s.ngam_total_backend = tot;
    if (s.raylog)
      *s.raylog << " Moments    : in flight " << static_cast<float>(tot / s.cfg.nrtcr) << " injected "
                << static_cast<float>(injected) << " absorbed " << static_cast<float>(absorbed) << " lost "
                << static_cast<float>(lost) << "\n";
    if (s.raylog && recomb)
      *s.raylog << " Recomb     : emitted " << static_cast<float>(launched) << " (included in injected)\n";
    // calc_rt_f
    const Units u = s.units;
    const bool he = has_e;
    double sig[3][5];
    std::memcpy(sig, sig_init, sizeof sig);
    for (ax::MFIter mfi(C); mfi.isValid(); ++mfi) {
      auto c = C.array(mfi);
      ax::ParallelForOMP(mfi.validbox(), [=, &sig](int i, int j, int k) { rrk::post_cell(c(i, j, k), u, sig, he); });
    }
    s.rays.Nray = 0;
    s.rays.save.nfree = 0;
    T.post += now() - t0;
    s.timers.post += now() - t0;
    if (state) {
      write_gas_back(s);
      return;
    }
    sync_fields(s);
  }

  // ---- checkpoint: the serial format, slab by slab (rt_backend.hpp) ---------

  // Components [scomp, scomp + ncomp) of src, kSlabPlanes planes at a time
  // on rank 0: fn(k0, nk, a), a the one-box FAB's Array4 (gas_state.hpp
  // gather_slabs with a component offset). Collective.
  template <class F>
  void gather_comps(const FFab& src, int scomp, int ncomp, F&& fn) const {
    const ax::DistributionMapping dm0(ax::Vector<int>{0});
    for (int k0 = 0; k0 < N; k0 += gs::kSlabPlanes) {
      const int nk = std::min(gs::kSlabPlanes, N - k0);
      const ax::Box slab(ax::IntVect(0, 0, k0), ax::IntVect(N - 1, N - 1, k0 + nk - 1));
      FFab s0(ax::BoxArray(slab), dm0, ncomp, 0);
      s0.ParallelCopy(src, scomp, 0, ncomp);
      fn(k0, nk, ax::ParallelDescriptor::IOProcessor() ? s0.const_array(0) : ax::Array4<float const>{});
    }
  }
  // The reverse: rank 0 fills (fill(k0, nk, a)), copied onto dst's
  // components [dcomp, dcomp + ncomp). Collective.
  template <class F>
  void scatter_comps(FFab& dst, int dcomp, int ncomp, F&& fill) {
    const ax::DistributionMapping dm0(ax::Vector<int>{0});
    for (int k0 = 0; k0 < N; k0 += gs::kSlabPlanes) {
      const int nk = std::min(gs::kSlabPlanes, N - k0);
      const ax::Box slab(ax::IntVect(0, 0, k0), ax::IntVect(N - 1, N - 1, k0 + nk - 1));
      FFab s0(ax::BoxArray(slab), dm0, ncomp, 0);
      fill(k0, nk, ax::ParallelDescriptor::IOProcessor() ? s0.array(0) : ax::Array4<float>{});
      dst.ParallelCopy(s0, 0, dcomp, ncomp);
    }
  }

  // The moment block (MomentState n, fx, fy, fz, e: [g][cell] each; the
  // budget) and the recombination block (RecombState ngam, egam: [cell][g];
  // ots, mode; launched, dropped, onspot)
  void ckpt_write(const std::function<void(const void*, std::size_t)>& sink) const override {
    const bool io = ax::ParallelDescriptor::IOProcessor();
    const std::int64_t ncell = static_cast<std::int64_t>(N) * N * N;
    const std::size_t npl = static_cast<std::size_t>(N) * N;
    std::vector<float> buf;
    const std::int64_t nv = kNfreq * ncell;
    for (int q = 0; q < 5; ++q) {
      sink(&nv, sizeof nv);
      for (int g = 0; g < kNfreq; ++g)
        gather_comps(M, 5 * g + q, 1, [&](int k0, int nk, ax::Array4<float const> const& a) {
          if (!io) return;
          buf.resize(npl * nk);
#pragma omp parallel for schedule(static)
          for (int k = k0; k < k0 + nk; ++k)
            for (int j = 0; j < N; ++j)
              for (int i = 0; i < N; ++i) buf[i + N * (j + static_cast<std::size_t>(N) * (k - k0))] = a(i, j, k);
          sink(buf.data(), buf.size() * sizeof(float));
        });
    }
    const double budget[3] = {injected, absorbed, lost};
    sink(budget, sizeof budget);
    if (!recomb) {
      const std::int64_t zero = 0;
      for (int q = 0; q < 4; ++q) sink(&zero, sizeof zero);
    } else {
      const std::int64_t nr = kNfreq * ncell;
      for (int which = 0; which < 2; ++which) {
        sink(&nr, sizeof nr);
        gather_comps(RB, which * kNfreq, kNfreq, [&](int k0, int nk, ax::Array4<float const> const& a) {
          if (!io) return;
          buf.resize(npl * nk * kNfreq);
#pragma omp parallel for schedule(static)
          for (int k = k0; k < k0 + nk; ++k)
            for (int j = 0; j < N; ++j)
              for (int i = 0; i < N; ++i)
                for (int g = 0; g < kNfreq; ++g)
                  buf[(i + N * (j + static_cast<std::size_t>(N) * (k - k0))) * kNfreq + g] = a(i, j, k, g);
          sink(buf.data(), buf.size() * sizeof(float));
        });
      }
      // ots (0: no cell is on the spot) and mode (1: every cell launches),
      // as the serial moment RT keeps them
      sink(&ncell, sizeof ncell);
      std::vector<float> z(npl, 0.0f);
      for (int k = 0; k < N; ++k) sink(z.data(), z.size() * sizeof(float));
      sink(&ncell, sizeof ncell);
      std::vector<std::uint8_t> one(npl, 1);
      for (int k = 0; k < N; ++k) sink(one.data(), one.size());
    }
    const double rbud[3] = {launched, 0, 0};
    sink(rbud, sizeof rbud);
  }

  void ckpt_read(const std::function<void(void*, std::size_t)>& rd, int version) override {
    const bool io = ax::ParallelDescriptor::IOProcessor();
    const std::int64_t ncell = static_cast<std::int64_t>(N) * N * N;
    const std::size_t npl = static_cast<std::size_t>(N) * N;
    std::vector<float> buf;
    for (int q = 0; q < 5; ++q) {
      std::int64_t n = 0;
      rd(&n, sizeof n);
      if (n == 0) continue;  // written before the first RT step
      if (n != kNfreq * ncell) throw std::runtime_error("checkpoint: moment RT state of the wrong size");
      for (int g = 0; g < kNfreq; ++g)
        scatter_comps(M, 5 * g + q, 1, [&](int k0, int nk, ax::Array4<float> const& a) {
          const std::size_t bytes = npl * nk * sizeof(float);
          if (!io) {
            rd(nullptr, bytes);
            return;
          }
          buf.resize(npl * nk);
          rd(buf.data(), bytes);
#pragma omp parallel for schedule(static)
          for (int k = k0; k < k0 + nk; ++k)
            for (int j = 0; j < N; ++j)
              for (int i = 0; i < N; ++i) a(i, j, k) = buf[i + N * (j + static_cast<std::size_t>(N) * (k - k0))];
        });
    }
    double budget[3];
    rd(budget, sizeof budget);
    injected = budget[0];
    absorbed = budget[1];
    lost = budget[2];
    if (version < 4) return;
    std::int64_t n = 0;
    for (int which = 0; which < 2; ++which) {
      rd(&n, sizeof n);
      if (n == 0) continue;
      if (!recomb) throw std::runtime_error("checkpoint and run differ in recomb_rad");
      if (n != kNfreq * ncell) throw std::runtime_error("checkpoint: recombination buffers of the wrong size");
      scatter_comps(RB, which * kNfreq, kNfreq, [&](int k0, int nk, ax::Array4<float> const& a) {
        const std::size_t bytes = npl * nk * kNfreq * sizeof(float);
        if (!io) {
          rd(nullptr, bytes);
          return;
        }
        buf.resize(npl * nk * kNfreq);
        rd(buf.data(), bytes);
#pragma omp parallel for schedule(static)
        for (int k = k0; k < k0 + nk; ++k)
          for (int j = 0; j < N; ++j)
            for (int i = 0; i < N; ++i)
              for (int g = 0; g < kNfreq; ++g)
                a(i, j, k, g) = buf[(i + N * (j + static_cast<std::size_t>(N) * (k - k0))) * kNfreq + g];
      });
    }
    // ots and mode: the moment RT keeps no cell on the spot
    rd(&n, sizeof n);
    rd(nullptr, static_cast<std::size_t>(n) * sizeof(float));
    rd(&n, sizeof n);
    rd(nullptr, static_cast<std::size_t>(n));
    double rbud[3];
    rd(rbud, sizeof rbud);
    launched = rbud[0];
  }

  // ---- no rays --------------------------------------------------------------

  void box_ray_counts(ax::Vector<ax::Real>& count) const override {
    if (static_cast<int>(count.size()) != static_cast<int>(ba.size())) count.resize(ba.size(), 0.0);
  }
  void remap(const ax::DistributionMapping&) override {
    throw std::runtime_error("amrex_balance = 1 does not support the moment RT");
  }
  std::int64_t nrays() const override { return 0; }
  int bytes_per_ray() const override { return 0; }
  void gather_rays(std::vector<RayP>& out) const override { out.clear(); }
  void gather_rays_chunks(const std::function<void(const std::vector<RayP>&)>&) const override {}
  void scatter_rays(const std::vector<RayP>& in, bool) override {
    if (!in.empty()) throw std::runtime_error("the moment level holds no rays");
  }

  std::string timers() const override {
    char buf[512];
    std::snprintf(buf, sizeof buf,
                  "AMReX moment timers [s]: prep %.3f calc_rt_1 %.3f inject %.3f transport %.3f (halo %.3f) "
                  "deposit %.3f phys %.3f step3 %.3f sums %.3f post %.3f sync %.3f tensor %.3f [flux %.3f glf %.3f copy %.3f]",
                  T.prep, T.rt1, T.inject, T.trans, T.halo, T.deposit, T.phys, T.absorb, T.sums, T.post, T.sync,
                  T.tensor, T.tflux, T.tglf, T.tcopy);
    return buf;
  }
};

}  // namespace

std::unique_ptr<LevelBase> make_moment_level(const Config& cfg, bool with_e) {
#ifdef AMREX_USE_GPU
  (void)cfg;
  (void)with_e;
  throw std::runtime_error("rt_method = m1 | otvet with rt_backend = amrex: CPU builds only");
#else
  return std::make_unique<MomentLevel>(cfg, with_e);
#endif
}

}  // namespace rtl
}  // namespace radray
