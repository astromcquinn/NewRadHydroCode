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
// The cone, recomb_rad = 1 and the checkpoint come in the next phases and
// are refused until then.
#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <AMReX_FArrayBox.H>

#include "exact_sum.hpp"
#include "radray/constants.hpp"
#include "radray/cosmology.hpp"
#include "radray/kernels/moment_cell.hpp"
#include "radray/otvet_closure.hpp"
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
  ax::FabArray<ax::BaseFab<float>> Pf;  // the transport fluxes of one group, 1 ghost
  ax::FabArray<ax::BaseFab<float>> Wn;  // the new values of one group
  std::vector<int> zorder;              // the local boxes (local indices) in z order
  std::vector<std::vector<CellPhys>> phys;  // per local box: the cell statistics
  // the local source cells: local box, cell (i, j, k), Luv, sfr (ascending
  // cell index)
  struct Src {
    int li, i, j, k;
    float Luv, sfr;
  };
  std::vector<Src> src;
  // OTVET: the tensor h (6) and the split's ndir (1), one ghost layer
  ax::FabArray<ax::BaseFab<float>> H, Nd;
  bool split = false;
  std::unique_ptr<OtvetConvolver> conv;  // the FFT tasks of this rank
  std::int64_t nsrc = 0;
  // the photon budget (MomentState's, rank 0's values on every rank)
  double injected = 0, absorbed = 0, lost = 0;
  double sig_init[3][5];

  MomentLevel(const Config& cfg, bool with_e) : CellLevel(cfg, with_e) {
    if (!cfg.amrex_slab)
      throw std::runtime_error("rt_method = m1 | otvet with rt_backend = amrex needs amrex_slab = 1 (z-slabs)");
    if (cfg.rt_method == RTMethod::OTVET && otvet_closure_params(cfg).closure == 2)
      throw std::runtime_error("rt_backend = amrex: otvet_closure = cone is not implemented yet");
    if (cfg.recomb_rad) throw std::runtime_error("rt_backend = amrex: recomb_rad = 1 with the moment RT is not implemented yet");
    stats_mode = cfg.amrex_moment_stats;
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
    Pf.define(ba, dm, 9, 1);
    Wn.define(ba, dm, 5, 0);
    const int nloc = static_cast<int>(C.local_size());
    zorder.resize(nloc);
    for (int l = 0; l < nloc; ++l) zorder[l] = l;
    std::sort(zorder.begin(), zorder.end(), [&](int a, int b) {
      return ba[C.IndexArray()[a]].smallEnd(2) < ba[C.IndexArray()[b]].smallEnd(2);
    });
    phys.resize(nloc);
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

  // The FFT convolution P[m] = (K_m * rho) / N^3 of the serial code onto the
  // slabs (comp m of P6), one component per task: rho = L / Lmax at the
  // sources. Collective.
  void fft_components(const std::vector<OtvetSource>& gs, double Lmax, OtvetKernelFill fill, const void* ctx,
                      const char* tag, double param, ax::FabArray<ax::BaseFab<float>>& P6) {
    const int np = ax::ParallelDescriptor::NProcs(), me = ax::ParallelDescriptor::MyProc();
    const std::int64_t nc = static_cast<std::int64_t>(N) * N * N;
    auto owner = [&](int m) { return static_cast<int>(static_cast<long long>(m) * np / 6); };
    bool mine = false;
    for (int m = 0; m < 6; ++m) mine = mine || owner(m) == me;
    if (mine) {
      if (!conv) conv = std::make_unique<OtvetConvolver>(N);
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
      ax::FabArray<ax::BaseFab<float>> one(ax::BoxArray(domain), dm1, 1, 0);
      if (owner(m) == me) {
        const void* ks = conv->kernel(m, fill, ctx, tag, param, conv->out());
        conv->convolve(ks, conv->out());
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

  // rt_moments otvet_update_tensor / compute_otvet_tensor /
  // otvet_closure_tensor (thin, split, screen) on the slabs
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
            for (int i = 0; i < N; ++i) {
              double w = 0;
              for (int g = 0; g < kNfreq; ++g) w += sigma18[g] * static_cast<double>(m(i, j, k, 5 * g));
              t += w;
            }
          part[static_cast<std::size_t>(k)] = t;
        }
      }
      target = plane_total(part);
    }
    double Lsum = 0, Lmax = 0;
    for (const OtvetSource& q2 : gs) Lsum += q2.L;
    for (const OtvetSource& q2 : gs) Lmax = std::max(Lmax, q2.L);
    double lam = 0, ratio = 0;

    auto iso = [&]() {
      H.setVal(0.0f);
      H.setVal(static_cast<float>(1.0 / 3), 0, 3, H.nGrow());
      if (split) Nd.setVal(0.0f);
    };
    // thin, or thin + split (compute_otvet_tensor / thin_tensor)
    if (p.closure == 0) {
      int method = s.cfg.otvet_tensor;
      if (method == 0)
        method = (!otvet_fft_available() || static_cast<double>(gs.size()) * static_cast<double>(nc) <= 2.7e8) ? 1 : 2;
      if (method == 2 && !otvet_fft_available()) throw std::runtime_error("otvet_tensor = fft: built without FFTW3f");
      if (gs.empty()) {
        iso();
      } else if (method == 1) {
        for (ax::MFIter mfi(H); mfi.isValid(); ++mfi) {
          auto h = H.array(mfi);
          ax::Array4<float> nd = split ? Nd.array(mfi) : ax::Array4<float>{};
          const bool sp = split;
          const int NN = N;
          const OtvetSource* sv = gs.data();
          const std::size_t ns = gs.size();
          ax::ParallelForOMP(mfi.validbox(), [=](int i, int j, int k) {
            double P[6] = {0, 0, 0, 0, 0, 0}, K[6];
            for (std::size_t n = 0; n < ns; ++n) {
              const OtvetSource& q2 = sv[n];
              auto mimg = [NN](int d) {
                d %= NN;
                if (d < 0) d += NN;
                if (2 * d > NN) d -= NN;
                return d;
              };
              otvet_kernel(NN, mimg(i - q2.i), mimg(j - q2.j), mimg(k - q2.k), K);
              for (int m = 0; m < 6; ++m) P[m] += q2.L * K[m];
            }
            float hc[6];
            otvet_normalise(P, hc);
            for (int m = 0; m < 6; ++m) h(i, j, k, m) = hc[m];
            if (sp) nd(i, j, k) = static_cast<float>(((P[0] + P[1]) + P[2]) / (4 * M_PI));
          });
        }
      } else {
        ax::FabArray<ax::BaseFab<float>> P6(ba, dm, 6, 0);
        fft_components(gs, Lmax, otvet_thin_fill, nullptr, "thin", 0.0, P6);
        const double ncd = static_cast<double>(nc);
        for (ax::MFIter mfi(H); mfi.isValid(); ++mfi) {
          auto h = H.array(mfi);
          auto pp = P6.const_array(mfi);
          ax::Array4<float> nd = split ? Nd.array(mfi) : ax::Array4<float>{};
          const bool sp = split;
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
    } else {
      // screen (otvet_closure_tensor closure 1)
      if (!otvet_fft_available()) throw std::runtime_error("otvet_closure = screen | cone need FFTW3f");
      if (gs.empty() || !(Lsum > 0)) {
        iso();
      } else {
        lam = p.screen_lambda;
        if (!(lam > 0)) lam = (target > 0 && qw > 0) ? target / (qw * Lsum) : 1.0;
        lam = std::max(lam, 1.0);
        const OtvetScreenCtx ctx{lam};
        ax::FabArray<ax::BaseFab<float>> P6(ba, dm, 6, 0);
        fft_components(gs, Lmax, otvet_screen_fill, &ctx, p.screen_lambda > 0 ? "screen" : nullptr, lam, P6);
        // trmin = 1e-6 sum rho, by plane in cell order (the empty cells add 0)
        double trmin = 0;
        {
          std::vector<double> part(static_cast<std::size_t>(N), 0.0);
          for (const OtvetSource& q2 : gs) {
            double rho = 0.0;
            rho += q2.L / Lmax;
            part[static_cast<std::size_t>(q2.k)] += rho;
          }
          double t = 0;
          for (double pk : part) t = t + pk;
          trmin = kOtvetNoiseFloor * t;
        }
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
      }
    }
    H.FillBoundary(period);
    if (split) Nd.FillBoundary(period);
    const double secs = now() - t0;
    T.tensor += secs;
    if (p.active() && s.raylog)
      *s.raylog << " OTVET      : " << otvet_closure_summary(s.cfg) << "; lambda " << static_cast<float>(lam)
                << " lambda_x " << 0.0f << " N_dir/N " << static_cast<float>(ratio) << " sources " << nsrc << " "
                << static_cast<float>(secs) << " s\n";
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
    T.inject += now() - t0;
  }

  // One GLF subcycle of group g (rt_moments transport_group)
  void transport_group(int g, mk::Closure cl, double C, double qg) {
    double t0 = now();
    M.FillBoundary(5 * g, 5, period);
    T.halo += now() - t0;
    for (ax::MFIter mfi(M); mfi.isValid(); ++mfi) {
      auto m = M.const_array(mfi, 5 * g);
      auto P = Pf.array(mfi);
      const ax::Box gbx = ax::grow(mfi.validbox(), 1);
      ax::Array4<float const> h = cl != mk::Closure::M1 ? H.const_array(mfi) : ax::Array4<float const>{};
      ax::Array4<float const> nd = cl == mk::Closure::Split ? Nd.const_array(mfi) : ax::Array4<float const>{};
      ax::ParallelForOMP(gbx, [=](int i, int j, int k) {
        float Pc[9], hc[6];
        if (cl != mk::Closure::M1)
          for (int q = 0; q < 6; ++q) hc[q] = h(i, j, k, q);
        const double ndc = cl == mk::Closure::Split ? qg * static_cast<double>(nd(i, j, k)) : 0.0;
        mk::flux_cell(cl, m(i, j, k, 0), m(i, j, k, 1), m(i, j, k, 2), m(i, j, k, 3), m(i, j, k, 4), hc, ndc, Pc);
        for (int q = 0; q < 9; ++q) P(i, j, k, q) = Pc[q];
      });
      auto w = Wn.array(mfi);
      auto Pc = Pf.const_array(mfi);
      ax::ParallelForOMP(mfi.validbox(), [=](int i, int j, int k) {
        const int di[7] = {0, -1, 1, 0, 0, 0, 0}, dj[7] = {0, 0, 0, -1, 1, 0, 0}, dk[7] = {0, 0, 0, 0, 0, -1, 1};
        float out[5];
        mk::glf_cell([&](int q, int p) { return m(i + di[p], j + dj[p], k + dk[p], q); },
                     [&](int p, int q) { return Pc(i + di[p], j + dj[p], k + dk[p], q); }, C, out);
        for (int q = 0; q < 5; ++q) w(i, j, k, q) = out[q];
      });
      auto mw = M.array(mfi, 5 * g);
      auto wc = Wn.const_array(mfi);
      ax::ParallelForOMP(mfi.validbox(), 5, [=](int i, int j, int k, int q) { mw(i, j, k, q) = wc(i, j, k, q); });
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
      const mk::Closure cl = !otvet ? mk::Closure::M1 : split ? mk::Closure::Split : mk::Closure::Otvet;
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
    X.want_nrab = false;
    const bool he = has_e;
    const Units u = s.units;
    const double dtsec_rt = s.clk.dtsec_rt;
    int err = 0;
    for (ax::MFIter mfi(C); mfi.isValid(); ++mfi) {
      const ax::Box bx = mfi.validbox();
      const int ncb = static_cast<int>(bx.numPts());
      const auto lo = ax::lbound(bx), len = ax::length(bx);
      auto c = C.array(mfi);
      std::vector<CellPhys>& ph = phys[mfi.LocalIndex()];
      ph.resize(ncb);
      CellPhys* pp = ph.data();
#pragma omp parallel for schedule(dynamic, 64)
      for (int q = 0; q < ncb; ++q) {
        const int i = lo.x + q % len.x, j = lo.y + (q / len.x) % len.y, k = lo.z + q / (len.x * len.y);
        RTCell& d = c(i, j, k);
        CellPhys& cp_ = pp[q];
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
      }
      for (int q = 0; q < ncb; ++q) err = std::max(err, pp[q].err);
    }
    ax::ParallelDescriptor::ReduceIntMax(err);
    if (err) throw std::runtime_error("radiative_physics: temperature index out of table");
    const double t1 = now();
    T.phys += t1 - t0;
    s.timers.physics += t1 - t0;

    // the sums
    Long cnt[2] = {0, 0};  // nsfill, ngfill
    for (const auto& ph : phys)
      for (const CellPhys& cp_ : ph) {
        if (!(cp_.flags & 1)) continue;
        if (cp_.flags & 2) cnt[1] += 1;
        if (cp_.flags & 4) cnt[0] += 1;
      }
    ax::ParallelDescriptor::ReduceLongSum(cnt, 2);
    const double Ncell = s.cfg.ncellr;
    const double nrtcr = s.cfg.nrtcr;
    const double navg[3] = {u.nH_avg, u.nHe_avg, u.nHe_avg};
    double v[kSNstat] = {};
    if (stats_mode == 0) {
      // exact: the addends of radiative_physics_plane, cell by cell in cell
      // order (each component its own sequential sum, as RPSink)
      chain(v, kSNstat, [&](double* a) {
        for (int li : zorder)
          for (const CellPhys& cp_ : phys[li]) {
            if (!(cp_.flags & 1)) continue;
            if (cp_.flags & 4)
              for (int inu = 0; inu < kNfreq; ++inu)
                for (int sp = 0; sp < 3; ++sp)
                  if (cp_.ngam_src[sp][inu] != 0) a[kSNgam + sp] = a[kSNgam + sp] + cp_.ngam_src[sp][inu] / nrtcr;
            const rrk::PhotoStats& st = cp_.st;
            for (int sp = 0; sp < 3; ++sp) {
              a[kSNion + sp] = a[kSNion + sp] + st.nion[sp] / navg[sp] / Ncell;
              a[kSNcol + sp] = a[kSNcol + sp] + st.ncol[sp] / navg[sp] / Ncell;
              a[kSNrec + sp] = a[kSNrec + sp] + st.nrec[sp] / navg[sp] / Ncell;
              a[kSNrab + sp] = a[kSNrab + sp] + st.nrab[sp] / navg[sp] / Ncell;
              a[kSGamma + sp] = a[kSGamma + sp] + st.gam[sp] / Ncell;
              a[kSLambda + sp] = a[kSLambda + sp] + st.lam[sp] / Ncell;
            }
            a[kSDtau] = a[kSDtau] + st.dtau / Ncell;
          }
      });
    } else {
      // fixed: 128-bit fixed point (associative), as the ray level
      i128 fx[kSNstat] = {};
      for (const auto& ph : phys)
        for (const CellPhys& cp_ : ph) {
          if (!(cp_.flags & 1)) continue;
          if (cp_.flags & 4)
            for (int inu = 0; inu < kNfreq; ++inu)
              for (int sp = 0; sp < 3; ++sp)
                if (cp_.ngam_src[sp][inu] != 0) fx[kSNgam + sp] += to_fx128(cp_.ngam_src[sp][inu] / nrtcr, kBitsRate);
          const rrk::PhotoStats& st = cp_.st;
          for (int sp = 0; sp < 3; ++sp) {
            fx[kSNion + sp] += to_fx128(st.nion[sp] / navg[sp] / Ncell, kBitsRate);
            fx[kSNcol + sp] += to_fx128(st.ncol[sp] / navg[sp] / Ncell, kBitsRate);
            fx[kSNrec + sp] += to_fx128(st.nrec[sp] / navg[sp] / Ncell, kBitsRate);
            fx[kSNrab + sp] += to_fx128(st.nrab[sp] / navg[sp] / Ncell, kBitsRate);
            fx[kSGamma + sp] += to_fx128(st.gam[sp] / Ncell, kBitsGamma);
            fx[kSLambda + sp] += to_fx128(st.lam[sp] / Ncell, kBitsLambda);
          }
          fx[kSDtau] += to_fx128(st.dtau / Ncell, kBitsTau);
        }
      allreduce_i128(fx, kSNstat);
      for (int a = 0; a < kSNstat; ++a) v[a] = from_fx128(fx[a], stat_bits(a));
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
                  "deposit %.3f phys %.3f step3 %.3f sums %.3f post %.3f sync %.3f tensor %.3f",
                  T.prep, T.rt1, T.inject, T.trans, T.halo, T.deposit, T.phys, T.absorb, T.sums, T.post, T.sync,
                  T.tensor);
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
