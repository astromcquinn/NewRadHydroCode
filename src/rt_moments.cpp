// Moment-based radiative transfer (rt_moments.hpp, docs/rt_moments.md).
#include "radray/rt_moments.hpp"
#include "radray/otvet_closure.hpp"
#include "radray/recomb_rad.hpp"
#include "radray/rt_tensor_diag.hpp"
#include "radray/kernels/moment_cell.hpp"
#include "radray/kernels/rt_cell.hpp"

#include <algorithm>
#include <cmath>
#include <ostream>
#include <stdexcept>

#include "radray/constants.hpp"
#include "radray/fortran.hpp"
#include "radray/parallel.hpp"
#include "radray/ray_tracing.hpp"
#include "radray/rt_rates.hpp"

namespace radray {

using namespace pc;
using namespace cp;

// ---------------------------------------------------------------------------
// Parameters
// ---------------------------------------------------------------------------

bool set_rt_method_param(Config& c, const std::string& key, const std::string& val) {
  if (set_recomb_param(c, key, val)) return true;         // recomb_rad.hpp
  if (set_otvet_closure_param(c, key, val)) return true;  // otvet_closure.hpp
  if (key == "rt_method") {
    if (val == "rays")
      c.rt_method = RTMethod::Rays;
    else if (val == "m1" || val == "M1")
      c.rt_method = RTMethod::M1;
    else if (val == "otvet" || val == "OTVET")
      c.rt_method = RTMethod::OTVET;
    else
      throw std::runtime_error("rt_method must be rays, m1 or otvet");
    return true;
  }
  if (key == "moment_nsub") {
    c.moment_nsub = std::stoi(val);
    if (c.moment_nsub < 3) throw std::runtime_error("moment_nsub must be >= 3 (3-D Courant limit)");
    return true;
  }
  if (key == "moment_inject") {
    c.moment_inject = std::stoi(val);
    if (c.moment_inject != 0 && c.moment_inject != 1) throw std::runtime_error("moment_inject must be 0 or 1");
    return true;
  }
  if (key == "otvet_tensor") {
    if (val == "auto" || val == "0")
      c.otvet_tensor = 0;
    else if (val == "direct" || val == "1")
      c.otvet_tensor = 1;
    else if (val == "fft" || val == "2")
      c.otvet_tensor = 2;
    else
      throw std::runtime_error("otvet_tensor must be auto, direct or fft");
    return true;
  }
  return false;
}

void copy_rt_method(const Config& from, Config& to) {
  to.rt_method = from.rt_method;
  to.moment_nsub = from.moment_nsub;
  to.otvet_tensor = from.otvet_tensor;
  to.moment_inject = from.moment_inject;
  to.otvet_closure = from.otvet_closure;
  to.otvet_split = from.otvet_split;
  to.otvet_screen_lambda = from.otvet_screen_lambda;
  to.otvet_cone_nside = from.otvet_cone_nside;
  to.otvet_cone_coarse = from.otvet_cone_coarse;
  to.otvet_cone_rc = from.otvet_cone_rc;
  to.otvet_cone_taucap = from.otvet_cone_taucap;
  to.otvet_cone_smax = from.otvet_cone_smax;
  to.otvet_cone_budget = from.otvet_cone_budget;
  to.otvet_cone_blend = from.otvet_cone_blend;
  copy_recomb_param(from, to);
  check_otvet_closure_param(to);
}

const char* rt_method_name(RTMethod m) {
  switch (m) {
    case RTMethod::Rays: return "rays";
    case RTMethod::M1: return "m1";
    case RTMethod::OTVET: return "otvet";
  }
  return "?";
}

std::string rt_method_summary(const Config& c) {
  std::string s = std::string("rt_method = ") + rt_method_name(c.rt_method);
  if (use_moments(c)) {
    s += " (moment_nsub = " + std::to_string(c.moment_nsub);
    if (c.rt_method == RTMethod::OTVET)
      s += std::string(", otvet_tensor = ") +
           (c.otvet_tensor == 0 ? "auto" : c.otvet_tensor == 1 ? "direct" : "fft");
    if (c.moment_inject) s += ", moment_inject = 1";
    const std::string cl = otvet_closure_summary(c);  // otvet_closure.hpp
    if (!cl.empty()) s += ", " + cl;
    s += ")";
  }
  // recomb_rad.hpp; printed only when not the original's case B
  if (c.recomb_case_A) s += ", recomb_case = A";
  if (c.recomb_rad_lines) s += ", recomb_rad_lines = 1";
  if (c.recomb_rad && use_moments(c))
    s += ", recomb_rad = 1 (every cell emits into the moments)";
  else if (c.recomb_rad)
    s += ", recomb_rad = 1 (nside " + std::to_string(c.recomb_rad_nside) + ", block " +
         std::to_string(c.recomb_rad_block) + ", pool " + std::to_string(c.recomb_rad_pool).substr(0, 4) + ")";
  return s;
}

// ---------------------------------------------------------------------------
// OTVET tensor
// ---------------------------------------------------------------------------

namespace {
// 1/r0^2 of the source cell itself (r0 = half a cell): an isotropic tensor
// that dominates the neighbours (weight 1 at one cell).
constexpr double kSelfWeight = 4.0;

inline int min_image(int d, int N) {
  d %= N;
  if (d < 0) d += N;
  if (2 * d > N) d -= N;
  return d;  // (-N/2, N/2]
}

struct SrcCell {
  int i, j, k;  // 0-based
  double L;
};

std::vector<SrcCell> source_cells(const RTGrid& G, const Units& u) {
  std::vector<SrcCell> v;
  for (int k = 0; k < G.N; ++k)
    for (int j = 0; j < G.N; ++j)
      for (int i = 0; i < G.N; ++i) {
        const std::int64_t c = G.cell(i + 1, j + 1, k + 1);
        if (static_cast<double>(G.sfr[c]) > u.sfr_min && G.Luv[c] > 0)
          v.push_back({i, j, k, static_cast<double>(G.Luv[c])});
      }
  return v;
}
}  // namespace

void otvet_kernel(int N, int dx, int dy, int dz, double K[6]) {
  if (dx == 0 && dy == 0 && dz == 0) {
    K[0] = K[1] = K[2] = kSelfWeight / 3;
    K[3] = K[4] = K[5] = 0;
    return;
  }
  const double x = dx, y = dy, z = dz;
  const double r2 = (x * x + y * y) + z * z;
  const double w = 1 / (r2 * r2);
  const bool tx = 2 * std::abs(dx) == N, ty = 2 * std::abs(dy) == N, tz = 2 * std::abs(dz) == N;
  K[0] = w * x * x;
  K[1] = w * y * y;
  K[2] = w * z * z;
  K[3] = (tx || ty) ? 0.0 : w * x * y;
  K[4] = (tx || tz) ? 0.0 : w * x * z;
  K[5] = (ty || tz) ? 0.0 : w * y * z;
}

int compute_otvet_tensor(const RTGrid& G, const Units& u, int method, std::vector<float>& h,
                         std::int64_t* nsrc_out) {
  const int N = G.N;
  const std::int64_t nc = G.ncell;
  h.assign(6 * nc, 0.0f);
  const std::vector<SrcCell> src = source_cells(G, u);
  if (nsrc_out) *nsrc_out = static_cast<std::int64_t>(src.size());
  if (method == 0)
    method = (!otvet_fft_available() ||
              static_cast<double>(src.size()) * static_cast<double>(nc) <= 2.7e8)
                 ? 1
                 : 2;
  if (method == 2 && !otvet_fft_available())
    throw std::runtime_error("otvet_tensor = fft: built without FFTW3f");

  auto store = [&](std::int64_t c, const double P[6]) {
    const double tr = (P[0] + P[1]) + P[2];
    float* hc = &h[6 * c];
    if (!(tr > 0)) {
      hc[0] = hc[1] = hc[2] = static_cast<float>(1.0 / 3);
      hc[3] = hc[4] = hc[5] = 0;
      return;
    }
    for (int m = 0; m < 6; ++m) hc[m] = static_cast<float>(P[m] / tr);
  };

  if (src.empty()) {
    for (std::int64_t c = 0; c < nc; ++c) {
      const double P[6] = {0, 0, 0, 0, 0, 0};
      store(c, P);
    }
    return method;
  }

  if (method == 1) {
#pragma omp parallel for schedule(static) if (par::worth(nc))
    for (int k = 0; k < N; ++k)
      for (int j = 0; j < N; ++j)
        for (int i = 0; i < N; ++i) {
          double P[6] = {0, 0, 0, 0, 0, 0}, K[6];
          for (const SrcCell& s : src) {
            otvet_kernel(N, min_image(i - s.i, N), min_image(j - s.j, N), min_image(k - s.k, N),
                         K);
            for (int m = 0; m < 6; ++m) P[m] += s.L * K[m];
          }
          store(G.cell(i + 1, j + 1, k + 1), P);
        }
    return 1;
  }

  // FFT: normalise the luminosities (the ratio P/trP does not depend on it)
  double Lmax = 0;
  for (const SrcCell& s : src) Lmax = std::max(Lmax, s.L);
  std::vector<double> rho(nc, 0.0);
  for (const SrcCell& s : src) rho[G.cell(s.i + 1, s.j + 1, s.k + 1)] += s.L / Lmax;
  std::array<std::vector<double>, 6> P;
  otvet_fft_convolve(N, rho, P);
#pragma omp parallel for schedule(static) if (par::worth(nc))
  for (std::int64_t c = 0; c < nc; ++c) {
    const double Pc[6] = {P[0][c], P[1][c], P[2][c], P[3][c], P[4][c], P[5][c]};
    store(c, Pc);
  }
  return 2;
}

// ---------------------------------------------------------------------------
// State
// ---------------------------------------------------------------------------

void moment_init(Simulation& s) {
  if (s.moments && s.moments->N == s.cfg.Nrt) return;
  auto m = std::make_shared<MomentState>();
  m->N = s.cfg.Nrt;
  m->ncell = s.cfg.Nrtc;
  const std::int64_t n = kNfreq * m->ncell;
  for (auto* v : {&m->n, &m->fx, &m->fy, &m->fz, &m->e}) v->assign(n, 0.0f);
  s.moments = m;
}

namespace {
// The rt_grid fields of cell c as the kernels' pointer view
rrk::CellRef grid_ref(RTGrid& G, std::int64_t c) {
  rrk::CellRef r;
  r.rhog = &G.rhog[c];
  r.fHI = &G.fHI[c];
  r.fHeI = &G.fHeI[c];
  r.fHeII = &G.fHeII[c];
  r.temp = &G.temp[c];
  r.Luv = &G.Luv[c];
  r.sfr = &G.sfr[c];
  r.dx = G.dx_p(c);
  r.gamma = G.gamma_p(c);
  r.lambda = G.lambda_p(c);
  r.nu = G.nu_p(c);
  r.ngam = G.ngam_p(c);
  r.fgam = G.fgam_p(c);
  r.tau = G.tau_p(c);
  return r;
}

// Photons injected into group g per unit Luv per substep, in the units of N
// (the injection of moment_step_1)
inline double inject_per_L(const Simulation& s, int g) {
  return fResc * pop2_init[g] * s.clk.dtsec_rt / s.units.vol_rt / s.units.nH_avg;
}

// The OTVET tensor of the hydro step: the original thin tensor, or an
// alternative closure (otvet_closure.hpp)
void otvet_update_tensor(Simulation& s) {
  MomentState& M = *s.moments;
  const OtvetClosureParams p = otvet_closure_params(s.cfg);
  if (!p.active()) {
    M.tensor_method = compute_otvet_tensor(s.rt, s.units, s.cfg.otvet_tensor, M.h, &M.nsrc);
    M.ndir.clear();
    return;
  }
  const RTGrid& G = s.rt;
  const Units& u = s.units;
  const std::int64_t nc = G.ncell;
  OtvetClosureIn in;
  in.N = G.N;
  for (int k = 0; k < G.N; ++k)
    for (int j = 0; j < G.N; ++j)
      for (int i = 0; i < G.N; ++i) {
        const std::int64_t c = G.cell(i + 1, j + 1, k + 1);
        if (static_cast<double>(G.sfr[c]) > u.sfr_min && G.Luv[c] > 0)
          in.src.push_back({i, j, k, static_cast<double>(G.Luv[c])});
      }
  M.nsrc = static_cast<std::int64_t>(in.src.size());
  const SigmaTable sig = photoion_sigma(freq_init);
  for (int g = 0; g < kNfreq; ++g) {
    in.q[g] = inject_per_L(s, g);
    in.sigma18[g] = sig[0][g] * 1e18;
  }
  // the solution's Gamma-weighted photons per cell
  std::vector<double> nw(nc, 0.0);
#pragma omp parallel for schedule(static) if (par::worth(nc))
  for (std::int64_t c = 0; c < nc; ++c) {
    double t = 0;
    for (int g = 0; g < kNfreq; ++g) t += in.sigma18[g] * static_cast<double>(M.n[g * nc + c]);
    nw[c] = t;
  }
  in.nw = nw.data();
  // HI column per cell (1e-18 cm^-2), as calc_rt_1's tau
  std::vector<float> col;
  if (p.closure == 2) {
    col.assign(nc, 0.0f);
#pragma omp parallel for schedule(static) if (par::worth(nc))
    for (std::int64_t c = 0; c < nc; ++c) {
      const double nH = static_cast<double>(G.rhog[c]) * u.rhob_avg * XHy / mH_cgs;
      col[c] = static_cast<float>(nH * static_cast<double>(G.fHI[c]) * u.len_rt * 1e-18);
    }
    in.col = col.data();
  }
  std::vector<float> xion;
  if (p.closure == 2) {
    xion.assign(nc, 0.0f);
#pragma omp parallel for schedule(static) if (par::worth(nc))
    for (std::int64_t c = 0; c < nc; ++c) xion[c] = 1.0f - G.fHI[c];
    in.xion = xion.data();
  }
  OtvetClosureOut out;
  otvet_closure_tensor(p, s.cfg.otvet_tensor, in, out);
  M.h = std::move(out.h);
  M.ndir = std::move(out.ndir);
  M.tensor_method = p.closure == 0 ? 0 : 2;
  M.closure_lambda = out.lambda;
  M.closure_lam_x = out.lam_x;
  M.closure_ratio = out.ratio;
  M.closure_seconds = out.seconds;
  if (s.raylog)
    *s.raylog << " OTVET      : " << otvet_closure_summary(s.cfg) << "; lambda "
              << static_cast<float>(out.lambda) << " lambda_x " << static_cast<float>(out.lam_x)
              << " N_dir/N " << static_cast<float>(out.ratio) << " sources " << M.nsrc << " "
              << static_cast<float>(out.seconds) << " s\n";
}
}  // namespace

void moment_prep(Simulation& s) {
  moment_init(s);
  MomentState& M = *s.moments;
  if (s.cfg.rt_method == RTMethod::OTVET) {
    otvet_update_tensor(s);
  } else {
    M.h.clear();
  }
}

double moment_photons_total(const Simulation& s) {
  if (!s.moments) return 0;
  const MomentState& M = *s.moments;
  double t = 0;
  for (std::int64_t c = 0; c < M.ncell; ++c) t = t + mk::photons_cell(&M.n[c], M.ncell);
  return t;
}

// ---------------------------------------------------------------------------
// Transport
// ---------------------------------------------------------------------------

namespace {

// One GLF subcycle of group g with Courant number C (periodic box; the
// per-cell arithmetic is mk::flux_cell and mk::glf_cell, moment_cell.hpp)
void transport_group(MomentState& M, int g, bool otvet, double C, const float* ndir, double qg) {
  const int N = M.N;
  const std::int64_t nc = M.ncell;
  float* n = &M.n[g * nc];
  float* fx = &M.fx[g * nc];
  float* fy = &M.fy[g * nc];
  float* fz = &M.fz[g * nc];
  float* e = &M.e[g * nc];
  // Per cell: P (xx, yy, zz, xy, xz, yz) and the E fluxes nu F (x, y, z)
  float* P = M.P.data();
  const mk::Closure cl = !otvet ? mk::Closure::M1 : ndir ? mk::Closure::Split : mk::Closure::Otvet;

#pragma omp parallel for schedule(static) if (par::worth(nc))
  for (std::int64_t c = 0; c < nc; ++c) {
    const double nd = cl == mk::Closure::Split ? qg * static_cast<double>(ndir[c]) : 0.0;
    mk::flux_cell(cl, n[c], fx[c], fy[c], fz[c], e[c], otvet ? &M.h[6 * c] : nullptr, nd, &P[9 * c]);
  }

  float* nn = M.nn.data();
  float* nfx = M.nfx.data();
  float* nfy = M.nfy.data();
  float* nfz = M.nfz.data();
  float* ne = M.ne.data();
  const std::int64_t sj = N, sk = static_cast<std::int64_t>(N) * N;

#pragma omp parallel for schedule(static) if (par::worth(nc))
  for (int k = 0; k < N; ++k) {
    const std::int64_t kl = (k == 0 ? N - 1 : k - 1) * sk, kr = (k == N - 1 ? 0 : k + 1) * sk;
    for (int j = 0; j < N; ++j) {
      const std::int64_t jl = (j == 0 ? N - 1 : j - 1) * sj, jr = (j == N - 1 ? 0 : j + 1) * sj;
      const std::int64_t row = j * sj + k * sk;
      for (int i = 0; i < N; ++i) {
        const std::int64_t c = row + i;
        const std::int64_t xl = row + (i == 0 ? N - 1 : i - 1), xr = row + (i == N - 1 ? 0 : i + 1);
        const std::int64_t yl = c - j * sj + jl, yr = c - j * sj + jr;
        const std::int64_t zl = c - k * sk + kl, zr = c - k * sk + kr;
        const std::int64_t st[7] = {c, xl, xr, yl, yr, zl, zr};
        const float* v[5] = {n, fx, fy, fz, e};
        float out[5];
        mk::glf_cell([&](int q, int p) { return v[q][st[p]]; }, [&](int p, int m) { return P[9 * st[p] + m]; }, C, out);
        nn[c] = out[0];
        nfx[c] = out[1];
        nfy[c] = out[2];
        nfz[c] = out[3];
        ne[c] = out[4];
      }
    }
  }
  // copy back in parallel (5 arrays of N^3: serial, this dominated at 256^3)
#pragma omp parallel for schedule(static) if (par::worth(nc))
  for (int k = 0; k < N; ++k) {
    const std::int64_t a = k * sk, b = a + sk;
    std::copy(nn + a, nn + b, n + a);
    std::copy(nfx + a, nfx + b, fx + a);
    std::copy(nfy + a, nfy + b, fy + a);
    std::copy(nfz + a, nfz + b, fz + a);
    std::copy(ne + a, ne + b, e + a);
  }
}

// rt_tensor_diag.hpp: the Gamma_HI-weighted moments N, F and D N of the
// state after the transport, D = the closure used (OTVET h or M1)
void moment_tensor_diag(const MomentState& M, TensorDiag& T, bool otvet, const float* ndir,
                        const double* qg) {
  const std::int64_t nc = M.ncell;
#pragma omp parallel for schedule(static) if (par::worth(nc))
  for (std::int64_t c = 0; c < nc; ++c) {
    double nu[kNfreq], sig[3][kNfreq];
    bool any = false;
    for (int g = 0; g < kNfreq; ++g) {
      const std::int64_t i = g * nc + c;
      nu[g] = M.n[i] > 0 ? static_cast<double>(M.e[i]) / M.n[i] : 0.0;
      any = any || M.n[i] > 0;
    }
    if (!any) continue;
    rrk::photoion_sigma(nu, sig);
    double* S = &T.acc[10 * c];
    for (int g = 0; g < kNfreq; ++g) {
      const std::int64_t i = g * nc + c;
      const double Ng = M.n[i];
      if (!(Ng > 0)) continue;
      const double w = sig[0][g] * 1e18;
      double Q[6];
      if (otvet && ndir) {
        mk::otvet_split_pressure(Ng, qg[g] * static_cast<double>(ndir[c]), &M.h[6 * c], Q);
      } else if (otvet) {
        for (int m = 0; m < 6; ++m) Q[m] = Ng * static_cast<double>(M.h[6 * c + m]);
      } else {
        mk::m1_pressure(Ng, M.fx[i], M.fy[i], M.fz[i], Q);
      }
      S[0] += w * Ng;
      S[1] += w * M.fx[i];
      S[2] += w * M.fy[i];
      S[3] += w * M.fz[i];
      for (int m = 0; m < 6; ++m) S[4 + m] += w * Q[m];
    }
  }
}

// Sum of per-plane partial sums in plane order (thread-count independent)
template <class F>
double plane_sum(int N, F&& plane) {
  std::vector<double> part(N, 0.0);
#pragma omp parallel for schedule(static) if (par::worth(static_cast<long long>(N) * N * N))
  for (int k = 0; k < N; ++k) part[k] = plane(k);
  double t = 0;
  for (int k = 0; k < N; ++k) t = t + part[k];
  return t;
}

}  // namespace

// ---------------------------------------------------------------------------
// Substep
// ---------------------------------------------------------------------------

void moment_step_1(Simulation& s) {
  double tw = wall_seconds();
  auto lap = [&](double& t) {
    const double now = wall_seconds();
    t += now - tw;
    tw = now;
  };
  s.timers.substeps += 1;
  moment_init(s);
  MomentState& M = *s.moments;
  RTGrid& G = s.rt;
  const Units& u = s.units;
  const int Nrt = s.cfg.Nrt;
  const std::int64_t nc = G.ncell;

  // Initialize RT and ray variables (rt_stat(1), ray_stat(1)) as rt_step_1
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

  // calc_rt_1 (with the B3/B4 fixes: dx reset and tau in every cell; tau is
  // a diagnostic here, the moment solver does not read it)
  const SigmaTable sig = photoion_sigma(freq_init);
  double sigd[3][kNfreq];
  for (int sp = 0; sp < 3; ++sp)
    for (int q = 0; q < kNfreq; ++q) sigd[sp][q] = sig[sp][q];
#pragma omp parallel for schedule(static) if (par::worth(nc))
  for (std::int64_t c = 0; c < nc; ++c) mk::calc_rt_1_cell(grid_ref(G, c), u, sigd);
  lap(s.timers.calc_rt1);

  // Sources and recombination photons (moment_inject = 0: before the
  // transport, the original order; 1: after it, so that photo_cell absorbs
  // them in the cell that emitted them before they move on)
  auto inject = [&]() {
    // Sources: the photons ray_source emits, into the source cell, isotropic
    double inj = 0;
    for (std::int64_t c = 0; c < nc; ++c) {
      if (!(static_cast<double>(G.sfr[c]) > u.sfr_min)) continue;
      mk::inject_cell(static_cast<double>(G.Luv[c]), s.clk.dtsec_rt, u, &M.n[c], &M.e[c], nc, inj);
    }
    M.injected += inj / s.cfg.nrtcr;

    // recomb_rad = 1 (recomb_rad.hpp): the ground-state recombination photons
    // photo_cell counted in the last substep, emitted isotropically where they
    // were made (every cell is in launch mode, i.e. case A; no blocks, no pool)
    if (s.recomb && s.recomb->active()) {
      RecombState& Rc = *s.recomb;
      const double rec = plane_sum(Nrt, [&](int k) {
        double t = 0;
        const std::int64_t c0 = static_cast<std::int64_t>(k) * Nrt * Nrt;
        for (std::int64_t c = c0; c < c0 + static_cast<std::int64_t>(Nrt) * Nrt; ++c)
          for (int g = 0; g < kNfreq; ++g) {
            const std::size_t kb = static_cast<std::size_t>(c) * kNfreq + g;
            if (!(Rc.ngam[kb] > 0)) continue;
            const std::int64_t i = g * nc + c;
            const float n0 = M.n[i];
            M.n[i] = static_cast<float>(static_cast<double>(n0) + static_cast<double>(Rc.ngam[kb]));
            M.e[i] = static_cast<float>(static_cast<double>(M.e[i]) + static_cast<double>(Rc.egam[kb]));
            t += static_cast<double>(M.n[i]) - static_cast<double>(n0);
            Rc.ngam[kb] = 0;
            Rc.egam[kb] = 0;
          }
        return t;
      });
      M.injected += rec / s.cfg.nrtcr;
      Rc.launched += rec / s.cfg.nrtcr;
    }
  };
  if (s.cfg.moment_inject == 0) inject();

  // Transport: moment_nsub GLF subcycles of Courant number 1/moment_nsub
  // (one cell of light travel per RT substep)
  if (M.P.size() != static_cast<std::size_t>(9 * nc)) {
    M.P.assign(9 * nc, 0.0f);
    for (auto* v : {&M.nn, &M.nfx, &M.nfy, &M.nfz, &M.ne}) v->assign(nc, 0.0f);
  }
  const bool otvet = s.cfg.rt_method == RTMethod::OTVET;
  if (otvet && M.h.size() != static_cast<std::size_t>(6 * nc)) otvet_update_tensor(s);
  // otvet_split = 1: N_dir,g = q_g ndir (otvet_closure.hpp)
  const float* ndir = (otvet && !M.ndir.empty()) ? M.ndir.data() : nullptr;
  std::array<double, kNfreq> qg{};
  for (int g = 0; g < kNfreq; ++g) qg[g] = ndir ? inject_per_L(s, g) : 0.0;
  const int nsub = s.cfg.moment_nsub;
  const double C = 1.0 / nsub;
  for (int g = 0; g < kNfreq; ++g) {
    const float* ng = &M.n[g * nc];
    bool any = false;
    for (std::int64_t c = 0; c < nc && !any; ++c) any = ng[c] > 0;
    if (!any) continue;
    for (int it = 0; it < nsub; ++it) transport_group(M, g, otvet, C, ndir, qg[g]);
  }
  if (s.cfg.moment_inject == 1) inject();
  lap(s.timers.trace1);
  if (s.tdiag && s.tdiag->acc) moment_tensor_diag(M, *s.tdiag, otvet, ndir, qg.data());  // diagnostic

  // Deposit for radiative_physics: ngam = N, nu = E (= nu*ngam, divided by
  // ngam in radiative_physics), dx = sum over groups of F
#pragma omp parallel for schedule(static) if (par::worth(nc))
  for (std::int64_t c = 0; c < nc; ++c)
    mk::deposit_cell(&M.n[c], &M.e[c], &M.fx[c], &M.fy[c], &M.fz[c], nc, grid_ref(G, c));
  (void)Nrt;
  lap(s.timers.sync);
}

void moment_step_3(Simulation& s) {
  const double tw0 = wall_seconds();
  MomentState& M = *s.moments;
  const RTGrid& G = s.rt;
  const int N = M.N;
  const std::int64_t nc = M.ncell;
  const double absorbed = plane_sum(N, [&](int k) {
    double t = 0;
    const std::int64_t c0 = static_cast<std::int64_t>(k) * N * N;
    for (std::int64_t c = c0; c < c0 + static_cast<std::int64_t>(N) * N; ++c)
      for (int g = 0; g < kNfreq; ++g) {
        const std::int64_t i = g * nc + c;
        if (M.n[i] == 0) continue;
        t += mk::absorb_group(M.n[i], M.e[i], M.fx[i], M.fy[i], M.fz[i], G.fgam[c * kNfreq + g]);
      }
    return t;
  });
  M.absorbed += absorbed / s.cfg.nrtcr;
  s.timers.trace3 += wall_seconds() - tw0;
}

void moment_redshift(Simulation& s) {
  if (!s.moments) return;
  MomentState& M = *s.moments;
  const int N = M.N;
  const std::int64_t nc = M.ncell;
  // (a - da)/a, or per distance with redshift_per_distance (cosmology.hpp)
  const double fscale = redshift_fscale(s.cfg, s.clk);
  const RedshiftTable T = make_redshift_table();
  const double ngam_min = s.units.ngam_min;
  const double lost = plane_sum(N, [&](int k) {
    double t = 0;
    const std::int64_t c0 = static_cast<std::int64_t>(k) * N * N;
    for (std::int64_t c = c0; c < c0 + static_cast<std::int64_t>(N) * N; ++c) {
      double l = 0;
      if (mk::redshift_cell(T, fscale, ngam_min, &M.n[c], &M.e[c], &M.fx[c], &M.fy[c], &M.fz[c], nc, l)) t += l;
    }
    return t;
  });
  M.lost += lost / s.cfg.nrtcr;
  if (s.raylog)
    *s.raylog << " Moments    : in flight " << static_cast<float>(moment_photons_total(s) / s.cfg.nrtcr)
              << " injected " << static_cast<float>(M.injected) << " absorbed "
              << static_cast<float>(M.absorbed) << " lost " << static_cast<float>(M.lost)
              << "\n";
  if (s.raylog && s.recomb && s.recomb->active())
    *s.raylog << " Recomb     : emitted " << static_cast<float>(s.recomb->launched)
              << " (included in injected)\n";
}

}  // namespace radray
