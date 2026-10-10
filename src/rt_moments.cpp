// Moment-based radiative transfer (rt_moments.hpp, docs/rt_moments.md).
#include "radray/rt_moments.hpp"
#include "radray/recomb_rad.hpp"

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
  if (set_recomb_param(c, key, val)) return true;  // recomb_rad.hpp
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
  copy_recomb_param(from, to);
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
    s += ")";
  }
  // recomb_rad.hpp; printed only when not the original's case B
  if (c.recomb_case_A) s += ", recomb_case = A";
  if (c.recomb_rad)
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

void moment_prep(Simulation& s) {
  moment_init(s);
  MomentState& M = *s.moments;
  if (s.cfg.rt_method == RTMethod::OTVET) {
    M.tensor_method = compute_otvet_tensor(s.rt, s.units, s.cfg.otvet_tensor, M.h, &M.nsrc);
  } else {
    M.h.clear();
  }
}

double moment_photons_total(const Simulation& s) {
  if (!s.moments) return 0;
  const MomentState& M = *s.moments;
  double t = 0;
  for (std::int64_t c = 0; c < M.ncell; ++c) {
    double sc = 0;
    for (int g = 0; g < kNfreq; ++g) sc += static_cast<double>(M.n[g * M.ncell + c]);
    t = t + sc;
  }
  return t;
}

// ---------------------------------------------------------------------------
// Transport
// ---------------------------------------------------------------------------

namespace {

// Pressure tensor P = D N (xx, yy, zz, xy, xz, yz) of the M1 closure
// (Levermore 1984): chi = (3 + 4f^2)/(5 + 2 sqrt(4 - 3f^2)), f = |F|/N,
// D = (1-chi)/2 I + (3chi-1)/2 u u^T with u = F/|F|.
inline void m1_pressure(double N, double fx, double fy, double fz, double P[6]) {
  if (!(N > 0)) {
    for (int m = 0; m < 6; ++m) P[m] = 0;
    return;
  }
  const double fn = std::sqrt((fx * fx + fy * fy) + fz * fz);
  const double f = std::min(fn / N, 1.0);
  const double chi = (3 + 4 * f * f) / (5 + 2 * std::sqrt(4 - 3 * f * f));
  const double a = (1 - chi) / 2, b = (3 * chi - 1) / 2;
  if (fn > 0) {
    const double ux = fx / fn, uy = fy / fn, uz = fz / fn;
    P[0] = N * (a + b * ux * ux);
    P[1] = N * (a + b * uy * uy);
    P[2] = N * (a + b * uz * uz);
    P[3] = N * b * ux * uy;
    P[4] = N * b * ux * uz;
    P[5] = N * b * uy * uz;
  } else {
    P[0] = P[1] = P[2] = N / 3;
    P[3] = P[4] = P[5] = 0;
  }
}

// One GLF subcycle of group g with Courant number C (periodic box).
//
// GLF (Rusanov) interface flux with wave speed c_hat: Phi_{i+1/2} =
// (f(U_i) + f(U_{i+1}))/2 - (U_{i+1} - U_i)/2 in cell units, so that
//   U_i' = U_i - C sum_d [ (f_d(U_R) - f_d(U_L))/2 - (U_R - 2 U_i + U_L)/2 ].
// The fluxes along d of (N, F_x, F_y, F_z, E) are (F_d, P_xd, P_yd, P_zd,
// nu F_d), P = D N. Written as the sum of the right- and left-going parts
// (N +- F_d)/2 >= 0 of each cell, the update is a positive combination for
// C <= 1/3 (N stays >= 0 and E/N stays within the neighbours' range).
void transport_group(MomentState& M, int g, bool otvet, double C) {
  const int N = M.N;
  const std::int64_t nc = M.ncell;
  float* n = &M.n[g * nc];
  float* fx = &M.fx[g * nc];
  float* fy = &M.fy[g * nc];
  float* fz = &M.fz[g * nc];
  float* e = &M.e[g * nc];
  // Per cell: P (xx, yy, zz, xy, xz, yz) and the E fluxes nu F (x, y, z)
  float* P = M.P.data();

#pragma omp parallel for schedule(static) if (par::worth(nc))
  for (std::int64_t c = 0; c < nc; ++c) {
    float* Pc = &P[9 * c];
    const double Nc = n[c];
    if (!(Nc > 0)) {
      for (int m = 0; m < 9; ++m) Pc[m] = 0;
      continue;
    }
    double Q[6];
    if (otvet) {
      const float* hc = &M.h[6 * c];
      for (int m = 0; m < 6; ++m) Q[m] = Nc * static_cast<double>(hc[m]);
    } else {
      m1_pressure(Nc, fx[c], fy[c], fz[c], Q);
    }
    for (int m = 0; m < 6; ++m) Pc[m] = static_cast<float>(Q[m]);
    const double nu = static_cast<double>(e[c]) / Nc;
    Pc[6] = static_cast<float>(nu * fx[c]);
    Pc[7] = static_cast<float>(nu * fy[c]);
    Pc[8] = static_cast<float>(nu * fz[c]);
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
        const float n0 = n[c];
        if (n0 == 0 && n[xl] == 0 && n[xr] == 0 && n[yl] == 0 && n[yr] == 0 && n[zl] == 0 &&
            n[zr] == 0) {
          nn[c] = nfx[c] = nfy[c] = nfz[c] = ne[c] = 0;
          continue;
        }
        // flux differences and second differences, axis by axis
        const double sN = ((double(n[xr]) + n[xl]) + (double(n[yr]) + n[yl])) +
                          (double(n[zr]) + n[zl]) - 6.0 * n0;
        const double dN = (double(fx[xr]) - fx[xl]) + (double(fy[yr]) - fy[yl]) +
                          (double(fz[zr]) - fz[zl]);
        const float *Pxl = &P[9 * xl], *Pxr = &P[9 * xr], *Pyl = &P[9 * yl], *Pyr = &P[9 * yr],
                    *Pzl = &P[9 * zl], *Pzr = &P[9 * zr];
        // F_x: P_xx along x, P_xy along y, P_xz along z (P_kd)
        const double dFx = (double(Pxr[0]) - Pxl[0]) + (double(Pyr[3]) - Pyl[3]) +
                           (double(Pzr[4]) - Pzl[4]);
        const double dFy = (double(Pxr[3]) - Pxl[3]) + (double(Pyr[1]) - Pyl[1]) +
                           (double(Pzr[5]) - Pzl[5]);
        const double dFz = (double(Pxr[4]) - Pxl[4]) + (double(Pyr[5]) - Pyl[5]) +
                           (double(Pzr[2]) - Pzl[2]);
        const double dE = (double(Pxr[6]) - Pxl[6]) + (double(Pyr[7]) - Pyl[7]) +
                          (double(Pzr[8]) - Pzl[8]);
        auto lap = [&](const float* v) {
          return ((double(v[xr]) + v[xl]) + (double(v[yr]) + v[yl])) + (double(v[zr]) + v[zl]) -
                 6.0 * v[c];
        };
        double U[5];
        U[0] = n0 - C * (0.5 * dN - 0.5 * sN);
        U[1] = fx[c] - C * (0.5 * dFx - 0.5 * lap(fx));
        U[2] = fy[c] - C * (0.5 * dFy - 0.5 * lap(fy));
        U[3] = fz[c] - C * (0.5 * dFz - 0.5 * lap(fz));
        U[4] = e[c] - C * (0.5 * dE - 0.5 * lap(e));
        // realizability: N >= 0, |F| <= N, E >= 0
        if (!(U[0] > 0)) {
          U[0] = U[1] = U[2] = U[3] = U[4] = 0;
        } else {
          const double fn = std::sqrt((U[1] * U[1] + U[2] * U[2]) + U[3] * U[3]);
          if (fn > U[0]) {
            const double r = U[0] / fn;
            U[1] *= r;
            U[2] *= r;
            U[3] *= r;
          }
          if (U[4] < 0) U[4] = 0;
        }
        nn[c] = static_cast<float>(U[0]);
        nfx[c] = static_cast<float>(U[1]);
        nfy[c] = static_cast<float>(U[2]);
        nfz[c] = static_cast<float>(U[3]);
        ne[c] = static_cast<float>(U[4]);
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
#pragma omp parallel for schedule(static) if (par::worth(nc))
  for (std::int64_t c = 0; c < nc; ++c) {
    for (int q = 0; q < kNfreq; ++q) {
      G.nu[c * kNfreq + q] = 0;
      G.ngam[c * kNfreq + q] = 0;
      G.fgam[c * kNfreq + q] = 0;
    }
    for (int d = 0; d < 3; ++d) {
      G.gamma[c * 3 + d] = 0;
      G.lambda[c * 3 + d] = 0;
      G.dx[c * 3 + d] = 0;
    }
    const double rhog = static_cast<double>(G.rhog[c]) * u.rhob_avg;
    const double nH = rhog * XHy / mH_cgs;
    const double nHe = rhog * YHe / mHe_cgs;
    const double nHI = nH * static_cast<double>(G.fHI[c]);
    const double nHeI = nHe * static_cast<double>(G.fHeI[c]);
    const double nHeII = nHe * static_cast<double>(G.fHeII[c]);
    float* tau = G.tau_p(c);
    for (int q = 0; q < kNfreq; ++q) {
      const double sm = ((0.0 + nHI * sig[0][q]) + nHeI * sig[1][q]) + nHeII * sig[2][q];
      tau[q] = static_cast<float>(sm * u.len_rt);
    }
  }
  lap(s.timers.calc_rt1);

  // Sources: the photons ray_source emits, into the source cell, isotropic
  double inj = 0;
  for (std::int64_t c = 0; c < nc; ++c) {
    if (!(static_cast<double>(G.sfr[c]) > u.sfr_min)) continue;
    const double Ldt = static_cast<double>(G.Luv[c]) * s.clk.dtsec_rt;
    for (int g = 0; g < kNfreq; ++g) {
      double ngam = fResc * pop2_init[g] * Ldt;
      ngam = ngam / u.vol_rt / u.nH_avg;
      const std::int64_t i = g * nc + c;
      const float n0 = M.n[i];
      M.n[i] = static_cast<float>(static_cast<double>(n0) + ngam);
      M.e[i] = static_cast<float>(static_cast<double>(M.e[i]) + ngam * freq_init[g]);
      inj += static_cast<double>(M.n[i]) - static_cast<double>(n0);
    }
  }
  M.injected += inj / s.cfg.nrtcr;

  // Transport: moment_nsub GLF subcycles of Courant number 1/moment_nsub
  // (one cell of light travel per RT substep)
  if (M.P.size() != static_cast<std::size_t>(9 * nc)) {
    M.P.assign(9 * nc, 0.0f);
    for (auto* v : {&M.nn, &M.nfx, &M.nfy, &M.nfz, &M.ne}) v->assign(nc, 0.0f);
  }
  const bool otvet = s.cfg.rt_method == RTMethod::OTVET;
  if (otvet && M.h.size() != static_cast<std::size_t>(6 * nc))
    M.tensor_method = compute_otvet_tensor(G, u, s.cfg.otvet_tensor, M.h, &M.nsrc);
  const int nsub = s.cfg.moment_nsub;
  const double C = 1.0 / nsub;
  for (int g = 0; g < kNfreq; ++g) {
    const float* ng = &M.n[g * nc];
    bool any = false;
    for (std::int64_t c = 0; c < nc && !any; ++c) any = ng[c] > 0;
    if (!any) continue;
    for (int it = 0; it < nsub; ++it) transport_group(M, g, otvet, C);
  }
  lap(s.timers.trace1);

  // Deposit for radiative_physics: ngam = N, nu = E (= nu*ngam, divided by
  // ngam in radiative_physics), dx = sum over groups of F
#pragma omp parallel for schedule(static) if (par::worth(nc))
  for (std::int64_t c = 0; c < nc; ++c) {
    double dx[3] = {0, 0, 0};
    for (int g = 0; g < kNfreq; ++g) {
      const std::int64_t i = g * nc + c;
      G.ngam[c * kNfreq + g] = M.n[i];
      G.nu[c * kNfreq + g] = M.e[i];
      dx[0] += static_cast<double>(M.fx[i]);
      dx[1] += static_cast<double>(M.fy[i]);
      dx[2] += static_cast<double>(M.fz[i]);
    }
    for (int d = 0; d < 3; ++d) G.dx[c * 3 + d] = static_cast<float>(dx[d]);
  }
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
        const float fg = std::min(std::max(G.fgam[c * kNfreq + g], 0.0f), 1.0f);
        const float n0 = M.n[i];
        M.n[i] = n0 * fg;
        M.e[i] = M.e[i] * fg;
        M.fx[i] = M.fx[i] * fg;
        M.fy[i] = M.fy[i] * fg;
        M.fz[i] = M.fz[i] * fg;
        t += static_cast<double>(n0) - static_cast<double>(M.n[i]);
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
  const double fscale = (s.clk.scalefactor - s.clk.da_sim) / s.clk.scalefactor;
  const RedshiftTable T = make_redshift_table();
  const double ngam_min = s.units.ngam_min;
  const double lost = plane_sum(N, [&](int k) {
    double t = 0;
    const std::int64_t c0 = static_cast<std::int64_t>(k) * N * N;
    for (std::int64_t c = c0; c < c0 + static_cast<std::int64_t>(N) * N; ++c) {
      double nu[kNfreq], ng[kNfreq], keep[kNfreq], down[kNfreq], F[kNfreq][3];
      double before = 0;
      bool any = false;
      for (int g = 0; g < kNfreq; ++g) {
        const std::int64_t i = g * nc + c;
        ng[g] = M.n[i];
        nu[g] = M.n[i] > 0 ? static_cast<double>(M.e[i]) / ng[g] : freq_init[g];
        F[g][0] = M.fx[i];
        F[g][1] = M.fy[i];
        F[g][2] = M.fz[i];
        before += ng[g];
        any = any || ng[g] > 0;
      }
      if (!any) continue;
      rrk::redshift_spectrum(T, fscale, ngam_min, nu, ng, keep, down);
      double after = 0;
      for (int g = 0; g < kNfreq; ++g) {
        const std::int64_t i = g * nc + c;
        double Fg[3];
        for (int d = 0; d < 3; ++d)
          Fg[d] = keep[g] * F[g][d] + (g + 1 < kNfreq ? down[g + 1] * F[g + 1][d] : 0.0);
        const float nf = static_cast<float>(ng[g]);
        if (!(nf > 0)) {
          M.n[i] = M.e[i] = M.fx[i] = M.fy[i] = M.fz[i] = 0;
          continue;
        }
        const double fn = std::sqrt((Fg[0] * Fg[0] + Fg[1] * Fg[1]) + Fg[2] * Fg[2]);
        const double r = fn > ng[g] ? ng[g] / fn : 1.0;
        M.n[i] = nf;
        M.e[i] = static_cast<float>(nu[g] * ng[g]);
        M.fx[i] = static_cast<float>(Fg[0] * r);
        M.fy[i] = static_cast<float>(Fg[1] * r);
        M.fz[i] = static_cast<float>(Fg[2] * r);
        after += static_cast<double>(nf);
      }
      t += before - after;
    }
    return t;
  });
  M.lost += lost / s.cfg.nrtcr;
  if (s.raylog)
    *s.raylog << " Moments    : in flight " << static_cast<float>(moment_photons_total(s) / s.cfg.nrtcr)
              << " injected " << static_cast<float>(M.injected) << " absorbed "
              << static_cast<float>(M.absorbed) << " lost " << static_cast<float>(M.lost)
              << "\n";
}

}  // namespace radray
