// Per-cell kernels of the moment RT (rt_moments.hpp, docs/rt_moments.md),
// written once for the serial solver (src/rt_moments.cpp, structure-of-
// arrays grids) and the AMReX level (src/amrex/moment_level.cpp, FabArrays
// on z-slabs; docs/stage_e_moments.md): the same arithmetic, so the two are
// bitwise identical. Host only (CPU; the moment level has no GPU path).
//
// The accessors keep the callers' layouts: a stencil accessor at(q, s)
// returns component q (0 n, 1 fx, 2 fy, 3 fz, 4 e) at stencil point s
// (kS0 the cell, kXL/kXR, kYL/kYR, kZL/kZR its face neighbours, periodic),
// and Pat(s, m) the transport flux m (0..8) of point s.
#pragma once

#include <algorithm>
#include <cmath>

#include "radray/constants.hpp"
#include "radray/kernels/ray_kernels.hpp"
#include "radray/kernels/rt_cell.hpp"

namespace radray {
namespace mk {

constexpr int kS0 = 0, kXL = 1, kXR = 2, kYL = 3, kYR = 4, kZL = 5, kZR = 6;

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

// otvet_split = 1 (otvet_closure.hpp): P = N I/3 + min(N, N_dir)(h - I/3)
inline void otvet_split_pressure(double Nc, double Nd, const float* hc, double Q[6]) {
  const double m = std::min(Nc, Nd);
  const double a = (Nc - m) / 3;
  for (int k = 0; k < 3; ++k) Q[k] = a + m * static_cast<double>(hc[k]);
  for (int k = 3; k < 6; ++k) Q[k] = m * static_cast<double>(hc[k]);
}

// The closure of a cell: M1, OTVET (P = N h), the split (h and N_dir =
// nd = q_g ndir) or isotropic (P = N I/3: the background field of
// otvet_background = ga01, docs/otvet_background.md; hc unused)
enum class Closure { M1, Otvet, Split, Iso };

inline void closure_pressure(Closure cl, double Nc, float fx, float fy, float fz, const float* hc, double nd,
                             double Q[6]) {
  if (cl == Closure::Iso) {
    Q[0] = Q[1] = Q[2] = Nc / 3;
    Q[3] = Q[4] = Q[5] = 0;
  } else if (cl == Closure::Split) {
    otvet_split_pressure(Nc, nd, hc, Q);
  } else if (cl == Closure::Otvet) {
    for (int m = 0; m < 6; ++m) Q[m] = Nc * static_cast<double>(hc[m]);
  } else {
    m1_pressure(Nc, fx, fy, fz, Q);
  }
}

// The 9 transport fluxes of a cell: P (xx, yy, zz, xy, xz, yz) and the E
// fluxes nu F (x, y, z)
inline void flux_cell(Closure cl, float n, float fx, float fy, float fz, float e, const float* hc, double nd,
                      float Pc[9]) {
  const double Nc = n;
  if (!(Nc > 0)) {
    for (int m = 0; m < 9; ++m) Pc[m] = 0;
    return;
  }
  double Q[6];
  closure_pressure(cl, Nc, fx, fy, fz, hc, nd, Q);
  for (int m = 0; m < 6; ++m) Pc[m] = static_cast<float>(Q[m]);
  const double nu = static_cast<double>(e) / Nc;
  Pc[6] = static_cast<float>(nu * fx);
  Pc[7] = static_cast<float>(nu * fy);
  Pc[8] = static_cast<float>(nu * fz);
}

// One GLF (Rusanov) subcycle of one cell with Courant number C (rt_moments
// transport_group): out = the new (n, fx, fy, fz, e).
//
// GLF interface flux with wave speed c_hat: Phi_{i+1/2} = (f(U_i) +
// f(U_{i+1}))/2 - (U_{i+1} - U_i)/2 in cell units, so that
//   U_i' = U_i - C sum_d [ (f_d(U_R) - f_d(U_L))/2 - (U_R - 2 U_i + U_L)/2 ].
// The fluxes along d of (N, F_x, F_y, F_z, E) are (F_d, P_xd, P_yd, P_zd,
// nu F_d), P = D N. Written as the sum of the right- and left-going parts
// (N +- F_d)/2 >= 0 of each cell, the update is a positive combination for
// C <= 1/3 (N stays >= 0 and E/N stays within the neighbours' range).
template <class At, class PAt>
inline void glf_cell(const At& at, const PAt& Pat, double C, float out[5]) {
  const float n0 = at(0, kS0);
  if (n0 == 0 && at(0, kXL) == 0 && at(0, kXR) == 0 && at(0, kYL) == 0 && at(0, kYR) == 0 &&
      at(0, kZL) == 0 && at(0, kZR) == 0) {
    for (int q = 0; q < 5; ++q) out[q] = 0;
    return;
  }
  // flux differences and second differences, axis by axis
  const double sN = ((double(at(0, kXR)) + at(0, kXL)) + (double(at(0, kYR)) + at(0, kYL))) +
                    (double(at(0, kZR)) + at(0, kZL)) - 6.0 * n0;
  const double dN = (double(at(1, kXR)) - at(1, kXL)) + (double(at(2, kYR)) - at(2, kYL)) +
                    (double(at(3, kZR)) - at(3, kZL));
  // F_x: P_xx along x, P_xy along y, P_xz along z (P_kd)
  auto dP = [&](int mx, int my, int mz) {
    return (double(Pat(kXR, mx)) - Pat(kXL, mx)) + (double(Pat(kYR, my)) - Pat(kYL, my)) +
           (double(Pat(kZR, mz)) - Pat(kZL, mz));
  };
  const double dFx = dP(0, 3, 4);
  const double dFy = dP(3, 1, 5);
  const double dFz = dP(4, 5, 2);
  const double dE = dP(6, 7, 8);
  auto lap = [&](int q) {
    return ((double(at(q, kXR)) + at(q, kXL)) + (double(at(q, kYR)) + at(q, kYL))) +
           (double(at(q, kZR)) + at(q, kZL)) - 6.0 * at(q, kS0);
  };
  double U[5];
  U[0] = n0 - C * (0.5 * dN - 0.5 * sN);
  U[1] = at(1, kS0) - C * (0.5 * dFx - 0.5 * lap(1));
  U[2] = at(2, kS0) - C * (0.5 * dFy - 0.5 * lap(2));
  U[3] = at(3, kS0) - C * (0.5 * dFz - 0.5 * lap(3));
  U[4] = at(4, kS0) - C * (0.5 * dE - 0.5 * lap(4));
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
  for (int q = 0; q < 5; ++q) out[q] = static_cast<float>(U[q]);
}

// moment_order = 2 (rt_moments.cpp moment_transport2): the GLF flux along
// axis d (0, 1, 2) at the face between cells b and c from the linearly
// reconstructed face states. U[q][p]: the 5 variables (n, fx, fy, fz, e) of
// the four cells a, b, c, e along the axis (p = 0..3: i-1, i, i+1, i+2); hb,
// hc and ndb, ndc: the closure data of cells b and c (flux_cell's hc, nd).
//
// Slopes: the monotonized central limiter, component by component, so a
// face value lies between the cell's and its neighbour's; then the face
// state is made realizable as glf_cell's result is (N > 0 or all zero,
// |F| <= N, E >= 0). Flux: Phi = (f_d(U_L) + f_d(U_R))/2 - (U_R - U_L)/2,
// f_d = (F_d, P_xd, P_yd, P_zd, nu F_d) with the closure of the face state
// (M1: of its N, F; OTVET and the split: h and N_dir of its own cell).
inline double mc_slope(double dl, double dr) {
  if (!(dl * dr > 0)) return 0.0;
  const double c = 0.5 * (dl + dr);
  const double m = 2.0 * (std::fabs(dl) < std::fabs(dr) ? dl : dr);
  return std::fabs(c) < std::fabs(m) ? c : m;
}

inline void face_flux(Closure cl, int d, double V[5], const float* hc, double nd, double f[5]) {
  // realizability of the face state
  if (!(V[0] > 0)) {
    for (int q = 0; q < 5; ++q) V[q] = f[q] = 0;
    return;
  }
  const double fn = std::sqrt((V[1] * V[1] + V[2] * V[2]) + V[3] * V[3]);
  if (fn > V[0]) {
    const double r = V[0] / fn;
    V[1] *= r;
    V[2] *= r;
    V[3] *= r;
  }
  if (V[4] < 0) V[4] = 0;
  double Q[6];
  closure_pressure(cl, V[0], static_cast<float>(V[1]), static_cast<float>(V[2]), static_cast<float>(V[3]), hc, nd,
                   Q);
  // P_kd: (xx, xy, xz), (xy, yy, yz), (xz, yz, zz)
  static const int col[3][3] = {{0, 3, 4}, {3, 1, 5}, {4, 5, 2}};
  f[0] = V[1 + d];
  f[1] = Q[col[d][0]];
  f[2] = Q[col[d][1]];
  f[3] = Q[col[d][2]];
  f[4] = V[4] / V[0] * V[1 + d];
}

inline void muscl_flux(Closure cl, int d, const float U[5][4], const float* hb, const float* hc, double ndb,
                       double ndc, float out[5]) {
  if (U[0][0] == 0 && U[0][1] == 0 && U[0][2] == 0 && U[0][3] == 0) {
    for (int q = 0; q < 5; ++q) out[q] = 0;
    return;
  }
  double L[5], R[5], fl[5], fr[5];
  for (int q = 0; q < 5; ++q) {
    const double a = U[q][0], b = U[q][1], c = U[q][2], e = U[q][3];
    L[q] = b + 0.5 * mc_slope(b - a, c - b);
    R[q] = c - 0.5 * mc_slope(c - b, e - c);
  }
  face_flux(cl, d, L, hb, ndb, fl);
  face_flux(cl, d, R, hc, ndc, fr);
  for (int q = 0; q < 5; ++q) out[q] = static_cast<float>(0.5 * (fl[q] + fr[q]) - 0.5 * (R[q] - L[q]));
}

// The realizability of a cell after a stage of moment_transport2, as at the
// end of glf_cell
inline void realizable_cell(float& n, float& fx, float& fy, float& fz, float& e) {
  if (!(n > 0)) {
    n = fx = fy = fz = e = 0;
    return;
  }
  const double fn = std::sqrt((double(fx) * fx + double(fy) * fy) + double(fz) * fz);
  if (fn > n) {
    const double r = n / fn;
    fx = static_cast<float>(fx * r);
    fy = static_cast<float>(fy * r);
    fz = static_cast<float>(fz * r);
  }
  if (e < 0) e = 0;
}

// calc_rt_1 of the moment solver (the B3/B4 fixes: dx reset and tau in every
// cell; tau is a diagnostic, the solver does not read it). sig =
// photoion_sigma(freq_init).
inline void calc_rt_1_cell(const rrk::CellRef& c, const radray::Units& u, const double sig[3][kNfreq]) {
  using namespace radray::pc;
  using namespace radray::cp;
  for (int q = 0; q < kNfreq; ++q) {
    c.nu[q] = 0;
    c.ngam[q] = 0;
    c.fgam[q] = 0;
  }
  for (int d = 0; d < 3; ++d) {
    c.gamma[d] = 0;
    c.lambda[d] = 0;
    c.dx[d] = 0;
  }
  const double rhog = static_cast<double>(*c.rhog) * u.rhob_avg;
  const double nH = rhog * XHy / mH_cgs;
  const double nHe = rhog * YHe / mHe_cgs;
  const double nHI = nH * static_cast<double>(*c.fHI);
  const double nHeI = nHe * static_cast<double>(*c.fHeI);
  const double nHeII = nHe * static_cast<double>(*c.fHeII);
  for (int q = 0; q < kNfreq; ++q) {
    const double sm = ((0.0 + nHI * sig[0][q]) + nHeI * sig[1][q]) + nHeII * sig[2][q];
    c.tau[q] = static_cast<float>(sm * u.len_rt);
  }
}

// The source photons of one cell (the number ray_source emits), into the
// cell, isotropic: group g at n[g * stride], e[g * stride]. inj += the
// photons added (float-rounded), group by group; rec (optional) receives
// the kNfreq addends.
inline void inject_cell(double Luv, double dtsec_rt, const radray::Units& u, float* n, float* e,
                        std::int64_t stride, double& inj, double* rec = nullptr) {
  using namespace radray::cp;
  const double Ldt = Luv * dtsec_rt;
  for (int g = 0; g < kNfreq; ++g) {
    double ngam = fResc * pop2_init[g] * Ldt;
    ngam = ngam / u.vol_rt / u.nH_avg;
    float& ng = n[g * stride];
    float& eg = e[g * stride];
    const float n0 = ng;
    ng = static_cast<float>(static_cast<double>(n0) + ngam);
    eg = static_cast<float>(static_cast<double>(eg) + ngam * freq_init[g]);
    const double add = static_cast<double>(ng) - static_cast<double>(n0);
    if (rec) rec[g] = add;
    inj += add;
  }
}

// The deposit for radiative_physics: ngam = N, nu = E (= nu*ngam, divided by
// ngam in radiative_physics), dx = the sum over groups of F
inline void deposit_cell(const float* n, const float* e, const float* fx, const float* fy, const float* fz,
                         std::int64_t stride, const rrk::CellRef& c) {
  double dx[3] = {0, 0, 0};
  for (int g = 0; g < kNfreq; ++g) {
    const std::int64_t i = g * stride;
    c.ngam[g] = n[i];
    c.nu[g] = e[i];
    dx[0] += static_cast<double>(fx[i]);
    dx[1] += static_cast<double>(fy[i]);
    dx[2] += static_cast<double>(fz[i]);
  }
  for (int d = 0; d < 3; ++d) c.dx[d] = static_cast<float>(dx[d]);
}

// otvet_background = ga01 (docs/otvet_background.md section 4a): move
// min(N_L, T) photons of one group of a cell from the local field L to the
// background B. They leave L in proportion to its (N, F, E) and arrive in B
// isotropically (F_B unchanged) with L's mean frequency. Returns the
// photons added to B (float-rounded).
inline double ga01_transfer(float& nL, float& fxL, float& fyL, float& fzL, float& eL, float& nB, float& eB,
                            double T) {
  const float n0 = nL;
  if (!(n0 > 0) || !(T > 0)) return 0;
  const float b0 = nB;
  if (T >= static_cast<double>(n0)) {
    nB = static_cast<float>(static_cast<double>(b0) + static_cast<double>(n0));
    eB = static_cast<float>(static_cast<double>(eB) + static_cast<double>(eL));
    nL = fxL = fyL = fzL = eL = 0;
  } else {
    const double r = 1 - T / static_cast<double>(n0);
    const float e0 = eL;
    nL = static_cast<float>(static_cast<double>(n0) * r);
    eL = static_cast<float>(static_cast<double>(e0) * r);
    fxL = static_cast<float>(static_cast<double>(fxL) * r);
    fyL = static_cast<float>(static_cast<double>(fyL) * r);
    fzL = static_cast<float>(static_cast<double>(fzL) * r);
    nB = static_cast<float>(static_cast<double>(b0) + (static_cast<double>(n0) - static_cast<double>(nL)));
    eB = static_cast<float>(static_cast<double>(eB) + (static_cast<double>(e0) - static_cast<double>(eL)));
  }
  return static_cast<double>(nB) - static_cast<double>(b0);
}

// The deposit of two fields (otvet_background = ga01): photo_cell sees the
// sum, ngam = N_L + N_B, nu = E_L + E_B, dx = the sum over groups and fields
// of F
inline void deposit_cell2(const float* n, const float* e, const float* fx, const float* fy, const float* fz,
                          const float* bn, const float* be, const float* bfx, const float* bfy, const float* bfz,
                          std::int64_t stride, const rrk::CellRef& c) {
  double dx[3] = {0, 0, 0};
  for (int g = 0; g < kNfreq; ++g) {
    const std::int64_t i = g * stride;
    c.ngam[g] = n[i] + bn[i];
    c.nu[g] = e[i] + be[i];
    dx[0] += static_cast<double>(fx[i]) + static_cast<double>(bfx[i]);
    dx[1] += static_cast<double>(fy[i]) + static_cast<double>(bfy[i]);
    dx[2] += static_cast<double>(fz[i]) + static_cast<double>(bfz[i]);
  }
  for (int d = 0; d < 3; ++d) c.dx[d] = static_cast<float>(dx[d]);
}

// moment_step_3 of one group of a cell with n != 0: N, F, E *= fgam
// (clamped to [0, 1]); returns the photons absorbed, n0 - n (doubles).
inline double absorb_group(float& n, float& e, float& fx, float& fy, float& fz, float fgam) {
  const float fg = std::min(std::max(fgam, 0.0f), 1.0f);
  const float n0 = n;
  n = n0 * fg;
  e = e * fg;
  fx = fx * fg;
  fy = fy * fg;
  fz = fz * fg;
  return static_cast<double>(n0) - static_cast<double>(n);
}

// moment_redshift of one cell: the group spectra over the hydro step (the
// redshift_spectrum of the rays; F follows the photons). Returns false for
// a cell without photons (untouched); otherwise lost = the photons before
// minus after.
inline bool redshift_cell(const rrk::RedshiftTable& T, double fscale, double ngam_min, float* n, float* e, float* fx,
                          float* fy, float* fz, std::int64_t stride, double& lost) {
  double nu[kNfreq], ng[kNfreq], keep[kNfreq], down[kNfreq], F[kNfreq][3];
  double before = 0;
  bool any = false;
  for (int g = 0; g < kNfreq; ++g) {
    const std::int64_t i = g * stride;
    ng[g] = n[i];
    nu[g] = n[i] > 0 ? static_cast<double>(e[i]) / ng[g] : radray::cp::freq_init[g];
    F[g][0] = fx[i];
    F[g][1] = fy[i];
    F[g][2] = fz[i];
    before += ng[g];
    any = any || ng[g] > 0;
  }
  if (!any) return false;
  rrk::redshift_spectrum(T, fscale, ngam_min, nu, ng, keep, down);
  double after = 0;
  for (int g = 0; g < kNfreq; ++g) {
    const std::int64_t i = g * stride;
    double Fg[3];
    for (int d = 0; d < 3; ++d) Fg[d] = keep[g] * F[g][d] + (g + 1 < kNfreq ? down[g + 1] * F[g + 1][d] : 0.0);
    const float nf = static_cast<float>(ng[g]);
    if (!(nf > 0)) {
      n[i] = e[i] = fx[i] = fy[i] = fz[i] = 0;
      continue;
    }
    const double fn = std::sqrt((Fg[0] * Fg[0] + Fg[1] * Fg[1]) + Fg[2] * Fg[2]);
    const double r = fn > ng[g] ? ng[g] / fn : 1.0;
    n[i] = nf;
    e[i] = static_cast<float>(nu[g] * ng[g]);
    fx[i] = static_cast<float>(Fg[0] * r);
    fy[i] = static_cast<float>(Fg[1] * r);
    fz[i] = static_cast<float>(Fg[2] * r);
    after += static_cast<double>(nf);
  }
  lost = before - after;
  return true;
}

// The photons of a cell over the groups (moment_photons_total's addend)
inline double photons_cell(const float* n, std::int64_t stride) {
  double sc = 0;
  for (int g = 0; g < kNfreq; ++g) sc += static_cast<double>(n[g * stride]);
  return sc;
}

// ---------------------------------------------------------------------------
// moment_inject = 2: the first cell of a source's light (docs/rt_moments.md
// section 10). The photons a steady source emits in one substep fill, at the
// end of it, the sphere of one cell of light travel around the source evenly
// in radius. The rays have exactly this: a ray emitted at a cell centre is
// cast one cell in the substep of its emission and shares its photons among
// the cells it crosses by its path length in each. The moment solver has no
// path: moment_inject = 1 counts the whole emission in the source cell for
// the substep, 0 lets the transport diffuse it away before anything absorbs.
// Here the emission enters as the cell averages of that sphere: in the source
// cell and its 26 neighbours the share w of the path r in [0, 1] cells from
// the centre that lies in the cell, averaged over the directions (0.6107 in
// the source cell, the mean distance from its centre to its faces; 0.3196 in
// the 6 face neighbours together, 0.0684 in the 12 edge and 0.0013 in the 8
// corner neighbours), with the flux m, the same average of the share times
// the direction: zero in the source cell, outward in the neighbours.
// ---------------------------------------------------------------------------
struct FirstCell {
  double w[27];     // offset (di, dj, dk) in {-1, 0, 1}^3 at (di + 1) + 3 (dj + 1) + 9 (dk + 1)
  double m[27][3];  // |m| <= w
};

// The weights by a midpoint rule over the directions through the faces of a
// cube (6 nq^2 directions p = (1, u, v) and permutations, of solid angle
// du dv / |p|^3). Only +, *, / and sqrt: the same bits on every platform.
inline FirstCell first_cell_weights(int nq = 512) {
  FirstCell W{};
  double sum = 0;
  const double h = 2.0 / nq;
  for (int face = 0; face < 6; ++face) {
    const int a = face / 2;
    const double sgn = face % 2 ? -1.0 : 1.0;
    for (int iu = 0; iu < nq; ++iu)
      for (int iv = 0; iv < nq; ++iv) {
        const double u = -1.0 + (iu + 0.5) * h, v = -1.0 + (iv + 0.5) * h;
        double p[3];
        p[a] = sgn;
        p[(a + 1) % 3] = u;
        p[(a + 2) % 3] = v;
        const double r2 = (p[0] * p[0] + p[1] * p[1]) + p[2] * p[2];
        const double r = std::sqrt(r2);
        const double wq = (h * h) / (r2 * r);
        const double d[3] = {p[0] / r, p[1] / r, p[2] / r};
        // the faces of the source cell that the path crosses before r = 1,
        // in the order of the crossings (0.5 / |d| < 1)
        double t[3];
        int ax[3], nt = 0;
        for (int q = 0; q < 3; ++q)
          if (std::fabs(d[q]) > 0.5) {
            t[nt] = 0.5 / std::fabs(d[q]);
            ax[nt] = q;
            ++nt;
          }
        for (int l1 = 1; l1 < nt; ++l1)
          for (int l = l1; l > 0 && t[l] < t[l - 1]; --l) {
            std::swap(t[l], t[l - 1]);
            std::swap(ax[l], ax[l - 1]);
          }
        int off[3] = {0, 0, 0};
        double t0 = 0;
        for (int l = 0; l <= nt; ++l) {
          const double len = (l < nt ? t[l] : 1.0) - t0;
          const int o = (off[0] + 1) + 3 * (off[1] + 1) + 9 * (off[2] + 1);
          W.w[o] += wq * len;
          for (int q = 0; q < 3; ++q) W.m[o][q] += wq * len * d[q];
          if (l < nt) {
            t0 = t[l];
            off[ax[l]] = d[ax[l]] > 0 ? 1 : -1;
          }
        }
        sum += wq;
      }
  }
  for (int o = 0; o < 27; ++o) {
    W.w[o] /= sum;
    for (int q = 0; q < 3; ++q) W.m[o][q] = o == 13 ? 0.0 : W.m[o][q] / sum;
  }
  return W;
}

// The photons of group g that a source cell of luminosity Luv emits in a
// substep (inject_cell's addend, the number ray_source emits)
inline double source_photons_group(double Luv, double dtsec_rt, const radray::Units& u, int g) {
  using namespace radray::cp;
  const double Ldt = Luv * dtsec_rt;
  double ngam = fResc * pop2_init[g] * Ldt;
  ngam = ngam / u.vol_rt / u.nH_avg;
  return ngam;
}

// The q photons of one group, of energy nu each, that the source in cell
// (i, j, k) of the periodic N^3 grid emits in a substep, added by W; n ..
// e: the group's fields. Returns the photons added (float-rounded).
inline double inject_first_cell(const FirstCell& W, int N, int i, int j, int k, double q, double nu, float* n,
                                float* fx, float* fy, float* fz, float* e) {
  double t = 0;
  for (int o = 0; o < 27; ++o) {
    const int ii = (i + (o % 3 - 1) + N) % N, jj = (j + ((o / 3) % 3 - 1) + N) % N, kk = (k + (o / 9 - 1) + N) % N;
    const std::int64_t c = ii + static_cast<std::int64_t>(N) * (jj + static_cast<std::int64_t>(N) * kk);
    const double a = W.w[o] * q;
    const float n0 = n[c];
    n[c] = static_cast<float>(static_cast<double>(n0) + a);
    e[c] = static_cast<float>(static_cast<double>(e[c]) + a * nu);
    if (o != 13) {
      fx[c] = static_cast<float>(static_cast<double>(fx[c]) + W.m[o][0] * q);
      fy[c] = static_cast<float>(static_cast<double>(fy[c]) + W.m[o][1] * q);
      fz[c] = static_cast<float>(static_cast<double>(fz[c]) + W.m[o][2] * q);
      realizable_cell(n[c], fx[c], fy[c], fz[c], e[c]);
    }
    t += static_cast<double>(n[c]) - static_cast<double>(n0);
  }
  return t;
}

}  // namespace mk
}  // namespace radray
