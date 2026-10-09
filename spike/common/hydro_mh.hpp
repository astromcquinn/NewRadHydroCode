// Unsplit MUSCL-Hancock hydro step for the framework spike (ideal gas,
// gamma = 5/3, double precision, periodic box, 2 ghost cells).
//
//   1. primitives W = (rho, vx, vy, vz, p) and MC-limited slopes in each
//      direction (cells i-1, i, i+1);
//   2. Hancock predictor: W_half = W - dt/(2 dx) sum_d A_d(W) dW_d;
//   3. face states W_half +- dW_d/2 and an HLLC flux at every face;
//   4. U_new = U - dt/dx sum_d (F_{d,+} - F_{d,-}).
//
// Three passes: predict_cell (interior + 1 ghost layer), face_flux (every
// face of the interior), update_cell. Each value depends only on U in a
// fixed stencil and is computed by the same function in every framework, so
// the update is bit-identical for any decomposition and thread count. This is the solver family the roadmap proposes for AMR
// (docs/ROADMAP.md section 5); it is NOT the original's split TVD scheme.
#pragma once

#include "rr_hd.hpp"

namespace rrh {

constexpr int NV = 5;  // rho, mx, my, mz, E
constexpr int NG = 2;  // ghost cells needed
constexpr double kGamma = 5.0 / 3.0;
constexpr double kSmallRho = 1e-10, kSmallP = 1e-12;

RR_HD void cons_to_prim(const double U[NV], double W[NV]) {
  const double rho = rrf::max(U[0], kSmallRho);
  W[0] = rho;
  W[1] = U[1] / rho;
  W[2] = U[2] / rho;
  W[3] = U[3] / rho;
  const double ke = 0.5 * rho * ((W[1] * W[1] + W[2] * W[2]) + W[3] * W[3]);
  W[4] = rrf::max((kGamma - 1) * (U[4] - ke), kSmallP);
}
RR_HD void prim_to_cons(const double W[NV], double U[NV]) {
  U[0] = W[0];
  U[1] = W[0] * W[1];
  U[2] = W[0] * W[2];
  U[3] = W[0] * W[3];
  U[4] = W[4] / (kGamma - 1) + 0.5 * W[0] * ((W[1] * W[1] + W[2] * W[2]) + W[3] * W[3]);
}
RR_HD double sound_speed(const double W[NV]) { return std::sqrt(kGamma * W[4] / W[0]); }

// Monotonised-central limited slope
RR_HD double mc_slope(double wm, double w0, double wp) {
  const double dl = w0 - wm, dr = wp - w0;
  if (dl * dr <= 0) return 0;
  const double dc = 0.5 * (dl + dr);
  const double s = dc > 0 ? 1.0 : -1.0;
  return s * rrf::min(std::fabs(dc), rrf::min(2 * std::fabs(dl), 2 * std::fabs(dr)));
}

// Physical flux along direction d (velocity index 1+d) of primitive W
RR_HD void phys_flux(const double W[NV], int d, double F[NV]) {
  const double vn = W[1 + d];
  double U[NV];
  prim_to_cons(W, U);
  F[0] = U[0] * vn;
  F[1] = U[1] * vn;
  F[2] = U[2] * vn;
  F[3] = U[3] * vn;
  F[1 + d] += W[4];
  F[4] = (U[4] + W[4]) * vn;
}

// HLLC (Toro 2009, section 10.4) with Davis wave-speed estimates
RR_HD void hllc(const double WL[NV], const double WR[NV], int d, double F[NV]) {
  const double cL = sound_speed(WL), cR = sound_speed(WR);
  const double uL = WL[1 + d], uR = WR[1 + d];
  const double SL = rrf::min(uL - cL, uR - cR);
  const double SR = rrf::max(uL + cL, uR + cR);
  const double rL = WL[0], rR = WR[0], pL = WL[4], pR = WR[4];
  const double SM = (pR - pL + rL * uL * (SL - uL) - rR * uR * (SR - uR)) /
                    (rL * (SL - uL) - rR * (SR - uR));
  if (SL >= 0) {
    phys_flux(WL, d, F);
    return;
  }
  if (SR <= 0) {
    phys_flux(WR, d, F);
    return;
  }
  const double* W = SM >= 0 ? WL : WR;
  const double S = SM >= 0 ? SL : SR;
  double Uk[NV], Fk[NV], Us[NV];
  prim_to_cons(W, Uk);
  phys_flux(W, d, Fk);
  const double uk = W[1 + d];
  const double fac = W[0] * (S - uk) / (S - SM);
  Us[0] = fac;
  Us[1] = fac * W[1];
  Us[2] = fac * W[2];
  Us[3] = fac * W[3];
  Us[1 + d] = fac * SM;
  Us[4] = fac * (Uk[4] / W[0] + (SM - uk) * (SM + W[4] / (W[0] * (S - uk))));
  for (int m = 0; m < NV; ++m) F[m] = Fk[m] + S * (Us[m] - Uk[m]);
}

// Hancock half-step state of the cell whose primitives are Wc, given the
// neighbours along each direction: Wn[d][0] = cell-1, Wn[d][1] = cell+1.
// Returns the slopes too.
RR_HD void hancock(const double Wc[NV], const double Wn[3][2][NV], double dtdx,
                   double Wh[NV], double dW[3][NV]) {
  for (int d = 0; d < 3; ++d)
    for (int m = 0; m < NV; ++m) dW[d][m] = mc_slope(Wn[d][0][m], Wc[m], Wn[d][1][m]);
  const double rho = Wc[0], p = Wc[4];
  for (int m = 0; m < NV; ++m) Wh[m] = Wc[m];
  for (int d = 0; d < 3; ++d) {
    const double un = Wc[1 + d];
    const double* s = dW[d];
    // A_d(W) dW for primitive Euler equations
    double AdW[NV];
    AdW[0] = un * s[0] + rho * s[1 + d];
    AdW[1] = un * s[1];
    AdW[2] = un * s[2];
    AdW[3] = un * s[3];
    AdW[1 + d] += s[4] / rho;
    AdW[4] = un * s[4] + kGamma * p * s[1 + d];
    for (int m = 0; m < NV; ++m) Wh[m] -= 0.5 * dtdx * AdW[m];
  }
  Wh[0] = rrf::max(Wh[0], kSmallRho);
  Wh[4] = rrf::max(Wh[4], kSmallP);
}

// Per-cell predictor: Q = (W_half(5), dW_x(5), dW_y(5), dW_z(5)), from U at
// the cell and its 6 face neighbours. U(i,j,k,m) is the caller's accessor
// (local index space, ghosts included). Computed on the interior plus one
// ghost layer, it needs 2 ghost cells of U.
constexpr int NQ = 4 * NV;
template <class UA>
RR_HD void predict_cell(const UA& U, int i, int j, int k, double dtdx, double Q[NQ]) {
  double Uc[NV], Wc[NV], Wn[3][2][NV];
  for (int m = 0; m < NV; ++m) Uc[m] = U(i, j, k, m);
  cons_to_prim(Uc, Wc);
  for (int e = 0; e < 3; ++e) {
    const int ei = e == 0, ej = e == 1, ek = e == 2;
    for (int s = 0; s < 2; ++s) {
      const int o = 2 * s - 1;
      double Un[NV];
      for (int m = 0; m < NV; ++m) Un[m] = U(i + o * ei, j + o * ej, k + o * ek, m);
      cons_to_prim(Un, Wn[e][s]);
    }
  }
  double dW[3][NV];
  hancock(Wc, Wn, dtdx, Q, dW);
  for (int d = 0; d < 3; ++d)
    for (int m = 0; m < NV; ++m) Q[NV * (1 + d) + m] = dW[d][m];
}

// HLLC flux through the face between cell L (below) and cell R (above) in
// direction d, from their predictor states.
RR_HD void face_flux(const double QL[NQ], const double QR[NQ], int d, double F[NV]) {
  double WL[NV], WR[NV];
  for (int m = 0; m < NV; ++m) {
    WL[m] = QL[m] + 0.5 * QL[NV * (1 + d) + m];
    WR[m] = QR[m] - 0.5 * QR[NV * (1 + d) + m];
  }
  WL[0] = rrf::max(WL[0], kSmallRho);
  WL[4] = rrf::max(WL[4], kSmallP);
  WR[0] = rrf::max(WR[0], kSmallRho);
  WR[4] = rrf::max(WR[4], kSmallP);
  hllc(WL, WR, d, F);
}

// Conservative update of one cell from its 6 face fluxes:
// F[d][0] = lower face, F[d][1] = upper face.
RR_HD void update_cell(double U[NV], const double* Fx0, const double* Fx1, const double* Fy0,
                       const double* Fy1, const double* Fz0, const double* Fz1, double dtdx) {
  for (int c = 0; c < NV; ++c)
    U[c] = U[c] - dtdx * (((Fx1[c] - Fx0[c]) + (Fy1[c] - Fy0[c])) + (Fz1[c] - Fz0[c]));
}

// Maximum signal speed of a cell, for the CFL time step
RR_HD double max_speed(const double U[NV]) {
  double W[NV];
  cons_to_prim(U, W);
  const double c = sound_speed(W);
  return rrf::max(rrf::max(std::fabs(W[1]), std::fabs(W[2])), std::fabs(W[3])) + c;
}

// Analytic initial condition "blast3d" at the centre of global cell (i,j,k)
// of an n^3 periodic unit box: a smooth density wave, a shear flow and a
// Gaussian over-pressure (a 3-D Sedov-like blast), so every flux path and
// all three directions are exercised. kind = 1: Sod tube along x.
RR_HD void initial_state(int kind, int i, int j, int k, int n, double U[NV]) {
  const double pi = 3.141592653589793238462643383279502884;
  const double x = (i + 0.5) / n, y = (j + 0.5) / n, z = (k + 0.5) / n;
  double W[NV];
  if (kind == 1) {
    const bool left = x < 0.25 || x >= 0.75;  // periodic: two interfaces
    W[0] = left ? 1.0 : 0.125;
    W[1] = W[2] = W[3] = 0;
    W[4] = left ? 1.0 : 0.1;
  } else {
    const double r2 = (x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5) + (z - 0.5) * (z - 0.5);
    W[0] = 1 + 0.3 * std::sin(2 * pi * x) * std::sin(4 * pi * y) + 0.2 * std::cos(2 * pi * z);
    W[1] = 0.5 * std::sin(2 * pi * y);
    W[2] = 0.3 * std::sin(2 * pi * z);
    W[3] = 0.2 * std::cos(2 * pi * x);
    W[4] = 0.1 + 10 * std::exp(-r2 / (0.08 * 0.08));
  }
  prim_to_cons(W, U);
}

}  // namespace rrh
