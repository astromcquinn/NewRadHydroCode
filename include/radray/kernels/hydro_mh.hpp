// Unsplit MUSCL-Hancock + HLLC hydro kernel (ideal gas, gamma = 5/3, double
// precision, periodic box, 2 ghost cells) with the dual-energy formulation
// of Bryan et al. (1995) and passively advected scalars. This is the
// `hydro_scheme = unsplit` solver of docs/hydro_unsplit.md, the solver family
// the roadmap proposes for AMR (docs/ROADMAP.md section 5); it is NOT the
// original's split moving-frame TVD scheme (src/hydro.cpp).
//
// One substep (dt, dx = 1 cell):
//   1. primitives W = (rho, vx, vy, vz, p) and MC-limited slopes in each
//      direction (cells i-1, i, i+1). The pressure comes from the thermal
//      energy selected by the dual-energy rule (thermal_energy below);
//   2. Hancock predictor: W_half = W - dt/(2 dx) sum_d A_d(W) dW_d;
//   3. face states W_half +- dW_d/2 and an HLLC flux at every face. The
//      face also returns the interface velocity v* and the upwind side,
//      which advect the scalars (the internal energy density and the four
//      species densities) with the HLLC mass flux;
//   4. U_new = U - dt/dx sum_d (F_{d,+} - F_{d,-}); the internal energy
//      gets the -p div(v) work from v* and the cell's half-step pressure;
//   5. finalize_cell: floors, the eta2 synchronisation of the internal
//      energy and the total-energy reset of cells whose pressure comes from
//      the internal energy (docs/hydro_unsplit.md section 2).
//
// Three passes: predict_cell (interior + 1 ghost layer), face_flux (every
// face of the interior), update_cell, then finalize_cell. Each value depends
// only on the state in a fixed stencil and is computed by the same function
// in every framework, so the update is bit-identical for any decomposition
// and thread count. Header-only and framework-independent on purpose (the
// AMReX port compiles the same header for its patches): no library state,
// host/device macros from rr_hd.hpp, fixed stencils. Gravity is not in the
// kernel: the host applies it as a kick between two half steps
// (docs/hydro_unsplit.md section 3).
//
// Provenance: spike/common/hydro_mh.hpp of the framework-spike branch
// (steps 1-4 without the scalars), extended here; this copy is the one
// maintained from now on.
#pragma once

#include "radray/kernels/rr_hd.hpp"

namespace rrh {

constexpr int NV = 5;  // rho, mx, my, mz, E
constexpr int NS = 5;  // scalars: internal energy density, then the 4 species w
constexpr int NG = 2;  // ghost cells needed
constexpr double kGamma = 5.0 / 3.0;
constexpr double kSmallRho = 1e-10, kSmallP = 1e-12;

// Dual-energy selection (Bryan et al. 1995, restated in Bryan et al. 2014
// section 4.1.1; eta1 = 1e-3 and eta2 = 0.1 are their values, the
// neighbourhood and the reset rule are this implementation's, see
// docs/hydro_unsplit.md section 2). eta1 <= 0 turns the formulation off
// (the thermal energy is always E - KE).
//   eta1: the pressure comes from E - KE when (E - KE)/E_ref > eta1, else
//         from the separately evolved internal energy. E_ref is the
//         maximum of E over the cell and its 6 face neighbours (the
//         truncation error of E comes from the fluxes, i.e. from the
//         neighbours' energies, so a stagnation cell between fast streams
//         is not fooled by its own small E), or the cell's E alone with
//         eta1_local (the cell-local form of Bryan et al. 1995, eq. 44 of
//         Bryan et al. 2014; see the doc);
//   eta2: the internal energy is synchronised to E - KE when
//         (E - KE)/max_neighbours(E) > eta2.
struct DualEnergy {
  double eta1 = 1e-3;
  double eta2 = 0.1;
  bool eta1_local = false;
  int nb_radius = 2;  // neighbourhood of the E maximum: +-1 (6 cells) or +-2 (12) along the axes
  // When the pressure comes from e, reset E := KE + e: 0 always; 1 only if
  // |E - KE - e| <= eta1 E (cell); 2 as 1, or if E - KE > E/2 (a thermal-
  // dominated stagnation cell); 3 only if E - KE <= eta1 max_{+-1}(E);
  // 4 never (Enzo: E is never modified)
  int reset_rule = 0;
  bool eta2_local = false;  // eta2 test against the cell's own E (Nyx's resync)
  // e_self: in cells whose pressure comes from E - KE the internal energy
  // is advected with its own specific value (reconstructed like the
  // species) and does the p div v work with its own pressure (e_ratio
  // below), so it stays on its own adiabat between eta2 synchronisations.
  // Off: with the face values and pressure of the selected energy E - KE,
  // so e is not evolved as itself there and drifts (docs/hydro_unsplit.md
  // section 7, U7).
  bool e_self = false;
};

RR_HD double kinetic_energy(const double U[NV]) {
  const double rho = rrf::max(U[0], kSmallRho);
  return 0.5 * ((U[1] * U[1] + U[2] * U[2]) + U[3] * U[3]) / rho;
}

// Thermal energy density of a cell: E - KE, or the evolved internal energy
// eD when E - KE is not a reliable fraction (eta1) of the reference energy
// Eref. Written so that a NaN or non-positive E - KE selects eD when the
// formulation is on.
RR_HD double select_thermal(const double U[NV], double eD, double Eref, double eta1) {
  const double eT = U[4] - kinetic_energy(U);
  if (eta1 > 0 && !(eT > eta1 * Eref)) return eD;
  return eT;
}
// Cell-local form (Eref = the cell's own E); used by the gravity kick
RR_HD double thermal_energy(const double U[NV], double eD, double eta1) {
  return select_thermal(U, eD, U[4], eta1);
}

// e_self: whether a cell's internal energy is evolved from its own value
// (its pressure comes from E - KE, es != eD), and the ratio of the evolved
// internal energy eD to the selected thermal energy es that scales the
// work term (exactly 1 where es = eD, so those cells are updated as
// without e_self, bit for bit).
RR_HD bool e_own(double eD, double es) { return es != eD; }
RR_HD double e_ratio(double eD, double es) { return es > 0 ? eD / es : 1.0; }

// Primitives from U and the (already selected) thermal energy density e
RR_HD void cons_to_prim(const double U[NV], double e, double W[NV]) {
  const double rho = rrf::max(U[0], kSmallRho);
  W[0] = rho;
  W[1] = U[1] / rho;
  W[2] = U[2] / rho;
  W[3] = U[3] / rho;
  W[4] = rrf::max((kGamma - 1) * e, kSmallP);
}
// Total-energy form (no dual energy), as in the spike
RR_HD void cons_to_prim(const double U[NV], double W[NV]) {
  cons_to_prim(U, U[4] - kinetic_energy(U), W);
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

// What the scalar advection and the dual-energy work term need from a face:
// the interface velocity v* (the contact speed S_M of the HLLC fan, or the
// velocity of the supersonic side) and the upwind side (0 = left, 1 =
// right), plus the specific internal energy of both face states.
struct FaceInfo {
  double v = 0;
  int up = 0;
  double eL = 0, eR = 0;
};

// HLLC (Toro 2009, section 10.4) with Davis wave-speed estimates
RR_HD void hllc(const double WL[NV], const double WR[NV], int d, double F[NV], FaceInfo& fi) {
  const double cL = sound_speed(WL), cR = sound_speed(WR);
  const double uL = WL[1 + d], uR = WR[1 + d];
  const double SL = rrf::min(uL - cL, uR - cR);
  const double SR = rrf::max(uL + cL, uR + cR);
  const double rL = WL[0], rR = WR[0], pL = WL[4], pR = WR[4];
  const double SM = (pR - pL + rL * uL * (SL - uL) - rR * uR * (SR - uR)) /
                    (rL * (SL - uL) - rR * (SR - uR));
  if (SL >= 0) {
    phys_flux(WL, d, F);
    fi.v = uL;
    fi.up = 0;
    return;
  }
  if (SR <= 0) {
    phys_flux(WR, d, F);
    fi.v = uR;
    fi.up = 1;
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
  fi.v = SM;
  fi.up = SM >= 0 ? 0 : 1;
}
RR_HD void hllc(const double WL[NV], const double WR[NV], int d, double F[NV]) {
  FaceInfo fi;
  hllc(WL, WR, d, F, fi);
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

// Per-cell predictor: Q = (W_half(5), dW_x(5), dW_y(5), dW_z(5)), from U
// and the selected thermal energy density (select_thermal) at the cell and
// its 6 face neighbours. U(i,j,k,m) and E(i,j,k) are the caller's
// accessors (local index space, ghosts included). Computed on the interior
// plus one ghost layer, it needs 2 ghost cells of U and E.
constexpr int NQ = 4 * NV;
template <class UA, class EA>
RR_HD void predict_cell(const UA& U, const EA& E, int i, int j, int k, double dtdx,
                        double Q[NQ]) {
  double Uc[NV], Wc[NV], Wn[3][2][NV];
  for (int m = 0; m < NV; ++m) Uc[m] = U(i, j, k, m);
  cons_to_prim(Uc, E(i, j, k), Wc);
  for (int e = 0; e < 3; ++e) {
    const int ei = e == 0, ej = e == 1, ek = e == 2;
    for (int s = 0; s < 2; ++s) {
      const int o = 2 * s - 1;
      double Un[NV];
      for (int m = 0; m < NV; ++m) Un[m] = U(i + o * ei, j + o * ej, k + o * ek, m);
      cons_to_prim(Un, E(i + o * ei, j + o * ej, k + o * ek), Wn[e][s]);
    }
  }
  double dW[3][NV];
  hancock(Wc, Wn, dtdx, Q, dW);
  for (int d = 0; d < 3; ++d)
    for (int m = 0; m < NV; ++m) Q[NV * (1 + d) + m] = dW[d][m];
}
// Total-energy form (E - KE from U), as in the spike
template <class UA>
struct TotalEnergyField {
  const UA& U;
  RR_HD double operator()(int i, int j, int k) const {
    double Uc[NV];
    for (int m = 0; m < NV; ++m) Uc[m] = U(i, j, k, m);
    return Uc[4] - kinetic_energy(Uc);
  }
};
template <class UA>
RR_HD void predict_cell(const UA& U, int i, int j, int k, double dtdx, double Q[NQ]) {
  predict_cell(U, TotalEnergyField<UA>{U}, i, j, k, dtdx, Q);
}

// HLLC flux through the face between cell L (below) and cell R (above) in
// direction d, from their predictor states.
RR_HD void face_flux(const double QL[NQ], const double QR[NQ], int d, double F[NV], FaceInfo& fi) {
  double WL[NV], WR[NV];
  for (int m = 0; m < NV; ++m) {
    WL[m] = QL[m] + 0.5 * QL[NV * (1 + d) + m];
    WR[m] = QR[m] - 0.5 * QR[NV * (1 + d) + m];
  }
  WL[0] = rrf::max(WL[0], kSmallRho);
  WL[4] = rrf::max(WL[4], kSmallP);
  WR[0] = rrf::max(WR[0], kSmallRho);
  WR[4] = rrf::max(WR[4], kSmallP);
  hllc(WL, WR, d, F, fi);
  fi.eL = WL[4] / ((kGamma - 1) * WL[0]);
  fi.eR = WR[4] / ((kGamma - 1) * WR[0]);
}
RR_HD void face_flux(const double QL[NQ], const double QR[NQ], int d, double F[NV]) {
  FaceInfo fi;
  face_flux(QL, QR, d, F, fi);
}

// Face value of a passively advected specific scalar (MUSCL-Hancock along
// the face direction only: s_half = s - dt/2 u ds/dx, then +- ds/2), for
// the cell's upper face (side = 1) or lower face (side = 0); un is the
// cell's velocity along that direction.
RR_HD double scalar_face(double sm, double s0, double sp, double un, double dtdx, int side) {
  const double ds = mc_slope(sm, s0, sp);
  return side ? s0 + 0.5 * ds * (1 - un * dtdx) : s0 - 0.5 * ds * (1 + un * dtdx);
}
// Scalar flux: the HLLC mass flux times the upwind face value
RR_HD double scalar_flux(const double F[NV], const FaceInfo& fi, double sL, double sR) {
  return F[0] * (fi.up ? sR : sL);
}

// Everything the update needs from one face: the HLLC flux F, the fluxes
// Fs of the NS scalars (the internal energy from the upwind face state's
// specific internal energy, the species from their specific abundances
// reconstructed along d) and the interface velocity v (the work term).
struct FaceData {
  double F[NV];
  double Fs[NS];
  double v;
};
// The face between cell a = (ia, ja, ka) and b = a + e_d, direction d, from
// the predictor states and the state. Accessors (local index space, ghosts
// included): Q(i, j, k) returns the NQ predictor values of a cell (const
// double*), U(i, j, k, m) the conserved state, S(i, j, k, s) the scalars.
// Reads Q of a and b, U and S of a - e_d .. b + e_d. own_a, own_b (e_self
// and e_own of the cell): the internal energy's face value on that side is
// its own specific value e/rho reconstructed along d, not that of the
// face state. The one implementation of the serial driver
// (src/hydro_unsplit.cpp) and the AMReX level (src/amrex/hydro_level.cpp).
template <class QA, class UA, class SA>
RR_HD void face_all(const QA& Q, const UA& U, const SA& S, int ia, int ja, int ka, int d,
                    double dtdx, FaceData& f, bool own_a = false, bool own_b = false) {
  const int ei = d == 0, ej = d == 1, ek = d == 2;
  const int ib = ia + ei, jb = ja + ej, kb = ka + ek;
  const double* Qa = Q(ia, ja, ka);
  const double* Qb = Q(ib, jb, kb);
  FaceInfo fi;
  face_flux(Qa, Qb, d, f.F, fi);
  f.v = fi.v;
  // species (and e with own_*): specific values s/rho reconstructed along d
  auto spec = [&](int i, int j, int k, int s) { return S(i, j, k, s) / rrf::max(U(i, j, k, 0), kSmallRho); };
  const double ua = Qa[1 + d], ub = Qb[1 + d];
  // internal energy: the specific internal energy of the upwind face state,
  // or with own_* the cell's own reconstructed e/rho
  const double eL = own_a ? scalar_face(spec(ia - ei, ja - ej, ka - ek, 0), spec(ia, ja, ka, 0),
                                        spec(ib, jb, kb, 0), ua, dtdx, 1)
                          : fi.eL;
  const double eR = own_b ? scalar_face(spec(ia, ja, ka, 0), spec(ib, jb, kb, 0),
                                        spec(ib + ei, jb + ej, kb + ek, 0), ub, dtdx, 0)
                          : fi.eR;
  f.Fs[0] = scalar_flux(f.F, fi, eL, eR);
  for (int s = 1; s < NS; ++s) {
    const double sL = scalar_face(spec(ia - ei, ja - ej, ka - ek, s), spec(ia, ja, ka, s), spec(ib, jb, kb, s),
                                  ua, dtdx, 1);
    const double sR = scalar_face(spec(ia, ja, ka, s), spec(ib, jb, kb, s), spec(ib + ei, jb + ej, kb + ek, s),
                                  ub, dtdx, 0);
    f.Fs[s] = scalar_flux(f.F, fi, sL, sR);
  }
}

// Maximum of a cell quantity E(i, j, k) over the cell and its face
// neighbours at distances 1..radius along the axes, in the fixed order
// (-x, +x, -y, +y, -z, +z) per distance (the dual-energy reference energy
// of the selection and of finalize_cell); m1 receives the radius-1 maximum.
template <class EA>
RR_HD double face_neighbour_max(const EA& E, int i, int j, int k, int radius, double& m1) {
  double m = E(i, j, k);
  m1 = m;
  for (int r = 1; r <= radius; ++r) {
    m = rrf::max(m, E(i - r, j, k));
    m = rrf::max(m, E(i + r, j, k));
    m = rrf::max(m, E(i, j - r, k));
    m = rrf::max(m, E(i, j + r, k));
    m = rrf::max(m, E(i, j, k - r));
    m = rrf::max(m, E(i, j, k + r));
    if (r == 1) m1 = m;
  }
  return m;
}

// The update of one cell from its six faces (x0 = lower x face, x1 = upper,
// ...) into Un and Sn, the serial driver's and the level's one
// implementation: update_cell for U, the internal energy with the work term
// p_half sum_d (v*_+ - v*_-) and the species without (update_scalar).
RR_HD void update_from_faces(const double U[NV], const double S[NS], const FaceData& x0, const FaceData& x1,
                             const FaceData& y0, const FaceData& y1, const FaceData& z0, const FaceData& z1,
                             double p_half, double dtdx, double Un[NV], double Sn[NS]);

// Conservative update of one cell from its 6 face fluxes:
// F[d][0] = lower face, F[d][1] = upper face.
RR_HD void update_cell(double U[NV], const double* Fx0, const double* Fx1, const double* Fy0,
                       const double* Fy1, const double* Fz0, const double* Fz1, double dtdx) {
  for (int c = 0; c < NV; ++c)
    U[c] = U[c] - dtdx * (((Fx1[c] - Fx0[c]) + (Fy1[c] - Fy0[c])) + (Fz1[c] - Fz0[c]));
}
// One scalar: s - dt/dx sum_d (F_+ - F_-), minus the work term w (the
// internal energy: w = p_half * sum_d (v*_+ - v*_-); 0 for the species).
RR_HD double update_scalar(double s, double Fx0, double Fx1, double Fy0, double Fy1, double Fz0,
                           double Fz1, double w, double dtdx) {
  return s - dtdx * ((((Fx1 - Fx0) + (Fy1 - Fy0)) + (Fz1 - Fz0)) + w);
}

RR_HD void update_from_faces(const double U[NV], const double S[NS], const FaceData& x0,
                                    const FaceData& x1, const FaceData& y0, const FaceData& y1,
                                    const FaceData& z0, const FaceData& z1, double p_half, double dtdx,
                                    double Un[NV], double Sn[NS]) {
  for (int m = 0; m < NV; ++m) Un[m] = U[m];
  update_cell(Un, x0.F, x1.F, y0.F, y1.F, z0.F, z1.F, dtdx);
  const double divv = ((x1.v - x0.v) + (y1.v - y0.v)) + (z1.v - z0.v);
  Sn[0] = update_scalar(S[0], x0.Fs[0], x1.Fs[0], y0.Fs[0], y1.Fs[0], z0.Fs[0], z1.Fs[0], p_half * divv, dtdx);
  for (int s = 1; s < NS; ++s)
    Sn[s] = update_scalar(S[s], x0.Fs[s], x1.Fs[s], y0.Fs[s], y1.Fs[s], z0.Fs[s], z1.Fs[s], 0.0, dtdx);
}

// End of a substep for one cell: density floor, the temperature floor
// (Tmin = thermal energy per unit density, applied to both energies), the
// eta2 synchronisation e := E - KE (Enb_max = max of E over the cell and
// its 6 face neighbours, after the update) and the reset of E to KE + e in
// cells whose pressure comes from the internal energy, when that changes E
// by at most eta1 E (so that the stored total energy carries the
// dual-energy temperature wherever that is cheap; docs/hydro_unsplit.md
// section 2). Counters: floor_rho, floor_T, synced, from_e, reset.
struct CellFlags {
  bool floor_rho = false, floor_T = false, synced = false, from_e = false, reset = false;
};
RR_HD CellFlags finalize_cell(double U[NV], double& eD, double Enb_max, double Enb1_max,
                              const DualEnergy& de, double rho_min, double Tmin) {
  CellFlags f;
  if (!(U[0] >= rho_min)) {
    U[0] = rho_min;
    f.floor_rho = true;
  }
  const double ke = kinetic_energy(U);
  const double emin = U[0] * Tmin;
  double eT = U[4] - ke;
  if (!(eT >= emin)) {
    // E - KE below the floor (in a cold supersonic cell this is the
    // truncation error of E, discarded below when the pressure comes
    // from e)
    eT = emin;
    U[4] = ke + emin;
    f.floor_T = true;
  }
  if (de.eta1 > 0) {
    if (eT > de.eta2 * (de.eta2_local ? U[4] : Enb_max)) {
      eD = eT;
      f.synced = true;
    }
    if (!(eD >= emin)) {
      eD = emin;
      f.floor_T = true;
    }
    if (!(eT > de.eta1 * (de.eta1_local ? U[4] : Enb_max))) {
      // pressure from e (the floor of E - KE above does not count, E - KE
      // is not used in this cell). The stored E is reset to KE + e when
      // the discarded difference is at most eta1 of the cell's own E, so
      // the reset changes the total energy by at most eta1 E per cell and
      // substep (a stagnation cell, whose E - KE is polluted by its
      // neighbours' kinetic-energy fluxes, and a pre-shock cell carrying
      // the numerical precursor of a strong shock keep their E).
      f.from_e = true;
      if (eT < emin) f.floor_T = false;
      bool reset = true;
      if (de.reset_rule == 1) reset = std::fabs(eT - eD) <= de.eta1 * U[4];
      else if (de.reset_rule == 2) reset = std::fabs(eT - eD) <= de.eta1 * U[4] || eT > 0.5 * U[4];
      else if (de.reset_rule == 3) reset = !(eT > de.eta1 * Enb1_max);
      else if (de.reset_rule == 4) reset = false;
      if (reset) {
        U[4] = ke + eD;
        f.reset = true;
      }
    }
  } else {
    eD = eT;
  }
  return f;
}

// Signal speed of a cell for the CFL condition: sum over the directions of
// |v_d| + c (RAMSES-style; the unsplit MUSCL-Hancock scheme without
// transverse corrections is stable for dt * sum_d(|v_d| + c)/dx < 1).
// e is the selected thermal energy density.
RR_HD double signal_sum(const double U[NV], double e) {
  double W[NV];
  cons_to_prim(U, e, W);
  const double c = sound_speed(W);
  return ((std::fabs(W[1]) + std::fabs(W[2])) + std::fabs(W[3])) + 3 * c;
}
// Maximum signal speed of a cell (max over directions), as in the spike
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
