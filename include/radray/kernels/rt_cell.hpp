// Per-cell RT kernels of rt_tools.f90 / rt_rates.f90, written once for
// every host: the gas -> RT grid average (prep_rt/calc_rt_i), the per-
// substep reset and optical depth (calc_rt_1), radiative_physics/photo_cell
// (the ionisation balance, photo-heating, cooling and the transmitted
// fraction fgam) and the final optical depth (post_rt/calc_rt_f).
//
// Provenance: spike/common/rt_cell.hpp of the framework-spike branch (a
// line-by-line copy of src/rt_tools.cpp and src/rt_rates.cpp for
// N_rt_gas = 1), extended here with everything the serial port has and the
// spike had left out: the legacy_rt switch of calc_rt_1 (B3, B4), the
// legacy_physics selection of the rates and of the cooling sum, the case A
// recombination (recomb_case) and the ground-state recombination count of
// the recombination radiation, the dual-energy field of hydro_scheme =
// unsplit, the rt_stat sums (nion, ncol, nrec, nrab, gamma, lambda, dtau)
// and the temperature-index error of the legacy rates (an error code: no
// exceptions here, the host throws). This copy is the maintained one; the
// serial port (src/rt_tools.cpp, src/rt_rates.cpp) calls these functions
// through CellRef views of its structure-of-arrays grid, the AMReX backend
// through CellData (docs/stage_d.md, milestone 1).
//
// Header-only, host/device (rr_hd.hpp). Precision: RT grid fields are
// float (real(4)), arithmetic is double (real(8)), as in the Fortran.
#pragma once

#include "radray/constants.hpp"
#include "radray/cosmology.hpp"
#include "radray/kernels/consts_hd.hpp"
#include "radray/kernels/rr_hd.hpp"

namespace rrk {

// One RT cell (rt_grid_type fields) plus its gas cell u_gas(5), w_gas(4)
// and the dual-energy field e of hydro_scheme = unsplit (used only when the
// host says has_e; N_rt_gas = 1), the layout of the AMReX backend's cell
// FabArray. 46 packed floats.
struct CellData {
  float rhog, fHI, fHeI, fHeII, temp, Luv, sfr;
  float dx[3], gamma[3], lambda[3];
  float nu[5], ngam[5], fgam[5], tau[5];
  float u[5], w[4], e;
};

// Pointers to the rt_grid fields of one cell, wherever they live: the
// structure-of-arrays RTGrid of the serial port (one pointer per field into
// its vectors) or a CellData.
struct CellRef {
  float *rhog = nullptr, *fHI = nullptr, *fHeI = nullptr, *fHeII = nullptr, *temp = nullptr;
  float *Luv = nullptr, *sfr = nullptr;
  float *dx = nullptr, *gamma = nullptr, *lambda = nullptr;  // 3 each
  float *nu = nullptr, *ngam = nullptr, *fgam = nullptr, *tau = nullptr;  // 5 each
};
RR_HD CellRef cell_ref(CellData& c) {
  CellRef r;
  r.rhog = &c.rhog;
  r.fHI = &c.fHI;
  r.fHeI = &c.fHeI;
  r.fHeII = &c.fHeII;
  r.temp = &c.temp;
  r.Luv = &c.Luv;
  r.sfr = &c.sfr;
  r.dx = c.dx;
  r.gamma = c.gamma;
  r.lambda = c.lambda;
  r.nu = c.nu;
  r.ngam = c.ngam;
  r.fgam = c.fgam;
  r.tau = c.tau;
  return r;
}

// rt_rates.f90: photoion_sigma (Verner et al. 1996), sigma[species][freq]
// in cm^2. Exponents and the literals are real(4) promoted to double;
// 2.136**2 and 2.039**2 are real(4) products.
RR_HD void photoion_sigma(const double nu[5], double sig[3][5]) {
  using namespace radray;
  using namespace radray::pc;
  auto s = [](float v) { return static_cast<double>(v); };
  for (int i = 0; i < 5; ++i) {
    const double fq = nu[i];
    double x, y;
    // HI
    if (fq < HI_eV) {
      sig[0][i] = 0;
    } else {
      x = fq / s(0.4298f);
      y = x;
      sig[0][i] = 5.475e-14 * rrf::pow2(x - 1) / std::pow(y, s(4.0185f)) /
                  std::pow(1 + std::sqrt(y / s(32.88f)), s(2.963f));
    }
    // HeI
    if (fq < HeI_eV) {
      sig[1][i] = 0;
    } else {
      x = fq / s(13.61f) - s(0.4434f);
      y = std::sqrt(rrf::pow2(x) + s(2.136f * 2.136f));
      sig[1][i] = 9.492e-16 * (rrf::pow2(x - 1) + s(2.039f * 2.039f)) /
                  std::pow(y, s(3.9060f)) / std::pow(1 + std::sqrt(y / s(1.469f)), s(3.188f));
    }
    // HeII
    if (fq < HeII_eV) {
      sig[2][i] = 0;
    } else {
      x = fq / s(1.720f);
      y = x;
      sig[2][i] = 1.369e-14 * rrf::pow2(x - 1) / std::pow(y, s(4.0185f)) /
                  std::pow(1 + std::sqrt(y / s(32.88f)), s(2.963f));
    }
  }
}

// Atomic-rate tables (AtomicRates::coef_cool_ etc., Fortran layout:
// coef_cool(16,Ntemp), alpha_rec(6,Ntemp), gamma_col(3,Ntemp))
struct RatesView {
  const double* cool = nullptr;   // (16, Ntemp)
  const double* alpha = nullptr;  // (6, Ntemp)
  const double* gcol = nullptr;   // (3, Ntemp)
  bool legacy_physics = true;     // nearest node (legacy) or interpolated (P9)
};
// All rates at one temperature. Index n-1 holds Fortran slot n.
struct RateSetP {
  double alpha[6], gcol[3], cool[16];
};
// Table entries at node iT (1..Ntemp): the legacy nearest-node lookup.
RR_HD RateSetP rates_node(const RatesView& R, int iT) {
  RateSetP r;
  for (int n = 0; n < 6; ++n) r.alpha[n] = R.alpha[n + 6 * (iT - 1)];
  for (int n = 0; n < 3; ++n) r.gcol[n] = R.gcol[n + 3 * (iT - 1)];
  for (int n = 0; n < 16; ++n) r.cool[n] = R.cool[n + 16 * (iT - 1)];
  return r;
}
// P9 (legacy_physics = 0; harness PATCHES.md 10, rt_rates.f90 rate_lerp).
// Linear interpolation of log R in log T: R = R_i exp(w log(R_{i+1}/R_i)).
// If either node is 0 (below the lambda/2 > 100 cut-offs), linear in R.
RR_HD double rate_lerp(double a, double b, double w) {
  if (a > 0 && b > 0) return a * std::exp(w * std::log(b / a));
  return a + w * (b - a);
}
// Between nodes i and i+1 (1 <= i < Ntemp) with weight w in [0,1)
RR_HD RateSetP rates_interp(const RatesView& R, int i, double w) {
  RateSetP r;
  for (int n = 0; n < 6; ++n) r.alpha[n] = rate_lerp(R.alpha[n + 6 * (i - 1)], R.alpha[n + 6 * i], w);
  for (int n = 0; n < 3; ++n) r.gcol[n] = rate_lerp(R.gcol[n + 3 * (i - 1)], R.gcol[n + 3 * i], w);
  for (int n = 0; n < 16; ++n) r.cool[n] = rate_lerp(R.cool[n + 16 * (i - 1)], R.cool[n + 16 * i], w);
  return r;
}
// Harness rt_rates.f90 rates_interp: u = (log10 T - lgT_min)/lgT_del,
// i = 1 + int(u), w = u - (i-1). Outside the table (u <= 0 or
// u >= Ntemp-1, including NaN) the edge node is used, no extrapolation.
RR_HD RateSetP rates_interp(const RatesView& R, double tmp) {
  using namespace radray;
  const double u = (std::log10(tmp) - cp::lgT_min) / cp::lgT_del;
  if (u > 0 && u < cp::Ntemp - 1) {
    const int i = 1 + static_cast<int>(u);
    const double w = u - (i - 1);
    return rates_interp(R, i, w);
  }
  return rates_node(R, u > 0 ? cp::Ntemp : 1);
}

// Average of the Eulerian cells inside an RT cell + derived gas state.
// Shared by calc_rt_i (prep_rt), calc_rt_f (post_rt) and the end of
// photo_cell, which contain the same code.
struct GasStateP {
  double rhog, nH, nHe, nHI, nHeI, nHeII, nHII, nHeIII, ne, neff, tmp;
};
// es: the averaged internal-energy density of hydro_scheme = unsplit (the
// dual-energy temperature, docs/hydro_unsplit.md section 4), or < 0 with
// the legacy scheme (the arithmetic below is then the original's).
RR_HD GasStateP gas_state_from_avg(const double us_in[5], const double ws[4],
                                   const radray::Units& u, double es = -1) {
  using namespace radray;
  using namespace radray::pc;
  using namespace radray::cp;
  GasStateP g;
  const double us1 = us_in[0];
  g.rhog = us1 * u.rho_unit;
  g.nH = g.rhog * XHy / mH_cgs;
  g.nHe = g.rhog * YHe / mHe_cgs;
  g.nHI = rrf::min(ws[0] * u.H_unit, g.nH);
  g.nHeI = rrf::min(ws[1] * u.He_unit, g.nHe);
  g.nHeII = rrf::min(ws[2] * u.He_unit, g.nHe - g.nHeI);
  g.nHII = g.nH - g.nHI;
  g.nHeIII = g.nHe - g.nHeI - g.nHeII;
  g.ne = g.nHII + g.nHeII + 2 * g.nHeIII;
  g.neff = g.nH + g.nHe + g.ne;
  // eth = (us(5)-sum(us(2:4)**2)/us(1)/2)*eden_unit   (real(8))
  const double ksum = ((0.0 + us_in[1] * us_in[1]) + us_in[2] * us_in[2]) + us_in[3] * us_in[3];
  double eth = (us_in[4] - ksum / us1 / 2) * u.eden_unit;
  if (es >= 0) eth = es * u.eden_unit;
  g.tmp = rrf::max((gamma_ideal - 1) * eth / k_cgs / g.neff, u.temp_min);
  return g;
}
RR_HD void store_gas_state(const CellRef& c, const GasStateP& g, const radray::Units& u) {
  *c.rhog = static_cast<float>(g.rhog / u.rhob_avg);
  *c.fHI = static_cast<float>(rrf::min(rrf::max(g.nHI / g.nH, 0.0), 1.0));
  *c.fHeI = static_cast<float>(rrf::min(rrf::max(g.nHeI / g.nHe, 0.0), 1.0));
  *c.fHeII = static_cast<float>(rrf::min(rrf::max(g.nHeII / g.nHe, 0.0), 1.0));
  *c.temp = static_cast<float>(g.tmp);
}
// post_rt / calc_rt_f: the optical depth from the gas state
RR_HD void tau_from_state(float* tau, const GasStateP& g, const double sig[3][5], double len_rt) {
  for (int i = 0; i < 5; ++i) {
    const double sm = ((0.0 + g.nHI * sig[0][i]) + g.nHeI * sig[1][i]) + g.nHeII * sig[2][i];
    tau[i] = static_cast<float>(sm * len_rt);
  }
}

// average_gas for N_rt_gas = 1 (normalised by N_rt_gas_c = 1, clamped)
RR_HD void average_gas_1(const float* u5, const float* w4, const radray::Units& u, double us[5],
                         double ws[4]) {
  for (int m = 0; m < 5; ++m) us[m] = (0.0 + static_cast<double>(u5[m])) / 1;
  for (int m = 0; m < 4; ++m) ws[m] = (0.0 + static_cast<double>(w4[m])) / 1;
  us[0] = rrf::max(us[0], u.Dgas_min);
}

// prep_rt: calc_rt_i (N_rt_gas = 1; has_e: the dual-energy temperature)
RR_HD void prep_cell(CellData& c, const radray::Units& u, bool has_e = false) {
  double us[5], ws[4];
  average_gas_1(c.u, c.w, u, us, ws);
  const CellRef r = cell_ref(c);
  store_gas_state(r, gas_state_from_avg(us, ws, u, has_e ? static_cast<double>(c.e) : -1.0), u);
}

// rt_step_1: calc_rt_1. Resets the radiation fields and updates the
// optical depth of (partially) ionised cells. legacy = true (legacy_rt = 1):
// rt_grid%dx is NOT reset (B3) and tau is only updated where fHI or
// fHeI < 0.99 (B4), as in the original; false: both fixed.
// sig = photoion_sigma(freq_init).
RR_HD void calc_rt_1_cell(const CellRef& c, const radray::Units& u, const double sig[3][5],
                          bool legacy) {
  using namespace radray::cp;
  using namespace radray::pc;
  for (int n = 0; n < 5; ++n) {
    c.nu[n] = 0;
    c.ngam[n] = 0;
    c.fgam[n] = 0;
  }
  for (int d = 0; d < 3; ++d) {
    c.gamma[d] = 0;
    c.lambda[d] = 0;
  }
  if (!legacy)
    for (int d = 0; d < 3; ++d) c.dx[d] = 0;
  // real(4) comparisons with the real(4) literal 0.99
  if (*c.fHI < 0.99f || *c.fHeI < 0.99f || !legacy) {
    const double rhog = static_cast<double>(*c.rhog) * u.rhob_avg;
    const double nH = rhog * XHy / mH_cgs;
    const double nHe = rhog * YHe / mHe_cgs;
    const double nHI = nH * static_cast<double>(*c.fHI);
    const double nHeI = nHe * static_cast<double>(*c.fHeI);
    const double nHeII = nHe * static_cast<double>(*c.fHeII);
    for (int i = 0; i < 5; ++i) {
      const double sm = ((0.0 + nHI * sig[0][i]) + nHeI * sig[1][i]) + nHeII * sig[2][i];
      c.tau[i] = static_cast<float>(sm * u.len_rt);
    }
  }
}

// post_rt: calc_rt_f (gas -> RT grid, then tau from the gas), N_rt_gas = 1
RR_HD void post_cell(CellData& c, const radray::Units& u, const double sig[3][5],
                     bool has_e = false) {
  double us[5], ws[4];
  average_gas_1(c.u, c.w, u, us, ws);
  const GasStateP g = gas_state_from_avg(us, ws, u, has_e ? static_cast<double>(c.e) : -1.0);
  const CellRef r = cell_ref(c);
  store_gas_state(r, g, u);
  tau_from_state(c.tau, g, sig, u.len_rt);
}

// ---------------------------------------------------------------------------
// radiative_physics (rt_tools.f90), in three parts per RT cell:
//   1. cell_radiation_prep: effective frequency and direction of the
//      deposited photons; false if the cell is skipped (no photons and
//      neutral);
//   2. photo_cell for each Eulerian gas cell of the RT cell: photo_gas_cell;
//   3. photo_cell_finish: fgam, the RT grid gas state, the transmitted
//      photons.
// The host loops over the gas cells (make_index_gas_rt for N_rt_gas > 1),
// accumulates the rt_stat sums (PhotoStats) in its own order and handles
// the recombination-radiation buffer (recomb_count).
// ---------------------------------------------------------------------------

constexpr int kPhotonIterMax = 1000;  // photon_conserve: iteration cap of photo_gas_cell

struct PhysCtx {
  radray::Units u;
  RatesView rates;
  double dtsec_rt = 0;
  // rt_cells_per_substep (docs/rt_multicell.md): a ray crosses this many RT
  // cells per substep, so the photons deposited in a cell (its share of the
  // ray's path) cross it at cells*len_rt/dtsec_rt. 1 is the original.
  int cells = 1;
  // B17 fix (photon_conserve, docs/rt_multicell.md): iterate the ionisation
  // balance until the photon limit holds in every frequency group (at most
  // kPhotonIterMax iterations), not only until ne converges (at most 5, the
  // original). false is the original.
  bool photon_conserve = false;
  // Config::cool_suppress_overdensity: the overdensity above which the
  // cooling and SN heating are suppressed ((x0/x)^2, (x0/x)^1.5); 200 is
  // the original, 0 switches the suppression off.
  double cool_suppress = 200.0;
  // recomb_case = A: ia = 0 (alpha_rec(1:3), coef_cool(1:3)); B (the
  // original): 3. The host lowers it to 3 for cells on the spot
  // (recomb_rad). want_nrab: count the ground-state recombinations
  // alpha_A - alpha_B (recomb_rad = 1).
  int ia = 3;
  bool want_nrab = false;
  // xray_heating (docs/xray_heating.md, xray_background.hpp): the uniform
  // X-ray background's primary photoionisation rates (1/s) and
  // photoelectron powers (erg/s) per HI, HeI, HeII atom, set by the host
  // for the hydro step. false (the default) is the original: no X-rays.
  bool xray = false;
  double xray_gamma[3] = {0, 0, 0};
  double xray_phi[3] = {0, 0, 0};
};

// xray_heating: the fractions of a fast photoelectron's energy that go into
// heat and into HI and HeI ionisations, at x = n_HII/n_H clamped to the
// fits' range [1e-4, 1]: Shull & van Steenberg (1985) as fitted in Ricotti,
// Gnedin & Shull (2002). Coefficients from memory, to be confirmed
// (docs/xray_heating.md section 5).
RR_HD void xray_secondary(double x, double& fheat, double& fion_HI, double& fion_HeI) {
  x = x < 1e-4 ? 1e-4 : (x > 1 ? 1 : x);
  fheat = 0.9971 * (1 - std::pow(1 - std::pow(x, 0.2663), 1.3163));
  fion_HI = 0.3908 * std::pow(1 - std::pow(x, 0.4092), 1.7592);
  fion_HeI = 0.0554 * std::pow(1 - std::pow(x, 0.4614), 1.6660);
}

// rt_stat addends of one gas cell (before the /navg/Ncell normalisation of
// the host) and the photo-heating/photoionisation rates of the RT grid
struct PhotoStats {
  double nion[3] = {0, 0, 0}, ncol[3] = {0, 0, 0}, nrec[3] = {0, 0, 0}, nrab[3] = {0, 0, 0};
  double nrB[3] = {0, 0, 0};  // case-B recombinations (recomb_rad_lines)
  double gam[3] = {0, 0, 0}, lam[3] = {0, 0, 0};
  double dtau = 0;   // sT*ne*c*dt (physical c), before /Ncell
  double T_in = 0;   // temperature before heating/cooling (for recomb_count)
};

// Error codes of photo_gas_cell
constexpr int kPhotoOK = 0, kPhotoBadTemperatureIndex = 1;

// Part 1. c.ngam/c.nu/c.dx hold the deposited photon sums. Returns whether
// the cell is processed; ngam_tot gets the photon total.
RR_HD bool cell_radiation_prep(const CellRef& c, const radray::Units& u, double& ngam_tot) {
  float* gng = c.ngam;
  float* gnu = c.nu;
  float* gdx = c.dx;
  // Effective frequency
  for (int n = 0; n < 5; ++n) {
    if (static_cast<double>(gng[n]) < u.ngam_min) {
      gnu[n] = static_cast<float>(rrk::freq_init(n));
      gng[n] = 0;
    } else {
      gnu[n] = gnu[n] / gng[n];
    }
  }
  ngam_tot = static_cast<double>(rrf::sum4<5>(gng));
  // Effective unit direction vector
  if (ngam_tot < u.ngam_min) {
    gdx[0] = gdx[1] = gdx[2] = 0;
  } else {
    const float sng = rrf::sum4<5>(gng);
    double dx[3];
    for (int d = 0; d < 3; ++d) dx[d] = static_cast<double>(gdx[d] / sng);
    const double dl = std::sqrt(((0.0 + dx[0] * dx[0]) + dx[1] * dx[1]) + dx[2] * dx[2]);
    if (dl > static_cast<double>(1e-9f)) {
      for (int d = 0; d < 3; ++d) gdx[d] = static_cast<float>(dx[d] / dl);
    } else {
      gdx[0] = gdx[1] = gdx[2] = 0;
    }
  }
  return ngam_tot > u.ngam_min || *c.fHI < 0.99f;
}

// photo_source of a stellar source cell: the rt_stat ngam addends per
// species and frequency group (before /nrtcr), 0 where the group is below
// the species' threshold. The host adds them group by group, in group
// order, as the original does.
RR_HD void source_stats(float Luv, double dtsec_rt, const radray::Units& u, double ngam_sp[3][5]) {
  using namespace radray;
  using namespace radray::pc;
  using namespace radray::cp;
  for (int inu = 0; inu < 5; ++inu) {
    const double nu = rrk::freq_init(inu);
    double ngam = fResc * rrk::pop2_init(inu) * (static_cast<double>(Luv) * dtsec_rt);
    ngam = ngam / u.vol_rt / u.nH_avg;
    ngam_sp[0][inu] = nu > HI_eV ? ngam : 0;
    ngam_sp[1][inu] = nu > HeI_eV ? ngam * (u.nH_avg / u.nHe_avg) : 0;
    ngam_sp[2][inu] = nu > HeII_eV ? ngam * (u.nH_avg / u.nHe_avg) : 0;
  }
}

// Part 2: photo_cell for one Eulerian gas cell uu(5), ww(4) (and the
// dual-energy field *ee if ee != nullptr). gng[5]: the RT cell's photon
// counts (rt_grid%ngam, per <nH>); nu[5]: its effective frequencies;
// sigma_pht = photoion_sigma(nu); nrad[5]: the photons of the RT cell
// (in/out: what remains after this gas cell). The returned gam/lam
// (PhotoStats) are the rates the host adds to rt_grid%gamma/lambda (/Nc).
RR_HD int photo_gas_cell(float* uu, float* ww, float* ee, const float gng[5], const double nu[5],
                         const double sigma_pht[3][5], double nrad[5], const PhysCtx& X,
                         PhotoStats& st) {
  using namespace radray;
  using namespace radray::pc;
  using namespace radray::cp;
  const Units& u = X.u;
  const double g = gamma_ideal;
  const double nu0[3] = {HI_eV, HeI_eV, HeII_eV};
  const double dtsec_rt = X.dtsec_rt;
  const int ia = X.ia;

  // Gas densities
  uu[0] = static_cast<float>(rrf::max(static_cast<double>(uu[0]), u.Dgas_min));
  const double rhog = static_cast<double>(uu[0]) * u.rho_unit;
  const double nH = rhog * XHy / mH_cgs;
  const double nHe = rhog * YHe / mHe_cgs;
  double nHI = rrf::min(static_cast<double>(ww[0]) * u.H_unit, nH);
  double nHeI = rrf::min(static_cast<double>(ww[1]) * u.He_unit, nHe);
  double nHeII = rrf::min(static_cast<double>(ww[2]) * u.He_unit, nHe - nHeI);
  double nHII = nH - nHI;
  double nHeIII = nHe - nHeI - nHeII;
  double ne = nHII + nHeII + 2 * nHeIII;
  double neff = nH + nHe + ne;
  double nff = nHII + nHeII + 4 * nHeIII;

  // Temperature. kin = sum(u(2:4)**2)/u(1)/2 is pure real(4).
  const float kin_f = (((0.0f + uu[1] * uu[1]) + uu[2] * uu[2]) + uu[3] * uu[3]) / uu[0] / 2.0f;
  const double kin = static_cast<double>(kin_f);
  double eth = (static_cast<double>(uu[4]) - kin) * u.eden_unit;
  // hydro_scheme = unsplit (docs/hydro_unsplit.md section 4): the thermal
  // energy is the dual-energy field, which gets the same heating/cooling as
  // u(5). Null with the legacy scheme: no-op.
  if (ee) eth = static_cast<double>(*ee) * u.eden_unit;
  double tmp = (g - 1) * eth / k_cgs / neff;
  if (tmp < u.temp_min) {
    tmp = u.temp_min;
    eth = neff * k_cgs * tmp / (g - 1);
  }
  st.T_in = tmp;

  // Temperature index for RT rates (legacy: the table node; out of the
  // table is an error the host reports)
  const int iT = 1 + rrf::int4(std::round((std::log10(tmp) - lgT_min) / lgT_del));
  if (X.rates.legacy_physics && (iT < 1 || iT > Ntemp)) return kPhotoBadTemperatureIndex;
  // P9 (legacy_physics = 0): rates interpolated in log T; with
  // legacy_physics the table values at the nearest node iT, as the original.
  const RateSetP rs = X.rates.legacy_physics ? rates_node(X.rates, iT) : rates_interp(X.rates, tmp);

  double ngam[5];
  for (int n = 0; n < 5; ++n) ngam[n] = static_cast<double>(gng[n]) * u.nH_avg;

  // Maximum photoionization rate in 1/s
  double gamma_pht[3][5], gamma_max[3][5];
  for (int inu = 0; inu < 5; ++inu)
    for (int sp = 0; sp < 3; ++sp)
      gamma_pht[sp][inu] = ngam[inu] * sigma_pht[sp][inu] * (u.len_rt * X.cells / dtsec_rt);
  for (int sp = 0; sp < 3; ++sp)
    for (int inu = 0; inu < 5; ++inu) gamma_max[sp][inu] = gamma_pht[sp][inu];
  double gam_pht[3];
  auto sum_dim2 = [&](double out[3]) {
    for (int sp = 0; sp < 3; ++sp) {
      double t = 0;
      for (int inu = 0; inu < 5; ++inu) t += gamma_pht[sp][inu];
      out[sp] = t;
    }
  };
  sum_dim2(gam_pht);

  // xray_heating: the X-ray background's ionisation rates per atom, primary
  // plus the secondary ionisations of HI and HeI by the photoelectrons, with
  // the densities and ionised fraction at the start of the substep. Not
  // limited by the RT's photons (a uniform background); zero without X-rays,
  // which leaves the sums below unchanged bit for bit.
  double gx[3] = {0, 0, 0};
  double fheat_x = 0;
  if (X.xray) {
    double fion_HI = 0, fion_HeI = 0;
    xray_secondary(nHII / nH, fheat_x, fion_HI, fion_HeI);
    const double edep = (X.xray_phi[0] * nHI + X.xray_phi[1] * nHeI) + X.xray_phi[2] * nHeII;
    gx[0] = X.xray_gamma[0] + (nHI > 0 ? fion_HI * edep / (HI_eV * eV2erg * nHI) : 0.0);
    gx[1] = X.xray_gamma[1] + (nHeI > 0 ? fion_HeI * edep / (HeI_eV * eV2erg * nHeI) : 0.0);
    gx[2] = X.xray_gamma[2];
  }

  // Implicit solution for ionisation balance. B17: the photon limit below
  // rescales gamma_pht after the densities were computed with the previous
  // rates, and the loop ends when ne has converged or after 5 iterations.
  // In a photon-starved cell with gamma*dt >> 1 every iteration ionises the
  // cell almost fully, so ne "converges" while the rates still exceed the
  // photons: the cell is ionised by photons it was not given (a few per cent
  // of the ionisations in dense fronts; the error grows with gamma_max*dt,
  // as k^2 with rt_cells_per_substep = k). photon_conserve = 1 iterates
  // until no group is limited by more than 1e-3 as well.
  int n = 0;
  double ne_prev = static_cast<double>(3.40282346638528859812e+38F);  // huge(0.)
  double nH1 = 0, nHe1 = 0, nHe2 = 0;
  bool limited = false;
  const int nit_max = X.photon_conserve ? kPhotonIterMax : 5;
  while ((std::fabs(ne / ne_prev - 1) > static_cast<double>(1e-3f) || limited) && n < nit_max) {
    n = n + 1;
    limited = false;
    ne_prev = ne;
    const double AH2 = rs.alpha[ia + 0] * ne * dtsec_rt;
    const double AHe2 = rs.alpha[ia + 1] * ne * dtsec_rt;
    const double AHe3 = rs.alpha[ia + 2] * ne * dtsec_rt;
    const double BH1 = (gam_pht[0] + gx[0] + rs.gcol[0] * ne) * dtsec_rt;
    const double BHe1 = (gam_pht[1] + gx[1] + rs.gcol[1] * ne) * dtsec_rt;
    const double BHe2 = (gam_pht[2] + gx[2] + rs.gcol[2] * ne) * dtsec_rt;
    const double sub1 = 1 + AHe3 + AHe2 + BHe2;
    const double sub2 = AHe2;
    const double sub3 = BHe1 - AHe3;
    const double sub4 = 1 + BHe1;
    const double sub5 = AHe3 * nHe;
    const double sub6 = sub1 * sub4 - sub2 * sub3;
    nH1 = (nHI + AH2 * nH) / (1 + AH2 + BH1);
    nHe1 = (nHeI * sub1 + (nHeII + sub5) * sub2) / sub6;
    nHe2 = (nHeI * sub3 + (nHeII + sub5) * sub4) / sub6;
    ne = (nH - nH1) + nHe2 + 2 * (nHe - nHe1 - nHe2);
    // Photoionisations can not exceed the number of photons
    for (int inu = 0; inu < 5; ++inu) {
      const double num1 = (((0.0 + gamma_pht[0][inu] * nH1) + gamma_pht[1][inu] * nHe1) +
                           gamma_pht[2][inu] * nHe2) *
                          dtsec_rt * u.vol_gas;
      if (num1 > u.ngam_min) {
        const double x = nrad[inu] / num1;
        if (X.photon_conserve && x < 1 - 1e-3) limited = true;
        for (int sp = 0; sp < 3; ++sp)
          gamma_pht[sp][inu] = rrf::min(x * gamma_pht[sp][inu], gamma_max[sp][inu]);
      }
    }
    sum_dim2(gam_pht);
  }
  // (n == 0 cannot happen: the first test is ~1 > 1e-3)

  // Update H, He number densities
  nHI = nH1;
  nHII = nH - nHI;
  nHeI = nHe1;
  nHeII = nHe2;
  nHeIII = nHe - nHeI - nHeII;
  ne = nHII + nHeII + 2 * nHeIII;
  neff = nH + nHe + ne;
  nff = nHII + nHeII + 4 * nHeIII;

  // Count photoionisations, photon conservation
  double nion[3] = {0, 0, 0};
  for (int inu = 0; inu < 5; ++inu) {
    const double d = rhog / u.rhob_avg;
    const double x = rrf::min(std::exp(-(d - rho_source) / rho_source), 1.0);
    const double num3[3] = {gamma_pht[0][inu] * nHI * dtsec_rt, gamma_pht[1][inu] * nHeI * dtsec_rt,
                            gamma_pht[2][inu] * nHeII * dtsec_rt};
    for (int sp = 0; sp < 3; ++sp) nion[sp] = nion[sp] + num3[sp];
    const double snum3 = ((0.0 + num3[0]) + num3[1]) + num3[2];
    nrad[inu] = rrf::max(nrad[inu] - x * snum3 * u.vol_gas, 0.0);
  }

  // Collisional ionisations and recombinations
  const double nsp[3] = {nHI, nHeI, nHeII};
  const double nsp_rec[3] = {nHII, nHeII, nHeIII};
  for (int sp = 0; sp < 3; ++sp) {
    st.nion[sp] = nion[sp];
    st.ncol[sp] = rs.gcol[sp] * nsp[sp] * ne * dtsec_rt;
    st.nrec[sp] = rs.alpha[sp + ia] * nsp_rec[sp] * ne * dtsec_rt;
    // ground-state recombinations: case A - case B (recomb_rad = 1)
    st.nrab[sp] = X.want_nrab ? (rs.alpha[sp] - rs.alpha[sp + 3]) * nsp_rec[sp] * ne * dtsec_rt : 0;
    // and, with recomb_rad_lines, the case-B ones (helium lines)
    st.nrB[sp] = X.want_nrab ? rs.alpha[sp + 3] * nsp_rec[sp] * ne * dtsec_rt : 0;
  }

  // Photoheating
  double lam_pht[3];
  for (int sp = 0; sp < 3; ++sp) {
    double t = 0;
    for (int inu = 0; inu < 5; ++inu) t += gamma_pht[sp][inu] * (nu[inu] - nu0[sp]) * eV2erg;
    lam_pht[sp] = t;
  }
  double heat_rate = ((0.0 + lam_pht[0] * nHI) + lam_pht[1] * nHeI) + lam_pht[2] * nHeII;
  // xray_heating: the heat share of the photoelectron power absorbed by the
  // updated densities (as the UV photo-heating)
  if (X.xray)
    heat_rate = heat_rate +
                fheat_x * ((X.xray_phi[0] * nHI + X.xray_phi[1] * nHeI) + X.xray_phi[2] * nHeII);

  // SN heating
  double sn_rate;
  if (tmp > 1e4 && rhog / u.rhob_avg > rho_source) {
    const double d = rhog / u.rhob_avg;
    const double x = rrf::max(1 - std::exp(-(d - rho_source) / rho_source), 0.0);
    const double tdyn = std::sqrt(3 * pi / (32 * G_cgs * rhog));
    sn_rate = x * eSN * eSF * rhog / tdyn * (c_cgs * c_cgs);
  } else {
    sn_rate = 0;
  }

  // Cooling
  const double Tcmb = u.Tcmb;
  double cool_rate;
  if (X.rates.legacy_physics) {
    cool_rate = rs.cool[ia + 0] * ne * nHII + rs.cool[ia + 1] * ne * nHeII +
                rs.cool[ia + 2] * ne * nHeIII + rs.cool[6] * ne * nHI + rs.cool[7] * ne * nHeI +
                rs.cool[8] * ne * nHeII + rs.cool[9] * ne * nHI + rs.cool[10] * ne * nHeI +
                rs.cool[11] * ne * nHeII + rs.cool[12] * ne * nff +
                rs.cool[13] * ne * (tmp - Tcmb) * rrf::pow4(Tcmb);
  } else {
    // P2/P3 (legacy_physics = 0): HeI collisional excitation (Cen 1992)
    // scales as ne^2*nHeII; slot 15 is HeII dielectronic recombination
    // cooling (Hui & Gnedin 1997), times ne*nHeII. Same term order as the
    // Fortran, left to right.
    cool_rate = rs.cool[ia + 0] * ne * nHII + rs.cool[ia + 1] * ne * nHeII +
                rs.cool[ia + 2] * ne * nHeIII + rs.cool[6] * ne * nHI + rs.cool[7] * ne * nHeI +
                rs.cool[8] * ne * nHeII + rs.cool[9] * ne * nHI + rs.cool[10] * ne * ne * nHeII +
                rs.cool[11] * ne * nHeII + rs.cool[12] * ne * nff +
                rs.cool[13] * ne * (tmp - Tcmb) * rrf::pow4(Tcmb) + rs.cool[14] * ne * nHeII +
                rs.cool[15] * ne * nHeI;
  }

  // Prevent overcooling/overheating (above 200 x the mean density in the
  // original; cool_suppress_overdensity, 0 = off)
  {
    const double x = rhog / u.rhob_avg;
    const double x0 = X.cool_suppress;
    if (x0 > 0 && x > x0) {
      cool_rate = cool_rate * rrf::pow2(x0 / x);
      sn_rate = sn_rate * std::pow(x0 / x, 1.5);
    }
  }
  const double net_rate = heat_rate + sn_rate - cool_rate;
  if (net_rate > 0) {
    eth = eth + net_rate * dtsec_rt;
  } else if (net_rate < 0) {
    eth = eth / (1 + std::log(1 - net_rate * dtsec_rt / eth));
  }
  tmp = (g - 1) * eth / k_cgs / neff;
  if (tmp < u.temp_min) {
    tmp = u.temp_min;
    eth = neff * k_cgs * tmp / (g - 1);
  }

  // Update u_gas and w_gas
  uu[4] = static_cast<float>(kin + eth / u.eden_unit);
  if (ee) *ee = static_cast<float>(eth / u.eden_unit);
  ww[0] = static_cast<float>(rrf::min(rrf::max(nHI, 0.0), nH) / u.H_unit);
  ww[1] = static_cast<float>(rrf::min(rrf::max(nHeI, 0.0), nHe) / u.He_unit);
  ww[2] = static_cast<float>(rrf::min(rrf::max(nHeII, 0.0), nHe) / u.He_unit);

  for (int sp = 0; sp < 3; ++sp) {
    st.gam[sp] = gam_pht[sp];
    st.lam[sp] = lam_pht[sp];
  }
  // Thomson optical depth: physical c (and physical dt), also with
  // c_reduced < 1
  st.dtau = sT_cgs * ne * (c_cgs * dtsec_rt);
  return kPhotoOK;
}

// Part 3: the end of photo_cell and of radiative_physics. us/ws: the sums
// of the updated gas cells (divided by Nc = N_rt_gas_c here); es: the sum
// of their dual-energy field, or < 0 with the legacy scheme; nrad: the
// photons that remain.
RR_HD void photo_cell_finish(const CellRef& c, const double nrad[5], double us[5], double ws[4],
                             double es, double Nc, const radray::Units& u) {
  float* gng = c.ngam;
  float* gfg = c.fgam;
  // Remaining photons -> fgam (numerator)
  for (int n = 0; n < 5; ++n)
    gfg[n] = static_cast<float>(static_cast<double>(gfg[n]) + nrad[n] / u.vol_rt / u.nH_avg);
  // Update gas on RT grid (no Dgas_min clamp here, as in the original)
  for (int m = 0; m < 5; ++m) us[m] = us[m] / Nc;
  for (int m = 0; m < 4; ++m) ws[m] = ws[m] / Nc;
  store_gas_state(c, gas_state_from_avg(us, ws, u, es < 0 ? -1.0 : es / Nc), u);
  // Update radiation on RT grid
  for (int m = 0; m < 5; ++m) {
    if (static_cast<double>(gng[m]) < u.ngam_min) {
      gfg[m] = 0;
      gng[m] = 0;
    } else {
      gfg[m] = gfg[m] / gng[m];
      gng[m] = gng[m] * gfg[m];
    }
  }
}

// radiative_physics for one RT cell with N_rt_gas = 1 (the AMReX backend):
// parts 1-3 on a CellData. st (optional) receives the rt_stat addends;
// processed (optional) whether the cell passed part 1 and ngam_tot its
// photon total; has_e: the dual-energy field c.e is the thermal energy.
// Returns the error code of photo_gas_cell (0 if the cell was skipped).
RR_HD int physics_cell(CellData& c, const PhysCtx& X, PhotoStats* st = nullptr,
                       bool has_e = false, bool* processed = nullptr, double* ngam_tot_out = nullptr) {
  const CellRef r = cell_ref(c);
  double ngam_tot = 0;
  const bool proc = cell_radiation_prep(r, X.u, ngam_tot);
  if (processed) *processed = proc;
  if (ngam_tot_out) *ngam_tot_out = ngam_tot;
  if (!proc) return kPhotoOK;
  double us[5] = {0, 0, 0, 0, 0}, ws[4] = {0, 0, 0, 0};
  double nu[5], nrad[5];
  for (int n = 0; n < 5; ++n) {
    nu[n] = c.nu[n];
    nrad[n] = static_cast<double>(c.ngam[n]) * X.u.nH_avg * X.u.vol_rt;
  }
  double sigma_pht[3][5];
  photoion_sigma(nu, sigma_pht);
  PhotoStats s;
  const int err = photo_gas_cell(c.u, c.w, has_e ? &c.e : nullptr, c.ngam, nu, sigma_pht, nrad, X, s);
  if (err) return err;
  for (int m = 0; m < 5; ++m) us[m] = us[m] + static_cast<double>(c.u[m]);
  for (int m = 0; m < 4; ++m) ws[m] = ws[m] + static_cast<double>(c.w[m]);
  const double Nc = 1;
  for (int sp = 0; sp < 3; ++sp) {
    c.gamma[sp] = static_cast<float>(static_cast<double>(c.gamma[sp]) + s.gam[sp] / Nc);
    c.lambda[sp] = static_cast<float>(static_cast<double>(c.lambda[sp]) + s.lam[sp] / Nc);
  }
  photo_cell_finish(r, nrad, us, ws, has_e ? static_cast<double>(c.e) : -1.0, Nc, X.u);
  if (st) *st = s;
  return kPhotoOK;
}

}  // namespace rrk
