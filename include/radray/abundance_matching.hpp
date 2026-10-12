// Abundance matching of halo UV luminosities to observed high-redshift UV
// luminosity functions (docs/abundance_matching.md).
//
//  1. UVLuminosityFunction: Schechter or double-power-law fits, one entry
//     per (redshift, reference), read from a text table (data/uv_lf/).
//  2. abundance_match: at the redshift of an entry, M_UV(M_halo) from
//     n_halo(>M) = n_gal(<M_UV), with an analytic halo mass function
//     (mass_function.hpp), optionally corrected for dust.
//  3. LuvModel: smooth parametric L_UV(M_halo, z), fitted to the matched
//     points of all entries (fit_luv_model). The source model
//     SfrModel::AbundanceMatched (halo_sources.hpp) evaluates it per halo.
#pragma once

#include <array>
#include <string>
#include <vector>

#include "radray/mass_function.hpp"

namespace radray {

// ---------------------------------------------------------------------------
// Magnitudes
// ---------------------------------------------------------------------------

// L_nu (erg/s/Hz) of an AB absolute magnitude: L = 4 pi (10 pc)^2
// 10^(-0.4 (M + 48.60)), i.e. M_UV = 51.595 - 2.5 log10 L.
double luv_from_muv(double M_uv);
double muv_from_luv(double L_uv);

// ---------------------------------------------------------------------------
// Observed UV luminosity functions
// ---------------------------------------------------------------------------

struct LFEntry {
  enum class Form { Schechter, DoublePowerLaw, SchechterCurved };
  double z = 0;
  Form form = Form::Schechter;
  // phi_star in Mpc^-3 mag^-1 (comoving, the observers' h = 0.7 units).
  //  Schechter: phi = 0.4 ln10 phi* x^(alpha+1) exp(-x), x = 10^(-0.4(M-M*))
  //  DPL:       phi = phi* / (10^(0.4(alpha+1)(M-M*)) + 10^(0.4(beta+1)(M-M*)))
  //  SchechterCurved (Bouwens et al. 2022, ApJ 940, 55, Sec. 3.2): Schechter
  //             times 10^(-0.4 delta (M+16)^2) faintward of M = -16, with
  //             delta stored in `beta`
  double Mstar = 0, phistar = 0, alpha = 0, beta = 0;
  // Magnitude range the fit is constrained by data (bright, faint)
  double M_bright = -23, M_faint = -17;
  double weight = 1;      // relative weight in the L_UV(M, z) fit
  std::string ref;        // citation key, e.g. Bouwens+21
  bool confirmed = false; // parameters checked against the paper by the user

  // phi(M_UV) in Mpc^-3 mag^-1
  double phi(double M_uv) const;
  // Number density brighter than M_UV, Mpc^-3: integral of phi from
  // M_UV = -30 to M_uv.
  double n_brighter(double M_uv) const;
};

// Text table, one entry per line, '#' starts a comment:
//   z  form(schechter|dpl)  Mstar  phistar  alpha  beta  M_bright  M_faint
//   weight  ref  status(confirmed|unconfirmed)
// beta is ignored (write 0) for Schechter entries.
std::vector<LFEntry> read_lf_table(const std::string& path);

// ---------------------------------------------------------------------------
// Abundance matching at one redshift
// ---------------------------------------------------------------------------

// Optional dust correction of the observed magnitudes: Meurer, Heckman &
// Calzetti (1999) A_UV = 4.43 + 1.99 beta_UV, beta_UV = beta0 + dbeta_dM
// (M_obs + 19.5), A_UV >= 0. Off by default (docs: to be decided).
struct DustModel {
  bool enabled = false;
  double beta0 = -2.0;     // beta_UV at M_UV = -19.5
  double dbeta_dM = -0.1;  // d beta_UV / d M_UV
  double irx_a = 4.43, irx_b = 1.99;
  double A_uv(double M_obs) const;
};

// The simulated FoF mass function as a smooth correction to an analytic one
// (user decision 2026-10-09: match to the simulated halos), fitted by
// tools/analysis/fit_fof_hmf.py:
//   dn/dlnM(sim) = dn/dlnM(base) * r_inf / (1 + (N/N_c)^-p),
//   N = M / mpart_msun, N_c = Nc0 ((1+z)/7)^q,
// and zero below nmin particles (no halos there in the catalogs).
struct SimHMF {
  bool enabled = false;
  double mpart_msun = 0;
  double r_inf = 1, p = 1, Nc0 = 1, q = 0;
  double nmin = 20;
  double factor(double M_msun, double z) const;
  // key = value file (fit_fof_hmf.py output)
  static SimHMF read(const std::string& path);
};

struct AMOptions {
  MassFunction::Kind hmf = MassFunction::Kind::Tinker08;  // analytic HMF (base of sim_hmf)
  SimHMF sim_hmf;  // when enabled: match to the simulated halos (base: hmf)
  MassFunctionCosmology cosmo;   // simulation cosmology (cosmo_parameters)
  double tinker_zmax_evol = 1e30;
  MassFunction::STParams st;  // Sheth-Tormen coefficients (--st-params)
  DustModel dust;
  // Halo masses (Msun, not Msun/h) of the matched points
  double logM_min = 8.0, logM_max = 13.0, dlogM = 0.1;
  // The shape fit also uses matched points down to this magnitude, beyond
  // the entry's M_faint (the LF's faint end extrapolated to the halos the
  // simulation resolves, user decision); 0 = only the observed range. The
  // high-z normalisation always uses the observed range only.
  double extend_faint_to = 0;
};

struct AMPoint {
  double z = 0;
  double logM = 0;    // log10 halo mass, Msun
  double n_halo = 0;  // cumulative halo density n(>M), Mpc^-3
  double M_obs = 0;   // matched observed magnitude
  double M_int = 0;   // intrinsic (dust-corrected) magnitude
  double logL = 0;    // log10 intrinsic L_UV, erg/s/Hz
  double weight = 1;
  bool in_range = false;  // inside the observed [M_bright, M_faint]
  bool in_shape = false;  // inside [M_bright, max(M_faint, extend_faint_to)]
};

// Cumulative halo number density above M (Msun), Mpc^-3, at redshift z
// (h = the hubble0 of the mass function's cosmology).
double halo_n_above(const MassFunction& mf, double h, double M_msun, double z,
                    MassFunction::Kind k);

// M_UV(M_halo) on the logM grid of opt for one LF entry.
std::vector<AMPoint> abundance_match(const LFEntry& e, const AMOptions& opt);

// ---------------------------------------------------------------------------
// Parametric L_UV(M_halo, z)
// ---------------------------------------------------------------------------

// L_UV = 2 L_p / ((M/M_p)^-a + (M/M_p)^-b),  x = z - z0,
//   log10 L_p = l0 + l1 x + l2 x^2,   a = a0 + a1 x,
//   log10 M_p = m0 + m1 x,            b fixed by default.
// M in Msun, L in erg/s/Hz. z is clamped to [z_lo, z_hi] before evaluation
// of the shape, so the shape is held fixed outside the redshifts it is
// fitted to (the HST LFs).
//
// Above z_hi only the normalisation changes (fixed shape, user decision
// 2026-10-08): log10 L += dlogL(z), piecewise linear through (z_hi, 0) and
// the knots (norm_z[k], norm_dlogL[k]) (increasing z, all > z_hi), held at
// the last knot's value beyond it. The knots are fitted to the luminous
// galaxies of the JWST LFs (fit_highz_norm).
struct LuvModel {
  static constexpr int kNpar = 8;
  // order: l0 l1 l2 a0 a1 m0 m1 b
  std::array<double, kNpar> p = {28.0, -0.1, 0.0, 1.3, 0.0, 11.5, 0.0, 0.3};
  double z0 = 8;
  double z_lo = 0, z_hi = 1e30;
  std::vector<double> norm_z, norm_dlogL;

  double dlogL_norm(double z) const;

  double log_luv(double logM, double z) const;        // shape + norm
  double log_luv_shape(double logM, double z) const;  // without dlogL_norm
  double luv(double M_msun, double z) const;

  static const char* name(int i);
  // `key = value` lines (am_l0 ... am_b, am_z0, am_zlo, am_zhi, and the
  // comma-separated lists am_norm_z, am_norm_dl)
  std::string to_params() const;
};

struct FitResult {
  LuvModel model;
  double rms_dex = 0;  // weighted rms of log L residuals over the used points
  int npoints = 0;
  int iterations = 0;
};

// Weighted least squares in log L of the shape (log_luv_shape) over the
// in_shape points with z <= start.z_hi (Nelder-Mead), with a penalty that
// keeps L(M) increasing (b >= 0, a(z) >= 0).
// free[i] = false keeps p[i] at its starting value.
FitResult fit_luv_model(const std::vector<AMPoint>& pts, const LuvModel& start,
                        const std::array<bool, LuvModel::kNpar>& free);

// High-redshift normalisation with the shape fixed: for each distinct
// redshift z > model.z_hi among the in-range points, the knot value is the
// weighted mean of logL_AM - log_luv(logM, z) over the points at that z
// (the exact least-squares solution for an offset). Returns the model with
// norm_z / norm_dlogL set (replacing any previous knots).
LuvModel fit_highz_norm(const std::vector<AMPoint>& pts, const LuvModel& shape);

// UV LF implied by a halo mass function and an L_UV(M) relation (no
// scatter), for checking a fit against the input LF: phi(M_UV) in
// Mpc^-3 mag^-1 at the magnitudes given (observed, i.e. dust re-applied).
std::vector<double> model_lf(const LuvModel& m, double z, const std::vector<double>& M_uv,
                             const AMOptions& opt);

}  // namespace radray
