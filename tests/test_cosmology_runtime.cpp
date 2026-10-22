// Unit tests of the runtime cosmology (Config::cosmo, cosmo_params.hpp) and
// the preheating floor (Config::T_min, z_floor; docs/production_run.md
// section 10).
//   * defaults: Cosmology{} is the original's float literals, and the
//     runtime formulas reproduce the original's folded constants (omegar,
//     temp_x) bit for bit; params.txt values are rounded to real(4), so the
//     original's values written explicitly are the defaults;
//   * a Planck cosmology changes exactly the cosmology-dependent units of
//     calc_sim_var (each checked against an independent closed form) and
//     leaves the others bit-identical; the clock (init_variables,
//     calc_dt_from_da, calc_da_from_dt), the RT substep count, the halo
//     units, the particle mass and the sources' cosmic time follow it;
//   * the floor: Units::temp_floor is temp_min unless T_min > temp_min at
//     z < z_floor; floor_gas_cell lifts a cold cell to temp_floor and leaves
//     a warm one bit-identical; photo_gas_cell stops cooling at the floor.
#include <cmath>
#include <cstring>
#include <stdexcept>
#include <string>

#include "check.hpp"
#include "radray/constants.hpp"
#include "radray/cosmo_params.hpp"
#include "radray/cosmology.hpp"
#include "radray/halo_sources.hpp"
#include "radray/halos.hpp"
#include "radray/kernels/rt_cell.hpp"
#include "radray/pm.hpp"
#include "radray/rt_rates.hpp"

using namespace radray;

namespace {

bool same(double a, double b) { return std::memcmp(&a, &b, sizeof a) == 0; }

// Planck 2018 (TT,TE,EE+lowE+lensing+BAO; from memory, see
// docs/production_run.md section 10): only a cosmology different from the
// original here, so the exact values do not matter.
Cosmology planck() {
  Cosmology c;
  CHECK(set_cosmology_key(c, "omegam", "0.3111"));
  CHECK(set_cosmology_key(c, "omegab", "0.04897"));
  CHECK(set_cosmology_key(c, "hubble0", "0.6766"));
  CHECK(set_cosmology_key(c, "sigma8", "0.8102"));
  CHECK(set_cosmology_key(c, "ns", "0.9665"));
  c.derive();
  return c;
}

void test_defaults() {
  const Cosmology c;
  CHECK(same(c.omegam, static_cast<double>(0.27f)));
  CHECK(same(c.omegab, static_cast<double>(0.045f)));
  CHECK(same(c.hubble0, static_cast<double>(0.70f)));
  CHECK(same(c.sigma8, static_cast<double>(0.80f)));
  CHECK(same(c.ns, static_cast<double>(0.96f)));
  CHECK(same(c.omegal, cp::omegal));
  CHECK(same(c.omegar, cp::omegar));
  CHECK(same(c.temp_x, folded::temp_x));
  CHECK(c.is_original());
  // the runtime formulas reproduce gfortran's compile-time folding
  CHECK(same(cp::omegar_c1, static_cast<double>(4.48e-7f * 1.69f)));
  CHECK(same(cosmology_omegar(cp::hubble0), cp::omegar));
  CHECK(same(cosmology_temp_x(cp::omegab, cp::hubble0), folded::temp_x));
  // explicit original values, as a params.txt would write them
  Cosmology e;
  CHECK(set_cosmology_key(e, "omegam", "0.27"));
  CHECK(set_cosmology_key(e, "omegab", "0.045"));
  CHECK(set_cosmology_key(e, "hubble0", "0.70"));
  CHECK(set_cosmology_key(e, "sigma8", "0.8"));
  CHECK(set_cosmology_key(e, "ns", "9.6d-1"));  // Fortran exponent
  e.derive();
  CHECK(std::memcmp(&e, &c, sizeof c) == 0);
  // other keys are not consumed; bad values throw
  CHECK(!set_cosmology_key(e, "box", "5"));
  bool threw = false;
  try {
    set_cosmology_key(e, "omegam", "0.3x");
  } catch (const std::runtime_error&) {
    threw = true;
  }
  CHECK(threw);
  threw = false;
  try {
    Cosmology b;
    b.omegab = 0.5;  // > omegam
    b.derive();
  } catch (const std::runtime_error&) {
    threw = true;
  }
  CHECK(threw);
  // the Config parser: cosmology and floor keys
  Config p;
  CHECK(set_cosmology_param(p, "hubble0", "0.6766"));
  CHECK(set_cosmology_param(p, "T_min", "100"));
  CHECK(!set_cosmology_param(p, "N", "32"));
  Config r;
  threw = false;
  try {
    copy_cosmology_param(p, r);  // T_min without z_floor
  } catch (const std::runtime_error&) {
    threw = true;
  }
  CHECK(threw);
  CHECK(set_cosmology_param(p, "z_floor", "15"));
  copy_cosmology_param(p, r);
  CHECK(r.T_min == 100 && r.z_floor == 15 && !r.cosmo.is_original());
  CHECK(cosmology_summary(Config::make(8)).empty());
  CHECK(!cosmology_summary(r).empty());
}

// Independent closed forms of the cosmology-dependent units of
// calc_sim_var (cosmo_tools.f90), from the definitions, not the code.
void test_units() {
  const double z = 9.5, a = 1 / (1 + z);
  SimClock clk;
  clk.redshift = z;
  clk.scalefactor = a;
  Config c0 = Config::make(16, 4, 8, 4.0), c1 = c0;
  c1.cosmo = planck();
  Units u0, u1;
  calc_sim_var(c0, clk, u0);
  calc_sim_var(c1, clk, u1);

  for (const Config* cf : {&c0, &c1}) {
    const Cosmology& c = cf->cosmo;
    const Units& u = cf == &c0 ? u0 : u1;
    const double h = c.hubble0, om = c.omegam, ob = c.omegab;
    const double ol = 1 - om, orad = 4.48e-7 * 1.69 * std::pow(2.725, 4) / (h * h);
    const double tol = 1e-12;
    CHECK_CLOSE(c.omegar, orad, 1e-6);  // real(4) literals in the code
    CHECK_CLOSE(u.Hubble, pc::H0_cgs * h * std::sqrt(om / (a * a * a) + ol + c.omegar / (a * a * a * a)), tol);
    CHECK_CLOSE(u.x_unit, a * cf->box / h / cf->Ngas * pc::Mpc2cm, tol);
    CHECK_CLOSE(u.rho_unit, pc::rhoc_cgs * om * h * h / (a * a * a), tol);
    CHECK_CLOSE(u.time_unit, 2 / (3 * pc::H0_cgs * h * std::sqrt(om)) * a * a, tol);
    CHECK_CLOSE(u.rhob_avg, u.rho_unit * ob / om, tol);
    CHECK_CLOSE(u.nH_avg, pc::rhoc_cgs * ob * h * h / (a * a * a) * cp::XHy / pc::mH_cgs, tol);
    CHECK_CLOSE(u.Dgas_avg, ob / om, tol);
    CHECK_CLOSE(u.Msun_unit, pc::rhoc_cgs * om * std::pow(cf->box * pc::Mpc2cm, 3) / h / cf->Npart / pc::Msun_cgs, tol);
    // velocity unit x_unit/time_unit = (a box/N Mpc) (3/2) H0 sqrt(om)/a^2
    CHECK_CLOSE(u.vel_unit, cf->box / cf->Ngas * pc::Mpc2cm * 1.5 * pc::H0_cgs * std::sqrt(om) / a, tol);
    const double x = 1e3 * std::pow(ob * h * h, 0.4);
    CHECK_CLOSE(c.temp_x, x, 1e-6);
    CHECK(1 + z < x);
    CHECK_CLOSE(u.temp_min, std::fmax(cp::Tcmb0 / c.temp_x / (a * a), 1.0), tol);
  }
  // Planck changes them ...
  CHECK(u1.Hubble != u0.Hubble && u1.x_unit != u0.x_unit && u1.rho_unit != u0.rho_unit);
  CHECK(u1.time_unit != u0.time_unit && u1.nH_avg != u0.nH_avg && u1.Dgas_avg != u0.Dgas_avg);
  CHECK(u1.temp1_unit != u0.temp1_unit && u1.temp_min != u0.temp_min);
  // ... and nothing that does not depend on the cosmology
  CHECK(same(u1.Tcmb, u0.Tcmb) && same(u1.temp_max, u0.temp_max));
  CHECK(same(u1.sfr_min, u0.sfr_min) && same(u1.ngam_min, u0.ngam_min));
  CHECK(same(u1.fgam_min, u0.fgam_min) && same(u1.seg_min, u0.seg_min));
  CHECK(same(u1.Tmf_min, u0.Tmf_min) && same(u1.Tmf_max, u0.Tmf_max) && same(u1.vgas_min, u0.vgas_min));
  CHECK(same(u1.len_rt / u1.len_gas, u0.len_rt / u0.len_gas));
  CHECK(same(u1.Dgas_min / u1.Dgas_avg, u0.Dgas_min / u0.Dgas_avg));
  CHECK(same(u1.H_unit / u1.He_unit, u0.H_unit / u0.He_unit));
  // the temperature conversion scales with the specific energy unit only
  CHECK_CLOSE(u1.temp1_unit / u0.temp1_unit,
              (u1.eden_unit / u1.rho_unit) / (u0.eden_unit / u0.rho_unit), 1e-14);
  // the floor is off by default
  CHECK(same(u0.temp_floor, u0.temp_min) && same(u1.temp_floor, u1.temp_min));
}

// The clock: calc_dt_from_da / calc_da_from_dt / init_variables against an
// independent integral of dt = da/(a H(a)) for the Planck cosmology.
void test_clock() {
  const Cosmology c = planck(), o;
  auto E = [&](double x) {
    return std::sqrt(c.omegam / (x * x * x) + c.omegal + c.omegar / (x * x * x * x));
  };
  const double a = 0.05, da = 0.005;
  double s = 0;
  const int n = 20000;  // Simpson
  for (int i = 0; i <= n; ++i) {
    const double x = a + da * i / n;
    const double w = (i == 0 || i == n) ? 1 : (i % 2 ? 4 : 2);
    s += w / (x * E(x));
  }
  const double yr = s * da / n / 3 / (pc::H0_cgs * c.hubble0) * pc::sec2yr;
  double dt = 0, dy = 0, dt0 = 0, dy0 = 0;
  calc_dt_from_da(c, a, da, dt, dy);
  calc_dt_from_da(o, a, da, dt0, dy0);
  CHECK_CLOSE(dy, yr, 1e-6);
  CHECK(std::fabs(dy / dy0 - 1) > 1e-3);  // and the original differs
  // the inverse: da from dt
  double da2 = 0, dy2 = 0;
  calc_da_from_dt(c, a, da2, dt, dy2);
  CHECK_CLOSE(da2, da, 1e-6);
  CHECK_CLOSE(dy2, dy, 1e-6);
  // init_variables: the analytic matter + Lambda age (no radiation, as the
  // original)
  SimClock k;
  init_variables(c, 9, k);
  const double ol = c.omegal, om = c.omegam, a9 = 0.1;
  const double age = 2 / (3 * pc::H0_cgs * c.hubble0 * std::sqrt(ol)) *
                     std::asinh(std::sqrt(ol / om * a9 * a9 * a9)) * pc::sec2yr;
  CHECK_CLOSE(k.year, age, 1e-12);
  CHECK(same(k.redshift, 9.0) && same(k.scalefactor, 0.1));
  // the RT substeps: nint(dl_c/dl_rt), dl_c/dl_rt = 2 c N/(H0 sqrt(om) box Mpc) (sqrt(a2) - sqrt(a1))
  Config cfg = Config::make(64, 4, 32, 10.0);
  cfg.cosmo = c;
  SimClock t;
  init_variables(c, 9, t);
  const int ns = calc_timestep(cfg, 0.002, t);
  const double cells = 2 * pc::c_cgs * 64 / (pc::H0_cgs * std::sqrt(c.omegam) * 10.0 * pc::Mpc2cm) *
                       (std::sqrt(0.102) - std::sqrt(0.1));
  CHECK(ns == static_cast<int>(std::lround(cells)));
}

// Halo units, particle mass, PM config, sources' cosmic time
void test_halos_sources() {
  const Cosmology c = planck(), o;
  Config cfg = Config::make(16, 4, 8, 4.0);
  cfg.cosmo = c;
  SimClock clk;
  init_variables(c, 8, clk);
  Units u;
  calc_sim_var(cfg, clk, u);
  const HaloUnits hu = make_halo_units(cfg, clk, u);
  CHECK(same(hu.hubble0, c.hubble0) && same(hu.cosmo.omegam, c.omegam));
  CHECK(same(hu.mpart_msun, particle_mass_msun(cfg.Npart, cfg.box, c)));
  CHECK_CLOSE(particle_mass_msun(512, 4.0, c) / particle_mass_msun(512, 4.0, o),
              (c.omegam / c.hubble0) / (o.omegam / o.hubble0), 1e-14);
#ifdef RADRAY_TEST_PM
  CHECK(same(pm::PMConfig::make(cfg).cosmo.omegab, c.omegab));
  CHECK(same(pm::PMConfig::make(16, 8).cosmo.omegam, cp::omegam));
#endif
  // cosmic time from z = 50 to 8 against an independent integral
  auto E = [&](double z) {
    const double a = 1 / (1 + z);
    return std::sqrt(c.omegam / (a * a * a) + c.omegal + c.omegar / (a * a * a * a));
  };
  double s = 0;
  const int n = 20000;
  for (int i = 0; i <= n; ++i) {
    const double z = 8 + 42.0 * i / n;
    const double w = (i == 0 || i == n) ? 1 : (i % 2 ? 4 : 2);
    s += w / ((1 + z) * E(z));
  }
  const double t = s * 42.0 / n / 3 / (pc::H0_cgs * c.hubble0) * pc::sec2yr;
  CHECK_CLOSE(cosmic_time_yr(8, 50, c), t, 1e-5);
  CHECK(std::fabs(cosmic_time_yr(8, 50, o) / t - 1) > 1e-3);
}

// The floor: temp_floor, floor_gas_cell, and photo_gas_cell's clamp
void test_floor() {
  Config cfg = Config::make(16);
  SimClock clk;
  init_variables(cfg.cosmo, 10, clk);
  Units u;
  calc_sim_var(cfg, clk, u);
  CHECK(same(u.temp_floor, u.temp_min));
  cfg.T_min = 1000;
  cfg.z_floor = 12;
  calc_sim_var(cfg, clk, u);
  CHECK(u.temp_floor == 1000 && rrk::temp_floor_active(u));
  cfg.z_floor = 10;  // z < z_floor only
  calc_sim_var(cfg, clk, u);
  CHECK(same(u.temp_floor, u.temp_min) && !rrk::temp_floor_active(u));
  cfg.z_floor = 12;
  cfg.T_min = 1e-3;  // below temp_min: max(temp_min, T_min)
  calc_sim_var(cfg, clk, u);
  CHECK(same(u.temp_floor, u.temp_min) && !rrk::temp_floor_active(u));
  cfg.T_min = 1000;
  calc_sim_var(cfg, clk, u);

  // A neutral gas cell at temperature T (grid units; w = XHy rho, YHe rho)
  auto cell = [&](double T, float uu[5], float ww[4], float& e) {
    const float rho = 0.3f;
    uu[0] = rho;
    uu[1] = 0.01f;
    uu[2] = uu[3] = 0;
    ww[0] = static_cast<float>(cp::XHy) * rho;
    ww[1] = static_cast<float>(cp::YHe) * rho;
    ww[2] = ww[3] = 0;
    const double rhog = rho * u.rho_unit;
    const double neff = rhog * cp::XHy / pc::mH_cgs + rhog * cp::YHe / pc::mHe_cgs;
    const double eth = neff * pc::k_cgs * T / (cp::gamma_ideal - 1) / u.eden_unit;
    const float kin = uu[1] * uu[1] / uu[0] / 2.0f;
    uu[4] = static_cast<float>(kin + eth);
    e = static_cast<float>(eth);
  };
  auto temp_of = [&](const float uu[5], double eth_grid) {
    const double rhog = uu[0] * u.rho_unit;
    const double neff = rhog * cp::XHy / pc::mH_cgs + rhog * cp::YHe / pc::mHe_cgs;
    return (cp::gamma_ideal - 1) * eth_grid * u.eden_unit / pc::k_cgs / neff;
  };
  float uu[5], ww[4], e;
  cell(50, uu, ww, e);
  rrk::floor_gas_cell(uu, ww, &e, u);
  CHECK_CLOSE(temp_of(uu, e), 1000, 1e-6);
  const float kin = uu[1] * uu[1] / uu[0] / 2.0f;
  CHECK_CLOSE(temp_of(uu, uu[4] - kin), 1000, 1e-4);
  float hot[5], hw[4], he;
  cell(5000, hot, hw, he);
  float hot0[5], he0 = he;
  std::memcpy(hot0, hot, sizeof hot);
  rrk::floor_gas_cell(hot, hw, &he, u);
  CHECK(std::memcmp(hot0, hot, sizeof hot) == 0 && he0 == he);

  // photo_gas_cell without photons: a half-ionised cell at 2000 K cools to
  // about 1060 K in one long call; with a floor of 1500 K it stops there
  AtomicRates rates;
  calc_atomic_rates(rates, true);
  rrk::PhysCtx X;
  X.u = u;
  X.rates.cool = rates.coef_cool_.data();
  X.rates.alpha = rates.alpha_rec_.data();
  X.rates.gcol = rates.gamma_col_.data();
  X.rates.legacy_physics = true;
  X.dtsec_rt = 1e17;  // long enough to cool far below 1000 K without the floor
  X.cool_suppress = 0;
  for (const bool floor_on : {false, true}) {
    X.u.temp_floor = floor_on ? 1500.0 : u.temp_min;
    float c5[5], c4[4], ce;
    cell(2000, c5, c4, ce);
    c4[0] = 0.5f * c4[0];  // half ionised: recombination and line cooling
    const float gng[5] = {0, 0, 0, 0, 0};
    double nu[5], nrad[5] = {0, 0, 0, 0, 0}, sig[3][5];
    for (int n = 0; n < 5; ++n) nu[n] = cp::freq_init[n];
    rrk::photoion_sigma(nu, sig);
    rrk::PhotoStats st;
    CHECK(rrk::photo_gas_cell(c5, c4, &ce, gng, nu, sig, nrad, X, st) == rrk::kPhotoOK);
    // T from the dual-energy field, with the cell's new ionisation state
    const double rhog = c5[0] * u.rho_unit;
    const double nH = rhog * cp::XHy / pc::mH_cgs, nHe = rhog * cp::YHe / pc::mHe_cgs;
    const double nHI = std::fmin(c4[0] * u.H_unit, nH);
    const double ne = nH - nHI;
    const double T = (cp::gamma_ideal - 1) * ce * u.eden_unit / pc::k_cgs / (nH + nHe + ne);
    if (floor_on)
      CHECK_CLOSE(T, 1500, 1e-5);
    else
      CHECK(T < 1400);
  }
}

}  // namespace

int main() {
  test_defaults();
  test_units();
  test_clock();
  test_halos_sources();
  test_floor();
  TEST_MAIN_END();
}
