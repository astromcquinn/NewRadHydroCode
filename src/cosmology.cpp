#include "radray/cosmology.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <ostream>

#include "radray/constants.hpp"
#include "radray/fortran.hpp"

namespace radray {

using namespace pc;
using namespace cp;

Config Config::make(int N, int nray_max_factor, int Ndm, double box_mpch) {
  Config c;
  c.Ncpu = 1;
  c.Ngas = N;
  c.Nrt = N;
  c.Hrt = c.Nrt / 2;
  c.Ndm = Ndm > 0 ? Ndm : c.Ngas / 2;
  c.Npart = static_cast<std::int64_t>(c.Ndm) * c.Ndm * c.Ndm;
  c.box = box_mpch;
  c.Ncell = static_cast<std::int64_t>(c.Ngas) * c.Ngas * c.Ngas;
  c.Nrtc = static_cast<std::int64_t>(c.Nrt) * c.Nrt * c.Nrt;
  c.ngasr = c.Ngas;
  c.nrtr = c.Nrt;
  c.hrtr = c.Hrt;
  c.nrtcr = static_cast<double>(c.Nrtc);
  c.ncellr = static_cast<double>(c.Ncell);
  c.N_rt_gas = c.Ngas / c.Nrt;
  c.N_rt_gas_c = c.N_rt_gas * c.N_rt_gas * c.N_rt_gas;
  c.x_rt_gas = c.N_rt_gas;
  c.x_gas_box = c.box / c.ngasr;
  c.x_dm_box = c.box / static_cast<double>(c.Ndm);
  c.x_rt_box = c.box / c.nrtr;
  c.nray_max_factor = nray_max_factor;
  c.Nray_max = nray_max_factor * c.Nrtc;
  c.Nfree_max = c.Nray_max / (4 * c.Ncpu);
  return c;
}

// cosmo_tools.f90: calc_sim_var
void calc_sim_var(const Config& cfg, const SimClock& clk, Units& u) {
  const double om = omegam, ol = omegal, ob = omegab, orr = omegar;
  const double h0 = hubble0, gam = gamma_ideal;

  const double z = clk.redshift;
  const double a = clk.scalefactor;
  // a**(3*(1+wde)) with wde = -1: exponent folds to 0.0 -> exactly 1.
  const double a_w = std::pow(a, 3 * (1 + wde));
  const double hsq = om / f::pow3(a) + ol / a_w + orr / f::pow4(a);
  u.Hubble = H0_cgs * h0 * std::sqrt(hsq);
  u.rho_crit = 3 * f::pow2(u.Hubble) / (8 * pi * G_cgs);
  u.Tcmb = Tcmb0 / a;

  // Grid units
  u.x_unit = (a * cfg.box / h0 / cfg.Ngas) * Mpc2cm;
  u.rho_unit = rhoc_cgs * om * f::pow2(h0) / f::pow3(a);
  u.mass_unit = u.rho_unit * f::pow3(u.x_unit);
  u.time_unit = 2 / (3 * H0_cgs * h0 * std::sqrt(om)) * f::pow2(a);
  u.year_unit = u.time_unit * sec2yr;
  u.vel_unit = u.x_unit / u.time_unit;
  u.acc_unit = u.x_unit / f::pow2(u.time_unit);
  u.eden_unit = u.rho_unit * f::pow2(u.vel_unit);
  u.cs_unit = u.vel_unit;
  // Particle units: mp_unit = rhoc_cgs*om*(box*Mpc2cm)**3/h0/Npart. In the
  // original box is a parameter and (box*Mpc2cm)**3 is folded at compile
  // time; only the halo finder uses these (C++ in both codes).
  u.mp_unit = rhoc_cgs * om * f::pow3(cfg.box * Mpc2cm) / h0 / static_cast<double>(cfg.Npart);
  u.Msun_unit = u.mp_unit / Msun_cgs;

  // Cell sizes in cgs
  u.len_gas = u.x_unit;
  u.vol_gas = f::pow3(u.len_gas);
  u.len_rt = u.x_unit * cfg.x_rt_gas;
  u.vol_rt = f::pow3(u.len_rt);

  // Average densities in cgs
  u.rhom_avg = u.rho_unit;
  u.rhob_avg = u.rhom_avg * ob / om;
  u.nH_avg = u.rhob_avg * XHy / mH_cgs;
  u.nHe_avg = u.rhob_avg * YHe / mHe_cgs;
  u.H_unit = u.rho_unit / mH_cgs;
  u.He_unit = u.rho_unit / mHe_cgs;

  // Temperature conversion (4. is a real(4) literal, exact)
  u.temp1_unit = u.eden_unit / u.rho_unit * (gam - 1) / k_cgs * mH_cgs /
                 (1 * XHy + 1 * YHe / 4.0);
  u.temp2_unit = u.eden_unit / u.rho_unit * (gam - 1) / k_cgs * mH_cgs /
                 (2 * XHy + 2 * YHe / 4.0);
  u.temp3_unit = u.eden_unit / u.rho_unit * (gam - 1) / k_cgs * mH_cgs /
                 (2 * XHy + 3 * YHe / 4.0);

  // Matter temperature in cgs; x = 1E3*(ob*h0**2)**0.4 (compile-time folded)
  const double x = folded::temp_x;
  if (1 + z > x) {
    u.temp_max = u.Tcmb;
    u.temp_min = u.Tcmb;
  } else {
    u.temp_max = folded::ten_lgT_max;
    u.temp_min = f::max(Tcmb0 / x / f::pow2(a), 1.0);
  }

  // Density limits in grid units
  u.Dgas_avg = ob / om;
  u.Dgas_min = u.Dgas_avg * folded::ten_lgD_min;
  u.Dgas_max = u.Dgas_avg * folded::ten_lgD_max;
  u.Dgas_frz = u.Dgas_avg * static_cast<double>(0.01f);

  // Velocity limits in grid units (3E8 is exact in real(4))
  u.vgas_min = 0;
  u.vgas_max = static_cast<double>(3e8f) / u.vel_unit;

  u.Tgas_min = u.temp_min / u.temp1_unit;
  u.Tgas_max = u.temp_max / u.temp1_unit;

  // MF temperature limits in grid units
  u.Tmf_min = 1;
  u.Tmf_max = 1;

  // Radiation density thresholds: single-precision literals!
  u.sfr_min = static_cast<double>(1e-6f);
  u.ngam_min = static_cast<double>(1e-6f);
  u.fgam_min = static_cast<double>(1e-4f);
  u.seg_min = static_cast<double>(1e-4f);
}

// cosmo_tools.f90: calc_dt_from_da
void calc_dt_from_da(double a, double da, double& dt, double& dy) {
  const double w = wde, om = omegam, ol = omegal, orr = omegar, h0 = hubble0;
  const std::int32_t n = f::max<std::int32_t>(100, f::int4(da / 1e-4));
  const double dx = da / n;
  dt = 0;
  dy = 0;
  for (std::int32_t i = 1; i <= n; ++i) {
    // (i-0.5): integer minus real(4) literal -> real(4), exact here
    const double x = a + static_cast<double>(static_cast<float>(i) - 0.5f) * dx;
    const double E = std::sqrt(om / f::pow3(x) + ol / std::pow(x, 3 * (1 + w)) +
                               orr / f::pow4(x));
    const double H = H0_cgs * h0 * E;
    dy = dy + 1 / (x * H) * dx;
    dt = dt + 3 * std::sqrt(om) / (2 * f::pow3(x) * E) * dx;
  }
  dy = dy * sec2yr;
}

// cosmo_tools.f90: calc_da_from_dt
void calc_da_from_dt(double a, double& da, double dt, double& dy) {
  const double w = wde, om = omegam, ol = omegal, orr = omegar, h0 = hubble0;
  const std::int32_t n = f::max<std::int32_t>(100, f::int4(dt / 1e-4));
  const double ds = dt / n;
  da = 0;
  dy = 0;
  for (std::int32_t i = 1; i <= n; ++i) {
    double x = a + da;
    double E = std::sqrt(om / f::pow3(x) + ol / std::pow(x, 3 * (1 + w)) +
                         orr / f::pow4(x));
    x = x + 2 * f::pow3(x) * E / (3 * std::sqrt(om)) * (ds / 2);
    E = std::sqrt(om / f::pow3(x) + ol / std::pow(x, 3 * (1 + w)) +
                  orr / f::pow4(x));
    const double H = H0_cgs * h0 * E;
    da = da + 2 * f::pow3(x) * E / (3 * std::sqrt(om)) * ds;
    dy = dy + 2 * f::pow2(x) * E / (3 * H * std::sqrt(om)) * ds;
  }
  dy = dy * sec2yr;
}

// radhydro.f90: init_variables
void init_variables(double z_init, SimClock& clk) {
  const double om = omegam, ol = omegal, h0 = hubble0;
  clk = SimClock{};
  clk.redshift = z_init;
  clk.scalefactor = 1 / (1 + z_init);  // ainit = 1/(1+zinit)
  clk.istep = 0;
  const double a = clk.scalefactor;
  if (ol > 0) {
    clk.year = 2 / (3 * H0_cgs * h0 * std::sqrt(ol)) *
               std::log(std::sqrt(ol / om * f::pow3(a)) +
                        std::sqrt(1 + ol / om * f::pow3(a))) *
               sec2yr;
  } else {
    clk.year = 2 / (3 * H0_cgs * h0) * std::sqrt(f::pow3(a) / om) * sec2yr;
  }
}

// radhydro.f90: calc_timestep, without the speed conversions and output
namespace {
TimestepInfo timestep_core(const Config& cfg, TimestepMode mode, double da_fixed, SimClock& clk) {
  using namespace stepp;
  TimestepInfo t;
  const double z = clk.redshift;
  const double a = clk.scalefactor;

  // Time steps. There are two hydro steps per PM step
  const bool first = (mode == TimestepMode::Fixed) ||
                     (mode == TimestepMode::Original && clk.istep == 0);
  if (first) {
    if (mode == TimestepMode::Fixed)
      clk.da_sim = da_fixed;
    else  // Restrict expansion timestep based on da/a ratio
      clk.da_sim = a * cfg.ra_max / (1 - cfg.ra_max);
    calc_dt_from_da(a, clk.da_sim, clk.dt_sim, clk.dt_yr);
    t.dt_vg = 0;
    t.dt_ap = 0;
    clk.dt_sec = clk.dt_yr * yr2sec;
  } else {
    // dt_ap = +Infinity without particles (amax_dm = 0)
    t.dt_vg = 2 * dt_gas_factor / clk.cmax_gas;
    t.dt_ap = dt_pm_factor / std::sqrt(clk.amax_dm);
    clk.dt_sim = f::gmin(t.dt_vg, t.dt_ap);
    calc_da_from_dt(a, clk.da_sim, clk.dt_sim, clk.dt_yr);
    clk.da_sim = f::gmin(clk.da_sim, a * cfg.ra_max / (1 - cfg.ra_max));
    calc_dt_from_da(a, clk.da_sim, clk.dt_sim, clk.dt_yr);
    t.dt_vg = clk.dt_yr * (t.dt_vg / clk.dt_sim);
    t.dt_ap = clk.dt_yr * (t.dt_ap / clk.dt_sim);
    clk.dt_sec = clk.dt_yr * yr2sec;
  }

  // Caps on the hydro step (config.hpp; docs/zoom_turbulence.md section
  // 4.4, docs/plane_parallel.md), all off by default:
  //   plane_z_re: the step before the turn-on redshift lands on it exactly,
  //               so the RT starts at z_re;
  //   dt_max_myr: the step is at most dt_max_myr Myr;
  //   rt_substeps_max: the step holds at most this many RT substeps (the
  //               light crossing of rt_substeps_max * rt_cells_per_substep
  //               RT cells at c_reduced).
  // Each one shortens da_sim; dt follows from da as above.
  {
    double da_cap = clk.da_sim;
    // plane_parallel = 1: the caps act from the turn-on (the RT is off
    // before plane_z_re; the step before it lands on it)
    bool caps = true;
    if (cfg.plane_parallel && cfg.plane_z_re >= 0) {
      const double a_re = 1 / (1 + cfg.plane_z_re);
      if (a < a_re * (1 - 1e-9)) {
        caps = false;
        if (a + da_cap > a_re) da_cap = a_re - a;
      }
    }
    if (caps && cfg.rt_substeps_max > 0) {
      // dl_c(da) = 2c/(H0 h sqrt(om)) (sqrt(a+da) - sqrt(a)) = nmax dl_sub / cscale
      const double dl_sub = cfg.box / hubble0 * Mpc2cm / cfg.Nrt * cfg.rt_cells_per_substep;
      const double x = cfg.rt_substeps_max * dl_sub / cfg.c_reduced * (H0_cgs * hubble0 * std::sqrt(omegam)) / (2 * c_cgs);
      const double sa = std::sqrt(a) + x;
      da_cap = f::gmin(da_cap, sa * sa - a);
    }
    if (caps && cfg.dt_max_myr > 0 && clk.dt_yr > cfg.dt_max_myr * 1e6) {
      // the da whose dt_yr is the cap: bisection (dt_yr grows with da)
      double lo = 0, hi = da_cap, dts = 0, dty = 0;
      for (int it = 0; it < 80; ++it) {
        const double mid = 0.5 * (lo + hi);
        calc_dt_from_da(a, mid, dts, dty);
        if (dty > cfg.dt_max_myr * 1e6)
          hi = mid;
        else
          lo = mid;
      }
      da_cap = f::gmin(da_cap, lo);
    }
    if (da_cap < clk.da_sim) {
      clk.da_sim = da_cap;
      calc_dt_from_da(a, clk.da_sim, clk.dt_sim, clk.dt_yr);
      clk.dt_sec = clk.dt_yr * yr2sec;
    }
  }

  // RT time step. Rays move one RT cell per substep, so the number of
  // substeps per hydro step sets the speed of light: cscale = c_reduced < 1
  // is the reduced speed of light approximation (the original hard-wires
  // cscale = 1). See docs/production_run.md. With rt_cells_per_substep = k
  // rays move k cells per substep (docs/rt_multicell.md): nstep_rt counts
  // the light crossings of k cells, and the carry and the >= 1 guard act on
  // that count.
  const double cscale = cfg.c_reduced;
  const double dl_c = 2 * c_cgs / (H0_cgs * hubble0 * std::sqrt(omegam)) *
                      (std::sqrt(a + clk.da_sim) - std::sqrt(a));
  const double dl_rt = cfg.box / hubble0 * Mpc2cm / cfg.Nrt;
  const double dl_sub = dl_rt * cfg.rt_cells_per_substep;  // = dl_rt for k = 1
  t.ns_rt = f::nint(dl_c / dl_rt);
  clk.nstep_rt = f::nint(dl_c / dl_sub * cscale);
  if (cfg.rt_substep_carry) {
    // light substeps of this step plus the remainder of the previous steps
    const double want = dl_c / dl_sub * cscale + clk.rt_carry;
    clk.nstep_rt = f::max(f::nint(want), 1);
    clk.rt_carry = want - clk.nstep_rt;
  }
  if (clk.nstep_rt < 1) {
    // Light (at c_reduced) moves less than half an RT cell in this hydro
    // step. At least one substep is needed; rays then move one cell per
    // hydro step, i.e. faster than c_reduced (same guard in the Fortran
    // harnesses; never reached with c_reduced = 1 at the validated sizes).
    std::fprintf(stderr,
                 " calc_timestep: warning: nstep_rt = %d < 1 (light moves %.3g RT cells in "
                 "this step at c_reduced = %g, rays %d per substep); using 1, so rays move "
                 "faster than c_reduced\n",
                 clk.nstep_rt, dl_c / dl_rt * cscale, cscale, cfg.rt_cells_per_substep);
    clk.nstep_rt = 1;
  }
  clk.dt_rt = clk.dt_sim / clk.nstep_rt;
  clk.dtyr_rt = clk.dt_yr * (clk.dt_rt / clk.dt_sim);
  clk.dtsec_rt = clk.dt_sec * (clk.dt_rt / clk.dt_sim);

  // Update timing variables
  clk.istep = clk.istep + 1;
  t.a1 = a;
  t.a2 = a + clk.da_sim;
  t.z1 = z;
  t.z2 = 1 / t.a2 - 1;
  t.y1 = clk.year;
  t.y2 = clk.year + clk.dt_yr;
  return t;
}
}  // namespace

TimestepInfo calc_timestep(const Config& cfg, const Units& u, TimestepMode mode, double da_fixed,
                           SimClock& clk, std::ostream* out) {
  const TimestepInfo t = timestep_core(cfg, mode, da_fixed, clk);

  // Velocity in km/s, acceleration in km/s/Gyr (in place, as the original;
  // 1E5 and 1E9 are real(4) literals, exact)
  const double c5 = static_cast<double>(1e5f), c9 = static_cast<double>(1e9f);
  clk.vmax_gas = clk.vmax_gas * u.vel_unit / c5;
  clk.vrms_gas = clk.vrms_gas * u.vel_unit / c5;
  clk.vmax_dm = clk.vmax_dm * u.vel_unit / c5;
  clk.vrms_dm = clk.vrms_dm * u.vel_unit / c5;
  clk.amax_dm = clk.amax_dm * u.acc_unit / c5 / sec2yr * c9;
  clk.arms_dm = clk.arms_dm * u.acc_unit / c5 / sec2yr * c9;

  if (out) {
    char buf[640];
    auto r = [](double x) { return static_cast<double>(static_cast<float>(x)); };
    std::snprintf(buf, sizeof buf,
                  "\n Step number  : %d\n Redshift     : %15.8g %15.8g %15.8g\n"
                  " Scale factor : %15.8g %15.8g %15.8g %15.8g\n"
                  " Year         : %15.8g %15.8g %15.8g %15.8g\n",
                  clk.istep, r(t.z1), r(t.z2), r(t.z2 - t.z1), r(t.a1), r(t.a2), r(t.a2 - t.a1),
                  r(clk.da_sim / t.a2), r(t.y1), r(t.y2), r(t.y2 - t.y1), r(clk.dt_yr / t.y2));
    *out << buf;
    std::snprintf(buf, sizeof buf,
                  " Velocity Gas : %15.8g %15.8g %15.8g\n Velocity DM  : %15.8g %15.8g\n"
                  " Acceleration : %15.8g %15.8g %15.8g\n RT           : %d %d %15.8g\n\n",
                  r(clk.vrms_gas), r(clk.vmax_gas), r(t.dt_vg), r(clk.vrms_dm), r(clk.vmax_dm),
                  r(clk.arms_dm), r(clk.amax_dm), r(t.dt_ap), clk.nstep_rt, t.ns_rt,
                  r(clk.dtyr_rt));
    *out << buf;
  }
  return t;
}

int calc_timestep(const Config& cfg, double da_fixed, SimClock& clk) {
  return timestep_core(cfg, TimestepMode::Fixed, da_fixed, clk).ns_rt;
}

double rt_path_fraction(const Config& cfg, const SimClock& clk) {
  // dl_c and dl_rt as in timestep_core (computed here, not stored in the
  // clock, so that SimClock and the checkpoint format are unchanged)
  const double a1 = clk.scalefactor - clk.da_sim;
  const double dl_c = 2 * c_cgs / (H0_cgs * hubble0 * std::sqrt(omegam)) *
                      (std::sqrt(a1 + clk.da_sim) - std::sqrt(a1));
  const double dl_rt = cfg.box / hubble0 * Mpc2cm / cfg.Nrt;
  // k = rt_cells_per_substep RT cells per substep
  return clk.nstep_rt * cfg.rt_cells_per_substep * dl_rt / dl_c;
}

double redshift_fscale(const Config& cfg, const SimClock& clk) {
  // the original: (a - da)/a per hydro step
  const double fscale = (clk.scalefactor - clk.da_sim) / clk.scalefactor;
  if (!cfg.redshift_per_distance || !(cfg.c_reduced < 1)) return fscale;
  // Per distance: the radiation moved s = nstep_rt*dl_rt/dl_c of the light
  // distance of the step, so it is redshifted by fscale^s (ln nu ~ -distance
  // travelled). s <= 1: a step whose substeps cover more than the light
  // distance (nint rounding, the nstep_rt >= 1 guard) keeps the full
  // per-time factor; with c_reduced >= 1 (above) the radiation moves at c
  // and the per-time factor is already per distance.
  const double s_path = std::min(rt_path_fraction(cfg, clk), 1.0);
  return std::pow(fscale, s_path);
}

// radhydro.f90: radhydro (clock part only)
void advance_clock_hydro_step(const Config& cfg, SimClock& clk, Units& u) {
  for (int half = 0; half < 2; ++half) {
    clk.year = clk.year + clk.dt_yr / 2;
    clk.scalefactor = clk.scalefactor + clk.da_sim / 2;
    clk.redshift = 1 / clk.scalefactor - 1;
    calc_sim_var(cfg, clk, u);
  }
}

}  // namespace radray
