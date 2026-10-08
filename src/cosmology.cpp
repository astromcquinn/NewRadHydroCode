#include "radray/cosmology.hpp"

#include <cmath>

#include "radray/constants.hpp"
#include "radray/fortran.hpp"

namespace radray {

using namespace pc;
using namespace cp;

Config Config::make(int N, int nray_max_factor) {
  Config c;
  c.Ncpu = 1;
  c.Ngas = N;
  c.Nrt = N;
  c.Hrt = c.Nrt / 2;
  c.Ndm = c.Ngas / 2;
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
  c.x_gas_box = box / c.ngasr;
  c.x_rt_box = box / c.nrtr;
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
  u.x_unit = (a * box / h0 / cfg.Ngas) * Mpc2cm;
  u.rho_unit = rhoc_cgs * om * f::pow2(h0) / f::pow3(a);
  u.mass_unit = u.rho_unit * f::pow3(u.x_unit);
  u.time_unit = 2 / (3 * H0_cgs * h0 * std::sqrt(om)) * f::pow2(a);
  u.year_unit = u.time_unit * sec2yr;
  u.vel_unit = u.x_unit / u.time_unit;
  u.acc_unit = u.x_unit / f::pow2(u.time_unit);
  u.eden_unit = u.rho_unit * f::pow2(u.vel_unit);
  u.cs_unit = u.vel_unit;
  // (mp_unit, Msun_unit: particle units, not used by RT -> omitted)

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

// radhydro.f90: calc_timestep (RT part; da_sim fixed by the spec)
int calc_timestep(const Config& cfg, double da_fixed, SimClock& clk) {
  const double a = clk.scalefactor;

  clk.da_sim = da_fixed;
  calc_dt_from_da(a, clk.da_sim, clk.dt_sim, clk.dt_yr);
  clk.dt_sec = clk.dt_yr * yr2sec;

  // RT time step
  const double cscale = 1;
  const double dl_c = 2 * c_cgs / (H0_cgs * hubble0 * std::sqrt(omegam)) *
                      (std::sqrt(a + clk.da_sim) - std::sqrt(a));
  const double dl_rt = box / hubble0 * Mpc2cm / cfg.Nrt;
  const int ns_rt = f::nint(dl_c / dl_rt);
  clk.nstep_rt = f::nint(dl_c / dl_rt * cscale);
  clk.dt_rt = clk.dt_sim / clk.nstep_rt;
  clk.dtyr_rt = clk.dt_yr * (clk.dt_rt / clk.dt_sim);
  clk.dtsec_rt = clk.dt_sec * (clk.dt_rt / clk.dt_sim);

  clk.istep = clk.istep + 1;
  return ns_rt;
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
