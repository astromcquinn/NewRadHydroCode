// Cosmology, unit conversions and time stepping.
// Ports of cosmo_tools.f90 (calc_sim_var, calc_dt_from_da, calc_da_from_dt),
// radhydro.f90 (init_variables: year; calc_timestep: RT part; radhydro:
// clock update) restricted to what the RT step needs.
#pragma once

#include <cstdint>
#include <iosfwd>

#include "radray/config.hpp"

namespace radray {

// radhydro_global.f90 time-stepping parameters (real(8) parameters; ra_max
// is initialised from the real(4) literal 0.1)
namespace stepp {
constexpr double ra_max = static_cast<double>(0.1f);
constexpr double dt_gas_factor = 0.5;
constexpr double dt_pm_factor = 0.5;
}  // namespace stepp

// Time-stepping state (radhydro_global.f90 "Sim" + "Cosmo" time variables).
struct SimClock {
  int istep = 0;
  int istep_rt = 0;
  int nstep_rt = 0;
  double da_sim = 0, dt_sim = 0;
  double dt_rt = 0, dtyr_rt = 0, dtsec_rt = 0;
  double dt_yr = 0, dt_sec = 0;
  double redshift = 0, scalefactor = 0, year = 0;
  // light-travel cells owed to / by the next hydro step (rt_substep_carry)
  double rt_carry = 0;

  // radhydro_global.f90 "Sim" speed variables. Set by summarize_gas (gas)
  // and pm_tools.f90 calc_kick (DM) in grid units; calc_timestep reads
  // cmax_gas and amax_dm and then converts all of them to km/s and
  // km/s/Gyr in place, as the original does (docs/port_notes_pm.md Q3).
  double cmax_gas = 0, vmax_gas = 0, vrms_gas = 0;
  double amax_dm = 0, arms_dm = 0, vmax_dm = 0, vrms_dm = 0;
};

// Quantities set by calc_sim_var (radhydro_global.f90 "Cosmo",
// "Unit conversions", "Gas", "RT" variables).
struct Units {
  double Hubble = 0, Tcmb = 0;
  double rho_crit = 0, rhob_avg = 0, rhom_avg = 0;
  double nH_avg = 0, nHe_avg = 0;

  double x_unit = 0, vel_unit = 0, acc_unit = 0;
  double time_unit = 0, year_unit = 0;
  double rho_unit = 0, mass_unit = 0, mp_unit = 0, Msun_unit = 0;
  double eden_unit = 0, cs_unit = 0;
  double temp1_unit = 0, temp2_unit = 0, temp3_unit = 0;
  double H_unit = 0, He_unit = 0;

  double len_gas = 0, vol_gas = 0;
  double Dgas_avg = 0, Dgas_min = 0, Dgas_max = 0, Dgas_frz = 0;
  double vgas_min = 0, vgas_max = 0;
  double Tgas_min = 0, Tgas_max = 0;
  double Tmf_min = 0, Tmf_max = 0;  // MF temperature limits (hydro smooth)
  double temp_min = 0, temp_max = 0;

  double len_rt = 0, vol_rt = 0;
  double sfr_min = 0, ngam_min = 0, fgam_min = 0, seg_min = 0;
};

// cosmo_tools.f90: calc_sim_var
void calc_sim_var(const Config& cfg, const SimClock& clk, Units& u);

// cosmo_tools.f90: calc_dt_from_da(a,da,dt,dy)
void calc_dt_from_da(double a, double da, double& dt, double& dy);

// cosmo_tools.f90: calc_da_from_dt(a,da,dt,dy)
void calc_da_from_dt(double a, double& da, double dt, double& dy);

// radhydro.f90: init_variables (redshift, scalefactor, year at z_init).
// The reference starts at zinit=300; the validation spec starts at z_init.
void init_variables(double z_init, SimClock& clk);

// radhydro.f90: calc_timestep. One implementation for every driver.
//   Original : the original: istep == 0 -> da = a*ra_max/(1-ra_max);
//              otherwise dt_sim = min(2*dt_gas_factor/cmax_gas,
//              dt_pm_factor/sqrt(amax_dm)), calc_da_from_dt, ra_max limit,
//              calc_dt_from_da.
//   Courant  : the istep > 0 branch from the first step on (hydro tests).
//   Fixed    : da_sim = da_fixed, dt from calc_dt_from_da (the form of the
//              istep == 0 branch; RT and PM tests).
// Then the RT step count (nstep_rt, dt_rt, dtyr_rt, dtsec_rt), istep + 1,
// and the in-place conversion of the speed variables of clk to km/s and
// km/s/Gyr. Prints the original's summary to *out if non-null.
enum class TimestepMode { Original, Courant, Fixed };
struct TimestepInfo {
  int ns_rt = 0;
  double dt_vg = 0, dt_ap = 0;  // years (0 in the istep == 0 / fixed branch)
  double z1 = 0, z2 = 0, a1 = 0, a2 = 0, y1 = 0, y2 = 0;
};
TimestepInfo calc_timestep(const Config& cfg, const Units& u, TimestepMode mode, double da_fixed,
                           SimClock& clk, std::ostream* out = nullptr);

// Fixed-da form used by the RT-only driver (= calc_timestep(..., Fixed,
// da_fixed, ...)). Returns ns_rt.
int calc_timestep(const Config& cfg, double da_fixed, SimClock& clk);

// Distance the radiation moves in the hydro step that just ended (nstep_rt
// RT substeps of k = rt_cells_per_substep RT cells, k dl_rt, each; 1 cell in
// the original) over the comoving light distance
// of the step, dl_c = 2c/(H0 sqrt(omegam)) (sqrt(a2) - sqrt(a1)), as in
// calc_timestep. Reads the clock at the end of the step (post_rt):
// a2 = scalefactor, a1 = scalefactor - da_sim, and the actual nstep_rt
// (after nint, the rt_substep_carry remainder and the nstep_rt >= 1 guard).
// About c_reduced (exactly so only up to the nint rounding).
double rt_path_fraction(const Config& cfg, const SimClock& clk);

// The factor by which the radiation's frequencies are scaled once per hydro
// step (ray_redshifting, moment_redshift), at the end of the step:
// fscale = (a - da)/a, the original (per time). With
// cfg.redshift_per_distance and cfg.c_reduced < 1: fscale^s with
// s = min(rt_path_fraction, 1), the redshift per distance travelled.
double redshift_fscale(const Config& cfg, const SimClock& clk);

// radhydro.f90: radhydro -- the two half-step clock updates (without hydro
// and PM), each followed by calc_sim_var.
void advance_clock_hydro_step(const Config& cfg, SimClock& clk, Units& u);

}  // namespace radray
