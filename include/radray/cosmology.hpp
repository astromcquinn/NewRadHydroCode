// Cosmology, unit conversions and time stepping.
// Ports of cosmo_tools.f90 (calc_sim_var, calc_dt_from_da, calc_da_from_dt),
// radhydro.f90 (init_variables: year; calc_timestep: RT part; radhydro:
// clock update) restricted to what the RT step needs.
#pragma once

#include <cstdint>

#include "radray/config.hpp"

namespace radray {

// Time-stepping state (radhydro_global.f90 "Sim" + "Cosmo" time variables).
struct SimClock {
  int istep = 0;
  int istep_rt = 0;
  int nstep_rt = 0;
  double da_sim = 0, dt_sim = 0;
  double dt_rt = 0, dtyr_rt = 0, dtsec_rt = 0;
  double dt_yr = 0, dt_sec = 0;
  double redshift = 0, scalefactor = 0, year = 0;
};

// Quantities set by calc_sim_var (radhydro_global.f90 "Cosmo",
// "Unit conversions", "Gas", "RT" variables).
struct Units {
  double Hubble = 0, Tcmb = 0;
  double rho_crit = 0, rhob_avg = 0, rhom_avg = 0;
  double nH_avg = 0, nHe_avg = 0;

  double x_unit = 0, vel_unit = 0, acc_unit = 0;
  double time_unit = 0, year_unit = 0;
  double rho_unit = 0, mass_unit = 0;
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

// radhydro.f90: calc_timestep, with the hydro/PM constraints replaced by a
// fixed da_sim (validation spec). Sets dt_sim, dt_yr, dt_sec, nstep_rt,
// dt_rt, dtyr_rt, dtsec_rt and increments istep. Returns ns_rt (log only).
int calc_timestep(const Config& cfg, double da_fixed, SimClock& clk);

// radhydro.f90: radhydro -- the two half-step clock updates (without hydro
// and PM), each followed by calc_sim_var.
void advance_clock_hydro_step(const Config& cfg, SimClock& clk, Units& u);

}  // namespace radray
