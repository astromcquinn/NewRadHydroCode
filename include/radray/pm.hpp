// Particle-mesh gravity for dark matter and gas: port of pm_tools.f90 (the
// -DMF branch by default, as the original Makefile sets; the non-MF branch is
// selectable), the DM/gravity parts of radhydro.f90 (calc_timestep) and
// hydro_tools.f90 (the time-step part of summarize_gas), calc_pk_3d from
// cosmo_tools.f90 and read_init_conditions from io_tools.f90.
//
// Conventions as in the rest of the port (docs/port_notes.md):
//   * gas cells (i,j,k) are 1-based Fortran indices, i fastest;
//   * xv_dm(6,Npart) is stored as xv[(ip-1)*6 + c], c = 0..5;
//   * particle indices ip are 1-based; 0 means "none" in the chains, as in
//     Fortran (ll[0] is a dummy);
//   * density(Ngas+2,Ngas,Ngas) keeps the in-place FFT padding.
//
// Requires FFTW3f (src/fft.cpp). See docs/port_notes_pm.md.
#pragma once

#include <array>
#include <cstdint>
#include <iosfwd>
#include <string>
#include <vector>

#include "radray/config.hpp"
#include "radray/cosmology.hpp"

namespace radray {
namespace pm {

// radhydro_global.f90 sizes used by the PM code.
struct PMConfig {
  int Ngas = 0, Hgas = 0, Ndm = 0;
  std::int64_t Npart = 0, Ncell = 0;
  double ngasr = 0, ndmr = 0;
  int N_dm_gas = 1, N_dm_gas_c = 1;
  double x_dm_gas = 1;
  // -DMF: gravity kicks the grid velocity v_gas; otherwise u_gas(2:5)
  bool mf = true;

  static PMConfig make(int Ngas, int Ndm, bool mf = true);
};

// domain_tools.f90: init_domain_gas. With one thread set_domain_gas hands
// out the domains in the order 1..Ndomain on every pass, so only the boxes
// are kept: box[n] = {i1,i2,j1,j2,k1,k2} (1-based, inclusive).
struct GasDomains {
  int Ndomain = 0;
  int Nx = 0;
  std::vector<std::array<int, 6>> box;
};
GasDomains init_domain_gas(int Ngas, int Ncpu = 1);

// DM state and PM work arrays (radhydro_global.f90 "DM" arrays, density,
// and the amax/arms/vmax/vrms globals).
struct PMState {
  PMConfig pc;
  GasDomains dom;
  std::int64_t dm_cpu[2] = {0, 0};      // dm_cpu(1:2,1)
  std::vector<float> xv;                // xv_dm(6,Npart)
  std::vector<std::int64_t> ll;         // ll_dm(Npart), index 1..Npart
  std::vector<std::int64_t> hoc;        // hoc_dm(Ngas,Ngas,Ngas)
  std::vector<std::int64_t> hoc_cpu;    // hoc_cpu(Ngas,Ngas,1)
  std::vector<float> density;           // density(Ngas+2,Ngas,Ngas)

  double amax_dm = 0, arms_dm = 0, vmax_dm = 0, vrms_dm = 0;

  float* xv_p(std::int64_t ip) { return &xv[(ip - 1) * 6]; }
  const float* xv_p(std::int64_t ip) const { return &xv[(ip - 1) * 6]; }
  std::int64_t gcell(int i, int j, int k) const {
    const std::int64_t N = pc.Ngas;
    return (i - 1) + N * ((j - 1) + N * (k - 1));
  }
  std::int64_t dcell(int i, int j, int k) const {  // density(i,j,k)
    const std::int64_t N = pc.Ngas;
    return (i - 1) + (N + 2) * ((j - 1) + N * (k - 1));
  }
};

// hydro_tools.f90 globals cmax_gas, vmax_gas, vrms_gas
struct GasSpeeds {
  double cmax_gas = 0, vmax_gas = 0, vrms_gas = 0;
};

// pm_tools.f90: init_pm (+ init_domain_gas). Zeroes all arrays.
void init_pm(PMState& s, const PMConfig& pc);

// pm_tools.f90 routines. u_gas is (5,Ngas^3), v_gas is (3,Ngas^3).
void calc_drift(PMState& s, double dt_sim);
void link_dm(PMState& s);
void cic_mass_assignment(PMState& s, const std::array<int, 6>& box);
void calc_density_field(PMState& s, const float* u_gas);
void calc_gravity_field(PMState& s, double scalefactor);
void calc_kick(PMState& s, float* u_gas, float* v_gas, double dt_sim,
               const Units& u);
// drift, link, density, gravity, kick, drift (uses clk.dt_sim and
// clk.scalefactor, which the caller has advanced by half a step).
void particle_mesh(PMState& s, float* u_gas, float* v_gas, const SimClock& clk,
                   const Units& u);

// pm_tools.f90: calc_power_spectrum. Fills pk_dm(5,Ngas) and pk_gas(5,Ngas)
// (Fortran layout pk[(k-1)*5 + c]) and leaves density/linked lists as the
// original does. write_power_spectrum writes the Fortran file format
// "(7es14.6)" for k = 2..Hgas+1.
void calc_power_spectrum(PMState& s, const float* u_gas, std::vector<double>& pk_dm,
                         std::vector<double>& pk_gas);
void write_power_spectrum(const std::string& path, const PMState& s,
                          const std::vector<double>& pk_dm,
                          const std::vector<double>& pk_gas);

// cosmo_tools.f90: calc_pk_3d(nc,field,Dk,deconvolve). field(nc+2,nc,nc)
// after a forward fft_3d; Dk(5,nc).
void calc_pk_3d(int nc, const float* field, std::vector<double>& Dk, bool deconvolve);

// pm_tools.f90: summarize_dm. Returns the reductions (x in DM cells, v in
// km/s) and prints the original's summary to *out if non-null.
struct DMSummary {
  std::array<double, 3> xavg{}, xmin{}, xmax{}, vavg{}, vmin{}, vmax{};
};
DMSummary summarize_dm(const PMState& s, const Units& u, std::ostream* out);

// The part of hydro_tools.f90:summarize_gas that feeds calc_timestep:
// clamps u_gas(1) >= Dgas_min and w_gas >= 0 (side effects of the original)
// and computes cmax_gas, vmax_gas, vrms_gas. To be replaced by the hydro
// port's summarize_gas when the full step is assembled.
GasSpeeds summarize_gas_speeds(int Ngas, float* u_gas, const float* v_gas, float* w_gas,
                               const Units& u);

// io_tools.f90: read_init_conditions. Reads <dir>/xv_dm.init and
// <dir>/u_gas.init, adds the minimum thermal energy, sets w_gas neutral.
// (The original, with Ncpu = 1, reads only xv_dm.init -- see
// docs/port_notes_pm.md Q1; this reads both, like the harness patch.)
void read_init_conditions(const std::string& dir, PMState& s, std::vector<float>& u_gas,
                          std::vector<float>& w_gas, const Units& u);

// radhydro.f90: calc_timestep, complete (PM/hydro constraints, ra_max
// limit, RT step count, unit conversions of the amax/vmax globals). If
// da_fixed > 0 the expansion step is fixed instead (as the istep == 0
// branch). Returns ns_rt. Prints the original's summary to *out.
int calc_timestep_full(const Config& cfg, PMState& s, GasSpeeds& g, SimClock& clk,
                       const Units& u, double da_fixed, std::ostream* out);

// cosmo_tools.f90: ztag (the redshift file tag set in calc_sim_var).
std::string ztag(double z);

}  // namespace pm
}  // namespace radray
