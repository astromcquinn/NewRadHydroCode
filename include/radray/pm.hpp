// Particle-mesh gravity for dark matter and gas: port of pm_tools.f90 (the
// -DMF branch by default, as the original Makefile sets; the non-MF branch is
// selectable), calc_pk_3d from cosmo_tools.f90 and read_init_conditions from
// io_tools.f90.
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
#include "radray/simulation.hpp"

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
  static PMConfig make(const Config& cfg, bool mf = true) { return make(cfg.Ngas, cfg.Ndm, mf); }
};

// DM state and PM work arrays (radhydro_global.f90 "DM" arrays, density,
// The amax/arms/vmax/vrms_dm globals live in SimClock.
struct PMState {
  PMConfig pc;
  // domain_tools.f90 gas domains (init_domain_gas in hydro.hpp). With one
  // thread set_domain_gas hands out the domains in the order 1..Ndomain on
  // every pass, so the PM loops over dom.box in that order.
  Domains dom;
  std::int64_t dm_cpu[2] = {0, 0};      // dm_cpu(1:2,1)
  std::vector<float> xv;                // xv_dm(6,Npart)
  std::vector<std::int64_t> ll;         // ll_dm(Npart), index 1..Npart
  std::vector<std::int64_t> hoc;        // hoc_dm(Ngas,Ngas,Ngas)
  std::vector<std::int64_t> hoc_cpu;    // hoc_cpu(Ngas,Ngas,1)
  std::vector<float> density;           // density(Ngas+2,Ngas,Ngas)

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

// pm_tools.f90: init_pm. Zeroes all arrays; takes the gas domains of
// init_domain_gas.
void init_pm(PMState& s, const PMConfig& pc, const Domains& dom);

// pm_tools.f90 routines. u_gas is (5,Ngas^3), v_gas is (3,Ngas^3).
void calc_drift(PMState& s, double dt_sim);
void link_dm(PMState& s);
void cic_mass_assignment(PMState& s, const std::array<int, 6>& box);
// cic_mass_assignment over all domains in order 1..Ndomain (parallel,
// bit-identical to the serial loop)
void cic_mass_assignment_all(PMState& s);
void calc_density_field(PMState& s, const float* u_gas);
void calc_gravity_field(PMState& s, double scalefactor);
// Sets clk.amax_dm, arms_dm, vmax_dm, vrms_dm (grid units).
void calc_kick(PMState& s, float* u_gas, float* v_gas, SimClock& clk, const Units& u);
// drift, link, density, gravity, kick, drift (uses clk.dt_sim and
// clk.scalefactor, which the caller has advanced by half a step).
void particle_mesh(PMState& s, float* u_gas, float* v_gas, SimClock& clk, const Units& u);
// hydro_scheme = unsplit (docs/hydro_unsplit.md section 3): the same, with
// the gas kick on the conserved lab-frame variables (momentum and total
// energy; the thermal energy, selected by the dual-energy rule with eta1
// from the internal-energy field e_gas, is unchanged). v_gas is not used.
void calc_kick_unsplit(PMState& s, float* u_gas, const float* e_gas, double eta1, SimClock& clk,
                       const Units& u);
void particle_mesh_unsplit(PMState& s, float* u_gas, const float* e_gas, double eta1,
                           SimClock& clk, const Units& u);

// pm_tools.f90: calc_power_spectrum. Fills pk_dm(5,Ngas) and pk_gas(5,Ngas)
// (Fortran layout pk[(k-1)*5 + c]) and leaves density/linked lists as the
// original does. write_power_spectrum writes the Fortran file format
// "(7es14.6)" for k = 2..Hgas+1.
void calc_power_spectrum(PMState& s, const float* u_gas, double box,
                         std::vector<double>& pk_dm, std::vector<double>& pk_gas);
void write_power_spectrum(const std::string& path, const PMState& s,
                          const std::vector<double>& pk_dm,
                          const std::vector<double>& pk_gas);

// cosmo_tools.f90: calc_pk_3d(nc,field,Dk,deconvolve). field(nc+2,nc,nc)
// after a forward fft_3d; Dk(5,nc); box in Mpc/h.
void calc_pk_3d(int nc, const float* field, std::vector<double>& Dk, bool deconvolve,
                double box);

// pm_tools.f90: summarize_dm. Returns the reductions (x in DM cells, v in
// km/s) and prints the original's summary to *out if non-null.
struct DMSummary {
  std::array<double, 3> xavg{}, xmin{}, xmax{}, vavg{}, vmin{}, vmax{};
};
DMSummary summarize_dm(const PMState& s, const Units& u, std::ostream* out);

// hydro_tools.f90 summarize_gas (the full routine, also for gravity-only
// runs) and radhydro.f90 calc_timestep are in hydro.hpp and cosmology.hpp.

// io_tools.f90: read_init_conditions. Reads <dir>/xv_dm.init and
// <dir>/u_gas.init, adds the minimum thermal energy, sets w_gas neutral.
// (The original, with Ncpu = 1, reads only xv_dm.init -- see
// docs/port_notes_pm.md Q1; this reads both, like the harness patch.)
void read_init_conditions(const std::string& dir, PMState& s, std::vector<float>& u_gas,
                          std::vector<float>& w_gas, const Units& u);

// cosmo_tools.f90: ztag (the redshift file tag set in calc_sim_var).
std::string ztag(double z);

}  // namespace pm
}  // namespace radray
