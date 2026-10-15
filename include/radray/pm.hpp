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
#include <memory>
#include <string>
#include <vector>

#include "radray/config.hpp"
#include "radray/cosmology.hpp"
#include "radray/simulation.hpp"

namespace radray {

struct AmrexPM;  // src/amrex/pm_level.cpp (pm_backend = amrex)

namespace pm {

// radhydro_global.f90 sizes used by the PM code.
//
// pm_refine = r (docs/pm_refine.md; 1 = the original): the PM mesh has
// Nm = r*Ngas cells per dimension. The gas density is spread onto the mesh
// (every mesh cell takes the density of the gas cell it lies in) and the
// gas is kicked by the average of the mesh forces over the r^3 mesh cells
// of its cell, divided by r (mesh to gas cell units): the two operations
// are adjoint, so the gas-DM forces stay pairwise antisymmetric and the
// total momentum is conserved. With r = 1 the mesh is the gas grid and
// every expression is the original's.
struct PMConfig {
  int Ngas = 0, Hgas = 0, Ndm = 0;
  int refine = 1;      // pm_refine
  int Nm = 0, Hm = 0;  // PM mesh: Nm = refine*Ngas, Hm = Nm/2
  std::int64_t Npart = 0, Ncell = 0, Nmcell = 0;  // Ndm^3, Ngas^3, Nm^3
  double ngasr = 0, ndmr = 0, nmr = 0;
  // Particles per mesh cell spacing: Nm/Ndm (an integer, >= 1). These are
  // the original's N_dm_gas, N_dm_gas_c and x_dm_gas (Ngas/Ndm) with r = 1.
  int N_dm_m = 1, N_dm_m_c = 1;
  double x_dm_m = 1;
  // Ngas/Ndm as a real (DM to gas grid velocity units; < 1 when Ndm > Ngas)
  double x_dm_gas = 1;
  // -DMF: gravity kicks the grid velocity v_gas; otherwise u_gas(2:5)
  bool mf = true;
  // pm_backend = amrex (pm_backend.hpp, docs/stage_e_pm.md): particle_mesh
  // and particle_mesh_unsplit run on the AMReX level instead of the host
  // code below (the drifts stay on the host). Serial by default.
  PMBackend backend = PMBackend::Serial;
  int amrex_max_grid_size = 16;  // the mesh's box size (amrex_pm_max_grid_size)
  int amrex_fft = 0;             // amrex_pm_fft: 0 pencil, 1 amrex::FFT
  // amrex_state = 1 (amrex_state.hpp): the gas on the shared decomposition,
  // the mesh boxes pm_refine times the gas boxes
  bool amrex_state = false;
  // pm_gas_levels with amrex_state = 1 (config.hpp): the refined levels of
  // the gas state enter the deposit and are kicked on their own cells
  bool gas_levels = false;

  static PMConfig make(int Ngas, int Ndm, bool mf = true, int refine = 1);
  static PMConfig make(const Config& cfg, bool mf = true) {
    PMConfig c = make(cfg.Ngas, cfg.Ndm, mf, cfg.pm_refine);
    c.backend = cfg.pm_backend;
    c.amrex_max_grid_size = cfg.amrex_pm_max_grid_size > 0 ? cfg.amrex_pm_max_grid_size
                                                           : cfg.amrex_max_grid_size;
    c.amrex_fft = cfg.amrex_pm_fft;
    c.amrex_state = cfg.amrex_state;
    c.gas_levels = cfg.amrex_state && cfg.pm_gas_levels != 0 && cfg.amrex_max_level > 0;
    if (cfg.amrex_state && cfg.amrex_pm_max_grid_size == 0)
      c.amrex_max_grid_size = cfg.pm_refine * cfg.amrex_max_grid_size;
    return c;
  }
};

// DM state and PM work arrays (radhydro_global.f90 "DM" arrays, density,
// The amax/arms/vmax/vrms_dm globals live in SimClock.
// A star particle of the star-formation sink (docs/stage_f.md section 2;
// pm_backend = amrex with amrex_state = 1): id > Npart, x and v in the
// particle's (DM grid) units, mass in mean-matter x level-0 gas-cell
// volumes, a_form the scale factor of its formation.
struct Star {
  std::int64_t id = 0;
  float x[3] = {0, 0, 0}, v[3] = {0, 0, 0};
  float mass = 0, a_form = 0;
};

struct PMState {
  PMConfig pc;
  // domain_tools.f90 gas domains (init_domain_gas in hydro.hpp). With one
  // thread set_domain_gas hands out the domains in the order 1..Ndomain on
  // every pass, so the PM loops over dom.box in that order.
  Domains dom;
  std::int64_t dm_cpu[2] = {0, 0};      // dm_cpu(1:2,1)
  std::vector<float> xv;                // xv_dm(6,Npart)
  std::vector<std::int64_t> ll;         // ll_dm(Npart), index 1..Npart
  std::vector<std::int64_t> hoc;        // hoc_dm(Nm,Nm,Nm)
  std::vector<std::int64_t> hoc_cpu;    // hoc_cpu(Nm,Nm,1)
  std::vector<float> density;           // density(Nm+2,Nm,Nm)
  // pm_refine > 1: the mesh force (fx,fy,fz per mesh cell, the float values
  // of calc_kick's force function), filled at the start of the kick, so the
  // kick reads it instead of re-evaluating the gradient for every particle
  // and gas cell. Identical values; empty and unused with pm_refine = 1.
  std::vector<float> fmesh;
  bool fmesh_valid = false;
  // The AMReX level of pm_backend = amrex (src/amrex/pm_level.cpp), built on
  // the first particle_mesh; null otherwise. It holds no state between
  // steps: xv above is the state on every rank.
  std::shared_ptr<AmrexPM> amrex;
  // pm_backend = amrex: gather the potential into density after the step
  // (pm_driver writes it); off by default, the host code always has it.
  bool gather_potential = false;
  // The star particles (rank 0's mirror, refreshed by
  // amrex_pm_gather_stars; every rank's after a checkpoint read); the
  // running count of stars ever made (ids = Npart + 1 .. Npart + nstar)
  std::vector<Star> stars;
  std::int64_t nstar = 0;

  float* xv_p(std::int64_t ip) { return &xv[(ip - 1) * 6]; }
  const float* xv_p(std::int64_t ip) const { return &xv[(ip - 1) * 6]; }
  std::int64_t gcell(int i, int j, int k) const {  // gas cell
    const std::int64_t N = pc.Ngas;
    return (i - 1) + N * ((j - 1) + N * (k - 1));
  }
  std::int64_t mcell(int i, int j, int k) const {  // mesh cell (hoc)
    const std::int64_t N = pc.Nm;
    return (i - 1) + N * ((j - 1) + N * (k - 1));
  }
  std::int64_t dcell(int i, int j, int k) const {  // density(i,j,k), mesh
    const std::int64_t N = pc.Nm;
    return (i - 1) + (N + 2) * ((j - 1) + N * (k - 1));
  }
};

// pm_tools.f90: init_pm. Zeroes all arrays; takes the gas domains of
// init_domain_gas (scaled to the mesh when pm_refine > 1). host_arrays =
// false leaves xv, ll, hoc, hoc_cpu and density unallocated: the ranks
// other than 0 with amrex_state = 1, whose particles and mesh live on the
// PM level (docs/stage_e_state.md), would free them right away.
void init_pm(PMState& s, const PMConfig& pc, const Domains& dom, bool host_arrays = true);

// pm_tools.f90 routines. u_gas is (5,Ngas^3), v_gas is (3,Ngas^3).
void calc_drift(PMState& s, double dt_sim);
void link_dm(PMState& s);
void cic_mass_assignment(PMState& s, const std::array<int, 6>& box);
// cic_mass_assignment over all domains in order 1..Ndomain (parallel,
// bit-identical to the serial loop)
void cic_mass_assignment_all(PMState& s);
void calc_density_field(PMState& s, const float* u_gas);
void calc_gravity_field(PMState& s, double scalefactor);
// Sets clk.amax_dm, arms_dm, vmax_dm, vrms_dm (PM mesh units; the gas grid
// with pm_refine = 1).
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

// The mesh force interpolated to a particle at xv (DM-grid units) with the
// CIC weights of calc_kick, in mesh units: what calc_kick adds to the
// particle's velocity times x_dm_m/dt_sim. For tests (docs/pm_refine.md).
void force_at(const PMState& s, const float* xv, double F[3]);

// pm_tools.f90: calc_power_spectrum. Fills pk_dm(5,Nm) and pk_gas(5,Ngas)
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
