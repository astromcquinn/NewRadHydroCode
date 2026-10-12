// On-the-fly friends-of-friends (FoF) halo finder for the dark-matter
// particles, with optional spherical-overdensity (SO) refinement.
//
// The original RadHydro does not find halos itself: read_halo
// (io_tools.f90) reads catalogs produced by an external halo finder. This
// module replaces that external step. See docs/halos_and_sources.md.
//
// Particle layout: xv_dm(6,Npart) real(4), Fortran order (component
// fastest), i.e. particle p has x,y,z at xv[6p+0..2] and vx,vy,vz at
// xv[6p+3..5]. In the original code (pm_tools.f90: calc_drift wraps with
// mod(...,ndmr)) the positions are in DM-grid units [0,Ndm) and the
// velocities in DM-grid units per time_unit (summarize_dm converts to km/s
// with x_dm_gas*vel_unit/1E5). The finder itself only needs a periodic box
// length L in the same units as the positions, so it works for any unit.
//
// Determinism: the FoF partition (connected components of the linking
// graph) is unique; groups are labelled by their smallest particle index
// and members are kept in ascending particle order, so every derived
// quantity is summed in a fixed order. The output is bit-identical for any
// number of OpenMP threads.
#pragma once

#include <cstdint>
#include <vector>

#include "radray/config.hpp"
#include "radray/cosmology.hpp"

namespace radray {

// ---------------------------------------------------------------------------
// Low level: FoF partition
// ---------------------------------------------------------------------------

struct FofGroups {
  // label[p] = smallest particle index in p's group (a particle that links
  // to nobody is its own group).
  std::vector<std::int64_t> label;
  // Groups with at least nmin members, ordered by their label (ascending).
  // members[offset[g] .. offset[g+1]) are the member indices, ascending.
  std::vector<std::int64_t> group_label;
  std::vector<std::int64_t> offset;
  std::vector<std::int64_t> members;
  std::int64_t ngroups() const { return static_cast<std::int64_t>(group_label.size()); }
  std::int64_t size(std::int64_t g) const { return offset[g + 1] - offset[g]; }
};

// Friends-of-friends with an absolute linking length `ll` in a periodic
// cube [0,L)^3. Two particles are friends when their minimum-image
// separation satisfies d^2 <= ll^2 (computed in double precision).
// `pos` points at x of particle 0; consecutive particles are `stride`
// floats apart (6 for xv_dm). Positions outside [0,L) are wrapped.
// `ncell_max` caps the linked-cell grid per dimension (0 = automatic).
FofGroups fof_partition(const float* pos, std::int64_t npart, int stride, double L,
                        double ll, std::int64_t nmin, int ncell_max = 0);

// Reference O(N^2) implementation (tests only).
FofGroups fof_partition_bruteforce(const float* pos, std::int64_t npart, int stride,
                                   double L, double ll, std::int64_t nmin);

// ---------------------------------------------------------------------------
// Halo catalog
// ---------------------------------------------------------------------------

// Conversion from particle units to the units of the catalog, and the
// cosmological quantities needed for virial definitions.
struct HaloUnits {
  double L = 0;           // box length in particle position units (Ndm)
  double box = 0;         // box length in comoving Mpc/h
  double pos_to_mpch = 0; // particle position unit -> comoving Mpc/h
  double vel_to_kms = 0;  // particle velocity unit -> peculiar km/s
  double mpart_msun = 0;  // particle mass in Msun (Msun_unit)
  double redshift = 0;
  double hubble0 = 0;     // little h
  double omegam_z = 0;    // Omega_m(z)
  double rho_crit = 0;    // proper critical density at z, g/cm^3
  double rhom_mean = 0;   // proper mean matter density at z, g/cm^3
};

// HaloUnits for the running simulation. Needs calc_sim_var to have been
// called for the current redshift. Particle mass = Msun_unit of
// cosmo_tools.f90 calc_sim_var: rhoc*omegam*(box Mpc)^3/h/Npart/Msun with
// Npart = Ndm^3 (it carries ALL the matter, Omega_m, not Omega_m-Omega_b).
HaloUnits make_halo_units(const Config& cfg, const SimClock& clk, const Units& u);

// Particle mass in Msun for Npart particles in a box of `box` Mpc/h
// (cosmo_tools.f90 calc_sim_var: Msun_unit = mp_unit/Msun_cgs).
double particle_mass_msun(std::int64_t npart_total, double box);

// Bryan & Norman (1998) virial overdensity relative to the CRITICAL
// density, flat universe: x = Omega_m(z)-1, Delta_c = 18pi^2 + 82x - 39x^2.
double delta_vir_bryan_norman(double omegam_z);

enum class HaloCentre {
  CentreOfMass,       // periodic-safe centre of mass of the FoF members
  PotentialMinimum,   // member with the lowest softened pairwise potential
                      // (O(N^2) per halo; for N > potmin_nmax the
                      // potential is summed over every k-th member only,
                      // k = ceil(N/potmin_nmax), see docs)
};

struct HaloFinderParams {
  double b = 0.2;                 // linking length / mean interparticle spacing
  std::int64_t min_npart = 20;    // minimum FoF members for a halo
  HaloCentre centre = HaloCentre::CentreOfMass;
  std::int64_t potmin_nmax = 5000;

  // Spherical-overdensity refinement around the FoF centre.
  bool so_refine = false;
  // SO threshold relative to rho_crit(z); <= 0 -> Bryan & Norman Delta_c.
  double so_delta_crit = 0;
  int so_iterations = 3;          // re-centring iterations (centre of mass
                                  // of the particles inside R_SO)

  // Mean molecular weight for T_vir (0.59: fully ionised primordial gas).
  double mu = 0.59;
  int ncell_max = 0;              // FoF grid cap per dimension (0 = auto)
};

// One halo in physical units. x in comoving Mpc/h (same unit as the
// legacy catalog's source_type%x, see halo_tools.f90 link_halo: x/x_rt_box),
// v in peculiar km/s, M in Msun, R (virial or SO radius) in comoving Mpc/h,
// T (virial temperature) in K.
struct Halo {
  std::int64_t id = 0;      // label of the FoF group (smallest member index)
  std::int64_t N = 0;       // number of particles (FoF, or inside R_SO)
  std::int64_t N_fof = 0;
  double M = 0;             // Msun
  double x[3] = {0, 0, 0};  // comoving Mpc/h, in [0, box)
  double v[3] = {0, 0, 0};  // km/s
  double R = 0;             // comoving Mpc/h
  double T = 0;             // K
};

// fof_dt_myr > 0 (docs/pm_refine.md, "FoF cadence"): the halos of the last
// FoF, reused by the on-the-fly sources until fof_dt_myr of cosmic time has
// passed. Masses, R and T stay those of the FoF; positions are drifted with
// the halo mean velocity (in gas-grid units at the FoF, vel_unit) over the
// code time dtau accumulated since. Checkpointed (version 5).
struct FofCache {
  std::vector<Halo> halos;  // as found by the last FoF
  bool valid = false;
  double year = 0;          // clk.year of the last FoF
  double vel_unit = 0;      // Units::vel_unit at the last FoF
  double dtau = 0;          // code time (sum of dt_sim) since the last FoF
  std::int64_t nfof = 0, nreuse = 0;  // FoF calls and reuses (log only)
};

// Find halos; output sorted by descending M, ties by descending N, then by
// ascending id (deterministic).
std::vector<Halo> find_halos(const float* xv, std::int64_t npart, const HaloUnits& hu,
                             const HaloFinderParams& p);

// Virial radius (comoving Mpc/h) and temperature (K) of a halo of mass M
// (Msun) at the redshift of `hu`, for overdensity delta_crit relative to
// rho_crit(z).
double virial_radius_mpch(double M_msun, double delta_crit, const HaloUnits& hu);
double virial_temperature(double M_msun, double R_mpch, double mu, const HaloUnits& hu);

}  // namespace radray
