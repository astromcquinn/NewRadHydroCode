// The unsplit hydro scheme (hydro_scheme = unsplit, docs/hydro_unsplit.md):
// host driver of the MUSCL-Hancock + HLLC kernel with the dual-energy
// formulation (include/radray/kernels/hydro_mh.hpp) on the uniform periodic
// gas grid, and the routing of the radhydro loop's two half steps to either
// scheme.
//
// State: GasGrid::u (float: rho, rho v, E in the LAB frame, i.e. the grid
// velocity HydroState::v is zero), GasGrid::w (species) and GasGrid::e (the
// internal energy density of the dual-energy formulation). A step converts
// the float state to double, sub-cycles the unsplit kernel on its CFL
// condition over the requested interval and rounds back to float once.
//
// Loop placement (docs/hydro_unsplit.md section 3): the original's six
// directional sweeps around particle_mesh become two unsplit steps over
// dt_sim/2 around the PM kick (Strang splitting with the gravity kick in
// the middle, the same call sites and the same calc_timestep); the PM kick
// of the gas acts on the conserved variables (pm::particle_mesh_unsplit).
#pragma once

#include <string>

#include "radray/config.hpp"
#include "radray/cosmology.hpp"
#include "radray/grid.hpp"
#include "radray/hydro.hpp"

namespace radray {

// params.txt keys hydro_scheme (legacy | unsplit), unsplit_cfl, dual_energy,
// dual_eta1, dual_eta2 (Config). Returns false for any other key. Throws on
// a bad value.
bool set_hydro_scheme_param(Config& cfg, const std::string& key, const std::string& val);
// Copies those fields (the drivers parse into a scratch Config first).
void copy_hydro_scheme(const Config& from, Config& to);
// One line for the drivers' banners
std::string hydro_scheme_summary(const Config& cfg);

// Prepares the state for the unsplit scheme (no-op with the legacy scheme):
// folds a non-zero grid velocity into u (lab frame) and zeroes it, and
// allocates and fills gas.e = max(E - KE, rho Tmin) if it is empty (initial
// conditions, or a restart from a legacy checkpoint). Call after the ICs or
// a checkpoint are read and calc_sim_var has run.
void init_unsplit_state(const Config& cfg, const Units& un, GasGrid& gas, HydroState& hs);

// Advances the gas by dt (code time) with the unsplit kernel, sub-cycled
// on dt_cfl = unsplit_cfl / max_cells(sum_d |v_d| + 3 c) (dx = 1 cell).
// Updates hs.ucounters and hs.dt_hy. Requires gas.e (init_unsplit_state).
// scalefactor: only the AMReX backend's star-formation sink uses it
void hydro_step_unsplit(const Config& cfg, const Units& un, double dt, GasGrid& gas,
                        HydroState& hs, double scalefactor = 1.0);

// The two halves of the hydro part of radhydro.f90's radhydro, for either
// scheme. Legacy: First = transposes + three forward sweeps, Second = three
// reverse sweeps + transposes (bit-identical to the original sequence).
// Unsplit: each half is hydro_step_unsplit over clk.dt_sim/2.
enum class HydroHalf { First, Second };
void hydro_half_step(const Config& cfg, const Units& un, const SimClock& clk, HydroHalf half,
                     GasGrid& gas, HydroState& hs);

}  // namespace radray
