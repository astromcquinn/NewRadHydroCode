// The PM's side of the sub-cycled gravity of the refined levels
// (pm_level_subcycle = 1 with pm_level_solve, docs/stage_g.md section 6),
// called by the hydro hierarchy (hydro_level.cpp) inside a hydro call;
// implemented in pm_level.cpp. Internal to the AMReX backends.
//
// Times are within the hydro call (0 at its start, the hydro's own t);
// level l's particles are kicked, drifted and solved for in its substeps
// (kick-drift-kick), level 0's are not (one kick per base step, in the
// PM step between the two hydro calls). U[l] is the hydro's state of
// level l (amrex_state = 1: the gas state's boxes and mapping), component
// 0 the density, 1-3 the momenta, 4 the total energy, 5 the internal
// energy.
#pragma once
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>

namespace radray::pmsub {

// A hydro call of a sub-cycled run is in progress (amrex_pm_subcycle_call
// was called for it): the hydro must call the functions below
bool active();
// The start of the call, after the hydro's push-down and regrid: the PM
// levels from the gas state's, the particles onto them (at the base step's
// middle with the corrections of level 0's kick), every level solved at
// t = 0. Collective.
void open(const amrex::Vector<amrex::MultiFab*>& U, int finest);
// The end of the call (the hydro's last substep done)
void done();
// Level l's potential (and l >= 1: its force) at time t, from the
// densities of the levels' particles at t (moved there virtually where
// they are at another time) and of U's gas; level 0's is the boundary
// values of level 1 only. Collective.
void solve(int l, double t, const amrex::MultiFab& U);
// Level l's gas (its valid cells and ng ghost layers) and particles kicked
// by its force over dt. Collective only through the momentum diagnostic.
void kick(int l, amrex::MultiFab& U, int ng, double dt, double eta1, double rho_min);
// Level l's particles drifted by h from time t0 (no reassignment); their
// time is then t0 + h (the time the hydro passes to the closing solve)
void drift(int l, double t0, double h);
// The end of a level-l substep: the particles of levels l..finest
// reassigned among them (l = 0: those of levels >= 1 that left level 1
// go to level 0). Collective.
void regroup(int l);

}  // namespace radray::pmsub
