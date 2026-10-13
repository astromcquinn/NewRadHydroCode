// The PM backend switch (docs/stage_e_pm.md): pm_backend = serial (the
// original's linked lists, float CIC and FFTW on the host, src/pm.cpp) |
// amrex (the dark matter in an AMReX ParticleContainer, int64 fixed-point
// CIC deposits, the Poisson solve on an AMReX level of the PM mesh and the
// kicks of the particles and the gas on that level, src/amrex/pm_level.cpp,
// MPI + OpenMP). This header is AMReX-free; with RADRAY_AMREX off the
// amrex_pm_* functions are stubs that throw (src/amrex/pm_level_stub.cpp).
// The process functions (backend_init, io_rank, ...) are those of
// rt_backend.hpp: one amrex::Initialize serves both backends.
#pragma once

#include <string>

#include "radray/config.hpp"
#include "radray/cosmology.hpp"

namespace radray {

namespace pm {
struct PMState;
}

// Parameters: pm_backend, amrex_pm_max_grid_size, amrex_pm_fft; false for
// any other key.
bool set_pm_backend_param(Config& c, const std::string& key, const std::string& val);
void copy_pm_backend_param(const Config& from, Config& to);
std::string pm_backend_summary(const Config& c);
// Throws for a combination the backend does not implement or a build
// without AMReX.
void check_pm_backend_param(const Config& c);
inline bool use_amrex_pm(const Config& c) { return c.pm_backend == PMBackend::AMReX; }

// One PM step on the level, in place of the host's link_dm /
// calc_density_field / calc_gravity_field / calc_kick (the two drifts of
// particle_mesh stay on the host, before and after): the particles of s.xv
// (after the first drift) and the gas of u_gas (and v_gas or e_gas) are
// copied onto the level at the start, the kicked velocities and gas are
// gathered back onto every rank at the end. s.density is left as the
// potential on the mesh (pm_driver writes it), as the host code does.
// calc_kick's amax/arms/vmax/vrms go into clk as on the host.
void amrex_pm_step(pm::PMState& s, float* u_gas, float* v_gas, SimClock& clk, const Units& u);
void amrex_pm_step_unsplit(pm::PMState& s, float* u_gas, const float* e_gas, double eta1,
                           SimClock& clk, const Units& u);
// The fields alone (tests): the density of s.xv (no drift) and u_gas, the
// potential and the mesh force on the level, gathered into s.density (the
// potential, as calc_gravity_field leaves it) and s.fmesh (fmesh_valid =
// true), so the host's force_at and gas_force read the backend's mesh.
void amrex_pm_fields(pm::PMState& s, const float* u_gas, double scalefactor);
// One line of the backend's phase timers ("" without the level)
std::string amrex_pm_timers(const pm::PMState& s);

}  // namespace radray
