// The hydro backend switch (docs/stage_e_hydro.md): hydro_backend = serial
// (the unsplit step on the host arrays, src/hydro_unsplit.cpp) | amrex (the
// same step on AMReX MultiFabs, src/amrex/hydro_level.cpp, MPI + OpenMP;
// hydro_scheme = unsplit only). This header is AMReX-free; with
// RADRAY_AMREX off the amrex_hydro_* functions are stubs that throw
// (src/amrex/hydro_level_stub.cpp). The process functions (backend_init,
// io_rank, ...) are those of rt_backend.hpp: one AMReX/MPI process set-up
// serves both backends.
//
// Stage E1 model: every rank holds the full host state (GasGrid); a step
// copies the rank's boxes in from it, sub-cycles the unsplit kernel on the
// level, and gathers the result back onto every rank (MPI_Allgatherv), so
// the host state stays complete on every rank until the PM, the halo finder
// and the I/O are distributed too.
#pragma once

#include <string>

#include "radray/config.hpp"
#include "radray/cosmology.hpp"
#include "radray/grid.hpp"
#include "radray/hydro.hpp"

namespace radray {

// Parameters: hydro_backend, amrex_hydro_max_grid_size; false for any other
// key.
bool set_hydro_backend_param(Config& c, const std::string& key, const std::string& val);
void copy_hydro_backend_param(const Config& from, Config& to);
std::string hydro_backend_summary(const Config& c);
// Throws for a combination the backend does not implement (hydro_scheme =
// legacy) or a build without AMReX.
void check_hydro_backend_param(const Config& c);
inline bool use_amrex_hydro(const Config& c) { return c.hydro_backend == HydroBackend::AMReX; }

// Advances the gas by dt (code time) with the unsplit kernel on the level,
// sub-cycled on its CFL condition exactly as hydro_step_unsplit, and
// updates hs.ucounters; called by hydro_step_unsplit after fold_frame.
// gas.u/w/e are the complete host arrays on every rank, in and out.
void amrex_hydro_step(const Config& cfg, const Units& un, double dt, GasGrid& gas, HydroState& hs);
// One line of the level's phase timers ("" without a level)
std::string amrex_hydro_timers();
// Frees the level (backend_finalize calls it before amrex::Finalize)
void amrex_hydro_release();

}  // namespace radray
