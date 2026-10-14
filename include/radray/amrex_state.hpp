// The state on the levels (docs/stage_e_state.md): amrex_state = 1 makes the
// AMReX levels the owners of the gas state between steps. The gas (u, w, e)
// lives in one MultiFab on a decomposition shared by the RT, hydro and PM
// levels, which read and write it locally instead of copying the full host
// grid in and gathering it back onto every rank; the host arrays of the
// gas grid and the RT grid exist only on rank 0, as a mirror refreshed
// from the levels for the statistics, the file outputs and the checkpoint
// (the file formats are unchanged), and nothing of the grid is allocated
// on the other ranks beyond their boxes.
//
// This header is AMReX-free; with RADRAY_AMREX off the amrex_state_*
// functions are stubs that throw (src/amrex/gas_state_stub.cpp). The
// process functions (backend_init, io_rank, ...) are those of
// rt_backend.hpp. Requires rt_backend = hydro_backend = pm_backend = amrex,
// hydro_scheme = unsplit, rt_method = rays, recomb_rad = 0, N_rt_gas = 1
// (check_amrex_state_param).
#pragma once
#include <cstdint>
#include <string>
#include <vector>

#include "radray/config.hpp"

namespace radray {

struct Simulation;
namespace pm {
struct PMState;
}

// Parameters: amrex_state (0 | 1), amrex_balance (0 | 1),
// amrex_balance_ray_weight, amrex_balance_min_gain, amrex_balance_method
// (sfc | knapsack); false for any other key.
bool set_amrex_state_param(Config& c, const std::string& key, const std::string& val);
void copy_amrex_state_param(const Config& from, Config& to);
std::string amrex_state_summary(const Config& c);
// Throws for a combination the state on the levels does not support or a
// build without AMReX.
void check_amrex_state_param(const Config& c);
inline bool use_amrex_state(const Config& c) { return c.amrex_state; }

// Builds the shared decomposition and the gas MultiFab (collective; once).
// The levels built afterwards use its BoxArray and DistributionMapping.
void amrex_state_init(Simulation& s);
// Rank 0's host gas arrays (gas.u, gas.w, gas.e) onto the level, slab by
// slab (collective): the initial conditions and a restart.
void amrex_state_scatter_gas(Simulation& s);
// The level's gas into rank 0's host arrays, slab by slab (collective):
// before the statistics, the file outputs and a checkpoint. The other
// ranks' host arrays stay empty.
void amrex_state_gather_gas(Simulation& s);
// summarize_gas's per-cell clamps (u(1) >= Dgas_min, w >= 0) on the level
// (collective; the serial function then finds nothing to clamp on rank 0's
// mirror, so its sums are bitwise those of the host run).
void amrex_state_clamp_gas(Simulation& s);
// add_halo_metals on the level: w(4) of cell c += dz[n] (in double, rounded
// to float once per entry, in list order) for every entry whose cell the
// rank owns. cells are 0-based linear gas cell indices.
void amrex_state_add_metals(Simulation& s, const std::vector<std::int64_t>& cells,
                            const std::vector<double>& dz);

// The RT level's cells (rt_level.cpp): all RT grid fields and the gas into
// rank 0's s.rt and s.gas (collective); stats_only: just rhog, fHI, fHeI,
// fHeII, temp (write_stats of every step). And the reverse (a restart):
// rank 0's s.rt fields onto the level's cells.
void amrex_state_gather_rt(Simulation& s, bool stats_only);
// A movie frame's fields only (rhog, fHI, fHeI, fHeII, temp, gamma(1)) into
// rank 0's s.rt (collective)
void amrex_state_gather_rt_movie(Simulation& s);
void amrex_state_scatter_rt(Simulation& s);

// amrex_balance = 1: rebalances the shared decomposition on the current
// per-box cost (cells + amrex_balance_ray_weight * live rays) before an RT
// step, moving the gas, the RT cells and rays and the particles when the
// efficiency (mean over max rank cost) would gain more than
// amrex_balance_min_gain. Collective. Returns whether it remapped; eff_old
// and eff_new are the efficiencies before and after (the proposed one when
// not remapped). The results are bitwise independent of the mapping.
bool amrex_balance_step(Simulation& s, pm::PMState& pm, double& eff_old, double& eff_new);

// Releases the state (backend_finalize)
void amrex_state_release();

}  // namespace radray
