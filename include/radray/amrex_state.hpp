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
#include <functional>
#include <string>
#include <vector>

#include "radray/config.hpp"

namespace radray {

struct Simulation;

// ---- streaming the level to rank 0 (docs/rank0_streaming.md) ---------------
// A run of consecutive k-planes of the level on rank 0: planes k0 .. k0 + nk
// - 1 (0-based), every i and j, in the serial cell order (i fastest, then j,
// then k), ncomp floats per cell with the components fastest:
// data[(i + N j + N^2 (k - k0)) * ncomp + m]. The cells of the slab are the
// cells c0() .. c0() + ncell() - 1 of the serial grid.
struct StreamSlab {
  int k0 = 0, nk = 0, N = 0, ncomp = 0;
  const float* data = nullptr;
  std::int64_t ncell() const { return static_cast<std::int64_t>(nk) * N * N; }
  std::int64_t c0() const { return static_cast<std::int64_t>(k0) * N * N; }
};
using StreamFn = std::function<void(const StreamSlab&)>;

// The components c0 .. c0 + ncomp - 1 of the RT level's cells (the float
// layout of rrk::CellData: rhog, fHI, fHeI, fHeII, temp, Luv, sfr, dx(3),
// gamma(3), lambda(3), nu(5), ngam(5), fgam(5), tau(5), u(5), w(4), e; the
// first 36 are the rt_grid record) onto rank 0, slab by slab in k order:
// fn(slab) on rank 0 only. Collective. Rank 0's transient is one slab
// (gs::kSlabPlanes planes). Makes the level if it does not exist yet.
void amrex_state_stream_rt(Simulation& s, int c0, int ncomp, const StreamFn& fn);
// The same for any list of components, streamed in the order listed
void amrex_state_stream_rt(Simulation& s, const std::vector<int>& comps, const StreamFn& fn);
// The same for the gas state's components (gas_state.hpp: u(5) at 0, e at
// 5, w(4) at 6)
void amrex_state_stream_gas(Simulation& s, int c0, int ncomp, const StreamFn& fn);
constexpr int kStreamGasU = 0, kStreamGasE = 5, kStreamGasW = 6;
// The reverse (a restart, the initial conditions): rank 0 fills each slab
// in k order, fill(slab, data) writing slab.ncell() * slab.ncomp floats in
// the StreamSlab layout, and the values are copied onto the level's
// components. Collective; fill runs on rank 0 only. Rank 0's transient is
// one slab.
using FillFn = std::function<void(const StreamSlab& slab, float* data)>;
void amrex_state_fill_gas(Simulation& s, int c0, int ncomp, const FillFn& fill);
// The RT cells' components comps (rrk::CellData floats); the level is made
// if it does not exist yet
void amrex_state_fill_rt(Simulation& s, const std::vector<int>& comps, const FillFn& fill);

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
// gammaHI_avg (rt_method = m1 | otvet on rt_backend = amrex): the
// substep-averaged Gamma_HI of every cell (gamma_avg_field, rt.hpp) into
// rank 0's out, N^3 floats; out is empty on the other ranks (collective)
void amrex_gamma_avg_field(Simulation& s, std::vector<float>& out);
// The same streamed onto rank 0 slab by slab (one component; collective)
void amrex_gamma_avg_stream(Simulation& s, const StreamFn& fn);
void amrex_state_scatter_rt(Simulation& s);
// fields_levels = 1 (radhydro): one plane of the RT cells at the finest
// resolution, normal to axis (0, 1, 2 = x, y, z) through base cell
// base_index (the middle fine plane of that cell: base_index 2^L + 2^L/2 in
// finest-level cells, L = the finest level). On rank 0 `out` holds
// (kLevelSliceComps, nf, nf) floats in Fortran order, the in-plane axes in
// increasing order (x,y for a z plane), nf = Nrt 2^L: rhog, fHI, fHeI,
// fHeII, temp and the level each pixel comes from (a coarse cell fills its
// 2^(L-l) x 2^(L-l) pixels); finest = L, plane = the fine index. Collective.
constexpr int kLevelSliceComps = 6;
void amrex_rt_level_slice(Simulation& s, int axis, int base_index, std::vector<float>& out, int& nf, int& finest,
                          int& plane);

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
