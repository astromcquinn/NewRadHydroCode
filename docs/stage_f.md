# Stage F: AMR for the gas, then the RT (`amrex_max_level`)

*Started 2026-10-14 (branch `stage-f-decisions`). The user's decisions are
in `docs/PLAN.md`, "Stage F decisions": one coarse RT substep across
levels, B17 under `legacy_physics = 0`, cooling on in dense gas with an
aggressive star-formation sink, the rest the implementer's choice, in the
order F1 (hierarchy, refinement off, bitwise), F2 (refined hydro, the
cooling flag and the sink), F3 (the RT across levels).*

## 0. Summary

| Item | State |
|---|---|
| F1: the gas hierarchy (`amrex_max_level`, AmrMesh regridding, FillPatch, refluxing, average-down) in `src/amrex/hydro_level.cpp`; `amrex_max_level = 0` is bitwise the Stage E level | in progress |
| F2: refined hydro validated on the hydro test suite; `cool_suppress_overdensity`; the star-formation sink | planned |
| F3: rays across levels with the coarse substep (`rt_cells_per_substep` = 2^l on level l) | planned |

## 1. The hierarchy (F1)

**Parameters** (`params.txt` / `--set`, AMReX-free in `hydro_backend.hpp`;
all backends refuse them unless `hydro_backend = amrex`):

| key | default | meaning |
|---|---|---|
| `amrex_max_level` | 0 | refinement levels of the gas hierarchy (0 = Stage E, one level) |
| `amrex_refine_mass` | 8 | a level-l cell is tagged when its gas mass exceeds `amrex_refine_mass` × the mean mass of a base cell, i.e. ρ/ρ̄_gas > `amrex_refine_mass` · 8^l (the Lagrangian criterion of `docs/ROADMAP.md` §5) |
| `amrex_regrid_every` | 1 | regrid every n hydro calls (two per `radhydro` step) |
| `amrex_fine_max_grid_size` | 32 | box size of the refined levels (level 0 keeps `amrex_hydro_max_grid_size` or the shared decomposition) |
| `amrex_blocking_factor` | 8 | AMReX `blocking_factor` (grids are multiples of it) |
| `amrex_n_error_buf` | 1 | cells of buffer around tagged cells |
| `amrex_grid_eff` | 0.7 | AMReX `grid_eff` |

**Data.** Level 0 is the Stage E level: with `amrex_state = 1` its
BoxArray and DistributionMapping are the shared decomposition and the gas
state (float, `gas_state.hpp`) owns the level-0 data between steps; the
refined levels (`U` in double, their own BoxArray from the tags and their
own SFC mapping) are owned by the hydro hierarchy between steps and are
written to the checkpoint as AMReX `VisMF` files next to it
(`<checkpoint>.amr/Level_l/`). The refinement ratio is 2 in every
direction; the geometry of level l is `refine(geom[l-1], 2)`, periodic.

**The step** (`advance(lev, t, dt)`, recursive; `dt` is the level-0 CFL
substep `dt_sub` of the serial driver, and because dx and dt halve
together the kernel's `dtdx` is `dt_sub` on every level):

1. ghosts of `U[lev]` (level 0: `FillBoundary`; finer: `FillPatchTwoLevels`
   from the parent's old and new states, interpolated in time, with
   `cell_cons_interp`, plus the level's own data), the dual-energy
   selection into `Es` (ghosts: `FillBoundary`, or the parent's current
   `Es` interpolated);
2. the substep of `hydro_mh.hpp` on the level (predictor, faces into three
   face-centred flux MultiFabs, the update, `En` ghosts, the floors and
   the dual-energy finalisation), as in Stage E1;
3. if `lev < finest`: the level's fluxes into the flux register of
   `lev+1` (`CrseInit`, factor −dtdx), two `advance(lev+1, ·, dt/2)`
   (each adds its fluxes with factor +dtdx/8: dt_f·A_f/V_c), then
   `Reflux` into `U[lev]` (scale 1) and `average_down` of `U[lev+1]` onto
   `U[lev]` (all 10 components: conserved variables, e, w).

The CFL step is the serial driver's arithmetic with the signal maximum
taken over every level (a level-l cell with signal s needs
dt/2^l ≤ cfl·2^-l/s, the same bound as a base cell), recomputed before
every substep as the serial driver does. With `amrex_max_level = 0` the
recursion is the Stage E1 loop, operation for operation: bitwise
identical (`ctest -L amrex_hydro`).

**Regridding.** `AmrMesh::MakeNewGrids` with `ErrorEst` tagging on the
level's density (above); new levels are filled by `FillPatchTwoLevels`
(the old fine data where it exists, the parent interpolated elsewhere),
at the start of a hydro call every `amrex_regrid_every` calls and at the
first call. Level 0's grids never change here (the balancing of Stage E4
moves them; the refined levels are then rebuilt at the next regrid).

**Coupling to the RT and PM before F3 and Stage G** (interim, approximate,
documented so it can be removed): the RT and the PM act on level 0 only.
Between hydro calls they change the level-0 gas (w, u(5), e by the RT;
the momenta and E by the PM kick; the metals). At the start of a hydro
call the hierarchy compares the level-0 state with what it wrote out last
(a float copy) and pushes every change down to the fine cells of that
coarse cell: the density and the species by the coarse ratio (additive
share where the old value is 0), the momenta by the coarse velocity
change (p_f += ρ_f Δv_c), E by the fine cell's kinetic-energy change plus
its thermal change, e by the coarse ratio. Unchanged coarse cells leave
their fine cells untouched exactly. After the step the fine levels are
averaged down, so level 0 is always consistent with them.

## 2. Refined hydro, cooling and the sink (F2)

Validation: the hydro test suite (`tests/hydro_suite`, `docs/hydro_tests.md`)
with `hydro_backend = amrex` and `amrex_max_level = 1, 2` against the
analytic solutions and the uniform runs at the matched finest resolution
(the composite of the finest data on a uniform grid of N·2^L cells,
`amrex_hydro_fine_grid`, test problems only); Sod, Sedov, Zel'dovich
pancake, Kelvin–Helmholtz. Plus the full-loop smoke problem with
refinement for restarts and thread/rank independence.

`cool_suppress_overdensity` (done with F1: `params.txt` / `--set`, both
drivers, `make_production.py --cool-suppress-overdensity`; default 200 =
the original: the cooling rate × (200/x)² and the SN heating × (200/x)^1.5
above x = ρ/ρ̄; 0 = no suppression; `rrk::PhysCtx::cool_suppress` in
`rt_cell.hpp`; tests `cool_suppress_default_identity` (200 is bitwise the
default) and `cool_suppress_off_smoke`). The sink (`sf_overdensity`, 0 =
off; `sf_efficiency` per free-fall time): gas above the threshold is
converted to collisionless star particles in the PM container, with
their mass.

**The sink (F2b, built 2026-10-14; tests label `sink`).** In every hydro
call, on every level's cells not covered by a finer level, a cell with
x = ρ_gas/ρ̄_gas > `sf_overdensity` loses
the mass dm = min(1, `sf_efficiency` · dt/t_ff) · ρ V, with
t_ff = √(3π/32Gρ_tot) and ρ_tot = ρ_gas Ω_m/Ω_b (the local total density
of the PM mesh when it is available), at constant specific quantities:
ρ, the momenta, the species densities and e scale by (1 − dm/ρV), E
loses dm (½v² + e/ρ)/V. The removed mass becomes one star particle per
cell and call at the cell centre with the cell's velocity, its mass in
units of the dark-matter particle mass, its formation time, and an id
above the dark-matter range. Stars live in the PM container
(`pm_backend = amrex`, `amrex_state = 1`: production) with a mass
component and a type flag; the CIC deposit weights particles by mass
(dark matter 1), the kick and drift treat them as dark matter, the
load balancing moves them with their boxes. They are not in rank 0's
`xv` mirror: the FoF and the `xv_dm` files see dark matter only (the
halo masses miss the stellar mass; documented), and the stars are
gathered to rank 0 for `DM/stars.<ztag>` (int64 n, then records
{int64 id; float x[3], v[3], mass, a_form}, sorted by formation time and
position) and `sink.txt` (the count and the total mass) at every data
output and at the end, and for the checkpoint (version 8, loaded at a
restart). Ids are `Npart + 1 + the running count`, numbered by a global
sort of each call's event positions, so they do not depend on the rank
count (the stars file and the checkpoint are byte-identical for 1, 2 and
3 ranks: `amrex_sink_ranks_smoke`; a restart is bit-identical:
`amrex_sink_restart_smoke`; `sink_check.py` checks the bookkeeping). The
hydro drivers without a PM (`hydro_driver`, the suite) apply the sink as
a mass loss (tests of the hydro part). The default `sf_overdensity = 0`
keeps every stored reference; `make_production.py --sf-overdensity` (and
`--amrex-max-level`, `--amrex-refine-mass`) write the production
parameters; the production value is a few thousand (user decision), to
be set with the first refined production run. Code: `hydro_level.cpp`
`sink`, `pm_level.cpp` `add_stars` / `gather_stars` / `load_stars`,
`pm_backend.hpp`, `radhydro.cpp` `write_stars` / `write_sink_summary`,
`checkpoint.cpp`.

## 3. The RT across levels (F3)

Rays live on the level of their cell and cross 2^l fine cells per coarse
substep in level-l patches; deposits per level; the exposure deposit and
`photon_conserve`; the angular resolution per level; no time
sub-cycling of the transport (user decision). Design when F2 is done.

## 4. Files and tests

| File | What |
|---|---|
| `include/radray/hydro_backend.hpp`, `src/hydro_backend.cpp` | the hierarchy's parameters (§1) and checks, `amrex_hydro_finest_level`, `amrex_hydro_write_levels` / `read_levels` (AMReX-free) |
| `src/amrex/hydro_level.cpp` | the hierarchy: `LevelData`, `Mesh` (AmrMesh, `ErrorEst` = `tag`), `regrid`, `fill_U` / `select` / `substep` per level, `advance` (recursive sub-cycling, flux register, average-down), `push_down`, `remap`, the VisMF files |
| `src/amrex/gas_state.hpp/.cpp` | the balancing calls `amrex_hydro_remap` (the refined levels survive a remap) |
| `apps/radhydro.cpp` | `<checkpoint>.amr/` written with every checkpoint and read at a restart; removed with its checkpoint |
| `include/radray/config.hpp`, `src/io.cpp`, `apps/rt_driver.cpp`, `apps/radhydro.cpp`, `include/radray/kernels/rt_cell.hpp` | `cool_suppress_overdensity` (§2) |
| `tests/hydro_suite/suite_driver.cpp`, `run_suite.py` | the backend keys and `backend_init`; `--amr L --refine-mass M` |
| `tools/hyak/hydro_suite_amr.slurm`, `tools/analysis/hydro_suite_amr.py` | the suite with refinement on klone and its tables |
| `tests/amr_conservation_check.py`, `tests/run_amrex.cmake` (`STRICT=cons`) | the conserved totals of a refined run against the uniform run (to 1e-7: the driver's float fields; the hierarchy's double level-0 mass is `mass0` in its timers line) |
| `tools/production/make_production.py` | `--cool-suppress-overdensity` |

Tests: `ctest -L amrex_hydro` (27, refinement off: bitwise the serial
step), `ctest -L amr` (the label matches `amrex` too; the new tests are
`amrex_hydro_amr_ref_<p>_L<L>` (conservation, the level count) and
`amrex_hydro_amr_bitwise_<p>_L<L>_r<np>t<nt>_mgs<m>` (layout independence)
for sod_x (1 level, refine_mass 0.5), sound_x (1, 0.99) and random3d (2
levels, 8)), `cool_suppress_default_identity`, `cool_suppress_off_smoke`
(label physics).
