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
| F1: the gas hierarchy (`amrex_max_level`, AmrMesh regridding, FillPatch, refluxing, average-down) in `src/amrex/hydro_level.cpp`; `amrex_max_level = 0` is bitwise the Stage E level | done 2026-10-14 |
| F2: refined hydro validated on the hydro test suite (§2: one level on 32³ matches the uniform 64³ run); `cool_suppress_overdensity`; the star-formation sink (F2b) | done 2026-10-14 (the pancake, KH and full-tier cases pending) |
| F3: the gas state on every level (F3a) and rays across levels with the coarse substep (F3b, `rt_cells_per_substep` = 2^l on level l) | done 2026-10-15: the I-front validated inside the refined region and across its boundary (§3) |
| F4: the PM on the refined gas levels (`pm_gas_levels`, default 1): the deposit from the finest level, every level kicked on its own cells; `make_production.py --ndm` (e.g. `pm_refine = 4` with 2× particles) | done 2026-10-15 (§4; on klone: pending) |

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

**Three traps met on the way** (all fixed, all with tests): (1) the
kernel's `dtdx` is the *level-0* substep on every level, since a level-l
cell of dx 2^-l takes dt 2^-l; with the level's own dt the fine levels
advanced at half speed per level (Sedov's shock radius came out 20–40 %
short while mass, momentum and energy were still conserved, the register
factors being consistent with the wrong update). (2) The time within a
hydro call restarts at zero, so every level's interpolation interval
(the parent's old and new states) must be reset at each call; with the
previous call's interval the first fine fill of a call extrapolated and
produced negative ghost densities. (3) AMReX's `cell_cons_interp` and
`lincc_interp` limit the slopes per direction, and in a field with
density contrasts of 1e6 (`random3d`) the three-direction sum made
negative fine densities next to small coarse values; the hierarchy
interpolates with `mf_linear_slope_minmax_interp` (the slope vector
limited so that no fine value leaves the range of the coarse neighbours),
which FillPatch dispatches only when passed as an `MFInterpolater*`.
`RADRAY_AMR_DEBUG=1` scans every level after each phase and reports the
first non-finite or non-positive cell with its parent column, and counts
the bad or untouched ghost cells of every fill.

**Regridding.** `AmrMesh::MakeNewGrids` with `ErrorEst` tagging on the
level's density (above); new levels are filled by `FillPatchTwoLevels`
(the old fine data where it exists, the parent interpolated elsewhere),
at the start of a hydro call every `amrex_regrid_every` calls and at the
first call. Level 0's grids never change here (the balancing of Stage E4
moves them; the refined levels are then rebuilt at the next regrid).

**Coupling to the RT and PM before F3 and Stage G** (interim, approximate,
documented so it can be removed; the RT acts on every level since F3, the
PM since §4 with `pm_gas_levels = 1`): the RT and the PM act on level 0 only.
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

**Validation on the hydro test suite** (`tests/hydro_suite`,
`docs/hydro_tests.md`; `run_suite.py --scheme unsplit --amr L
--refine-mass 0.2`, `tools/hyak/hydro_suite_amr.slurm`, one klone node,
2026-10-14; JSON in `docs/data/stage_f/`). With `refine_mass = 0.2` every
cell above 0.2× the mean gas density is refined: the whole Sedov box (the
ambient medium is at the mean density) and everything but the
undisturbed low-density state of a tube, so a refined run on N³ is the
uniform (2^L N)³ solution averaged down to N³, and the analysis (the
level-0 fields) should sit between the uniform N³ and (2N)³ runs. L1
errors against the exact solutions (ρ, u, p for the tubes; ρ, p, v and
the shock radius for Sedov at t = 10), uniform / 1 level / 2 levels:

| case | ρ | u or v | p | R error | wall [s] |
|---|---|---|---|---|---|
| Sod x, N = 32 | 0.0227 / 0.0098 / 0.0098 | 0.0242 / 0.0099 / 0.0099 | 0.0196 / 0.0080 / 0.0080 | | 3 / 11 / 13 |
| Sod x, N = 64 | 0.0120 / 0.0057 / 0.0057 | 0.0142 / 0.0066 / 0.0066 | 0.0090 / 0.0042 / 0.0042 | | 11 / 148 / 168 |
| Sod xyz (diagonal), N = 32 | 0.0248 / 0.0189 / 0.0189 | 0.0262 / 0.0186 / 0.0186 | 0.0228 / 0.0170 / 0.0170 | | 3 / 6 / 7 |
| strong shock, N = 32 | 0.159 / 0.086 / 0.036 | 0.0257 / 0.0128 / 0.0091 | 0.0249 / 0.0104 / 0.0096 | | 1 / 11 / 50 |
| strong shock, N = 64 | 0.114 / 0.066 / 0.034 | 0.0142 / 0.0068 / 0.0054 | 0.0137 / 0.0063 / 0.0047 | | 8 / 154 / 508 |
| 123 (double rarefaction), N = 32 | 0.0327 / 0.0210 / 0.0221 | 0.0649 / 0.0196 / 0.0184 | 0.279 / 0.057 / 0.069 | | 1 / 8 / 32 |
| Sedov, N = 32 | 0.0608 / 0.0427 / 0.0345 | 0.0987 / 0.0683 / 0.0604 | 0.0765 / 0.0523 / 0.0436 | +5.9 % / +4.3 % / +3.5 % | 2 / 39 / 113 |
| Sedov, N = 64 | 0.0406 / 0.0279 / 0.0202 | 0.0648 / 0.0418 / 0.0332 | 0.0467 / 0.0318 / 0.0250 | +3.0 % / +1.9 % / +1.5 % | 43 / 936 / 2587 |

One level on N = 32 matches or beats the uniform N = 64 run (Sod: 0.0098
against 0.0120; Sedov: 0.0427 against 0.0406, R error 4.3 % against
3.0 %), two levels halve the errors again where the second level is
tagged (the strong shock, Sedov; Sod's second level never is: the
threshold 0.2 · 8 is above its densities, so amr1 = amr2 there). The
Sedov anisotropy falls from 2.7 % to 1.3 % to 1.0 % (N = 32). The cost
is the whole-box refinement's: 8× the cells and 2× the substeps per level
(the N = 64 two-level Sedov run is a 256³ run, 60× the uniform wall), with
the flux-register and FillPatch overheads on top (`reflux` and `fill` in
the timers line; to profile). Before the two fixes of §1 ("three traps")
the refined runs were *worse* than the uniform ones: Sedov's radius
20–40 % short, the tube errors 4–10× larger.

Still to do in F2: the pancake and Kelvin–Helmholtz cases of the suite,
the full tier, and the composite of the finest data for a comparison at
the matched finest resolution; `cool_suppress_overdensity` and the sink
are built (below).

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

The user's decision: one coarse substep everywhere, no time sub-cycling
of the transport (the I-fronts in refined, dense gas move much slower
than light; quasars are not modelled). Design (2026-10-14), in two steps.

**F3a: the gas state on every level.** `GasState` grows one float
MultiFab per refined level (`levels[l].G`, the hydro hierarchy's BoxArray
and mapping of that level), so the gas lives in one place: the hydro
copies every level in and out (as level 0 today), a regrid rebuilds the
state's fine levels with the hydro's, the checkpoint writes them (the
`.amr` directory moves from the hydro to the state), and the push-down of
the level-0 changes the RT and PM make (§1) moves to the state, where it
belongs. Nothing changes for `amrex_max_level = 0` (bitwise), and the
level-0-only RT keeps working on refined runs as today.

**Built (2026-10-14/15, branch `stage-f3`).** F3a as designed. F3b in
`src/amrex/rt_level.cpp`: `Level` carries its refinement index, cells
per coarse substep (`rt_cells_per_substep` 2^l), per-level units and RT
parameters, a coverage mask (the cells under the finer level: no
emission, no physics, no statistics), an in-patch mask (where a ray
leaves the level) and the rays leaving it; the substep is split into
four phases (tau, cast, phys, advance), and `Hierarchy` (one `Level`
per level of the gas state, rebuilt at a regrid with its rays parked on
the parent) interleaves the levels' phases: after the tau phase a
refined level's ghost cells outside its grids take the parent's tau
(halved) and after the physics the parent's fgam (1 − √(1 − f)); after
the casts the deposits cross the boundary both ways (a level's ghost
deposits beyond its grids summed into the parent's cells, the parent's
deposits under the level shared among the 8 children exactly in the
int64 fixed point); after the advance the rays migrate (`leaving` up,
`collect_down` down, scaled to the other level's units with dr = rmax(nside)
− r); the statistics are summed over the levels and published once; after
`post` the finer levels' cells and gas are averaged down onto the
parents. Sources: `apply_sources` records the catalogue's positions
(`Simulation::rt_sources`) and the hierarchy puts a covered base cell's
Luv/sfr on the finest cell holding each source. The checkpoint carries
every level's rays in base units (a restart scatters them onto the base
level and the first prep moves them down). `amrex_refine_box` refines a
fixed region (tests). With `amrex_max_level = 0` the plain `Level` runs:
bitwise Stage E (176/176 on klone), and the refined smoke runs are
byte-identical across 1, 2 and 3 ranks and through a restart (the `sink`
tests run with one level).

**Validation, first result (2026-10-15, `tools/hyak/rt_amr_ifront.slurm`,
`tests/problems_full/make_ifront.py`, `tools/analysis/ifront_full.py`).**
One source at the centre of a uniform 4 Mpc/h box, the fixed physics,
ĉ = 0.1c, RT from z = 30 on 8 ranks; the I-front radius from the ionised
volume, R = (3V/4π)^{1/3}, in base cells. The refined 32³ run (one level
over the central half of the box) against the uniform 32³ run and the
uniform 64³ run as the reference:

| step | R (64³) | R (32³) | err | R (32³ + 1 level) | err |
|---|---|---|---|---|---|
| 3 | 4.138 | 4.159 | 0.022 | 4.142 | 0.004 |
| 4 | 4.868 | 4.897 | 0.029 | 4.854 | 0.014 |
| 5 | 5.448 | 5.479 | 0.031 | 5.430 | 0.018 |
| 6 | 6.030 | 6.039 | 0.008 | 5.986 | 0.044 |

The refined run tracks the 64³ front to a few hundredths of a cell, as
well as the base run does; here the front stayed inside the refined
region (radius 8). With the whole box refined (the RT entirely on level
1, the hydro and PM at 32³) the front reaches 8.60 base cells at step 10
against the 64³ run's 8.74 and the 32³ run's 8.65, with the same photons
in flight (526 against 552): the 1.5 % lag is the two-cell coarse
substep's (the k = 2 result of `docs/rt_multicell.md`, the decision's
price), not a defect. The first run with the refined region of radius 4,
so that the front crosses the level boundary, overshot by 0.57 cells
right after the crossing and had 8× the photons in flight: the kernels
keep a ray's photons per RT-cell volume (`ray_source`: ngam / vol_rt /
nH_avg), so a ray and the exchanged deposits must be rescaled by 8 across
a level (now done). The rerun (`rt_amr_ifront5`): the front crosses the
level boundary between steps 8 and 10 without a jump and continues on
the base level; refined 32³ against the 64³ reference, errors in base
cells: 0.004, 0.014, 0.018, 0.044, 0.061, 0.080, 0.108, 0.140 (steps
3–10; the uniform 32³ run: 0.022, 0.029, 0.031, 0.008, 0.010, 0.029,
0.060, 0.087), i.e. the fine level's lag carried across, nothing from the
crossing. CTest `amrex_rt_amr_ifront` (label `amr`, nightly) repeats the
three runs on 2 ranks. Three more traps on the way: the
refined levels' RT parameters must carry their own grid size for the
periodic wrap (`nrtr`, `hrtr`: a source ray at x = 32.5 on level 1
wrapped to cell 0); after the RT hierarchy's step the hydro's push-down
reference must be the current level 0 (the RT changed the fine cells
itself and averaged them down; the hydro pushed the same change down
again, and a species that had grown from ~0 overflowed in the ratio rule,
now guarded); and `fgam` is the transmitted fraction, so a half cell's is
√fgam.

**Next: the production configuration (scripts ready 2026-10-15, not yet
run).** `tools/hyak/amr_validate_runs.sh` sets up and submits three
10 Mpc/h, 0.1c runs with the abundance-matched sources on the backends:
`u128` (gas 128³, pm_refine = 2), `amr128` (the same ICs, two levels,
`amrex_refine_mass` 8) and `u256` (gas 256³, pm_refine = 1, the DM ICs of
`u128` byte for byte), so all three share the dark matter and the 256³ PM
mesh; `tools/analysis/amr_validate_compare.py` compares x_HII(z), τ_es,
T_M and the zdata1 fields (the reference averaged onto 128³) with `u256`.
`tools/hyak/amr_cost_probe.slurm` restarts the 512³ checkpoint of
`reion512_40mpc_f0.23` on 5 nodes for two steps with 0 levels, 2 levels,
and 2 levels on the 2048³ mesh (`pm_refine = 4` at the restart) and
reports the per-step times, the cells per level (`RADRAY_AMR_DEBUG`) and
the memory: the cost of the hierarchy at the target size before a
refined 512³ production run (whose PM at pm_refine = 2 adds FoF and P(k)
on rank 0, `docs/stage_e_pm.md` §4.1). The 128³ runs default to the
512³ target's PM, `pm_refine = 4` with 2× particles (the 512³ mesh, 256³
particles), and the refined gas in the PM (`pm_gas_levels = 1`, §4).

**F3b, as designed: rays on the levels.** Each level has its own RT cells (`RTCell`
FabArray), deposit and optical-depth FabArrays with their ghost layers,
and ray container. A ray lives on the finest level covering its position;
positions stay in level-0 cell units, so a migration between levels is a
re-indexing (the cell of the new level) and a `Redistribute` of that
level's container. Per coarse substep a level-l ray crosses 2^l of its
cells (`RTP::cells = 2^l`, the k-cell path of `docs/rt_multicell.md`, the
rate's crossing speed `k len_rt/dtsec_rt`, the exposure deposit and
`photon_conserve`, which the decision requires), then the level's
deposits are summed in its ghosts and the cell physics runs on its cells
with the level's `len_rt`. Coarse cells under fine grids see no rays and
no physics: after the substep's physics the gas and the RT diagnostic
fields of a refined level are averaged down onto the parent (as the hydro
does), so the level-0 fields are always the composite. Splitting and
merging stay per cell of the ray's level, so the angular resolution
follows the local cell size (the split radius is measured in local
cells); a ray entering a finer level splits at its next substep if the
local criterion says so, and one entering a coarser level merges there.
Sources emit on the finest level covering the source cell with the
emission nside of that level (`rrk::emit_nside` for k = 2^l). The
ray-weighted load balancing counts the rays of every level per level-0
box (the fine levels' mapping is their own; the balancing moves level 0
and the state's fine levels follow the hydro's regrid). The serial code
and the moment solvers stay single-level (refused with
`amrex_max_level > 0`).

**Validation.** Refinement off: bitwise Stage E (every existing test).
With refinement: `uniform_1src` and `uniform_8src` with a refined region
(the I-front crossing a level boundary) against the uniform run at the
matched finest resolution (the composite of §2) and the k = 2, 4 results
of `docs/rt_multicell.md` §4; the Iliev tests of the suite with
refinement; the full-loop smoke problem for rank/thread independence and
restarts; a 128³ production-configuration run with two levels against the
uniform 128³ and 256³ runs (x_HII(z), τ_es, T_M, the front positions).

## 4. The PM on the refined levels (`pm_gas_levels`)

*2026-10-15, branch `claude/pensive-goodall-kvpwgz`; the cheaper step
before multigrid gravity (Stage G), on the user's request, with
`pm_refine = 4` as the target.*

**Before** (§1, "Coupling"): the PM saw only level 0. The mesh took the
density of the level-0 gas cell each mesh cell lies in, the level-0 cell
was kicked by the mean mesh force over its r³ mesh cells, and the hydro's
push-down gave every fine cell its parent's velocity change. A finer PM
mesh (`pm_refine`) therefore sharpened the dark matter only: the refined
gas felt, and sourced, gravity at the level-0 resolution.

**Now** (`pm_gas_levels = 1`, the default; needs `pm_backend = amrex` and
`amrex_state = 1`, no effect without refinement, `pm_refine` a power of
2), in `src/amrex/pm_level.cpp`:

* **Deposit** (`composite_density`): every mesh cell takes the gas density
  of the finest level l ≤ log₂ r covering it (each level's cells spread
  over their q³ = (r/2^l)³ mesh cells; ParallelCopy of the levels in
  order, the finer overwriting). A level finer than the mesh is already
  in its average-down onto level log₂ r (the hydro, the sink, the RT and
  this kick all average down).
* **Kick** (`kick_fine`, after level 0's own kick): every refined cell is
  kicked (the unsplit kick of `kick_gas`) by the mean mesh force over its
  q³ mesh cells, or, for a level finer than the mesh, by the force of the
  mesh cell holding it; divided by r (velocities of every level are in
  level-0 cell units). Spread and average are transposes, so the gas-DM
  forces stay pairwise antisymmetric.
* **Average-down**: the kicked u(2:5) of every level onto its parent's
  covered cells (finest first, the mean of the 8 children in double), so
  level 0 is the composite; then `amrex_hydro_level0_is_current()`, so
  the hydro does not push the kick down again.
* `pm_gas_levels = 0` keeps the level-0 coupling (for comparison).
* `RADRAY_PM_MOMENTUM_CHECK=1` prints the DM + gas momentum change of every
  PM step relative to the summed |change| of the gas cells (a diagnostic).
* `make_production.py --ndm` sets the particle count apart from the mesh:
  `--N 512 --pm-refine 4 --ndm 1024` is the 2048³ mesh with 1024³
  particles (the ICs of the `pm_refine = 2` runs), and `--pm-gas-levels`.

**Tests** (`tests/pm_gas_levels_check.py`, CTest
`amrex_pm_gas_levels_check`, label `amr`, ~30 s on 2 ranks; run in a cloud
container with AMReX 26.10 fetched, OpenMPI 4.1.6):

1. *Equivalence.* Uniform gas, the same 16³ particles, z = 100 → 4.9,
   against the uniform 32³ run on the same 32³ mesh (B), the composite
   averaged onto 16³:

   | run | gas ρ rms | gas momentum rms | particle positions rms [cells] |
   |---|---|---|---|
   | 16³, no refinement, `pm_refine = 2` | 1.7e-1 | 1.7e-1 | 6.0e-3 |
   | 16³ + level 1 everywhere, `pm_gas_levels = 0` | 8.8e-2 | 1.7e-1 | 7.7e-3 |
   | 16³ + level 1 everywhere, `pm_gas_levels = 1` | **2.7e-3** | **7.5e-3** | **1.0e-4** |

   With the level coupling the refined run reproduces the uniform run at
   its resolution (what is left is the hierarchy's time stepping); with
   the level-0 coupling the particles were no better than without
   refinement.
2. *Momentum.* The smoke problem with two levels over the whole box:
   |ΔP_DM + ΔP_gas| / Σ|Δp_gas| per PM step ≤ 1e-7 at `pm_refine` = 1
   (both levels finer than the mesh), 2 (mixed) and 4 (both spanning mesh
   cells); the test requires < 1e-5.
3. `amrex_pm_gas_levels_ranks_smoke` (1, 2, 3 ranks byte-identical) and
   `amrex_pm_gas_levels_restart_smoke` (bit-identical restart), two levels,
   `pm_refine = 4`.

**The first refined production run** (`tools/hyak/amr512_runs.sh`, set up
2026-10-15, to launch on klone after the cost probe): `fesc512_b2.5_amr`
(40 Mpc/h, gas 512³ + two levels, 1024³ particles on the 2048³ mesh, 0.1c,
case B, f_esc(z) = 0.135 ((1+z)/7)^2.5 from the f_esc thread, the exposure
deposit, cooling in dense gas with the sink at 3000× the mean, FoF every
30 Myr, the ICs of `hmf512_40mpc_r2`) and its control `fesc512_b2.5_uniform`
(the same without refinement), 5 nodes each. Estimated 4–8 h for the
refined run (to be replaced by the probe's numbers).

**Not done (Stage G):** the particles still see one uniform mesh, so the
force resolution is the mesh's everywhere; a level finer than the mesh
gets the mesh resolution; the cost and memory of the uniform 2048³ mesh
at 512³ gas (`pm_refine = 4`) are to be measured on klone. Multigrid
Poisson on the hierarchy remains the long-term route.

## 5. Files and tests

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
| `tools/hyak/amr_validate_runs.sh`, `tools/analysis/amr_validate_compare.py` | the 128³ two-level production configuration against uniform 128³ and 256³ (§3, "Next") |
| `tools/hyak/amr_cost_probe.slurm` | the cost of 0 and 2 levels, and of 2 levels on the 2048³ mesh, on a 512³ checkpoint, 5 nodes (§3, "Next") |
| `src/amrex/pm_level.cpp` (`composite_density`, `kick_fine`, `momentum`), `src/pm_backend.cpp`, `include/radray/pm.hpp`, `config.hpp` | `pm_gas_levels` (§4) |
| `tests/pm_gas_levels_check.py` | the equivalence and momentum tests of §4 |
| `tools/hyak/amr512_runs.sh` | the refined 512³ production run and its unrefined control (§4) |

Tests: `ctest -L amrex_hydro` (27, refinement off: bitwise the serial
step), `ctest -L amr` (the label matches `amrex` too; the new tests are
`amrex_hydro_amr_ref_<p>_L<L>` (conservation, the level count) and
`amrex_hydro_amr_bitwise_<p>_L<L>_r<np>t<nt>_mgs<m>` (layout independence)
for sod_x (1 level, refine_mass 0.5), sound_x (1, 0.99) and random3d (2
levels, 8)), `cool_suppress_default_identity`, `cool_suppress_off_smoke`
(label physics).
