# Stage G: gravity on the hierarchy (`pm_level_solve`)

*Written 2026-10-22 on branch `claude/modest-mendel-8k29sb`. WP3 of the
zoom-in plan (`docs/zoom_turbulence.md` §4.3) and item 9 of
`docs/STATUS.md`'s "Next": the particles on the finest level covering
them, a CIC deposit per level, level 0 by the FFT (bitwise the Stage E
level) and every refined level by its own Poisson solve with AMReX's MLMG
(`MLPoisson`, Dirichlet values from the parent's potential), the 4-point
gradient and the kicks per level. Behind `pm_level_solve = 1` (default 0);
`pm_gas_levels` (the uniform `pm_refine` mesh fed by the levels,
`docs/stage_f.md` §4) stays as the comparison.*

## 0. Summary

| Item | State | Validation |
|---|---|---|
| The flag and its parameters (`pm_level_solve`, `pm_mlmg_iters`, `pm_mlmg_tol`, `pm_mlmg_bottom_sweeps`, `pm_mlmg_bottom_scale`, `pm_mlmg_verbose`; `include/radray/config.hpp`, `pm_backend.hpp`) | done | without refinement `pm_level_solve = 1` is byte-identical to 0 (`amrex_pm_level_solve_identity_smoke`) |
| The particles on their levels (one `ParticleContainer` with the hierarchy's levels, `Redistribute` by level after every drift and regrid) | done | the particle count and masses add up on every level (`RADRAY_PM_LEVEL_DEBUG=1`) |
| The deposit per level, the [1, 3, 3, 1] restriction onto the parents, the spill of the fine CIC clouds outside a level | done | §3: the far field of a refined particle is the base mesh's CIC force (0.9 % at 8 cells) |
| The level solve (MLMG, fixed V-cycles, smoother bottom, the multigrid depth fixed by the blocking factor) | done | level 1 over the whole box reproduces the uniform `pm_refine = 2` force to 3.7e-7; bitwise identical on 1, 2 and 3 ranks for `pm_refine` 1 and 2, bit-identical restarts |
| The forces and kicks per level, the average-down of the kicked momenta | done | momentum conserved to 1e-9 with levels over the whole box, 5e-4 with a region (the Dirichlet boundary) |
| Tests `ctest -L level_solve` (10 tests) | done | all pass here (4-core VM, OpenMPI 4.1.6, AMReX 26.10 with the linear solvers) |
| T1 of the zoom ladder (32³ base, two levels over the central eighth, against uniform 128³) | §4 | nightly `amrex_pm_zoom_t1` |
| Sub-cycled gravity of the refined levels (`pm_level_subcycle`, §6) | done 2026-10-31 | §6.9: the orbit energy error 77× smaller (64× expected), momentum as before, 1/2/3 ranks byte-identical, bit-identical restarts, byte-identical to 0 without a refined level; klone profile to run |
| Timing at 128³ with two levels | §5: klone, one node (2026-10-24, re-run 2026-10-28): the level solve's PM 0.16–0.17 s per step against 0.27 s for the uniform 512³ mesh; the step 1.5 s in both (the hydro hierarchy, 1.2 s); the VM numbers as well | the profile found the bottom-solve stall of a 144-cell level, fixed at no cost (width² sweeps, a floor of 8: 0.122 s per step against 0.124 s, 1.5e-11 against 4e-7); the 2026-10-28 step slowdown was `/gscratch` metadata latency, §5.1 |

Ground rules kept: `reference/` untouched; every new behaviour behind
`pm_level_solve = 1`, which only acts with `amrex_max_level > 0` (and
needs `amrex_state = 1`, `pm_backend = amrex`); the Fortran harnesses need
no mirror (the 2026-10-08 rule); results do not depend on the rank or
thread count (rule 5), the box layout or the restart.

## 1. The method

Level 0 is the Stage E level: the PM mesh of `Nm = pm_refine · Ngas`
cells, the int64 fixed-point CIC deposit, the FFT with the original's
Green's function (the 7-point Laplacian's), the float potential, the
4-point gradient, the CIC kick. Unchanged, bitwise, when no refined level
exists.

**Level l ≥ 1** (`FineMesh` in `src/amrex/pm_level.cpp`): the gas
hierarchy's level l (`gs::GasState::fine[l-1]`, its BoxArray and mapping)
gets a PM mesh of its gas boxes refined by `pm_refine`: `Nm 2^l` cells
per dimension in the periodic box, dx = 2^-l level-0 mesh cells.

1. **Particles.** The particle container has `amrex_max_level + 1`
   levels (level l's geometry: the mesh refined 2^l times in the same
   box; the refined levels' BoxArrays are set from the gas state after
   every hydro regrid, `sync_levels`). A particle's AMReX position is the
   centre of its home cell at the finest level's resolution, so that
   `floor(pos 2^l)` is its home cell on every level (exact in float while
   `Nm 2^max_level ≤ 2^23`); the exact float position in DM-grid units
   stays in the SoA, as before. The drift's `Redistribute` moves every
   particle to the finest level whose boxes hold its home cell (AMReX's
   `AmrAssignGrid`: finest first). A restart loads the particles onto
   level 0 and the first drift re-levels them; `amrex_balance` remaps
   level 0 only.
2. **Deposit.** Every level deposits its own particles (and stars) with
   the CIC at its own resolution, in density units (mass per level-0
   mesh-cell volume: a particle of mass m in a level-l cell is a density
   m 8^l) with `bits` fractional bits: 38 to level 2 and 3 fewer per
   level below (26 at most), so the range 2^(63 − bits) per cell holds
   the deep levels' densities (a particle's mass times 8^l must fit, else
   the step aborts). `SumBoundary` folds the ghost deposits into the
   level's cells; the CIC mass that fell outside the level stays in the
   uncovered ghost cells (the *spill*). `own(l)` = this deposit;
   `tot(l)` = own(l) + the finer levels' totals and spills restricted
   onto it, finest first.
3. **Restriction.** Not the mean of the 8 children: a fine CIC cloud
   restricted by the mean is a sharper source than the parent's own CIC
   (NGP-like at the parent's resolution), and the 7-point Laplacian's
   lattice anisotropy, which the CIC smoothing cancels, then gives 5–8 %
   force errors 8 parent cells away (measured, §3). The parent CIC's
   spectrum is the fine one's times cos²(k h_f/2), i.e. the fine field
   filtered with [1, 2, 1]/4; sampling it at the parent centres (half a
   fine cell off the fine centres) with the two-point mean gives the
   [1, 3, 3, 1]/8 stencil per dimension over the parent cell's two
   children and their outer neighbours. The parent then sees every finer
   level's mass as its own CIC would, times one more half-fine-cell
   top-hat (a slightly softer source: the far field is the base mesh's).
   Mass is conserved; the spill goes through the same stencil into the
   parent cells outside the level. Integer arithmetic throughout.
4. **Density.** Level l's density is its gas (the level's own cells,
   spread over the `pm_refine`³ mesh cells of each; level 0's gas is the
   composite the hydro averages down) + tot(l) + own(l − 1) prolonged
   piecewise-constantly (the parent-level particles whose CIC cloud
   reaches into the level), rounded to float once. The mean of level 0's
   density (an exact fixed-point sum) is the k = 0 mode the FFT drops;
   it is subtracted on the refined levels.
5. **Level solve.** lap φ = (2a/3)(ρ − mean) on the level's mesh, the
   same 7-point operator the FFT inverts on level 0, with
   `amrex::MLPoisson` + `MLMG`: the Dirichlet values on the level's
   boundary are the parent's potential (`setCoarseFineBC`, AMReX's
   coarse-fine interpolation), the domain is periodic, the parent's
   potential interpolated trilinearly is the initial guess, `pm_mlmg_iters`
   V-cycles (default 20) or until the residual max-norm is below
   `pm_mlmg_tol` (default 1e-10) of the initial one, red-black
   Gauss–Seidel (2 + 2 sweeps), the smoother as the bottom solver with
   at least `pm_mlmg_bottom_sweeps` sweeps (default 8, AMReX's) and, on a
   rectangular level, at least the square of its bottom grid's width
   (⌈`pm_mlmg_bottom_scale` × width²⌉, scale 1 by default; AMReX's 8 left a level whose width has an odd factor, 144 = 16 × 9
   cells, with a 9³ bottom grid that the V-cycles could not converge,
   §5.1); a solve that ends at `pm_mlmg_iters` above 1e-6 of the initial
   residual is reported (no
   dot products: every operation is per cell and the norms are maxima,
   so the solve is bitwise rank-, thread- and layout-independent). The
   multigrid depth is fixed by `sync_levels`: a level whose boxes fill
   their bounding box lets MLMG agglomerate and coarsen to a width of 2
   (layout-independent), any other level is coarsened log2(blocking
   factor × `pm_refine`) − 1 times (the boxes are multiples of the
   blocking factor however the rank count chops them). `setEnforceSingularSolvable(false)`:
   a level covering the whole box is a periodic (singular) problem that
   the fixed iterations handle without the mean-subtraction's reduction.
   The solution is stored as a float potential; its 3 ghost layers are
   the level's own values where covered and the parent's potential
   interpolated tri-quadratically elsewhere (the 4-point gradient of the
   boundary cells and the force in the ghost layer the kicks read). A
   level must be nested 3 parent cells inside its parent (AmrMesh nests by
   blocking_factor / 2 = 4; checked).
6. **Forces and kicks.** The 4-point gradient on each level, times 2^l
   (an exact float scaling: the force is in level-0 mesh units on every
   level), on the boxes grown by one cell; the DM kick per level (the
   CIC at the level's resolution; `amax`, `arms`, `vmax`, `vrms`
   accumulated across the levels in 128-bit fixed point and as maxima);
   the gas of level 0 kicked as before, the gas of every refined level by
   its own mesh force (the unsplit kick on the conserved variables), then
   the average-down of the kicked momenta onto the parents, finest first,
   and `amrex_hydro_level0_is_current()` as in `pm_gas_levels`.

With `pm_level_solve = 1` and no refined level in use the code path is
the Stage E one (the deposit of level 0 through the same routine with the
same arithmetic): byte-identical.

### 1.1 Costs and limits

- Each particle is deposited and kicked once, on its level; the
  restriction is one 64-point stencil per parent cell under a level.
- A level solve: `pm_mlmg_iters` V-cycles of 4 sweeps on the level's
  cells (and the coarser multigrid levels); 8–20 cycles reach 1e-8 to
  1e-11 (§3, §5).
- The Dirichlet level solve is not the composite solve: the forces across
  a level boundary are not pairwise antisymmetric (momentum changes by
  5e-4 of the summed gas |Δp| per step with a region, §3) and a particle
  near its level's boundary feels the mismatch between its fine and coarse
  potentials (the self-force of a refined particle 6 cells inside a level
  is 4e-6 of the force at one cell, §3). Nyx's composite multigrid would
  remove both; the level solve was chosen for its bitwise validation and
  its match with the sub-cycled hydro (`docs/zoom_turbulence.md` §4.3).
- `pm_refine` is allowed (every level's mesh is its gas grid times
  `pm_refine`) but the zoom does not need it: the levels refine the mesh.
  More particles than base-mesh cells are now allowed (`PMConfig`:
  `Ndm` a multiple of the mesh; the particle mass `dm_mass_ratio` is a
  double), which the T1 test uses (128³ particles on a 32³ base).
- The hydro's `amrex_refine_box` region is grown by `amrex_n_error_buf`
  and snapped to the blocking factor; `amrex_n_error_buf = 0` keeps a
  region exact.

## 2. Where things live

| File | What |
|---|---|
| `include/radray/config.hpp` | `pm_level_solve`, `pm_mlmg_iters`, `pm_mlmg_tol`, `pm_mlmg_bottom_sweeps`, `pm_mlmg_bottom_scale`, `pm_mlmg_verbose` |
| `include/radray/pm_backend.hpp`, `src/pm_backend.cpp` | the parameters, `check_pm_backend_param` (needs `amrex_state = 1` with refinement), `amrex_pm_force_at` (tests) |
| `include/radray/pm.hpp`, `src/pm.cpp` | `PMConfig::level_solve`, `max_level`, `blocking_factor`, `mlmg_*`; `dm_mass_ratio`; `Ndm` a multiple of the mesh |
| `src/amrex/pm_level.cpp` | `FineMesh`, `bits_of`, `interp_parent(_quad)`; `PMLevel::sync_levels`, `deposit_level`, `restrict_and_spill`, `rho_level`, `solve_level`, `force_level`, `kick_dm_level`, `kick_levels`, `average_down_momenta`, `force_at`; `fpos`/`home_l` (the positions by level), `for_tiles` |
| `src/amrex/pm_level_stub.cpp` | the no-AMReX stub of `amrex_pm_force_at` |
| `src/amrex/pm_subcycle.hpp` | `pm_level_subcycle` (§6): the PM's functions the hydro hierarchy calls (`pmsub::open`, `solve`, `kick`, `drift`, `regroup`, `done`); in `pm_level.cpp`: `sub_open`, `sub_solve` (`sub_density`, `parent_at`), `sub_kick`, `sub_drift`, `sub_regroup` (`to_level0`), `step_level0`, `sub_finish`; `amrex_pm_subcycle_call`, `amrex_pm_subcycle_finish` (`pm_backend.hpp`) |
| `src/amrex/hydro_level.cpp` | `HydroLevel::advance`/`step`: the calls of §6.1 |
| `tests/test_amrex_pm_subcycle.cpp`, `tests/pm_subcycle_check.py`, `tests/pm_subcycle_collapse_check.py` | §6.9 |
| `CMakeLists.txt` | the `LSOLVERS` component of `find_package(AMReX)`, `AMReX_LINEAR_SOLVERS ON` for `RADRAY_FETCH_AMREX` |
| `tools/hyak/build_amrex.sh` | AMReX with the linear solvers into `$RADRAY_SW/amrex-26.10-mlmg` |
| `tests/test_amrex_pm_levels.cpp` | the unit tests (§3) on a gas state with refined levels built by the test |
| `tests/pm_level_solve_check.py`, `tests/flag_identity_check.py` | the equivalence, momentum and identity checks |
| `tests/pm_zoom_t1_check.py` | T1 (§4) |
| `tools/hyak/pm_level_solve_profile.slurm`, `pm_level_solve_table.py` | the klone profile (§5) |
| `docs/stage_g.md` | this document |

## 3. Validation (`ctest -L level_solve`)

Registered with `RADRAY_AMREX=ON`; AMReX built with its linear solvers
(`AMReX_LINEAR_SOLVERS=ON`, the default; `tools/hyak/build_amrex.sh` builds
it so). All pass on the cloud VM (4 cores, OpenMPI 4.1.6, AMReX 26.10,
FFTW 3.3.10).

- `amrex_pm_levels_unit_r{1,2,3}` (`test_amrex_pm_levels`): a 32³ gas
  state with refined levels built by the test (no hydro), `pm_refine = 1`,
  `pm_mlmg_iters = 30`, `pm_mlmg_tol = 1e-11`:
  - *Ewald.* The force of a single particle at random points at distance
    d, rms relative error against the Ewald sum, for the uniform backend
    (`pm_refine` 1 and 2) and the level solve with level 1 over the
    whole box, level 1 over 16³ level-0 cells around the particle (6–8
    cells from the faces) and levels 1 and 2 over 16³ and 8³ cells:

    | d [cells] | uniform r = 1 | uniform r = 2 | level 1 everywhere | level 1 region | levels 1 + 2 |
    |---|---|---|---|---|---|
    | 0.5 | 0.916 | 0.521 | 0.521 | 0.521 | 0.206 |
    | 1 | 0.555 | 0.181 | 0.181 | 0.181 | 0.036 |
    | 2 | 0.177 | 0.028 | 0.028 | 0.028 | 0.011 |
    | 4 | 0.026 | 0.0091 | 0.0091 | 0.0092 | 0.0057 |
    | 5.7 | 0.016 | 0.0051 | 0.0051 | 0.0054 | 0.0052 |
    | 8 | 0.0094 | 0.0026 | 0.0026 | 0.0060 | 0.0060 |

    Level 1 over the whole box is the uniform r = 2 mesh at every
    distance; inside a region the force is the finer level's (the two-level
    case at d ≤ 2 matches the uniform r = 4 criteria of
    `test_amrex_pm.cpp`), and at 8 cells, where the points straddle the
    level boundary, the error (0.6 %) lies between the r = 2 (0.26 %) and
    r = 1 (0.94 %) meshes: the base level's force on the far side of the
    boundary is its own CIC force (0.92 % on the level-0 points). The
    self-force is 1.5e-7 (whole box) and 3.7e-6 (regions) of the force at
    one cell (< 1e-5 required).
  - *Against the uniform mesh.* Level 1 over the whole box on 16³ random
    particles: the forces at 2000 random points against the uniform
    `pm_refine = 2` backend's, max |ΔF| / max |F| = 3.7e-7 (the MLMG
    solution to 1e-11 against the FFT's; < 1e-4 required): the same
    discrete problem.
  - *Momentum.* One unsplit step from rest with 16³ random particles and
    random gas: |ΔP_DM + ΔP_gas| / Σ|p| = 5e-10 (level 1 everywhere),
    2e-9 (levels 1 and 2 everywhere; < 1e-6 required), 5.3e-4 and 6.4e-4
    with the regions (< 1e-2 required; the Dirichlet boundary, §1.1).
- `amrex_pm_level_solve_check` (`pm_level_solve_check.py`, 2 ranks,
  34 s): the equivalence test of `pm_gas_levels_check.py` with the level
  solve: uniform gas, the same 16³ particles, z = 100 → 5, the 16³ run
  with level 1 over the whole box against the uniform 32³ run (B):

  | run | gas ρ rms | gas momentum rms | particle positions rms [cells] |
  |---|---|---|---|
  | A1: `pm_gas_levels = 1`, `pm_refine = 2` (the uniform 32³ mesh) | 2.7e-3 | 7.5e-3 | 1.0e-4 |
  | G1: `pm_level_solve = 1`, `pm_refine = 1` (level 0 on 16³ by the FFT, level 1 on 32³ by MLMG) | 2.7e-3 | 7.5e-3 | 1.0e-4 |

  (within 1e-2, 2e-2, 1e-3 and no worse than 1.5× A1.) The momentum of
  the smoke problem with two levels and `RADRAY_PM_MOMENTUM_CHECK=1`:
  |ΔP| / Σ|Δp_gas| per PM step ≤ 7e-8 with the levels over the whole box
  for `pm_refine` 1 and 2 (< 1e-5 required), 2–3e-4 with them over half
  the box (< 2e-2).
- `amrex_pm_level_solve_ranks_smoke_r{1,2}`: the smoke problem with two
  levels over half the box, `pm_refine` 1 and 2, to z = 28 (12 steps, the
  gravity-only steps): 1, 2 and 3 ranks byte-identical (every field,
  movie frame and checkpoint); `amrex_pm_level_solve_restart_smoke_r{1,2}`:
  bit-identical restarts from step 9;
  `amrex_pm_level_solve_identity_smoke`: without refinement
  `pm_level_solve = 1` is byte-identical to 0. Nightly, through the RT
  steps 13–15 (`pm_refine = 1`): `amrex_pm_level_solve_rt_ranks` (1 and 3
  ranks) and `amrex_pm_level_solve_rt_restart` (from step 13), after the
  fixes below.

**Found on the way, and fixed (2026-10-29; the RT's, not the
gravity's):** with the RT on and a refined region that does not cover the
whole box (the rank and restart tests above used to stop before the smoke
problem's sources appear, z = 28), three things broke the bitwise rank
independence of the checkpoints and the exactness of a restart through an
RT step. None changes the results of an uninterrupted run except where
the run read undefined memory (3).

1. *The checkpoint's ray order.* AMReX chops a refined level's grids by
   the rank count, and the checkpoint gathered a refined level's rays box
   by box in the level's box order: the checkpoints of 1 and 3 ranks
   listed the same rays in different orders (`rank_check.py`, steps
   13–15: `restart/checkpoint_00{13,14,15}.bin` differed, every field,
   frame and statistic identical). The refined levels' rays are now
   gathered by order boxes independent of the layout (the base boxes
   refined to the level, chopped to 32 cells, where they meet the level;
   `gather_rays_ordered`), each box's rays sorted by home cell and bytes.
2. *A ray's `dr` through the checkpoint.* The kernel advances `dr` by the
   path length, so it drifts from `rmax(nside) − r` by round-off. The
   checkpoint holds every ray in base units; `rescaled` recomputed `dr` on
   the way up and again when the restart moved the rays down, which
   shifted a split by a substep (one photon cell fewer in the first
   substep after the restart). The refined levels' rays now keep their
   `dr` (checkpoint version 9), and `amrex_rt_restore_levels`, called right
   after the hydro levels are read, builds the RT levels of the checkpoint
   and moves the rays down with it, so they are those of the uninterrupted
   run; the first prep then regrids, if the hydro did, as the
   uninterrupted run does. Checkpoints of version ≤ 8 restore with `dr`
   recomputed, as before.
3. *`TF` read before it was set.* `exchange_deposits` weighs the deposits
   a refined level's rays make beyond its grids by the parent's `fgam` of
   the previous substep (`amrex_ghost_rule = 3`, the `TF` components 5–9).
   `TF` was never initialised, so the first substep after a level was
   built (the first RT step, a regrid, a `remap`, and every level after a
   restart) read undefined memory there; after a restart, where the
   uninterrupted run had the previous step's values, level 1 drifted
   while level 2's rays stayed identical. `TF` is now 0 when defined (the
   thick-gas weight, what freshly mapped memory held) and follows its
   cells through a `remap`; the checkpoint's level directory holds every
   level's `fgam` (`RT_fgam_<l>`, `amrex_rt_write_levels`), read back with
   the ghosts filled as at the end of a substep.

Also saved, though no test here needs it: the hydro hierarchy's push-down
reference (`Level_0_prev`, level 0 as the fine levels last saw it), so
that a change to level 0 between the RT's post and the checkpoint would
still be pushed down after a restart (a restart from step 13 is
bit-identical with and without it here).

Checked: `rank_check.py` with 1 and 3 ranks through steps 13–15, every
file byte-identical (checkpoints included); `restart_check.py` from step
13 (the first RT step) bit-identical with `pm_level_solve = 1` and with
`pm_gas_levels = 1`; the nightly tests `amrex_pm_level_solve_rt_ranks`
and `amrex_pm_level_solve_rt_restart` run these (~10 min each). Before
the fixes, on this branch, the `pm_gas_levels` restart from step 13 came
out identical (the `dr` and `fgam` differences happened not to change an
output bit; an earlier version of the branch saw it differ), and with the
level solve every RT field of step 14 differed.

Also run here: the `amrex_pm` label (the Stage E tests, `pm_level_solve`
off: 76 tests, all passing, the stored references bitwise) after the
particle-mass change (`dm_mass_ratio`).

## 4. T1: the zoom ladder's first rung

`tests/pm_zoom_t1_check.py` (nightly `amrex_pm_zoom_t1`, 2 ranks, 3.2
min here): the same Zel'dovich particles (128³, seed 7, 10 Mpc/h, z =
100 → 5, uniform gas, no RT) in three runs:

| run | gas | PM |
|---|---|---|
| A | 128³ | the uniform 128³ mesh (the reference) |
| B | 32³ + levels 1 and 2 over the central eighth (`amrex_refine_box=8,8,8,23,23,23`, `amrex_n_error_buf=0`: level 2 at the fine resolution there) | `pm_level_solve = 1`, `pm_refine = 1`: the 32³ mesh outside the region by the FFT, the levels by MLMG |
| C | the same hierarchy | `pm_gas_levels = 1`, `pm_refine = 4`: the uniform 128³ mesh everywhere (the Stage F coupling, the comparison) |

Inside the region (64³ fine cells; 238 000 particles at z = 5):

| | B (level solve) | C (uniform fine mesh) |
|---|---|---|
| particle displacement from A, rms [fine cells]: inner half of the region / whole region / outside | 0.14 / 0.30 / 0.84 | 0.002 / 0.006 / 0.07 |
| by distance from the centre (half-widths 0–0.25 / 0.5–0.75 / 0.75–1 / 1–1.25) | 0.17 / 0.14 / 0.38 / 0.49 | 0.001 / 0.002 / 0.009 / 0.03 |
| velocities, rms of the difference / rms velocity | 5.0e-2 | 2.0e-3 |
| cells above overdensity 1 / 10 (A: 0.0669 / 0.0023) | 0.0739 / 0.0024 | 0.0669 / 0.0023 |
| density on the fine grid, rms difference / rms | 0.33 | 0.011 |
| A's 12 halos with ≥ 128 particles in the region: found within 2 fine cells / within 40 % in mass | 12 / 11 | 12 / 12 |
| their mass ratios | 0.64–1.17 (rms Δlog₁₀ M 0.09) | 0.88–1.01 (0.018) |
| FoF halos ≥ 32 particles in the whole box / in the region | 99 / 47 | 394 / 41 |

Reading: C has the fine mesh everywhere, so its particles are A's to
1e-2 cells; B's base is 4× coarser outside the region, and the matter
that flows into the region carries the base's errors (0.8 cells rms
outside, 0.4 at the region's edge, 0.14 in its inner half). That is the
zoom's resolution, not the solver's: the solver's own accuracy is the
refined-everywhere test of §3 (1e-4 cells). Every halo of A above 128
particles is found at its position in B, with masses within 25 % but for
two (0.64, 0.75: halos at the region's edge, fed from the coarse base).
The test passes when B's inner-half particles stay within 0.3 fine cells
rms of A's, 90 % of those halos are found within 2 cells and 75 % of
them within 40 % in mass, and the density PDF fractions agree with A's to
25 % where A has ≥ 100 cells above the threshold.

## 5. Timing

The 4-core VM (2 ranks × 2 threads, OpenMPI 4.1.6): 128³ gas and 128³
particles from Zel'dovich ICs (10 Mpc/h, z = 100), no sources, 4 hydro
steps, `amrex_max_grid_size = 32`, `amrex_fine_max_grid_size = 64`; the
hierarchy with two levels over the central quarter width
(`amrex_refine_box=48,48,48,79,79,79`: level 1 80³ = 512 000 cells,
level 2 128³ = 2.1 M cells; the central eighth of
`tools/hyak/pm_level_solve_profile.slurm`, level 2 of 256³ = 17 M cells,
exceeds this VM's 15 GB through the hydro hierarchy). Mean per step
(`timing.txt`; the backend's phases from the cumulative timers):

| configuration | step [s] | PM [s] | hydro [s] | PM phases per step [s] |
|---|---|---|---|---|
| uniform 128³, no refinement | 1.68 | 0.40 | 0.88 | drift 0.10, deposit 0.10, fft 0.035, force 0.03, kick 0.12 |
| two levels, `pm_level_solve = 1` (`pm_refine = 1`) | 9.53 | 1.40 | 7.69 | level 0 as above; levels: deposit 0.07, **solve 0.56**, force 0.035, kick 0.04, sync 0.004 |
| two levels, `pm_gas_levels = 1`, `pm_refine = 4` (the uniform 512³ mesh) | 17.2 | 9.52 | 7.26 | deposit 2.23, **fft 4.44**, force 2.30, kick 0.14; gas levels: deposit 0.58, kick 0.08 |

The level solves: level 1 in 16 V-cycles to a residual of 4e-11 of the
initial one, level 2 in 7 (the finer level's initial guess, the parent's
potential, is better), 0.56 s for 2.6 M cells on 4 cores, i.e. about 0.2
µs per cell per V-cycle. The level solve costs 1.0 s per step over the
base PM here, the uniform 512³ mesh 9.1 s (its FFT, deposit and force
on 134 M cells); the hydro hierarchy (7.3–7.7 s, FillPatch and the
reflux on 64 boxes) dominates both refined runs, as in the 512³
production run (`docs/stage_f.md` §4).

### 5.1 klone (2026-10-24)

`tools/hyak/pm_level_solve_profile.slurm` on one cpu-g2 node, 24 ranks ×
8 threads, the branch at 7faad08 (before the bottom-sweep fix below),
10 hydro steps of the same problem with the levels over the central
eighth (`amrex_refine_box=32,32,32,95,95,95`: level 1 144³ = 3.0 M
cells, the region plus the nesting band of level 2; level 2 256³ =
16.8 M cells). The build job ran the 311 AMReX tests, all passing. Mean
per step (`tools/hyak/pm_level_solve_table.py`; the phases from the
cumulative timers):

| configuration | step [s] | PM [s] | hydro [s] | PM phases per step [s] |
|---|---|---|---|---|
| uniform 128³, no refinement | 0.23 | 0.012 | 0.03 | drift 0.003, deposit 0.003, fft 0.004, kick 0.002 |
| two levels, `pm_level_solve = 1` | 1.52 | 0.169 | 1.17 | level 0: drift 0.014, deposit 0.004, fft 0.003, kick 0.004; levels: deposit 0.007, **solve 0.124**, force 0.002, kick 0.003 |
| two levels, `pm_gas_levels = 1`, `pm_refine = 4` (the uniform 512³ mesh) | 1.58 | 0.267 | 1.13 | deposit 0.102, **fft 0.097**, force 0.024, kick 0.009; gas levels: deposit 0.006, kick 0.022 |

- The level solve's PM is 0.17 s per step against 0.27 s for the uniform
  512³ mesh, which costs four times the base mesh's resolution in every
  cell of the box. Both are small next to the hydro hierarchy (1.1–1.2 s:
  FillPatch 0.5 s, the reflux 0.12 s per step), which makes the refined
  step 6.6× the unrefined one.
- **The bottom solve.** Level 2 converged in 7 V-cycles to 5e-11, but
  level 1 ran to the 20-cycle cap at 4e-7: its 144 cells per side (16 × 9)
  coarsen only four times, to a 9³ bottom grid, and the 8 Gauss–Seidel
  sweeps of AMReX's smoother bottom solve reduce that grid's lowest mode
  by only ~0.5 per V-cycle. Reproduced here with one level of 144³ cells
  (`amrex_refine_box=28,28,28,99,99,99` on the 128³ problem, 2 ranks):
  20 cycles, residual 5.6e-7, 1.4 s per step; with
  `pm_mlmg_bottom_sweeps = 64` (now the default) 8 cycles, 3.5e-11,
  0.65 s per step. The sweeps run on the tiny bottom grid, cost nothing
  measurable on a level that coarsens to 2 cells, and keep the solve
  deterministic (no dot products). With the fix the klone level solve
  should drop from 0.12 s to about 0.05 s per step; the profile was not
  rerun. A level-1 residual of 4e-7 is still far below the float
  potential's precision, so the klone run's results were not affected.
- **Confirmed cycle by cycle** (2026-10-27, `pm_mlmg_verbose = 2`, one
  144³ level as above, 2 ranks): with 8 bottom sweeps the residual falls
  as on level 2 for three cycles (2.3e-2, 1.5e-3, 2.3e-4) and then only
  by a factor 0.70 per cycle (1.6e-4, 1.2e-4, 8.4e-5, … 1.4e-5 at cycle
  11): the 9³ bottom grid is never solved. With 64 sweeps every cycle
  gains about 20× (7.9e-5 at cycle 3, 1.7e-7 at 5, 3.4e-11 at 8).
- **The sweeps now scale with the bottom grid** (branch
  `claude/modest-mendel-8k29sb` after PR #7): the sweeps one V-cycle
  needs grow as the bottom grid's width squared, and the nesting pads of
  the zoom levels give widths with larger odd factors (528 = 16 × 33,
  1040 = 16 × 65). A rectangular level (one whose boxes fill their
  bounding box, so that MLMG agglomerates it) has a bottom grid of the odd
  part of its widths, known from its geometry alone, so it gets
  max(`pm_mlmg_bottom_sweeps`, width²) sweeps: 81 for 144³ (8 cycles to
  1.2e-11), about 4200 for a 1040-wide level (on 65³ cells, a fifth of
  that level's own smoothing work). Layout-independent, deterministic,
  and no Krylov bottom solver. A level that is not rectangular is
  coarsened only log2(blocking factor × `pm_refine`) − 1 times, and a wide
  irregular level would leave a wide bottom grid that no fixed sweep
  count solves; for those the new warning (`pm_level_solve WARNING: level
  l … reached pm_mlmg_iters … at relative residual … > 1e-06`, the first
  time per level and every 100th after; the timers line shows `bottom`,
  the bottom grid's width, and `capped N×`) is the guard. The zoom's
  levels are boxes and agglomerate; the density-tagged production levels
  are small patches whose bottom grids stay small. If a zoom level ever
  warns, padding its PM mesh to a box that coarsens to 2–4 cells is the
  next step.
- **The klone re-run** (2026-10-28, job 42496915, commit a28d87c: the
  width² rule with the 64-sweep floor of PR #7; the build job's 9
  `level_solve` tests pass; the same problem, nodes and layout as above):

  | configuration | step [s] | step without I/O [s] | PM [s] | hydro [s] | level solve per step [s] |
  |---|---|---|---|---|---|
  | uniform 128³ | 1.83 (before: 0.23) | 0.20 (0.21) | 0.013 | 0.03 | |
  | two levels, `pm_level_solve = 1` | 3.26 (1.52) | 1.51 (1.50) | 0.191 (0.169) | 1.16 | 0.146 (0.124) |
  | two levels, `pm_gas_levels = 1` | 3.22 (1.58) | 1.55 (1.55) | 0.269 (0.267) | 1.13 | |

  (The whole-step column is the filesystem's, see the last item; the
  step without `t_io` did not change.)

  Level 1 now converges: 8 cycles to 1.5e-11 (bottom grid 9, 81 sweeps),
  against 20 cycles stalled at 4.0e-7; level 2 as before (7 cycles,
  5.3e-11, bottom grid 2); no capped solve, no warning. But the level
  solve got **slower**, 0.146 s against 0.124 s per step, where on the
  4-core VM it was 2.2× faster. The hypothesis (confirmed below): the
  sweeps cost nothing to compute on a 9³ or 2³ grid, so the time is in
  what each sweep does besides the arithmetic (its boundary fill and
  kernel launches), and the 64-sweep floor made level 2 pay 8× AMReX's
  sweeps for a bottom grid of width 2 that needs none of them. Hence the
  floor is now AMReX's 8 (`pm_mlmg_bottom_sweeps` default), and only a
  wide bottom grid gets more (width²). To measure
  rather than infer, the profile job takes `EXTRA="--set
  pm_mlmg_verbose=1"`, which prints MLMG's own Solve / Iter / Bottom
  timers per level solve.
- **Measured (2026-10-28, 777cccf, the floor of 8; MLMG's timers, mean
  per solve; one node, 24 × 8):** level 1 (bottom 9³, 81 sweeps, 8
  cycles to 1.5e-11) 0.070 s, of which the bottom solve 0.045 s (81 × 8
  sweeps at about 70 µs each on a 729-cell grid: two boundary fills and
  two kernel launches per sweep, not arithmetic); level 2 (bottom 2³, 8
  sweeps, 7 cycles to 5.3e-11) 0.054 s, its bottom 0.0008 s. The level
  solve is 0.122 s per step (job 42554825; 0.136 s in job 42554081 with
  the same code on the same node: the latency-bound bottom varies by
  ~25 % between runs), against 0.124 s at 7faad08 where level 1 stalled
  at 4e-7 and 0.146 s at a28d87c, whose 64-sweep floor cost level 2
  about 64 × 7 sweeps ≈ 27 ms per solve for a 2³ grid. So the
  convergence fix is now free: same cost as before it, 4e-7 → 1.5e-11.
- **Fewer sweeps, more cycles** (`pm_mlmg_bottom_scale`, the width² rule
  as ⌈scale × width²⌉, default 1; the same job, level 1 per solve):

  | scale | sweeps | cycles | final residual | level-1 solve [s] | its bottom [s] | level solve per step [s] |
  |---|---|---|---|---|---|---|
  | 1 (default) | 81 | 8 | 1.5e-11 | 0.070 | 0.045 | 0.122 |
  | 0.5 | 41 | 11 | 3.3e-11 | 0.071 | 0.032 | 0.123 |
  | 0.33 | 27 | 15–16 | 7.5e-11 | 0.074 | 0.029 | 0.130 |
  | 0.25 | 21 | 19 | 5.7e-11 | 0.095 | 0.028 | 0.143 |
  | 0.15 | 13 | 20 (cap) | 9.2e-9 | 0.077 | 0.018 | 0.128 |

  Each cycle away from the bottom costs ~3 ms on level 1, so the sweeps
  saved come back as cycles: width² is the cheapest setting that
  converges, and the default stays. The key is there for wide bottom
  grids (a 1040-wide zoom level: 65², 4225 sweeps, ~0.3 s per cycle at
  this per-sweep cost) where the balance may differ; measure it there
  with `pm_mlmg_verbose = 1` before changing it. Padding the level's PM
  mesh to a width that coarsens to 2–4 (144 → 160 = 32 × 5 or 192 =
  64 × 3) would remove most of the bottom at +37 % or +137 % cells on
  the level; for 0.045 s of a 1.5 s step it is not worth building now.
- **The whole step 1.6 s slower in the second run: the filesystem, not
  the build.** All of it is `t_io` in `timing.txt` (0.024 → 1.6–1.7 s
  per step in every configuration), and on a step without output that
  is the gather of five fields onto rank 0 (the RT `sync` timer:
  0.124–0.139 s for the whole run in both jobs, unchanged), then
  `create_directories(step_NNNN)` and a 1 KB `stats.txt` on
  `/gscratch/astro`. No code on that path changed between the runs, and
  the binary links only the module's gcc 13 runtime and
  `/sw/ompi/4.1.6-2` (RUNPATH), FFTW statically, no conda library.
  Measured on compute nodes (2026-10-28, `runs/iotest_20261028`):
  - the same binary, node and rank count (job 42554783, 24 ranks):
    1.3–3.5 s per step of `t_io` writing to `/gscratch/astro`, twice,
    and 0.029–0.034 s writing to the node-local `$TMPDIR` in between;
  - 4 against 24 ranks alternating on one node (job 42556904): 0.04–0.07
    s per step for both, so the rank count does not matter;
  - a single process creating 20 directories with a 1 KB file each:
    `mkdir` 0.4 ms to 70 ms on average and up to 0.72 s, the file up to
    1.4 s, on `/gscratch/astro` (93 % full); 5 ms and 0.8 ms on
    `/gscratch/scrubbed`; 0.1 ms on `$TMPDIR`; good and bad minutes
    alternate within half an hour.
  The 2026-10-28 05:25 run started one minute after the 1024³ run ended
  and alongside its post-processing job on the same fileset. Without
  the I/O the three runs agree: level_solve 1.46–1.52 s, gas_levels
  1.55 s, uniform 0.20–0.21 s per step. Hence
  `pm_level_solve_table.py` now prints `io` and the step without it, and
  `LOCAL=1` makes the profile job run in `$TMPDIR` and copy back only
  the logs and `timing.txt`. Production runs pay the same per-step
  metadata cost (one directory and one small file per step): negligible
  against their 30–60 s steps, but a short benchmark on `/gscratch`
  should not be read without the `io` column.

## 6. Sub-cycled gravity on the refined levels (`pm_level_subcycle`)

*Design written 2026-10-31 (branch `stage-g-subcycle`), the user's decision
of 2026-10-29 (`docs/PLAN.md`, "Stage G"): the fine levels need gravity on
their own timescale. Behind `pm_level_subcycle = 1` (default 0), which acts
only with `pm_level_solve = 1`, `amrex_max_level > 0` and `pm_refine = 1`;
with 0 every result is byte-identical (§6.8). Results in §6.9.*

### 6.1 The time-stepping order per level

A base step [0, Δ] is two hydro calls of Δ/2 around the PM step. Without the
flag, every level is drifted, solved and kicked once, at Δ/2 (§1). With it:

- **Level 0** keeps its one kick per base step (drift–kick–drift, DKD): its
  particles drift Δ/2, level 0's density (every level's particles and gas at
  Δ/2, the finer levels restricted onto it) is FFT-solved, its particles and
  its gas are kicked by g₀(Δ/2)·Δ. The second drift of level 0's particles
  is deferred to the end of the second hydro call (the same arithmetic, so
  without refinement the step is bitwise unchanged). The refined levels are
  neither solved nor kicked here.
- **Level l ≥ 1**, in every one of its hydro substeps [T, T + h] (h is
  whatever the hydro takes; today h_{l−1}/2), kick–drift–kick (KDK):
  1. its state at T is the children's interpolation source (`Uold`), as
     before;
  2. K1: its gas (the valid cells and the two ghost layers the substep
     reads) and its particles kicked by g_l(T)·h/2;
  3. the hydro substep, unchanged;
  4. its particles drift by h;
  5. the closing solve at T + h: the deposit (§6.2), the MLMG level solve,
     the force g_l(T + h);
  6. K2: gas (valid cells) and particles kicked by g_l(T + h)·h/2;
  7. the children's substeps (each the same, recursively), the reflux and
     the average-down, as before;
  8. the particles of levels l and finer reassigned (§6.3).
  
  The next substep's K1 reuses g_l(T + h): one solve per substep.
  K(h/2) H(h) K(h/2) is Strang splitting of hydro and gravity, second order
  like the base step's H(Δ/2) K(Δ) H(Δ/2); for the particles it is the
  leapfrog, second order and symplectic within a level. The base step's
  form (half hydro, gravity, half hydro) per fine substep would need the
  hydro substep split in two, i.e. twice the hydro work and a new reflux
  pattern; KDK has the same order with one force per substep.
- **Level 0 for level 1:** after each level-0 hydro substep, before level 1's
  substeps, level 0 is solved at the substep's end (deposit and FFT) for
  level 1's boundary values only; level 0 is never kicked by it.
- **Opening:** at the start of each hydro call, after the hydro's push-down
  and regrid, the PM levels follow the gas state's, every particle goes to
  the finest level holding it, and every level is solved at the call's start
  (g_l for the first K1). Nothing is carried from one call to the next
  except the particles and the gas.

Per base step, level l ≥ 1 does N_l + 2 solves (N_l its substeps, 2·nsub·2^l
today, nsub the level-0 CFL substeps per half step; the 2 the openings) and
level 0 does 1 FFT for its kick, 2 for the openings and 2·nsub for the
boundary values. **This is more than "level 0 keeps the FFT once per base
step":** level 0's gravity (its kick) is still once per step, but level 1
needs level 0's potential at level 1's times, which only a level-0 solve at
the end of each level-0 substep provides (§6.2). The FFT is ~1 % of a step
(§5), so the extra 2–4 cost nothing measurable; the alternative is in §6.7.

### 6.2 What is interpolated in time

- **The boundary values.** A level-l solve at time t takes its Dirichlet
  values and the potential of its ghost layers from the parent's potential
  interpolated linearly in time between the parent's last two solves, i.e.
  between the start and the end of the parent's substep that contains t
  (the parent is ahead: its closing solve precedes its children's
  substeps). The time-interpolated parent potential replaces the parent's
  potential of §1 step 5, cell by cell (exact for t at either end).
- **The density of a solve at t** is level l's gas (the hydro's double
  state, rounded once), its own particles at t, the finer levels' particles
  (still at T, the parent's substep start) and the parent level's particles
  (`own(l−1)` prolonged, already at the parent's substep end) each moved
  **virtually** to t: deposited at x + v (t − t_q) with the drift's
  arithmetic, nothing stored (level 0's particles: from the call's start,
  their DKD positions, with the velocity of the drift they are in). This is
  a second-order prediction (Enzo deposits predicted positions the same
  way). The parent's closing solve therefore sees its children's particles
  predicted and its covered cells' gas from its own coarse substep (the
  coarse prediction of what the average-down will give); it is not redone
  after the children (a corrector solve would double the solves of every
  level that has children; §6.7). A virtual drift of more than 2 of the
  level's cells in a direction is dropped for that particle (it is
  deposited where it is; counted: `unpredicted` in the timers line). The
  rule depends on the particle alone: a first version tested the predicted
  cell against the particle's box, whose size depends on the rank count
  (AmrMesh chops the refined levels' grids by it), and 1 and 3 ranks
  diverged at the smoke problem's step 15 (`amrex_pm_subcycle_rt_ranks`).
- **The scale factor** in lap φ = (2a/3)(ρ − ρ̄) is linear in time across
  the call, as the base step's a at Δ/2 is.
- **The hydro's ghosts** interpolate the parent's states before its K1 and
  after its K2: both synchronous (gravity included up to that time), so
  the children's coarse–fine ghosts carry the parent's gravity in time.
  Level 0's gravity is not in its substeps (its kick is at Δ/2), so level
  1's ghosts from level 0 lag or lead by up to g₀Δ/2 in velocity: the same
  first-order mismatch the boundary of the level solve already has, at the
  base step's resolution.

### 6.3 Particles between levels

- A particle stays on one level for a whole substep of that level, and for
  the whole substep of its parent: the reassignment (AMReX `Redistribute`
  over levels l..finest) happens at the end of every level-l substep.
- A particle of level l that left its boxes during the parent's substep
  stays on level l (an *orphan*, at most 2 of its cells outside: one cell
  per substep, two substeps per parent substep), deposits its CIC mass into
  level l's ghost layers (the spill restricted onto the parent, §1) and is
  kicked by the force of level l's ghost layers (the parent's potential
  interpolated, §1); the forces of every level carry 4 ghost layers, the
  deposits 6 (the prediction above), the potentials 6 (the 4-point
  gradient), so a level must be nested
  4 parent cells inside its parent (AmrMesh's blocking_factor/2 with the
  default 8; checked). An orphan farther out, or a finer-level particle that
  left its parent's boxes within one substep, stops the run with a message
  (a particle faster than a cell per substep).
- **Level 0 ↔ level 1.** Level 0's particles sit at the call's start with
  the velocity of their drift (DKD), so they enter the refined levels only
  at the call openings: at the base step's start (both schemes synchronous)
  without a correction, at its middle minus the half of level 0's kick
  that the refined level will apply again, v −= g₀(Δ/2)·Δ/2. Particles of
  the refined levels that left level 1 become level-0 particles at the end
  of every level-0 substep (call time t_k, their level's force g): in the
  first call v' = v − g t_k, x' = x − v t_k + g t_k Δ/2, so that level 0's
  kick and first drift give v + g(Δ − t_k) and x + v(Δ/2 − t_k) to O(gΔ²);
  in the second call v' = v + g(Δ/2 − t_k), x' = x − v' t_k, which the
  deferred drift takes to the step's end. These few particles get level 0's
  first-order force placement; inside the levels the scheme is KDK.
- Stars of the sink are created as level-0 particles (as now: those of the
  first call are drifted by the PM step's first drift) and enter their level
  at the next opening.

### 6.4 Momentum: average-down and reflux

The fine levels' gas is kicked in the hydro's double state of the level. The
hydro's average-down after the children (all ten components) carries the
children's kicks into the parent's covered cells; the reflux corrects the
hydro fluxes only (gravity is a cell source, not a flux, and nothing of it
enters the flux registers). Level 0's kick at Δ/2 is applied to all of
level 0 and then its covered cells are overwritten by the average of the
fine levels' momenta and energies (`average_down_momenta`, as with
`pm_gas_levels`), so the fine gas gets no level-0 kick. The total momentum
is not conserved exactly: the Dirichlet boundary (as in §1.1) and, new, the
time offset between a level's mass and the other levels' kicks (the parent
feels its children's predicted mass, the children their parent's
interpolated potential; O(h²)). `RADRAY_PM_MOMENTUM_CHECK=1` prints the
momentum the kicks of a whole base step add against their summed |Δp|.

### 6.5 The push-down and `level0_is_current`

- Within a hydro call the fine levels' kicks change only the hydro's own
  state; the push-down compares level 0 with what the last call wrote out,
  so nothing new enters it.
- The PM step's level-0 kick changes level 0 between the calls; the
  average-down above restores the covered cells and then declares level 0
  current (`amrex_hydro_level0_is_current`), so the second call pushes
  nothing down, as with `pm_gas_levels` and `pm_level_solve` before.
- The RT acts between base steps on every level, averages down and declares
  level 0 current (unchanged); its gravity-free changes need no push-down.
  The push-down keeps carrying only what it carried before (the summaries'
  clamps, the halo metals).

### 6.6 Restarts, ranks and the timestep

- Checkpoints are written at the ends of base steps, where every particle is
  synchronous (level 0's deferred drift done, the refined levels at their
  last K2) and on its level after a full reassignment; no potential or force
  is carried (each call's opening recomputes them). The checkpoint is
  unchanged; a restart reproduces the uninterrupted run bit for bit (§6.9).
- Every new operation is per cell or per particle, the deposits are fixed
  point, the solves the fixed-iteration MLMG, the statistics 128-bit fixed
  point and maxima: no floating-point reduction enters the physics, so the
  results do not depend on the rank or thread count (the particles' order
  in a tile does not matter).
- The kick statistics: amax and arms over level 0's kick and every refined
  level's last kick of the step, vmax and vrms over every particle's
  velocity at the step's end; the timestep criteria are unchanged (in
  level-0 units), so the base step is as long as without sub-cycling.

### 6.7 Choices, and the options not taken

1. **Level 1's boundary values** (the one deviation from the request's
   "level 0 keeps the FFT once per base step"): a level-0 solve after each
   level-0 hydro substep (built; time-interpolated like every other level).
   The alternative, level 0's potential held constant over each half step
   (the opening's, then the PM step's), saves 2·nsub FFTs and their
   deposits per step and is first order in time at level 1's boundary.
2. **No corrector solve** after the children (built: the parent's closing
   solve uses its children's predicted mass). A corrector would make the
   parent's K2 and its children's last boundary values exact at the cost of
   one more solve per substep on every level that has children (~2× on the
   intermediate levels of a zoom).
3. **Level 0's particles enter the hierarchy only at the call openings**
   (built); fine particles leave to level 0 at every level-0 substep end.
4. **The substep ratio:** everything follows the hydro's substeps and their
   times (the parent interpolation uses times, not the ratio), so phase 2
   (each level's substep count from its own CFL limit) needs no change on
   the gravity side.

### 6.8 Without the flag

`pm_level_subcycle = 0` leaves every code path as it was: the new ghost
widths, the deposit's prediction and the per-tile skip are taken only with
the flag (the shared routines reduce to the old arithmetic: one ghost
layer, offset 0). The fast suite (its stored references bitwise) and the
`level_solve` tests pass unchanged (§6.9). With the flag and
no refined level in use, the step is byte-identical to the flag off
(`amrex_pm_subcycle_identity_smoke`: level 0's deferred drift is the same
arithmetic).

### 6.9 Validation (`ctest -L subcycle`)

All on the 4-core VM (OpenMPI 4.1.6, AMReX 26.10 with the linear solvers,
2026-10-31).

1. **An orbit inside level 2** (`amrex_pm_subcycle_orbit`,
   `tests/test_amrex_pm_subcycle.cpp`, 95 s on 2 ranks × 2 threads): a light
   particle on a circular orbit 6 level-2 cells around a clump of 511
   particles at one point (the clump at the centre of level 2 of a 32³
   base; a static universe; a cold, light gas), radhydro's base step at a
   fixed base step of 1/64 of the period, one orbit; the energy relative
   to the clump against a sub-cycled run at half the base step (the
   energy's swing of a few per cent around the orbit, the mesh's angular
   anisotropy against the direction-averaged Φ, is common to all runs and
   cancels):

   | | max \|E − E_ref\| / (v²/2) |
   |---|---|
   | `pm_level_subcycle = 0` (level 2 kicked once per base step) | 3.63e-3 |
   | `pm_level_subcycle = 1` (8 times: 2 hydro calls × 4 level-2 substeps) | 4.71e-5 |
   | ratio | 77 (expected (2·2²)² = 64; the request's (2^l)² = 16 counted one hydro call per base step) |

   At 16 base steps per orbit the cold gas falling onto the clump makes the
   hydro take 4 substeps per call (both runs), so the sub-cycled run has
   32 level-2 kicks per base step and the error ratio is 217; the test uses
   64 steps per orbit, where the hydro takes one substep per call. Checked:
   ratio ≥ 16.
2. **A collapsing overdensity** (nightly `amrex_pm_subcycle_collapse`,
   `tests/pm_subcycle_collapse_check.py`, ~6 min on 2 × 2): 64³ particles
   in 4 Mpc/h with one spherical Gaussian Zel'dovich displacement (σ = 5
   cells, linear central overdensity 1.686 at z = 15), uniform gas, z = 100
   → 6; B0/B1: 16³ base with levels 1 and 2 over the central half (level 2
   at 64³ resolution), A: uniform 64³ with its step capped at 1/8 of B1's
   median base step (level 2's substep; 348 steps against 26). The halo
   at z = 6 (A: M200 = 544 particles):

   | | centre offset [fine cells] | M200 / A's − 1 | profile, rms log difference (0.5–12 cells) |
   |---|---|---|---|
   | B0 (`pm_level_subcycle = 0`) | 0.025 | +0.075 | 0.088 |
   | B1 (`pm_level_subcycle = 1`) | 0.008 | 0.000 | 0.008 |

   With gravity once per base step the halo's core is 11 % too dense and
   the next shell 20 % too thin; sub-cycled, the profile is A's to 1 % in
   every shell. (A smaller halo, σ = 3, z_c = 12 → 7, M200 = 54: profile
   0.044 against 0.028, M200 −9 % against −6 %.)
3. **Ranks and restarts:** the smoke problem with two levels over half the
   box to z = 28 (12 gravity-only steps; `amrex_pm_subcycle_ranks_smoke`,
   `amrex_pm_subcycle_restart_smoke`): 1, 2 and 3 ranks byte-identical
   (every field, movie frame and checkpoint), a restart from step 9
   bit-identical. Nightly, through the RT steps: a restart from step 13
   (the first RT step) to 15 bit-identical (`amrex_pm_subcycle_rt_restart`,
   7 min), and 1, 2 and 3 ranks byte-identical to step 14
   (`amrex_pm_subcycle_rt_ranks`). Two things found on the way:
   - *The restriction's rounding* (`restrict_and_spill`, §1 item 3): a
     parent cell whose [1, 3, 3, 1] stencil spans several of the level's
     boxes gets one partial sum per box, and each was divided by 512 and
     truncated before the partials were added, so the parent's total
     depended on how the level's boxes are cut, i.e. on the rank count
     (AmrMesh chops the refined grids by it): 1 and 3 ranks differed by
     ~1e-13 of TOT(0) at the smoke problem's step 15 and the particles
     diverged. With the flag the partials carry their floor quotient and
     remainder through the sum and are divided once (exact). Without the
     flag the truncation is kept (byte-identity); it is a latent rank
     dependence of `pm_level_solve` itself, harmless in the tests so far
     (the float rounding of the density absorbs it) and left for a
     separate fix.
   - *The RT at step 15*: with the gravity fixed, the gas of every level and
     the particles entering step 15's RT are identical on 1 and 3 ranks
     (checksums, `RADRAY_PM_SUB_DEBUG=1`), and the RT's fields still
     differ (photons in flight in a few cells, 9e-23 against 0). The same
     happens with `pm_level_subcycle = 0` restarted from the same
     checkpoint (with the exact restriction, so that the particles agree):
     a rank dependence of the RT across two levels that this state
     triggers and the `pm_level_solve` nightly test's states do not, not
     the gravity's. Hence the rank test stops after step 14; the RT issue
     is reported for a separate fix.
4. **Momentum** (`amrex_pm_subcycle_momentum`,
   `tests/pm_subcycle_check.py`, `RADRAY_PM_MOMENTUM_CHECK=1`, the smoke
   problem, 6 steps): the momentum all kicks of a base step add against
   their summed gas |Δp|, with `pm_level_subcycle = 0`'s per-PM-step value
   next to it:

   | levels over | `= 0` | `= 1` |
   |---|---|---|
   | the whole box (periodic solves) | 2.0e-8 | 2.4e-8 (< 1e-5 required) |
   | half the box (Dirichlet boundary) | 2.4e-5 | 2.4e-5 (< 2e-2 required) |

   The level crossing's time offsets (§6.4) add nothing measurable here.
5. **T1** (`tests/pm_zoom_t1_check.py --subcycle`, nightly
   `amrex_pm_zoom_t1_subcycle`; §4's problem with run D = B with the
   flag), inside the region:

   | | B (§4) | D (`pm_level_subcycle = 1`) | C (uniform fine mesh) |
   |---|---|---|---|
   | displacement from A, rms [fine cells]: inner half / whole region / outside | 0.137 / 0.298 / 0.84 | 0.074 / 0.262 / 0.84 | 0.002 / 0.006 / 0.07 |
   | by distance from the centre (half-widths 0–0.25 / 0.5–0.75 / 0.75–1 / 1–1.25) | 0.165 / 0.135 / 0.379 / 0.487 | 0.068 / 0.073 / 0.343 / 0.454 | 0.001 / 0.002 / 0.009 / 0.029 |
   | velocities, rms difference / rms | 5.0e-2 | 4.4e-2 | 2.0e-3 |
   | density, rms difference / rms | 0.329 | 0.297 | 0.011 |
   | A's 12 halos ≥ 128 particles: within 2 cells / within 40 % in mass | 12 / 11 | 12 / 10 | 12 / 12 |
   | their mass ratios | 0.64–1.17 | 0.67–1.59 | 0.88–1.01 |

   Sub-cycling halves the inner region's displacement error (0.14 → 0.07
   fine cells); what is left there and at the region's edge is the coarse
   base's error carried in by the inflow (§4), which the flag does not
   touch. One halo at the region's edge came out 1.59× A's mass (0.79× in
   B: fed from the base, it scatters either way). Passes B's criteria.
6. **Cost.** The extra work is particle-bound (every solve deposits the
   particles of the levels l − 1 to the finest, the reassignments) and
   cell-bound (the solves). Per base step, two levels, 2 ranks × 2 threads:

   | problem | `= 0` step [s] | `= 1` step [s] | of which: deposits (level 0's particles) / solves / openings / reassignments / kicks [s] | solves, FFTs per step |
   |---|---|---|---|---|
   | smoke, 16³ base, 512 particles | 0.75 | 2.0 | 0.16 / 1.17 / 0.30 / 0.002 / 0.03 | 16, 4 (against 2, 1) |
   | collapse, 16³ base, 262 K particles | 0.27 | 1.21 | 0.32 / 0.32 / 0.23 / 0.18 / 0.12 | 16, 4 |
   | T1, 32³ base, 2.1 M particles | 1.40 | 6.2 | 2.08 (0.80) / 0.90 / 1.28 / 0.88 / 0.57 | 16, 4 |

   (The openings include their own deposits and solves.) In these small
   tests the hydro is cheap and most particles sit on level 0, whose
   extra deposits (the openings and the boundary solves, §6.7 item 1) are
   0.8 s of T1's 4.8 s; in a zoom the finest level holds most cells and
   particles, its own substeps are the cost (one solve, one deposit of
   its particles, two kicks and one reassignment per substep, against
   the hydro's two fills, predictor, fluxes and update), and the coarser
   levels add about one more deposit of the finest particles per finest
   substep. The klone profile of §5 with the flag on measures it at 128³:

   ```
   CONFIGS="level_solve level_subcycle" LOCAL=1 sbatch tools/hyak/pm_level_solve_profile.slurm
   python3 tools/hyak/pm_level_solve_table.py $RUNS/pm_level_solve_128 level_solve level_subcycle
   ```
7. **Not affected by the flag** (§6.8): the fast suite
   (`ctest -LE "long|nightly"`, 488 tests with the new ones, all passing,
   every stored reference bitwise; the live Fortran comparisons are not
   built here, the private submodule being absent) and `-L level_solve`.

## 7. How to run

```
# AMReX with the linear solvers (the default): -DAMReX_LINEAR_SOLVERS=ON in the
# recipe of docs/stage_e_pm.md section 7; klone: tools/hyak/build_amrex.sh
mpirun -np 4 build-amrex/radhydro RUN out --set rt_backend=amrex --set hydro_backend=amrex \
  --set pm_backend=amrex --set amrex_state=1 --set hydro_scheme=unsplit \
  --set amrex_max_level=2 --set amrex_refine_box=8,8,8,23,23,23 --set amrex_n_error_buf=0 \
  --set pm_level_solve=1 [--set pm_mlmg_iters=20 --set pm_mlmg_tol=1e-10 --set pm_mlmg_bottom_sweeps=8 --set pm_mlmg_bottom_scale=1 --set pm_mlmg_verbose=1] \
  --set step_timers=1      # the MLMG iterations and residual per level in the timers line
ctest --test-dir build-amrex -L level_solve
RADRAY_PM_LEVEL_DEBUG=1 ...   # the mass of every level's deposits after each PM deposit
# the refined levels' gravity in their own substeps (section 6; pm_refine = 1)
mpirun -np 4 build-amrex/radhydro RUN out ... --set pm_level_solve=1 --set pm_level_subcycle=1
ctest --test-dir build-amrex -L subcycle     # (the nightly ones with -DRADRAY_NIGHTLY_TESTS=ON)
RADRAY_PM_MOMENTUM_CHECK=1 ...   # the momentum the kicks of each base step add
# klone: the 128^3 two-level profile with and without the flag (section 6.9)
CONFIGS="level_solve level_subcycle" LOCAL=1 sbatch tools/hyak/pm_level_solve_profile.slurm
```
