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
| The flag and its parameters (`pm_level_solve`, `pm_mlmg_iters`, `pm_mlmg_tol`, `pm_mlmg_verbose`; `include/radray/config.hpp`, `pm_backend.hpp`) | done | without refinement `pm_level_solve = 1` is byte-identical to 0 (`amrex_pm_level_solve_identity_smoke`) |
| The particles on their levels (one `ParticleContainer` with the hierarchy's levels, `Redistribute` by level after every drift and regrid) | done | the particle count and masses add up on every level (`RADRAY_PM_LEVEL_DEBUG=1`) |
| The deposit per level, the [1, 3, 3, 1] restriction onto the parents, the spill of the fine CIC clouds outside a level | done | §3: the far field of a refined particle is the base mesh's CIC force (0.9 % at 8 cells) |
| The level solve (MLMG, fixed V-cycles, smoother bottom, the multigrid depth fixed by the blocking factor) | done | level 1 over the whole box reproduces the uniform `pm_refine = 2` force to 3.7e-7; bitwise identical on 1, 2 and 3 ranks for `pm_refine` 1 and 2, bit-identical restarts |
| The forces and kicks per level, the average-down of the kicked momenta | done | momentum conserved to 1e-9 with levels over the whole box, 5e-4 with a region (the Dirichlet boundary) |
| Tests `ctest -L level_solve` (10 tests) | done | all pass here (4-core VM, OpenMPI 4.1.6, AMReX 26.10 with the linear solvers) |
| T1 of the zoom ladder (32³ base, two levels over the central eighth, against uniform 128³) | §4 | nightly `amrex_pm_zoom_t1` |
| Timing at 128³ with two levels | §5 (this VM); klone: `tools/hyak/pm_level_solve_profile.slurm`, to run | |

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
   `pm_mlmg_tol` (default 1e-8) of the initial one, red-black
   Gauss–Seidel (2 + 2 sweeps), the smoother as the bottom solver (no
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
| `include/radray/config.hpp` | `pm_level_solve`, `pm_mlmg_iters`, `pm_mlmg_tol`, `pm_mlmg_verbose` |
| `include/radray/pm_backend.hpp`, `src/pm_backend.cpp` | the parameters, `check_pm_backend_param` (needs `amrex_state = 1` with refinement), `amrex_pm_force_at` (tests) |
| `include/radray/pm.hpp`, `src/pm.cpp` | `PMConfig::level_solve`, `max_level`, `blocking_factor`, `mlmg_*`; `dm_mass_ratio`; `Ndm` a multiple of the mesh |
| `src/amrex/pm_level.cpp` | `FineMesh`, `bits_of`, `interp_parent(_quad)`; `PMLevel::sync_levels`, `deposit_level`, `restrict_and_spill`, `rho_level`, `solve_level`, `force_level`, `kick_dm_level`, `kick_levels`, `average_down_momenta`, `force_at`; `fpos`/`home_l` (the positions by level), `for_tiles` |
| `src/amrex/pm_level_stub.cpp` | the no-AMReX stub of `amrex_pm_force_at` |
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
  levels over half the box, `pm_refine` 1 and 2, 12 steps: 1, 2 and 3
  ranks byte-identical (every field, movie frame and checkpoint);
  `amrex_pm_level_solve_restart_smoke_r{1,2}`: bit-identical restarts;
  `amrex_pm_level_solve_identity_smoke`: without refinement
  `pm_level_solve = 1` is byte-identical to 0.

**Found on the way (not the gravity's):** with the RT on and a refined
region that does not cover the whole box, the hydro hierarchy's fine-level
boxes are chopped by the rank count, and (a) the checkpoints of 1 and 3
ranks differ in their ray section from the first RT step on (the backend
gathers the rays box by box in the level's box order, sorted within a box;
`rank_check.py` on the smoke problem with the levels over half the box,
steps 13–15: `restart/checkpoint_00{13,14,15}.bin` differ on 3 ranks,
every field, frame and statistic identical) and (b) a restart from the
first RT step's checkpoint is not bit-identical (the rays re-enter the
tiles in another order). Both reproduce with `pm_gas_levels = 1`
(the same `rank_check.py` run: checkpoints 13 and 14 differ on 3 ranks; `restart_check.py` from step 13: every field of steps 14–15 differs); the existing rank and restart tests use a
region covering the whole box, whose layout no rank count chops. The Stage
G rank and restart tests therefore stop before the smoke problem's sources
appear (12 steps; the restart from step 9); the RT's checkpoint order on
a chopped layout is an item for the RT (`docs/stage_f.md` §3).

Also run here: the `amrex_pm` label (the Stage E tests, `pm_level_solve`
off) and the fast suite of the default build.

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

[filled below after the runs]

## 6. How to run

```
# AMReX with the linear solvers (the default): -DAMReX_LINEAR_SOLVERS=ON in the
# recipe of docs/stage_e_pm.md section 7; klone: tools/hyak/build_amrex.sh
mpirun -np 4 build-amrex/radhydro RUN out --set rt_backend=amrex --set hydro_backend=amrex \
  --set pm_backend=amrex --set amrex_state=1 --set hydro_scheme=unsplit \
  --set amrex_max_level=2 --set amrex_refine_box=8,8,8,23,23,23 --set amrex_n_error_buf=0 \
  --set pm_level_solve=1 [--set pm_mlmg_iters=20 --set pm_mlmg_tol=1e-8 --set pm_mlmg_verbose=1] \
  --set step_timers=1      # the MLMG iterations and residual per level in the timers line
ctest --test-dir build-amrex -L level_solve
RADRAY_PM_LEVEL_DEBUG=1 ...   # the mass of every level's deposits after each PM deposit
```
