# Stage E2: the PM gravity on AMReX (`pm_backend = amrex`)

*Written 2026-10-12 on branch `stage-e-pm`. Item 2 of Stage E in
`docs/framework_decision.md` §6: the dark matter in an AMReX
`ParticleContainer`, the CIC deposit onto a MultiFab of the PM mesh, the
Poisson solve, the force interpolation and the kicks on one AMReX level,
behind the new flag `pm_backend = amrex`. The default, `serial`, is the
original's PM (`src/pm.cpp`) and stays bit-for-bit. Stage E1 (the hydro on
AMReX) is a separate thread; the hydro files are untouched here.*

## 0. Summary

| Item | State | Validation |
|---|---|---|
| The shared PM kernels (`include/radray/kernels/pm_kernels.hpp`: CIC weights and corners, the 4-point force stencil, the Green's function) | done; `src/pm.cpp` calls them | the `pm` label unchanged: every stored gfortran reference byte-identical, the pm_refine and thread tests pass |
| The switch (`pm_backend`, `amrex_pm_max_grid_size`, `amrex_pm_fft`; `include/radray/pm_backend.hpp`) | done, AMReX-free; stub for builds without AMReX | the default build and the fast suite need neither AMReX nor MPI |
| The level `src/amrex/pm_level.cpp` in `pm_driver` and `radhydro` (`calc_kick` and `calc_kick_unsplit`, `pm_refine`) | done | §3: round-off against `pm_driver` on the 8 PM problems; bitwise identical for 1/2/4/8 ranks × 1/2 threads × box sizes 32/16/8, both FFTs; Ewald force, no self-force, momentum, the r = 2 pancake and its twin on the level; the full loop with `rt_backend = amrex` rank- and thread-independent, bit-identical restarts (legacy, unsplit, pm_refine = 2) |
| Timing on the 4-core VM, 128³ and 256³ DM | §4 | |
| klone (2026-10-13): the 64³ production configuration, three backends + `amrex_state = 1` vs serial; `pm_refine = 2` with 1024³ particles on 1/2/4 nodes | §3.1 (nightly `amrex_state_full64_vs_serial`), §4.1 | round-off against the serial code (x_HII 1.5e-6, T 3.2e-5, τ_es 2.2e-7 at 64³; gas ≤ 5e-6 rms, particles 3e-6 cells at 1024³); bitwise across 1/2/4 nodes and box sizes 32/64 for each FFT; PM 2.7× the serial on 1 node, 1.3 s per step on 4 |
| Tests `ctest -L amrex_pm` (65 tests, 73 with the serial fixtures) | done | all pass here |

Ground rules kept: `reference/` and the serial PM are untouched in
behaviour (the kernels the serial code now calls are its own arithmetic
moved into a header); every new behaviour is behind `pm_backend = amrex`;
the Fortran harnesses need no mirror (the 2026-10-08 rule for new flags).

## 1. What runs where

`pm::particle_mesh` and `pm::particle_mesh_unsplit` keep their call sites
in `radhydro` and `pm_driver`. With `pm_backend = amrex` they do

```
calc_drift (host)  ->  amrex_pm_step (the level)  ->  calc_drift (host)
```

in place of the host's `link_dm`, `calc_density_field`,
`calc_gravity_field` and `calc_kick`. The two drifts stay on the host:
they are per-particle float arithmetic that every rank does identically on
its full copy of `xv_dm`, and keeping them there saves a `Redistribute`
(the particles never move while they are on the level).

**MPI model of this stage.** As in Stage D, the hydro and the RT hold the
full state on every rank. So, like `amrex_rt_prep`/`amrex_rt_post`, the
level copies the host state in at the start of the PM step and gathers
the results back onto every rank at the end:

- in: every rank walks the full `xv_dm` (after the first drift), computes
  each particle's home mesh cell with `link_dm`'s formula
  (`int(x·x_dm_m) mod Nm`), asks an AMReX `ParticleLocator` which box owns
  the cell, and creates the particles of its own boxes in its own tiles
  (no `Redistribute`; the global count is checked); every rank fills its
  gas boxes from `u_gas` and `v_gas` (or the dual-energy field `e`);
- out: each rank writes its kicked particle velocities (by the particle
  id = the host index) and its kicked gas cells into zeroed host arrays,
  then one MPI bitwise-OR all-reduce per array (exact, since the pieces
  are disjoint and the unchanged entries are identical on every rank;
  chunked at 2²⁸ elements). `pm_driver` also gathers the potential (its
  `phi.bin` output; `PMState::gather_potential`, off in `radhydro`).

The level therefore holds no state between steps: checkpoints, restarts
and the FoF/halo code are the host's, unchanged, and a restart is
bit-identical by construction.

The process functions are those of the RT backend (`backend_init`,
`io_rank`, `backend_nranks`, … in `include/radray/rt_backend.hpp`): one
`amrex::Initialize` serves both backends; `pm_driver` now calls them and
rank 0 writes.

### 1.1 The level

One uniform level of the PM mesh (`Nm = pm_refine·Ngas` cells per
dimension, periodic `Geometry` with dx = 1, `BoxArray::maxSize(mgs)` with
`mgs = amrex_pm_max_grid_size` (0 = `amrex_max_grid_size`, a multiple of
`pm_refine` that divides `Nm`), the default `DistributionMapping`). The gas
grid is the same `BoxArray` coarsened by `pm_refine`, so each gas cell's
r³ mesh cells lie in one box (the spread and the average of `pm_refine`
need no communication).

| Data | AMReX object | Ghosts |
|---|---|---|
| particles | `ParticleContainerPureSoA<9, 0>`: AMReX position = the centre of the home mesh cell, id = the host index ip, SoA reals = position in DM-grid units and velocity (the original's float `xv_dm(6,ip)`) | |
| `acc` | `FabArray<BaseFab<Long>>`, the CIC mass in int64 fixed point (2⁻³⁸: resolution 3.6e-12, range 3.4e7 per cell) | 1 (the CIC reach) |
| `pot` | `FabArray<BaseFab<float>>`, the density, then the potential | 3 (the 4-point gradient on the box plus one ghost layer) |
| `frc` | the mesh force, 3 components (the float values of `calc_kick`'s `force(i,j,k)`) | 1 (the CIC reach of the kick) |
| `gas` | `u(5)` and `v(3)` (or `e`) on the gas grid | 0 |

One step:

1. **deposit** (`cic_mass_assignment`): per box, each thread adds
   `to_fx(w·mass)` into its own int64 buffer of the grown box (no atomics;
   `rpk::cic_corners` with the Fortran's weights), the buffers are summed
   into `acc`, then `SumBoundary` folds the ghost deposits into their
   owners (integers: exact, associative). The density is the gas density of
   the gas cell the mesh cell lies in (`calc_density_field`'s spread) plus
   the DM mass, rounded to float once (the original rounds after every
   particle);
2. **Poisson** (`calc_gravity_field`): forward FFT, every spectral
   coefficient multiplied by `rpk::green_wphi` (the original's
   `Aphi/(ka²+kb²+kc²)` with the k = 0 mode zeroed, double product stored as
   float), backward FFT, the real output scaled by the float
   `1/(Nm³)` as the MKL/FFTW shim does; `FillBoundary` of the potential;
3. **force**: `rpk::force_stencil` on each box grown by one cell;
4. **DM kick** (`accelerate_dm`): per particle the CIC weights from the
   SoA position, `rpk::particle_force` from `frc` (indices mapped to the
   periodic image nearest the home cell, range-checked), `v += (f/x_dm_m)·dt`;
   `amax`/`vmax` as maxima, `arms`/`vrms` summed in 128-bit fixed point
   (2⁻⁵⁶) and reduced as 32-bit limbs, so `clk.arms_dm` etc. are rank- and
   thread-independent;
5. **gas kick** (`accelerate_gas`): the three variants of `src/pm.cpp`
   (MF `v_gas += F dt`, non-MF on `u_gas(2:5)`, unsplit on the conserved
   lab-frame variables with `rrh::thermal_energy`), with `gas_force` = the
   mesh force (r = 1) or the k, j, i-ordered double average over the r³
   mesh cells divided by r (r > 1), as on the host;
6. **gather** (§1).

### 1.2 The FFT

`amrex_pm_fft = pencil` (default) is a distributed 3-D real FFT built from
1-D FFTW transforms: the real data are copied (`ParallelCopy`) onto
z-slabs (full x and y per box, one slab per rank), transformed along x
(r2c) and y (c2c) line by line, copied onto y-slabs and transformed along
z; the backward transform is the reverse. Every line uses the same four
plans (`FFTW_ESTIMATE | FFTW_UNALIGNED`, one transform per call,
single-threaded, made once), executed on thread-local buffers with the
new-array interface. So every line gets exactly the same arithmetic
whatever the rank and thread layout: the result is bitwise rank- and
thread-independent by construction.

`amrex_pm_fft = amrex` uses `amrex::FFT::R2C<float>` (AMReX built with
`AMReX_FFT=ON`, which on the CPU links FFTW in both precisions). It plans
`fftwf_plan_many` transforms over each rank's local pencils with FFTW's
planner set to the OpenMP thread count. In the matrix of §3 it is also
bitwise identical for 1–8 ranks and 1–2 threads, but that rests on FFTW
choosing the same algorithm for the same 1-D size whatever `howmany` and
the thread count, which is not guaranteed; the pencil transform is the
default because its determinism does not depend on the planner. The two
differ from each other, and from the serial 3-D FFTW plan of `src/fft.cpp`,
at round-off. Where AMReX uses the threaded FFTW, the amrex path resets
FFTW's planner thread count to 1 after each use, so the host code's own
plans (the hydro's 1-D transforms, the serial PM) stay the serial ones. On
klone AMReX is built with `AMReX_FFTW_IGNORE_OMP` (`tools/hyak/
build_amrex.sh`): it never initialises FFTW's threads, its plans are
single-threaded, and the reset is skipped (`fftwf_planner_nthreads() > 1`
guards it; an unguarded `fftwf_plan_with_nthreads` before
`fftwf_init_threads` runs `fftwf_cleanup`, which left the solver's plans
undefined and crashed their destruction at exit at 1024³, §4.1).

### 1.3 Determinism

Every cross-cell or cross-rank sum is an integer (the int64 deposits, the
128-bit statistics) or a maximum, every line of the FFT is transformed
with the same plan, and the particle and gas kicks are per element. The
order of the particles in a tile (which depends on the thread count) does
not enter any result. Hence the bitwise identity across ranks, threads and
box sizes of §3.

Against the serial PM the differences are round-off: the serial CIC adds
each particle's mass to a float density in linked-list order; the serial
3-D FFTW plan and the per-line transforms round differently. Both are
~1e-7 per cell of the density and the potential; the 4-point gradient and
30 steps of the leapfrog amplify that to a few 1e-6 of the maximum force
where the force is small (the growth test: 5e-6), see §3.

## 2. Where things live

| File | What |
|---|---|
| `include/radray/kernels/pm_kernels.hpp` | `rpk::MeshGeom`, `cic_weights`, `cic_corners`, `force_stencil`, `particle_force`, `green_wphi`: the per-particle and per-cell arithmetic of `pm_tools.f90`, shared by `src/pm.cpp` and the level |
| `include/radray/pm_backend.hpp`, `src/pm_backend.cpp` | the switch and its parameters (AMReX-free): `pm_backend = serial \| amrex`, `amrex_pm_max_grid_size` (0 = `amrex_max_grid_size`), `amrex_pm_fft = pencil \| amrex`; `check_pm_backend_param`; the level entry points `amrex_pm_step`, `amrex_pm_step_unsplit`, `amrex_pm_fields` (tests), `amrex_pm_timers` |
| `include/radray/config.hpp` | `PMBackend`, the three fields |
| `include/radray/pm.hpp` | `PMConfig::backend/amrex_max_grid_size/amrex_fft` (copied from `Config` by `PMConfig::make(cfg)`), `PMState::amrex` (the level), `gather_potential` |
| `src/pm.cpp` | calls the kernels; `particle_mesh`/`particle_mesh_unsplit` dispatch to the level between the drifts |
| `src/amrex/pm_level.cpp` | the level (C++20, target `radray_amrex`, links `AMReX::amrex` and FFTW3f) |
| `src/amrex/pm_level_stub.cpp` | the no-AMReX build: the level functions throw |
| `apps/pm_driver.cpp` | `backend_init`, `--set key=value` (`pm_backend`, `amrex_pm_*`, `amrex_max_grid_size`), rank 0 writes, the backend's timers |
| `apps/radhydro.cpp` | the `pm_backend` parameters, the summary line, the timers |
| `tests/test_amrex_pm.cpp`, `tests/pm_ewald.hpp`, `tests/amrex_pm_check.py` | the unit tests, the Ewald force (shared with `test_pm_refine.cpp`), the run/compare script |
| `tools/hyak/build_amrex.sh`, `build_amrex_gpu.sh` | AMReX with `AMReX_FFT=ON` (and a user-space double-precision FFTW, which AMReX's FFT module links although the backend uses only the float transforms); the `amrex\|amrex_pm` labels |

The CMake edits are additive: `src/pm_backend.cpp` in `radray`,
`src/amrex/pm_level.cpp` in `radray_amrex` (which now links FFTW3f), the
stub in the default build and the sanitizer source lists, the `FFT`
component in `find_package(AMReX)` and `AMReX_FFT ON` for
`RADRAY_FETCH_AMREX`. `simulation.hpp` is untouched (the level hangs off
`PMState`), the hydro files are untouched.

## 3. Validation (`ctest -L amrex_pm`, 65 tests)

Registered only with `RADRAY_AMREX=ON`; the MPI launcher is
`RADRAY_MPIEXEC` as for the `amrex` label. All pass on the cloud VM
(4 cores, OpenMPI 4.1.6, AMReX 26.10, FFTW 3.3.10).

- `amrex_pm_unit_r{1,2,4}_{pencil,amrex}` (`test_amrex_pm`, the
  counterparts of `test_pm_refine.cpp` on the level, 1/2/4 ranks): the
  force of a single particle against the Ewald sum for pm_refine = 1, 2, 4
  (the backend's mesh force gathered by `amrex_pm_fields`, interpolated with
  the host's `force_at`; the same criteria as the host test: 4 % at 4 mesh
  cells, 1.2 % at 8, r = 2 as accurate at d/2); no self-force (< 1e-5 of
  the force at one cell); the gas kick against Ewald at the gas cell
  centres (2 %/1 % rms at 3–8 cells for r = 1/2); the backend's potential
  and particle forces against the serial code on random particles and gas
  (within 1e-5 of the maximum; measured 2e-7 and 1e-6); total momentum of
  DM + gas conserved by one step to 2e-6 for MF, non-MF and unsplit kicks
  with pm_refine = 1, 2, 4 (measured ≤ 1e-7, as on the host);
- `amrex_pm_roundoff_<problem>` for the 8 PM problems (`zeldovich_3d`,
  `growth_3d`, `pancake_x`, `pancake_z_nomf`, `zeldovich_3d_r2`,
  `pancake_x_r2`, `pancake_z_nomf_r2`, `pancake_x_n32`): `pm_driver --set
  pm_backend=amrex` with 2 ranks × 2 threads and 8³ boxes against the serial
  `pm_driver` run of the `pm` label, every value of `xv_dm`, `u_gas`,
  `v_gas` and `phi` within 1e-4 relative or 1e-5 of the field's maximum at
  every step, the stats within 1e-4 (the DM velocity and position
  statistics relative to the velocity scale and the box);
- `amrex_pm_momentum_<problem>`, `amrex_pm_pancake_{pancake_x,
  pancake_x_r2, pancake_x_n32, pancake_z_nomf, pancake_z_nomf_r2}`,
  `amrex_pm_growth_3d`: the physics checks of the `pm` label
  (`pm_physics_check.py`: momentum conservation over the run, the analytic
  Zel'dovich pancakes, the linear growth) on the backend's runs;
  `amrex_pm_r2_twin_pancake_x`: the pm_refine = 2 run's dark matter is
  bit-identical to its pm_refine = 1 twin on a gas grid as fine as the mesh;
- `amrex_pm_bitwise_<problem>_<fft>_r<ranks>t<threads>_mgs<box>`: for
  `zeldovich_3d` (pm_refine 1) and `zeldovich_3d_r2` (pm_refine 2), both
  FFTs, (ranks, threads, box) ∈ {(1,2,16), (2,1,16), (2,2,8), (4,1,16),
  (4,2,8), (8,1,8), (8,2,16)} byte-identical to the 1 × 1 run with one 32³
  box (every `.bin` and `stats.txt`);
- the full loop with `rt_backend = amrex` and `pm_backend = amrex` on
  `problems_full/smoke`: `amrex_pm_full_ranks_smoke` (1, 2 and 4 ranks
  byte-identical: fields, movie frames, checkpoints),
  `amrex_pm_thread_smoke` (1, 2, 4 threads, unsplit hydro),
  `amrex_pm_restart_smoke`, `amrex_pm_restart_smoke_unsplit`,
  `amrex_pm_restart_smoke_pm_refine` (`onthefly_r2`: Ndm = 2N through the
  checkpoint).

Also run here: the same bitwise matrix at 8 ranks on `pancake_z_nomf_r2`
(non-MF kick); the whole `amrex` label (54 tests) and the `pm` label in
the AMReX build; the fast suite in the default build.

**The full loop against the serial PM.** On the smoke problem with the
RT backend, swapping the serial PM for the backend changes the histories
by 2e-5 (T_V, T_M) and single cells by up to 70 % by step 15, while
swapping the serial RT for the RT backend (Stage D) changes them by 1e-8
to 3e-7. That is the sensitivity of the 16³ smoke problem to the force
round-off (1e-6 of the maximum force, §1.3) through the hydro, not an
error of the backend: the gravity-only problems above follow the serial
trajectories to 1e-5 over 30 steps, and the unit tests pin the forces. The
64³ production-configuration comparison is §3.1.

### 3.1 The 64³ production configuration (klone, 2026-10-13)

Nightly test `amrex_state_full64_vs_serial` (labels `amrex;amrex_state;
full;nightly`; `tests/amrex_full64_check.py ... state`): `make_production.py
--N 64 --recomb-rad 0` with today's production defaults (10 Mpc/h, 0.01c,
factor 16, on-the-fly sources, `legacy_rt = legacy_physics = 0`,
`hydro_scheme = unsplit`, `redshift_per_distance = 1`, case B), z = 100 → 5
in 27 hydro steps. The serial code (4 threads, `deterministic_parallel`)
against the three backends with the state on the levels (`rt_backend =
hydro_backend = pm_backend = amrex`, `amrex_state = 1`, 2 ranks × 2
threads), compared as `amrex_full64_vs_serial` does for the RT: the
histories of `stats.txt` (x_HII and T, volume and mass weighted, every
step), τ_es from the movie frames, and the per-cell fields of the last
step. Measured on klone (job 41832828, the build of this branch; the rows
other than the test are runs of the same problem with fewer backends,
the RT + hydro run bitwise the RT-only one):

| 64³, serial vs backends | x_HII,V | x_HII,M | T_V | T_M | τ_es | fHI rms (z = 5) | T rms | ρ_gas rms |
|---|---|---|---|---|---|---|---|---|
| RT backend only (2 × 2) | 3.2e-8 | 3.7e-7 | 2.1e-7 | 1.1e-6 | 1.1e-8 | 1.6e-4 | 6.7e-6 | 3.2e-6 |
| RT + hydro backends | 3.2e-8 | 3.7e-7 | 2.1e-7 | 1.1e-6 | 1.1e-8 | 1.6e-4 | 6.7e-6 | 3.2e-6 |
| PM backend only (serial RT and hydro) | 1.2e-6 | 1.6e-6 | 3.2e-5 | 1.9e-5 | 1.5e-7 | 5.4e-2 | 6.4e-4 | 5.2e-5 |
| **the three backends, `amrex_state = 1` (the test)** | 1.2e-6 | 1.5e-6 | 3.2e-5 | 1.9e-5 | 2.2e-7 (0.05728) | 5.4e-2 | 7.1e-4 | 5.2e-5 |
| tolerance of the test | 1e-3 | 1e-3 | 1e-3 | 1e-3 | 1e-3 | 0.15 | 1e-2 | 5e-4 |

(history columns: the largest relative difference over the steps; fields:
rms difference over rms value.) The three backends with `amrex_state = 1`
are bitwise identical to the same three with the host state
(`amrex_state = 0`), as at 256³ and 512³ (`docs/stage_e_state.md` §4).

Reading: the histories and τ_es follow the serial code to 3e-5, the gas
density to 5e-5 and the temperature to 7e-4 per cell, all from the PM
backend's round-off (its fixed-point CIC and per-line FFT, §1.3; the hydro
backend is bitwise the serial step and the RT backend contributes 1e-4).
The 5 % rms in fHI at z = 5, after overlap, is not a density or
temperature effect: per cell, the fHI difference is the Γ_HI difference
with the sign reversed (correlation −0.9998, slope −0.99; median |ΔΓ/Γ|
1e-3, rms of Δln Γ 4.5e-2, against 1.3e-4 in Δln ρ). The 1e-4 density
round-off moves rays across the split, merge and removal thresholds, which
reshuffles the rays' per-cell sampling noise in Γ (the noise of a single
substep that `docs/otvet_closure.md` finds in the rays' snapshot Γ); the
RT backend alone, whose densities are bitwise the serial ones, stays at
1.6e-4. The test takes 95 s on 4 cores (serial 38 s, backends 53 s), well
inside the nightly budget.

The existing nightly `amrex_full64_vs_serial` had stopped working: it read
τ_es with the frame layout of before `be49827` (the Γ columns of
`frames.txt`) and wrote into a work directory it did not create. Both are
fixed here; it passes on klone (x_HII 4e-6, T 8e-5, τ_es 8e-7, fHI 6e-3,
T 1e-3 per cell, 105 s).

## 4. Timing

`pm_driver`, 3 steps of a Zel'dovich problem (`tools/ics/make_ics.py
--kind zeldovich`, Ngas = Ndm, pm_refine = 1), the `particle_mesh` wall
time per step on the 4-core VM (OpenMPI 4.1.6, `--bind-to none`), serial
PM (`src/pm.cpp`, OpenMP) against the backend:

| `particle_mesh` per step [s] | 128³ (2.1 M particles) | 256³ (16.8 M) |
|---|---|---|
| serial PM, 4 threads | 0.39 | 2.74 |
| serial PM, 1 thread | 0.83 | 7.16 |
| backend, 1 rank × 4 threads, 32³ boxes | 0.44 | 2.74 |
| backend, 2 ranks × 2 threads, 32³ boxes | 0.61 | 4.55 |
| backend, 4 ranks × 1 thread, 32³ boxes | 1.02 | 8.14 |
| backend, 4 × 1, 16³ boxes | 1.04 | 7.60 |
| backend, 4 × 1, 64³ boxes | 1.02 | 7.43 |
| backend, 1 rank × 1 thread, 32³ boxes | 0.85 | 7.07 |
| backend, 4 × 1, 32³ boxes, `amrex_pm_fft = amrex` | 0.91 | 7.14 |

(`pm_driver`'s outputs are 1.3 GB per step at 256³; the runs were timed
with the disk free, after a first set that filled it and doubled the
copy-out times.)

On one rank the backend matches the serial OpenMP PM at 256³ and costs
10 % more at 128³; on one thread it matches the serial code. With more
ranks it gets slower, not faster: every rank walks the full host particle
array to pick its own particles (`prep`, about 1.5 s per step at 256³) and
the gathers all-reduce 340 MB per step (the DM velocities and the gas),
so at 4 ranks the copy in and out costs more than the whole serial step.
The compute phases themselves (deposit, FFT, force, kicks) scale with the
ranks: the backend's own timers (3 steps, 32³ boxes; a repeat of
the 1 × 4 run at 256³ took 11.1 s against the 8.2 s of the table, so the
numbers on this shared VM scatter by 30 %):

| per 3 steps [s] | prep (locate, fill) | deposit | FFT | force | DM kick | gas kick | gather (v write, reduce, copy; gas; potential) |
|---|---|---|---|---|---|---|---|
| 128³, 1 × 4 | 0.19 (0.10, 0.06) | 0.37 | 0.16 | 0.13 | 0.31 | 0.02 | 0.07 (0.03, 0, 0.01; 0.02; 0.01) |
| 128³, 4 × 1 | 0.35 (0.27, 0.05) | 0.28 | 0.08 | 0.08 | 0.28 | 0.01 | 0.22 (0.04, 0.06, 0.03; 0.07; 0.02) |
| 256³, 1 × 4 | 1.63 (0.87, 0.43) | 2.77 | 1.14 | 0.99 | 2.66 | 0.17 | 0.69 (0.39, 0, 0.09; 0.16; 0.05) |
| 256³, 4 × 1 | 3.35 (2.24, 0.72) | 2.15 | 0.80 | 0.73 | 1.92 | 0.10 | 3.38 (0.96, 1.21, 0.27; 0.70; 0.24) |

The deposit and the DM kick are the particle loops (one CIC per particle
each); the FFT and the force are 10–15 % of a step. The potential gather
is `pm_driver`'s only (`gather_potential`).

This is the price of Stage E2's MPI model (the full state on every rank),
not of the level. When the particles stay on the level (the drifts and a
`Redistribute` there, the halo finder and the checkpoint reading the
container), `prep` disappears and the gather shrinks to the gas kick; the
distributed hydro of E1 removes the latter too. `amrex_pm_fft = amrex` is
10 % faster than the pencil transform at 4 ranks (its `plan_many`
transforms and fewer copies).

### 4.1 klone: pm_refine = 2 with 1024³ particles (2026-10-13)

The `hmf512_40mpc_r2` configuration (`docs/pm_refine.md`: 40 Mpc/h, 512³
gas, a 1024³ PM mesh and 1024³ particles, `pm_refine = 2`, unsplit hydro,
no sources, `tvir_min = 1e30`, FoF catalogues every step) restarted from
its last checkpoint (step 55, z = 5.50) for 3 hydro steps (56–58, z = 5.41
→ 5.26), `zdata1 = "15 5.5 5.35 4"` so that step 57 writes P(k) and the
particles, the gas fields (`u_gas`, `e_gas`, `rhog`) every step, and
`--set step_timers=1` (the backends' cumulative phase timers and the FoF
and P(k) times after every step, rank 0). The three backends with the
state on the levels (`rt_backend = hydro_backend = pm_backend = amrex`,
`amrex_state = 1`) on 1, 2 and 4 cpu-g2 nodes, 24 ranks × 8 threads per
node, the UCX settings of `env.sh`, for `amrex_pm_fft = pencil | amrex`
and gas boxes `amrex_max_grid_size = 32 | 64` (the PM mesh boxes are
`pm_refine` times larger, 64³ and 128³); the serial code on one node with
192 threads. Script `tools/hyak/amrex_pm_scaling.slurm`, tables
`tools/hyak/amrex_pm_scaling_table.py`, comparison with the serial code
`tools/hyak/amrex_pm_vs_serial.py`; outputs in
`$RUNS/hmf512_40mpc_r2/amrex_pm_scaling` (jobs 41832848 serial and the
1-node pencil runs, 41835129 the 1-node amrex runs, 41835130 2 nodes,
41835131 4 nodes; logs in `logs/stage_e_pm/`).

**Bitwise.** On each FFT, every run is identical to the 1-node, 32³-box
run of that FFT: every file of the three steps (`tests/amrex_pm_check.py
bitwise`: `u_gas`, `e_gas`, `rhog`, `stats.txt`), the particles, P(k) and
the three halo catalogues (sha256). That holds for 1, 2 and 4 nodes and
both box sizes: 5 comparisons per FFT, all identical. The two FFTs differ from each other
at round-off, as designed (§1.2). The amrex FFT itself is
rank-, node- and box-independent here as well (it also was in every test
of §3). `pencil` stays the default, since its independence holds by
construction.

**Against the serial code** (step 56–58; `vs_serial_n1_*.txt`; the pencil
and amrex runs give the same magnitudes):

| | pencil, 1 node | amrex, 1 node |
|---|---|---|
| `u_gas` rms diff / rms (ρ, ρv, E), step 58 | 8e-7 – 3.4e-6 | 7e-7 – 3.4e-6 |
| `u_gas` max diff / max, all steps | ≤ 4.1e-5 | ≤ 4.1e-5 |
| `e_gas` rms, step 58 | 3.5e-6 | 3.4e-6 |
| `stats.txt` (DM a/v max and rms, gas c/v) | ≤ 1.4e-6 | ≤ 1.3e-6 |
| N_halo, steps 56/57/58 (serial 452785 / 462732 / 471867) | 452786 / 462733 / 471863 | 452785 / 462734 / 471860 |
| particles at z = 5.33: positions rms / max [DM cells] | 3.2e-6 / 1.2e-4 | 3.2e-6 / 1.8e-4 |
| velocities rms / max [rms velocity] | 4.4e-7 / 3.3e-5 | 4.4e-7 / 3.4e-5 |
| P(k), largest relative difference | 1.0e-6 | 1.0e-6 |

Round-off, as in §3: the fixed-point CIC and the FFT against the serial
float CIC and the 3-D FFTW plan. A few halos out of 4.6·10⁵ near the
linking length change their membership.

**Per step** (mean over the 3 steps; PM = the `pm` column of `timing.txt`:
drift, deposit, FFT, force, kicks; FoF = the halo finder on rank 0,
including the gather of the particles into the host mirror it reads;
P(k) once, at step 57; summaries + io = the statistics and the field
output of every step, with the state gathered from the levels):

| [s per step] | PM | hydro | FoF | P(k) (step 57) | summaries + io | step |
|---|---|---|---|---|---|---|
| serial, 1 node, 192 threads | 11.18 | 15.89 | 32.1 | 20.0 | 7.6 | 73.5 |
| backends, 1 node, pencil, 32 | 4.13 | 3.12 | 188.0 | 48.4 | 20.4 | 231.8 |
| backends, 2 nodes, pencil, 32 | 2.97 | 2.08 | 183.8 | 48.6 | 21.4 | 226.4 |
| backends, 4 nodes, pencil, 32 | 1.32 | 0.87 | 182.0 | 49.9 | 22.8 | 223.6 |

The PM backend's phases (mean per step, s; `gather` = the particle mirror
for the FoF, counted in the FoF column above):

| nodes, FFT, box | PM | drift | deposit | (SumBoundary) | FFT | force | kick_dm | kick_gas | gather |
|---|---|---|---|---|---|---|---|---|---|
| 1, pencil, 32 | 4.13 | 1.07 | 1.04 | 0.29 | 0.78 | 0.23 | 0.95 | 0.04 | 8.20 |
| 1, pencil, 64 | 4.30 | 1.08 | 1.20 | 0.01 | 0.75 | 0.22 | 0.98 | 0.04 | 7.88 |
| 1, amrex, 32 | 4.50 | 1.07 | 1.01 | 0.19 | 1.20 | 0.22 | 0.92 | 0.04 | 7.62 |
| 1, amrex, 64 | 4.55 | 1.06 | 1.11 | 0.04 | 1.15 | 0.21 | 0.96 | 0.04 | 7.34 |
| 2, pencil, 32 | 2.97 | 0.50 | 0.59 | 0.12 | 0.92 | 0.11 | 0.81 | 0.03 | 6.14 |
| 2, pencil, 64 | 3.06 | 0.61 | 0.66 | 0.10 | 1.08 | 0.10 | 0.58 | 0.02 | 5.87 |
| 2, amrex, 32 | 2.97 | 0.53 | 0.57 | 0.21 | 1.22 | 0.11 | 0.50 | 0.02 | 6.13 |
| 2, amrex, 64 | 3.24 | 0.60 | 0.70 | 0.15 | 1.27 | 0.10 | 0.53 | 0.02 | 5.87 |
| 4, pencil, 32 | 1.32 | 0.35 | 0.27 | 0.07 | 0.31 | 0.06 | 0.33 | 0.01 | 4.92 |
| 4, pencil, 64 | 1.54 | 0.40 | 0.40 | 0.04 | 0.30 | 0.06 | 0.37 | 0.01 | 4.70 |
| 4, amrex, 32 | 1.41 | 0.36 | 0.27 | 0.08 | 0.39 | 0.06 | 0.31 | 0.01 | 4.91 |
| 4, amrex, 64 | 1.63 | 0.41 | 0.40 | 0.06 | 0.39 | 0.06 | 0.36 | 0.01 | 4.70 |

Peak RSS (`peak RSS` lines of the log and `sacct`): the serial code 168 GB.
The backends: rank 0 120–126 GB on every node count (its host state: the
mirror of the particles and the gas, and the FoF's and P(k)'s work
arrays, all of which the serial code also holds); the other ranks 12 GB
on 1 node and 7 GB on 2 and 4. The peaks of the job steps (sacct MaxRSS,
always rank 0) are 126.0–127.5 GB on 2 nodes and 128.5–131.4 GB on 4.

**Reading.**
- The PM backend is 2.7× faster than the serial PM on one node (4.1 s
  against 11.2 s), and 3.1× faster again on 4 nodes (1.3 s; 3.0 s on 2
  nodes, where the FFT's transposes first cross the network and it is
  slower than on one node). The deposit and the kicks scale with the
  ranks. The drift (positions on the level, `amrex_state = 1`) costs as
  much as the deposit. The FFT is 0.3–0.4 s on 4 nodes. `SumBoundary` is
  under 0.3 s, and 0.01–0.15 s with 128³ mesh boxes.
  The larger boxes are 4 % slower on one node and 15–17 % slower on 4
  nodes (512 mesh boxes for 96 ranks: worse balance in the deposit and the
  drift).
- `amrex_pm_fft = amrex` costs 0.4 s more than the pencil on one node (1.2
  s against 0.8 s: its plans are single-threaded on klone, §1.2), 0.2–0.3 s
  on 2 and 0.1 s on 4. It brings nothing here, so `pencil` stays the default.
- The hydro backend is 5× faster than the serial hydro on one node and 18×
  on 4.
- The step is 3× *slower* than the serial code: the FoF and P(k) run on
  rank 0 alone, with that rank's 8 threads, on a host mirror of all 1024³
  particles (188 s against 32 s on 192 threads, and 48 s against 20 s for
  P(k)). They take 90 % of the step and do not change with the node count.
  The gather into the mirror also scales poorly (8.2 s on 1 node, 4.9 s on
  4). This is the open item of `docs/stage_e_state.md` (the halo finder on
  the distributed particles, E3). Until then, a run with FoF every step
  should use `fof_dt_myr` or run the halo finder less often. Rank 0 with
  more threads would also help (a heterogeneous launch).
- Rank 0's 124 GB is the host state that the FoF, P(k) and the output
  read. It does not shrink with the node count. The other ranks hold
  their part of the levels.

**Two fixes found by these runs** (neither changes a result):
- *Start-up memory.* With `amrex_state = 1`, every rank allocated the
  host PM arrays and ran the host initial FFT, then freed them: 45 GB per
  rank at 1024³, so 24 ranks went OOM at 700 GB on one node.
  `pm::init_pm(..., host_arrays)` now allocates them on the ranks that
  need them (rank 0, which reads the checkpoint and does the I/O), and the
  others skip the host FFT. The state after start-up is the one the ranks
  had after freeing the arrays.
- *The amrex FFT's exit crash* (§5, the planner reset): every
  `amrex_pm_fft = amrex` run at 1024³ finished its steps and then crashed
  destroying its plans. The step outputs of the crashed run are
  byte-identical to the fixed one.

## 5. Findings and traps

- The serial `link_dm` cell of a particle (`int(x·x_dm_m)`) and the CIC
  corners (`1 + int(x·x_dm_m − 0.5)`, …) differ by at most one cell, so one
  ghost layer suffices for the deposit and the kick; every ghost access is
  range-checked and the step aborts if a corner falls outside.
- The AMReX position of a particle is the centre of its home cell, not the
  physical position: `Redistribute` is never called, and the float
  position in DM-grid units stays exact in the SoA components (a
  conversion to mesh units would not be exact when Nm/Ndm is not a power
  of two).
- `amrex::FFT` sets FFTW's planner to the OpenMP thread count when its
  first solver is built (`fftwf_plan_with_nthreads`), a global that the
  host's later plans would inherit; the backend resets it to 1 after each
  solve, but only when it is above 1: in FFTW 3.3.10,
  `fftwf_plan_with_nthreads` before `fftwf_init_threads` calls
  `fftwf_cleanup`, and an AMReX built with `AMReX_FFTW_IGNORE_OMP` (klone)
  never initialises the threads. The unguarded reset left the R2C plans
  undefined; they still ran, but destroying them segfaulted at the end of
  the 1024³ runs (§4.1; the 32³ tests use transforms without twiddle
  tables and never showed it).
- The bitwise-OR all-reduce gathers disjoint float pieces exactly; a
  summing all-reduce would turn −0.0 into +0.0.
- `pm_driver`'s `phi.bin` is the potential on the PM mesh: the backend
  gathers it only when `gather_potential` is set (`radhydro` never needs
  it).
- The particle order within a tile depends on the thread count; nothing
  depends on it (integer sums, maxima, per-particle kicks, writes by id).

## 6. What is left for Stage E

- The distributed hydro (E1, the separate thread) and then the end of the
  redundant host state: the level would keep the particles between steps
  (drift and `Redistribute` on the level, the gather replaced by the
  checkpoint and the halo finder reading the container).
- The halo finder and sources on the distributed particles (E3), load
  balancing (E4), AMReX checkpoints (E5), the klone scaling (E6; the PM
  part is §4.1).
- The halo finder and P(k) on the distributed particles: on rank 0 they
  are 90 % of a step at 1024³ and make the full loop 3× slower than the
  serial code (§4.1).
- The GPU build (the loops are plain OpenMP here; the level's structure
  follows `rt_level.cpp`, whose device path is in place).

## 7. How to run

```
# AMReX 26.10 with the FFT module (needs FFTW3 double and float on the CPU)
cmake -S amrex -B build-amrex-lib -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$HOME/ext/inst-amrex \
  -DAMReX_SPACEDIM=3 -DAMReX_MPI=ON -DAMReX_OMP=ON -DAMReX_PARTICLES=ON \
  -DAMReX_PARTICLES_PRECISION=SINGLE -DAMReX_FORTRAN=OFF -DAMReX_FASTMATH=OFF -DAMReX_FFT=ON
cmake --build build-amrex-lib -j4 && cmake --install build-amrex-lib
cmake -S . -B build-amrex -DCMAKE_BUILD_TYPE=Release -DRADRAY_AMREX=ON \
  -DAMReX_DIR=$HOME/ext/inst-amrex/lib/cmake/AMReX
cmake --build build-amrex -j4
ctest --test-dir build-amrex -L amrex_pm -j4
# runs
mpirun -np 4 -x OMP_NUM_THREADS=1 build-amrex/pm_driver tests/problems_pm/zeldovich_3d out \
  --set pm_backend=amrex [--set amrex_pm_max_grid_size=16] [--set amrex_pm_fft=amrex]
mpirun -np 4 build-amrex/radhydro RUN out --set legacy_rt=0 --set rt_backend=amrex --set pm_backend=amrex \
  --set deterministic_parallel=0
```
