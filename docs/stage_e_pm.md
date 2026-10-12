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
at round-off. The amrex path resets FFTW's planner thread count to 1 after
each use, so the host code's own plans (the hydro's 1-D transforms, the
serial PM) stay the serial ones.

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
trajectories to 1e-5 over 30 steps, and the unit tests pin the forces. A
64³ production-configuration comparison (as `amrex_full64_vs_serial` does
for the RT) is the right nightly check and is not added here.

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
  solve.
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
  balancing (E4), AMReX checkpoints (E5), the klone scaling (E6).
- A 64³ production-configuration nightly comparison of the full loop with
  both backends against the serial code.
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
