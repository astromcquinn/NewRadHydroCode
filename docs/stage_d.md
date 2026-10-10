# Stage D: the ray-tracing RT on AMReX (uniform grid, CPU first)

*Written 2026-10-10 on branch `stage-d-amrex`. Stage D of
`docs/ROADMAP.md` §7 and of the migration plan in
`docs/framework_decision.md` §6: the sorted-ray-array RT of the framework
spike (`spike/`) turned into a library component on an AMReX level, with
the ray kernels shared with the serial port and the compact R-struct ray
encoding of `docs/ray_compression.md`. The serial port stays the bit-for-bit
legacy reference; the AMReX path is a new option, `rt_backend = amrex`,
implemented for `legacy_rt = 0` only.*

## 0. Summary

| Milestone | State | Validation |
|---|---|---|
| 1. One physics implementation (`include/radray/kernels/`) | done | fast suite 100/100 and the nightly label unchanged, every flag combination byte-identical to the stored and live gfortran references |
| 2. Library component `src/amrex/rt_level.cpp`, `rt_backend = amrex` in `rt_driver` and `radhydro` | done | bitwise identical to `spike/serial rt_serial` for 1/2/4/8 ranks × 1/2 threads × box sizes 32/16/8; round-off against `rt_driver --set legacy_rt=0`; rank- and thread-independent full loop; bit-identical restarts (checkpoint version 5) |
| 3. R-struct encoding, `ray_encoding = rstruct` | done | 28 B of ray data per ray (80 B float); x_HII, τ_es and T errors in §3 |
| 4. Tests (`ctest -L amrex`), timing, klone job | done here; the klone runs are the user's | §4 |
| 5. GPU | scaffolding only: CUDA CMake target, device rate tables, build and `ckpt` job scripts; **not compiled** (no nvcc ≥ 12.2 in the cloud container) | §5 |

Ground rules kept: the serial code is untouched in behaviour (the kernels
it now calls are its own arithmetic moved into headers); every stored-
reference, live-Fortran, restart and thread test passes unchanged; AMReX is
optional (`-DRADRAY_AMREX=ON`, `find_package` or `-DRADRAY_FETCH_AMREX=ON`);
the default build and the fast suite need neither AMReX nor MPI.

## 1. Milestone 1: one implementation of the ray and RT-cell physics

`spike/common`'s kernels moved to `include/radray/kernels/` (as PR #15 did
for `hydro_mh.hpp`) and became the single implementation:

| Header | Contents | Serial callers |
|---|---|---|
| `ray_kernels.hpp` | `rrk::RayP` (the 80-byte payload), `RTP` (substep parameters), `ray_segment`, `cast_1`, `cast_3`, `split`, `merge_add`/`merge_finish`, source emission (`SourceSpec`, `source_ray`, `source_thick_deposit`), `classify`, `merge_cone`, the photon-count histogram bins, the sorted-array `step1_count`/`step1_write`, the redshift table and `redshift_spectrum`/`redshift_ray`, `interp_lin` | `src/ray_tracing.cpp` (`ray_source`, `ray_cast_split`, `ray_cast_merge`, `ray_cast_1/3`, `ray_segment`, the nside/radius helpers, `ray_redshifting`, `interpolate_lin`, `make_redshift_table`), `src/recomb_rad.cpp` (`ray_source_recomb`), `src/rt_moments.cpp` (`redshift_spectrum`) |
| `rt_cell.hpp` | `CellData` (46 floats: rt_grid fields, u_gas, w_gas, e) and `CellRef` (pointer view, so the serial SoA grid and the backend's AoS cell share the code), `photoion_sigma`, the rate lookups (`rates_node`, `rate_lerp`, `rates_interp`), `gas_state_from_avg`/`store_gas_state`/`tau_from_state`, `calc_rt_1_cell`, and `radiative_physics` in three parts (`cell_radiation_prep`, `photo_gas_cell`, `photo_cell_finish`) plus `physics_cell` for N_rt_gas = 1 | `src/rt_tools.cpp` (`prep_rt`, `rt_step_1`'s `calc_rt_1`, `radiative_physics_plane`, `post_rt`), `src/rt_rates.cpp` |
| `healpix_hd.hpp` | table-free NESTED `pix2vec`/`vec2pix` (same libm calls in the same order as HEALPix 3.11) | `src/healpix.cpp` (adds the range checks) |
| `consts_hd.hpp`, `rr_hd.hpp` | device-usable constants, Fortran-semantics helpers | |

What the kernels gained over the spike's copies, so that the serial code
could call them without changing a bit: the `legacy_rt` switch of
`calc_rt_1` (B3, B4) and the segment count (B8, `RTP::nseg`), the merge-
cone nside that `recomb_prepare` lowers (`RTP::nside_merge`), the rt_stat
sums (`PhotoStats`: nion, ncol, nrec, nrab, gamma, lambda, dtau, each
returned per gas cell and added by the host in its own order), case A
recombination (`PhysCtx::ia`) and the ground-state recombination count,
the dual-energy field of `hydro_scheme = unsplit`, the temperature-index
error of the legacy rates as a code (the host throws), the clamp flag of a
split (ray-dump provenance), source emission with a given spectrum
(recombination rays). `radray::Ray` is now `rrk::RayP` plus the two chain
links (the links moved to the end; checkpoint versions ≤ 4 are converted
on reading). The ngam_ray thresholds and the F3 pool guard became shared
host functions (`calc_ngam_ray`, `pool_split_candidates`, `pool_guard` in
`rt.hpp`).

One trap found on the way: the serial `photo_source` statistic adds each
frequency group's photons to the `rt_stat%ngam` accumulator separately;
summing the groups first changes the last bits of a diagnostic sum.
`source_stats` therefore returns the per-group addends.

Not carried: `ray_rotate` (branch `claude/exciting-clarke-zn62ih` is not on
the remote; nothing to merge). The moment solvers (`rt_method = m1 |
otvet`) stay on the serial path.

`spike/common`'s copies of the kernels are deleted. `spike/serial
rt_serial` (the bitwise reference of the backend) builds against the
library kernels and now applies the pool guard; `spike/amrex` and
`spike/parthenon` are frozen prototypes that no longer compile against the
moved kernels (their maintained successor is `src/amrex/rt_level.cpp`).

## 2. Milestone 2: the library component

### 2.1 Where things live

| File | What |
|---|---|
| `include/radray/rt_backend.hpp`, `src/rt_backend.cpp` | the switch and its parameters (AMReX-free): `rt_backend = serial \| amrex`, `ray_encoding = float \| rstruct`, `amrex_max_grid_size` (16), `amrex_ghosts` (2), `amrex_redist_local` (3; 0 = global Redistribute); `check_rt_backend_param` refuses `legacy_rt = 1`, `rt_method ≠ rays`, `recomb_rad`, `deterministic_parallel`; the process functions (`backend_init/finalize`, `io_rank`, …) |
| `src/amrex/rt_level.cpp` | the level (C++20; target `radray_amrex`, links `AMReX::amrex`) |
| `src/amrex/ray_encoding.hpp` | `FloatEnc` and `RStructEnc`: ray ↔ particle components |
| `src/amrex/rt_level_stub.cpp` | the no-AMReX build: the level functions throw |
| `src/rt_tools.cpp` | `radiative_transfer` dispatches: `prep_rt` → `amrex_rt_prep`, each substep → `amrex_rt_substep`, `post_rt` → `amrex_rt_post`; `accumulate_rt_stats` (the end of `rt_step_3`) is shared |
| `src/checkpoint.cpp` | version 5: the backend's rays gathered on rank 0 |
| `apps/rt_driver.cpp`, `apps/radhydro.cpp` | `backend_init` at start, rank 0 writes the files, the movie hook syncs the level's fields before a frame |

### 2.2 The level

One uniform AMReX level of the RT grid (periodic `Geometry` with dx = 1,
`BoxArray::maxSize(amrex_max_grid_size)`, the default `DistributionMapping`):

- cells: `FabArray<BaseFab<rrk::CellData>>` (46 floats per cell, AoS, no
  ghosts: the per-cell kernels take the struct; the spike measured the SoA
  layout 3–5× slower on the CPU), `FabArray<BaseFab<float>>` of tau(5) and
  fgam(5) with `amrex_ghosts` ghost layers (`FillBoundary` before `cast_1`
  and `cast_3`);
- deposits: `FabArray<BaseFab<Long>>` of the 13 int64 fixed-point
  accumulators with the same ghosts, atomic adds into the box and its ghosts,
  then `SumBoundary` (exact, associative);
- rays: `ParticleContainerPureSoA<kNR, kNI>` (float build of AMReX;
  position = the centre of the home cell, as the periodic image nearest to
  the emitting cell so the local `Redistribute` sees a small move); gathered
  per box into a cell-sorted `rrk::RayP` scratch array with `DenseBins`,
  migrated by `Redistribute` after every substep;
- the reach assertion: every ghost access is range-checked and the run
  aborts with "a ray reached outside the ghost region" if a target is
  beyond `amrex_ghosts` (2 is the reach of one substep, design doc §3.2; 4
  also validated).

One substep, in the order of the serial `rt_step_1/2/3`:
1. `calc_rt_1_cell` on every cell (legacy_rt = 0: B3, B4), tau into the
   ghost array;
2. the photon-count histograms of the live rays (`ngam_pdf`, and `spl_pdf`
   of the rays that may split), the live count, the exact photon sum
   (128-bit fixed point) → `ReduceLongSum` → `calc_ngam_ray` and
   `pool_guard` (the F3 guard, with the number of stellar source cells);
3. per box: `DenseBins` sort by home cell, then per cell `step1_count`
   (emission, classification, cones), an exclusive scan, `step1_write`
   (source rays, casts, 4 daughters per split, merges in the canonical
   order), the thick-source deposit; per-cell maxima of casts/splits/merges
   for `ray_stat%ncst2` etc.;
4. `cast_1` of every output ray (deposits, ≤ 4 segments), `SumBoundary`;
5. `physics_cell` on every cell (the three parts of `radiative_physics`
   with the dual-energy field when present), the per-cell `PhotoStats`
   folded serially per box into 128-bit fixed-point sums, fgam into the
   ghost array;
6. `cast_3` (advance, attenuate, delete), compaction, the survivors
   written back as particles, `Redistribute`;
7. the counters and sums reduced (`ReduceLongSum`, 128-bit sums as four
   32-bit limbs) into `s.rays.stat`, `s.rays.Nray` (the live count),
   `s.rt_stat`, then `accumulate_rt_stats`, so `summarize_rt`, the ray.log
   lines and the history outputs are the serial ones.

`amrex_rt_prep` copies `gas.u/w/e` and `rt.Luv/sfr` of the host grids into
the cells (every rank holds the full host state) and computes the RT grid
gas state; `amrex_rt_post` redshifts the rays (`redshift_ray`), computes
`calc_rt_f`, the live count and the exact photon total
(`ngam_rays_total`), and copies the cells back: `ParallelCopy` to a
one-box FabArray on rank 0, then `Bcast` of the flat 46 × N³ floats. The
movie hook calls the same copy (`amrex_rt_sync_fields`) only when a frame
is due.

**MPI model of Stage D.** `radhydro`'s hydro, PM, halo finder and source
model run redundantly on every rank on the full serial state (they are
deterministic, so every rank has the same state); only the RT is
distributed. Rank 0 writes the files. This costs one broadcast of the cell
fields per hydro step (0.4 GB at 128³, 3 GB at 256³) and the memory of the
full state per rank (about 3.5 GB at 256³ with Ndm = N); it is what
Stage E replaces.

### 2.3 Determinism

Every cross-cell or cross-rank sum is an integer: int64 deposits,
integer histograms and counters, and 128-bit fixed-point sums of the
double statistics (per quantity a power-of-two scale: 2⁻¹⁰⁰ for the
per-atom rates, 2⁻¹²⁰ for Γ, 2⁻¹⁶⁰ for Λ, 2⁻¹¹⁰ for dτ, 2⁻⁴⁰ for photons per
ray; reduced with `ReduceLongSum` as 32-bit limbs). Merges use the
canonical order (cone, payload hash, payload bytes). The result is
therefore bitwise identical for any rank count, box size and thread
count, and identical to `rt_serial`; the statistics (printed as floats)
agree with the serial port to the printed precision.

### 2.4 What the backend does not do (yet)

- `recomb_rad = 1` (the recombination-ray emission exists in the kernels,
  `source_spec_spectrum`; the block bookkeeping of `recomb_prepare` is not
  on the level);
- the ray dumps and `ray_quant` research instruments;
- `legacy_rt = 1` (impossible by design, design doc §0.1);
- distributed hydro/PM (Stage E).

## 3. Milestone 3: the R-struct encoding

`ray_encoding = rstruct` stores each ray in 7 int32 particle components
(28 B of ray data; the float layout uses 18 float + 2 int components, 80 B).
AMReX adds its 3 float positions and the 8-byte id/cpu to either: 48 B
against 100 B per particle in the container. Bytes per ray are reported by
`amrex_rt_bytes_per_ray`.

| Field | Bits | Encoding |
|---|---|---|
| position | 3 × 16 | offset in the home cell (1/65536 cell); the cell is the particle position |
| nside, ipix | 30 + 1 | the unique id 4 nside² + ipix; bit 30: on-centre (never merged: the direction is the pixel centre, no offset stored) |
| direction offset | 8 + 8 | (θ, φ) tangent-plane offset from the pixel centre in units of pixel radius / 128, merged rays only |
| r | 16 | 1/8 cell; dr = rmax(nside) − r |
| photon total | 32 | float |
| spectral shape | 5 × 8 | −16 log₂(ngam_g / total), 1/16 octave over 16 octaves (255 = zero); decoded as total × f_g / Σf |
| nu | 5 × 8 | linear in the group's band |

Every lossy field is rounded with deterministic dithering (`floor(v + u)`,
u from a hash of the substep, the ray's float payload and the field),
which is unbiased and independent of the thread and rank layout. The
on-centre test is |dx − float(centre)| ≤ 2·10⁻⁶ (float rounding of
never-merged directions).

Accuracy against the float layout (`tools/ray_compression/precision_metrics.py`,
targets of `docs/ROADMAP.md`: x_HII within 1–2 %, τ_es within 1 %, T within
a few percent, I-front within half a cell):

| Run | x_HII,V | x_HII,M | T_M | T per ionised cell (rms) | τ_es | front shift [cells] |
|---|---|---|---|---|---|---|
| `lognormal_16src`, 10 steps, 32³ | 1.3e-4 | 5.0e-6 | 2.8e-5 | 1.2e-3 | 3.4e-6 | 4.7e-4 |
| `lognormal_16src`, 100 steps, factor 16 (the study's setting; x_HII,V = 0.95 at the end) | 5.2e-4 | 9.2e-5 | 5.9e-4 | 1.3e-2 | 4.4e-6 | 2.7e-2 |
| `uniform_8src`, 100 steps, factor 16 | 2.9e-5 | 2.9e-5 | 1.0e-3 | 7.8e-3 | 2.9e-6 | 8.3e-3 |
| 64³ production configuration (10 Mpc/h, 0.01c, z = 100 → 5.5; §4.2) | 1.2e-4 | 1.2e-4 | 2.9e-4 | 1.1e-2 | < 1e-4 (0.05224 in both) | 2.8e-2 |

All within the acceptance metrics, with a margin of 10–100 on x_HII and
τ_es; the per-cell temperature scatter of 1 % is the 1/16-octave shape
quantisation averaging down over the rays of a cell, as the compression
study predicted for dithered 8-bit photons. The decode (`pix2vec_nest`
per ray) made the first R-struct run of the 64³ configuration about twice
as slow as the float one; the substep now decodes each ray once (the
histograms read the gathered scratch array), see §4.3.

## 4. Validation, tests and timing

### 4.1 `ctest -L amrex`

Registered only with `RADRAY_AMREX=ON` (43 tests; the MPI launcher is
`RADRAY_MPIEXEC` with `RADRAY_MPIEXEC_FLAGS`, OpenMPI defaults; `srun` and
empty flags on a cluster):

- `amrex_bitwise_<case>_r<ranks>t<threads>_mgs<box>`: `rt_driver --set
  rt_backend=amrex` against `rt_serial` with `compare.py --rtol 0
  --atol-frac 0` (every field of every step and all stats), for the three
  cases of `spike/amrex/validate.sh` (`uniform_8src`,
  `lognormal_16src_smallbuf`, the latter with `legacy_physics = 0`) and
  (ranks, threads, max_grid_size) ∈ {(1,1,32), (1,2,16), (2,1,16), (2,2,8),
  (4,1,16), (4,2,8), (8,1,8), (8,2,16)};
- `amrex_roundoff_<problem>`: the five RT problems against `rt_driver
  --set legacy_rt=0` (rtol 1e-5; dx of source cells and fgam of I-front
  cells excepted, as explained in `docs/framework_decision.md` §2);
- `amrex_thread_lognormal_16src`: 1, 2 and 4 threads byte-identical;
- `amrex_rstruct_lognormal_16src`: the R-struct run stays within 10 % of
  the serial fields;
- `amrex_full_ranks_smoke`: `radhydro` on the full-loop smoke problem with
  1 and 2 ranks, every file (fields, movie frames, checkpoints)
  byte-identical;
- `amrex_restart_smoke`, `amrex_restart_smoke_unsplit`: bit-identical
  restarts from the version-5 checkpoint, legacy and unsplit hydro;
- nightly: `amrex_full64_vs_serial`, the 64³ production configuration
  against the serial code (history, τ_es, per-cell fields).

Results in this session (cloud VM, 4 cores, OpenMPI 4.1.6): every test
passes (the 27 bitwise tests, the 5 round-off, threads, R-struct, the
full-loop rank check and the three restarts); the serial code's fast suite
(100 tests) and the whole nightly label (20 tests) pass unchanged with the
moved kernels. The three round-off tests needed the 1e-3 field criterion
(single I-front or removal-threshold cells at 1e-4; `fHeII` 2.5e-4 in one
cell of `lognormal_16src_smallbuf`, where the pool guard and merges are
active) and 1e-6 on the float stats (`ngam_rays_total`, a double sum in a
different order). The bitwise block: see the run log in the PR.

### 4.2 The 64³ production configuration

`make_production.py --N 64` (10 Mpc/h, 0.01c, factor 16, on-the-fly
sources, legacy_rt = legacy_physics = 0): serial (deterministic_parallel),
backend with float rays and backend with R-struct rays, 2 threads each.

| Run (2 threads, 1 rank, loaded machine) | wall [s] | RT [s] | peak RSS [MB] | x_HII,V vs serial | T_M | T per ionised cell (rms) | front shift | τ_es |
|---|---|---|---|---|---|---|---|---|
| serial (`deterministic_parallel`) | 128 | 99 | 405 | | | | | 0.05224 |
| backend, float rays | 359 | 323 | 1396 | 7.3e-6 | 4.0e-5 | 1.3e-3 | 1.2e-4 cells | 0.05224 |
| backend, R-struct rays | 744 | 699 | 1252 | 1.2e-4 (vs float) | 2.9e-4 | 1.1e-2 | 2.8e-2 cells | 0.05224 |

(`tools/ray_compression/precision_metrics.py` on the 10 zdata1 outputs;
τ_es from `tools/analysis/reion.py` on the movie frames.) The backend's
history is the serial one to 1e-5: the differences are the round-off of
the deposit and merge order (design doc §2.2), which a reionisation run
amplifies only to 1e-3 per ionised cell over 1479 substeps. The wall
times here were taken with the nightly suite and the `amrex` label
running on the same 4 cores; the clean comparison is §4.3. The peak memory
of the backend is the serial state plus the level (cells 46 floats, tau
and fgam with ghosts, 13 int64 deposits with ghosts, the particle tiles
and per-box scratch) plus the flat gather buffer of `sync_fields`; the
R-struct run has the smaller ray tiles (48 B against 100 B per particle).

### 4.3 Timing

`rt_backend = amrex` against the serial OpenMP code on the 4-core VM at 64³
and 128³ (`tools/production/timing_table.py` constants from the measured
runs): TBD.

Klone: `tools/hyak/amrex_scaling.slurm` runs one hydro step of a 256³
production checkpoint with 24 ranks × 8 threads (one rank per CCD) and the
serial code on 192 threads for the comparison with `docs/cluster.md` §6
(the user runs it; `tools/hyak/build_amrex.sh` builds AMReX and the code on
a node).

## 5. GPU (not compiled here)

Prepared, compile-only pending nvcc ≥ 12.2 (the cloud container has none):

- CMake: with an AMReX built with `AMReX_GPU_BACKEND=CUDA`, `rt_level.cpp`
  is compiled as CUDA (`--fmad=false`, C++20, separable compilation);
- the rate tables are copied once into device memory
  (`Gpu::DeviceVector`, `rates_view`); the cell scratch, histograms and
  counters are device vectors with `HostDevice::Atomic`;
- every per-cell and per-ray loop is an `AMREX_GPU_DEVICE` lambda in
  `ParallelForOMP` (OpenMP inside a box on the CPU, a GPU kernel otherwise);
- `tools/hyak/build_amrex_gpu.sh` builds the AMReX CUDA library and the
  code on a `ckpt` GPU node (A100 sm_80, H200 sm_90);
  `tools/hyak/amrex_gpu.slurm` runs one hydro step of the 128³ benchmark
  for the comparison with `docs/framework_decision.md` §4.

Still host-only, to be done when a GPU build is possible: the per-box
serial folds (the 128-bit statistics, the exact photon sum, the per-cell
maxima) and the host-pointer copies of `prep`/`sync_fields` (pinned
buffers or managed memory); a team-level `step1` for source cells (192
emissions and a heap sort of merges per thread serialise); the AoS cell
FAB and ray scratch are uncoalesced.

## 6. Findings and traps

- AMReX's `ParallelCopy`/`FillBoundary` need `operator+=` on the element
  type: the AoS cell FAB is copied through a 46-component float FAB for the
  gather (spike gotcha 5, unchanged).
- `ParallelFor` is serial on the CPU; threads come from `ParallelForOMP`
  inside a box. The per-box serial folds are cheap (a few ms per substep at
  32³ per box) but are the first thing to move for a GPU.
- Particle precision is build-wide (SINGLE): the float build is the one
  that reproduces the original's `real(4)` rays and dark matter.
- The merge must use the canonical order, and every reduction an integer:
  with that, the backend was bitwise identical to `rt_serial` on its first
  complete run, for every rank layout tried.
- The nearest-image particle position (spike gotcha 2) is still needed for
  the local `Redistribute`; `amrex_redist_local = 0` (global) also passes.
- The multi-node hydro slowdown of `docs/framework_decision.md` §3.3 is
  untouched here (Stage E); nothing new learned about it on one machine.

## 7. How to build and run

```
# AMReX 26.10 (once; 3 min at -j4)
git clone --branch 26.10 --depth 1 https://github.com/AMReX-Codes/amrex.git amrex
cmake -S amrex -B build-amrex-lib -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$HOME/ext/inst-amrex \
  -DAMReX_SPACEDIM=3 -DAMReX_MPI=ON -DAMReX_OMP=ON -DAMReX_PARTICLES=ON \
  -DAMReX_PARTICLES_PRECISION=SINGLE -DAMReX_FORTRAN=OFF -DAMReX_FASTMATH=OFF
cmake --build build-amrex-lib -j4 && cmake --install build-amrex-lib
# the code with the backend (or -DRADRAY_FETCH_AMREX=ON instead of AMReX_DIR)
cmake -S . -B build-amrex -DCMAKE_BUILD_TYPE=Release -DRADRAY_AMREX=ON \
  -DAMReX_DIR=$HOME/ext/inst-amrex/lib/cmake/AMReX
cmake --build build-amrex -j4
ctest --test-dir build-amrex -L amrex -LE nightly -j4
# runs
mpirun -np 4 -x OMP_NUM_THREADS=2 build-amrex/rt_driver tests/problems/lognormal_16src out \
  --set legacy_rt=0 --set rt_backend=amrex [--set ray_encoding=rstruct] [--set amrex_max_grid_size=16]
mpirun -np 4 build-amrex/radhydro RUN out --set legacy_rt=0 --set rt_backend=amrex --set deterministic_parallel=0
```
