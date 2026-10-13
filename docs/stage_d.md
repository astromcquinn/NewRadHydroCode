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
| 4. Tests (`ctest -L amrex`), timing, klone job | done; klone 256³ on 1, 2 and 4 nodes, 512³ on 2 and 4 nodes (2026-10-12) | §4, §4.4 |
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

Carried at the merge with `main` (2026-10-11): the B16 fix `ray_rotate`
(the frame rotation is `rrk::frame_rotation`, keyed by `RTP::ray_rotate`,
`rot_seed`, `rot_abits`, `rot_istep`; `SourceSpec` holds the rotation and
`source_ray` applies it, so the serial `ray_source`, the recombination rays
and the backend's emission share it; the backend rotates with the 1-based
cell like the serial code), `redshift_per_distance` (`redshift_fscale` in
`ray_redshifting` and the backend's `post_rt`) and `recomb_rad_lines`
(`PhotoStats::nrB`, the case-B recombinations). The moment solvers
(`rt_method = m1 | otvet`) stay on the serial path.

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

Known limit of the AMReX build: the live-Fortran full-loop tests
(`full_vs_fortran_live_*`) fail to link there, because the Fortran harness
links `libradray.a`, whose `rt_tools.cpp` now references the backend's
`amrex_rt_*` symbols from the separate AMReX library. Run the fast suite in
the default build and only `ctest -L amrex` in the AMReX build (as CI does);
all seven live tests pass in the default build of the merged branch.

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

`rt_backend = amrex` against the serial OpenMP code on the 4-core VM
(idle; OpenMPI 4.1.6, `--bind-to none`): two hydro steps restarted from
the last checkpoint of a production run (`make_production.py`, 10 Mpc/h,
0.01c, factor 16; 64³: z = 4.9, 17 substeps, 3.2 M live rays), the RT
time of `timing.txt` summed over the two steps. The serial code uses
`deterministic_parallel`; the backend `amrex_max_grid_size = 16` unless
noted. The serial checkpoint restarts on the backend directly (its live
rays go to the level).

| 64³, 2 hydro steps | RT [s] | peak RSS of rank 0 [MB] |
|---|---|---|
| serial, 4 threads | 23.5 | 401 |
| serial, 1 thread | 86.8 | 401 |
| backend, 1 rank × 4 threads | 26.9 | 1983 |
| backend, 2 ranks × 2 threads | 24.3 | 1363 |
| backend, 4 ranks × 1 thread | 24.4 | 1300 |
| backend, 4 × 1, `amrex_max_grid_size = 32` | 24.7 | 1206 |
| backend, 1 rank × 1 thread | 73.7 | 1818 |
| backend, 4 × 1, R-struct rays | 31.8 | 839 |

So on 4 cores the backend is within 4–15 % of the serial OpenMP code
whatever the rank/thread split, and on one core it is 15 % faster than the
serial code (the sorted arrays against the chains; the spike measured
1.3–1.5× slower, before the contention fix of §6 and with the AoS cell
FAB). The R-struct run costs 30 % more (the decode, one `pix2vec` per ray
per substep, plus a second one for the on-centre test at encode) and
halves the memory of the rays. The backend's memory is dominated by the
redundant serial state plus the level (§2.2); the 1-rank case also holds
the flat gather buffer.

| 128³, 2 hydro steps (z = 5.0, 13 substeps, 22.5 M live rays) | RT [s] | peak RSS of rank 0 [MB] |
|---|---|---|
| serial, 4 threads | 148.9 | 2817 (the whole run) |
| serial, 1 thread | 598.7 | |
| backend, 1 rank × 4 threads | 152.2 | 11480 |
| backend, 2 ranks × 2 threads | 136.8 | 10260 |
| backend, 4 ranks × 1 thread | 136.4 | 9787 |
| backend, 4 × 1, `amrex_max_grid_size = 32` | 128.9 | 9677 |
| backend, 1 rank × 1 thread | 390.9 | 11455 |
| backend, 4 × 1, R-struct rays | 172.1 | 6233 |

At 128³ the backend is 8–13 % faster than the serial OpenMP code on the
same 4 cores (4 ranks with 32³ boxes is best), and 1.5× faster on one
core; R-struct costs 26 % (decode) and saves 36 % of rank 0's peak memory.

The first attempt at the 128³ backend runs was killed by the container's
15 GB memory limit for every multi-rank case: each rank read the whole
checkpoint (the serial ray pool and the ray list, 4 GB) and the level held
three copies of the rays (the particle tile and two 80-byte scratch
arrays). Now ranks other than 0 skip the ray blocks of a checkpoint, and
the level releases each copy as soon as it is used (at most two exist at
once; AMReX's arena reuses the blocks). The peak RSS above is rank 0's and
is reached during the restart (the serial pool, the ray list and the
particle tile before the first `Redistribute`), not in the substeps; the
steady-state memory per ray of the level is the tile (100 B float, 48 B
R-struct) plus the transient scratch (80 B), against 88 B in the serial
pool. Reducing the restart peak (streaming the rays to their owners) and
the redundant serial state per rank are open (Stage E).

The RT cost model of `tools/production/timing_table.py` (`--a-cell`,
`--b-ray` per substep) is unchanged for the backend to within these
ratios; its constants are those of `docs/production_run.md` §5.

Klone: `tools/hyak/amrex_scaling.slurm` runs one hydro step of a
production checkpoint with 24 ranks × 8 threads per node (one rank per
CCD) and, on one node, the serial code on 192 threads for the comparison
with `docs/cluster.md` §6; `tools/hyak/build_amrex.sh` builds AMReX and
the code on a node. The results are in §4.4.

### 4.4 klone: 256³ and 512³ on 1, 2 and 4 nodes (2026-10-12)

Stage E prerequisites on `cpu-g2` (AMD EPYC 9654, 2 × 96 cores, 1.5 TB;
gcc 13.2, OpenMPI 4.1.6, AMReX 26.10 built by `tools/hyak/build_amrex.sh`
into `$RADRAY_SW/amrex-26.10`; `ctest -L amrex` 54/54 (81/81 with the 27 `amrex_hydro` tests of Stage E1) on a 32-core
allocation, launched with `mpirun --oversubscribe` because a one-task
allocation gives OpenMPI a single slot). One hydro step restarted from the
last checkpoint of a production run, `run_steps = 1`, no checkpoint or
movie output, the fields of the last step written and compared
(`tools/hyak/amrex_scaling_compare.py`). The backend runs 24 ranks × 8
threads per node at 256³ and 12 ranks × 16 threads per node at 512³
(memory, below), `amrex_max_grid_size = 32`, `OMP_PLACES=cores
OMP_PROC_BIND=close`, `--cpu-bind=cores`; the serial code 192 threads with
`deterministic_parallel`. "Peak RSS" is VmHWM of the process, reported per
rank by `radhydro` (rank 0, and the min/max over the ranks) and by
`sacct` (the mean over the ranks).

**256³** (`reion256_f0.2`, `docs/production_256.md`: 10 Mpc/h, 0.1c,
factor 16, `unsplit`, `ray_rotate = 3`; checkpoint of step 48, step 49 to
z = 5.32, 92 RT substeps, 144 M live rays, 5437 source cells). Run twice:
first with the build of commit `d927abd` on one node (job 41796867: RT
159.0 s, serial 144.0 s) and again, after Stage E1 was merged, with the
build of commit `5c0c974` on 1, 2 and 4 nodes (jobs 41812550–41812552,
2026-10-12 evening), which the table shows. The fields of the three
node counts are bitwise identical to each other (every field and
`stats.txt`), as are those of the two builds.

| 256³, step 49 | RT [s] | step [s] | RT speed-up vs 1 node | peak RSS rank 0 [GB] | min / mean / max over ranks [GB] |
|---|---|---|---|---|---|
| serial, 192 threads (1 node) | 144.4 | 173.2 | | 30.6 | |
| backend, 1 node, 24 ranks × 8 threads | 156.8 | 172.8 | 1 | 37.3 | 12.8 / 14.7 / 37.3 |
| backend, 2 nodes, 48 ranks × 8 threads | 100.2 | 116.2 | 1.56 | 37.4 | 12.6 / 13.6 / 37.4 |
| backend, 4 nodes, 96 ranks × 8 threads | 76.8 | 92.9 | 2.04 | 37.4 | 12.5 / 13.4 / 37.4 |

RT phase timers of the step, serial (`ray_tracing_1` 48.7,
`radiative_physics` 41.6, `ray_tracing_3` 44.9, `calc_ray_1` 6.0,
`post_rt` 1.1 s) against the backend on 1 / 2 / 4 nodes (`step1` 28.0 /
17.2 / 12.4, `cast1` 40.3 / 22.6 / 16.6, `phys` 22.5 / 13.2 / 9.5,
`cast3` 32.5 / 19.7 / 13.9, `sort` 13.6 / 8.6 / 5.9, `migrate` 8.9 / 7.4
/ 6.5, `hist` 3.9 / 2.2 / 1.7, `cell` 1.2 / 3.9 / 4.3, `sync` 3.9 / 3.9 /
4.7, `post_rt` 6.5 / 5.5 / 5.7 s; inside the phases `SumBoundary` 10.5 /
5.4 / 6.2 and `FillBoundary` 0.8 / 5.4 / 7.4 s). So on one node the
backend's RT is 9 % slower than the serial OpenMP code on the same 192
cores: the physics is 1.8× faster (each rank's cells are local to its
CCD), the casts are the same, and the backend pays for the sort, the
migration, the boundary sums and the gather/broadcast of `post_rt` that
the serial code does not have. Across nodes the RT scales 1.56× on 2 and
2.04× on 4 (against the spike's 1.42× and 2.13× for its RT at 256³,
`docs/framework_decision.md` §3.3): the casts, the physics and the sort
halve per doubling, while `migrate`, `post_rt` (the gather of the rays
to rank 0 and the broadcast of the fields), `sync` and the boundary
exchanges do not shrink, and `cell` (the per-rank serial part of
`calc_rt_1`) grows. 256³ is small for 4 nodes: 512 boxes over 96 ranks
is 5 boxes and 1.5 M rays per rank. The non-RT parts of the step run
redundantly on every rank with 8 threads (hydro 9.9 s against 4.9 s on
192 threads), which the step total hides only because the serial run's
`io` includes writing the fields (21.8 s against 2.1 s on the backend).
Memory: rank 0's peak is the restart (the serial ray pool read from the
checkpoint and the ray list, before the first `Redistribute`; it was
61.7 GB before the slab broadcast and the chunked scatter of commit
`d927abd`); the other ranks (12.5 GB at least, 13–15 GB on average) hold
the redundant serial state plus their share of the level (6 M rays and
0.7 M cells per rank on one node, about 1.5 GB), i.e. the redundant
state costs about 12 GB per rank at 256³, four times the §2.2 estimate
of the field arrays alone, and it does not shrink with more nodes: what
else every rank keeps (the serial RT grid, the halo and source lists,
the hydro work arrays) is to measure in Stage E, where it goes away.

**The UCX eager setting** (`docs/stage_e_hydro.md` §5: the klone
module's UCX makes AMReX's `FillBoundary` 10× slower across nodes;
`tools/hyak/env.sh` now exports `UCX_RNDV_THRESH=inf UCX_ZCOPY_THRESH=inf`
for every AMReX job). The 256³ step on 2 nodes with it (job 41813041):
RT 95.1 s against 100.2 s with the module's defaults, the fields and
`stats.txt` bitwise identical. Phases: `FillBoundary` 5.4 → 1.6 s,
`migrate` 7.4 → 5.7, `cast3` 19.7 → 18.7, `phys` 13.2 → 11.7, `step1`
17.2 → 14.2, `SumBoundary` 5.4 → 5.2 (unchanged: it is the reduction,
not the transport), while the large broadcasts get slower, `sync` 3.9 →
8.9 s (the fields to every rank) and `post_rt` 5.5 → 10.5 s (the ray
gather and the field broadcast), which eager copying through bounce
buffers does not suit. A net gain of 5 % at 256³, but **a net loss of
8 % at 512³ on 4 nodes** (job 41813111: RT 844.5 s against 782.1 s,
step 930.1 against 876.2, fields bitwise identical): `FillBoundary`
16.3 → 11.0 s, `sort` 77.3 → 62.7 and `migrate` 29.6 → 27.6 improve,
but `sync` 34.9 → 99.4 s (the 25 GB field broadcast), `post_rt` 41.6 →
106.0 s (the ray gather to rank 0 and the broadcast back) and the
deposit `SumBoundary` 50.4 → 100.3 s double. So "eager for everything"
is right for the hydro's many medium-sized messages and wrong for the
RT's large transfers; a finite threshold (eager below a few MB,
rendezvous above) is the compromise. With `UCX_RNDV_THRESH=8m
UCX_ZCOPY_THRESH=inf` (job 41815516) the 512³ step on 4 nodes is the
fastest of the three: **RT 667.1 s** (15 % under the module defaults'
782.1 s), step 746.6 s, fields bitwise identical; `cast1` 243.4 →
178.3 s with `SumBoundary` 50.4 → 56.5, `sort` 77.3 → 57.5, `step1`
136.4 → 128.2, `sync` 34.9 → 25.2 and `post_rt` 41.6 → 31.8 s (the
field broadcast's 8 MB slabs and the ray gather's chunks stay on the
rendezvous path, the many smaller messages go eager). At 256³ on 2
nodes the same setting gives RT 89.7 s (job 41815515; 100.2 with the
module's defaults, 95.1 with eager for everything): `FillBoundary`
0.5 s, `sync` 3.9, `post_rt` 5.4, `migrate` 5.7, `step1` 14.3 s, the
fields bitwise identical. So `tools/hyak/env.sh` exports
`UCX_RNDV_THRESH=8m UCX_ZCOPY_THRESH=inf` for every AMReX job
(`radray_ucx_eager`, `docs/stage_e_hydro.md` §5), and the Stage E
removal of the gather and the broadcast changes the balance again.

Fields of the last step, backend against serial (rms relative difference
over the grid, the fraction of cells beyond 1e-3, and the largest
relative difference in a single cell): fHI 1.2e-4, 3.3e-4 of the cells,
0.28; fHeI 2.2e-4, 3.8e-4, 0.22; fHeII 5e-7, 0, 2.8e-4; Γ 4.2e-4, 1.2e-4,
0.47; T 5.9e-7, 0, 3.6e-4; ρ identical. The global statistics agree to
1e-9 (x_HII,V 0.99994530 in both, T_M to 6e-9); the ray counts differ by
40 of 144 M (merges and splits by 2e-6, threshold decisions on rays whose
photon counts differ in the last bits). This is the round-off of the
deposit and merge order of §2.3 amplified in single I-front and
removal-threshold cells, the same pattern as the 64³ nightly comparison
(§4.2); nothing larger.

**512³** (`reion512_40mpc_f0.23`: 40 Mpc/h, 0.3c, factor 15, `unsplit`,
`ray_rotate = 3`, case B, no `recomb_rad`; checkpoint of step 52 at
z = 4.93, step 53 to z = 4.83, 151 RT substeps, 1.12 G live rays, 59072
source cells, 6.8 G ray-substeps). The first attempt with 24 ranks × 8
threads per node ran out of memory on 1 and 2 nodes (rank 0 at about
490 GB during the restart, the other ranks at 36 GB plus the 25 GB
broadcast buffer of `sync_fields`, over 1.4 TB per node), after which
`sync_fields` broadcasts the fields in slabs of 8 planes and the restart
scatters the rays in chunks of 32 M (commit `d927abd`). The 256³ step
re-run with these changes (job 41799553) is bitwise identical to the
first backend run (every field, `stats.txt` and the ray counts), with
rank 0's peak down from 61.7 to 37.3 GB and the others' from 14.9–18 to
12.8–14.7 GB. The 512³ runs below use 12 ranks × 16 threads per node
(jobs 41799554, 41799555, 41799556; `--mem=1400G`):

| 512³, step 53 | RT [s] | step [s] | peak RSS rank 0 [GB] | min / mean / max over ranks [GB] |
|---|---|---|---|---|
| serial, 192 threads (1 node) | 2593.4 | 2878.7 | 203.8 | |
| backend, 1 node, 12 ranks × 16 threads | OOM | | 300.5 at the kill | 140.6 mean |
| backend, 2 nodes, 24 ranks × 16 threads | 1460.7 | 1540.9 | 292.9 | 103.0 / 119.3 / 292.9 |
| backend, 4 nodes, 48 ranks × 16 threads | 782.1 | 876.2 | 292.9 | 99.7 / 109.7 / 292.9 |

The 1-node run was killed by the cgroup after 24 min, during the RT
(`srun --kill-on-bad-exit=1` ends the step; the serial run of the same
job completed): the 12 ranks summed to about 1.7 TB, 300 GB on rank 0
and 126 GB on each of the others. On 2 nodes, with half the level per
rank, the other ranks peak at 103–112 GB and the node with rank 0 holds
about 1.45 TB, which just fits. So every rank's redundant serial state is
about 100 GB at 512³ (8 × the 12.8 GB of 256³, i.e. it scales with the
grid, about 200 B per cell: the fields, the serial RT grid, the PM mesh,
the halo and source lists and the hydro work arrays), and rank 0 adds
190 GB for the restart (the serial ray pool of 1.12 G × 80 B read from
the checkpoint and the ray list it is converted to before the scatter,
both freed afterwards) and the gathered ray list of `post_rt`. The
serial code itself peaks at 204 GB (the state plus the live ray pool).
With this Stage D model (every rank a full copy of the serial state)
512³ needs 2 nodes or more at 12 ranks per node; fewer, fatter ranks
(8 × 24 threads) would fit on one node but were not run. Stage E removes
the redundant state, which is the point of it.

Timing on 2 nodes (384 cores against the serial code's 192): the RT is
1.78× faster, the step 1.87×. Phases, serial (`ray_tracing_1` 814.2,
`radiative_physics` 852.0, `ray_tracing_3` 767.2, `calc_ray_1` 125.2,
`calc_rt_1` 13.2, `post_rt` 10.5 s) against the backend (`step1` 143.4,
`cast1` 624.1 of which `SumBoundary` 368.0, `phys` 214.8, `cast3` 200.1,
`sort` 120.7, `migrate` 60.8, `hist` 33.1, `cell` 24.0, `sync` 27.8,
`post_rt` 40.1 s). The physics is 4.0× faster (the serial code's
`radiative_physics` does not scale to 192 threads on this 2-socket node
while every rank's cells are local), the casts with sort and migration
1.5× faster, and a quarter of the backend's RT is the `SumBoundary` of
the deposits over the 4096 boxes of 32³ (25 % at 512³ against 7 % at
256³): the first target for the next optimisation of the backend, with
the gather of `post_rt`. The non-RT parts of the step (hydro 42.1 s on
16 threads against 33.1 s on 192, PM 7.2 against 6.1, halos 10.8 against
9.4) run redundantly on every rank; the serial run's `io` of 226.8 s
(the same fields written by the same code, against 14.7 s on the
backend) is Lustre, not investigated.

Fields of step 53, backend (2 nodes) against serial (rms relative
difference, fraction of cells beyond 1e-3, largest single-cell relative
difference): fHI 5.8e-4, 1.1e-4, 0.61; fHeI 5.3e-4, 1.2e-4, 0.53; fHeII
1.1e-6, 4e-7, 3.6e-3; Γ 4.6e-4, 6.2e-5, 1.15; T 4.7e-7, 1.5e-8, 1.9e-3;
ρ identical. The global statistics agree to 1e-9 or better (x_HII,V
0.999995006130 in both, x_HII,M to 7e-10, T_M to 1e-9), the ray counts
differ by 56 of 1.12 G (merges by 4, splits by 36): the same round-off
of the deposit and merge order as at 256³, in single I-front cells, at
the same level, i.e. it does not grow with the problem size or the
number of ranks. The 4-node run's fields, `stats.txt` and ray counts are
bitwise identical to the 2-node run's (rank independence, as in the
tests).

On 4 nodes (768 cores, job 41799556) the RT is 1.87× faster than on 2
nodes and 3.3× faster than the serial code on one node (the step 3.3×).
Phases, 2 → 4 nodes: `step1` 143.4 → 136.4 (the emission does not
scale: it is serial work per rank), `cast1` 624.1 → 243.4, of which
`SumBoundary` 368.0 → 50.4, `phys` 214.8 → 108.1, `cast3` 200.1 →
112.0, `sort` 120.7 → 77.3, `migrate` 60.8 → 29.6, `hist` 33.1 → 17.4,
`cell` 24.0 → 13.3, `post_rt` 40.1 → 41.6 s. The physics, the sort and
the migration halve as the ranks double; the `SumBoundary` of the
deposits falls 7×, far more than the rank count, so its 368 s on 2 nodes
is not the communication volume (which hardly changes) but contention on
the node that holds rank 0 (about 1.45 TB in use of 1.5 TB, every rank
at 103–112 GB): with 4 nodes the same ranks use half the memory and the
boundary sum costs what it did at 256³. So the backend scales from 2 to
4 nodes at 512³, and its one-node cost is memory, not time. Peak RSS:
rank 0 292.9 GB again (the restart), the others 99.7 GB at least, 109.7
on average: the redundant state does not shrink with more ranks, as
expected.

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
- **Atomics on the CPU serialise the threads.** The first clean timing of
  the backend at 64³ (1 rank × 4 threads, 2 hydro steps restarted from a
  production checkpoint, 3.2 M rays) took 137 s of RT against 23.5 s for
  the serial OpenMP code: 11 s in the photon histogram (one shared
  histogram with an atomic per ray), 22 s in `cast_1` (13 int64 atomics
  per segment, rays sorted by cell so adjacent threads hit the same
  cells), and the contended cache lines slowed everything else (the gather
  phase took 40 s). With per-thread histograms and per-thread deposit
  buffers of the grown box, summed into the accumulator FAB afterwards
  (integers, so any order), and dynamic scheduling of the per-cell loops
  (source cells do 192 emissions), the same run takes 27 s, still bitwise
  identical. The spike did not see this at 32³ with 7 k rays. The GPU
  build keeps the atomics (`kThreadBuffers` is off there).
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
