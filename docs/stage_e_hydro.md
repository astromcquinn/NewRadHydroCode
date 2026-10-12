# Stage E1: the unsplit hydro step on AMReX (`hydro_backend = amrex`)

*Written 2026-10-12 on branch `claude/clever-fermi-aa4km4`. Stage E of
`docs/ROADMAP.md` §7 and of the migration plan in
`docs/framework_decision.md` §6 ("full uniform-grid code on AMReX with
MPI"), part 1: the hydro. The RT was distributed in Stage D
(`docs/stage_d.md`); the PM, the halo finder and the I/O are the other
parts of Stage E.*

## 0. Summary

| Item | State |
|---|---|
| `hydro_backend = amrex`: the unsplit MUSCL-Hancock + HLLC step with the dual-energy formulation (`hydro_scheme = unsplit`) on AMReX MultiFabs, MPI + OpenMP (`src/amrex/hydro_level.cpp`) | done |
| One implementation of the face and update arithmetic (`rrh::face_all`, `face_neighbour_max`, `update_from_faces` in `include/radray/kernels/hydro_mh.hpp`), called by the serial driver and the level | done; the serial unsplit outputs are bitwise unchanged |
| **Bitwise identical to `hydro_backend = serial`** for 1–8 ranks × 1–2 threads × box sizes 32/16/8 on three hydro problems, and in the full loop (with `rt_backend = amrex` and `serial`) | done, `ctest -L amrex_hydro` (27 tests) |
| Timing on the 4-core VM at 128³ and 256³ against the serial OpenMP code | §4 |
| The multi-node hydro slowdown of the spike (`docs/framework_decision.md` §3.3) | a hypothesis and the job that measures it (`tools/hyak/amrex_hydro_scaling.slurm`), §5; **the klone runs are the user's** |
| The legacy split scheme (`hydro_scheme = legacy`) | stays on the serial code, by design (its pencil FFT and transposes do not decompose) |

Ground rules kept: the serial code is untouched in behaviour (the kernel
code it now calls is its own arithmetic moved into the header; every
stored-reference, live-Fortran, restart and thread test passes unchanged);
AMReX stays optional; the default build and the fast suite need neither
AMReX nor MPI.

## 1. The level

`src/amrex/hydro_level.cpp` holds one uniform periodic AMReX level of the
gas grid (`BoxArray::maxSize(amrex_hydro_max_grid_size)`, default 32, the
default `DistributionMapping`):

| FabArray | Components | Ghosts | Communication |
|---|---|---|---|
| `U` (MultiFab, double) | 10: the 5 conserved variables, then the scalars of `hydro_mh.hpp` (the internal energy density `e`, the 4 species densities `w`) | 2 | `FillBoundary` at the start of every substep |
| `Es` (MultiFab) | 1: the dual-energy selected thermal energy density | 2 | `FillBoundary` after the selection |
| `En` (MultiFab) | 1: the total energy after the update, before `finalize_cell` | 2 | `FillBoundary` before the finalisation (its neighbour maximum is the reference of the η₁/η₂ rules) |
| `Q` (`FabArray<BaseFab<QCell>>`, AoS: the 20 predictor doubles of a cell) | 1 | 1 | none: computed on the valid box grown by one, from the ghost-filled `U` and `Es` |
| faces | per-box `FArrayBox` temporaries, 11 components (the HLLC flux, the 5 scalar fluxes, the interface velocity v*) on the face-centred box of each direction | 0 | none: a box computes every face of its own cells from `Q`, so two boxes sharing a face compute the same values |

One substep, in the order of the serial `hydro_step_unsplit`:

1. `U.FillBoundary`; the selection (`rrh::face_neighbour_max` of E over
   the ±`dual_nb_radius` face neighbours, `select_thermal`) on the valid
   cells; `Es.FillBoundary`;
2. the CFL speed: the maximum of `signal_sum` over the valid cells
   (`reduction(max)` inside a box, `ReduceRealMax` across ranks: a maximum
   is exact and order-independent), then the sub-cycling arithmetic of the
   serial driver, copied (`dt_cfl`, `nsub`, `dt_sub`);
3. the predictor (`predict_cell`) on the valid box grown by one;
4. the faces (`face_all`) of the three face-centred boxes, then the update
   (`update_from_faces`) in place on `U`'s valid cells, with the updated E
   also into `En`;
5. `En.FillBoundary`; `finalize_cell` on the valid cells with the
   neighbour maxima from `En`; the five cell counters as integer atomics,
   `ReduceLongSum` at the end of the step.

Threads: `ParallelForOMP` inside each box (the per-cell loops), an OpenMP
`reduction(max)` for the CFL speed; every cell is computed by the same
function from the same stencil whatever the tiling, so the thread count
does not matter.

**Host ↔ level.** In Stage E1 every rank still holds the full `GasGrid`
(the PM, the halo finder, the summaries and the I/O read it). A step copies
the rank's own boxes in from the host arrays (float → double, no
communication), sub-cycles on the level, and gathers the whole level back
onto every rank: each rank packs its boxes in global box order as floats
(the one rounding of the serial driver) and broadcasts its block, and
every rank unpacks each block as it arrives into `gas.u`, `gas.e`,
`gas.w` (so the extra memory is one block, not the level). The traffic is
40 B per cell per half step (0.67 GB at 256³), comparable to Stage D's
gather and broadcast of the RT cells, and goes away when the rest of
Stage E puts the state on the level for good. Between half steps the level holds nothing, so
the checkpoint format is unchanged and `fold_frame`, `init_unsplit_state`
and the counters of `HydroState::ucounters` are the serial ones (the
counters identical on every rank).

**Why it is bitwise identical to the serial driver.** The serial
`hydro_step_unsplit` computes each cell and face from a fixed stencil with
one function per pass (`select_thermal`, `predict_cell`, `face_all`,
`update_from_faces`, `finalize_cell`, all in `hydro_mh.hpp` now), and its
only reductions are a maximum, a minimum and integer counters. The level
calls the same functions on the same values (ghost cells are periodic
copies, as the serial padded arrays are), so every intermediate double is
the same, for any decomposition. The one trap was the finalisation: it
needs the neighbour maximum of the updated, not yet finalised E, which the
serial code reads from its separate `Un` array; the level keeps that copy
in `En` and fills its ghosts, instead of reading the in-place `U`, which
other cells' finalisation would already be changing.

## 2. Parameters and where things live

| File | What |
|---|---|
| `include/radray/hydro_backend.hpp`, `src/hydro_backend.cpp` | the switch and its parameters (AMReX-free): `hydro_backend = serial \| amrex`, `amrex_hydro_max_grid_size` (32); `check_hydro_backend_param` refuses `hydro_scheme = legacy` and a build without AMReX; the process functions are those of `rt_backend.hpp` (one `backend_init` serves both backends) |
| `src/amrex/hydro_level.cpp` | the level (target `radray_amrex`, C++20) |
| `src/amrex/hydro_level_stub.cpp` | the no-AMReX build: the level functions throw |
| `src/hydro_unsplit.cpp` | `hydro_step_unsplit` dispatches to `amrex_hydro_step` after `fold_frame` |
| `include/radray/kernels/hydro_mh.hpp` | `FaceData`, `face_all`, `face_neighbour_max`, `update_from_faces` (new), used by both drivers |
| `apps/hydro_driver.cpp` | `backend_init` at start, rank 0 writes and prints, `--set key=value` (new; `nsteps`, the hydro scheme and backend keys) |
| `apps/radhydro.cpp` | parses, checks and prints `hydro_backend`; the level's timers in the per-step output and `timing.txt` |

Both backends combine freely: `rt_backend` and `hydro_backend` are
independent switches (`radhydro --set rt_backend=amrex --set
hydro_scheme=unsplit --set hydro_backend=amrex` is the Stage E1
configuration; with MPI every rank runs the whole driver as in Stage D).

```
mpirun -np 4 -x OMP_NUM_THREADS=2 build-amrex/hydro_driver PROB out --set hydro_scheme=unsplit \
  --set hydro_backend=amrex [--set amrex_hydro_max_grid_size=32]
mpirun -np 4 build-amrex/radhydro RUN out --set legacy_rt=0 --set rt_backend=amrex \
  --set hydro_scheme=unsplit --set hydro_backend=amrex --set deterministic_parallel=0
```

## 3. Validation: `ctest -L amrex_hydro`

Registered with `RADRAY_AMREX=ON` (inside the `amrex` label; 27 tests):

- `amrex_hydro_bitwise_<problem>_r<ranks>t<threads>_mgs<box>`:
  `hydro_driver --set hydro_scheme=unsplit --set hydro_backend=amrex`
  against the serial unsplit run of the same problem with `compare.py
  --rtol 0 --atol-frac 0` (every field of every step, every stat), for
  `sod_x`, `random3d` and `sound_x` (N = 32) and (ranks, threads, box) ∈
  {(1,1,32), (1,2,16), (2,1,16), (2,2,8), (4,1,16), (4,2,8), (8,1,8)};
- `amrex_hydro_full_ranks_smoke`: `radhydro` on the full-loop smoke
  problem with `rt_backend = amrex` and `hydro_backend = amrex`, 1, 2 and 3
  ranks, every file (fields, movie frames, checkpoints) byte-identical
  (`rank_check.py`);
- `amrex_hydro_full_vs_serial_hydro`: the same run with `hydro_backend =
  serial` and `amrex`, byte-identical (`backend_pair_check.py`, new);
- `amrex_hydro_restart_smoke`: a bit-identical restart through both
  backends.

Checked by hand in this session before the tests were written: the
serial-RT full loop with `hydro_backend = amrex` is byte-identical to the
all-serial run (95 files), the AMReX-RT loop with the AMReX hydro at 1, 2
and 3 ranks is byte-identical to the AMReX-RT loop with the serial hydro,
and `random3d` on 2 ranks is identical with 1 and 4 threads. The default
build's fast suite is unchanged (the stored hydro references are the
legacy scheme; the unsplit scheme's outputs on four hydro problems and
the smoke run were recorded before the kernel refactor and are bitwise
the same after it).

## 4. Timing (4-core VM, OpenMPI 4.1.6, gcc 13)

`random3d` at 128³ and 256³ (`tools/hyak/make_hydro_problem.py`, the
random3d field of the hydro tests at that size: Mach 20–40 compressive
flow, log-normal density, fixed da = 0.002, 2 hydro steps = 4 unsplit
calls), the hydro time of `hydro_driver` ("hydro" of the final line;
the serial code with OpenMP, the backend with `mpirun --bind-to none`).

| 128³, 2 steps (4 unsplit calls, 8 substeps) | hydro [s] | peak RSS, all processes [MB] | `fill` [s] | `copyout` [s] (of which the broadcast) |
|---|---|---|---|---|
| serial, 4 threads | 12.2 | 819 | | |
| serial, 1 thread | 33.8 | 803 | | |
| backend, 1 rank × 4 threads, boxes 32³ | 14.7 | 1494 | 0.36 | 0.56 (0) |
| backend, 2 × 2, 32³ | 10.8 | 1701 | 1.03 | 0.33 (0.13) |
| backend, 4 × 1, 32³ | 10.1 | 1965 | 1.33 | 0.46 (0.25) |
| backend, 4 × 1, 16³ | 11.5 | 2418 | 1.34 | 0.76 (0.53) |
| backend, 4 × 1, 64³ | 10.3 | 2126 | 0.99 | 0.74 (0.50) |
| backend, 1 × 1, 32³ | 32.2 | 1438 | 1.13 | 0.37 (0) |

| 256³, 2 steps | hydro [s] | peak RSS, all processes [MB] | `fill` [s] | `copyout` [s] (broadcast) |
|---|---|---|---|---|
| serial, 4 threads | 114.5 | 6349 | | |
| backend, 2 × 2, 64³ | 99.4 | 12128 | 4.9 | 7.2 (4.1) |
| backend, 4 × 1, 32³ and 64³ | killed (out of memory at 13.6 GB on the 15 GB VM) | | | |

(The 16³ and 32³ rows of the 128³ table with 4 × 1 were taken with the
first, `MPI_Allgatherv` version of the gather; the 2 × 2 and 4 × 1, 32³
rows with the per-rank broadcast that replaced it, which is why their
`copyout` is smaller.)

Reading:
- **Speed.** With MPI ranks the backend beats the serial OpenMP code on the
  same 4 cores: 4 × 1 is 17 % faster at 128³, 2 × 2 is 13 % faster at
  256³; on one core the two are equal (32.2 against 33.8 s). One rank with
  4 threads is 20 % slower than the serial code: `ParallelForOMP` threads
  each box's loops separately (fork/join per box and per pass), where the
  serial driver threads whole k-planes. Ranks are the better use of cores,
  as the spike found (its best layout was one rank per CCD).
- **Box size.** At 128³ on one machine the box size matters little
  (16³ costs 10 % more than 32³ and 64³ from the ghost volume, `fill`
  1.33 against 0.99 s). The multi-node question is §5.
- **Communication.** `fill` (three `FillBoundary` per substep) is 5 % of
  the time and the host gather 7 % at 256³ on 2 ranks; both are shared
  memory here. The gather (`copyout`) is the Stage E1 price and goes away
  with the state on the level.
- **Memory.** The Stage E1 model costs the full host state per rank (the
  driver's u, w, v, e: 0.8 GB per rank at 256³) plus the level (about
  330 B per cell with the ghosts of 32³ boxes: U 114, Es and En 23, Q 192)
  and one block of the gather per rank. That is why 4 ranks at 256³ do
  not fit in this VM's 15 GB while 2 ranks (12.1 GB) do, and why the
  serial driver needs half as much (6.3 GB, its 41 doubles per cell
  without ghosts). On a klone node (1.5 TB) the same model at 512³ with
  24 ranks is 24 × 6.4 GB of host state plus 45 GB of level, which fits,
  but the per-rank host state is the first thing Stage E2–E4 should
  remove. The level's own largest item, `Q` (20 doubles per cell), could
  be halved by computing the predictor per face pair instead of storing
  it, at the cost of recomputing it 6× per cell.

## 5. The multi-node hydro slowdown: hypothesis and the measurement

`docs/framework_decision.md` §3.3: the spike's AMReX hydro cost
6.8–7.6 µs per cell-update·core on 2 and 4 nodes against 0.69 on one
node, with 16³ boxes, while the RT on the same `FillBoundary` and
Parthenon's hydro did not slow down. It could not be reproduced here (one
machine). What the level is built to test:

- **Ghost volume.** A 16³ box with 2 ghost layers exchanges 3904 ghost
  cells per 4096 valid cells, i.e. about one ghost cell per cell per
  `FillBoundary`; at 5 doubles that is 40 B per cell-update across the
  network for the fraction of faces that cross nodes, per `FillBoundary`.
  At 0.69 µs per cell-update·core a node updates 2.8×10⁸ cells/s, so the
  ghost traffic of the spike's hydro is of order 10 GB/s, at the limit of
  one node's injection bandwidth, whereas inside a node it is a memory
  copy. The RT does not show it because its cost per cell is 100× larger.
  The spike's prototype also did one `FillBoundary` per substep of 5
  components, with no overlap. The level here needs three exchanges per
  substep (`U`: 10 components, `Es` and `En`: 1 each), so the hypothesis
  predicts the slowdown appears here too with 16³ boxes and shrinks with
  32³ (0.42 ghost cells per cell) and 64³ (0.2).
- **The measurement** is `tools/hyak/amrex_hydro_scaling.slurm`: 24 ranks
  × 8 threads per node on 1, 2 and 4 nodes with box sizes 16, 32 and 64,
  on a 256³ random3d problem, against the serial OpenMP code on one node;
  it prints the level's timers (`fill` separately from the compute
  phases) and checks that every layout gives byte-identical outputs. If
  `fill` dominates on 2 nodes with 16³ boxes and not with 64³, the
  hypothesis holds and the production layout uses large boxes (one or two
  boxes per rank) and, later, overlapped communication. If the slowdown
  is independent of the box size, the next suspect is the MPI transport of
  the klone OpenMPI module across nodes (UCX/verbs settings), which the
  RT would hide.

## 6. What this stage does not do (yet)

- The state on the level between steps: the host gather per half step
  (40 B per cell) and the redundant host state per rank (0.8 GB per rank
  at 256³, the reason 4 ranks at 256³ do not fit in a 15 GB VM, §4) stay
  until the PM (E2), the halo finder and sources (E3) and the I/O (E4)
  are distributed.
- The legacy split scheme on the backend (never: its moving-frame pencil
  sweeps with the FFT smoothing do not decompose; `check_hydro_backend_param`
  refuses it).
- GPU: the level's loops are `ParallelForOMP` lambdas with host
  accessors; `hydro_level.cpp` is not a CUDA target yet (the AoS `QCell`
  FAB and the per-box face temporaries are the parts to revisit, as for
  the RT level's serial folds).
- Load balancing: the hydro level uses the default `DistributionMapping`;
  the ray-weighted mapping of Stage E (E5) applies to the RT level, and
  the two levels need not share a decomposition (`docs/design_parallel_and_encoding.md`
  §3.1).

## 7. Findings and traps

- A `FabArray` with a static lifetime must be released before
  `amrex::Finalize`: the hydro level lives in a static in
  `hydro_level.cpp` (as the serial driver's work arrays do), and
  `backend_finalize` releases it first (`amrex_hydro_release`); otherwise
  the destructor runs after the arenas are gone ("pure virtual method
  called" at exit, the first run of every test).
- `hydro_driver` had no `--set`; it has one now, so the tests can switch
  the scheme and the backend without editing `params.txt`.
- The flux temporaries are per box, not per level: a face-centred
  MultiFab for each direction would cost 33 doubles per cell for data that
  no other box ever reads. With per-box temporaries the level costs about
  32 doubles per cell plus ghosts (U 10, Es, En, Q 20), close to the
  serial driver's 41.
