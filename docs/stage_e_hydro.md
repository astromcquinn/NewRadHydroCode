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
| Timing on the 4-core VM at 128³ and 256³ against the serial OpenMP code; klone, 256³ on 1, 2 and 4 `cpu-g2` nodes | §4, §4.1 |
| The multi-node hydro slowdown of the spike (`docs/framework_decision.md` §3.3) | measured on klone (2026-10-12): it is the UCX rendezvous (zero-copy) path of the klone OpenMPI, not the ghost volume; `UCX_RNDV_THRESH=8m UCX_ZCOPY_THRESH=inf` (eager below 8 MB) removes it, §5 |
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

### 4.1 klone: 256³ on 1, 2 and 4 nodes (2026-10-12)

`tools/hyak/amrex_hydro_scaling.slurm` on `cpu-g2` (2 × 96-core EPYC 9654
per node, InfiniBand, OpenMPI 4.1.6 with UCX 1.16 from the klone `ompi`
module): the 256³ random3d problem of `tools/hyak/make_hydro_problem.py`,
2 steps (4 unsplit calls, 22 substeps, 3.7×10⁸ cell-updates), 24 ranks ×
8 threads per node, `--cpu-bind=cores`, against the serial OpenMP code on
192 threads of one node (jobs 41812553–41812555 with the module's
default MPI settings, 41813037–41813039 with `UCX_RNDV_THRESH=inf
UCX_ZCOPY_THRESH=inf`, "UCX eager" below; the production setting of §5
is the 8 MB threshold, whose hydro `fill` lies between the two).
Every layout, box size and transport setting below is byte-identical to
the serial run (`diff -r` of every step's fields and stats).

| 256³, 2 steps, hydro [s] (`fill` [s]) | boxes 16³ | 32³ | 64³ |
|---|---|---|---|
| serial, 192 threads | 38.0 | | |
| backend, 1 node, 24 × 8, module defaults | 8.20 (0.71) | 8.08 (0.44) | 8.46 (0.42) |
| backend, 2 nodes, 48 × 8, module defaults | 8.63 (4.47) | 7.85 (3.92) | 8.27 (3.80) |
| backend, 4 nodes, 96 × 8, module defaults | 9.15 (4.15) | 12.34 (6.84) | 8.65 (3.21) |
| backend, 1 node, 24 × 8, UCX eager | 7.47 (0.75) | 6.99 (0.41) | 7.31 (0.47) |
| backend, 2 nodes, 48 × 8, UCX eager | 5.51 (0.64) | 5.32 (0.47) | 5.93 (0.39) |
| backend, 4 nodes, 96 × 8, UCX eager | 5.22 (0.43) | 5.39 (0.29) | 6.34 (0.32) |

The other timers at 32³ boxes, 1 / 2 / 4 nodes (UCX eager): `predict`
0.45 / 0.21 / 0.11, `flux` 0.80 / 0.42 / 0.21, `update` + `finalize`
0.7 / 0.3 / 0.1, `copyout` 3.2 / 2.4 / 3.0 (of which the per-rank
broadcasts 3.0 / 2.2 / 2.8), and the rest of the driver's step outside
the level (`fold_frame`, the serial per-rank work of
`radhydro_hydro_step`) about 2 s on every layout. Reading:

- On one node the backend is 4.7× faster than the serial OpenMP code
  (8.1 against 38.0 s; the serial unsplit driver does not scale to 192
  threads across the two sockets, while every rank's cells are local to
  its CCD).
- The compute phases scale: they halve per doubling of the nodes
  (`predict` + `flux` + `update` + `finalize` 1.95 → 0.93 → 0.43 s).
  What does not is the Stage E1 model: the gather of the whole level to
  every rank (`copyout`, 0.67 GB broadcast per half step to every rank,
  2.4–3.0 s whatever the node count) and the redundant host work (2 s),
  which together are 5 of the 5.3 s on 2 nodes. They go away with the
  state on the level (E2–E4), after which the hydro step at 256³ is
  about 1 s on 2 nodes.
- The box size matters little: 32³ is the best or within 3 % of it on
  every layout; 16³ costs 5–10 % more (ghost volume in `fill` and
  per-box overhead), 64³ gives 64 boxes, fewer than the 96 ranks of 4
  nodes (a third idle) and 1–4 % slower on 2 nodes. **Production keeps
  32³** (`amrex_hydro_max_grid_size = 32`, the default, the same as the
  RT level's `amrex_max_grid_size`).

## 5. The multi-node hydro slowdown: measured, and it is the UCX transport

`docs/framework_decision.md` §3.3: the spike's AMReX hydro cost
6.8–7.6 µs per cell-update·core on 2 and 4 nodes against 0.69 on one
node, with 16³ boxes, while the RT on the same `FillBoundary` and
Parthenon's hydro did not slow down. The hypothesis written before the
measurement was the ghost volume of 16³ boxes against the node's
injection bandwidth (about one ghost cell per cell per `FillBoundary`
at 16³, 0.42 at 32³, 0.2 at 64³, so the slowdown should shrink with the
box size). The measurement (§4.1, jobs of 2026-10-12) says otherwise:

1. **The slowdown appears.** With the klone module's default MPI
   settings the step takes 7.9–8.6 s on 2 nodes and 8.7–12.3 s on 4,
   against 8.1 s on one: no speed-up at all, while the compute phases
   halve per doubling. All of it is `fill` (the three `FillBoundary` per
   substep, 66 exchanges): 0.44 s on one node (shared memory), 3.9 s on
   2 nodes, 6.8 s on 4, i.e. 60–100 ms per exchange.
2. **It does not scale with the ghost volume.** On 2 nodes `fill` is
   4.47, 3.92 and 3.80 s for 16³, 32³ and 64³ boxes (ghost volume 5×
   apart); on 4 nodes 4.15, 6.84 and 3.21. The volume across the network
   is small: at 32³ about 12 MB per rank per exchange, 0.8 GB per rank
   per run, under a second at any sensible bandwidth. So the hypothesis
   is wrong and the per-exchange cost is per message, not per byte.
3. **It is the MPI transport.** The diagnostics (job 41812852 and the
   2-node runs at 32³ below, every one byte-identical to the serial run):
   the nodes have InfiniBand (mlx5_0, port active), the module's OpenMPI
   selects the UCX PML (priority 51 over ob1's 20, devices `mlx*`, UCX
   picks `dc_mlx5`; the module sets `UCX_TLS=all`), and forcing it
   (`OMPI_MCA_pml=ucx`) gives the same `fill` (4.9 s), so the default
   runs were already on UCX over InfiniBand, not on TCP.

   | 2 nodes, 48 × 8, 32³ boxes | `fill` [s] | hydro [s] |
   |---|---|---|
   | module defaults (UCX, `dc_mlx5`) | 3.92 | 7.85 |
   | `OMPI_MCA_pml=ucx` (forced, verbose: `dc_mlx5/mlx5_0:1`) | 4.89 | 10.24 |
   | `UCX_TLS=rc_x,sm,self` | 3.85 | 10.38 |
   | `UCX_TLS=tcp,sm,self` (UCX over TCP) | 4.96 | 9.03 |
   | `UCX_RNDV_SCHEME=put_zcopy` | 4.13 | 10.02 |
   | `UCX_PROTO_ENABLE=n` (the pre-1.15 protocol selection) | 2.54 | 6.31 |
   | `MALLOC_MMAP_THRESHOLD_=32M MALLOC_TRIM_THRESHOLD_=2G MALLOC_TOP_PAD_=512M` | 3.86 | 7.69 |
   | `OMPI_MCA_pml_ucx_opal_mem_hooks=1` | 4.03 | 7.87 |
   | `OMP_WAIT_POLICY=passive` | 3.48 | 7.14 |
   | 24 ranks × 1 thread per node | 2.81 | 17.0 |
   | 8 ranks × 24 threads per node | 1.59 | 5.29 |
   | 96 ranks × 2 threads per node | 10.7 | 19.1 |
   | `UCX_RNDV_THRESH=inf` (no rendezvous) | 1.25 | 6.11 |
   | `UCX_TLS=rc_x,sm,self UCX_RNDV_THRESH=inf` | 1.40 | 6.43 |
   | **`UCX_RNDV_THRESH=inf UCX_ZCOPY_THRESH=inf`** (eager, copied) | **0.76** | **5.63** |
   | `OMPI_MCA_pml=ob1 OMPI_MCA_btl=tcp,self,vader` (no UCX, TCP over IPoIB) | 0.37 | 5.36 (bus error at exit) |

   The cost is UCX's rendezvous zero-copy path for the medium-sized
   messages of `FillBoundary` (a few hundred KB to a few MB per
   neighbour pair, posted all at once and waited on together): it grows
   with the number of ranks per node (1.6 → 3.9 → 10.7 s for 8, 24 and 96
   ranks per node) and goes away when UCX is told to send everything
   eagerly through its bounce buffers (`UCX_RNDV_THRESH=inf
   UCX_ZCOPY_THRESH=inf`: `fill` 0.76 s, and 0.47 / 0.29 s in the full
   runs of §4.1 on 2 / 4 nodes, the same as one node's shared memory),
   or when UCX is bypassed altogether (ob1 over TCP: 0.37 s, but the run
   died with a bus error in the shared-memory BTL at `MPI_Finalize`, so
   not a usable configuration without more work). It is not the memory
   registration cache (the glibc thresholds that keep the buffers mapped
   change nothing, nor OpenMPI's memory hooks), not the InfiniBand
   transport (UCX over TCP is as slow), not the OpenMP threads (1 thread
   per rank is still 2.8 s, `passive` waiting 3.5 s). AMReX's
   communication buffers on a CPU build come from the plain-malloc
   `BArena` (`AMReX_Arena.cpp`, `the_comms_arena = The_BArena()` without
   GPU), allocated and freed per exchange, which is what the zero-copy
   path dislikes; a caching comms arena is not selectable at run time on
   CPU. The RT did not show the slowdown in the spike because its
   `FillBoundary` traffic is small against its compute and its ray
   migration is a few large messages.
4. **What production uses.** Eager for everything is wrong for the RT:
   its large transfers (the 25 GB field broadcast and the ray gather of
   `post_rt` at 512³, the deposit `SumBoundary`) slow down by 2× on the
   bounce-buffer path and the 512³ RT on 4 nodes is 8 % slower overall
   (`docs/stage_d.md` §4.4). A finite threshold keeps those on the
   rendezvous path and sends the many medium-sized messages eagerly:

   | `UCX_RNDV_THRESH` (`UCX_ZCOPY_THRESH=inf`) | hydro 2 nodes, 32³: `fill` [s], two runs | RT 256³, 2 nodes [s] | RT 512³, 4 nodes [s] |
   |---|---|---|---|
   | module defaults | 3.92, 4.89 | 100.2 | 782.1 |
   | `8m` | 1.47, 0.60 | **89.7** | **667.1** |
   | `16m` | 2.14 | | |
   | `32m` | 2.05, 1.82 | | 667.8 |
   | `inf` | 0.47, 0.76, 0.51 | 95.1 | 844.5 |

   The single jobs above ran on different node pairs and scatter by up
   to 2×; `tools/hyak/amrex_hydro_ucx_sweep.slurm` (job 41820425) runs
   the settings interleaved three times on one fixed pair, with every
   run byte-identical:

   | fixed node pair, 48 × 8, 32³ boxes: `fill` [s], three runs | `copyout` [s] | hydro [s] |
   |---|---|---|
   | module defaults: 4.37, 3.94, 4.00 | 1.2 | 7.8–8.3 |
   | `8m`: **0.38, 0.34, 0.37** | 1.1–1.2 | **4.0–4.1** |
   | `16m`: 1.32, 0.47, 0.57 | 2.4 | 5.3–6.2 |
   | `32m`: 1.31, 1.49, 1.33 | 2.4 | 6.1–6.3 |
   | `inf`: 1.29, 0.53, 0.47 | 2.4 | 5.4–6.1 |

   OpenMPI's `pml_monitoring` histogram of the same run explains the
   order: the exchange's messages are 64 KB to 4 MB (66 per neighbour
   pair per run, i.e. one per `FillBoundary`), while the Stage E1
   gather's per-rank broadcast blocks are 16–32 MB. With the 8 MB
   threshold the former go eager and the latter stay on the rendezvous
   path; from 16 MB up the broadcast blocks go eager too, which doubles
   `copyout` (1.2 → 2.4 s) and, by occupying UCX's eager queues, slows
   the exchange that follows. For the RT at 512³ `8m` and `32m` are the
   same (667.1 and 667.8 s), so 8 MB is the setting.
   `tools/hyak/env.sh` exports **`UCX_RNDV_THRESH=8m UCX_ZCOPY_THRESH=inf`**
   (`radray_ucx_eager`, called by the AMReX job scripts;
   `RADRAY_UCX_RNDV_THRESH` changes the threshold, `RADRAY_UCX_EAGER=0`
   keeps the module's defaults) and the box size stays 32³ (§4.1). With
   it the hydro step on 2 nodes is 6.8–7.8 s (the gather and the
   redundant host work of Stage E1 are the rest, and E2–E4 remove them),
   the 256³ RT on 2 nodes is 10 % faster than with the module's defaults
   and the 512³ RT on 4 nodes 15 % faster, all outputs bitwise identical
   to the module-default runs.

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
