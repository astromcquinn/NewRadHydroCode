# Parthenon prototype of the framework spike

Two programs built on [Parthenon](https://github.com/parthenon-hpc-lab/parthenon)
v25.12 (bundled Kokkos 4.6, OpenMP backend). They reproduce the serial
references of `spike/serial` **bit for bit**:

| program | reproduces | what Parthenon provides |
|---|---|---|
| `hydro_parthenon` | `hydro_serial` (MUSCL-Hancock + HLLC, periodic n^3) | meshblocks, `cons` variable with `Metadata::FillGhost`, ghost exchange tasks, `par_for`/`par_reduce` |
| `rt_parthenon` | `rt_serial` (sorted-ray-array RT substep, `legacy_rt = 0`) | meshblocks, RT cell variables, FillGhost exchange of `tau`/`fgam`, a `Swarm` for the rays with swarm boundary communication for migration |

Both use the shared kernels of `spike/common` unchanged and keep the kernel
bodies in `KOKKOS_LAMBDA` form.

Files: `CMakeLists.txt`, `parth_common.hpp` (ghost-exchange task
collection, block geometry, rank-0 gather, argv filter),
`hydro_parthenon.cpp`, `rt_parthenon.cpp`, input decks `hydro.in`, `rt.in`.
A third executable, `rt_parthenon_aos`, is `rt_parthenon.cpp` built with
`-DRT_CELLS_AOS`: the 45-float `CellData` lives in a per-block float AoS
`Kokkos::View` and only tau/fgam are Parthenon variables. It exists only to
measure what the Parthenon variable layout costs, and it is bitwise
identical as well.

## Build

```
cmake -S spike/parthenon -B build-spike-parthenon \
      -DPARTHENON_SOURCE_DIR=/home/user/ext/parthenon     # a v25.12 checkout (git clone --recursive)
cmake --build build-spike-parthenon -j2
```

* Parthenon is consumed with `add_subdirectory` (as AthenaPK does). Without
  `PARTHENON_SOURCE_DIR` the project uses `FetchContent` at
  `PARTHENON_TAG` (default `v25.12`); `-DFETCHCONTENT_SOURCE_DIR_PARTHENON=<dir>`
  points that at a local checkout instead.
* The CMake file forces `Kokkos_ENABLE_OPENMP=ON`, `PARTHENON_DISABLE_HDF5=ON`,
  examples/tests off, `Real = double`, and defaults
  `PAR_LOOP_LAYOUT=MANUAL1D_LOOP`. Parthenon's own host default
  `SIMDFOR_LOOP` runs every `par_for` serially (see the notes). Our targets get
  `-ffp-contract=off -fno-fast-math` and `-O2` (as `spike/serial`);
  `../cmake/radray_core.cmake` provides `radray_core`/`spike_common`.
* Needs MPI (OpenMPI here); Parthenon enables it automatically when found.

Build times on this 4-core VM (shared with another agent's builds, so
upper bounds):

| what | wall | CPU |
|---|---|---|
| Parthenon + Kokkos library only, inside this project, `-j2` | 5 min 23 s | 9 min 27 s user |
| everything from scratch (Parthenon + Kokkos + our 3 executables), `-j2`, machine shared | 10 min 39 s | |
| Parthenon + Kokkos + all examples, standalone, `-j4` (`/home/user/ext/build-parthenon`) | 7 min 27 s | 24 min 19 s user |
| `rt_parthenon` (compile + link) | 40 s | |
| `hydro_parthenon` (compile + link) | 31 s | |
| compile-only Kokkos **CUDA** build of both TUs (nvcc 12.0, sm_80, configure 24 s) | 144–304 s + 115 s (machine load) | |

## Run

Both programs need a Parthenon input deck (`-i`); the mesh-related values in
it are overwritten by the program, so the decks here only make Parthenon
happy. Any `block/par=value` argument is passed on to Parthenon as an
override; everything else is our own options.

```
# hydro: same options as hydro_serial, plus --mb (meshblock size)
mpirun -np 4 build-spike-parthenon/hydro_parthenon -i spike/parthenon/hydro.in \
       --n 64 --steps 20 --ic blast --mb 16 --out U.bin

# RT: same arguments as rt_serial, plus --mb (meshblock size, >= 8)
mpirun -np 4 build-spike-parthenon/rt_parthenon -i spike/parthenon/rt.in \
       tests/problems/lognormal_16src_smallbuf out_dir --set legacy_rt=0 --mb 16
```

`OMP_NUM_THREADS` sets the Kokkos OpenMP thread count. With several ranks
per node use `mpirun --bind-to none` (OpenMPI binds each rank to one core by
default, so 2 threads per rank would share a core), and when
oversubscribing set `OMP_WAIT_POLICY=passive`.

rt_parthenon prints the per-step line of rt_serial and at the end

```
 ranks R  threads T  blocks B (mb M)  wall W s  ray-substeps S  ns/ray-substep*threads*ranks X
 timers (max over ranks): cell .. halo .. hist .. step1 .. cast1 .. fold .. phys .. cast3 .. sort .. migrate .. post ..
```

The timers are those of rt_serial plus `halo` (tau/fgam ghost exchange and
the float mirror), `fold` (ghost -> owner deposit sums), `migrate` (swarm
communication); `sort` is the swarm -> cell-sorted array gather plus the
write-back of the surviving rays into the swarm.

## Design

### hydro_parthenon

* Variable `cons` (5 components, `Metadata::Cell | Independent | FillGhost`),
  nghost = 2, periodic mesh `[0,1]^3`.
* Per step: ghost exchange as a `TaskCollection` (one task list per default
  block partition: `StartReceiveBoundBufs` then
  `AddBoundaryExchangeTasks`), then per block: `par_reduce(Max)` of
  `max_speed` (+ `MPI_Allreduce(MAX)`, exact), `par_for` over
  interior + 1 layer for `rrh::predict_cell` (the `U(i,j,k,m)` accessor reads
  the Parthenon variable directly), `par_for` over faces for
  `rrh::face_flux`, `par_for` for `rrh::update_cell`. Q and F live in
  per-block scratch `Kokkos::View`s (component fastest, as the kernels want
  contiguous 5/20-vectors). Our own dt and step count; no Parthenon
  driver/integrator.
* Output: rank 0 gathers (`MPI_Gatherv` of (cell index, 5 doubles)) and
  writes raw float64 (5,n,n,n), identical to hydro_serial. The mass/energy
  sums use hydro_serial's order.

### rt_parthenon

Same phase order as rt_serial, per substep:

| phase | implementation |
|---|---|
| gather (`sort`) | swarm -> cell-sorted `RayP` array per block: counting sort by home cell (atomic counts, `parallel_scan`, atomic scatter), then all swarm particles are removed |
| `cell` | `calc_rt_1_cell` per interior cell |
| `halo` | FillGhost exchange of `tau` (task collection), float mirror `tf(k,j,i,0..4)` incl. ghosts |
| `hist` | per-block atomic int64 histogram + npos -> `MPI_Allreduce(SUM)` -> `RTRun::thresholds` |
| `step1` | `step1_count` per cell, `parallel_scan`, `step1_write` per cell (scratch idx/key per ray, as rt_serial), thick-source deposits |
| `cast1` | `cast_1` per output ray; deposits `Kokkos::atomic_add` into per-block int64 accumulators `acc(k,j,i,13)` that include a 4-cell ghost zone |
| `fold` | hand-written reverse exchange: every ghost slab is added into the owning neighbour's interior (same rank: direct add; else `MPI_Isend/Irecv`, tag = receiver gid*27 + direction) |
| `phys` | `physics_cell` per interior cell from the folded accumulators |
| `halo` | FillGhost exchange of `fgam`, mirror `tf(k,j,i,5..9)` |
| `cast3` | `cast_3` per output ray, reading tau/fgam (incl. ghosts) |
| `sort` | surviving rays -> `AddEmptyParticles` + write; swarm position = centre of the image of the home cell nearest the block (may lie outside `[0,N)`; Parthenon's periodic swarm BC wraps it on receipt) |
| `migrate` | swarm boundary communication (`ResetCommunication` in a first task region, then `Send` -> `Receive` per block); one hop is enough |
| `post` (per hydro step) | `redshift_ray` on the swarm, `post_cell`, fixed-point ngam total and ray count via `MPI_Allreduce` |

* Cells: `rtc` (35 components: `CellData` without `tau`, `fgam`), `tau` (5),
  `fgam` (5), all `Real` = double holding the float values (float -> double
  -> float is exact). `tau` and `fgam` have `FillGhost`, nghost = 4. Each
  kernel loads/stores only the `CellData` fields it uses (a compile-time bit
  mask per phase; access through one base pointer + component stride per
  variable), because every component of a Parthenon variable is its own
  memory stream.
* Rays: `Swarm "rays"` with `int ri[2]` (nside, ipix) and `Real rf[18]` (the
  18 floats of `rrk::RayP`); positions x,y,z are only used for routing; the
  authoritative position is the float `x` in `rf`. `NoPersistentParticleIds`.
* cast_1/cast_3 get wrapped global cell indices: `local_img` maps them to the
  nearest periodic image in the block's ghost-grown range and
  `Kokkos::abort`s if none is in range (never triggered in the tests).
* Determinism: within a cell the ray order is arbitrary (the atomic scatter
  is thread-order dependent), which is harmless because merges sort
  canonically and deposits/statistics are integer sums. Hence results do
  not depend on rank count, block size or thread count.
* Output: rank 0 gathers `CellData` in global order and calls
  `RTRun::write_step` (same files as rt_serial).

## Validation

References: `spike/serial` built with `cmake -S spike/serial -B build-serial`,
run with `OMP_NUM_THREADS=1`:

```
build-serial/rt_serial tests/problems/uniform_8src            ref/uniform_8src            --set legacy_rt=0
build-serial/rt_serial tests/problems/lognormal_16src_smallbuf ref/lognormal_16src_smallbuf --set legacy_rt=0
build-serial/rt_serial tests/problems/lognormal_16src_smallbuf ref/lognormal_lp0 --set legacy_rt=0 --set legacy_physics=0
build-serial/hydro_serial --n 32 --steps 10 --ic blast --out ref_32_blast.bin      # etc.
```

Each Parthenon run:

```
OMP_NUM_THREADS=$t OMP_WAIT_POLICY=passive mpirun --bind-to none --oversubscribe -np $np \
   rt_parthenon -i spike/parthenon/rt.in tests/problems/$p out --mb $mb --set legacy_rt=0 [--set legacy_physics=0]
python3 tools/compare/compare.py ref/$ref out --rtol 0 --atol-frac 0      # all fields + stats, 10 steps

OMP_NUM_THREADS=$t mpirun -np $np hydro_parthenon -i spike/parthenon/hydro.in --n $n --steps $s --ic $ic --mb $mb --out U.bin
cmp ref.bin U.bin
```

Results (final code, `PAR_LOOP_LAYOUT=MANUAL1D_LOOP` build):

| program | problem / flags | ranks | meshblock | OMP threads | result |
|---|---|---|---|---|---|
| rt_parthenon | uniform_8src, `legacy_rt=0` | 1, 2, 4, 8 | 16^3, 8^3 | 1, 2 | **16/16 PASS** |
| rt_parthenon | uniform_8src, `legacy_rt=0` | 1 | 32^3 (one block, its own periodic neighbour) | 1, 2 | **2/2 PASS** |
| rt_parthenon | lognormal_16src_smallbuf, `legacy_rt=0` | 1, 2, 4, 8 | 16^3, 8^3 | 1, 2 | **16/16 PASS** |
| rt_parthenon | lognormal_16src_smallbuf, `legacy_rt=0 legacy_physics=0` | 1, 4, 8 | 16^3, 8^3 | 1, 2 | **12/12 PASS** |
| hydro_parthenon | blast, n=32, 10 steps | 1, 2, 4 | 16^3, 8^3 | 1, 2 | **12/12 PASS** (`cmp`) |
| hydro_parthenon | sod, n=32, 10 steps | 1, 2, 4 | 16^3, 8^3 | 1, 2 | **12/12 PASS** |
| hydro_parthenon | blast, n=64, 20 steps | 1, 2, 4 | 16^3, 8^3 | 1, 2 | **12/12 PASS** |
| hydro_parthenon | blast, n=32, one 32^3 block | 1 | 32^3 | 1, 2 | **2/2 PASS** |
| rt_parthenon_aos | uniform_8src and lognormal_16src_smallbuf, `legacy_rt=0` | 1, 4, 8 | 16^3, 8^3 | 1, 2 | **24/24 PASS** |

The same matrix also passed (46/46 RT, 38/38 hydro) with Parthenon's default
`SIMDFOR_LOOP` build before the loop-layout switch, so the results do not
depend on the loop pattern either.

"PASS" for RT is `compare.py ... --rtol 0 --atol-frac 0` reporting
`RESULT: PASS` over all 10 steps (14 fields and all integer/real stats per
step), i.e. bitwise equality with rt_serial. The 8-rank runs with 8^3 blocks
have 8 blocks per rank; the 8-rank 16^3 runs have exactly one block per
rank. The `Kokkos::abort` guards (ray outside its home block after
migration, cast reads/deposits outside the 4-cell ghost zone, ray moved
beyond the ghost zone) never fired.

## Timings

4-core VM, shared with another agent: every run below waited until no
foreign compile/benchmark process existed and was repeated if one appeared
meanwhile; best of 2 runs. `MANUAL1D_LOOP` build, `--bind-to none`,
`OMP_PROC_BIND=false`, no output. These are small problems (32³), so fixed
per-substep costs weigh heavily; treat the numbers as indicative.

**RT, `lognormal_16src_smallbuf --set legacy_rt=0`** (10 steps, 120
substeps, 1 125 901 ray-substeps; wall = time loop only)

| run | wall [s] | ns/ray-substep x threads x ranks |
|---|---|---|
| rt_serial, 1 thread | **1.93** | 1716 |
| rt_serial, 4 threads | 0.72 | 2566 |
| rt_parthenon 1 rank x 1 thread, 16³ blocks (8 blocks) | **7.21** | 6402 |
| rt_parthenon 1 rank x 1 thread, one 32³ block | 4.82 | 4283 |
| rt_parthenon 4 ranks x 1 thread, 16³ blocks | **2.70** | 9608 |
| rt_parthenon 4 ranks x 1 thread, 8³ blocks (64 blocks) | 6.01 | 21358 |
| rt_parthenon 1 rank x 4 threads, 16³ blocks | **3.40** | 12087 |
| rt_parthenon 1 rank x 4 threads, one 32³ block | 2.27 | 8062 |
| rt_parthenon 2 ranks x 2 threads, 16³ blocks | 3.11 | 11043 |
| rt_parthenon_aos 1 x 1, 16³ / 32³ | 6.93 / 4.59 | |
| rt_parthenon_aos 4 ranks x 1, 16³ | 2.71 | |
| rt_parthenon_aos 1 x 4 threads, 32³ | 1.98 | |

Phase timers [s], 1 rank x 1 thread, 16³ blocks, vs rt_serial 1 thread:

| | cell | halo | hist | step1 | cast1 | fold | phys | cast3 | sort | migrate | post |
|---|---|---|---|---|---|---|---|---|---|---|---|
| rt_serial | 0.07 | - | 0.02 | 0.36 | 0.62 | - | 0.26 | 0.51 | 0.02 | - | 0.05 |
| rt_parthenon | 0.51 | 2.36 | 0.04 | 1.05 | 0.64 | 0.79 | 0.94 | 0.46 | 0.18 | 0.12 | 0.09 |

Reading: the ray kernels themselves (cast1, cast3) cost the same as in
rt_serial. The overhead is in (a) Parthenon's ghost exchange of tau/fgam
with nghost = 4 (`halo`, 1/3 of the wall time; with one 32³ block it drops
to 0.85 s), (b) our hand-written deposit fold (`fold`), which grows with the
number of blocks, and (c) the per-cell phases (cell/step1/phys), which are
2.5–7× slower than in rt_serial. The AoS variant recovers only part of (c)
(phys 0.94 -> 0.70 s, cell 0.51 -> 0.43 s). The rest is the extra memory
traffic of the 10 double tau/fgam streams (with 4 KiB-aliased strides), the
float mirror and the ghost-grown accumulators, plus per-block launches; I did
not fully pin it down (no `perf` on this VM; callgrind instruction counts
were consistent with this). With 8³ blocks (64 blocks) the halo and fold
costs double: the ghost volume is 7× the interior. Each substep also does
a few scalar device->host copies per block (ray counts for buffer sizes) and
3 allreduces; on CPU that is noise, on a GPU it would not be.

**Hydro, n = 64, 20 steps, blast**

| run | wall [s] |
|---|---|
| hydro_serial, 1 thread | **1.63** |
| hydro_serial, 4 threads | 0.49 |
| hydro_parthenon 1 x 1, one 64³ block / 32³ blocks / 16³ blocks | **1.81** / 2.04 / 2.57 |
| hydro_parthenon 4 ranks x 1, 32³ / 16³ blocks | **0.72** / 0.92 |
| hydro_parthenon 1 rank x 4 threads, 64³ / 32³ / 16³ blocks | 0.69 / **0.79** / 1.11 |
| same, Parthenon's default `SIMDFOR_LOOP` build, 64³ / 32³ blocks | 1.73 / 1.70 (no thread speedup) |

The hydro update is within ~10 % of the serial code on one block; the
difference grows with the block count through the ghost exchange (halo
0.16 s for one 64³ block, 0.82 s for 64 16³ blocks, single thread).

## Assessment notes

### What was used, and how well it fit

| feature | used for | fit |
|---|---|---|
| Meshblocks, `StateDescriptor::AddField`, `Metadata::FillGhost` | hydro `cons`, RT `tau`/`fgam` | Good. Worked first time, bitwise, including a single periodic block that is its own neighbour (26 times). |
| Ghost exchange: `StartReceiveBoundBufs` + `AddBoundaryExchangeTasks` in a `TaskCollection` | hydro every step; `tau` before cast_1, `fgam` before cast_3 | Good API. A subset of variables is exchanged by exchanging a shallow `MeshData` (`mesh_data.AddShallow(label, base, {"tau"})`). **Slow on CPU**: `SetBounds` does a per-element logical-coordinate inverse transform (forest topology) plus an `IsActive` test, and `SendBoundBufs` a sparse-threshold test (`bvals/comms/boundary_communication.cpp`). Measured about 25–30 ns per ghost value on one core (callgrind: `SetBounds` is 28 % of all RT instructions). For RT, nghost = 4 makes the ghost volume of a 16³ block 2.4× its interior, so `halo` is about 1/3 of the 1-thread RT wall time with 16³ blocks. |
| Variables (`ParArrayND`, layout `(comp,k,j,i)`, `Real` = double) | all RT cell fields (35 + 5 + 5 components) | **Poor fit for a 45-field float record.** Every component is a separate memory stream with a (n+2g)³·8-byte stride that is a multiple of 4 KiB, and there is no per-variable float storage (`PARTHENON_SINGLE_PRECISION` switches everything, swarms included). The kernels need contiguous `float*` (tau/fgam pointers in cast_1/cast_3), so a float mirror `tf(k,j,i,10)` is rebuilt after each exchange. Loading/storing only the fields a phase uses (compile-time masks) helped; `rt_parthenon_aos` (cells in a float AoS `Kokkos::View`) quantifies the rest (see Timings). |
| `Swarm` (`AddSwarm`, vector-valued `AddSwarmValue`, `AddEmptyParticles`, `MarkParticleForRemoval`/`RemoveMarkedParticles`, `SwarmDeviceContext`) | ray storage between substeps | Fine as a container. Positions x,y,z are mandatory doubles (we store the home-cell centre for routing; the authoritative float x is payload). Persistent uint64 ids are on by default (`Metadata::NoPersistentParticleIds` turns them off). |
| Swarm boundary communication (`ResetCommunication`, `Send`, `Receive` in task lists) | ray migration after cast_3 | Works, including periodic wrap (swarm periodic BC on receipt) and the self-neighbour case. Gotchas: (1) same-rank `Send` sets the target's flag directly (`bvals/bvals.cpp:107-125`), so `ResetCommunication` of *all* blocks must run in an earlier task region than any `Send`, or a message can be wiped; (2) `int` payload goes through `Real` buffers (cast, fine below 2⁵³); (3) `Swarm::LoadBuffers_` sorts **all** active particles of the block by destination with `parthenon::sort`, which on CPU is a serial `std::sort` (`utils/sort.hpp`; thrust on CUDA; host round trip on HIP) — O(N log N) per block per migration even if one ray leaves; (4) a particle may be at most half a block outside its block (`GetNeighborBlockIndex` aborts otherwise), so with our ≤3-cell moves blocks must be ≥ 8³ for one hop. |
| `SortParticlesByCell` | **not used** | Keys every slot `0..max_active_index` without an `IsActive` test (`interface/swarm.cpp:511-517`), so stale slots are counted unless the swarm is `Defrag`ed first ("not fully parallelized", docs); serial `std::sort` on CPU; finds per-cell starts with a guess-and-walk search; gives index lists, not a packed payload. Our counting sort (atomic counts, `Kokkos::parallel_scan`, atomic scatter into a `View<RayP>`) is about 20 lines and leaves the payload contiguous for `step1_*`. |
| `pmb->par_for` / `par_reduce`, plain Kokkos (`parallel_scan`, `atomic_add`, `atomic_fetch_add`, `View<RayP>`, `Kokkos::abort`) | all kernels | Good, with one trap: on host backends Parthenon's default loop pattern is `SIMDFOR_LOOP` (`src/CMakeLists.txt:41-46`, `doc/sphinx/src/par_for.rst`) -- plain C loops with an `omp simd` inner loop, **not threaded**. With it, hydro_parthenon got no speedup from 4 threads (1.70 s vs 1.59 s) and only the raw `Kokkos::parallel_*` parts of rt_parthenon scaled. This project sets `PAR_LOOP_LAYOUT=MANUAL1D_LOOP` (flat `RangePolicy`) before `add_subdirectory`; then hydro runs in 0.69–0.79 s on 4 threads. Parthenon's CPU model is clearly "one MPI rank per core". Parthenon wrappers and raw Kokkos mix without friction; `par_for` bounds are inclusive. |
| Task lists | ghost exchange, migration | Usable, but our substep is a strict sequence of global phases with two allreduces, so tasking buys little; phases are driven from a plain host loop over blocks. No `EvolutionDriver`, integrator or outputs are used. |
| Packs (`SparsePack`, `MeshBlockPack`, `SwarmPack`) | **not used** | Kernels are launched per block (8–64 blocks × ~15 kernels per substep). For a GPU build they should be packs over `MeshData` (one launch for all blocks), which needs the per-block scratch (rays, CSR offsets, accumulators) to become pack-indexed as well. |
| Reverse ("ghost → owner, summed") exchange | deposits | **Absent.** Checked `bvals/comms/*` and `BoundaryType` (`basic_types.hpp:70`: `local, nonlocal, any, flxcor_send/recv, gmg_*`): flux correction is face-only and coarse–fine; nothing sums ghost data into the owner. Written by hand (`fold_deposits` + slab helpers, ~100 lines) from `pmb->neighbors` (`rank`, `gid`, `offsets`): same-rank neighbours by direct add, others by `MPI_Isend/Irecv` with tag = receiver gid·27 + direction. Integer adds, so the order does not matter. Same-level only; on AMR it would need restriction/prolongation of deposits. |
| Global reductions | histogram, statistics, ray totals, dt | Plain `MPI_Allreduce` on integers (exact), and `MAX` on doubles for hydro dt. Parthenon's `AllReduce<T>` (task-based `MPI_Iallreduce`) offers no reproducible floating-point sum. |

### Gotchas met

* Default host loop pattern `SIMDFOR_LOOP` is serial (see the table).
* `ParameterInput::ModifyFromCmdline` aborts on any argument that contains
  `=` but no `/` (`parameter_input.cpp:415`), e.g. our `--set legacy_rt=0`.
  We pass Parthenon a filtered argv (`sp::parthenon_args`).
* `Kokkos::view_alloc(WithoutInitializing, const char*)` treats the
  `const char*` as a pointer to user memory; use `std::string`.
* The global cell offset of a block is computed from
  `block_size.xmin` (the forest `LogicalLocation` is per tree).
* OpenMPI binds ranks to cores by default; with OpenMP threads use
  `--bind-to none` (2 threads per rank on one core ran 8–20× slower), and
  `OMP_WAIT_POLICY=passive` when oversubscribing.
* The spike/common of the base commit does not compile for CUDA
  (`radray::cp::freq_init` in device code); the version on
  `framework-spike` (consts_hd.hpp) does.

### GPU plausibility

A compile-only Kokkos CUDA build (nvcc 12.0 + g++-12, sm_80, Kokkos
CUDA + OpenMP, Parthenon headers configured for CUDA) of both TUs succeeds
against the `framework-spike` spike/common, with no warnings in our code. It
would not run well as is: (a) `fold_deposits` hands `View::data()` to MPI
(needs GPU-aware MPI or host staging); (b) several scalar device→host copies
per block per substep (ray counts for buffer sizing) are synchronisation
points; (c) per-block launches instead of packs; (d) `step1_*` runs one
thread per cell with a per-cell heap sort (divergent; fine on CPU);
(e) `parthenon::sort` inside swarm `Send` is thrust on CUDA but a host
round trip on HIP and unsupported on SYCL.

### Docs and sources consulted

`doc/sphinx/src/{particles,boundary_communication,tasks,reductions,
load_balancing,solvers}.rst`; `example/particles` (swarm creation,
transport, iterative Send/Receive), `example/advection` (driver, partitions,
`AddBoundaryExchangeTasks`); sources `src/interface/swarm{,_comms,_container,
_device_context}.{hpp,cpp}`, `src/bvals/bvals.cpp`,
`src/bvals/boundary_conditions{,_generic}.{hpp,cpp}`,
`src/bvals/comms/{bvals_in_one.hpp,boundary_communication.cpp}`,
`src/bvals/neighbor_block.hpp`, `src/mesh/mesh-gmg.cpp` (neighbour lists),
`src/interface/data_collection.hpp`, `src/parameter_input.cpp`,
`src/utils/sort.hpp`, `src/utils/cell_center_offsets.hpp`.

### Lines of code (non-blank, non-comment)

| file | lines |
|---|---|
| `rt_parthenon.cpp` | ~840 (952 incl. comments/blank) |
| `hydro_parthenon.cpp` | ~205 (235) |
| `parth_common.hpp` | ~116 (155) |
| `CMakeLists.txt` | ~39 (61) |
| for comparison: `spike/serial/rt_serial.cpp`, `hydro_serial.cpp` | 269, 110 |

About 100 lines of rt_parthenon are the hand-written reverse exchange,
~60 the swarm <-> sorted-array plumbing, ~80 the CellData <-> variables
mapping (incl. the AoS variant), ~130 setup/argument handling, ~70 the
output gather/timers; the rest is the phase sequence of rt_serial.

### What Parthenon offers for the later plans

* **Octree AMR with refinement tags**: `parthenon/mesh/refinement =
  adaptive|static`; per-package `CheckRefinementBlock`/`CheckRefinementMesh`
  returning `AmrTag` plus built-in `derivative_order_1/2` criteria
  (`src/amr_criteria/`), per-variable prolongation/restriction, a forest of
  octrees (`src/mesh/forest/`). Mature (used by AthenaPK, KHARMA, Phoebus,
  Artemis — from memory, not checked here).
* **Swarms on AMR**: the docs say *"AMR is currently not supported, but
  support will be added in the future"* (`doc/sphinx/src/particles.rst`,
  section Communication), and `src/mesh/mesh-amr_loadbalance.cpp` contains no
  swarm code, so particles are presumably not carried along when blocks are
  refined or load-balanced **[not tested]**. `SetNeighborIndices_` already
  tests neighbour regions at half-block resolution, so finer neighbours are
  anticipated in the routing. For adaptive RT with rays in a swarm this is
  the main gap: we would have to move rays on refinement/rebalancing
  ourselves.
* **Load balancing**: per-block cost via `MeshBlock::SetCostForLoadBalancing`
  with `<parthenon/loadbalancing> balancer = manual`
  (`doc/sphinx/src/load_balancing.rst`); matches the `w_brick` cost model of
  design §3.1, but see the swarm caveat above. No timer-based balancing yet.
* **Gravity / multigrid**: `src/solvers/{mg_solver,cg_solver,bicgstab_solver}.hpp`
  (geometric MG, CG and BiCGSTAB with MG preconditioning), demonstrated in
  `example/poisson_gmg` with static refinement and `multigrid = true`
  (two-level composite grids, `src/mesh/mesh-gmg.cpp`). Their dot products
  are `AllReduce<Real>` (`MPI_Iallreduce`), so iteration counts/results are
  not bitwise independent of the decomposition.
* **Determinism**: ghost exchange is a pure copy (×1) and hence bitwise.
  There is no reproducible floating-point reduction; Kokkos `Real`
  reductions depend on thread order. Bitwise reproducibility (our rule) has
  to come from our integer fixed-point design, as here.
* **Kokkos backends**: bundled Kokkos 4.6 (Serial, OpenMP, CUDA, HIP, SYCL);
  Parthenon documents CUDA machine configs and CUDA-aware MPI
  (`doc/sphinx/src/building.rst`); `parthenon::sort` is CPU/CUDA/HIP only.

