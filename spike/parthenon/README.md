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
  examples/tests off, `Real = double`. Our targets get
  `-ffp-contract=off -fno-fast-math` and `-O2` (as `spike/serial`);
  `../cmake/radray_core.cmake` provides `radray_core`/`spike_common`.
* Needs MPI (OpenMPI here); Parthenon enables it automatically when found.

Build times on this 4-core VM (shared with another agent's builds, so
upper bounds):

| what | wall | CPU |
|---|---|---|
| Parthenon + Kokkos library only, inside this project, `-j2` | 5 min 23 s | 9 min 27 s user |
| Parthenon + Kokkos + all examples, standalone, `-j4` (`/home/user/ext/build-parthenon`) | 7 min 27 s | 24 min 19 s user |
| `rt_parthenon` (compile + link) | 40 s | |
| `hydro_parthenon` (compile + link) | 31 s | |
| compile-only Kokkos **CUDA** build of both TUs (nvcc 12.0, sm_80, configure 24 s) | 144 s + 115 s | |

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
  kernel loads/stores only the `CellData` fields it uses (compile-time bit
  masks), because every component of a Parthenon variable is its own memory
  stream.
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

Results (all on the final code, `/tmp/.../final.txt` of this session):

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

"PASS" for RT is `compare.py ... --rtol 0 --atol-frac 0` reporting
`RESULT: PASS` over all 10 steps (14 fields and all integer/real stats per
step), i.e. bitwise equality with rt_serial. The 8-rank runs with 8^3 blocks
have 8 blocks per rank; the 8-rank 16^3 runs have exactly one block per
rank. The `Kokkos::abort` guards (ray outside its home block after
migration, cast reads/deposits outside the 4-cell ghost zone, ray moved
beyond the ghost zone) never fired.

## Timings

@TIMINGS@

## Assessment notes

@NOTES@
