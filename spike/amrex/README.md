# AMReX prototype of the spike (hydro_amrex, rt_amrex)

This is the AMReX half of the framework spike. It runs the framework-independent kernels in
`spike/common` (unchanged) on AMReX data structures, distributed over MPI
ranks and OpenMP threads. It must reproduce `spike/serial` **bit for bit**.

| executable | reproduces | AMReX pieces |
|---|---|---|
| `hydro_amrex` | `spike/serial/hydro_serial` (unsplit MUSCL-Hancock + HLLC, `hydro_mh.hpp`) | `Geometry` (periodic), `MultiFab` U (5 comps, 2 ghosts), predictor states Q in `FabArray<BaseFab<QCell>>` (20 doubles per cell, 1 ghost; `-DHYDRO_Q_SOA`: 20-comp `MultiFab`), 3 face-centred flux `MultiFab`s; `FillBoundary`; `MFIter` tiling + OpenMP; `ParallelFor`; `ReduceOps<ReduceOpMax>` + `ParallelDescriptor::ReduceRealMax`; `ParallelCopy` to rank 0 for output |
| `rt_amrex` | `spike/serial/rt_serial` (sorted-ray-array RT substep, `ray_kernels.hpp`, `rt_cell.hpp`, `rrhost::RTRun`) | `FabArray<BaseFab<rrk::CellData>>` (cells), `FabArray<BaseFab<float>>` (tau/fgam, 4 ghosts, `FillBoundary`), `FabArray<BaseFab<Long>>` (13 int64 deposit accumulators, 4 ghosts, `HostDevice::Atomic::Add` + `SumBoundary`), `ParticleContainerPureSoA<21,2>` (rays) + `DenseBins` (cell sort) + local `Redistribute` (migration), `Scan::ExclusiveSum`, `ParallelForOMP`, `Reduce::Sum`, `ParallelDescriptor::ReduceLongSum` |

## Build

### AMReX (once)

The prototype needs AMReX 26.10 or later, built 3-D with MPI, OpenMP and particles in **single**
precision. Positions and the ray payload are `float`, as in `rrk::RayP`. On this machine
it is installed at `/home/user/ext/inst-amrex`. The from-scratch build took
3 min 3 s wall (10.5 min CPU) with `-j4`:

```
git clone --branch 26.10 --depth 1 https://github.com/AMReX-Codes/amrex.git amrex
cmake -S amrex -B build-amrex -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX=$HOME/ext/inst-amrex \
      -DAMReX_SPACEDIM=3 -DAMReX_MPI=ON -DAMReX_OMP=ON \
      -DAMReX_PARTICLES=ON -DAMReX_PARTICLES_PRECISION=SINGLE \
      -DAMReX_LINEAR_SOLVERS=ON -DAMReX_FORTRAN=OFF -DAMReX_FASTMATH=OFF
cmake --build build-amrex -j4 && cmake --install build-amrex
```

`AMReX_LINEAR_SOLVERS` is not needed here. It was built because it is wanted later
(MLMG). AMReX 26.x requires **C++20**; our two targets set `CXX_STANDARD 20`, while
`radray_core` stays C++17.

### The prototype

```
cmake -S spike/amrex -B build-spike-amrex -DAMReX_DIR=/home/user/ext/inst-amrex/lib/cmake/AMReX
cmake --build build-spike-amrex -j2
# reference:
cmake -S spike/serial -B build-spike-serial && cmake --build build-spike-serial -j2
```

* `-DSPIKE_FETCH_AMREX=ON` (with no `AMReX_DIR`) instead downloads and builds AMReX
  26.10 with the options above via `FetchContent`. That path was not tested,
  since the installed AMReX was used.
* Our targets get `-O2 -ffp-contract=off -fno-fast-math` (same as `spike/serial`).
  They also link `AMReX::Flags_INLINE` (`-finline-limit=43210`); see gotchas.
* `-DRT_CELL_SOA` (add to `CMAKE_CXX_FLAGS`) switches the RT cell storage to a
  45-component float SoA `FabArray`, which is the GPU-friendly layout. `-DHYDRO_Q_SOA` does the same for the hydro predictor states.
* A clean build of the prototype takes 21 s at -j2, including `radray_core`.
  Compiling `rt_amrex.cpp` takes 16.5 s, `hydro_amrex.cpp` 7.0 s, and
  `rt_serial.cpp` 2.8 s for comparison. `rt_amrex.cpp` preprocesses to 208 k lines;
  `rt_serial.cpp` preprocesses to 75 k.

## Run

```
hydro_amrex --n 64 --steps 20 --ic blast|sod --cfl 0.4 --out U.bin --max_grid_size 32
rt_amrex <problem_dir> <outdir> [--nsteps n] [--set key=value]... [--no-output]
         [--max_grid_size 16] [--redist_local 3] [--ghosts 4]
mpirun --oversubscribe --allow-run-as-root -np 4 -x OMP_NUM_THREADS=1 rt_amrex ...
```

* Output formats are those of the serial codes. hydro writes raw float64 (5,n,n,n).
  RT writes `step_NNNN/*.bin` + `stats.txt` via `RTRun::write_step` on rank 0.
* `--redist_local L` sets the maximum number of cells a ray may move per substep, which
  enables AMReX's neighbour-only `Redistribute`. `0` means a global `Redistribute`.
* `--ghosts G` sets the ghost width of tau/fgam and of the deposit accumulators. The
  default is 4 (as specified). The design's reach estimate is 2 (§3.2), and 2 also
  passes. Every ghost access is range-checked, and the run aborts if a target falls
  outside the ghost region.

## Validation

`spike/amrex/validate.sh <serial_build> <amrex_build> <workdir>` runs the whole matrix:

It ran on 4 cores shared with another agent's builds. All results are bitwise:
`cmp` for hydro, and `compare.py --rtol 0 --atol-frac 0` for RT (every field of every
step, plus all stats).

```
spike/amrex/validate.sh build-spike-serial build-spike-amrex /tmp/val       # -> SUMMARY: 75 PASS, 0 FAIL
spike/amrex/validate.sh build-spike-serial build-spike-amrex /tmp/val hydro # re-run after the AoS-Q change: 24 PASS, 0 FAIL
```

| case | ranks | max_grid_size | OMP threads | result |
|---|---|---|---|---|
| hydro blast, n=64, 20 steps (`cmp`) | 1, 2, 4 | 32, 16 | 1, 2 | 12/12 PASS |
| hydro sod, n=64, 20 steps (`cmp`) | 1, 2, 4 | 32, 16 | 1, 2 | 12/12 PASS |
| RT `uniform_8src`, `legacy_rt=0` | 1 (mgs 32), 1, 2, 4, 8 | 32; 16, 8 | 1, 2 | 17/17 PASS |
| RT `lognormal_16src_smallbuf`, `legacy_rt=0` | 1 (mgs 32), 1, 2, 4, 8 | 32; 16, 8 | 1, 2 | 17/17 PASS |
| RT `lognormal_16src_smallbuf`, `legacy_rt=0 legacy_physics=0` | 1 (mgs 32), 1, 2, 4, 8 | 32; 16, 8 | 1, 2 | 17/17 PASS |

Extra spot checks, all on `lognormal_16src_smallbuf` and all PASS:
* `--ghosts 2` with 4 ranks / mgs 8 and with 1 rank x 2 threads / mgs 16;
* global Redistribute (`--redist_local 0`) with 4 ranks;
* the SoA build (`-DRT_CELL_SOA -DHYDRO_Q_SOA`) with 1 rank x 2 threads and 4 ranks
  (RT), and 2 ranks / mgs 16 (hydro `cmp`);
* 8 ranks with `--max_grid_size 4`, i.e. 512 boxes of 4^3, where the 4 ghost layers
  are as wide as the box;
* 1 rank x 4 threads with mgs 32.

For the record, the first runs of `rt_amrex` passed `compare.py` on the first try
with a global Redistribute. Every bug found later was a performance or
Redistribute-placement issue (see gotchas), never a numerical difference. The
order-independent design (fixed-point deposits, canonical merge order) did its job.


## Timings

Best of 5 runs, `--no-output`, on 4 vCPUs shared with another agent. Single-thread numbers
are stable to about 5%. Multi-thread and multi-rank numbers vary up to 3x between runs
because of host contention, so treat them as indicative. "ns/r-s" is wall
x threads x ranks per ray-substep. Both codes do 1 125 901 ray-substeps on this problem.

**RT, `lognormal_16src_smallbuf --set legacy_rt=0` (N=32, 10 steps)**

| run | wall [s] | ns/r-s | cell | hist | step1 | cast1 | phys | cast3 | sort | migrate | post |
|---|---|---|---|---|---|---|---|---|---|---|---|
| rt_serial, 1 thread | **1.80** | 1598 | 0.062 | 0.019 | 0.343 | 0.587 | 0.242 | 0.471 | 0.016 | - | 0.053 |
| rt_serial, 4 threads | **0.71** | 2514 | 0.027 | 0.014 | 0.139 | 0.226 | 0.091 | 0.167 | 0.019 | - | 0.019 |
| rt_amrex 1 rank x 1 thread, mgs 16 (8 boxes) | **3.26** | 2899 | 0.200 | 0.032 | 0.527 | 1.481 | 0.387 | 0.469 | 0.059 | 0.035 | 0.058 |
| rt_amrex 1x1, mgs 32 (1 box) | 2.68 | 2378 | 0.155 | 0.025 | 0.553 | 1.004 | 0.329 | 0.463 | 0.060 | 0.015 | 0.058 |
| rt_amrex 1x1, mgs 16, `--ghosts 2` | 2.53 | 2251 | 0.158 | 0.026 | 0.443 | 0.923 | 0.341 | 0.471 | 0.062 | 0.035 | 0.060 |
| rt_amrex 4 ranks x 1 thread, mgs 16 | **1.52** | 5392 | 0.087 | 0.031 | 0.197 | 0.678 | 0.205 | 0.230 | 0.030 | 0.030 | 0.024 |
| rt_amrex 4x1, mgs 16, `--ghosts 2` | 1.12 | 3961 | 0.066 | 0.026 | 0.176 | 0.399 | 0.154 | 0.213 | 0.022 | 0.027 | 0.025 |
| rt_amrex 1 rank x 4 threads, mgs 16 | **1.51** | 5362 | 0.090 | 0.051 | 0.242 | 0.653 | 0.192 | 0.173 | 0.060 | 0.025 | 0.019 |
| rt_amrex 1x4, mgs 32 | 1.17 | 4161 | 0.056 | 0.047 | 0.203 | 0.456 | 0.157 | 0.156 | 0.055 | 0.015 | 0.018 |

Communication inside those phases, from `rt_amrex`'s third timer line:

| run | FillBoundary(tau, fgam) | SumBoundary(deposits) | setVal(deposits incl. ghosts) |
|---|---|---|---|
| 1x1 mgs 16 | 0.121 | 0.785 | 0.126 |
| 1x1 mgs 32 | 0.046 | 0.287 | 0.061 |
| 1x1 mgs 16 ghosts 2 | 0.052 | 0.284 | 0.067 |
| 4x1 mgs 16 | 0.160 | 0.357 | 0.032 |
| 4x1 mgs 16 ghosts 2 | 0.115 | 0.150 | 0.015 |
| 1x4 mgs 16 | 0.055 | 0.403 | 0.044 |

Reading the numbers:
* `rt_amrex` is **1.5-1.8x slower than `rt_serial` on one core**. The ray kernels
  themselves cost the same (cast3 0.47 vs 0.47 s; step1 is +30-50% for the per-box
  CSR/scan bookkeeping and the strided cell reads). The excess is mesh-side:
  * `SumBoundary` of the 13 int64 accumulators. Its cost scales with ghost volume:
    a 16^3 box with 4 ghosts has 2.4x as many ghost cells as valid cells.
  * Zeroing the ghost-grown accumulators.
  * Cell load/store between the FAB and the kernels' `CellData` struct.
* With the reach-2 ghost width of the design (`--ghosts 2`), the gap shrinks to 1.4x.
* `hist` does not scale with threads: one shared int64 histogram with OpenMP atomics.
  A per-thread or per-box histogram would fix it.
* Scaling 1 to 4 ranks (or threads) is 2.1-2.3x at mgs 16 (2.3x with 2 ghosts),
  against 2.5x for `rt_serial`'s OpenMP. N=32 is tiny: per substep there are only
  ~7 k rays and 32 k cells, so fixed costs dominate.

**hydro, n=64, 20 steps, blast**

| run | wall [s] | ns/cell-update | FillBoundary | maxspeed | predict | flux | update |
|---|---|---|---|---|---|---|---|
| hydro_serial, 1 thread | **1.41** | 269 | | | | | |
| hydro_serial, 4 threads | **0.46** | 87 | | | | | |
| hydro_amrex 1 rank x 1 thread, mgs 32 (8 boxes) | **1.56** | 298 | 0.034 | 0.056 | 0.750 | 0.568 | 0.155 |
| hydro_amrex 1x1, mgs 64 (1 box) | 1.50 | 286 | 0.017 | 0.057 | 0.698 | 0.585 | 0.143 |
| hydro_amrex 4 ranks x 1 thread, mgs 32 | **0.46** | 88 | 0.060 | 0.018 | 0.211 | 0.136 | 0.033 |
| hydro_amrex 1 rank x 4 threads, mgs 32 | **0.61** | 117 | 0.015 | 0.051 | 0.277 | 0.210 | 0.058 |

`hydro_amrex` runs within 10% of `hydro_serial` once the predictor states are AoS
(`FabArray<BaseFab<QCell>>`). With AMReX's native SoA `MultiFab` for Q
(`-DHYDRO_Q_SOA`), the flux pass alone took 1.43 s instead of 0.57 s, and the code
was about 1.7x slower in total (best of 3: 2.44 s). The max-speed reduction loop is
not threaded (a `ReduceOps` loop outside an OpenMP region), which partly explains
why 1x4 is slower than 4x1.


## Notes for the framework assessment

### AMReX features used, and how well they fit

| feature | used for | fit |
|---|---|---|
| `Geometry`, `BoxArray::maxSize`, `DistributionMapping` | periodic N^3 domain, decomposition | Trivial. The default DM is SFC/knapsack by box count. `DistributionMapping::makeKnapSack(weights)` / `makeSFC` (`AMReX_DistributionMapping.H:307`) exist for the cost-weighted bricks of design §3.1 but were not tried. |
| `MultiFab`; `FabArray<BaseFab<T>>` with T = `float`, `Long`, our structs | hydro state; RT tau/fgam (float), deposits (int64), cells (`BaseFab<rrk::CellData>`), hydro Q (`BaseFab<QCell>`) | Good. Any trivially copyable T works for storage, `Array4` access and `ParallelFor`. Communication (`ParallelCopy`, `FillBoundary`) needs an arithmetic T (gotcha 5). |
| `FillBoundary(scomp, ncomp, period)` | tau before cast_1, fgam before cast_3, hydro U | Exactly what was needed. Component subsets avoid moving all 10 components twice. |
| `SumBoundary(period)` on `FabArray<BaseFab<Long>>` | fold ghost deposits into their owners | Correct out of the box for int64, so exact and decomposition-independent. It is the single biggest CPU overhead (see timings): every ghost tag is copied into temporary FABs and then added (`FB_local_add_cpu`, `AMReX_FBI.H:455`), so cost scales with ghost volume. The `deterministic` flag (`AMReX_FBI.H:650, 1243`) matters on GPU for floats and is irrelevant for our integers. |
| `HostDevice::Atomic::Add` | int64 deposits (incl. ghost cells), histogram, stats | Fine: `#pragma omp atomic` on the host, device atomics on the GPU. Plain `Gpu::Atomic::AddNoRet` is *not* atomic between host OpenMP threads, so `HostDevice::` is the right call. |
| `ParticleContainerPureSoA<21,2>` | rays: 3 position + 18 payload floats; nside, ipix ints | A good fit for "a ray has a home cell". Positions and real components are `ParticleReal`, a *global* build option (gotcha 8). Ray IDs are not needed: every ray gets id 1, since `NextID()` would add a serial counter. |
| `Redistribute(0,0,0,local=3)` | migration after cast_3 | Works and is cheap: `migrate` is 0.015-0.035 s of 1.1-3.3 s. The local (neighbour-only) algorithm needs *unwrapped* positions (gotcha 2). The global algorithm also passes (`--redist_local 0`). |
| `DenseBins<int>` (`BinPolicy::OpenMP` on CPU) | counting sort of each box's rays by local cell to CSR offsets and permutation | Fits directly: offsets = our `start[]`, permutation = gather order. We still gather into an AoS `RayP` scratch every substep (80 B/ray). `ReorderParticles` / `SortParticlesByCell` (`AMReX_ParticleContainer.H:548-577`) could keep the tile itself sorted, but the kernels want AoS `RayP` anyway. |
| `MFIter` tiling + `#pragma omp parallel` + `ParallelFor` | hydro | The standard pattern; bitwise identical for any tiling and thread count. `ParallelFor` on the CPU is a serial loop: threads come only from the OpenMP region around a tiled `MFIter` (gotcha 6). |
| `ParallelForOMP` (26.x, `AMReX_GpuLaunch.H:310-390`) | all RT loops (boxes in sequence, threads inside a box) | Useful: OpenMP inside a box on the CPU, `ParallelFor` on GPU. Static schedule only (the serial reference uses `dynamic` for step1). Not mentioned in the sphinx docs. |
| `Scan::ExclusiveSum`, `Reduce::Sum`, `ReduceOps<ReduceOpMax>` | output offsets, survivor compaction, ngam_tot, hydro max speed | Fine. Integer scans and sums are exact. The CPU `Scan` is serial. |
| `ParallelDescriptor::ReduceLongSum / ReduceRealMax / ReduceIntMax` | histogram (2101 int64), ray stats, ngam_tot, dt, error flag | Thin `MPI_Allreduce` wrappers (`AMReX_ParallelDescriptor.H:1321`). Integer sum and max are exact, so decomposition-independent. `ReduceRealSum` is an ordinary MPI float sum (not reproducible across rank counts); it is never used here. |

**Lines of code** (non-blank, non-comment):

| file | lines | serial counterpart |
|---|---|---|
| `rt_amrex.cpp` | 635 | `rt_serial.cpp`: 269 |
| `hydro_amrex.cpp` | 176 | `hydro_serial.cpp`: 110 |

Of the ~370 extra lines in `rt_amrex`:
* ~60 are the AoS/SoA cell-layout switch;
* ~45 are RayP <-> SoA particle load/store;
* ~40 are the periodic-image index mapping for ghost reads and deposits;
* ~50 are per-box scratch, CSR and survivor compaction;
* the rest are timers, error checks and output gathering.

No kernel had to change; `spike/common` is used untouched. Build: AMReX takes 3 min
at -j4, the prototype 21 s at -j2.

### Gotchas (in the order we hit them)

1. **Name clash.** `amrex::CellData<T>` exists: the per-cell component view returned
   by `Array4::cellData(i,j,k)` (`AMReX_Array4.H:46`). With `using namespace amrex;`,
   `rrk::CellData` becomes ambiguous. We alias it as `RTCell`.
2. **Local `Redistribute` silently misplaces particles whose positions were wrapped
   by hand.** A ray crossing the periodic boundary from cell 31 to cell 0 looks like a
   31-cell jump. With `local=L` smaller than the box size, AMReX left such rays in a
   box that does not contain them, *without any error*. Only our own "home cell inside
   this box" check in the gather caught it. Fix: store the *unwrapped* image of the new
   home cell nearest to the emitting cell (it may be just outside the domain), and let
   `Redistribute` apply the periodic shift. That shift is exact in float.
3. **AMReX does not pass its inlining flag to consumers.** AMReX itself is built with
   `-finline-limit=43210` (GCC), but `AMReX::amrex` does not propagate it. The
   communication copy/add loops are header templates instantiated in *our*
   translation units, and without the flag they ran ~2x slower (`SumBoundary`
   1.9-2.0 s vs 0.86 s in one A/B run). We link `AMReX::Flags_INLINE`; it does not
   change results.
4. **SoA FABs vs struct-per-cell kernels.** AMReX FABs are component-slowest
   (Fortran order). Our per-cell kernels (`calc_rt_1_cell`, `physics_cell`,
   `predict_cell`, `face_flux`) take a contiguous struct or array per cell. With 45
   float components, every cell load/store touches 45 streams 16-128 kB apart:
   * the RT cell phases were 3-5x slower and the hydro flux pass 2.5x slower than
     the AoS serial code;
   * `FabArray<BaseFab<struct>>` (1 component) restores AoS and works with
     `ParallelFor`/`Array4`, so it is the default here;
   * on a GPU the SoA layout (`-DRT_CELL_SOA`, `-DHYDRO_Q_SOA`; both validated) is
     the coalesced one, unless the kernels are restructured to touch only what they
     use.
5. **`ParallelCopy`/`FillBoundary` do not compile for `BaseFab<struct>`.** They
   instantiate the `FabArrayBase::ADD` path, which needs `operator+=` on the element
   type. Output therefore goes AoS FAB, then a 45-component float FAB, then
   `ParallelCopy` to rank 0.
6. **`ParallelFor` is serial on the CPU.** OpenMP needs `#pragma omp parallel` plus
   a tiled `MFIter`, or `ParallelForOMP`; never nest the two. A `ReduceOps` loop
   outside an OpenMP region is single-threaded (the hydro max-speed loop).
7. **`DenseBins::build(n, v, Box, f)` clamps out-of-box indices** to the box edge
   (`AMReX_DenseBins.H:430-447`) instead of failing. We use the 1-D `nbins` overload
   with our own range check.
8. **Particle precision is global.** `AMReX_PARTICLES_PRECISION` fixes `ParticleReal`
   for *all* particle containers in the build (`AMReX_REAL.H:91`). Rays want float
   payloads; dark-matter particles for PM would want double positions. In one build,
   one of the two must store its data in int/runtime components or pay 2x memory.
9. AMReX 26.x **requires C++20** (`std::integral` and `requires` in public headers)
   and **nvcc >= 12.2**. Our targets are C++20; `radray_core` stays C++17.
10. `amrex::Initialize(argc, argv, false, MPI_COMM_WORLD, f)` (no ParmParse of argv)
    keeps our own command line. Set `amrex.verbose=0` in the callback to silence the
    banner.
11. `mpirun` here needs `--oversubscribe --allow-run-as-root`. The validation uses
    `--bind-to none` so a rank's OpenMP threads are not pinned to a single core.

### GPU readiness (honest)

**Not compiled for GPU.** AMReX 26.10 refuses the installed nvcc 12.0 at configure
time (minimum 12.2), so even a compile-only check was impossible here.

As written:
* **Device-style already.** Every per-cell and per-ray loop is an
  `AMREX_GPU_DEVICE` lambda through `ParallelForOMP`/`ParallelFor`. Deposits use
  device atomics on `Long`; `DenseBins` falls back to its GPU policy; `Scan` and
  `Reduce::Sum` are used throughout. Scratch is `Gpu::DeviceVector`, and the error
  flag and histogram come back via `Gpu::copy`. `spike/common` compiles as CUDA on
  the `framework-spike` branch.
* **Would not compile for GPU / is host-only.**
  * The initial fill: `LoopOnCpu` over a host `radray::GasGrid`, and sources written
    through host `Array4`s.
  * The rank-0 output loop, which reads a device FAB on the host.

  Both need a pinned host FAB plus a copy, or `The_Managed_Arena`.
* **Would compile but fail at run time.** `rrk::PhysCtx` holds *host* pointers to
  the rate tables (`RTRun::phys()`); they must be copied to device memory once.
* **GPU performance issues to expect.**
  * `step1_count`/`step1_write` run one thread per cell, so source cells (192
    emitted rays, a heap sort of merges) serialize. A team or warp per cell is
    needed.
  * The AoS cell FAB and the AoS `RayP` scratch are uncoalesced.
  * The int64 deposit atomics from adjacent threads hit the same few cells
    (design §2.4).
  * `S.nout`/`nkeep` from `Scan` mean one device-to-host sync per box per phase.

### What AMReX offers for the later plans (pointers into `/home/user/ext/amrex/Src`)

* **AMR with refinement on overdensity.**
  * `AmrCore` (`AmrCore/AMReX_AmrCore.H`, `AMReX_AmrMesh.H`) manages the hierarchy,
    regridding and Berger-Rigoutsos clustering (`AMReX_Cluster.H`).
  * Tagging uses `ErrorList`/`TagBox` (`AMReX_ErrorList.H`, `AMReX_TagBox.H`); an
    overdensity criterion is a few lines.
  * `FillPatchUtil` and `Interpolater` fill coarse-fine ghosts; `FluxRegister`
    handles refluxing.
* **Subcycling.** Time subcycling is built into the older `Amr`/`AmrLevel` framework
  (`Amr/AMReX_Amr.H:117-124`: `subCycle()`, `subcyclingMode()`, per-level cycle
  counts). With `AmrCore` the application writes its own (short) recursive timestep.
* **Particles on AMR levels.**
  * `AmrParticleContainer` (`AmrCore/AMReX_AmrParticles.H:294`) assigns particles to
    the finest level containing them. `Redistribute(lev_min, lev_max, nGrow)` has
    the `nGrow` slack meant for subcycling.
  * `Particle/AMReX_ParticleMesh.H` provides generic `ParticleToMesh`/`MeshToParticle`
    (CIC deposit/interpolate). `NeighborParticles` and `DenseBins`/`SparseBins` also
    exist.
  * An AMR version of the rays would keep them in an `AmrParticleContainer`, with
    per-level ghost deposits plus `SumBoundary`. The coarse/fine sync of an
    arbitrary deposit is *not* provided; we would write it (cf. `ParticleToMesh` and
    `average_down`).
* **Gravity.**
  * MLMG (`LinearSolvers/MLMG/AMReX_MLMG.H`, `AMReX_MLPoisson.H`,
    `AMReX_MLNodeLaplacian.H`) is geometric multigrid with composite multi-level
    solves, periodic BCs and GPU support: the natural multigrid gravity for the
    roadmap. Free-space BCs are in `LinearSolvers/OpenBC` and
    `FFT/AMReX_FFT_OpenBCSolver.H`.
  * `FFT/AMReX_FFT_Poisson.H` (needs `AMReX_FFT=ON` and FFTW/cuFFT; off in this
    build) could replace our PM FFT on the uniform level.
* **Nyx heritage.** Nyx is AMReX's cosmology code: comoving hydro, dark-matter
  particles, MLMG gravity, AMR with subcycling. It is the closest existing analogue
  of an AMR RadHydro. It is *not* in this source tree and was not inspected, so
  anything said about its internals (time-stepping, PM/MLMG split, halo finding) is
  from memory and **unverified**. It is the first thing to read before an AMR port.
* **Determinism of its reductions.**
  * Integer `FillBoundary`, `SumBoundary` and `ReduceLongSum` are exact.
  * Float reductions (`MultiFab::sum`, `ReduceRealSum` = `MPI_Allreduce`) are
    reproducible only for a fixed rank count and decomposition.
  * On GPU, float ADD-communication uses atomics unless `deterministic=true`
    (`AMReX_FBI.H:650`).
  * Our fixed-point int64 design sidesteps all of this: bitwise decomposition
    independence needed no AMReX-side special handling.

### Docs consulted

* Sphinx docs: `/home/user/ext/amrex/Docs/sphinx_documentation/source/Particle.rst`
  (PureSoA, particle tiling, local `Redistribute`).
* Header doc comments:
  * `Src/Base/AMReX_FabArray.H` (`SumBoundary`, `OverrideSync`);
  * `AMReX_GpuLaunch.H` (`ParallelForOMP`), `AMReX_GpuAtomic.H`, `AMReX_Reduce.H`,
    `AMReX_Scan.H`, `AMReX_FBI.H`;
  * `Src/Particle/AMReX_DenseBins.H`, `AMReX_ParticleContainer.H`,
    `AMReX_ParticleTile.H`.
* Example code: `Tests/Particles/RedistributeSOA/main.cpp`, the most useful PureSoA
  example.

The sphinx docs lag the code: `ParallelForOMP` appears in no `.rst` file.

### Verdict (for the comparison with Parthenon)

* **Coverage.** AMReX covered every distributed piece of the sorted-ray-array design
  with stock features: ghost exchange, integer ghost-to-owner sums, ray migration,
  binning and reductions. Each took a few lines, and the result was bitwise
  identical to the serial reference from the first working run.
* **Costs.**
  * ~1.5x single-core overhead on the CPU at N=32, mostly `SumBoundary` over 4 ghost
    layers and the SoA-vs-AoS mismatch. It drops to 1.4x with 2 ghosts and is timed
    per phase.
  * C++20 and nvcc >= 12.2.
  * Global particle precision.
  * Two silent-failure modes that need our own assertions: local `Redistribute`
    with wrapped positions, and `DenseBins` clamping.
* **Later plans.** For AMR on overdensity, particles on levels, MLMG gravity and Nyx
  as a template, AMReX is a complete toolbox.
