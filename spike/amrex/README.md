# AMReX prototype of the spike (hydro_amrex, rt_amrex)

This is the AMReX half of the framework spike. It runs the framework-independent kernels in
`spike/common` (unchanged) on AMReX data structures, distributed over MPI
ranks and OpenMP threads. It must reproduce `spike/serial` **bit for bit**.

| executable | reproduces | AMReX pieces |
|---|---|---|
| `hydro_amrex` | `spike/serial/hydro_serial` (unsplit MUSCL-Hancock + HLLC, `hydro_mh.hpp`) | `Geometry` (periodic), `MultiFab` U (5 comps, 2 ghosts), Q (20 comps, 1 ghost), 3 face-centred flux `MultiFab`s; `FillBoundary`; `MFIter` tiling + OpenMP; `ParallelFor`; `ReduceOps<ReduceOpMax>` + `ParallelDescriptor::ReduceRealMax`; `ParallelCopy` to rank 0 for output |
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
  26.10 with the options above via `FetchContent`. That path was not tested here
  because the sandbox has no GitHub access.
* Our targets get `-O2 -ffp-contract=off -fno-fast-math` (same as `spike/serial`).
  They also link `AMReX::Flags_INLINE` (`-finline-limit=43210`); see gotchas.
* `-DRT_CELL_SOA` (add to `CMAKE_CXX_FLAGS`) switches the RT cell storage to a
  45-component float SoA `FabArray`, which is the GPU-friendly layout.
* Building the prototype takes about 40 s for `rt_amrex` and 25 s for `hydro_amrex` at
  -j2. That is mostly AMReX headers: `rt_amrex.cpp` pulls in ~400 k lines after
  preprocessing.

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

VALIDATION_RESULTS

## Timings

TIMINGS

## Notes for the framework assessment

NOTES
