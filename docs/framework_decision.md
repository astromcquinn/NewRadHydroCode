# Framework decision: AMReX vs Parthenon (Stage C spike)

*Written 2026-10-08 on branch `framework-spike`. Prototypes in `spike/`,
benchmarks in `spike/bench/results/`, Hyak job scripts in `spike/hyak/`.*

## Recommendation

**Adopt AMReX (26.10) for Stages D–G.** Keep the spike's
framework-independent kernel layer (`spike/common`) as the place where all
physics lives, so the framework only owns the mesh, the particles and the
communication.

The two prototypes are both correct, bitwise, so the decision comes down to
the later stages and to CPU cost:

1. **AMR with particles.** Parthenon's swarms do not support AMR. Its own
   docs say "AMR is currently not supported" (`doc/sphinx/src/particles.rst`),
   and the remesh/load-balance code (`src/mesh/mesh-amr_loadbalance.cpp`)
   never touches swarms; only static refinement works with swarms. Both our
   rays and our dark matter are particles, so this blocks Stages F and G.
   AMReX has `AmrParticleContainer`: particles on AMR levels, redistributed
   to the finest level that covers them, with ghost/virtual particles for
   gravity.
2. **Sub-cycling.** RT and hydro on fine levels need finer time steps.
   AMReX supports sub-cycling in time (`Amr`/`AmrLevel`, or by hand on
   `AmrCore`, documented in `AmrCore.rst` §"subcycling"). Parthenon advances
   all blocks with one global time step.
3. **Gravity.** AMReX has MLMG multigrid Poisson on AMR, an FFT module
   (`Src/FFT`, including open-boundary Poisson) for the uniform-grid PM, and
   particle–mesh deposit helpers. Nyx, cosmological hydro plus N-body on
   AMR, is built on exactly these. Parthenon has a geometric multigrid and
   MG-preconditioned CG/BiCGSTAB, which are good, but no FFT and no
   particle–mesh layer.
4. **CPU cost**, the production platform. On this 4-core VM, AMReX ran the
   RT about 1.3× slower than the frameworkless reference on one core.
   Parthenon ran it about 1.8× slower. Both match the reference on the ray
   kernels themselves; the difference is mesh bookkeeping (§3).
5. **Ease of use for a newcomer.** AMReX has the largest documentation set
   (a sphinx manual of about 30 chapters, many tutorials) and a large user
   community: Nyx, Castro, MAESTROeX, WarpX, PeleC/PeleLM, AMR-Wind,
   Quokka. It also has a stock reverse ghost exchange (`SumBoundary`) for
   particle deposits. The Parthenon prototype needed more workarounds: a
   hand-written reverse exchange, its own particle sort, argv filtering, and
   a non-default loop layout to get OpenMP threading at all.

What Parthenon would have given us, and what we give up:
- Kokkos, with one source tree for CUDA, HIP, SYCL and OpenMP. A
  compile-only CUDA build of the Parthenon prototype works today.
- A lighter octree-of-meshblocks AMR, a natural fit for GPUs.

AMReX is GPU-native too, through CUDA, HIP and SYCL backends. That is
enough for our secondary GPU target. Its constraints are C++20 and
nvcc ≥ 12.2, which is why it could not be compiled for GPU here (§4).

**Main risks of AMReX:**
- Particle precision is one build-wide switch: we want float rays and
  double DM.
- Its SoA FABs fit badly with per-cell struct kernels on the CPU.
- Keeping RT deposits consistent between coarse and fine levels is ours to
  write.

None of these is a blocker; each has a mitigation in the plan below.

## 1. What was built

All three implementations share `spike/common` (about 1500 lines,
header-only, host/device). It holds line-by-line copies of the serial port's
kernels:
- the ray kernels: cast, advance, split, merge (as a reduction), source
  emission, classification and redshifting;
- the per-cell RT physics (`radiative_physics` for N_rt_gas = 1);
- a new unsplit MUSCL-Hancock + HLLC hydro kernel.

The RT substep follows the sorted-ray-array design
(`docs/design_parallel_and_encoding.md` §2):
- rays sorted by home cell;
- a two-pass count/scan/write per cell for emission, split and merge;
- canonical merge order (cone, payload hash, payload bytes);
- photon deposits into **int64 fixed-point accumulators**, which are exact
  and associative;
- integer histograms for the `ngam_ray` thresholds;
- migration of rays to their owning rank after every substep.

| code | what | lines (no blanks/comments) |
|---|---|---|
| `spike/serial/rt_serial.cpp` | reference, OpenMP, no framework (the "in-house" option) | 269 |
| `spike/amrex/rt_amrex.cpp` | `FabArray` cells, `ParticleContainerPureSoA` rays, `FillBoundary`, `SumBoundary`, `Redistribute`, `DenseBins` | 635 |
| `spike/parthenon/rt_parthenon.cpp` (+ `parth_common.hpp`) | Parthenon variables, Swarm rays, task-list ghost exchange and swarm comm, hand-written reverse exchange | 854 + 129 |
| `spike/*/hydro_*.cpp` | MUSCL-Hancock on `MultiFab` / meshblock variables | 110 / 208 / 207 |

Kernel (a) is a new unsplit MUSCL-Hancock + HLLC scheme, not the original's
split TVD scheme. The original smooths the moving-frame velocity with an FFT
along whole pencils and sub-cycles each pencil, which does not decompose
into blocks. ROADMAP §5 already plans the move to an unsplit Godunov scheme
for AMR.

## 2. Validation

**Serial sorted-array reference vs the serial port**
(`rt_serial` vs `rt_driver --set legacy_rt=0`), on all five RT problems:
- **Ray statistics are identical every step:** cast, split, merge, delete,
  live rays. That includes the merge path of `lognormal_16src_smallbuf`
  (1477 merges in step 10) and a 64³ benchmark with 400k rays.
- **Fields agree to round-off**, with differences only from the summation
  order of deposits and merges. The exceptions are explained:
  - `dx` in source cells: it is the normalised sum of 192 symmetric
    directions. The float sum leaves noise of order 0.1 after
    normalisation; the integer sum cancels exactly and gives 0.
  - `fgam` in a few I-front cells: up to 1.5e-3 relative where
    fgam ≈ 1e-4 (absolute ≤ 1e-5). fgam is a remaining fraction, and
    round-off is amplified there.
  - A build with 256× coarser deposits (`rt_serial_fx36`) moves the same
    cells by proportionally more, which confirms these are summation-order
    effects.
- `rt_serial` is **bitwise identical across OpenMP thread counts**.

**Prototypes vs the serial reference: bitwise (`compare.py --rtol 0
--atol-frac 0`, `cmp` for hydro)**, re-run after merging onto the final
`spike/common`:

| prototype | RT cases | hydro cases |
|---|---|---|
| AMReX | uniform_8src, lognormal_16src_smallbuf, lognormal with legacy_physics=0 × 1/2/4/8 ranks × box size 16/8 (and 32) × 1/2 threads: 75/75 PASS in the agent's full run, 49/49 in the quick re-run after the merge | blast and sod × 1/2/4 ranks × 32/16 × 1/2 threads, all PASS |
| Parthenon | the same problems × 1/2/4/8 ranks × 16³/8³ blocks × 1/2 threads: 46/46 PASS, plus 12/12 in the re-run after the merge | 38/38 PASS |

So the sorted-array design with integer deposits gives
**decomposition-independent, thread-count-independent results in both
frameworks**. That is the property the design doc (§2.3, §3.4) asked for.
Neither framework helped or hindered: neither offers reproducible float
reductions, and both handle integer data exactly.

**Hydro sanity check:** Sod in a periodic box (`spike/bench/check_sod.py`)
gives an L1 density error of 2.8e-2 / 1.5e-2 / 7.8e-3 at n = 32 / 64 / 128.
That is first-order convergence, as expected with discontinuities. Mass is
conserved to round-off.

## 3. Benchmarks

### 3.1 Here (4-core Intel Xeon VM, 2.1 GHz, OpenMPI, gcc 13, -O2)

RT problem: 64³ log-normal box, 64 sources, 3 steps with 24 substeps each,
400k live rays at the end and 1.36e7 ray-substeps in total
(`spike/bench/results/local_4core_N64.csv`). Best of 2 runs. "Block" is the
box or meshblock edge length.

| RT | serial | AMReX 32³ | AMReX 16³ | Parthenon 32³ | Parthenon 16³ |
|---|---|---|---|---|---|
| 1 rank × 1 thread | **19.6 s** | 26.6 s | 29.6 s | 35.7 s | 48.3 s |
| 1 rank × 4 threads | 7.5 s | 11.6 s | 13.7 s | 15.8 s | 21.2 s |
| 2 ranks × 2 threads | — | 11.3 s | 13.8 s | 15.2 s | 25.0 s |
| 4 ranks × 1 thread | — | 13.3 s | 14.5 s | 16.8 s | 19.6 s |

Hydro: 64³ blast, 10 steps.

| hydro | serial | AMReX 32³ | AMReX 16³ | Parthenon 32³ | Parthenon 16³ |
|---|---|---|---|---|---|
| 1 rank × 1 thread | 0.75 s | 0.81 s | 0.92 s | 1.04 s | 1.35 s |
| 1 rank × 4 threads | 0.21 s | 0.34 s | 0.38 s | 0.33 s | 0.47 s |
| 4 ranks × 1 thread | — | 0.35 s | 0.36 s | 0.31 s | 0.48 s |

Reading:
- The ray kernels cost the same in all three codes: cast_3 takes
  0.46–0.51 s in each on the merge test problem. The spike adds nothing
  per ray.
- The framework overhead is on the mesh:
  - **AMReX:** mostly `SumBoundary` of the 13 int64 deposit components
    over 4 ghost layers, plus zeroing the ghost-grown accumulators. It
    drops with 2 ghost layers (about 20% faster), which suffice for the
    measured ray reach.
  - **Parthenon:** mostly its ghost exchange of tau and fgam. The unpack
    does a per-element index transform, and callgrind puts `SetBounds` at
    28% of all instructions. Then the hand-written deposit fold, and
    per-cell phases 2.5–7× slower than serial on double, component-slowest
    variables.
- Both overheads scale with surface/volume, so they shrink with bigger
  blocks: 32³ is better than 16³ everywhere. Production blocks would be
  32³–64³.
- The 4-core VM scales poorly even for the serial code (2.6× on 4
  threads), so the parallel numbers here say little about Genoa. **The
  Hyak runs below are the ones that count.**

### 3.2 Hyak (to run)

`spike/hyak/README.md` has the commands. `spike/hyak/bench_node.sbatch`
measures one `cpu-g2` node:
- 128³ (512 sources) and 256³ RT boxes, and a 256³ hydro blast;
- layouts 1:192, 2:96, 8:24 (one rank per NUMA domain), 24:8 (one per CCD),
  48:4, 96:2 and 192:1, plus a thread-scaling curve.

`spike/hyak/bench_scale.sbatch` runs 1, 2 and 4 nodes:
- strong scaling at 256³;
- weak scaling at 192³, 240³ and 304³.

Copy the CSVs into `spike/bench/results/` and update this section. The
build script uses `-march=znver4`, and the projects keep
`-ffp-contract=off`. Step 3 of the README re-runs the bitwise validation on
the node.

What would change the recommendation: if Parthenon came within about 10%
of AMReX at the best Genoa layout with 32³–64³ blocks, the AMR/particle
argument still decides. The CPU numbers only confirm or refute the
secondary point.

## 4. GPU

- `spike/common` compiles as CUDA device code (nvcc 12.0, sm_80,
  `spike/gpu/kernels_compile_test.cu`). Array constants had to move into
  function-local `constexpr` arrays.
- The Parthenon prototype builds with Kokkos CUDA (nvcc 12.0, sm_80),
  compile only.
- AMReX 26.10 refuses nvcc < 12.2 at configure time, and only 12.0 is
  installable here. The AMReX prototype's loops are all device lambdas;
  its CMake has no CUDA target yet.
- Nothing has been run on a GPU. Known blockers before a GPU run:
  - `PhysCtx` holds host pointers to the rate tables; they must be copied
    to device memory;
  - host-side setup and output loops;
  - per-cell `step1` serialises in source cells (192 rays each), so it
    needs a team-level kernel there;
  - AoS layouts are uncoalesced.
- `spike/hyak/build_gpu.sh` builds AMReX with CUDA and the Parthenon
  prototype with Kokkos CUDA on a `ckpt` node, needing a CUDA module
  ≥ 12.2. This is the Stage D GPU benchmark from PLAN.md (≥5× of a node
  ⇒ apply to ACCESS).

## 5. Assessment by criterion

| criterion | AMReX | Parthenon |
|---|---|---|
| AMR on overdensity | `AmrCore` + `ErrorList`/`TagBox` tagging (e.g. a mass threshold per cell), `FillPatch`, `FluxRegister` for conservative coarse/fine fluxes; Nyx and Castro refine on density | mature octree AMR, refinement tags per block (`AmrTag`), prolongation/restriction per variable, very good for GPUs |
| Star-formation sink | particles can be created anywhere on any level (`AmrParticleContainer`); Nyx/Castro have star/sink-particle precedents **[from memory, not checked in source]** | swarm creation on a block is easy, but swarms do not survive AMR remeshing |
| RT sub-cycling on fine levels | sub-cycling in time supported; we would run our RT substeps per level ourselves anyway. Coarse/fine deposit and ray hand-off between levels is ours to write (no stock "particle deposit FluxRegister" for rays) | single global dt; sub-cycling would be ours to write |
| Particles (DM, rays) | `ParticleContainer`/`AmrParticleContainer`, `Redistribute` (local and global), `DenseBins`, neighbour particles, particle–mesh deposit/interpolation; one build-wide particle precision | Swarm: pool + mask + defrag, one-hop neighbour communication in a blocking iterative task region, **no AMR support**, `SortParticlesByCell` unsuitable (includes inactive slots, serial `std::sort`) |
| Multigrid gravity | MLMG (`MLPoisson`, `MLABecLaplacian`) on AMR; FFT module for the uniform-grid PM; open-BC Poisson | GMG and MG-preconditioned CG/BiCGSTAB on the AMR hierarchy (`src/solvers`, `example/poisson_gmg`); no FFT |
| Determinism | integer collectives exact; float `ReduceRealSum`/`MultiFab::sum` are plain MPI reductions (fixed decomposition only); GPU float ADD comm deterministic on request (`deterministic=true`, `AMReX_FBI.H`). Our integer design is bitwise on any layout | integer collectives exact; no reproducible float reductions; solvers use float `MPI_Iallreduce`. Same conclusion: our integer design is bitwise on any layout |
| Reverse ghost exchange (deposits) | `SumBoundary` / `ParallelAdd`, stock | **missing**, about 100 lines written by hand |
| CPU threading | MFIter tiling + OpenMP (`ParallelFor` alone is serial on CPU; `ParallelForOMP` exists but is not in the docs) | Kokkos OpenMP, but the default `SIMDFOR_LOOP` layout is not threaded; needs `PAR_LOOP_LAYOUT=MANUAL1D_LOOP`; design assumes one rank per core |
| Data types | `FabArray<BaseFab<T>>` for float, int64 and structs; `ParallelCopy`/`FillBoundary` need `operator+=` so not on struct FABs | variables are `Real` (double unless the whole build is single); component-slowest |
| GPU | CUDA/HIP/SYCL native; needs C++20, nvcc ≥ 12.2 | Kokkos (CUDA/HIP/SYCL/OpenMP); compile-only CUDA build done; per-block launches would need packs |
| Build complexity | CMake, library build 3 min (-j4), prototype 21 s; install + `find_package` | CMake `add_subdirectory` of Parthenon + Kokkos, 5–7 min (-j2/-j4), 4 min with our code at -j4; rebuilt per application |
| Documentation | very large sphinx manual, guided tutorials, Doxygen; plus the gotchas below | good sphinx manual (about 30 pages) and examples; particle and boundary-comm docs thin; several traps found only by reading source |
| Community / heritage | ECP/DOE-funded, large: Nyx (cosmology), Castro, WarpX, Pele, AMR-Wind, Quokka (RHD) **[project list from memory]** | LANL-led: AthenaPK, KHARMA, Phoebus, Artemis **[from memory]**; smaller |
| Code we wrote | 635 lines RT, 208 hydro | 854 + 129 lines RT, 207 hydro |

Gotchas found (details in `spike/amrex/README.md` and
`spike/parthenon/README.md`):

AMReX:
- `amrex::CellData` clashes with ours under `using namespace amrex`.
- Local `Redistribute` silently leaves a particle in the wrong box if we
  wrap its position ourselves; store the unwrapped image and let AMReX
  apply the periodic shift.
- Consumers must link `AMReX::Flags_INLINE`, or the header loops run 2×
  slower.
- `DenseBins`' box overload silently clamps out-of-range cells.

Parthenon:
- Any argv token with `=` but no `/` aborts.
- `ResetCommunication` must run in an earlier task region than `Send`.
- Swarm `Send` sorts every particle of the block.

## 6. Migration plan for Stages D–G (AMReX)

Ground rules carry over: `reference/` and the serial port stay the bitwise
legacy reference, and every new behaviour is behind a flag.

**Stage D: ray data structure in AMReX on a uniform grid, CPU first, then GPU.**
1. Promote `spike/common` to `include/radray/kernels/` and make the serial
   port call the same per-cell physics and per-ray functions, so each piece
   of physics has a single implementation. The serial port's bit-for-bit
   tests must stay green.
2. Turn `rt_amrex` into a library component (`src/amrex/rt_level.cpp`).
   It should work on one `AmrLevel`-like object: a ray container per level,
   cell `FabArray`s, deposits as int64 with `SumBoundary`, 2 ghost layers
   with an assertion on the ray reach.
3. Choose the ray encoding from the compression study (Stage B). Store rays
   as SoA float components, so AMReX particle precision is `SINGLE` for
   rays. Dark matter goes in a second container in double precision:
   either build AMReX with double particles and store ray payloads in
   float runtime components, or keep DM positions as double runtime
   components. Decide this at the start of D.
4. GPU: add the CUDA target, move the rate tables to device memory,
   implement a team-level step-1 for source cells, then benchmark on
   `ckpt` A100/H200 against a node (PLAN.md: ≥5× ⇒ ACCESS proposal).
5. Validation: bitwise against `rt_serial` (the spike's test matrix
   becomes CTest), and round-off against `rt_driver --set legacy_rt=0`, as
   in §2.

**Stage E: full uniform-grid code on AMReX with MPI (production platform).**
1. Hydro: the MUSCL-Hancock kernel, validated with the hydro test suite
   (Sod, Sedov, Zel'dovich pancake, sound waves). It is a new scheme behind
   a flag; the legacy split scheme stays in the serial code.
2. PM gravity: DM particles in an AMReX `ParticleContainer`, CIC deposit,
   and the Poisson solve with the `amrex::FFT` module (or MLMG). Validate
   against `pm_driver` to round-off on the PM problems.
3. Halo finder and sources: gather DM on the existing FoF (or a
   distributed FoF later), with sources as per-cell `Luv`/`sfr`.
4. Load balancing: AMReX `DistributionMapping` with weights
   `a·cells + b·rays` (design doc §3.1), repartitioned every hydro step.
5. Restarts (required for `ckpt`): AMReX plotfile/checkpoint for fields
   and particles.
6. Scaling on the 5 owned nodes; choose the layout from
   `bench_node.sbatch`.

**Stage F: AMR for gas, then RT.**
1. `AmrCore` with tagging on overdensity (a Lagrangian mass threshold, e.g.
   8× the base mass per cell). Use `FillPatch` for ghosts and
   `FluxRegister` for conservative refluxing, with sub-cycling in time.
2. The star-formation sink behind a flag: a fixed efficiency per
   free-fall time above a density threshold, creating star particles. The
   original's no-cooling-above-200× switch stays available.
3. RT on AMR:
   - rays live on the level of their home cell and move one local cell per
     substep, with fine levels sub-cycled within the coarse RT substep;
   - deposits stay per level, with ghost sums and a coarse/fine rule for
     deposits that cross a level boundary;
   - rays are handed off between levels when they cross a coarse/fine
     boundary (split on entry to a finer level, merge on exit).
   This is new research. Validate on the Iliev tests with refinement and
   against uniform-grid runs at the matched finest resolution.

**Stage G: gravity on AMR.** MLMG Poisson on the AMR hierarchy, with
`AmrParticleContainer` deposit for DM and virtual particles across levels,
following Nyx. Check the force accuracy against the uniform-grid PM and the
halo mass function. TreePM stays an option.

**Fallback.** If a blocker appears in Stage D that the spike did not
show, the in-house option is the spike's own `rt_serial` plus a thin
layer. That layer would need our own MPI (bricks, a 2-cell halo, ray
migration) and Kokkos for GPU, and would be a few months of
infrastructure (ROADMAP §1). Parthenon becomes attractive again if its
swarms gain AMR support. Watch its CHANGELOG.

## 7. Open points for the user

- Approve AMReX as the framework, or wait for the Hyak numbers first.
- Precision layout for particles (Stage D item 3).
- Moving `spike/common` into the main library (Stage D item 1) touches
  `src/rt_tools.cpp` and `src/ray_tracing.cpp`, which the production-
  readiness thread (OpenMP) may also be editing. Do it after that PR is
  merged.
