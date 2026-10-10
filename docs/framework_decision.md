# Framework decision: AMReX vs Parthenon (Stage C spike)

*Written 2026-10-08 on branch `framework-spike`. Prototypes in `spike/`,
benchmarks in `spike/bench/results/`, Hyak job scripts in `spike/hyak/`.
Hyak CPU and GPU results (§3.2–3.3, §4) added 2026-10-10 on branch
`spike-hyak-bench`.*

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
   Parthenon ran it about 1.8× slower. On a Hyak `cpu-g2` node (§3.2) the
   two are within 13% of each other at their best layout (AMReX 4% faster
   at 256³), and AMReX beats the OpenMP-only reference by 1.34× at 256³.
   Both match the reference on the ray
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
nvcc ≥ 12.2. Both are available on Hyak (gcc 13.2, CUDA 12.4.1), where
the AMReX CUDA library builds (§4).

**Main risks of AMReX:**
- Its SoA FABs fit badly with per-cell struct kernels on the CPU.
- Keeping RT deposits consistent between coarse and fine levels is ours to
  write.
- It needs C++20 (gcc ≥ 11) and, for GPU, nvcc ≥ 12.2. Checked on Hyak:
  gcc 13.2 and CUDA 12.4.1 work.

None of these is a blocker; each has a mitigation in the plan below.

**Particle precision (reviewed 2026-10-09).** AMReX's `ParticleReal` is
one build-wide type (`AMReX_PARTICLES_PRECISION`), and SoA integer
components are 32-bit. The spike flagged "float rays but double dark
matter" as a conflict. It is not one: the original stores the dark matter
as `real(4)` (`xv_dm` in `radhydro_global.f90`), the port keeps it as
`float` (`PMState::xv`), and every PM validation, bit-for-bit against the
gfortran harness, already runs with float positions and velocities in
cell units. So the **single-precision particle build is the one that
reproduces the original**, and it halves the ray memory relative to a
double build. Quantitatively, a float position in [0, N) cells has a
resolution of N·2⁻²⁴, i.e. 7.6e-6 cells at N = 128, 3e-5 at 512 and
6e-5 at 1024; even a random walk of that error over 10⁴ steps stays
below 10⁻² cells, against a PM force resolution of one to two cells. With
AMR at 3–4 levels over a 512³–1024³ base the resolution in fine cells is
still below 10⁻³. The only cases where float positions become marginal
are effective resolutions of order 10⁵ cells per side, which are not in
the plan. The halo-mass-function gate (`STATUS.md`, Next 1) is the test
that would reveal a precision problem if one appeared.

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

### 3.2 Hyak klone, one `cpu-g2` node (AMD Genoa, 2 × 96 cores)

Built with `spike/hyak/build.sh` (gcc 13.2, OpenMPI 4.1.6, `-march=znver4`,
`-ffp-contract=off`) and checked bitwise on the node first: AMReX quick
matrix 49/49 PASS, Parthenon 52/52 PASS. Jobs 41688290 (128³) and 41688291
(256³), 2026-10-09: `spike/bench/results/hyak_node_N128_41688290.csv`,
`hyak_node_N256_41688291.csv`. The same matrix was run on 2026-10-08
(`hyak_20261008_*.csv`); RT wall times agree to a median 0.4–1.6% between
the two days, except the pure-OpenMP 1 × 192 layouts, which are unstable
(up to 30%) and in any case the worst.

RT, log-normal box, 4 steps. 128³ with 512 sources (2.84e9 ray-substeps);
256³ with 512 sources (3.83e10 ray-substeps). Wall seconds; "b" is the box
or meshblock edge. `rt_serial` is OpenMP only (one rank).

| layout (ranks × threads) | serial | AMReX b32 | AMReX b16 | Parthenon b32 | Parthenon b16 | 256³ AMReX b32 | 256³ Parthenon b32 | 256³ serial |
|---|---|---|---|---|---|---|---|---|
| 1 × 12 | 267 | 457 | 487 | 484 | 567 | | | |
| 1 × 192 | **55.6** | 255 | 395 | 304 | 809 | 3382 | 4251 | 738 |
| 2 × 96 | | 143 | 172 | 177 | 314 | 1684 | 2150 | |
| 8 × 24 (1/NUMA) | | 83.9 | 84.7 | 85.6 | 107 | 771 | 1171 | |
| 24 × 8 (1/CCD) | | 69.1 | **65.1** | **56.7** | 59.6 | **550** | **572** | |
| 48 × 4 | | 97.3 | 68.9 | 57.2 | 68.8 | 647 | 580 | |
| 96 × 2 | | 96.3 | 69.9 | — ¹ | 64.5 | 770 | 800 | |
| 192 × 1 | | 193 | 83.6 | — ¹ | 81.1 | 1073 | 977 | |

¹ 128³ has only 64 blocks of 32³, fewer than the ranks; Parthenon aborts.

Hydro, 256³ blast, 10 steps (1.68e8 cell updates): best AMReX 0.43 s (48 × 4,
b32), best Parthenon 0.33 s (24 × 8 and 48 × 4, b32), serial OpenMP 0.75 s.

Reading:
- **24 ranks × 8 threads (one rank per CCD) is the best layout** for both
  frameworks and both sizes. It is the default `LAYOUT` of the scaling job.
- At **256³, the production-like size, AMReX is 4% faster than Parthenon**
  (550 vs 572 s) and **1.34× faster than the frameworkless OpenMP code**
  (738 s), which cannot use MPI. At 128³ Parthenon is 13% faster (56.7 vs
  65.1 s) and the serial code ties with it. Either way the two frameworks
  are within the 10% band of §3.1's criterion, so the CPU numbers do not
  change the recommendation; the AMR/particle argument decides.
- Pure OpenMP across the whole node is poor for both frameworks (3.9–7.4×
  slower than 24 × 8), so production should run several ranks per node.
- Per-core efficiency still drops from 12 to 192 cores: AMReX at 1 × 12 needs
  1.9 µs per ray-substep·core, at 24 × 8 on 192 cores 4.4 µs. Memory
  bandwidth per core on Genoa is the likely limit; it is worth a profile in
  Stage D.

### 3.3 Hyak, 1, 2 and 4 nodes

Job 41688292 (`hyak_scale_41688292.csv`), 24 × 8 per node, 2026-10-09.

Strong scaling, RT 256³, b32 (3.83e10 ray-substeps):

| nodes | AMReX wall | speed-up | efficiency | Parthenon wall | speed-up | efficiency |
|---|---|---|---|---|---|---|
| 1 | 550 s | 1 | 100% | 574 s | 1 | 100% |
| 2 | 389 s | 1.42 | 71% | 345 s | 1.67 | 83% |
| 4 | 259 s | 2.13 | 53% | 268 s | 2.15 | 54% |

Weak scaling, RT at 192³ / 240³ / 304³ on 1 / 2 / 4 nodes (sources ∝
volume), b16. The ray work per node is not constant (8.9e9, 13.6e9,
20.5e9 ray-substeps per node), so the table compares the cost per unit of
work, ns per ray-substep·core:

| nodes | AMReX | efficiency | Parthenon | efficiency |
|---|---|---|---|---|
| 1 | 4669 | 100% | 5073 | 100% |
| 2 | 6006 | 78% | 7455 | 68% |
| 4 | 6190 | 75% | 6763 | 75% |

Reading:
- 256³ is small for 4 nodes: 512 boxes over 96 ranks is about 5 boxes
  per rank, and the sources are clustered, so the 4-node strong-scaling
  efficiency (53–54%) is mostly load imbalance. The ray-weighted
  `DistributionMapping` planned for Stage E (§6, item E4) targets exactly
  this. The weak-scaling efficiency of 75% at 4 nodes is the more relevant
  production number.
- **AMReX hydro across nodes is anomalous**: 6.8–7.6 µs per cell-update·core
  on 2 and 4 nodes against 0.69 on one node, while Parthenon goes from 0.63
  to 1.0–1.5. The RT, which uses the same `FillBoundary` across nodes, does
  not show it. Not investigated yet. Likely candidates are the many small
  16³ boxes per message or the ghost-cell fill pattern in
  `spike/amrex/hydro_amrex.cpp`. Look at it first in Stage E, before using
  the AMReX hydro at scale.
- The CPU results confirm the recommendation: AMReX and Parthenon are within
  ~10% on RT at the best layout and on multi-node RT; AMReX stays the
  choice for AMR with particles.

## 4. GPU

### 4.1 What runs

- `spike/common` compiles as CUDA device code. Array constants had to move
  into function-local `constexpr` arrays.
- **The Parthenon prototype runs on Hyak GPUs** (Kokkos 4.6 CUDA, CUDA
  12.4.1, gcc 13.2; A100 80GB PCIe as sm_80, H200 as sm_90; one rank, one
  GPU). `spike/hyak/build_gpu.sh` builds it on a `ckpt` node. Two fixes
  were needed:
  - `--fmad=false` for the CUDA targets (`spike/parthenon/CMakeLists.txt`).
    nvcc contracts a·b+c into FMAs in device code by default, which
    `-ffp-contract=off` does not reach.
  - Native 64-bit integer atomics (`ll_add`/`ll_fetch_add` in
    `rt_parthenon.cpp`). In this Kokkos build `Kokkos::atomic_add` on
    `long long` compiles to a compare-and-swap loop (no `ATOMG.ADD.64`
    in `cuobjdump`). On the photon histogram, where every ray adds to one
    counter, that made a 128³ substep take about 200 s. Device code now
    calls CUDA's `atomicAdd(unsigned long long*)`, which gives the same
    bits; host code is unchanged.
  The CPU builds of the final source still pass the quick validation
  bitwise (AMReX 49/49, Parthenon 52/52, job 41703917).
- **The AMReX prototype still has no CUDA target.** The AMReX 26.10 CUDA
  library itself builds on klone (nvcc 12.4). Adding the target and running
  it is Stage D item 4. **All GPU numbers below are therefore from
  Parthenon**, on the same shared kernels.

### 4.2 Validation on the GPU

Run against `rt_serial`/`hydro_serial` with `compare.py --rtol 0
--atol-frac 0` / `cmp`. A100 and H200 give identical results.
- Hydro Sod (no transcendental functions): **bitwise PASS**, mb 32 and 16.
- Hydro blast: round-off differences. Its initial Gaussian is computed with
  `exp` on the device.
- RT (`uniform_8src`, `lognormal_16src_smallbuf` with `legacy_physics` 1
  and 0, `legacy_rt=0`):
  - **all ray statistics are identical in every step** (casts, splits,
    merges, deletions, live rays);
  - fHI, fHeI, fHeII, temp, rhog, tau and the gas velocities are bitwise
    identical;
  - ngam and fgam differ by at most 1e-16 to 5e-16 in 10–450 cells per
    step, where one side is 0 or a few fixed-point quanta, and 1–2 cells of
    nu, dx, gamma and lambda differ at the float-ulp level.
  The cause is the absorbed fraction `1 - exp(-tau·dl)`
  (`ray_kernels.hpp:158, 217`). CUDA's double `exp` is accurate to 1 ulp
  while glibc's is essentially correctly rounded, and a 1-ulp difference
  can flip the 2⁻⁴⁴ fixed-point rounding of a deposit. On the 128³ bench,
  over 196 substeps, this changed 1 ray-substep in 2.8e9. Bitwise GPU = CPU
  results would need the same correctly rounded `exp` (e.g. CORE-MATH or a
  table) in both builds. That is a Stage D design decision; the integer
  bookkeeping itself is already exact on the GPU.

### 4.3 Timing: 128³ bench problem on one GPU vs one `cpu-g2` node

The same problem as the node benchmark (`bench_N128`, 512 sources, 4 steps
× 49 RT substeps = 196 substeps, 2.84e9 ray-substeps).
`spike/bench/results/hyak_gpu_N128.csv`, notes and per-phase timers in
`hyak_gpu_N128_notes.md`.

| GPU | mb (blocks) | wall, 4 steps | median RT substep | ns per ray-substep |
|---|---|---|---|---|
| A100 | 32 (64) | 56.8 s | 283 ms | 20.0 |
| A100 | 64 (8) | 24.1 s | 122 ms | 8.5 |
| A100 | 128 (1) | 16.9 s | 94 ms | 6.0 |
| H200 | 32 (64) | 45.4 s | 234 ms | 16.0 |
| H200 | 64 (8) | 16.0 s | 86 ms | 5.6 |
| H200 | 128 (1) | **11.0 s** | **62 ms** | 3.9 |

Small blocks are slow because each block is a separate set of kernel
launches. The largest phase on the H200 at mb 128 is `step1` (3.9 of
11.0 s): one thread per cell, with a serial sort in source cells. Then sort
1.7 s, halo 1.2 s, cast_1 1.2 s, migrate 0.9 s.

Against one `cpu-g2` node (2026-10-09, §3.2, same problem):

| one GPU (mb 128) vs one node | best CPU code: `rt_serial` 1 × 192, 55.6 s | Parthenon 24 × 8, 56.7 s | AMReX 24 × 8, 65.1 s |
|---|---|---|---|
| H200, 11.0 s | **5.0×** | **5.1×** | **5.9×** |
| A100, 16.9 s | 3.3× | 3.4× | 3.8× |

### 4.4 Is one GPU worth ≥ 5 cpu-g2 nodes? (PLAN.md, NSF ACCESS rule)

**H200: just, at 5.0–5.9×. A100: no, at 3.3–3.8×.** This holds on the RT
alone, which dominates the cost.

The ratio understates what a GPU code could do:
- the GPU code is untuned: one thread per cell in `step1`, per-block
  kernel launches, host synchronisation every substep, AoS ray layout;
- 128³ is small for a GPU.

It also flatters the GPU in two ways:
- at mb 128 the whole box is one block, so the GPU pays no exchange
  between blocks, while the CPU runs do; a multi-GPU run would pay
  exchange between GPUs;
- the comparison covers the RT only, not hydro or gravity.

The ratio is measured with Parthenon, not AMReX, the chosen framework.

**Reading:** the ACCESS threshold is met on H200-class GPUs, without
tuning, on the dominant kernel. An ACCESS proposal is justified if it
targets H100/H200-class nodes, not A100s. Before writing the
performance section of a proposal, confirm the ratio with the AMReX CUDA
build and a team-level `step1` (Stage D item 4), at 256³ and on more than
one GPU.

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
| GPU | CUDA/HIP/SYCL native; needs C++20, nvcc ≥ 12.2 | Kokkos (CUDA/HIP/SYCL/OpenMP); CUDA build runs on A100/H200 (§4); per-block launches would need packs |
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
   tests must stay green. **Partly done by PR #15:** the hydro kernel and
   `rr_hd.hpp` already live in `include/radray/kernels/`, extended with the
   dual-energy formulation, and are maintained there; `spike/common/`
   keeps the frozen copies the prototypes were validated against. The ray
   and RT-cell kernels still have to move.
2. Turn `rt_amrex` into a library component (`src/amrex/rt_level.cpp`).
   It should work on one `AmrLevel`-like object: a ray container per level,
   cell `FabArray`s, deposits as int64 with `SumBoundary`, 2 ghost layers
   with an assertion on the ray reach.
3. Choose the ray encoding from the compression study (Stage B). Build
   AMReX with `AMReX_PARTICLES_PRECISION=SINGLE`: rays as SoA float
   components, and dark matter as float positions and velocities in a
   second container, the precision of the original (see "Particle
   precision" above). Use 2 ghost layers for the deposit accumulators
   (validated in the spike; `--ghosts 2`), with the assertion on the ray
   reach kept.
4. GPU: add the CUDA target, move the rate tables to device memory,
   implement a team-level step-1 for source cells, then benchmark on
   `ckpt` A100/H200 against a node (PLAN.md: ≥5× ⇒ ACCESS proposal).
5. Validation: bitwise against `rt_serial` (the spike's test matrix
   becomes CTest), and round-off against `rt_driver --set legacy_rt=0`, as
   in §2.

**Stage E: full uniform-grid code on AMReX with MPI (production platform).**
1. Hydro: the library's unsplit MUSCL-Hancock kernel
   (`include/radray/kernels/hydro_mh.hpp`, `hydro_scheme = unsplit`,
   PR #15, `docs/hydro_unsplit.md`), already validated with the hydro
   test suite and coupled to RT and PM on the uniform grid. On AMReX
   patches it needs only the `FillBoundary`/flux pattern of
   `spike/amrex/hydro_amrex.cpp`. The legacy split scheme stays in the
   serial code.
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

- Approve AMReX as the framework. The Hyak numbers (§3.2–3.3) are in and
  do not change the recommendation.
- ACCESS: one H200 is worth 5.0–5.9 `cpu-g2` nodes on the RT, an A100
  3.3–3.8 (§4.4). Decide whether to prepare a proposal for H100/H200-class
  GPUs now or after the AMReX CUDA port confirms the ratio.
- Moving `spike/common` into the main library (Stage D item 1) touches
  `src/rt_tools.cpp` and `src/ray_tracing.cpp`. PR #10 (production
  readiness, OpenMP) is merged, so this can start.

## 8. Review (2026-10-09, after merging `main` with PR #10)

- The decisive claims were re-checked in the sources: Parthenon's
  `particles.rst` states that swarms do not support AMR, and
  `src/mesh/mesh-amr_loadbalance.cpp` and `mesh_refinement.cpp` contain no
  swarm handling; Parthenon's driver sets one global time step
  (`SetGlobalTimeStep`). AMReX's `ParticleReal` is build-wide and its SoA
  integer components are 32-bit.
- `main` (production readiness, PR #10) was merged into the branch. With
  the merged library: the fast suite passes (80/80); `rt_serial` is
  bitwise unchanged and still matches the new `rt_driver --set
  legacy_rt=0` with identical ray statistics on all five problems; the
  AMReX quick matrix passes 49/49; Parthenon passes at 1, 4 and 8 ranks.
- The particle-precision risk is withdrawn (see "Particle precision"
  above): the original and the port already use float dark matter.
- The recommendation stands: AMReX.
