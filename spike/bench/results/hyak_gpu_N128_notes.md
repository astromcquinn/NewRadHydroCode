# GPU RT timing, 128^3 bench problem (2026-10-10)

Code: `rt_parthenon` (Parthenon v25.12, bundled Kokkos 4.6, CUDA backend), one
rank, one GPU, built with CUDA 12.4.1, gcc 13.2.0, `--fmad=false`, from the
`spike-hyak-bench` worktree. Problem: `$SPIKE_WORK/problems/bench_N128` (512
sources, 4 steps x 49 RT substeps = 196 substeps), the same as the CPU node bench.

Valid results: `gpu_N128_summary.csv` (= gpu_41703754/41703755 for mb 32/64,
gpu_41703836/41703837 for mb 128). Earlier gpu_*.csv are from the slow
pre-fix build (CAS atomics) or killed runs; ignore them.

| GPU | mb (blocks) | wall 4 steps | median substep | ns/ray-substep |
|---|---|---|---|---|
| A100 80GB PCIe (sm_80) | 32 (64) | 56.8 s | 0.283 s | 20.0 |
| A100 | 64 (8) | 24.1 s | 0.122 s | 8.5 |
| A100 | 128 (1) | 16.9 s | 0.094 s | 6.0 |
| H200 (sm_90 build) | 32 (64) | 45.4 s | 0.234 s | 16.0 |
| H200 | 64 (8) | 16.0 s | 0.086 s | 5.6 |
| H200 | 128 (1) | 11.0 s | 0.062 s | 3.9 |

work = 2 835 088 334 ray-substeps (CPU codes: ...335, see validation).

Phase timers (s, whole run), H200 mb 128: cell 0.10, halo 1.22, hist 0.53,
step1 3.87, cast1 1.22, fold 0.10, phys 0.52, cast3 0.85, sort 1.66,
migrate 0.91, post 0.03. A100 mb 128: cell 0.25, halo 1.61, hist 0.67,
step1 5.59, cast1 2.06, fold 0.16, phys 0.87, cast3 1.29, sort 2.85,
migrate 1.53. step1 (per-cell classify/emit/split/merge, one thread per
cell, serial heap sort in source cells) is the largest phase.

Same problem on one cpu-g2 node (node_41630446.csv, 2026-10-08): best
serial 1:192 53.4 s, Parthenon 24:8 mb32 57.7 s, AMReX 24:8 mb16 64.8 s.
=> H200 is 4.8x (vs best CPU code) / 5.2x (Parthenon vs Parthenon);
   A100 is 3.2x / 3.4x.

## Fixes needed to get here
1. `--fmad=false` for the CUDA targets (spike/parthenon/CMakeLists.txt):
   nvcc contracts FMAs in device code by default.
2. 64-bit integer atomics: Kokkos::atomic_add on long long (and on
   unsigned long long) compiled to CAS loops in this build (cuobjdump: no
   RED/ATOMG .ADD.64). Under contention (the photon histogram: every ray
   adds to one counter) a 128^3 substep took ~200 s, the 4-step run >1.4 h.
   Device code now calls native atomicAdd(unsigned long long*)
   (`ll_add`/`ll_fetch_add` in rt_parthenon.cpp): 5042 s for step 1 alone ->
   16.9 s for all 4 steps (A100, mb 128).

## Validation on GPU (A100 and H200 identical)
- hydro sod (no transcendentals): bitwise PASS (cmp), mb 32 and 16.
- hydro blast: FAIL at round-off; its IC is a Gaussian (std::exp,
  hydro_mh.hpp:200) evaluated on the device.
- RT uniform_8src, lognormal_16src_smallbuf (legacy_physics 1 and 0),
  legacy_rt=0, vs rt_serial, compare.py --rtol 0 --atol-frac 0: FAIL, but
  all 10 integer/stat lines are bitwise identical in all 10 steps (ray
  counts, casts, splits, merges, deletions, ngam_rays_total), and fHI,
  fHeI, fHeII, temp, rhog, tau, u/w_gas are bitwise identical. Differences:
  ngam (max abs 9.5e-17) and fgam (4.6e-16) in 10-450 cells per step where
  one side is 0 or a few fixed-point quanta, plus 1-2 cells of nu (1.9e-6),
  dx (3e-8), gamma (2.6e-23), lambda (1.2e-32) at float-ulp level.
  Cause: the absorbed fraction 1 - exp(-tau*dl) (ray_kernels.hpp:158, 217)
  in double; CUDA's exp is within 1 ulp, glibc's is ~correctly rounded, and
  a 1-ulp difference flips the 2^-44 fixed-point rounding of a deposit.
  At 128^3 over 196 substeps this changed one ray-substep in 2.8e9.
  Bitwise GPU = CPU would need a shared correctly-rounded exp (e.g. CORE-MATH
  or a table) in both builds.
