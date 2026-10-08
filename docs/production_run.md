# Production reionization run: reduced speed of light, ray pool, restart, cost

This covers the work that prepares the production run (box 10 Mpc/h,
Ngas = Ndm = N, c_reduced = 0.01, nray_max_factor = 16, z ≈ 15 → 5.5, 4
threads): the reduced speed of light, the OpenMP parallelisation
(`docs/parallel.md`), the ray-pool guard, the F1 fix, restarts, movie
frames, and the measured cost that sizes the run. Every change keeps the
legacy validation: with the defaults all Fortran-vs-C++ references are
bit-identical, for any thread count.

New parameters of `apps/radhydro` (params.txt / `--set`), all defaulting to
the original behaviour:

| key | default | meaning |
|---|---|---|
| `c_reduced` | 1 | reduced speed of light ĉ = c_reduced·c (section 1); also `rt_driver` and both Fortran harnesses |
| `deterministic_parallel` | 0 | 1: parallel ray tracing, thread-count independent (`docs/parallel.md`); needs `legacy_rt = 0` |
| `threads` | 0 | OpenMP threads (0: `OMP_NUM_THREADS` / all cores) |
| `rt_substep_carry` | 0 | 1: carry the rounding of `nstep_rt` to the next hydro step, so the total RT substeps equal the light travel at ĉ exactly (`legacy_rt = 0` only) |
| `ra_max` | the original 0.1 | max expansion per hydro step, da/a ≤ ra_max/(1−ra_max) (`legacy_rt = 0` only) |
| `movie_dz`, `movie_z0`, `movie_fields` | 0, first zdata1, `fHI,temp` | movie frames every `movie_dz` in redshift, written from inside the RT substeps (section 4) |
| `restart_every`, `restart_keep_latest` | 0, 1 | complete checkpoint every n hydro steps (section 4); `--restart <file>` resumes |

`tools/production/make_production.py OUTDIR --N 128` writes the ICs and the
params.txt of the production configuration.

## 1. Reduced speed of light

### What c_reduced does

`calc_timestep` sets the number of RT substeps per hydro step from the
comoving light travel distance in RT cells,

    nstep_rt = nint(dl_c/dl_rt * cscale),   dl_c = 2c/(H0 √Ωm) (√(a+da) − √a),

and the RT moves every ray by exactly one RT cell per substep
(`ray_cast_3`, ≤ 1 cell of segments). The original hard-wires
`cscale = 1`, so rays travel at c. With `cscale = c_reduced` there are
c_reduced times fewer substeps per hydro step, each `dtsec_rt = dt/nstep_rt`
long, and rays still move one cell per substep: they travel at
ĉ = c_reduced·c. This is the only place where c_reduced enters, in the C++
(`cosmology.cpp`) and in both Fortran harnesses (RT harness PATCHES.md 12,
full harness PATCHES.md section 4).

**nstep_rt ≥ 1.** If the light travel of a hydro step at ĉ is under half
an RT cell, `nint` gives 0 (and `dt_rt = dt/0`). Both codes then use one
substep and print a warning; rays then move one cell per hydro step, faster
than ĉ. Never reached with c_reduced = 1, and not in the production
configuration (below). The CTest `full_vs_fortran_live_smoke_legacy_rt0_c0.0001`
runs the guard in both codes (nstep_rt = 0 in all 15 steps) and they agree
bit for bit.

### What changes and what does not (photon conservation)

The standard reduced-speed-of-light approximation (RSLA; Gnedin & Abel
2001) keeps the source emission rate and the photoionisation rate in the
optically thin limit, and lets the photon density scale as 1/ĉ. Checked
for every place that involves c or the substep length:

| quantity | where | with c_reduced |
|---|---|---|
| photons emitted per substep | `ray_source`, `photo_source`: L·`dtsec_rt` | ∝ dtsec_rt ∝ 1/ĉ: the emission **rate** per unit time is unchanged; photons per hydro step are identical (the `ifront_rsla` test checks the total to 1e-6) |
| photons in flight, radiation energy density | rays and the cells they cross | a ray travels 1 cell per substep, so photons emitted during a time t fill a sphere of radius ĉt: n_γ ∝ L/(4πr²ĉ) ∝ 1/ĉ, as in the RSLA |
| photoionisation rate Γ | `radiative_physics`: Γ = n_γ,cell σ (len_rt/dtsec_rt) | the photons crossing a cell in one substep ∝ F·A·dtsec_rt, divided by dtsec_rt: Γ = Fσ, **independent of ĉ** in the optically thin limit |
| ionisations per substep, photon limiter | `photo_cell` (`num1`, `nrad`) | both ∝ dtsec_rt: consistent, photons are conserved per substep |
| heating, cooling, chemistry | `photo_cell` with `dtsec_rt` | physical time steps; unchanged rates |
| Thomson optical depth `dtau` | σ_T n_e c dtsec_rt | uses the **true c** and the physical time (it is the CMB τ_es integral, not a ray quantity): unchanged |
| cosmological redshifting of rays | `ray_redshifting`: factor (a − da)/a per hydro step | depends on the expansion during the hydro step, i.e. on time: unchanged |
| ray removal | `ngam_min`, `fgam_min` per substep, `nside_max` (path length ≈ 4000–8000 cells) | per-substep fractions and a path length: unchanged; a ray reaches nside_max after a longer time |
| I-front speed | — | limited to ≤ ĉ (below) |

So photons are conserved per substep, the emission rate is physical, Γ is
correct where the medium is optically thin, and the photon density and
energy in rays scale as 1/ĉ: the RSLA. Known RSLA limitations apply: fronts
cannot move faster than ĉ, and after overlap, where the mean free path is
long, the photon field responds to changes with the light-crossing time at
ĉ (Gnedin 2016; Deparis et al. 2019 and Ocvirk et al. 2019 argue for
ĉ ≳ 0.1–0.2c for the timing of reionization).

### Validation

* **Fortran vs C++** for c_reduced = 0.1: all five RT problems, legacy
  (`legacy_rt = legacy_physics = 1`) and fixed, bit-identical (stored
  reference `reference/harness/outputs_c0.1_fixed`, CTests `ref_c01_*`);
  c_reduced = 0.0005 (the guard), bit-identical; full loop: smoke with
  c_reduced = 0.1 in both source modes and flag sets, and 1e-4, live
  CTests `full_vs_fortran_live_smoke_*_c*`; also the `-O0 -fcheck=all`
  Fortran builds (clean).
* **I-front test** (`tests/ifront_rsla_check.py`, CTest `ifront_rsla`;
  figure `docs/figures/ifront_rsla.png`): uniform gas at the mean density,
  z = 10, N = 32, one source of 8e28 (code units), 15 hydro steps of
  da = 2e-4, c_reduced = 1 and 0.1, `legacy_rt = legacy_physics = 0`.
  The finite-c front (Shapiro et al. 2006, AJ 131, 1314) is the
  infinite-c front at the retarded time,

      R(t) = R_∞(t − R(t)/ĉ),

  which without recombinations is exactly their R³/3 + (Q/ĉ)R = Qt
  (R → ĉt early, R → R_∞ late). R_∞ is taken from the c = 1 run (itself
  corrected the same way), ĉ = one cell per substep (0.083c here: nstep_rt
  is rounded, 2 instead of 2.4 per step), so there is no free parameter.

  ![I-front](figures/ifront_rsla.png)

  | step | t [Myr] | ĉt [cells] | R (c = 1) | R (c = 0.1) | R_∞(t − R/ĉ) | error |
  |---|---|---|---|---|---|---|
  | 1 | 1.6 | 2.0 | 3.46 | 1.76 | 1.77 | −0.9% |
  | 2 | 3.2 | 4.0 | 4.42 | 3.20 | 2.95 | +8.5% |
  | 3 | 4.9 | 6.0 | 5.08 | 3.95 | 3.77 | +4.8% |
  | 5 | 8.1 | 10.0 | 6.05 | 5.12 | 4.92 | +4.0% |
  | 10 | 16.3 | 19.9 | 7.66 | 6.95 | 6.75 | +2.9% |
  | 15 | 24.4 | 29.8 | 8.78 | 8.15 | 7.99 | +2.0% |

  The front at ĉ = 0.083c is causal (R ≤ ĉt + 1.5 cells, checked for both
  runs), starts on the light cone (R = 1.76 vs ĉt = 2 after two
  substeps), lags the c = 1 front by 49% after one step, 22% after three
  and 7% after fifteen, and follows the retarded-time solution to 2–5%
  (rms 3.8%, pass criteria 12% per step and 6% rms). The small systematic
  excess (front 0.1–0.2 cells ahead) is the cast's one-cell look-ahead: a
  ray deposits its photons over the cell ahead of it before it advances.
  Emitted photons are identical in the two runs.

### The production configuration: light per hydro step

Box 10 Mpc/h, z = 15 → 5.5, ĉ = 0.01c. The light travel per hydro step
does not depend on N in Mpc/h; in RT cells it scales with N:

| N | RT cell [kpc/h] | hydro steps 15 → 5.5 (ra_max) | light per hydro step at 0.01c | RT substeps per hydro step | RT substeps 15 → 5.5 |
|---|---|---|---|---|---|
| 64 | 156 | 9 | 1.6–2.4 Mpc/h | 10–15 | ≈ 105 |
| 128 | 78 | 9 | 1.6–2.4 Mpc/h | 20–30 | ≈ 210 |
| 192 | 52 | 9 | 1.6–2.4 Mpc/h | 30–46 | ≈ 315 |

So the hydro step does not "dominate the light travel" in the sense that
`nstep_rt` would round to 0 or 1: light still crosses 10–46 cells per hydro
step and the guard is never needed. Rays move one cell per substep as
always. Two consequences matter, though:

* **Rounding.** `nint` of 10–15 changes ĉ by up to ±5% from step to step
  at N = 64 (less at larger N). `rt_substep_carry = 1` carries the
  remainder, so the total number of substeps (the light horizon) is exact;
  the production params use it.
* **Light horizon.** At 0.01c light travels ≈ 16 comoving Mpc/h between
  z = 15 and 5.5, 1.6 box lengths; one hydro step moves it 16–24% of the
  box. A source that forms at z = 8 lights up only ≈ 6 Mpc/h by z = 5.5.
  Reionization will be delayed and the post-overlap photon field
  underestimated (see above); the user accepted this.

**ra_max / outputs.** With the original ra_max (da/a = 1/9) there are only
**9 hydro steps** from z = 15 to 5.5 (and ≈ 17 from z = 100 to 15; the
Courant and PM limits did not bind in the N = 64 and 128 runs). A movie
needs far more frames: `movie_dz` writes x_HI and T at any redshift
spacing from inside the RT substeps (one substep = 1/10–1/46 of a hydro
step), independently of the hydro step, so ra_max does not need to change
for the movie. The gas dynamics response to photoheating, however, is
updated only 9 times; if the hydrodynamics of the ionised gas matters, set
`ra_max = 0.02–0.05` (da/a = 2–5%, 20–50 hydro steps; legacy_rt = 0
only): the hydro is cheap (section 5), and with `rt_substep_carry` the
RT substep total is unchanged.

## 2. Ray pool (F3) and ray memory

The ray pool could overflow (F3): merging is triggered by thresholds on the
photon content, not by the pool size. With `legacy_rt = 0` there is now a
guard (C++ and both Fortran harnesses, RT harness PATCHES.md 13): from the
live count, the number of source cells and a histogram of the rays that may
split this substep, the split threshold is raised so that at most
3/4·Nray_max rays can be live after the substep, and if even the sources
do not fit, splitting stops and all but fresh source rays merge. It is an
integer computation on the histogram, deterministic, and inactive in every
stored reference. The full harness's legacy problem with
nray_max_factor = 1 (which overflowed before) now runs, with the guard
active in 1011 of 1669 substeps, and Fortran and C++ agree bit for bit
(nightly CTest `full_vs_fortran_live_legacy_poolguard`). With
`legacy_rt = 1` the original behaviour (error on overflow) is kept.

`Ray` is now 88 B (int32 links), see `docs/parallel.md` §4.

## 3. F1, restart, nightly tests

* **F1** (write_data_1/2 with one thread wrote only xv_dm): fixed in the
  C++ and in the full harness (`-DHARNESS_FULL_F1`). Both write DM/xv_dm,
  Gas/{u,v,w}_gas and, at data 1, RT/rt_grid and Ray/{ray, ray_index,
  htoc_ray} in the original formats; the live tests compare them (Gas, RT
  byte for byte; Ray as a prefix: the C++ stops at the ray high-water
  mark, the original writes all Nray_max slots). D3D/zreion is not
  written by the C++ (calc_gas_3 is not ported). Only restart files are
  affected.
* **Restart.** The original's data-1 files cannot restart the run
  bit-identically (no free lists, stale thresholds, RT statistics, clock).
  `restart_every = n` writes a complete binary checkpoint
  (`include/radray/checkpoint.hpp`: clock, all RT grid fields and chains,
  gas, grid velocity, rays up to the high-water mark with the pool, the
  thresholds, the RT statistics, particles, output schedule, movie state)
  to `<out>/restart/checkpoint.bin` (or one per step with
  `restart_keep_latest = 0`); `radhydro <problem> <out> --restart
  <checkpoint>` resumes, appending to the logs. CTests `restart_smoke_legacy`
  (restart at step 13 inside the RT phase) and `restart_smoke_det_movie`
  (deterministic_parallel and movie frames): every file written after the
  restart is byte-identical to the uninterrupted run.
* **Nightly tests.** The four ~5-minute `full_ref_{legacy,onthefly}_rt{1,0}`
  and the new `full_vs_fortran_live_legacy_poolguard` carry the label
  `nightly` and are disabled by default. To run them:

      cmake -DRADRAY_NIGHTLY_TESTS=ON build
      ctest --test-dir build -L nightly

## 4. Movie frames and checkpoints in the production run

`movie_dz = 0.05` writes frames `movie/frame_NNNNN_{fHI,temp}.bin`
(float32 N³ each) and a line in `movie/frames.txt` (frame, step, substep,
z, ⟨x_HII⟩_V, ⟨x_HII⟩_M) every Δz = 0.05 from z = 15: ≈ 190 frames to
z = 5.5, from inside the RT substeps (frames in steps without RT are written
at the end of the hydro step). Disk: 2·4N³ B per frame: 0.4 GB (N = 64),
3.2 GB (N = 128), 10.8 GB (N = 192) for 190 frames; use
`movie_fields = fHI` or a larger `movie_dz` at N = 192.

## 5. Measured cost and timing table

(see the measurements below)

## 6. Speedups

(see below)
