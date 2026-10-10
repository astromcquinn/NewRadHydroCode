# Project status

## Production run at 256³ (branch `production-256`, started 2026-10-10)

The production thread of the 2026-10-09 list ("Threads to launch", item 2),
run on klone from the worktree `/gscratch/astro/mcquinn/radhydro/prod256-worktree`
with runs under `$RUNS = /gscratch/astro/mcquinn/radhydro/runs` and Slurm logs
in `/gscratch/astro/mcquinn/radhydro/logs/prod256/`. Settings (user decisions):
256³, 10 Mpc/h, `c_reduced = 0.1`, `nray_max_factor = 16`,
`legacy_rt = legacy_physics = 0`, `deterministic_parallel`, `rt_substep_carry`,
`ray_rotate = 3`, `sfr_model = abundance`, log ξ_ion = 25.38, f_esc calibrated
to τ_es = 0.054, and **`hydro_scheme = unsplit`** (the recommendation of
`docs/hydro_unsplit.md` §10; the decision placeholder in the thread brief was
not filled in, so this thread took the documented recommendation; every
earlier 256³ run used `legacy`). `make_production.py` gained `--hydro-scheme`.

**Step 0, recovery of the earlier attempt (2026-10-09 16:35–16:46).** `sacct`
shows no production, calibration or no-source run from that attempt: its
build jobs 41686559/60 were cancelled, 41686777 is the fast suite (the 10
known stored-reference failures, hence "FAILED") and 41686778 the nightly
suite (19/19). The run directories it left: `$RUNS/cal128_c01/` (only the
`ics -> ../reion128/ics` link, no runs) and `$RUNS/hmf256/` (ICs and a
`params.txt` with `tvir_min = 1e30`, never run; its ICs are byte-identical
to `reion256/ics`, seed 7). Its draft `--c-reduced` option for
`calibrate_fesc.py` sits uncommitted in the main checkout; this thread uses
`--extra c_reduced=0.1` instead. Nothing was reusable except the ICs: the
finished 256³ runs of `docs/cluster.md` §7 and the rt_methods / recomb
siblings all use the global source model and `ray_rotate = 0`, so the
no-source, calibration and production runs below are all new. The older
root `/gscratch/astro/$USER/runs` no longer exists (moved to
`radhydro/runs` on 2026-10-09).

**Step 1.** `ray_rotate` was not on `main` at the start of the thread; while
this thread's merge of `main` into `claude/exciting-clarke-zn62ih` was in
its fast suite (job 41705532: 94/104, the known 10), the B16 thread merged
the same branch as PR #20. The duplicate merge was dropped and
`production-256` branched from `main` at 6dab3d1.

**Progress (updated as the thread goes):**

| step | state |
|---|---|
| 2. build and test `main` on a node | done: fast suite 94/104 (job 41705708; the 10 known stored-RT-reference failures), nightly 20/20 (job 41705709, 58 min); `docs/cluster.md` §3 |
| 3. 256³ mass function | done: `docs/hmf_validation.md` (256³ section), `data/hmf/fof_256_10mpc.txt`, `data/uv_lf/am_params_256_10mpc.txt`; run `$RUNS/hmf256` (2 min 15 s on 96 threads) |
| 4. f_esc at 0.1c | done: **f_esc = 0.28** (`docs/production_run.md` §8, runs `$RUNS/cal128/fesc_*`, 9–10 min each on 96 threads); the 256³ τ_es check is the production run itself |
| 5. production run | `$RUNS/reion256_prod` running on a whole node (job 41706017); `docs/production_256.md` when it finishes |
| 6. c convergence | 128³ runs at f_esc = 0.28: 0.1c (`$RUNS/cal128/fesc_0.28`, done, τ_es = 0.0541), 0.3c (`$RUNS/cal128_c03`, job 41706019) and 1.0c (`$RUNS/cal128_c10`, job 41706020) running |
| 7. ray-precision check | done: `docs/ray_compression.md` §6.2 (`$RUNS/p128/prec`, job array 41705715) |
| 8. STATUS, PLAN, PR | pending |


*Last updated: 2026-10-10: B16 fix `ray_rotate` merged into `main` (branch `claude/exciting-clarke-zn62ih`; production uses `ray_rotate = 3`). The same day: the framework-spike Hyak CPU and GPU benchmarks (PR #19, branch `spike-hyak-bench`: AMReX approved, no ACCESS proposal for now), case A recombination and recombination radiation (branch `recomb-radiation`), and moment RT (M1, OTVET) with its 256³ klone comparison (branch `claude/elegant-bell-rqd3l1`). Before that, 2026-10-09: the survey (PR #16), `hyak-setup` (PR #17), `movies` (PR #12) and `production-run` (PR #7).*

## Done, merged into `main` 2026-10-10 (branch `recomb-radiation`)

**Case A recombination and recombination radiation transported by the rays
(`docs/recomb_radiation.md`).**
- New flags, all off by default:
  - `recomb_case = A|B`;
  - `recomb_rad = 1`, which needs case A, the rays and `legacy_rt = 0`;
  - `recomb_rad_nside`, `recomb_rad_block`, `recomb_rad_pool`.
- How it works:
  - Each ionised cell's ground-state recombinations (α_A − α_B, for HII,
    HeII and HeIII) are buffered.
  - Each b³ block launches them as 12·nside² rays from one rotating cell,
    while there is room in the ray pool. The stellar splits are reserved
    first.
  - Blocks that are optically thick, too faint, or have no room fall back
    to case B on the spot, so no photons are lost.
  - Checkpoint version 4.
- With the defaults the code is bit-identical. The fast suite gives
  92/102: the same 10 known 1-ulp failures, with mismatch details
  byte-identical to `main`.
- New tests:
  - fast: `thread_rt_recomb_lognormal_16src`, `restart_smoke_recomb`;
  - nightly: `recomb_stromgren`, run by `tools/analysis/recomb_stromgren.py`.
- Strömgren check: case A gives R_A/R_B to 1% of the temperature-corrected
  value, and case A + RR recovers R_B (1.002 with block 1).
- 256³ at 0.1c (block 4), runs
  `/gscratch/astro/mcquinn/radhydro/runs/reion256_c01_caseA{,_rr}`, comparison
  in `.../recomb_reion256_c01/`. Movie: `docs/movies/recomb_256_c01.mp4`.
  - Case A alone delays reionization by Δz ≈ 0.25 and lowers τ_es from
    0.0904 to 0.0880.
  - Case A + RR follows case B to Δz ≤ 0.013 up to overlap (τ_es 0.0905),
    with ionised gas 1–3% hotter.
  - After overlap it moves toward case A, as it physically should, since
    the ionised IGM is thin to its own recombination photons. At z = 6 the
    ionised gas is 14% hotter than case B (case A: 25%), and the residual
    x_HI is about α_A/α_B higher.
  - The pool is full during overlap, so only 53–56% of the recombinations
    up to then are launched; 82% are by the end.
  - Cost: 1.46× wall (2:16:39 against 1:33:51), with 13% more memory.
- Tools: `tools/hyak/recomb_runs.sh` and `tools/hyak/recomb_post.slurm`.
- Found, not caused by this work: with `deterministic_parallel`,
  `lognormal_16src_smallbuf` (`nray_max_factor = 1`) overflows the ray pool
  on `main` too.
- Open:
  - transport HeII Lyα and the other ionising recombination lines;
  - a larger pool, or coarser emission during overlap, so that more of the
    emission is launched;
  - tag recombination rays to measure their losses separately.

## Done, merged into `main` 2026-10-10 (branch `claude/elegant-bell-rqd3l1`)

**Moment RT: M1 and OTVET as alternatives to the ray tracing
(`docs/rt_moments.md`).**
- `rt_method = rays (default) | m1 | otvet` in `rt_driver` and `radhydro`.
  A two-moment photon solver per frequency group on the RT grid
  (`src/rt_moments.cpp`): explicit first-order GLF transport at ĉ,
  `moment_nsub = 4` subcycles per RT substep, and the M1 closure or the
  optically thin variable Eddington tensor (direct sum or FFT, every hydro
  step).
- Sources, `photo_cell` (absorption through `rt_grid%fgam`) and the redshift
  are the rays' physics, unchanged. `ray_redshifting`'s spectrum update is
  factored out (`redshift_spectrum`) with identical arithmetic.
- The defaults are untouched: the full fast suite, 93 tests including the
  live gfortran comparisons, passes. New: `ctest -L moments`, 7 tests (unit,
  I-front vs finite c, threads, restarts). Checkpoint version 3.
- Results:
  - Photons are conserved to float precision.
  - Thread-count independent; bit-identical restarts.
  - I-fronts converge to the finite-c solution: +9% (M1) and +7% (OTVET)
    at 10 cells of light travel, +2.5% and +1.6% at 30. The first-order
    diffusion makes small fronts run ahead (+63% at 2 cells).
  - Shadows behind a dense clump: OTVET sharp, M1 slightly filled in
    (x_HII 0.06).
  - 16-source lognormal box: the moment methods end 4–8% above the rays in
    x_HII, with about 1.6× thicker fronts.
  - 5.5–6.5× faster than the rays on that box (1 thread, 32³), at about
    0.7 µs per cell-substep.
- **256³ production comparison on klone, done** (`docs/rt_moments.md` §6).
  M1 and OTVET siblings of `reion256_c01` (0.1c) and `reion256` (0.01c),
  made by `tools/hyak/rt_methods_runs.sh` and post-processed by
  `tools/hyak/rt_methods_post.slurm`. Runs are in
  `/gscratch/astro/mcquinn/radhydro/runs/reion256{,_c01}_{m1,otvet}`, the
  comparison in `.../rt_methods_{c01,c001}/`.
  - History: |Δx_HII,V| ≤ 0.018 from the rays at every z, and τ_es within
    1% (0.1c: 0.0904 / 0.0905 / 0.0906). The end of reionization moves by
    Δz ≤ 0.12.
  - Morphology: OTVET is much closer to the rays than M1. It has 2–3× fewer
    cells off by more than 0.1 in x_HII, and corr(log T) ≥ 0.87, against
    0.08–0.6 for M1 around overlap. The last neutral islands sit in
    different places with M1.
  - The temperature of ionised gas agrees to 1–2%.
  - At 0.01c only, the small early fronts run ahead (+20% x_V at z = 13).
  - Cost: 1.8× less RT time than the rays at 0.1c (49:35 and 51:39 wall,
    against 1:32:42), and 2.4–3.3× less at 0.01c. Half the memory.
    - The transport, 0.29 s per substep on 192 threads, meets the
      projection.
    - The unchanged `radiative_physics` is now 51–55% of the moment RT time.
    - The cost per substep follows the ionised volume, not N_source.
  - Side-by-side movies (rays | M1 | OTVET): `docs/movies/rt_methods_256_{c01,c001}.mp4`,
    made with the new `tools/movies/compare_movie.py` (`docs/movies.md`).
  - Fast suite on klone: 90/100. The 10 failures are the known 1-ulp
    stored-RT-reference platform differences (`docs/cluster.md` §3, plus
    the two ray-dump tests on the same references).
  - `rt_methods_runs.py` now compares at matched movie frames. The zdata1
    output redshifts drift apart once the hydro steps differ.
- Open: second-order transport; Iliev-test validation; M1 crossing beams;
  the OTVET formulation (time-dependent two-moment here; the published
  implementations differ and could not be checked, no journal access).

## Done (merged into `main`)

**Port of the original, bit-identical to gfortran builds (single-threaded, FFTW in place of MKL):**

| Component | C++ | Fortran reference | Notes |
|---|---|---|---|
| RT: ray tracing with HEALPix split/merge, photo-chemistry/heating | `apps/rt_driver` | `reference/harness` | `docs/port_notes.md` |
| Hydro (`hydro_tools.f90`, -DMF) | `apps/hydro_driver` | `reference/harness_hydro` | `docs/port_notes_hydro.md` |
| PM gravity (`pm_tools.f90`, `fft_tools.f90`) and Zel'dovich ICs | `apps/pm_driver`, `tools/ics` | `reference/harness_pm` | `docs/port_notes_pm.md` |
| Full `radhydro.f90` loop | `apps/radhydro` | `reference/harness_full` | `docs/full_run.md` |
| On-the-fly FoF halo finder, source model, legacy catalog reader/writer | `src/halos.cpp`, `src/halo_sources.cpp` | Fortran calls the same finder via a C API | `docs/halos_and_sources.md` |

The full loop works with legacy catalogs and with on-the-fly sources.

**Flags.** In both codes, `legacy_rt=1` and `legacy_physics=1` (the defaults)
give the original code.

`legacy_rt=0` turns on these RT bookkeeping fixes:
- B1: freed rays no longer steer the split/merge thresholds;
- B2: the thresholds are reset every substep;
- B3: the direction field is reset;
- B4: τ is updated in all cells;
- B8: up to 4 segments per step;
- B9: a 32-bit overflow is avoided.

`legacy_physics=0` turns on these physics fixes:

| Fix | Change |
|---|---|
| P1 | case-B recombination cooling |
| P2 | HeI 2³S excitation |
| P3 | HeII dielectronic recombination cooling |
| P4 | HeIII recombination typo |
| P5 | collisional-ionisation cooling uses the ionisation potential |
| P6 | HeII dielectronic recombination in the ionisation balance |
| P7 | HeII case A (Hui & Gnedin 1997) |
| P8 | T-dependent Gaunt factor |
| P9 | interpolated rate tables |
| P10 | HeI ground-state excitation (Glover & Jappsen 2007) |

Recombination radiation stays off. Deliberate model choices are unchanged:
overcooling suppression, the rho_source photon consumption, SN heating.

**Tools:** the comparison tool, quick-look figures
(`tools/analysis/test_gallery.py` → `docs/figures/tests_*.png`), and the cooling
figure.

**Ray compression study (ROADMAP §4, Stage B): done** (PR #13, branch
`ray-compression`); see `docs/ray_compression.md`.

- Optional ray snapshots (`ray_dump_*`) and a research-only lossy-precision
  instrument (`ray_quant_*`). Both are off by default, and the C++ stays
  bit-identical to the Fortran (`ctest -L raycomp`).
- Never-merged rays are determined by (source, age, nside, ipix) to
  float32 precision.
- Ray size:
  - fixed-width "R-struct": 29–31 B per ray (80 B of ray data, 88 B per slot);
  - rANS prototype (`build/ray_codec`): 5–9 B before overlap and 12.7 B
    after, with a 6-bit dithered photon total and 10-bit offsets. It costs
    230–310 ns to decode and 290–450 ns to encode per ray, against ~2 µs of
    RT per ray-substep.
- Dithered photon quantisation: 8 bits is at the 10⁻⁴ level in every
  test including the 64³ production configuration (0.01c); 1–2 bits pass
  the 32³ tests but reach 1 % and 15 % per-cell T scatter at 0.01c.
  Round-to-nearest is biased (≈1 % error at 6 bits).
- Many sources (128³, up to 1506 sources at the dumps): the source id
  costs ~3 bits per ray; the spectra diversify and the models must be
  refit every hydro step (redshifting), after which the size is the same
  12–14 B per ray as with 17 sources.
- Recommendation: R-struct with 8 dithered photon bits as the Stage D
  in-memory format, and rANS for checkpoints and MPI.
- The 128³ production-configuration precision check
  (`tools/production/hyak_precision.sbatch`) was run on 2026-10-10
  (`docs/ray_compression.md` §6.2): the 64³ conclusion holds.

**Design doc for Phases 2–3:** `docs/design_parallel_and_encoding.md`
covers Kokkos, MPI, sorted ray arrays, compact ray encoding and adaptivity.

**Production readiness (PR #10; `docs/production_run.md`, `docs/parallel.md`):**
- `c_reduced`: reduced speed of light via the original's `cscale`. Bit-identical
  between C++ and both Fortran harnesses. Includes an nstep_rt ≥ 1 guard,
  `rt_substep_carry` for exact total light travel, and an I-front test against
  the finite-c solution.
- OpenMP in hydro, PM, summaries, FoF and all RT grid loops. These are
  bit-identical to serial for any thread count.
- Ray tracing runs in parallel only with `deterministic_parallel = 1`
  (needs `legacy_rt = 0`). It uses 8³ domains swept by colour and is
  thread-count independent. It differs from the legacy order at about 1e-7.
  The 4-core speedup was not measured on an idle machine; re-measure on Hyak.
- Ray-pool guard (F3, `legacy_rt = 0`). Rays are now 88 B.
- F1 fix; bit-identical restarts (`restart_every`, `--restart`).
- `movie_dz` frames, written from inside the RT substeps; `ra_max` parameter
  (`legacy_rt = 0` only); `fields` and `data_files` output options.
- `nightly` CTest label (`-DRADRAY_NIGHTLY_TESTS=ON`).
- `tools/production/`: `make_production.py` (parameters), the timing table
  and the speed-up scripts.

**Hydro test suite (PR #10; `docs/hydro_tests.md`, `tests/hydro_suite/`):**
- Covers Sod (axes and diagonals, strong shock, 123), Galilean invariance,
  Sedov, sound waves, the Zel'dovich pancake with gas, an adiabatic
  cosmological box, and Kelvin–Helmholtz.
- **Key finding: the original scheme spuriously heats cold supersonic gas.**
  The temperature of unshocked cold gas is wrong by factors of 10–1000:
  underdense gas sits near 10⁴ K instead of 0.02–2 K, and pancake gas before
  the caustic is too hot. The error falls steeply with resolution only for
  smooth flows. There is no entropy fallback. **Don't trust pre-reionisation
  gas temperatures** from this scheme.
- Densities, velocities and shocks are accurate (Sod L1 order ≈ 0.8, Sedov
  radius to 1% at 128³, sound-wave order ≈ 2).
- Other issues: energy is not conserved (up to 5% in Sedov); contacts smear
  over about 5 cells; orientation dependence (H2); the Courant step shrinks
  with the bulk speed.
- These matter for the plane-parallel / DM20-type runs, whose cold
  pre-reionisation gas structure is the point. A fix, e.g. a dual-energy /
  entropy formulation for cold supersonic gas behind a flag, is a candidate
  next task.

**Unsplit hydro scheme with dual energy (PR #15, branch `hydro-unsplit`;
`docs/hydro_unsplit.md`, `docs/hydro_tests.md` §12):**
- `hydro_scheme = legacy | unsplit` (default `legacy`, bit-identical to the
  Fortran; every stored-reference, live-Fortran, restart and thread test
  passes unchanged). `unsplit` is the spike's MUSCL-Hancock + HLLC kernel
  (`include/radray/kernels/hydro_mh.hpp`, header-only, the AMReX thread
  takes it from here) with the dual-energy formulation of Bryan et al.
  1995 (`GasGrid::e`, η₁/η₂ against the ±2-cell maximum, a total-energy
  reset) and the species advected with the mass flux. Two sub-cycled
  unsplit steps around the PM kick, same call sites, same `calc_timestep`.
- **The spurious heating of cold supersonic gas is fixed**: pancake infall
  median T/T_adiabatic 1.000 at every N (legacy 7–20 at N = 64), box
  median 2–9 K at z = 10 (legacy 6500 K), coupled 64³ run 2.7 K at z = 14
  (legacy 133 K). Also: energy conserved to round-off except ahead of
  strong shocks in Mach ≫ 1 gas (U4, ≤ 0.6 %), contacts 3.4 cells, exact
  symmetry and orientation independence for 1-D problems, 2–3× less
  dissipative sound waves, KH grows at 1 cell per layer, half the cost per
  hydro step. Worse: Sedov at 128³ (radius +1.5 % vs +0.9 %, slower
  convergence), two hot stagnation cells before a caustic (U2), the
  cold-gas temperature depends on the hydro step (use `ra_max` ≈ 0.02).
- The RT's `photo_cell` reads and writes the field (the one physics
  coupling point); the RT grid temperature and `summarize_gas` read it for
  diagnostics. Checkpoint format version 2 (version 1 still readable).
- **Constants confirmed** against the Enzo (Bryan et al. 2014 §4.1.1) and
  Nyx (Almgren et al. 2013 §3.2) method papers: η₁ = 10⁻³, η₂ = 0.1, the
  total-energy reset, and tested: the exact Enzo form reproduces the
  cold-gas results (energy to round-off, but 28× stagnation cells at z = 10
  in the pancake); the Nyx resync rule does not work in this kernel (box
  underdense gas at 5000 K). Deliberate deviations (the ±2 all-axis
  neighbourhood, the time-centred pressure in the work term) are measured
  and documented in `docs/hydro_unsplit.md` §2; all are parameters.
- Suite: `run_suite.py --scheme unsplit`, tolerances per scheme,
  `hydro_suite_unsplit_*` (fast) and `hydro_suite_full_unsplit_*`
  (nightly), `restart_smoke_unsplit`, `thread_full_unsplit_hydro_pm`.

**Hyak klone setup (PR #17, branch `hyak-setup`, 2026-10-08; `docs/cluster.md`):**

- **Environment:** gcc/g++/gfortran 13.2.0 (module `gcc/13.2.0`), FFTW
  3.3.10 float built SSE2-only in user space (AVX2 / no-SIMD builds break 5–6
  stored PM/hydro references), CMake 4.4.4, Python 3.12, numpy, matplotlib
  and ffmpeg from a conda env, all under `/gscratch/astro/mcquinn/sw`.
  `tools/hyak/` holds env, setup, build, production (auto-restart) and
  scaling jobs. `radhydro` has a new `run_steps` parameter (stop after n steps
  of an invocation; for timing).
- **Tests:** fast suite 72/80, nightly 12/12. The 8 failures are RT stored
  references (`ref_*lognormal_16src*`, `ref_c01_*`): 1-ulp differences in
  double-precision stats and 2–10 cells of ray photons at the removal
  threshold. The live gfortran RT harness built on klone shows exactly the
  same differences and agrees with the C++ byte for byte (17/17), so the
  stored references are platform-dependent (glibc 2.28 vs 2.39 is the likely
  cause). They were not regenerated.
- **OpenMP scaling on an idle node:** byte-identical outputs for 1–192
  threads. RT 50× on 96 threads at 128³, flat beyond that (512 domains per
  colour); 94× on 192 threads at 256³ (ray tracing 116×).
- **Runs** (global source model): 128³ 0.01c in 3 min, 0.1c in 16 min;
  256³ 0.01c in 26 min (shared nodes); 256³ 0.1c with movie in 1 h 32 min on
  a whole node, 27 GB peak RSS. **0.01c delays reionization by Δz ≈ 2–2.5**
  against 0.1c and lowers τ_es by 20% (same ICs at each N). The 128³ and 256³ ICs are
  different realisations (make_ics.py phases depend on N), so N = 128 vs 256
  is not a clean resolution test. Both c reionize too early with the uncalibrated global model
  (τ_es = 0.07–0.09).
- **Findings for later:** movie frames written before the RT starts log
  x_HII,V = 1 and x_HII,M = nan. (The other two findings of this thread,
  uninitialised checkpoint padding bytes and the `ctest -L long` label that
  selected no test, were fixed in the 2026-10-09 survey.)
- **Open PRs touched:** #7 (`production-run`, draft): its `tools/hyak/`
  scripts were copied into `hyak-setup` and adapted (env, time limits,
  restart from the newest checkpoint); the `hyak-setup` versions are the
  ones to keep when #7 is merged.

**Movie tooling (PR #12, branch `movies`; `docs/movies.md`):**
- `tools/movies/make_movie.py` renders an `apps/radhydro` run directory to an
  MP4 (ffmpeg, piped raw frames) and a few PNG stills.
- Style `cosmo` shows slices or thin projections of x_HII (or log x_HI), T
  and ρ_gas, with fixed scales and perceptually uniform colour maps. The time
  panel has x_HII,V/M(z), τ_es accumulated, the total τ_es if the IGM were
  ionised below z, and a Planck band. Labels give z, cosmic time and the box.
- Style `planar` shows a slice along the ray direction and the I-front
  position against time, for the future plane-parallel mode.
- Options for zoom, slice axis, position and thickness. Memory is bounded:
  frames are memory-mapped one slab at a time.
- `make_production.py` now writes `rhog` in the movie frames as well.
- Test: CTest `movie_tool` (synthetic run).
- Test movie: `docs/movies/reion64_test.mp4` (0.9 MB), from a 64³ production
  run (90 s run, 16 s render).
- Finding (B16, `docs/port_notes.md`): the original ray tracer leaves a
  photon deficit along ±x/±y from every source. Its HEALPix frame is fixed
  to the grid, so those axes are pixel vertices and the rays around them
  are symmetric. In the movies this shows as grid-aligned neutral streaks
  and "+" crosses on sources. The deficit is up to 0.25 in x_HII near the
  front in the uniform test, and present in the legacy (= Fortran) runs.
  Fixed behind a flag on 2026-10-10 (`ray_rotate`, below).
- Finding: frames written before the first RT step have an uninitialised
  rt_grid (fHI = rhog = 0). The tool skips them; a C++ fix that skips them at
  write time is optional.

**B16 fix: `ray_rotate` (branch `claude/exciting-clarke-zn62ih`, 2026-10-10;
`docs/port_notes.md` B16):**
- `ray_rotate = 1 | 2 | 3` rotates the 192 HEALPix directions a source
  emits: one fixed rotation, a new one every RT substep, or a new one every
  substep and source cell (default 0 = the original, so every reference
  test is unchanged). `ipix` becomes the grid-frame pixel that contains the
  rotated direction, and the existing split deviation handles the rest. So
  there is no frame id and no extra ray memory. The rotation is a hash of
  `ray_rotate_seed`, the substep (restart-safe) and the cell, so it is
  deterministic and independent of the thread count. Works in
  `rt_driver` and `radhydro`, with any `legacy_*` flags. It also rotates
  the recombination-radiation rays (`recomb_rad`, merged from `main`; own
  hash stream; checked thread-count independent). No effect with
  `rt_method = m1 | otvet`.
- Uniform test: the worst axis cell goes from −0.46 to −0.018 in x_HII
  (also −0.47 → −0.024 at `c_reduced = 0.1`). A fixed rotation (1) leaves
  −0.10 on one axis, so use 2 or 3. In the 240-step run with production
  ray settings, the 12-fold ring of neutral, cold cells inside the sphere
  (`docs/figures/raycomp_slices_uniform_1src.png`) is gone: 36 → 0
  interior cells below x_HII = 0.9, and the scatter of log T in the
  interior falls from 0.52 to 0.03 dex. So the ring was B16 too.
- 64³ cosmological run: the streaks, the "+" crosses and the faint grid
  lines through bubbles are gone, and x_HII,V(z) is unchanged. In exchange,
  cell-to-cell noise of log x_HI in ionised gas rises by 0.0–0.03 dex
  (`docs/figures/ray_rotate_stills64.png`).
- Cost: ray tracing takes 3–10 % longer single-threaded
  (`rt_driver`, 30 steps, fixed flags). Part of it is extra work (on
  `lognormal_16src` about 4 % more casts and 3 % more splits; callgrind:
  5–6 % more instructions in the casts, same mispredicts), cause not
  pinned down. A fixed rotation costs the same, so it comes from rays
  that are no longer grid-aligned, not from changing frames. Not visible
  in the 64³ `radhydro` run at 4 threads (RT 40.8 s → 39.1 s, within the
  noise).
- Tests: `rays` (unit: proper rotations, the keys, the emitted rays) and
  `ray_rotate` (physics label, 2 s: the deficit is seen without and gone
  with the flag, same front and photons).
- Tools: `tools/analysis/ray_rotate_compare.py`
  (`docs/figures/ray_rotate_uniform.png`).
- **Production uses `ray_rotate = 3`** (user decision 2026-10-10);
  `make_production.py` writes it (`--ray-rotate 0` for the original). The
  user notes that the single-source tests exaggerate the effect for a full
  cosmological run, where many overlapping sources already fill each
  other's gaps. The f_esc calibration (`calibrate_fesc.py` calls
  `make_production.py`) now includes the rotation; the 64³ history did not
  change with it.

**Production-run thread (PR #7, branch `production-run`; merged 2026-10-09):**

*Abundance-matched sources* (`docs/abundance_matching.md`):
- `sfr_model = abundance`, a new source model; the existing models are
  unchanged and bit-identical.
- L_UV(M_halo, z) from abundance matching to the observed UV LFs, with all
  LF parameters read from the papers the user supplied: Bouwens+22
  (z ≈ 6, 8), GLIMPSE/Atek+26 (z ≈ 7), Donnan+24 (z = 9–12.5). The shape
  is fitted at z ≤ 9 with the faint end extrapolated to M_UV = −12; at
  z = 10–12.5 only the normalisation is fitted; L(M) is held fixed beyond
  the data. L(M) is matched to the simulated FoF mass function
  (`data/hmf/fof_128_10mpc.txt`; user decision), so the box's LF matches
  the observations. The fit is `data/uv_lf/am_params_128_10mpc.txt`.
- Ionizing rate = f_esc ξ_ion L_UV with log ξ_ion = 25.38 (Pahl et al.
  2025). The RT code is untouched.
- Tools: `apps/am_fit`, `apps/hmf_table` (Tinker in `MassFunction`),
  `tools/analysis/plot_am.py`, `box_lf.py`, `fit_fof_hmf.py`,
  `hmf_compare.py`, `tests/test_abundance_matching.cpp`.

*Halo mass function* at 128³ / 10 Mpc/h (`docs/hmf_validation.md`):
- FoF agrees with Sheth–Tormen above about 100 particles and is 1.3–2.2×
  above Tinker; a 2–4× deficit at 20–60 particles grows with z. The first
  halo forms at z ≈ 11.3. Not yet checked at 256³.

*f_esc calibration* (`docs/production_run.md` §8, folded in from the
thread's `reionization_production.md`):
- At 128³ and **0.01c**, f_esc = 0.5 gives τ_es = 0.0516 with the final
  LFs (Planck 2018: 0.054 ± 0.007), midpoint z = 6.72, x_M = 0.967 at
  z ≈ 5.3. τ is a weak constraint (f_esc 0.2–0.8 gives 0.046–0.058); the
  end of reionization is the sharper one.
- The production run was to start at f_esc = 0.5 (user decision), but the
  calibration **must be redone at 0.1c** (`docs/cluster.md` §7: 0.01c
  lowers τ_es by 20 %) and at 256³; `tools/production/calibrate_fesc.py`.
- Movie, test movie `docs/figures/reion_test_128.mp4` (the thread's own
  movie tool was dropped in favour of `tools/movies/make_movie.py` from
  PR #12; `reion_production.slurm` now runs `reion.py`, the stills and the
  movie when a run finishes). `make_production.py` gained `--am-params`,
  `--xi-ion`, `--f-esc` and `--restart-every`.

## Measured production costs (PR #10, `docs/production_run.md`)

Box 10 Mpc/h, Ngas = Ndm, factor 16, z = 100 → 5.5, fixed flags,
`deterministic_parallel`:

| N | c_reduced | Measured / projected time | Memory |
|---|---|---|---|
| 128 | 0.01 | 69 min measured on about 2 cores; about 30 min projected on 4 | 2.8 GB |
| 192 | 0.01 | about 2.5–3 h on 4 cores | 8–11 GB; checkpoints of 6–9 GB each |

At 0.01c the light horizon is only about 16 Mpc/h over z = 15 → 5.5.

## Framework spike (Stage C), PR #11 (branch `framework-spike`)

- `spike/`: AMReX and Parthenon prototypes of (a) an unsplit MUSCL-Hancock
  hydro step and (b) the sorted-ray-array RT substep with HEALPix split and
  merge, int64 fixed-point deposits and ray migration between ranks, all
  built on shared host/device kernels (`spike/common`).
- Validation:
  - the serial sorted-array reference `rt_serial` matches `rt_driver
    --set legacy_rt=0` with identical ray statistics on all five RT
    problems, fields to round-off (the explained exceptions are in
    `docs/framework_decision.md` §2);
  - both prototypes are bitwise identical to `rt_serial`/`hydro_serial`
    for 1–8 ranks, several block sizes and 1–2 threads.
- Local 4-core timings: the RT on one core costs serial 19.6 s, AMReX
  26.6 s, Parthenon 35.7 s (64³, 400k rays, 32³ blocks).
- **Recommendation: AMReX** (`docs/framework_decision.md`). Parthenon swarms
  do not support AMR, there is no sub-cycling, and the mesh overhead was
  higher. **Approved by the user on 2026-10-10.**
- Hyak CPU (2026-10-09, branch `spike-hyak-bench`, `docs/framework_decision.md`
  §3.2–3.3; built and validated bitwise on the node, AMReX 49/49, Parthenon
  52/52):
  - best layout 24 ranks × 8 threads per node (one rank per CCD);
  - RT 256³ on one node: AMReX 550 s, Parthenon 572 s, OpenMP-only
    `rt_serial` 738 s; at 128³ Parthenon is 13% faster than AMReX;
  - 1 → 4 nodes: weak scaling 75% (both), strong scaling at 256³ 53–54%
    (too small a problem for 96 ranks);
  - open: AMReX *hydro* across nodes is ~10× slower per cell than on one
    node (Parthenon is not); investigate before Stage E.
  The numbers do not change the recommendation. **AMReX approved by the
  user on 2026-10-10** (`docs/PLAN.md`).
- GPU (2026-10-10, §4): the Parthenon prototype runs on A100 and H200
  (Kokkos CUDA, CUDA 12.4.1) after two fixes: `--fmad=false` and native
  64-bit `atomicAdd` (Kokkos used CAS loops). Ray statistics are identical
  to `rt_serial`; deposits differ at the 1e-16 level from CUDA's `exp`
  (explained in §4.2). One GPU vs one `cpu-g2` node on the 128³ bench:
  **H200 5.0–5.9×, A100 3.3–3.8×**, so the ACCESS ≥ 5× rule is just met
  on H200-class GPUs, untuned. The AMReX prototype still has no CUDA
  target. **No ACCESS proposal for now** (user, 2026-10-10).

## Branches and PRs (merged on 2026-10-09)

The survey found three branches off `main`; they were merged in the order
`hyak-setup` (PR #17), `movies` (PR #12), `production-run` (PR #7), with
these resolutions: the `tools/hyak/` scripts of `hyak-setup` were kept
over the older copies on `production-run` (plus the post-run `reion.py`,
stills and movie calls); `tools/movie/make_movie.py` of `production-run`
was dropped for `tools/movies/make_movie.py` of PR #12;
`docs/reionization_production.md` was folded into `docs/production_run.md`
§8. Remote branches with no commits beyond `main` after the merges:
`framework-spike`, `hydro-unsplit`, `ray-compression`, `hyak-setup`,
`movies`, `production-run`, `claude/hopeful-thompson-w5m60t`,
`claude/relaxed-brahmagupta-9v7xb5`, `claude/lucid-brahmagupta-85vs76`.
All nine were deleted by the user on 2026-10-09. Each merge was
preceded by a merge of `main` into the branch and a full test run in the
cloud container (gfortran 13.3, FFTW 3): fast suite (including the live
Fortran, restart and thread tests) and nightly label, all passing.

**Facts from these branches that change the plan:**
- One klone node runs 128³ at 0.01c in 3 min and 256³ at 0.1c with a
  movie in 1.5 h (192 threads, RT 94× over one core). Resolution is no
  longer the constraint; 256³ is the production size.
- **c_reduced = 0.01 delays reionisation by Δz ≈ 2–2.5 and lowers τ_es by
  20 %** against 0.1c (same ICs). Use 0.1c; check convergence with one
  0.3c or 1c run at 128³.
- The f_esc = 0.5 calibration of PR #7 was done at 128³ and
  0.01c; it has to be redone at 0.1c (and at 256³, whose HMF has not been
  checked).
- On klone 8 RT stored-reference tests fail by 1 ulp in double sums and by
  rays at the removal threshold; the live gfortran harness built there shows
  the same differences and agrees with the C++ byte for byte. The stored RT
  references are platform-dependent (libm). Decision needed: regenerate
  them per platform, or compare the double sums with a tolerance.
- B16 (PR #12, `docs/port_notes.md`): the fixed HEALPix frame leaves a photon deficit
  along ±x/±y from every source. Fixed behind `ray_rotate` (2026-10-10).
  A fixed global rotation was not enough; the frame has to change every
  substep.
- `make_ics.py` draws different phases for different N with the same seed,
  so 128³ against 256³ mixes resolution with cosmic variance.
- The linear growth of the fundamental mode is 2 % low in power by z = 5.4
  at both N and both c (not resolution); cause unchecked.

**Repository cleanup done in this survey:** `tools/compare/compare.py.orig`
removed; `ctest -L long` in CLAUDE.md replaced by the nightly label (no test
carried `long`); README status, layout and dependencies brought up to date;
`.gitignore` covers `out*/` and `runs/`; checkpoint files no longer contain
uninitialised padding bytes (`DriverState` is written from a zeroed copy;
format unchanged).

## Threads to launch (2026-10-09)

On Hyak (`claude remote-control` in tmux, or batch jobs):

1. ~~Merge and validate~~ Done 2026-10-09 (PRs #16, #17, #12, #7). Still
   to do on klone: the fast and nightly suites on the merged `main`
   (expect the 8 known RT stored-reference failures of `docs/cluster.md` §3).
2. **Production run at 256³, 0.1c, abundance-matched sources.** Recalibrate
   f_esc at 0.1c (four 128³ runs of about 16 min, `calibrate_fesc.py`),
   check the 256³ HMF (`hmf_compare.py`), then the 256³ 0.1c run with
   movie on a whole node (1.5 h), and one 0.3c or 1c run at 128³ for the
   c convergence. Decide `hydro_scheme = unsplit` for it (the legacy scheme
   overheats cold gas; `make_production.py` does not set it yet).
3. ~~Framework-spike benchmarks~~ Done 2026-10-10 (branch
   `spike-hyak-bench`; `docs/framework_decision.md` §3.2–4). AMReX
   approved and no ACCESS proposal for now (user, 2026-10-10).
4. **Ray-compression 128³ precision check**
   (`tools/production/hyak_precision.sbatch`), the one unrun item of the
   compression study.

In the cloud (branch + PR each):

5. **B16 fix behind a flag** (frame rotation; `tools/analysis/axis_deficit.py`
   and the 64³ movie as the checks). Blocks science runs.
6. **Stage D: sorted-ray-array RT in AMReX** on a uniform grid, with the
   R-struct encoding (8 dithered photon bits), validated against
   `rt_driver --set legacy_rt=0` as the spike was. Needs the user's AMReX
   confirmation first.
7. **Plane-parallel / RT-domain mode** (`docs/plane_parallel.md`).
8. Smaller: the stored-RT-reference platform policy; the 2 % growth
   deficit; the unsplit open items (Sedov radius convergence, cold-gas
   time-step sensitivity at 128³).

**Decisions needed from the user:** production
`c_reduced` (0.1 recommended) and N (256); `hydro_scheme` for production;
the stored-reference policy on klone; deleting the merged branches.

## Next

1. ~~Halo mass function at production resolution~~ Done at 128³
   (`docs/hmf_validation.md`); redo at 256³ before the production run.
2. ~~Abundance-matched sources~~ Done (`sfr_model = abundance`,
   `docs/abundance_matching.md`). **Open: recalibrate f_esc at 0.1c and
   256³** (`docs/production_run.md` §8, `calibrate_fesc.py`); the 0.01c
   value is 0.5.
3. **Production reionisation run and movie** on one klone node
   (`docs/production_run.md` §7–8, `docs/cluster.md` §5, `docs/movies.md`).

   | Setting | Value |
   |---|---|
   | Box | 10 Mpc/h |
   | c_reduced | 0.1 (0.01 delays reionisation by Δz ≈ 2–2.5, `docs/cluster.md` §7); convergence check with one 0.3c or 1c run at 128³ |
   | nray_max_factor | 16 |
   | Ngas = Ndm | 256³ (1.5 h at 0.1c with a movie on a whole node) |
   | Flags | legacy_rt = legacy_physics = 0, deterministic_parallel, rt_substep_carry, ray_rotate = 3 |
   | Sources | sfr_model = abundance, log ξ_ion = 25.38, f_esc recalibrated at 0.1c |
   | hydro_scheme | user decision pending (`unsplit` fixes the cold-gas overheating; `make_production.py` does not set it yet) |
   | Before the run | 256³ HMF and f_esc |

4. **User's cluster: UW Hyak `klone`, account `astro`.**
   - Owned: about 1000 `cpu-g2` cores (about 5 nodes of AMD Genoa, 192 cores
     and 1.5 TB per node) with 7.8 TB of memory, a 64-core,
     1 TB `cpu-g2-mem2x` slice, and 40 `compute-bigmem` cores.
   - No owned GPUs: GPUs only via preemptible `ckpt` partitions. No time
     limits. See `docs/ROADMAP.md` (answers section).
   - Access: Claude Code installed in user space, `claude remote-control`
     inside tmux, batch jobs with restarts.
5. ~~B16 fix before any science run~~ Done behind `ray_rotate`
   (2026-10-10, see above): streaks, crosses and the 12-fold ring are gone.
   Production uses `ray_rotate = 3` (user decision 2026-10-10), set by
   `make_production.py`.
6. ~~Hydro fix for cold supersonic gas~~ Done (`hydro_scheme = unsplit`,
   PR #15; η₁/η₂ confirmed against the Enzo and Nyx
   papers). Open: convergence check of the Sedov radius bias and the
   cold-gas time-step sensitivity at 128³.
7. **Plane-parallel / RT-domain mode** (`docs/plane_parallel.md`).
8. After the production run, follow `docs/ROADMAP.md`:
   - ray compression study: done (`docs/ray_compression.md`); 64³
     production-configuration check done, 128³ job script ready for Hyak;
   - framework spike (AMReX vs Parthenon): done, AMReX recommended
     (`docs/framework_decision.md`), **approved 2026-10-10**; Hyak CPU and GPU
     benchmarks done (H200 ≈ 5 nodes on the RT; no ACCESS proposal yet);
   - ray data-structure rewrite;
   - MPI;
   - AMR for gas, then RT;
   - multigrid/tree gravity;
   - adaptive RT.

## Open questions for the user

- Moment RT (`docs/rt_moments.md` §5–§6): whether to invest in second-order
  transport or first in faster photo-chemistry (now half the moment RT
  time), and whether the moment methods should become the
  transparent-region hand-over of ROADMAP §4.5 or stay a standalone
  alternative. Which OTVET formulation the user has in mind (the paper
  could not be read here).

- Original `Halo/` catalogs, `sfr.txt`, or the Trac & Cen halo finder. These
  would pin down the original source model; it is reconstructed, see
  `docs/halos_and_sources.md`.
- Whether to keep the original's default `nray_max_factor = 64` in legacy
  mode. Production uses 16.

## Known limits

- **OpenMP only.** MPI/GPU is Phase 2. On one Genoa node the RT reaches 94× on 192 threads at 256³ and about 50× at 128³ (`docs/cluster.md` §6).
- **RT dominates the cost.** Hydro+PM take under 1%. At full c, z=12→6 at
  128³ would take days on one core, hence `c_reduced`.
- **Network access from cloud sessions.** arXiv, ADS and journal sites are
  unreachable.
