# Project status

*Last updated: 2026-10-09, merge of `hyak-setup` (PR #17) and `movies` (PR #12) into `main` (after the repository survey of branch `claude/lucid-brahmagupta-85vs76`): branch and PR state, merge order, Hyak threads and next steps are in the two sections before "Next".*

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

**Ray compression study (ROADMAP §4, Stage B): done** on branch
`ray-compression`; see `docs/ray_compression.md`.

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
- Not yet run: the 128³ production-configuration precision check
  (`tools/production/hyak_precision.sbatch`).

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

**Unsplit hydro scheme with dual energy (branch `hydro-unsplit`;
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

**Hyak klone setup (branch `hyak-setup`, 2026-10-08; `docs/cluster.md`):**

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

**Movie tooling (branch `movies`; `docs/movies.md`):**
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
  **Look into it and fix it behind a flag (frame rotation) before any
  science run.**
- Finding: frames written before the first RT step have an uninitialised
  rt_grid (fHI = rhog = 0). The tool skips them; a C++ fix that skips them at
  write time is optional.

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

## Framework spike (Stage C), branch `framework-spike`

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
  higher. **Pending the user's decision** and the Hyak runs
  (`spike/hyak/README.md`: single-node layouts, 1–4 nodes).
- GPU: the shared kernels and the Parthenon prototype compile with CUDA
  (nvcc 12.0, compile only). AMReX needs nvcc ≥ 12.2: build on Hyak
  `ckpt` with `spike/hyak/build_gpu.sh`. Nothing has been run on a GPU yet.

## Branches and PRs (survey, 2026-10-09; all merged the same day)

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

**Facts from these branches that change the plan:**
- One klone node runs 128³ at 0.01c in 3 min and 256³ at 0.1c with a
  movie in 1.5 h (192 threads, RT 94× over one core). Resolution is no
  longer the constraint; 256³ is the production size.
- **c_reduced = 0.01 delays reionisation by Δz ≈ 2–2.5 and lowers τ_es by
  20 %** against 0.1c (same ICs). Use 0.1c; check convergence with one
  0.3c or 1c run at 128³.
- The f_esc = 0.5 calibration on `production-run` was done at 128³ and
  0.01c; it has to be redone at 0.1c (and at 256³, whose HMF has not been
  checked).
- On klone 8 RT stored-reference tests fail by 1 ulp in double sums and by
  rays at the removal threshold; the live gfortran harness built there shows
  the same differences and agrees with the C++ byte for byte. The stored RT
  references are platform-dependent (libm). Decision needed: regenerate
  them per platform, or compare the double sums with a tolerance.
- B16 (`movies` branch): the fixed HEALPix frame leaves a photon deficit
  along ±x/±y from every source, up to 0.25 in x_HII near the front; it is
  the "+" crosses and streaks in the stills. Fix behind a flag (a global
  frame rotation first) before any science run.
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
3. **Framework-spike benchmarks** (`spike/hyak/README.md`): `bench_node`
   (MPI × OpenMP layouts on one node), `bench_scale` on 1–4 nodes, and the
   GPU build on `ckpt` (`build_gpu.sh`, nvcc ≥ 12.2). Copy the CSVs into
   `spike/bench/results/` and update `docs/framework_decision.md`. These
   numbers confirm the AMReX choice and feed the NSF ACCESS decision rule.
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

**Decisions needed from the user:** confirm AMReX (Stage D); production
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
   | Flags | legacy_rt = legacy_physics = 0, deterministic_parallel, rt_substep_carry |
   | Sources | sfr_model = abundance, log ξ_ion = 25.38, f_esc recalibrated at 0.1c |
   | hydro_scheme | user decision pending (`unsplit` fixes the cold-gas overheating; `make_production.py` does not set it yet) |
   | Before the run | B16 fix (item 5), 256³ HMF and f_esc |

4. **User's cluster: UW Hyak `klone`, account `astro`.**
   - Owned: about 1000 `cpu-g2` cores (about 5 nodes of AMD Genoa, 192 cores
     and 1.5 TB per node) with 7.8 TB of memory, a 64-core,
     1 TB `cpu-g2-mem2x` slice, and 40 `compute-bigmem` cores.
   - No owned GPUs: GPUs only via preemptible `ckpt` partitions. No time
     limits. See `docs/ROADMAP.md` (answers section).
   - Access: Claude Code installed in user space, `claude remote-control`
     inside tmux, batch jobs with restarts.
5. **Investigate grid-aligned RT artifacts.** The 64³ stills
   (`docs/figures/stills64_*.png`, `tools/movies/stills.py`) show straight
   neutral streaks and "+"-shaped crosses centred on sources, aligned with the
   grid axes, plus faint radial spokes inside bubbles. A suspect is rays moving
   exactly along grid planes (B10, positions on cell faces). Check this before
   science runs, and see whether it also happens with legacy_rt = 1.
   Possibly related: the 32³ `uniform_1src` baseline shows a 12-fold ring of
   isolated neutral cells inside the ionised sphere
   (`docs/figures/raycomp_slices_uniform_1src.png`).
6. ~~Hydro fix for cold supersonic gas~~ Done (`hydro_scheme = unsplit`,
   branch `hydro-unsplit`; η₁/η₂ confirmed against the Enzo and Nyx
   papers). Open: convergence check of the Sedov radius bias and the
   cold-gas time-step sensitivity at 128³.
7. **Plane-parallel / RT-domain mode** (`docs/plane_parallel.md`).
8. After the production run, follow `docs/ROADMAP.md`:
   - ray compression study: done (`docs/ray_compression.md`); 64³
     production-configuration check done, 128³ job script ready for Hyak;
   - framework spike (AMReX vs Parthenon): done, AMReX chosen
     (`docs/framework_decision.md`); Hyak benchmarks to run
     (`spike/hyak/README.md`);
   - ray data-structure rewrite;
   - MPI;
   - AMR for gas, then RT;
   - multigrid/tree gravity;
   - adaptive RT.

## Open questions for the user

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
