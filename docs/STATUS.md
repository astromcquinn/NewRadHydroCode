# Project status

*Last updated: 2026-10-26: the Fortran split for a public repository (below). Before that, 2026-10-23: MUSIC installed on klone (below). Before that, 2026-10-22: the uniform 512³ suite of every RT method run and reported (`docs/rt_suite_512.md`; the Sheth–Tormen sources, f6 = 0.11, the continuous ĉ schedule `c_reduced_schedule_interp`; branch `rt-suite-512`, below). Before that, 2026-10-22: the refined 512³ run against its control (branch `amr512-analysis`, below); the zoom-in turbulence plan and a PR audit (below). Before that, 2026-10-18: the substep-averaged Γ_HI (`gammaHI_avg`), the launch script of the uniform 512³ suite (`tools/hyak/rt_suite512.sh`) and the suite's decisions; before that the assessment and plan for the unified RT-methods suite (`docs/PLAN.md`, "Unified RT-methods suite", with the thread prompts; branch `claude/sleepy-feynman-dx6ic2`, below). Before that, 2026-10-17: ĉ by ionisation state, `c_reduced_schedule` (branch `c-schedule`, below); a falling escape fraction (branch `fesc-evolution`, PR #45, below). Before that, 2026-10-14 (evening): Stage F3, the 8 sources of `uniform_8src` with one and two refined levels, three fixes of the RT across levels (branch `stage-f3-8src`, below). Before that, 2026-10-15: Stage F3, the RT across refinement levels (branch `stage-f3`, below). Before that, 2026-10-14: Stage F decisions, B17 under `legacy_physics = 0`, Stage F1 (the gas hierarchy) and F2 (refined hydro validated, `cool_suppress_overdensity`, the star-formation sink), `srun --mpi=pmix` (branch `stage-f-decisions`, PR #40, below). Also 2026-10-13: Stage E2 on klone, the PM backend validated and scaled (PR #38, below). Before that, 2026-10-13: Stage E4, load balancing (`amrex_balance = 1`) and production runs on the backend (`make_production.py --backend amrex`, `reion_production.slurm` over several nodes), 512³ on 5 nodes, the SIGTERM backtraces (branch `stage-e-balance`, PR #37, below); Stage E3, the state on the levels, `amrex_state = 1` (PR #34 merged, `docs/stage_e_state.md`, below): 512³ runs on one node, bitwise identical to the 2-node host-state run; PR #32 (Stage E1 + the klone runs) and PR #33 (Stage E2) merged into `main`. Also 2026-10-13: OTVET closures merged (PR #35); faster OTVET tensor (branch `otvet-fft-speed`); OTVET closures implemented behind flags (`otvet_closure`, `otvet_split`) and validated at 256³ and 512³: the split and the cone remove OTVET's after-overlap P_Γ excess, and the cone also matches the rays' overlap era best; 512³ runs done (branch `otvet-closure`, PR #35, below, `docs/otvet_closure.md`). Also 2026-10-13: Stage E3, the state on the levels, `amrex_state = 1` (branch `stage-e-state`, `docs/stage_e_state.md`, below): 512³ runs on one node; PR #32 (Stage E1 + the klone runs) and PR #33 (Stage E2) merged into `main`. Before that, 2026-10-12 (evening): Stage E2, the PM gravity on AMReX behind `pm_backend = amrex` (branch `stage-e-pm`, below, `docs/stage_e_pm.md`); Stage E prerequisites on klone, the AMReX backend built, tested (81/81) and timed at 256³ and 512³ on 1–4 nodes, and the multi-node hydro slowdown found to be the UCX transport (`UCX_RNDV_THRESH=8m UCX_ZCOPY_THRESH=inf` in `env.sh`, also 10–15 % off the RT), box size 32³ (branch `claude/clever-fermi-aa4km4`, `docs/stage_d.md` §4.4, `docs/stage_e_hydro.md` §4.1/§5, below); Stage E1, the unsplit hydro on AMReX, `hydro_backend = amrex`, bitwise identical to the serial step (`docs/stage_e_hydro.md`); M1's neutral source cells, `moment_inject` (branch `m1-clumps`, below); what's next (Stage E, MPI, then AMR), the production default back to case B, the Hyak module-load fix (branch `claude/clever-fermi-aa4km4`, below); `pm_refine` (PM mesh finer than the gas grid), FoF cadence option and the Sheth–Tormen abundance matching (branch `pm-refine`, below); the nightly workflow fixed (branch `claude/new-session-knktw3`, below). Also 2026-10-12: rays vs M1 vs OTVET at 40 Mpc/h, 512³ with case A + recombination radiation, and the Γ_HI comparison (branch `rt-methods-512`, `docs/rt_methods_512.md`). Also 2026-10-12: the stored RT references that differ on klone are skipped there (label `platform_ref`); OTVET kept as implemented. Also 2026-10-12: rays crossing k cells per substep, `rt_cells_per_substep` with `photon_conserve` (B17) and `ray_deposit = exposure` (branch `rt-multicell`, PR #27, `docs/rt_multicell.md`, merged into `main`). Before that, 2026-10-11: production runs use `hydro_scheme = unsplit` and `redshift_per_distance = 1` (user decisions, now the `make_production.py` defaults); Stage D, the ray-tracing RT on AMReX (branch `stage-d-amrex`, PR #24, `docs/stage_d.md`, merged into `main`); the 256³ production run (branch `production-256`, PR #23, below) and slab-only movie frames `movie_slab`, the movie-reader speed fix and `POST=sbatch` (branch `movie-slab`, merged into `main`). Before that, 2026-10-10: `redshift_per_distance` (branch `redshift-per-distance`, merged); IGM diagnostics (branch `igm-diagnostics`, merged); recombination radiation with the moment RT and helium lines (branch `recomb-moments`, merged into `main`); status assessment and the Hyak production and Stage D prompts ("Threads to launch"); B16 fix `ray_rotate` merged into `main` (branch `claude/exciting-clarke-zn62ih`; production uses `ray_rotate = 3`). The same day: the framework-spike Hyak CPU and GPU benchmarks (PR #19, branch `spike-hyak-bench`: AMReX approved, no ACCESS proposal for now), case A recombination and recombination radiation (branch `recomb-radiation`), and moment RT (M1, OTVET) with its 256³ klone comparison (branch `claude/elegant-bell-rqd3l1`). Before that, 2026-10-09: the survey (PR #16), `hyak-setup` (PR #17), `movies` (PR #12) and `production-run` (PR #7).*

## Done 2026-10-22: the uniform 512³ suite on klone (branch `rt-suite-512`)

The launch of `tools/hyak/rt_suite512.sh` (the Thread D prompt of
`docs/PLAN.md`): base rays (case B), the case A + RR rays, M1 and the four
OTVET variants, one cpu-g2 node each. Report: `docs/rt_suite_512.md` (when
the set is done).

- **Disk before the launch** (`du -sh`, 2026-10-19): the 2026-10-11 set
  without its checkpoints (deleted earlier): `reion512_rays` 84 GB,
  `_caseB` 80, `_m1` 80, `_otvet` 83, `_otvet_split` 137, `_otvet_screen`
  137, `_otvet_cone` 137 (738 GB; the six names the suite overwrites
  hold 658 GB; `reion512_rays_caseB` is not overwritten and stays).
  `/gscratch/astro` was at 97 % (38.8 of 39.9 TB, 1.1 TB free), so the
  suite's ≈ 1.3 TB while running (one checkpoint kept per run) fits only
  with the old set removed and leaves ≈ 0.4 TB for the group.
- **Build** (2026-10-18, `tools/hyak/build_hyak.sh` on a cpu-g2 node,
  gcc 13.2.0): `ctest -R "gamma_avg|c_schedule|movie_slab"` 6/6
  (`gamma_avg_rays`, `gamma_avg_m1`, `c_schedule`,
  `restart_smoke_c_schedule`, `movie_slab`, `restart_smoke_det_movie_slab`).
- **Probe** (`rt_suite512.sh probe`, job 42013527): the old
  `reion512_rays` removed; three steps of the base at the 2048³ PM mesh
  for the peak memory and the PM time per step: **peak RSS 311 GB**
  (326178764 kB; the estimate was ≈ 275 GB of PM arrays plus the fields,
  no RT yet), **PM 220–235 s per step** (hydro 6–13 s, summaries 5.5 s;
  13.7 min for the three steps on node n3440, 192 threads). Both are under
  the stop limits (700 GB, 10 min); the PM adds ≈ 4 min × ≈ 55 steps ≈ 3.7 h
  to each run against ≈ 1 h with `pm_refine = 1` at 512³.
- **Base** `reion512_rays` submitted 2026-10-18 20:55 (job 42014611, one
  cpu-g2 node, 192 threads, 1000 GB, 4 days); started 20:51 on n3440: 23
  hydro steps of PM alone (≈ 240 s each) to z = 15 in 1.6 h, the first
  sources (10191 halos) at z = 13.7 in step 24, steps of 380–760 s with
  the RT at z = 13.7–12.4.
- **Siblings** submitted 2026-10-18 23:17 after the base had run four RT
  steps with sources (f6 is fixed, so the base's check is a report, not a
  gate): the old `_m1`, `_otvet`, `_otvet_{split,screen,cone}` removed
  (80–88 GB each by then); `reion512_rays_caseA` 42020712, `_m1`
  42020713, `_otvet` 42020714, `_otvet_split` 42020715, `_otvet_screen`
  42020716, `_otvet_cone` 42020717. Each sibling's `params.txt` differs
  from the base's only in `rt_method` and the closure keys, or the
  recombination keys (checked with diff). The old `reion512_rays_caseB`
  (80 GB) is not part of the suite and stays.
- **STOPPED 2026-10-19 08:30: the suite's sources are 20× too bright.**
  The base reads x_HII,V = 0.29 at z = 12.1, 0.83 at z = 10.8 (the ĉ
  schedule reached 0.5c at step 34) and 0.99999 at z = 9.2 with
  Γ_HI ≈ 5.6×10⁻¹² s⁻¹, against z50 = 7.1 and z99 = 5.5 of the
  calibration run `fesc512_b2.5`. Cause (from `step_NNNN/stats.txt`,
  `Luv` = the source luminosity density): with 1024³ particles the
  20-particle FoF floor (`fof_min_npart = 20`) is 1.3×10⁸ M☉ instead of
  the 1.0×10⁹ M☉ of the 128³/512³ calibration, and the fitted
  L_UV(M, z) of `am_params_128_10mpc.txt` (L ∝ M^1.19 below the break,
  with the +0.4–0.9 dex high-z normalisation) is applied to those far more
  numerous halos: 10191 sources at z = 13.7 against 15 at z = 13.2 in
  `fesc512_b2.5`.

  | z | `reion512_rays` (this suite) Luv | `fesc512_b2.5` (512³ particles) | `fesc512_b2.5_uniform` (1024³, ST L(M)) |
  |---|---|---|---|
  | 13.7 | 2.1×10²⁶ | 3.5×10²⁴ (z = 13.2) | 3.0×10²⁵ |
  | 12.1 | 5.8×10²⁶ | — | 4.8×10²⁵ |
  | 11.5 | 6.8×10²⁶ | 3.1×10²⁵ (z = 11.6) | 9.1×10²⁵ |
  | 11.0 | 7.5×10²⁶ | — | 1.4×10²⁶ |
  | 10.5 | 7.0×10²⁶ | 7.1×10²⁵ | 9.7×10²⁵ |
  | 9.4 | 5.0×10²⁶ | 8.8×10²⁵ (z = 9.5) | 9.5×10²⁵ |

  (erg s⁻¹ Hz⁻¹ per comoving volume unit of `stats.txt`; the calibration
  run peaks at 2.1×10²⁶ at z = 6.9.) The decision of 2026-10-18 ("with
  1024³ particles the halo counts shift a little") does not hold: the
  emissivity at z = 11–13 is 20–60× the calibration's, and the RT itself
  is consistent with it (the ionised volume per photon matches the other
  runs). The six siblings carry the same sources. The runs were left to
  the user to cancel (the session's `scancel` was blocked): base
  42014611, siblings 42020712–42020717. Options for the relaunch (user
  decision): a halo-mass floor equal to the calibration's, `m_min =
  1.02e9` M☉ (or `fof_min_npart = 160`; neither is a `make_production.py`
  option yet, `fof_min_npart = 20` is written by it), which keeps the
  512³-calibrated luminosities literally; or a refit of L(M) to the 1024³
  FoF mass function (`am_fit`, the route the Sheth–Tormen table took,
  whose uniform run ended at z = 5.27 with Γ at half the data at 0.1c
  with the sink).
- **User decision 2026-10-19 (morning): the Sheth–Tormen table**
  `data/uv_lf/am_params_st_fof.txt` with the 1024³ particles and
  `fof_min_npart = 20`, the source population of `fesc512_b2.5_uniform`
  and `fesc512_b2.5_amr` (the 2026-10-11 decision to calibrate to an
  analytic mass function); no precheck, f6 = 0.135 and β = 2.5 kept. The
  seven jobs cancelled at 09:05 after 9–11.5 h; the base relaunched with
  `AM_PARAMS=data/uv_lf/am_params_st_fof.txt OVERWRITE=1 rt_suite512.sh
  base` (job 42035621; **every later phase needs the same `AM_PARAMS`**,
  the script's default is still the 128-box table). Plan: let the base
  reach x_V ≈ 0.3 and compare its x_V(z) with the two earlier runs
  (0.1c, the sink, cooling suppression 0, `ray_deposit = exposure`,
  `fof_dt_myr = 30`, the AMReX backend; identical sources otherwise):
  they had x_V = 0.023 at z = 12.05, 0.059–0.062 at 11.05, 0.113–0.118 at
  10.05, 0.170–0.174 at 9.05, 0.204–0.206 at 8.55, 0.263–0.265 at 8.05,
  0.35–0.36 at 7.55, 0.48–0.49 at 7.05; the schedule stays at 0.1c until
  x_V = 0.25, so the base should track them to z ≈ 8.1 and then run
  ahead. Then the siblings.
- **The base's check, 2026-10-19 12:45 (z = 9.67, 4.2 h, still at
  0.1c):** the sources are the reference's (`stats.txt`: at the steps
  where `fesc512_b2.5_uniform` refreshed its FoF catalogue, N and Luv
  agree to 3 digits: 10191 / 2.964×10²⁵ at z = 13.69, 79664 / 1.412×10²⁶
  at 11.00, 171906 / 1.354×10²⁶ at 9.73; between refreshes the reference
  lagged, e.g. 21601 against 36516 sources at z = 12.11). The base runs
  ≈ 1.3× more ionised at fixed z from z = 13 on (x_V 0.0107 / 0.0085 at
  z = 13, 0.0563 / 0.0419 at 11.5, 0.116 / 0.086 at 10.5, 0.172 / 0.133 at
  9.7 against the uniform run; the AMR run is within 0.005 of the
  uniform one): the reference's 30 Myr FoF cadence (fewer photons), the
  sink and cooling suppression 0 against none and 200 (more
  recombinations in dense gas there), and `ray_deposit = exposure`. Same
  shape, a constant ratio, the sources verified: taken as the rough match
  asked for.
- **Siblings relaunched 2026-10-19 12:48** with the same `AM_PARAMS`
  (the first launch's partial outputs removed, 204 + 104 + 102 GB):
  `reion512_rays_caseA` 42044005, `_m1` 42044006, `_otvet` 42044007,
  `_otvet_split` 42044008, `_otvet_screen` 42044011, `_otvet_cone`
  42044012; each `params.txt` differs from the base's only in the method
  or recombination keys (diff).
- **Stopped again 2026-10-19 15:30 (user): tune f6 so that reionization
  ends at z = 5.5.** The base had reached z = 8.02, x_V = 0.345 (7 h; the
  schedule at 0.2c from z = 8.59); the siblings 2–2.6 h. Two estimates of
  the base's end disagreed: the lead over `fesc512_b2.5_uniform` at
  fixed x_V (Δz = 0.41–0.47 from x_V = 0.05 to 0.3) plus the schedule's
  faster end (≈ 0.4 from the ĉ study) gave z99 ≈ 6.1; a one-zone model
  (dQ/dt = ṅ/n̄_H − C α_B n̄_e Q with the run's own emissivity, the
  reference's refresh-step `Luv` below z = 8, C = 4.75 fitted to the
  base's x_V(z) to 25 %; `scratchpad/onezone.py`) gave z50 = 7.13,
  z90 = 5.76, z99 = 5.53 at f6 = 0.135, and 5.28 / 5.10 / 4.92 at f6 =
  0.12 / 0.11 / 0.10 (β = 2.5), 5.46 / 5.42 at β = 2.0 / 1.5 with f6 =
  0.135. So the calibration is measured instead: **`rt_suite512.sh
  precheck`** with `AM_PARAMS=data/uv_lf/am_params_st_fof.txt` and
  `PRECHECK_F6="0.06 0.08 0.1 0.135"` (0.135 anchors the 20 → 40 Mpc/h
  offset against the cancelled base's history to z = 8): the ICs job
  42051860 (256³ gas, 512³ particles, 20 Mpc/h, seed 7), the runs
  `$RUNS/rt512cal256_f{0.06,0.08,0.1,0.135}` 42051864, 42051915–42051917
  (one node and ≈ 3 h each, the suite's configuration at the suite's
  cell size and particle mass, `pm_refine = 2`). Then `precheck-post`;
  pick the f6 whose z99 at 256³ in 20 Mpc/h is ≈ 5.3 (the 40 Mpc/h box
  ends ≈ 0.2 later, `docs/fesc_evolution.md` §4), β = 2.5 kept; then the
  base with `F6=<that>` and the siblings.
- **Precheck results (2026-10-19 20:46; ICs 3 min, the runs 1.6–2.9 h
  each):**

  | f6 (β = 2.5) | z10 | z25 | z50 | z75 | z90 | z99 | z99.9 | τ_es | x_V(z = 4.95) |
  |---|---|---|---|---|---|---|---|---|---|
  | 0.135 | 10.30 | 8.33 | 7.37 | 6.80 | 6.51 | 6.21 | 6.09 | 0.0614 | 1 |
  | 0.1 | 9.55 | 7.67 | 6.72 | 6.03 | 5.67 | 5.36 | 5.26 | | 1 |
  | 0.08 | 8.67 | 7.18 | 6.09 | 5.33 | 5.05 | — | — | | 0.944 |
  | 0.06 | 7.89 | 6.42 | 5.13 | — | — | — | — | | 0.574 |

  z99 moves by 2.8 per e-fold of f6 between 0.1 and 0.135 (the earlier
  "0.07 per 11 %" reading was wrong; the one-zone model, which gave
  z99 = 5.53 at f6 = 0.135, is 0.6–0.7 too early and is dropped; the
  lead-based estimate of ≈ 6.1 for the 40 Mpc/h base was right). The
  anchor: the 40 Mpc/h base at f6 = 0.135 was *ahead* of the 20 Mpc/h
  256³ run early (z at x_V = 0.05, 0.1, 0.2, 0.3, 0.34: 11.64 / 11.28,
  10.74 / 10.30, 9.29 / 8.79, 8.28 / 8.07, 8.05 / 7.93), the gap
  shrinking from 0.4–0.5 to 0.1 by x_V = 0.34, unlike the 128-table
  calibration where the big box was later throughout (0.3 at z10, 0.2 at
  z99). Taking the end offset as ≈ +0.1 (the big box a little later), the
  target is z99 = 5.6 at 20 Mpc/h, i.e. **f6 = 0.11** by interpolation in
  ln f6 (5.36 + 2.8 ln(0.11/0.1) = 5.63); the expected 40 Mpc/h end is
  z99 = 5.5 ± 0.15. **Base relaunched 2026-10-19 20:52 with `F6=0.11`**
  (job 42061340; the cancelled f6 = 0.135 base removed). `precheck-post`
  running (igm_diag 42061305, summary 42061306:
  `$RUNS/rt512cal256/igm/summary.md`, done 21:05): the IGM of the grid
  against the data (ratio in brackets):

  | f6 | τ_es | ṅ(5.4) [log] | Γ(5.4) | Γ(5.1) | ⟨F⟩(5.4) | ⟨F⟩(5.1) | λ(5.4) | λ(5.1) |
  |---|---|---|---|---|---|---|---|---|
  | 0.1 | 0.0548 (+0.1σ) | 50.74 | 0.190 (0.51×, x_V = 0.98) | 0.403 (0.79×) | 0.043 (0.5×) | 0.168 (1.2×) | 5.9 (1.2×) | 21.0 (2.3×) |
  | 0.135 | 0.0614 (+1.1σ) | 50.87 | 1.072 (2.9×) | 1.339 (2.6×) | 0.265 (3.3×) | 0.353 (2.5×) | 39.2 (8.1×) | 60.5 (6.7×) |

  (f6 = 0.08 and 0.06 are unfinished at z = 5.1: Γ 0.16× and 0.07×.)
  f6 = 0.11 sits between: at z = 5.1 Γ ≈ 1.0–1.2× and ⟨F⟩ ≈ 1.4× the
  data, ṅ_ion ≈ 0.4 dex below Becker+21 as in `docs/fesc_evolution.md`;
  λ_mfp stays 2–3× the data after overlap (the known resolution excess,
  §5–6 there).
- **The f6 = 0.11 base's check (2026-10-20 04:05, z = 7.75, x_V = 0.300,
  7.2 h; 0.2c from z = 7.94):** z at x_V = 0.05, 0.1, 0.2, 0.25, 0.3 is
  11.32, 10.28, 8.42, 7.98, 7.75 against 10.94, 9.79, 8.26, 7.88, 7.66
  interpolated (in ln f6) between the 256³ f6 = 0.1 and 0.135 runs: the
  40 Mpc/h box ahead by 0.4–0.5 early, 0.1 by x_V = 0.3, as the 0.135
  anchor was. On track for z99 ≈ 5.5. **Siblings launched 04:08** from
  this base (the cancelled second launch's partial outputs removed, 73 +
  90 + 88 + 74 GB): `reion512_rays_caseA` 42070978, `_m1` 42070979,
  `_otvet` 42070980, `_otvet_split` 42070981, `_otvet_screen` 42070982,
  `_otvet_cone` 42070983; `params.txt` checked (f_esc = 0.11, the ST
  table; the method or recombination keys the only differences).
- **Base done 2026-10-21 01:00** (job 42061340, n3476; 28 h 10 min wall,
  99 hydro steps, 5225 RT substeps, peak RSS 475 GB; `DONE` at z = 4.99;
  the ray pool guard never fired; 281 GB with the last checkpoint).
  Cost (`timing.txt`): RT 16.3 h (ray_tracing_1 6.4 h + ray_tracing_3
  5.8 h + physics 3.1 h), PM 6.3 h (the 2048³ mesh, 230 s per step), I/O
  4.4 h, halos 0.6 h, P(k) 0.5 h, hydro 0.2 h. Against the 2026-10-11
  case-B rays (8.6 h): the schedule's 0.5c end, the PM and twice the
  steps. **History** (`reion_history.py`, volume-weighted): z10 = 10.28,
  z25 = 7.98, **z50 = 6.92**, z75 = 6.27, **z90 = 5.95**, **z99 = 5.59**,
  z99.9 = 5.46, **τ_es = 0.059** (+0.7σ); mass-weighted z50 = 7.33, z99 =
  5.62. The end is 0.09 above the z = 5.5 target (the 256³ interpolation
  put the 20 Mpc/h box at 5.63; the 40 Mpc/h box ended 0.04 earlier, not
  0.1 later). Against `fesc512_b2.5` (z50 = 7.11, z99 = 5.51, τ = 0.056):
  a later start and midpoint, the same end to 0.1. IGM statistics: job
  42100957 (`$RUNS/rt_methods_512/igm_base/`). The automatic post job
  42100693 (stills, movie).
- **The checkout moved (2026-10-21 01:10):** another session checked out
  `split-on-entry` in `NewRadHydroCode/`. The suite's jobs resolve
  `build/radhydro` there at start, so the suite's binary (built
  2026-10-18 20:28 from `rt-suite-512`) was copied to
  `$RUNS/rt_suite512_bin/radhydro` (with `COMMIT`, `README`) and the
  three queued closure jobs were resubmitted with `EXE=` that path:
  `_otvet_split` 42101054, `_otvet_screen` 42101055, `_otvet_cone`
  42101056 (the running three keep their loaded binary; a resubmission
  after a failure must use the copy). This thread's docs work continues
  in the worktree `$RADHYDRO/worktrees/rt-suite-512`.
- **The base's IGM (2026-10-21, `igm_diag.slurm` 22 min + `fesc_summary.py`;
  `$RUNS/rt_methods_512/igm_base/`), next to `fesc512_b2.5` (128-box
  table, 512³ particles, 0.2c, case A + RR; `docs/fesc_evolution.md` §4);
  ratio to the data in brackets; Γ is the mean over x_HI < 10⁻² as in that
  table, the PR #48 median estimator (G_medx, the neutral volume as zero)
  in the last column:**

  | z | run | x_V | Γ_HI [10⁻¹²] | λ_mfp [pMpc] | ⟨F⟩ | A | Γ_medx/obs |
  |---|---|---|---|---|---|---|---|
  | 5.9 | `reion512_rays` (f6 = 0.11, this suite) | 0.918 | 0.205 (1.36×) | 3.9 (3.6×) | 0.023 (2.0×) | 1.41 | 1.05 (z = 5.98) |
  | | `fesc512_b2.5` | 0.924 | 0.237 (1.6×) | 3.2 (2.9×) | 0.022 (1.9×) | 1.40 | |
  | 5.4 | `reion512_rays` | 0.9999 | 0.462 (1.24×) | 19.0 (3.9×) | 0.149 (1.9×) | 1.82 | 1.10 (z = 5.49) |
  | | `fesc512_b2.5` | 0.997 | 0.371 (1.00×) | 9.4 (1.9×) | 0.090 (1.1×) | 1.12 | |
  | 5.1 | `reion512_rays` | 1 | 0.721 (1.42×) | 37.0 (4.1×) | 0.259 (1.8×) | 2.09 | 1.28–1.45 (z = 5.25–4.99) |
  | | `fesc512_b2.5` | 1 | 0.493 (0.97×) | 16.3 (1.8×) | 0.155 (1.1×) | 1.14 | |

  ṅ_ion = 10^50.73–50.82 at z = 5.9–5.1 (0.4–0.75 dex below Becker+21).
  Reported only (f6 not retuned): Γ is 1.2–1.4× the data by the mean
  estimator (1.05–1.1× by the median one at z ≥ 5.5, rising to 1.45× at
  z = 5), ⟨F⟩ 1.8–2.0× and λ_mfp 3.6–4.1× the data, i.e. about 2× the λ
  excess of `fesc512_b2.5` at the same Γ: the Sheth–Tormen sources with
  the schedule's 0.5c end give a more transparent post-overlap IGM at
  this resolution (the faint 1.3×10⁸ M☉ halos ionise the voids early;
  the 1024³ particles on the 2048³ mesh do not add Lyman-limit opacity,
  the open item of `docs/fesc_evolution.md` §6).
- **M1 done 2026-10-21 04:52** (job 42070979, 14 h 32 min, peak RSS
  367 GB; RT 6.7 h, PM 6.2 h) and **OTVET thin done 05:08** (42070980,
  14 h 45 min, 373 GB; RT 6.8 h). Histories (V-weighted): M1 z50 = 6.93,
  z90 = 5.98, z99 = 5.67, τ = 0.0593; OTVET 6.92 / 5.99 / 5.64 / 0.0594;
  the rays 6.92 / 5.95 / 5.59 / 0.0591: the same midpoint, the moment
  methods end 0.05–0.08 later. Running: the case A + RR rays (z = 5.60 at
  22 h), `_otvet_split` (z = 9.9), `_otvet_screen` (started 05:10);
  `_otvet_cone` queued.
- **The Γ spikes and the continuous schedule (2026-10-21 morning).** The
  median Γ_HI of the ionised cells (and Γ_V) jumps at every switch of the
  ĉ schedule, in the rays and M1 alike: ×1.9 at step 52 (0.1 → 0.2c,
  z = 7.94), ×1.46 at step 63 (0.2 → 0.3c, z = 6.89), ×1.71 at step 74
  (0.3 → 0.5c, z = 6.25), then a relaxation over 3–5 steps: the photon
  reservoir n_γ ≈ ṅ λ/ĉ of the old speed is drained at the new one
  (Γ = σ ĉ n_γ; the steady state σ ṅ λ is ĉ-independent; no photons are
  created). Built: `c_reduced_schedule_interp = 1` (ĉ linear in x_HII,V
  between the plateau midpoints, `docs/production_run.md` §1;
  `make_production.py --c-schedule-interp`; the `c_schedule` test's
  third check; 4/4 tests on klone, commit df3f26a; the binary
  `$RUNS/rt_suite512_bin/radhydro_interp`). **User decision: rerun the
  whole suite with the continuous schedule, overwriting the step-schedule
  runs; the base first, checked, then the six siblings.** The step-schedule
  set's record: the histories and IGM tables above, the partial movie
  `$RUNS/rt_methods_512/compare_partial.mp4` (rays, M1, OTVET, case A + RR
  to z = 5.2) and `igm_base/`; its runs are overwritten. The four running
  jobs (case A + RR at z = 5.2 after 27 h, split, screen, cone) cancelled
  at 09:55; **base relaunched 09:56** (job 42114332; `MAKE_EXTRA="--c-schedule-interp 1"
  EXE=$RUNS/rt_suite512_bin/radhydro_interp`, submitted from the
  worktree so the job's scripts are the branch's). Every later phase:
  the same `AM_PARAMS`, `F6=0.11`, `MAKE_EXTRA`, `EXE`.
- **The continuous base's check (2026-10-21 17:56, z = 7.70, x_V = 0.305,
  5.1 h):** ĉ rose from 0.1c at x_V = 0.125 by ≈ 0.0045c per step to 0.169c
  at x_V = 0.30 (the step schedule would have jumped to 0.2c at 0.25); no
  spike: Γ_V 2.92 → 3.07 → 3.36 ×10⁻¹⁴ over steps 50–52 (the step base:
  2.86 → 5.43); z at x_V = 0.05, 0.1, 0.2, 0.25, 0.3 = 11.32, 10.28, 8.46,
  8.02, 7.72 against 11.32, 10.28, 8.42, 7.98, 7.75 of the step base.
  **Siblings launched 18:00** (the step-schedule runs removed, 985 GB):
  `_caseA` 42132367, `_m1` 42132368, `_otvet` 42132369, `_otvet_split`
  42132370, `_otvet_screen` 42132374, `_otvet_cone` 42132375, all with
  `c_reduced_schedule_interp = 1`, f6 = 0.11 and the `radhydro_interp`
  binary (checked). Another session's `conserve-test` job holds one of
  the account's five nodes.
- **Continuous set, M1 and OTVET thin done 2026-10-22 03:5x** (9 h 47 min
  and 9 h 56 min; peak RSS 367 / 373 GB). Histories: M1 z50 = 6.95, z90 =
  5.97, z99 = 5.66, τ = 0.0594; OTVET 6.94 / 5.98 / 5.64 / 0.0594, within
  0.02 of the step-schedule runs (6.93 / 5.98 / 5.67 and 6.92 / 5.99 /
  5.64): the continuous ĉ changes the history by nothing measurable. Cost:
  the RT the same (6.6–6.8 h) but the PM 1.65 h against 6.2 h in the step
  set (60 against 230 s per step, M1 on the same node n3484 both times):
  the step set's binary (`rt_suite512_bin/radhydro`, built 2026-10-18 in
  the main checkout) was linked against the single-threaded FFTW because
  that checkout's `build/CMakeCache.txt` of 2026-10-12 still named
  `fftw-3.3.10-float` (`FFTW3F_OMP_LIBRARY-NOTFOUND`), so the 2048³ PM
  FFTs ran on one thread; the worktree's fresh build (`radhydro_interp`)
  has `libfftw3f_omp`. The probe's 220–235 s per PM step was this too.
  The continuous set is consistent within itself (one binary); the
  cost table of the report uses it.
- **Continuous set: OTVET split done 2026-10-22 05:00** (10 h 38 min;
  z50 = 6.94, z90 = 5.98, z99 = 5.64, τ = 0.0594; RT 7.5 h) and **the base
  done 12:41** (job 42114332, 23 h 52 min, peak RSS 475 GB; RT 16.4 h, PM
  1.9 h, I/O 4.4 h; `DONE` at z = 4.99). Base history: z10 = 10.28, z25 =
  8.02, **z50 = 6.94, z90 = 5.94, z99 = 5.58**, z99.9 = 5.46, **τ_es =
  0.0591**, within 0.02 of the step-schedule base. Across the former
  switch points Γ_V now moves by 4–5 % per step (7.62 → 7.96 ×10⁻¹⁴ at
  step 63, 1.24 → 1.29 ×10⁻¹³ at step 74) instead of ×1.46 and ×1.71. IGM
  diagnostics: job 42163385 (`$RUNS/rt_methods_512/igm_base_interp/`).
  Running: case A + RR (z = 5.54 at 19 h), screen (z = 5.24), cone
  (z = 5.38). **Its IGM** (`igm_base_interp/`, done 13:05): the same as
  the step base's at every snapshot after z = 5.98 (Γ_medx/obs 1.06, 1.06,
  1.09, 1.28, 1.45 at z = 5.98, 5.70, 5.49, 5.25, 4.99; ⟨F⟩ 1.76–2.03×;
  λ 2.4–4.3×), but at the z = 6.20 snapshot Γ_medx is 0.118 against 0.153
  in the step base, whose snapshot there sat inside the 0.3 → 0.5c spike:
  the continuous schedule removes a 30 % artefact from that snapshot and
  changes nothing after relaxation.
- **OTVET screen done 2026-10-22 13:59** (10 h 08 min, peak 370 GB; z50 =
  6.96, z90 = 6.01, z99 = 5.68, τ = 0.0596) and **OTVET cone done 14:34**
  (10 h 35 min, 376 GB; z50 = 6.94, z90 = 5.96, z99 = 5.60, τ = 0.0594;
  RT 7.4 h of which the cone tensor 0.5 h). **Case A + RR rays done
  18:50** (24 h 51 min, peak 483 GB; z50 = 6.98, z90 = 5.97, z99 = 5.61,
  τ = 0.0595; RT 17.7 h, only 1.04× the case-B rays' wall time against
  1.67× in the 2026-10-11 set: with `recomb_rad_cell_onspot = 1` most
  recombination photons are absorbed on the spot and the pool guard never
  fired). **The continuous set is complete** (7 runs, 1.40 TB with one
  checkpoint each: 281 + 284 GB the ray runs, 167 GB each moment run;
  /gscratch/astro at 96 %, 2.0 TB free). Analysis jobs submitted 19:05:
  `rt_suite512.sh post` (42179001: `rt_methods_runs.py`, the seven-run
  movie, `igm_stats.py` on every run into `$RUNS/rt_methods_512/`),
  `gamma_power.py --field gammaHI_avg` (plain and `--clip 3`), the
  last-substep field for comparison, and `rt_methods_runs.py --gamma-field
  gammaHI_avg` (42179002, `logs/rt512_gamma_avg.slurm`), and the paper
  figures `tools/analysis/rt_suite_figures.py` (42179003: maps of x_HII and
  Γ_HI at x_V = 0.5, 0.8, 0.97 (and 1 for Γ), the 3D P(k) of both fields
  and r(k) against the case-B rays; `figures/` and `figures_clip3/`).
- **Report written 2026-10-22: `docs/rt_suite_512.md`** (configuration,
  the first launch and the table change, the probe, cost, histories, Γ
  maps and PDFs, P(k) and r(k) of x_HII and Γ against the rays, the IGM
  table, the movie `docs/movies/rt_suite_512.mp4`, conclusions and what
  differs from the first set; figures `docs/figures/rt_suite_512_*`).
  Headline: histories method-independent to 0.1 in z (z99 = 5.58–5.68,
  τ = 0.059); M1 the outlier (fronts decorrelated at k ≈ 1 h/Mpc, no
  proximity zones, 25× too little large-scale P_Γ after overlap); OTVET
  thin 11–14× too much P_Γ after overlap, removed by split/screen/cone
  (0.7–1.9× the rays with r = 0.92–0.96); the cone the best match to the
  rays overall; case A + RR within 5 % of case B in Γ and r ≥ 0.99 at
  k ≤ 2 h/Mpc; Γ within 10 % of Gaikwad at z ≥ 5.5 for every method,
  λ_mfp 3.5–5× and ⟨F⟩ 1.9–2.1× the data for all. New tool
  `tools/analysis/rt_suite_figures.py`. Ordering changes from the first
  set: the moment methods now end slightly *later* than the rays (the
  schedule's 0.5c end), and recombination radiation costs 4 % instead of
  67 % (`recomb_rad_cell_onspot = 1`).
- `origin/main` merged into the branch 2026-10-22 (one STATUS conflict);
  the fast suite on klone 145/145 (job 42184776, 38 min). PR opened.
- **2026-10-26, after the repository move:** `origin` re-pointed, the
  branch reset to the re-created `origin/rt-suite-512` (no unpushed
  commits), the new `main` (the Fortran split) merged, the submodule
  `reference/fortran` initialised; PR #2 on the public repository. **The
  live Fortran tests had been broken on klone since 2026-10-23:** the
  conda-forge gcc/gfortran 15 compilers installed into the shared `radray`
  env (for MUSIC) are exported by conda's activation as `CC`/`CXX`/`FC`,
  and the harness Makefiles' `$(FC)` then built PIE objects that cannot
  link the static FFTW ("relocation R_X86_64_32S ... recompile with
  -fPIE"). Fix: `tools/hyak/env.sh` exports the module's gcc 13 `CC`,
  `CXX`, `FC` after the activation (commit 5c4c6c5), **and** the stale
  `reference/fortran/harness*/build` directories must be removed once
  (objects from gfortran 15 with compressed debug sections the gcc 13
  linker rejects: "file not recognized"). With both, the fast suite is
  145/145 (137 in job 42408423 plus the 8 `fortran_live` tests rerun on a
  node, job 42409410). Other checkouts built after 2026-10-23 need the
  same cleanup.

## Done 2026-10-18: `gammaHI_avg` and the 512³ suite's launch script (branch `claude/sleepy-feynman-dx6ic2`)

- **`gammaHI_avg`** (`fields`, `movie_fields`; `docs/igm_diagnostics.md`,
  `docs/movies.md`): Γ_HI averaged over the hydro step's RT substeps,
  weighted by the substep length, accumulated in double after every
  substep's physics (`Simulation::gamma_avg`, `gamma_avg_begin/add/field`
  in `rt.hpp`, hooked in `radiative_transfer`); the step average in the
  snapshots, the running average in frames written inside a step. Output
  only (every other file byte-identical), thread-independent, serial RT
  only (rays and the moment methods; `rt_backend = amrex` refuses it).
  `movie_every_substep = 1` writes a frame after every RT substep (tests).
  Tools: `gamma_power.py --field gammaHI_avg`, `rt_methods_runs.py
  --gamma-field gammaHI_avg`. Tests `gamma_avg_rays`, `gamma_avg_m1`
  (`tests/gamma_avg_check.py`: the mean of the per-substep frames, the
  step field against the last frame, 1 against 2 threads, with against
  without the field).
- **`tools/hyak/rt_suite512.sh`**: the uniform 512³ rerun of every RT
  method with the user's decisions of 2026-10-18 (`docs/PLAN.md`): phases
  `probe` (a few steps at the 2048³ PM mesh, peak memory and PM time),
  `base` (the ray run, case B), `siblings` (the case A + RR rays, M1,
  OTVET thin, split, screen, cone, through the existing sibling scripts)
  and `post`; one checkpoint kept (`RESTART_KEEP`);
  `OVERWRITE=1` removes the old directories of the same names; never
  reuses an old `params.txt`; the disk estimate in its header.
- Decisions recorded in `docs/PLAN.md`: the 512³-calibrated luminosities
  kept, cooling suppression 200 and no sink (a serial-path sink later),
  overwrite, no dark matter in the outputs, case B as the suite's default
  with one case A + RR ray run; the disk estimate.

## 2026-10-17: the unified RT-methods suite, assessed (branch `claude/sleepy-feynman-dx6ic2`)

**`docs/PLAN.md`, "Unified RT-methods suite".** The user's target: every
RT method at the same motivated ĉ, 512³, the same f_esc, ideally with
AMR, plus flux-limited diffusion. Against it: the moment methods exist
only in the serial code (no MPI, no AMR); the AMReX backend has no
recombination radiation, so the AMR production runs are case B while the
method comparison was case A + RR; the 512³ sets so far differ in ĉ (0.2
against 0.1) and f_esc (constant 0.23 against the falling law of
`fesc-evolution`, since confirmed at 512³ and merged, PR #45); the ĉ study of
`docs/production_run.md` §9 covers the history but not the post-overlap
Γ; FLD does not exist and is not an isotropic Eddington tensor (the flux
is slaved to ∇N; it needs an implicit solver). **Decided the same day
(user): the default ĉ is a schedule on x_HII,V**, implemented by another
thread as `c_reduced_schedule` (PR #46) and revised on 2026-10-18 to 0.1c
below 0.25, 0.2c to 0.5, 0.3c to 0.75, 0.5c after. Proposed order:
the schedule's check against full c (thread A, two runs; #42, #44, #45
and #46 are in `main` as of 2026-10-18, and #43, AMR gravity, is not
needed for the uniform suite); FLD in
the serial moment framework (thread B) and the uniform 512³ suite of all
methods at one ĉ, f_esc and recombination treatment; then the moment
methods, FLD and `recomb_rad` on the AMReX hierarchy (thread C) and the
suite with two levels. Five decisions are listed for the user there.

## Done 2026-10-26: the Fortran split for a public repository (branch `fortran-split`)

Hy Trac agreed to this repository being public with his code removed.
`docs/fortran_split.md` is the plan and the record. The tree: his original
(`reference/Radhydro`), the harness sources, Makefiles, run scripts and
`PATCHES.md`, `harness_common`, and the two Fortran test helpers moved to
the **private** repository `astromcquinn/RadHydro-fortran` (with their
history, 20 commits), carried here as the submodule `reference/fortran`;
the stored `outputs*` stay. CMake: `RADRAY_FORTRAN_DIR` and
`RADRAY_HAVE_FORTRAN_REF`; the live gfortran comparisons exist only with
the submodule, the stored-reference comparisons always. The workflows
populate the submodule from the Actions secret `FORTRAN_TOKEN`. CLAUDE.md
rule 1, README, `reference/README.md` and the path mentions updated.
Verified in the cloud container: the fast suite with the reference absent
and the `fortran` tests with it present. **Left for the user:** create the
private repository and the token, merge, then the history rewrite
(`tools/git/scrub_fortran_history.sh`, §3 of the plan: a new public
repository from the scrubbed history, recommended, or an in-place
force-push) at a quiet moment after the thread PRs merge.

## Done 2026-10-23: MUSIC on klone (user)

`tools/hyak/setup_music.sh` ran on a klone login node: `$RADRAY_SW/music/MUSIC`
built (cosmo-sims/MUSIC, the radray conda env with gsl, hdf5, the double
FFTW and camb), `tools/ics/music/test32.conf` ran, and the CAMB table of
`camb_transfer.py` (the Cain et al. 2026 cosmology at z = 300) is at
`$RUNS/ics_test/camb_z300.dat`. Thread D of `docs/zoom_turbulence.md` §6
(the MUSIC converter) is unblocked; the configs for the C26 box and the
Z3 zoom are in `tools/ics/music/`.

## Done 2026-10-22: zoom-in turbulence plan; PR audit (branch `claude/epic-cannon-5laltb`)

* **Plan:** `docs/zoom_turbulence.md`: refine one region of the 40 h⁻¹ Mpc
  base to ≥ 1024³ cells at 0.3–1.2 kpc/h (C26's requirement for
  reionization-driven turbulence), with hierarchical zoom ICs, levels
  tagged by particle mass class, multigrid gravity on the hierarchy
  (Stage G, now a prerequisite: the uniform `pm_refine` mesh stops at
  20 kpc/h), rays into the zoom (the HEALPix cap `lside_max = 13` blocks
  splitting below r/8192, 1.2 kpc/h at 10 Mpc/h), a preheating floor,
  diagnostics, the ladder Z1–Z4, cost estimates (Z3 ≈ 10⁵ core-hours,
  days on the owned klone cores; the suite a few 10⁶ for a proposal) and
  MHD afterwards (kinematic induction with constrained transport first).
  Decisions needed from the user: its §9. Request recorded in
  `docs/PLAN.md`.
* **Open PRs and branches (2026-10-22, against `main` at 1132b74):**
  - PR #47 `moments-amrex` (moment RT on the backend, 25 commits, +4243):
    no CI check recorded; conflicts with `main` in `apps/radhydro.cpp`,
    `config.hpp`, `rt_backend.hpp`, `src/amrex/rt_level.cpp`,
    `rt_level_stub.cpp`, `src/rt_backend.cpp` (PRs #48 and #50 landed after
    its last merge of `main`). Needs a conflict resolution and a CI run
    before it can merge.
  - PR #49 `claude/sleepy-feynman-dx6ic2` (gammaHI_avg, `rt_suite512.sh`,
    the unified-suite plan; CI green): conflicts only in `docs/PLAN.md`.
  - Branch `rt-suite-512` (no PR, 33 commits, contains PR #49 plus
    `c_reduced_schedule_interp`, `docs/rt_suite_512.md`, the suite's
    launch and status): conflicts only in `docs/PLAN.md`. It supersedes
    PR #49: one PR from `rt-suite-512` (and closing #49) is the simplest.
  - Branch `claude/kind-wozniak-4fiw74` (docs only, 2026-10-11, the
    "thread 10" planning that `rt-methods-512` carried out): stale, can be
    deleted.
  - Every other remote branch is merged.

## Done 2026-10-22: two refined levels rank-independent (branch `rank-independent-deposits`)

The open item of 2026-10-16 is fixed. Two refined levels differed between
rank counts from about step 8 (relative 1e-10). It was a dependence on the
fine levels' box layout, which AMReX chops by the rank count (level 1 of
the 8-source test: 1, 2, 4 boxes for 1, 2, 3 ranks; on one rank, fine boxes
of 16 instead of 32 cells also changed the result from step 6): the
deposits a fine level's rays leave in ghost cells outside the level were
summed per fine box, and each box's partial sum was divided (by 8, or the
`conserve` weight) and rounded before the boxes' sums were added in the
parent. Now the raw int64 sums are added in the parent and divided or
weighted once per parent cell. 1 rank (fine boxes 32 or 16) and 2 and 3
ranks give byte-identical outputs (all 51 files of the 8-source two-level
run). New nightly CTest `amrex_rt_amr_layout` (16³, two nested levels; 1 rank
with two box sizes, 2 and 3 ranks; ~18 min). `ctest -L amrex -LE nightly`
183/183 on klone. Refined results change in the last digits.

## Done 2026-10-22: `fesc512_b2.5_amr` against `fesc512_b2.5_uniform` (branch `amr512-analysis`)

`docs/stage_f.md` §4 (results, cost, figures `docs/figures/amr512_*`,
movie `docs/movies/amr512.mp4`); the cost also in `docs/zoom_turbulence.md`
§7 (on the zoom-plan branch `claude/epic-cannon-5laltb`, PR #52).
Analysis jobs `igm_diag.slurm` (13 min) and `fesc_post.slurm` (13 min;
now takes `REF`/`REFLAB`/`REFTITLE`/`NEWTITLE`, `SPREAD=` for no third
run); `fesc_summary.py` for the tables.
* **IGM:** two levels above 8 ρ̄ / 64 ρ̄ (19.5 kpc/h cells, ≥ 6 % more
  cells at z = 5) end reionization 0.1 later (z99 5.16 against 5.27),
  lower Γ_HI by 10–17 % and λ_mfp by 20–33 % at fixed z (≈ 10–20 % at
  fixed Γ), raise the relative P_Γ by 1.3–1.8× at k ≤ 5 h/Mpc with the
  same pattern, and keep self-shielded specks in the ionized gas. The sink
  takes 1.0 % of the gas (0.012 % unrefined).
* **Cost:** 2.6× the step time (11.3 h against 4.3 h on 5 nodes; 3–4× per
  step after z = 6): the RT 3.7× (2.9× the rays, born at the finest level
  where the sources are), the hydro 46× (FillPatch, refluxing and the CFL
  reduction: per-box overhead of small scattered boxes on 119 ranks), a
  typical rank 39 GB against 17 GB. Checkpoints at `restart_every = 1`
  cost 760 s each (290 s unrefined) and were 62 % of the last job's 23.2 h.
  Wall-clock of the refined run 37.5 h over three jobs (two crashes, a
  9.5 h hang); core-hours 10 900 (steps) / 36 000 (wall) against
  4 100 / 12 900.
* **Open:** the λ excess at fixed Γ is only partly removed by this
  refinement (`docs/fesc_evolution.md` §6 remains to discuss); the
  per-box hydro overhead (a larger `amrex_max_grid_size` / blocking
  factor, or fewer ranks for the fine levels) before the next refined
  production run; the refined cell count per level is not logged
  (add it to the timers line).

## Done 2026-10-21: rays at the refinement boundaries (branch `split-on-entry`)

`docs/stage_f.md` §3, "Rays at the level faces". Two new defaults, refined
runs only: `amrex_split_on_entry = 1` (a ray entering a finer level splits
at once, without the kernel's radius clamp, which pulled its children back
toward the source); `amrex_ghost_rule = conserve` (a fine ray's deposit in
a ghost cell outside its level goes to the parent with weight
1/(1 + √fgam), photon conserving; the previous rule created photons just
outside the level: up to +0.008 in x_HII there). With both, the cells
outside a refined level are as accurate as a uniform base-grid run. The
`parent` and `half` rules stay as options (not recommended: `half` loses
photons in thick cells). Two refined levels were not rank-independent at the 1e-10 level
(fixed 2026-10-22, above).

## Done 2026-10-17: the reduced speed of light by ionisation state (branch `c-schedule`)

`docs/production_run.md` §1, "By ionisation state". Decision: `docs/PLAN.md`,
"Reduced speed of light by ionisation state".
* **`c_reduced_schedule = c0 x1 c1 x2 c2 ...`** (apps/radhydro, default
  off = the original): each hydro step's ĉ from the box's x_HII,V at the
  end of the previous step (x = 0 until the RT has written its grid).
  Stateless, so restarts are exact with no checkpoint change; on the AMReX
  backends rank 0 computes x and broadcasts it.
* **Production default** (`make_production.py --c-schedule`, user
  decision, revised 2026-10-18): 0.1c for x < 0.25, 0.2c to 0.5, 0.3c to
  0.75, 0.5c after (was 0.1c / 0.2c from 0.3 / 0.4c from 0.8). An explicit
  `--c-reduced` alone keeps a constant speed; `calibrate_fesc.py` and the
  stored-tolerance AMReX test are unchanged.
* **Tests:** `c_schedule` (a one-speed schedule is the constant
  `c_reduced` byte for byte; the switching rule), `restart_smoke_c_schedule`,
  and on the backends `amrex_state_restart_smoke_c_schedule` and
  `amrex_state_c_schedule_r2` (2 ranks). Fast suite 141/141, `amrex_state`
  9/9 on klone, nightly `amrex_state_full64_vs_serial` passes (114 s).
* **The nightly 64³ backend test is pinned to 0.01c** (`--c-schedule none`,
  the speed its tolerances were measured at). With the schedule it failed on
  per-cell fHI (0.57 rms against 0.15): both runs switched ĉ at the same
  steps, but at 0.1–0.4c the snapshot Γ is ray-sampling noise, independent
  between the serial and backend runs, and the run took 70 min.
* **Next:** (f6, β) of `docs/fesc_evolution.md` were calibrated at a
  constant 0.2c and need a recalibration with the schedule.


## Done 2026-10-17: a falling escape fraction (branch `fesc-evolution`)

Report: `docs/fesc_evolution.md`. Decisions: `docs/PLAN.md`, "f_esc
evolution thread". Parameters only; no code change.
* **Model:** f_esc(z) = f6 ((1+z)/7)^β (`f_esc`, `f_esc_beta`,
  `f_esc_z0 = 6`, the pivot redshift), already in the sources.
* **256³ calibration** (20 Mpc/h, the 512³ cell size; `otvcl256` with only
  the source lines changed; `tools/hyak/fesc_runs.sh`): β = 0, 2, 2.5, 3, 5,
  7 with f6 = 0.23, 0.15, 0.135, 0.12, 0.075, 0.045 (f6 for z50 ≈ 7 from a
  one-zone model; z50 came out 7.04–8.04). β ≥ 5 does not finish by
  z = 5. β = 2 and 3 bracket Γ_HI and ⟨F⟩; λ_mfp ∝ Γ^≈1 at fixed z.
* **512³ run `fesc512_b2.5`** (β = 2.5, f6 = 0.135; otherwise
  `reion512_rays` + `recomb_rad_cell_onspot = 1`): z99 = 5.51 (predicted
  ≈ 5.5), z50 = 7.11, τ_es = 0.056. At z = 5.1–5.4 Γ_HI = 0.97–1.00×
  Gaikwad+23, ⟨F⟩ 1.1× Bosman+22 (A = 1.1); λ_mfp 1.6–1.9× the observed
  at z = 5.1–5.7 (from 4.2–5.1×), 2.9× at 5.9. P_Γ: the rays' pattern
  (r ≥ 0.99), relative power 1.5–3× during the end of reionization,
  within 10–20 % by z = 5. Movie `docs/movies/fesc_evolution_512.mp4`.
* **Left for discussion with the user:** the remaining ≈ 1.8× λ excess at
  fixed Γ (resolution: 1024³ or AMR; a flag-gated subgrid Lyman-limit
  opacity; or accept it), §6; and the recalibration of (f6, β) with the new
  ĉ schedule (§7).
* Tools: `fesc_summary.py`, `fesc_post.slurm`. Census of all reionization
  runs: `$RUNS/reion_census.md`. The finished runs' checkpoints (459 GB)
  were deleted with the user's OK (quota at 99 %).

## Done 2026-10-16: refined-run crash fix, abort on error, `fields_levels` (branch `fields-levels`)

- **Crash of the refined runs** (`fesc512_b2.5_amr` at step 23, `amr_val128/amr128` at step 20: "a rank with rays to add but no box"): a refined level with fewer boxes than ranks; rays parked or leaving a finer level arrived on a rank without a box of that level. `add_rays` now hands such rays to the owner of box 0 (commit a90c2ec on `claude/pensive-goodall-kvpwgz`). `fesc512_amr` resumed from its checkpoint with the fix and PR #42's RT fixes and passed step 23.
- **An error on one rank ends the job** (`backend_abort`: `MPI_Abort` with more than one rank, from the error handlers of radhydro, rt_driver, pm_driver, hydro_driver; 72f7731). Before, the failing rank hung in the collective `amrex::Finalize` and the others waited: `fesc512_amr` held 5 nodes for 9.5 h after its crash. Tested: an error on rank 0 alone ends a 2-rank run in about 30 s.
- **`fields_levels = 1`**: a plane at the finest resolution of the RT levels with every field output (`levels_slice.bin/.txt`; `docs/full_run.md`, `docs/stage_f.md` §3); CTest `amrex_fields_levels`; `tools/analysis/amr_figures.py` (the AMR figures of `docs/stage_f.md`).
- **Open: two refined levels are not rank-independent.** The 8-source problem with two levels differs between 2, 3 and 8 ranks from step 8 on (relative 1e-10 in x_HII, the photons in flight); also with the binary before the `add_rays` fix, so it predates it. One level was checked byte-identical across ranks in F3; two levels never were. To find (the RT hierarchy with nested levels, or the hydro's fine layouts, which AMReX chops by the rank count).
- Jobs after the 2026-10-16 cluster restart: `fesc512_b2.5_uniform` finished (z = 4.98, 13.4 h); `fesc512_amr`, `amr_val128/{amr128,u256}` resubmitted on `build-amrex-fix2` (resuming from their checkpoints); the cost probe has no usable checkpoint (`reion512_40mpc_f0.23`'s was removed; `fesc512_b2.5`'s has recomb_rad = 1).

## Done 2026-10-14 (evening): F3 validation with 8 sources (branch `stage-f3-8src`)

**`docs/stage_f.md` §3, "Eight sources".** `make_ifront.py --sources 8src`
puts the 8 sources of `tests/problems/uniform_8src` in the full-loop
problem; `tools/hyak/rt_amr_ifront8.slurm` runs uniform 32³, 32³ with one
level over a box where fronts cross its faces and overlap there, uniform
64³, then two nested levels against uniform 128³;
`tools/analysis/ifront_multi.py` compares x_HII,V/M, the photons in
flight and the ionised volume around each source, and checks for jumps.
Both pass: with one level, x_HII,V/M and the photons in flight are better
than the 32³ run's at steps 3–7 and equal within 0.0006 after; per-source
radii within 0.013 cells of the base run's error (sources 3 and 6, whose
fronts meet at the level's faces, are 0.03 cells worse late); no jump at
any crossing (≤ 0.026 cells). Two levels: errors as the base run's (per
source 0.13 against 0.10 cells). With the whole box refined the hierarchy
reproduces uniform 64³ with `rt_cells_per_substep = 2` to 1e-4. Nightly
CTest `amrex_rt_amr_ifront8`. On the way:
- **Fix: the covered cells' fgam.** A parent ray crossing a cell under a
  finer level lost 1 − fgam of that cell, which was one step stale (no
  physics there; only `average_down` refreshed it); rays entering the
  level through newly ionised gas lost photons (4.2 % of the photons in
  flight at step 10, 14 % near the level's faces). Now the mean of the
  children's fgam of the same substep (`fill_covered_fgam`).
- **Fix: each level's ray pool** is `Nray_max × cells(l)/N³` (the split /
  merge thresholds and the pool guard used the base pool on every level).
- **Fix: `ngam_rays_total` in base units** (it summed the levels' rays in
  their own units).
- **Correction:** the 1.5 % lag of the whole-box one-source run reported
  on 2026-10-15 was not the coarse substep's (k = 2 is 0.013 cells *ahead*
  of k = 1 at 64³) but gravity: `make_ifront.py` wrote Zel'dovich dark
  matter, and the refined runs (PM on level 0 only) grew different gas
  structure. The problem now has the dark matter at rest; the one-source
  refined run is then 3–4× closer to 64³ than the 32³ run at every step.
- `amrex_refine_box` takes one box per level (nested regions);
  `RADRAY_AMR_DEBUG=1` prints each level's bounding box.
- Tests: `ctest -L amrex -LE nightly` 176/176 on klone; the two nightly
  I-front tests pass. `amrex_max_level = 0` is bitwise as before.
- Open: a residual of a few hundredths of a cell at the level's faces
  (more ionisation just outside), likely the fine-ghost deposit sent up
  as photons over the whole coarse cell (absorbed 1 − fgam there, the ray
  losing 1 − √fgam); unverified.

## Done 2026-10-15: Stage F3, the RT across refinement levels (branch `stage-f3`)

**`docs/stage_f.md` §3.** F3a: the gas state holds every refined level
(the hydro copies them in and out per call, the checkpoint writes the
state's floats). F3b: the RT hierarchy in `src/amrex/rt_level.cpp`, one
`Level` per level of the state with its own cells, deposits, rays, units
and RT parameters; the substep in four phases interleaved across the
levels; a refined level's outside ghosts take the parent's tau (halved)
and transmitted fraction (√fgam); the deposits cross the level boundary
both ways and the rays migrate after every substep, with the photons
rescaled by the cell volume; a level's cells under the finer level emit
nothing and have no physics; the finer levels are averaged down after
the step; the statistics are summed over the levels; the sources go to
the finest cell holding them; the checkpoint carries every level's rays
in base units. `amrex_max_level = 0` is bitwise Stage E (176/176 on
klone); the refined smoke runs are byte-identical across ranks and
through a restart. Validation: one source in uniform gas (`make_ifront.py`,
`rt_amr_ifront.slurm`, `ifront_full.py`): the refined 32³ run tracks the
uniform 64³ front to 0.004–0.14 base cells inside the refined region and
across its boundary (the lag is the two-cell coarse substep's, the user's
choice); CTest `amrex_rt_amr_ifront` (nightly). Traps met: the per-level
periodic wrap, the hydro's push-down after the RT, the photon
normalisation across levels, the transmitted fraction of a half cell.
Left for later: the PM on the refined levels (Stage G), FoF/P(k) on
rank 0, the pancake/KH suite cases with refinement, the flux-register
and FillPatch overheads, a refined production run (and the sink's
production threshold).

## Done 2026-10-14: Stage F decisions; B17 under `legacy_physics = 0`; Stage F1 and F2 (branch `stage-f-decisions`, PR #40)

- **Stage F decisions taken** (`docs/PLAN.md`, "Stage F decisions"): one
  coarse RT substep across levels (rays cross 2^l fine cells, no time
  sub-cycling of the transport); B17 (`photon_conserve`) part of
  `legacy_physics = 0`; cooling on in dense gas with an aggressive
  star-formation sink (overdensity of a few thousand); the rest the
  implementer's choice; order: gas AMR with refinement off (bitwise),
  refined hydro on the test suite with the cooling flag and the sink,
  then the RT across levels.
- **B17 folded into `legacy_physics = 0`:** the key `photon_conserve`
  defaults to −1 = `!legacy_physics`, 0 | 1 override
  (`Config::photon_conserve_effective`, both drivers, `make_production.py
  --photon-conserve -1|0|1`). The Fortran harnesses mirror it (RT harness
  PATCHES.md 14: `photo_cell` iterates while any group is photon-limited,
  `nit_max = 1000`; the full harness reads the key too), so every
  combination of the mirrored flags stays byte-identical. The stored
  `legacy_physics = 0` references (`outputs_legacy_rt0_phys0`,
  `outputs_c0.1_fixed`, harness_full `*_rt0`) are regenerated on the CI
  platform by the new on-demand workflow
  `.github/workflows/regen_references.yml` (artifact `references`), with
  `reference/harness/store_outputs.py` for the RT subset; the legacy sets
  are untouched.
- **Stage F2, refined hydro validated** (`docs/stage_f.md` §2): the hydro
  suite with refinement (`run_suite.py --amr L`, one klone node): one
  level on 32³ matches or beats the uniform 64³ run (Sod L1 ρ 0.0098
  against 0.0120; Sedov's shock radius 4.3 % off against 3.0 %, two levels
  3.5 %; the strong shock 0.16 → 0.086 → 0.036). Three traps on the way,
  all fixed with tests (§1): the fine levels' dtdx (half speed per level
  before), the interpolation interval reset per call, and the
  positivity of the interpolater (`mf_linear_slope_minmax_interp`).
- **The star-formation sink (F2b)** (`docs/stage_f.md` §2):
  `sf_overdensity`, `sf_efficiency`; star particles on the PM level (the
  dark matter deposits bitwise as before), `DM/stars.<ztag>`, `sink.txt`,
  checkpoint version 8, rank-independent ids; tests label `sink`.
- **klone: `srun --mpi=pmix`** (`docs/cluster.md` §3): the backend's
  intermittent start-up crash is a use-after-free inside `ompi/4.1.6-2`'s
  PMI-1 `MPI_Init` (found with an address-sanitizer build); `pmix` avoids
  it; the test launcher and the job scripts pass it.
- **Stage F1, the gas hierarchy** (`docs/stage_f.md` §1):
  `amrex_max_level` (0 = Stage E, bitwise: `ctest -L amrex_hydro` 27/27
  on klone), AmrMesh regridding on the gas overdensity
  (`amrex_refine_mass` · 8^l), FillPatch ghosts, two substeps of dt/2 per
  parent substep, refluxing and average-down, the CFL over every level,
  the push-down of the level-0 changes the RT and PM make between hydro
  calls, the refined levels through the checkpoint (`<checkpoint>.amr/`)
  and across a rebalancing. Tests (label `amr`): conservation against
  the uniform run and layout independence at 1 and 2 levels.
- `cool_suppress_overdensity` (§2; 200 = the original, 0 = dense gas
  cools), `make_production.py --cool-suppress-overdensity`.
- Next: Stage F3, the RT across levels (`docs/stage_f.md` §3); the
  pancake, KH and full-tier suite cases with refinement; the flux-register
  and FillPatch overheads (to profile).
- Storage: the Stage D/E one-step timing outputs (581 GB) were deleted on
  the user's instruction.
## Done 2026-10-13: Stage E2 on klone, the PM backend's validation and scaling (branch `stage-e-pm-klone`, PR #38)

Report: `docs/stage_e_pm.md` §3.1 (64³) and §4.1 (1024³). Decisions and
choices: `docs/PLAN.md`, "Stage E2 on klone". Nothing changes a result;
the serial PM is untouched.
* **Tests on klone** (`build_amrex.sh`, AMReX 26.10 + FFT): the labels
  `amrex_pm`, `amrex`, `amrex_state`, `pm` 185/185 (at the start and again on
  the final branch); the fast suite of the default build 137/137.
* **64³ production configuration** (`make_production.py --N 64`,
  `legacy_rt = legacy_physics = 0`, unsplit): the three backends with
  `amrex_state = 1` against the serial code. x_HII,V 1.2e-6, x_HII,M
  1.5e-6, T_V 3.2e-5, T_M 1.9e-5, τ_es 2.2e-7; per cell at z = 5, fHI rms
  5.4e-2, T 7.1e-4, ρ 5.2e-5. The fHI difference is the rays' Γ noise
  reshuffled by the PM round-off (the RT-only and RT + hydro runs give
  1.6e-4). New nightly test `amrex_state_full64_vs_serial` (95 s).
  Also fixed the nightly `amrex_full64_vs_serial`, which crashed on the
  frame layout of `be49827`.
* **`pm_refine = 2`, 1024³ particles** (`hmf512_40mpc_r2` restarted for 3
  steps, 24 ranks × 8 threads per node):
  - Bitwise identical across 1/2/4 nodes and box sizes 32/64 for each FFT
    (every step file, the particles, P(k), the halo catalogues). The
    amrex FFT is bitwise as well; `pencil` stays the default.
  - Against the serial code: round-off (gas ≤ 5e-6 rms, particles
    3e-6 cells rms, P(k) 1e-6, a few halos out of 4.6·10⁵).
  - PM per step: serial 11.2 s, backend 4.1 s (1 node), 3.0 s (2), 1.3 s
    (4); hydro 15.9 / 3.1 / 2.1 / 0.9 s.
  - But the step is 3× slower than the serial code (224–232 s against
    74 s): FoF and P(k) run on rank 0 with 8 threads (188 s and 48 s
    against 32 s and 20 s). Peak RSS: rank 0 120–126 GB, other ranks
    7–12 GB, serial 168 GB.
* **Fixes:**
  - `init_pm(..., host_arrays)`: with `amrex_state = 1` the non-I/O ranks
    no longer allocate the host PM arrays (45 GB per rank, OOM on one
    node).
  - The amrex FFT's planner reset is guarded (klone's AMReX never
    initialises FFTW threads; the reset ran `fftwf_cleanup` and the plans'
    destruction segfaulted at exit).
* **Tools:** `--set step_timers=1` (output only),
  `tools/hyak/amrex_pm_scaling.slurm`, `amrex_pm_scaling_table.py`,
  `amrex_pm_vs_serial.py`.
* **Open:** the distributed halo finder and P(k) (E3's open item), now the
  bottleneck of the PM-heavy runs; until then, FoF less often
  (`fof_dt_myr`).

## In progress 2026-10-14: the moment RT on AMReX (branch `moments-amrex`)

Design: `docs/stage_e_moments.md` (approved 2026-10-14, PLAN
"Moments-on-AMReX thread"). `rt_backend = amrex` with `rt_method = m1 |
otvet` on z-slabs (`amrex_slab = 1`), every sum in the serial order.

* **Phase 0, done** (`f6f971c`): `amrex_slab = 1` (z-slabs in z order, also
  the shared decomposition of `amrex_state = 1`); sparse sources from
  `calc_halo` (`RTGrid::src_*`, byte-identical; no N³ `Luv`/`sfr`/`per_cell`
  on the ranks other than 0); `src/amrex/exact_sum.hpp`.
* **Phase 1, M1, done** (`f6f971c`): `src/amrex/moment_level.cpp`, the
  serial per-cell arithmetic shared through
  `include/radray/kernels/moment_cell.hpp`; `amrex_moment_stats = exact |
  fixed`.
* **Phase 2, OTVET thin / split / screen, done** (`63af7ab`): the direct sum
  per cell, or the serial FFT one component per task (`OtvetConvolver`).
* **Phase 3, the cone, done** (`68233e5`): 12 FFT tasks, mip-maps per slab
  with the finest level used gathered, the march per slab, the budget's
  receiver sum as a chain, the Gaussian's z pass on y-slabs.
* **Phase 3b, `recomb_rad = 1`, done** (`d25c694`). Needed for the 1024³
  comparison: the 512³ cone run used `recomb_case = A`, `recomb_rad = 1`.
* **Phase 4, checkpoints, done** (`098a9f9`): the format unchanged, streamed
  through rank 0; a serial moment checkpoint restarts on the level
  byte-identically (so the serial 256³/512³ checkpoints seed phase 5).
* **Validation:** label `amrex_moments` (84 tests): `rt_driver` against the
  serial moment RT for M1, thin (direct, FFT), split, screen (fixed and
  budget λ), cone (λ_x from the budget, no blend), recomb_rad, on 1-4
  ranks × 1-2 threads × slab planes; the full loop with the state on the
  slabs against the serial moment RT; restarts. All bitwise.
  The fast suite: only the 10 `platform_ref` tests fail, as on `main`
  (outputs byte-identical to `main`'s).
* **Threaded OTVET FFT plans, done** (`fd1d8f0`): `amrex_moment_fft_threads`,
  the guard at the first convolution (tests `otvet_fft_threads_{64,128}`).
* **Phase 5, 256³ on klone (2026-10-14):** the serial production runs
  `otvcl256_*` restarted from their step-40 checkpoints on the level with
  the serial PM (`tools/hyak/moments_validate.slurm`, MODE = host; the
  AMReX PM agrees with the serial PM to round-off only, so the production
  configuration is compared with the serial RT by the full-loop pair
  tests instead): M1, thin, split, screen + split bitwise the production
  runs' step 41 (19 files each); the cone bitwise the current serial code
  on 1 and 2 nodes (7 files, 12 movie frames). The cone production run
  itself predates the final cone code (its step 41 differs from today's
  serial code, λ_x 760.865 vs 760.866).
* **Transport:** rewritten on raw indices (`5a4f5d4`), GLF 3× faster;
  the FFT tasks of all kernels computed at once (`9637d0e`; the cone tensor
  107 → 38 s at 512³ on one node); a movie frame gathers only its 6 fields
  (`44e431d`).
* **Phase 5, 512³ cone on klone:** bitwise the production run (steps 53-54,
  8 files, 12 movie frames) on 1 and 2 nodes (host mode; the current serial
  code reproduces those steps too). Production configuration on 1 node: RT
  293 s against the serial 495 s for the two steps, run 415 against 594 s
  (`docs/stage_e_moments.md` section 10). Pending (ckpt-g2 queue, frequent
  preemptions): 4 nodes bitwise; 2 and 4 nodes timed.
* **Fast suite on the final code (`44e431d`):** 411/422; the 10
  `platform_ref` tests as on `main`, and `end_to_end_legacy0`, which hung
  in MPI_Finalize (a Slurm PMI barrier of a singleton run inside an srun
  step) and passes alone (13.5 s).
* **Phase 6:** `runs/reion1024_otvet_cone` (ICs with the 512³ run's large
  scales), job 41889599 queued on ckpt-g2 (4 nodes × 16 × 12, binary
  `scratch/moments_bin/radhydro_44e431d`). Disk: ~1.4 TB expected
  (0.85 TB of checkpoints), 2.0 TB free on 2026-10-14.
* **Then:** the comparison with the 512³ cone (`gamma_power.py`,
  `reion_power.py`, `rt_methods_runs.py`, `igm_stats.py`).

## Done 2026-10-13: a better Eddington tensor for OTVET (branch `otvet-closure`, PR #35)

Report: `docs/otvet_closure.md`. Decisions and choices: `docs/PLAN.md`,
"OTVET-closure thread". `main` merged in after PR #31.

* **Step 1, ground truth: done.**
  - `tensor_dump = 1` (radhydro): the rays' Γ-weighted S0, S1, S2 per
    cell at the zdata1 snapshots, with odd substeps for the noise, the
    thin OTVET tensor and the sources. The moment methods dump the same
    moments of their own N, F, D N.
  - Off by default and byte-identical. Tests `tensor_dump_rays`,
    `tensor_dump_otvet` and `thread_full_tensor_dump`; fast suite 127/127
    after the merge.
  - 256³, 20 Mpc/h runs: `otvcl256` (rays), `_m1`, `_otvet`.
  - Analysis tools: `tensor_compare.py`, `cone_tensor.py` (numpy
    prototypes) and `tensor_summary.py`.
* **Findings.**
  - Before overlap, OTVET is good in the bulk but too isotropic in shadowed
    and low-Γ cells.
  - During overlap it is too isotropic on large scales: it adds sources
    behind neutral gas.
  - After overlap it is too anisotropic by λ_eff/0.61 L, Wu+21's
    normalisation, which is measured here. Its large-scale force is
    anti-correlated with the truth.
  - The 256³ runs reproduce the 512³ symptom: OTVET's P_Γ excess grows from
    1.1–1.4× at z = 6.4 to 15–19× at z = 4.9, with r ≥ 0.9.
* **Step 2, candidates: done.** Ranking: cone + budget + blend, then the
  split, then the screened kernel. The user approved implementing all three
  and running them at 512³.
* **Step 3, implementation: done** (`include/radray/otvet_closure.hpp`,
  `src/otvet_closure.cpp`, `docs/otvet_closure.md` §5).
  - `otvet_closure = thin | screen | cone`, `otvet_split`,
    `otvet_screen_lambda`, `otvet_cone_*`. The defaults run the original
    code path, byte-identical (smoke and lognormal runs against e55def7).
  - `apps/otvet_tensor` checks the C++ against the numpy prototypes on the
    rays' snapshots: ≤ 2.5·10⁻⁵ S0-weighted at 256³.
  - New unit and ctest tests (threads, restart, tensor dump); fast suite
    green.
* **256³ full runs** (`otvcl256_otvet_{split,screen,screensplit,cone}`,
  45–47 min each; `docs/otvet_closure.md` §6). P_Γ / rays on the hydro-step
  means after overlap, at k ≤ 7 h/Mpc:

  | run | P_Γ / rays |
  |---|---|
  | split | 0.85–1.09 |
  | cone | 0.69–1.29 |
  | screen | up to 1.5 |
  | OTVET | 3–17 |
  | M1 | 0.1–0.5 |

  - r(k) ≥ 0.92 for all closures.
  - The cone also has the best overlap-era morphology and the closest end
    of reionization (x_V = 0.99 at z = 6.05, against the rays' 6.00).
  - The split and the cone match the upper Γ quantiles. Their lower tail is
    too narrow by z = 4.9.
  - `gamma_power.py --s0` uses the tensor dump's step-averaged Γ. The rays'
    snapshot Γ is dominated by single-substep sampling noise on cell
    scales.
* **The cone's photon budget counts the ionised cells only** (`c53a503`).
  Before overlap the unattenuated near field had driven λ_x to 1 cell, and
  the first 256³ cone run's own tensor at z = 9.6 was worse than OTVET's.
  The rerun (`otvcl256_otvet_cone2`) fixes that; nothing changes after
  overlap.
* **512³ runs: done** (`reion512_rays_otvet_{split,cone,screen}`;
  `docs/otvet_closure.md` §7; analysis in `$RUNS/otvet_closure_512`).
  - **Wall time:** split 4.90 h, screen 4.77 h, cone 5.22 h; thin OTVET
    4.80 h.
  - **After overlap:** OTVET's 2.3–9× P_Γ excess is gone with every
    closure. The split keeps a 1.8× excess on the box scale; the cone does
    not.
  - **The cone is the best overall.** End of reionization at z = 5.838
    against the rays' 5.816 (OTVET and split 5.868). It has the best
    overlap-era morphology, x_HII power and temperatures, and the rays' P_Γ
    during overlap (r = 0.98–1.00).
  - **Still pending:** the IGM statistics and the comparison movie
    (job 41825950) run after the Oct 13–14 maintenance.
* **2026-10-13:** PR #35 merged into `main`.
  - **Movie:** `$RUNS/otvet_closure_512/compare.mp4` (rays, M1, OTVET,
    split, cone, screen).
  - **After-overlap P_Γ figure** (`docs/otvet_closure.md` §7) at the epochs
    of `rt_methods_512.md` §6. The case-B rays differ from the rays by
    4–16 %; the split is within 0.77–1.06 except 1.8× on the box scale; the
    cone within 0.72–1.16 except 0.5× at k ≈ 0.5 h/Mpc.
  - **Faster tensor** (branch `otvet-fft-speed`): cached kernel spectra
    and component-parallel FFTs, byte-identical. At 256³: thin 5.1×,
    screen 3.1×, cone 2.8×. Fast suite 137/137.
* **Open question for the user:** make the cone the recommended OTVET
  setting?
* **Open:** the user referred to a "previous message" with numbered
  instructions that never reached the session. Only its point 1, the
  merge of `main`, is done.

## Done 2026-10-13: Stage E4, load balancing and production runs on the backend (branch `stage-e-balance`, PR #37)

**`docs/stage_e_state.md` §5.** `amrex_balance = 1` (with `amrex_state =
1`): before every RT step the shared decomposition is rebuilt on the cost
cells + `amrex_balance_ray_weight` × live rays per box (space-filling
curve or knapsack) and, when the balance efficiency gains more than
`amrex_balance_min_gain`, the gas, the RT cells and rays and the
particles move onto it. The levels are mapping-independent, so the run
is byte-identical to the fixed mapping: `amrex_balance_full_pair_r2/r3`
(remapped at every RT step, knapsack and sfc) in `ctest -L amrex_state`,
26/26 with the restart label on klone. The checkpoint's rays are now
gathered box by box in a canonical order (sorted by home cell and bytes),
so the file does not depend on the mapping or the rank count. Also
fixed: `radray_amrex` links the threaded FFTW itself (its objects reach
the PM level now), and `env.sh` prefers the threaded float FFTW over a
stale `RADRAY_FFTW` in the login environment.
- **klone:** the RT step at 512³ is 4 % faster on 2 nodes, 13 % on
  4 and 10 % on 5 nodes, the whole allocation (936, 470 and 398 s
  against 974, 538 and 444 s), at 256³ 11 % and 12 % (74 and 42.5 s),
  all bitwise identical to the fixed mapping; 512³ on 5 nodes: the step
  506 s, 5.7× the serial code on one node (§4). One fixed-mapping run
  at 256³ on 4 nodes hung inside the RT step with identical data up to
  the hang; eight further runs of the case were clean and bitwise
  identical (§5). For the next one the backend writes every rank's call
  stack to `Backtrace.<rank>` on SIGTERM (`scancel`, the time limit),
  `tools/hyak/backtrace_lines.py` resolves them, the job scripts run in
  the run directory and the klone build keeps line tables.
- **Stage E is complete** (hydro, PM, the state on the levels, load
  balancing, scaling on the whole allocation, production runs). Left
  over, none blocking: rank 0's restart transient for serial
  checkpoints (the serial ray pool, ~100 GB at 512³; backend
  checkpoints stream), the per-step mirror gathers and the host P(k)
  and FoF on rank 0 (deferred by decision until round-off there is
  accepted; with FoF every step at 1024³ particles they are 90 % of the
  step and make it 3× slower than the serial code, `docs/stage_e_pm.md`
  §4.1), `Luv`/`sfr` on every rank. Next: Stage F (AMR); the
  user's decisions (rays across levels, `photon_conserve` under
  `legacy_physics = 0`, cooling in dense gas) were taken on 2026-10-14
  (`docs/PLAN.md`, "Stage F decisions"); and the group's
  storage: `/gscratch/astro` is 96 % full, the timing runs of this and
  the Stage D/E threads hold about 450 GB under
  `runs/reion512_40mpc_f0.23/amrex_*` and 40 GB under
  `runs/reion256_f0.2/amrex_*`, all one-step outputs that can go.
- **Production runs on the backend** (`docs/stage_e_state.md` §6):
  `make_production.py --backend amrex` writes the backend configuration
  into `params.txt`, and `reion_production.slurm` launches MPI ranks
  when it sees it there (24 × 8 per node, the module and UCX settings,
  restart/`DONE`/post-processing unchanged). Validated at 128³ to z = 5
  against the serial code (histories to 1 in the fifth digit, 20 %
  less wall time on one node) and through a time-limit kill and
  resubmission (bitwise the uninterrupted run).

## Done 2026-10-13: Stage E3, the state on the levels (branch `stage-e-state`, PR #34, merged)

**`docs/stage_e_state.md`.** `amrex_state = 1`: the gas lives in one
MultiFab on a decomposition shared by the RT, hydro and PM levels, which
read and write it locally; the particles live in the PM container
(drifts and `Redistribute` on the level); the host grids and particle
arrays exist only on rank 0, as a mirror refreshed by slab gathers for
the statistics, the outputs and the checkpoint (formats unchanged); the
halo finder runs on rank 0's mirror and the catalogue is broadcast.
Needs the three backends, `hydro_scheme = unsplit`, `rt_method = rays`.
- **Byte-identical** to the same three backends with the host state on
  every rank (`amrex_state = 0`): the smoke problem at 1/2/3 ranks and 2
  threads, every file (fields, stats, frames, checkpoints), and a
  bit-identical restart (`ctest -L amrex_state`, 5 tests); the whole
  `amrex` label 159/159 on klone with the phase-2 build.
- **klone, 256³ on one node:** bitwise identical to the host-state
  three-backend run; the ranks other than 0 hold 4.5 GB (8 GB with the
  host state, 12.8 GB in Stage D); RT 151 s. **512³ on one node runs**
  (24 ranks × 8 threads, rank 0 293 GB in the restart, the others 35 GB
  with the phase-1 build): RT 1756 s, step 1849 s, 1.48× faster than the
  serial code on the same node (Stage D needed 2 nodes, or 8 fat ranks
  at 3258 s), and **bitwise identical at 512³** to the 2-node run with
  the same backends and the host state; with phase 2 the other ranks
  hold 24 GB. On 2 and 4 nodes: RT 974 s and 538 s (1.83× per doubling,
  bitwise identical to the one-node run), the other ranks at 13 GB; the
  step on 4 nodes 598 s, 4.8× the serial code on one node.
- Costs to remove next: the per-step mirror gathers (60 B per cell, for
  bitwise statistics; a 128-bit on-level path where round-off is fine),
  rank 0's restart transient (the serial ray pool and list of the
  checkpoint), the host P(k) and FoF on rank 0, and `Luv`/`sfr` on every
  rank.
- Also: the checkpoint's backend rays are streamed in chunks (written
  while gathered, scattered while read: rank 0 holds one chunk, and the
  2 GB `MPI_Gatherv` limit is gone); the live-Fortran full-loop tests
  link again in an AMReX build (`radray_backend_stubs`, which they
  could not since Stage D); `build_amrex.sh` builds AMReX with the FFT module on
  klone (`AMReX_FFTW_IGNORE_OMP`); `amrex_hydro_max_grid_size = 0` means
  `amrex_max_grid_size`.

## Done 2026-10-12: Stage E2, the PM gravity on AMReX (branch `stage-e-pm`)

**`docs/stage_e_pm.md`.** Item 2 of Stage E (`docs/framework_decision.md`
§6): `pm_backend = amrex` (default `serial`, byte-identical). The dark
matter in an AMReX `ParticleContainerPureSoA`, the CIC deposit in int64
fixed point with `SumBoundary`, the Poisson solve, the mesh force, the DM
kick and the gas kick (`calc_kick`, `calc_kick_unsplit`, `pm_refine`) on
one AMReX level of the PM mesh (`src/amrex/pm_level.cpp`; stub without
AMReX). Stage E2's MPI model: the hydro and the RT hold the full state on
every rank, so the level copies the host state in at the start of the PM
step (each rank creates the particles of its own boxes from `xv_dm`) and
gathers the kicked velocities and gas back onto every rank (bitwise-OR
all-reduce of disjoint pieces); the two drifts stay on the host; the level
holds no state between steps (checkpoints unchanged).
- **Shared kernels:** `include/radray/kernels/pm_kernels.hpp` (CIC, force
  stencil, Green's function), called by `src/pm.cpp`: the `pm` label
  unchanged (stored gfortran references byte-identical).
- **FFT:** `amrex_pm_fft = pencil` (default; per-line 1-D FFTW transforms
  on slab layouts, rank- and thread-independent by construction) or
  `amrex` (`amrex::FFT::R2C`; also bitwise rank/thread-independent in
  every test here, not guaranteed by construction). AMReX needs
  `-DAMReX_FFT=ON` (FFTW double and float on the CPU; `build_amrex.sh`
  builds the double library on klone).
- **Tests, `ctest -L amrex_pm` (65, all pass here; the `amrex` and `pm` labels and the fast suite unchanged):** Ewald single-particle
  force, no self-force, gas kick, momentum (MF/non-MF/unsplit, r = 1, 2, 4)
  and backend-vs-serial fields on the level (1/2/4 ranks, both FFTs);
  round-off against `pm_driver` on the 8 PM problems (1e-4 relative or
  1e-5 of the maximum; growth_3d reaches 5e-6 of the maximum force after
  30 steps); the physics checks (momentum, pancakes incl. r = 2, growth,
  the r = 2 twin bit-identical) on the backend's runs; bitwise identical
  for (ranks, threads, box) ∈ {1,2,4,8} × {1,2} × {32,16,8} on pm_refine
  1 and 2, both FFTs; the full loop with `rt_backend = amrex` +
  `pm_backend = amrex` rank-independent (1, 2, 4), thread-independent,
  bit-identical restarts (legacy, unsplit, pm_refine = 2).
- **Full loop vs the serial PM:** the 16³ smoke problem amplifies the
  force round-off (1e-6 of the maximum) to 2e-5 in T_V by step 15 (the RT
  backend's round-off: 1e-8–3e-7); the 64³ production-configuration
  nightly comparison is `amrex_state_full64_vs_serial` (2026-10-13, above).
- **Timing (4-core VM, `pm_driver`, 3 steps, Ngas = Ndm):** per step, serial PM on 4 threads
  0.39 s (128³) / 2.74 s (256³); backend 1 rank × 4 threads 0.44 / 2.74 s,
  1 × 1 0.85 / 7.07 s (serial 1 thread 0.83 / 7.16 s), 4 ranks × 1 thread
  1.02 / 8.1 s. On one rank the backend matches the serial OpenMP PM; more
  ranks cost more because every rank walks the full host particle array
  and the gathers all-reduce 340 MB per step at 256³ (the Stage E2 model
  of a full host state, not the level).
- **Not done / open:** the particles staying on the level between steps
  (needs E1 and the end of the redundant host state); E3–E6; the GPU
  build (plain OpenMP loops here).

## Done 2026-10-12: Stage E prerequisites on klone (branch `claude/clever-fermi-aa4km4`)

**`docs/stage_d.md` §4.4 and `docs/stage_e_hydro.md` §4.1, §5.** The AMReX
backend built, tested and timed on `cpu-g2`, 1 to 4 nodes.
- **Build and tests** (build of commit `5c0c974`, Stage E1 included):
  `tools/hyak/build_amrex.sh` on a 32-core allocation, `ctest -L amrex`
  **81/81** (54 RT + the 27 `amrex_hydro` tests: the hydro backend
  byte-identical to the serial unsplit step for every rank, thread and
  box-size combination), launched with `mpirun --oversubscribe`; the
  default build's fast suite **125/125** outside `platform_ref`. Fixed in
  the repo earlier the same day: the module loading under `set -u`
  (`radray_module_load`), the test launcher flags, the slab-wise field
  broadcast and the chunked restart scatter (`d927abd`).
- **256³ RT** (`reion256_f0.2`, step 49, 144 M rays; 24 ranks × 8
  threads per node): RT 156.8 s on one node against the serial code's
  144.4 s on 192 threads (**9 % slower**: physics 1.8× faster, casts
  equal, plus sort, migration, boundary sums and the gather), 100.2 s on
  2 nodes (1.56×) and 76.8 s on 4 (2.04×; the spike's RT scaled 1.42×
  and 2.13×). Peak RSS 37.3 GB on rank 0 (the restart), 12.5–15 GB on
  the others (the redundant serial state, 12 GB per rank, does not
  shrink with more nodes), serial 30.6 GB. The fields are bitwise
  identical across 1, 2 and 4 nodes and agree with the serial run to
  round-off (rms 1e-4 fHI, 4e-4 Γ, 6e-7 T, ρ identical, single I-front
  cells up to 0.3–0.5; global statistics to 1e-9, ray counts to 3e-7).
- **512³ RT** (`reion512_40mpc_f0.23`, the only case B `legacy_rt = 0`
  checkpoint; `reion512_rays_caseB` kept none and `reion256_c01` has no
  `run/` at all, so the 256³ run is `reion256_f0.2`): step 53, 1.12 G
  rays, 12 ranks × 16 threads per node (24 × 8 ran out of memory in the
  restart). **One node does not fit** (12 ranks sum to 1.7 TB: 300 GB on
  rank 0, 126 GB on each other rank, the redundant serial state about
  100 GB per rank, 200 B per cell); 2 nodes RT 1460.7 s, 4 nodes 782.1 s
  (1.87×) against the serial 2593.4 s on 192 threads (1.78× and 3.3×);
  peak RSS 292.9 GB on rank 0, 100–112 GB on the others, serial 203.8 GB.
  The 2-node `SumBoundary` of 368 s (a quarter of the RT) was memory
  contention on the node holding rank 0, 50 s on 4 nodes. Fields to
  round-off as at 256³, bitwise identical between 2 and 4 nodes. These
  512³ runs are from the earlier build (`d927abd`); the RT level is
  unchanged since. Stage D needs ≥ 2 nodes at 512³ until Stage E removes
  the redundant state; 8 ranks × 24 threads on one node (2026-10-13) fits
  (rank 0 293 GB, the others 171–196 GB) but takes 3258 s for the RT,
  1.26× the serial code, fields bitwise identical.
- **The multi-node hydro slowdown is the UCX transport, not the ghost
  volume** (`docs/stage_e_hydro.md` §5; 256³ random3d, 2 steps, 24 × 8
  per node, every layout byte-identical to the serial run). With the
  klone module's defaults the hydro step is 8.1 s on 1 node, 7.9 on 2,
  12.3 on 4: no speed-up, all of it `FillBoundary` (0.44 → 3.9 → 6.8 s),
  independent of the box size (16³/32³/64³ within 15 %). UCX's
  rendezvous zero-copy path is the cost (per message, growing with the
  ranks per node; UCX over TCP as slow; registration-cache knobs,
  memory hooks and OpenMP waiting policy change nothing); **`UCX_RNDV_THRESH=inf
  UCX_ZCOPY_THRESH=inf`** brings `fill` to 0.3–0.5 s on 2 and 4 nodes
  (as on one) and the step to 5.3 s on 2 and 4 nodes against 7.0 on one,
  the rest being the Stage E1 gather (2.4–3 s) and the redundant host
  work (2 s). But eager for everything slows the RT's large transfers
  (the field broadcast, the ray gather, the deposit sums): 256³ RT on 2
  nodes 95.1 s against 100.2 with the module defaults, 512³ on 4 nodes
  844.5 against 782.1. **`UCX_RNDV_THRESH=8m UCX_ZCOPY_THRESH=inf`**
  (eager below 8 MB) is the compromise and the best of all for the RT:
  256³ 2 nodes 89.7 s (−10 %), 512³ 4 nodes 667.1 s (−15 %), hydro
  `fill` 0.34–0.38 s on 2 nodes in a fixed-node-pair sweep (the module
  defaults 4.0–4.4 s, 16 MB and above 0.5–1.5 s: the gather's 16–32 MB
  broadcast blocks then go eager too, `amrex_hydro_ucx_sweep.slurm`);
  every run bitwise identical to the module-default one. `tools/hyak/env.sh` exports it for every AMReX
  job (`radray_ucx_eager`; `RADRAY_UCX_RNDV_THRESH`, `RADRAY_UCX_EAGER=0`).
  Box size for production: **32³** (16³ costs 5–10 % more, 64³ leaves
  ranks idle at 4 nodes). The backend on one node is 4.7× faster than
  the serial OpenMP hydro (8.1 against 38.0 s).
- Substitutions: `reion256_c01` has no checkpoint, so the 256³ run is
  `reion256_f0.2`; the 512³ checkpoint without `recomb_rad` is the case
  B run `reion512_40mpc_f0.23`.
- `radhydro` reports the peak RSS of every rank (min/max in
  `timing.txt`); `amrex_scaling.slurm` runs on several nodes and compares
  the fields (`tools/hyak/amrex_scaling_compare.py`).

## Done 2026-10-12: Stage E1, the unsplit hydro on AMReX (branch `claude/clever-fermi-aa4km4`)

**`docs/stage_e_hydro.md`.** `hydro_backend = amrex`: the unsplit
MUSCL-Hancock + HLLC step with the dual-energy formulation on AMReX
MultiFabs, MPI + OpenMP (`src/amrex/hydro_level.cpp`), `hydro_scheme =
unsplit` only; the legacy split scheme stays serial by design.
- The face and update arithmetic moved into `include/radray/kernels/hydro_mh.hpp`
  (`face_all`, `face_neighbour_max`, `update_from_faces`), so the serial
  driver and the level share one implementation; the serial unsplit
  outputs are bitwise unchanged.
- **Bitwise identical to `hydro_backend = serial`** for 1–8 ranks × 1–2
  threads × box sizes 32/16/8 on `sod_x`, `random3d`, `sound_x`, and in
  the full loop with either RT backend at 1–3 ranks; bit-identical
  restarts. `ctest -L amrex_hydro`: 27/27 here; the whole `amrex` label
  54/54 + 27 and the default build's fast suite 135/135 (gfortran live
  tests included) after the change.
- Stage E1 model: every rank still holds the full gas grid; a step copies
  its boxes in and gathers the level back with per-rank broadcasts (40 B
  per cell per half step). `hydro_driver` gained `--set` and rank-0 I/O.
- Timing on the 4-core VM (`docs/stage_e_hydro.md` §4): with 4 ranks × 1
  thread the backend is 17 % faster than the serial OpenMP code at 128³
  (10.1 vs 12.2 s), 2 × 2 is 13 % faster at 256³ (99 vs 114 s); one rank
  × 4 threads is 20 % slower (per-box fork/join); `FillBoundary` 5 % and
  the host gather 7 % at 256³. Memory: the full host state per rank (0.8
  GB at 256³) plus the level (330 B per cell): 4 ranks at 256³ exceed the
  15 GB VM, 2 ranks need 12 GB against the serial 6.3 GB.
- The spike's multi-node hydro slowdown: hypothesis (ghost volume of 16³
  boxes against the node's injection bandwidth) and the job that measures
  it on 1/2/4 nodes with box sizes 16/32/64,
  `tools/hyak/amrex_hydro_scaling.slurm` (+ `make_hydro_problem.py`);
  **the user runs it** (the Hyak thread).
- Open: the state on the level between steps (after E2–E4), GPU target
  for the hydro level, the hydro level's decomposition vs the RT level's.

## 2026-10-12, branch `claude/clever-fermi-aa4km4`: what's next, Stage E prerequisites

* Assessment for the user (`docs/PLAN.md`, "Stage E prerequisites"): MPI
  next (Stage E: hydro, PM, halos/sources, I/O distributed on AMReX; the
  RT already is), then AMR (Stage F). Before the port: the klone 256³
  timing of the Stage D backend (`tools/hyak/amrex_scaling.slurm`, the
  user's Hyak thread), the multi-node hydro slowdown of
  `docs/framework_decision.md` §3.3, and the open Stage F decisions (RT on
  AMR mixture, `photon_conserve` in `legacy_physics = 0`, cooling in dense
  gas).
* **Production default back to case B** (user): `make_production.py`
  defaults `--recomb-case B --recomb-rad 0`, so the production
  configuration runs on `rt_backend = amrex`. Case A + RR with
  `--recomb-case A --recomb-rad 1`.
* Hyak scripts: `module load` after `env.sh` failed under `set -u`
  (`build_amrex.sh` on klone: "mpicxx: command not found"); every script now
  loads modules through `radray_module_load` (`tools/hyak/env.sh`).

## Done 2026-10-12: M1's neutral source cells, `moment_inject` (branch `m1-clumps`)

**`docs/rt_moments.md` §5a.** After overlap the 512³ M1 run kept 1.5 × 10⁻³
of its hydrogen neutral (rays 1.3 × 10⁻⁴), 96 % of it in 1590 cells at
100–500 ρ̄ that the rays and OTVET ionise. At 128³ with the source map
written (`fields = ...,Luv`), most of that HI sits in the brightest source
cells. Cause: the moment RT injected the source and recombination photons
before the transport, which moves 3/4 of a cell's photons to its
neighbours in the first subcycle, so dense source cells saw too few of
their own photons. **Fix:** `moment_inject = 1` injects after the
transport (C++ default 0; `make_production.py --moment-inject` default 1).
At 128³ M1's post-overlap neutral fraction falls 30–60× to the level of
OTVET and case B rays; the moment proximity peaks rise to 2–3× the rays';
the diffuse Γ distribution does not change. New test
`restart_smoke_m1_inject`. Runs `$RUNS/m1_clumps/`. The 512³ moment runs
used 0 and are not rerun; `docs/rt_methods_512.md` §6b, §6c and §9 carry
the correction.

## Done 2026-10-11/12, branch `pm-refine` (merged with `main` 2026-10-12; PR pending)

**`pm_refine`: a PM mesh of pm_refine·N cells for Ndm up to pm_refine·N
particles (`docs/pm_refine.md`).** Goal: 1.3×10⁸ M☉ halos (20 particles)
in the 40 Mpc/h, 512³ runs.
- How it works:
  - The gas density is spread onto the fine mesh and the mesh forces are
    averaged back (adjoint) in `calc_kick` and `calc_kick_unsplit`.
  - CIC, `link_dm` and the particle mass are on the mesh.
  - `x_dm_gas = Ngas/Ndm` is now a real number (`summarize_dm`,
    `halos.cpp` halo velocities).
  - Threaded FFTW: `fftw-3.3.10-float-omp`, the same SSE2 build with
    `--enable-openmp`, which `env.sh` now selects; capped at 32 threads.
  - `make_ics.py` handles Ndm > Ngas and has `--phases-from-ndm`.
  - `make_production.py --pm-refine 2 [--data-files dm]`.
  - The checkpoint format is unchanged (the mesh has no state).
- **pm_refine = 1 is byte-identical.** Fast suite 115/125 (the known 10
  stored-RT failures); every PM and hydro stored reference passes;
  nightly 22/22 with the OpenMP FFTW library. **After the merge with
  `main` (2026-10-12, b56d4ed, and again after PR #29, c7a5f3f): fast suite 125/125 on klone** (the
  platform-dependent RT references are skipped there since `main`'s
  `platform_ref` label), including every pm_refine / FoF-cadence test;
  nightly 22/22 (jobs 41787154 and 41791863, 58–59 min).
- Tests (all pass), 16 new:
  - single-particle force against Ewald: the error is a function of the
    distance in mesh cells, so pm_refine = 2 is as good at d/2;
  - no self-force;
  - momentum conserved to 3×10⁻⁹ for every kick variant;
  - an r = 2 pancake is bit-identical to its r = 1 twin on a gas grid as
    fine as the mesh;
  - pancakes with r = 2;
  - FFT and full-loop thread-count independence;
  - restart with Ndm > N.
- Gas+DM Zel'dovich pancake (both kicks): the DM errors of pm_refine = 2
  at N equal those of the original at 2N, and the gas is as good or
  better.
- **HMF, 40 Mpc/h with 1024³ DM** (`$RUNS/hmf512_40mpc_r2`, 55 steps,
  173 GB), against `hmf256`:
  - 0.73–0.89 down to 20 particles (1.3×10⁸ M☉), with no extra falloff
    near the limit;
  - pm_refine = 1 (512³ DM) has none of these halos and 10–45% at
    10⁹–10^9.5 M☉;
  - 10 Mpc/h controls with the same DM ICs attribute the 10–25% deficit
    to the gas grid being 2× coarser than the mesh;
  - against the same-configuration 10 Mpc/h control the big box is flat
    in particle number (0.84–0.94, sample variance).
- Cost per step (no RT): ~60 s against 29 s, mostly the on-the-fly FoF on
  10⁹ particles (32 s). PM 13 s (30 s before the r > 1 speed-ups:
  parallel `link_dm`, stored mesh force, FFT cap); 174 GB peak.
- **FoF cadence** (user: "FoF should definitely run less often than
  every step"): `fof_dt_myr`.
  - The FoF runs every Δt Myr; in between the cached halos are reused,
    drifted, and the sources re-assigned at the current z.
  - The cache is in the checkpoint (version 7).
  - Default 0 = every hydro step (byte-identical).
  - Tests: `thread_full_fof_cadence`, `restart_smoke_fof_cadence`.
  - In the 512³ production schedule (39 source steps, 20–30 Myr apart
    until overlap): 20 Myr → 33 FoF calls, 30 Myr → 22, 50 → 16. That is
    only 3–12 min saved at 1024³; 20 Myr suggested.
- **Abundance matching (user decision 2026-10-11):** Sheth–Tormen refit to
  the FoF halos ≥ 1000 particles of the 40 Mpc/h 1024³ run (A = 0.3458,
  a = 0.6164, p = 0.3), `data/uv_lf/am_params_st_fof.txt`, for every box
  and pm_refine. Data/fit 0.92–1.07 at z = 5.5–8; the 512³ and 10 Mpc/h
  boxes agree; the 1024³ box LF matches the observed LFs to 0.9–1.1
  wherever halos are resolved. `am_fit`/`hmf_table --st-params`,
  `tools/analysis/fit_st_fof.py`. FoF stays every hydro step (user).
  **f_esc needs recalibrating with the new parameters.**
- Decided (user, 2026-10-11): FoF every hydro step (`fof_dt_myr` stays
  as an option, default 0); abundance matching on the ST refit (the
  R(N, z) correction is superseded).
- **Open (user):** a production run with pm_refine = 2 (f_esc to
  recalibrate with `am_params_st_fof.txt`).

## Done 2026-10-12: nightly workflow out of disk (branch `claude/new-session-knktw3`)

The first scheduled `nightly.yml` run (on `main` at d8d6029) died after
2 h 50 min with "No space left on device", leaving no log.
- **Cause:** the full hydro tier keeps every run in `build/`. A dry run
  of the generated problems gives 6.9 GB for the legacy scheme and 7.6 GB
  for unsplit (KH and the pancake ~2 GB each). The repository is private,
  so the job gets the standard private-repo runner: 2 cores and 14 GB of
  disk (from memory, not checked against GitHub's docs; the workflow now
  prints `nproc`, `free` and `df`).
- **Fix:** `run_suite.py --clean` removes each test's problem and run
  directories once it is analysed (the results JSON stays); the nightly
  full-tier CTests pass it. `nightly.yml` runs the label as three jobs,
  each on its own runner: `rt-full` (full-loop references, pool guard,
  Stromgren), `hydro-legacy` and `hydro-unsplit`. It uses ctest
  `-j$(nproc)` and a 350-minute limit per job.
- **Measured here** (4-core container; each hydro scheme pinned to 2
  cores, as on a runner): `rt-full` 8/8 in 30 min on 4 cores, peak 2.5 GB
  (`recomb_stromgren` 2.2 GB). `hydro-legacy` 7/7 in 1 h 41 min, peak
  3.1 GB. `hydro-unsplit` 7/7 in 1 h 39 min, peak 2.7 GB. Not yet run on
  GitHub.
- **Open, for the user:** the daily cadence uses roughly 4–5 runner-hours
  per night, beyond the included Actions minutes of a private repository
  on the Free and Pro plans. The hydro jobs could run weekly instead.

## Done 2026-10-12: rays, M1 and OTVET at 40 Mpc/h, 512³ (branch `rt-methods-512`)

**`docs/rt_methods_512.md`.** Three runs differing only in `rt_method`
(case A + transported recombination radiation, `c_reduced = 0.2`,
f_esc = 0.23, `redshift_per_distance`) plus rays with case B; runs
`$RUNS/reion512_rays{,_m1,_otvet,_caseB}`, comparison
`$RUNS/rt_methods_512/`, movie `docs/movies/rt_methods_512.mp4`.
- **Code:** whole-box Γ_HI columns in `movie/frames.txt` (volume mean,
  median and mean over x_HII > 0.5, f_V(Γ > 10⁻¹⁴)); field `gammaHI`
  (HI rate only) for movie frames and snapshots; `restart_keep = n`;
  `data_files = gas`; `make_production.py --fields --restart-keep
  --data-files --recomb-case --recomb-rad --recomb-rad-block`, **case A +
  RR is now its default** (user); readers of `frames.txt` accept both
  formats. Fast suite on klone 107/107 after merging `main`; `movie_slab` checks
  the Γ columns against the cubes.
- **Tools:** Γ row in `compare_movie.py`, Γ column in `stills.py`, Γ
  history/maps/PDFs in `rt_methods_runs.py`, `gamma_power.py` (P_Γ, r(k)),
  `igm_stats.py` reads `gammaHI.bin`, `rt_methods_post.slurm`
  (`SETS`, `EXTRA`, `IGM=1`), `field_cubes.slurm` (field cubes at every
  step from a checkpoint or a rerun, bit-identical).
- **Cost (one node):** rays A + RR 14.3 h / 246 GB, rays B 8.6 h / 207 GB,
  M1 4.2 h / 143 GB, OTVET 4.8 h / 146 GB; pool guard never fired.
  Recombination rays cost the rays 1.67× (up to 20× per step at z ≈ 11).
- **Results:** histories and τ_es (0.053) agree to |Δx_V| ≤ 0.017; the
  mean Γ after overlap agrees to 1–2 %; the Γ distribution does not
  (16–84 %: rays ±35 %, OTVET ±22 %, M1 ±13 → ±6 %). P_Γ at k ≲ 6 h/Mpc:
  M1 0.07–0.4× the rays and decorrelated (r = 0.1–0.6 at k ≲ 1.5), OTVET
  2–9× with r ≥ 0.95; both grow with time after overlap. ⟨F_Lyα⟩ at fixed
  ⟨Γ⟩ differs by 20–35 % between methods. RR: Γ +10 %, λ_mfp −30 %.
- **Found and fixed (2026-10-12):** the rays' block emission of recombination
  radiation moves a dense cell's own recombination photons to the block's
  emitting cell, so case A + RR makes spurious self-shielded clumps (≈ 90
  cells holding 62 % of the HI at z = 5.44; case B none). New option
  `recomb_rad_cell_onspot = τ` (default 0; `make_production.py` 1): cells
  thicker than τ at the HI edge keep their photons on the spot. 128³ test:
  residual HI −4×, history and τ_es unchanged (`docs/rt_methods_512.md` §6c).
  The 512³ case A + RR ray run predates it and is not rerun (user); its
  post-overlap residual HI, ρ_HI power and λ_mfp are flagged in the doc.
- **Open:** compare P_Γ with the user's analytic prediction; a
  substep-averaged Γ to remove the rays' sampling noise; disk — the runs'
  checkpoints (≈ 0.9 TB, incl. a saved 166 GB step-48 checkpoint of the
  ray run) can be deleted once no restart is needed.

## Done 2026-10-12: rays crossing k cells per substep (branch `rt-multicell`, PR #27)

**`docs/rt_multicell.md`.** This tests a coarse RT substep for the RT on
AMR (Stage F). Instead of sub-cycling level l 2^l times, rays cross 2^l
fine cells per coarse substep. On the uniform grid that is k cells per
substep with a k times longer substep at the same ĉ.
- **Flag:** `rt_cells_per_substep = k` (1–16, default 1 = the original,
  byte-identical; `legacy_rt = 0` and `rt_method = rays` only; `rt_driver`
  and `radhydro`; checkpoint version 6, a restart may change k). In the
  shared kernels, so the AMReX backend supports it (ghosts and local
  `Redistribute` ≥ k + 1; round-off tests at k = 2, 4).
  `deterministic_parallel` domains d ≥ 2k + 2.
- **Generalisations k forced:** the rate uses the crossing speed
  k·len_rt/dtsec_rt; the split clamp is [rmin, rmin + k]; for k > 4
  sources emit at the nside that covers the first path.
- **New options, off by default:** `photon_conserve = 1` fixes **B17**
  (new in `docs/port_notes.md`: the photon limit of `photo_gas_cell` stops
  before converging, so photon-starved front cells ionise without the
  photons; 1.8 % of the ionisations at k = 1, growing as k²).
  `ray_deposit = exposure` keeps a ray's photons for the cells further
  along its path unless absorbed. With the share deposit they jump up to
  k − 1 cells past I-fronts.
- **k = 1:** byte-identical; fast suite 117/117 (5 new tests), nightly label 22/22, AMReX label 54/54.
- **Results** (exposure + pc against k = 1): the uniform I-front at c = 1
  is within 0.04 cells for k ≤ 4 and 0.15 for k = 8. At 0.1c, k = 4 lags
  2.3 cells while the front is light-limited. `uniform_8src` fronts are
  within 0.02 cells (k ≤ 4); `lognormal_16src` x_HII,M is −1.8/−2.9/−14 %
  at k = 2/4/8. In the 64³ production configuration at 0.1c, τ_es is off
  by −0.12/−0.52/−1.6 % and T_M by 5/10/20 % early. At c = 1 (to z = 9), k = 4 matches k = 1 to 3e-5 in τ_es, Δz = 0.003 and 0.07 cells, with T_M 3 % high. The RT per hydro step falls only 1.3× at k = 2 and 4 (more live rays, k cells per ray-substep). k = 8 is 2× slower at 64³, because the deterministic domains leave one per colour.
- **What fails first:** the lagged τ at fronts moving more than one cell
  per substep (v_front > ĉ/k). The transport over k cells is accurate
  (k = 8 at 8ĉ matches k = 1 at ĉ to 1 %); the error is the longer
  substep. So on the uniform grid k > 1 equals fewer substeps and is not a
  speed option.
- **Recommendation for the AMR RT (§9, for the user to decide):** a
  mixture. Use the coarse substep with 2^l cells per substep, the exposure
  deposit and `photon_conserve`. Add photo-chemistry sub-cycling in fine
  cells, and transport sub-cycling only in patches where a front outruns
  the substep.
- **Tests:** `ctest -L multicell` (`multicell`, `multicell_exposure`,
  `photon_conserve`, the thread and restart tests, AMReX round-off).
  `tools/hyak/rt_multicell.slurm` is a 128³ repeat of the production
  comparison for the user to run.
- **Follow-ups (not done):**
  - an exposure field for Γ and a corrector pass of step 1 (§6);
  - 27 deterministic colours so that k = 8 parallelises at 64³;
  - `amrex_full64_vs_serial` (nightly) was not re-run on this branch;
  - the 128³ repeat on klone.

## Done 2026-10-11: Stage D, the ray-tracing RT on AMReX (branch `stage-d-amrex`, PR #24)

**`docs/stage_d.md`.** Uniform grid, CPU first; the serial port stays the
bit-for-bit legacy reference and the AMReX path is the new option
`rt_backend = amrex` (`legacy_rt = 0` only; AMReX optional,
`-DRADRAY_AMREX=ON`).
- **Milestone 1, one physics implementation:** the spike's ray and RT-cell
  kernels moved to `include/radray/kernels/` (`ray_kernels.hpp`,
  `rt_cell.hpp`, `healpix_hd.hpp`, `consts_hd.hpp`) and the serial port
  calls them (`Ray` is `rrk::RayP` + links; cast, advance, split, merge,
  emission, redshifting, `photo_cell`, rates, HEALPix exist once). Fast
  suite 100/100 and the nightly label unchanged, every flag combination
  byte-identical. Checkpoint version 5 (the old ray layout is read).
- **Milestone 2, the library component** `src/amrex/rt_level.cpp`: rays in a
  `ParticleContainerPureSoA`, cell FabArrays, int64 deposits with
  `SumBoundary`, 2 ghost layers with the reach assertion, `Redistribute`
  every substep, 128-bit fixed-point RT statistics; the full source model,
  thick-source deposits, the pool guard, redshifting and `post_rt` in the
  hydro loop, movie frames, checkpoints; wired into `rt_driver` and
  `radhydro` (rank 0 writes; the hydro/PM run redundantly on every rank in
  Stage D). Bitwise identical to `spike/serial rt_serial` for 1/2/4/8
  ranks × 1/2 threads × box sizes 32/16/8 (and 4 ghosts, global
  Redistribute); round-off against `rt_driver --set legacy_rt=0` on the
  five RT problems; the full loop rank- and thread-independent with
  bit-identical restarts.
- **Milestone 3, R-struct:** `ray_encoding = rstruct`, 28 B of ray data per
  ray (80 B float; 48 vs 100 B per particle with AMReX's position and id),
  dithered. Errors against the float layout: x_HII ≤ 5e-4, τ_es ≤ 5e-6,
  T_M ≤ 1e-3, 1 % rms per ionised cell, I-front ≤ 0.03 cell on
  `lognormal_16src`/`uniform_8src` (100 steps) and the 64³ production
  configuration, i.e. 10–100× inside the acceptance metrics.
- **Milestone 4, tests:** `ctest -L amrex` (43 tests: the bitwise matrix,
  the five round-off comparisons, threads, R-struct, 1 vs 2 ranks of the
  full loop, restarts; nightly `amrex_full64_vs_serial`). The 64³
  production run on the backend reproduces the serial history to 1e-5
  (τ_es = 0.05224 in both). Timing (4-core VM, 64³, two hydro steps
  from a production checkpoint): serial OpenMP 23.5 s of RT on 4 threads,
  backend 24–27 s for any of 1×4, 2×2, 4×1 ranks×threads, 74 s on one
  core (serial 87 s), R-struct +30 %; after removing the CPU atomics
  (first measurement 137 s: a shared histogram and 13 atomics per
  deposit serialised the threads). 128³ (22.5 M rays): serial 149 s,
  backend 129–152 s on 4 cores (4 ranks × 1, 32³ boxes best), 391 s on one
  core (serial 599 s), R-struct +26 % time, −36 % memory. Rank 0's peak
  memory at restart (10–11 GB at 128³) is open. The klone
  256³ job (24 ranks × 8 threads) is `tools/hyak/amrex_scaling.slurm`,
  **the user runs it**.
- **Milestone 5, GPU:** scaffolding only (CUDA CMake target, device rate
  tables, `tools/hyak/build_amrex_gpu.sh`, `amrex_gpu.slurm`); **not
  compiled** (no nvcc ≥ 12.2 here).
- Open: `recomb_rad` on the backend; the per-box serial folds and host
  copies for a GPU build; a team-level step 1 for source cells; Stage E
  (distributed hydro/PM, the §3.3 multi-node hydro slowdown).

*Last updated: 2026-10-11: slab-only movie frames `movie_slab`, the movie-reader speed fix and `POST=sbatch` (branch `movie-slab`, not merged). Before that, 2026-10-10: `redshift_per_distance` (branch `redshift-per-distance`, not merged); IGM diagnostics (branch `igm-diagnostics`); recombination radiation with the moment RT and helium lines (branch `recomb-moments`, merged into `main`); status assessment and the Hyak production and Stage D prompts ("Threads to launch"); B16 fix `ray_rotate` merged into `main` (branch `claude/exciting-clarke-zn62ih`; production uses `ray_rotate = 3`). The same day: the framework-spike Hyak CPU and GPU benchmarks (PR #19, branch `spike-hyak-bench`: AMReX approved, no ACCESS proposal for now), case A recombination and recombination radiation (branch `recomb-radiation`), and moment RT (M1, OTVET) with its 256³ klone comparison (branch `claude/elegant-bell-rqd3l1`). Before that, 2026-10-09: the survey (PR #16), `hyak-setup` (PR #17), `movies` (PR #12) and `production-run` (PR #7).*

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
| 5. production run | **done: `$RUNS/reion256_f0.2`, f_esc = 0.2** (`docs/production_256.md`): 1 h 30 min on a whole node, reionization ends at z = 5.6–5.7, τ_es = 0.057; movie `docs/movies/production_256_c01.mp4`. The f_esc = 0.28 run (`reion256_prod`, ends z = 6.44, τ_es = 0.063) 0.15 (`reion256_f0.15`, x_V = 0.87 at z = 5.5) and 0.1 (`reion256_f0.1`, does not finish) bracket it |
| 6. c convergence | done: `docs/production_run.md` §9; 0.1c delays the end of reionization by Δz = 0.4 against full c (0.3c: 0.14), the midpoint by 0.15, τ_es 2.5 % low |
| 7. ray-precision check | done: `docs/ray_compression.md` §6.2 (`$RUNS/p128/prec`, job array 41705715) |
| 8. STATUS, PLAN, PR | this entry; PR opened for `production-256` |


## Done 2026-10-11, branch `movie-slab` (merged into `main`)

**Slab-only movie frames, faster movie tools and post-processing as a
separate job (`docs/movies.md`, "Slab frames").**
- **Problem.** A movie frame was a full cube per field: about 160 GB per
  512³ run at Δz = 0.1. The tools read only a slab of it.
- **New option.** `movie_slab = t` (default 0 = the full cube, byte-identical
  to before) and `movie_slab_axes` (default `xyz`).
  - Each frame field becomes `movie/frame_NNNNN_<field>.slab`: one block
    per axis, holding cells N/2 − t/2 … N/2 − t/2 + t − 1 along that axis,
    float32, Fortran order.
  - The geometry is in `movie/slab.txt`. A restart with a different
    geometry stops with an error.
  - `frames.txt` is unchanged: x_HII,V and x_HII,M are whole-box sums.
  - Disk with t = 8, at 512³: 25 MB per field instead of 537 MB, so 75 MB
    instead of 1.6 GB per frame of fHI, temp and rhog.
  - `make_production.py --movie-slab 8` sets it, recommended for N ≥ 256;
    the default 0 writes no line.
- **Tools.** A shared reader, `tools/movies/movie_frames.py`, serves
  `make_movie.py`, `compare_movie.py`, `stills.py` and `compare_stills.py`.
  They read both kinds of frame.
  - With slabs, the slice must lie in the stored slab; otherwise the tool
    says which cells are stored.
  - `reion.py` and the other tools that read only `frames.txt` need nothing.
  - `rt_methods_runs.py` needs whole cubes for its field statistics. With
    slab frames it falls back to the zdata1 field outputs, and says so.
- **Speed: whole-cube reads fixed.** The old slab reader called `np.take` on
  the memory-mapped cube. That copies the whole file, so each 512³ field was
  read and transposed in full. This made `make_movie.py` take 1–2.5 h after
  every 512³ leg.
  - The new reader slices instead. It reads the same values, and stills are
    byte-identical to before.
  - Synthetic 512³, 6 frames: 106 s → 1.6 s (slabs) for files not yet in
    the page cache, and 9.5 s → 5.4 s for cached cubes. `stills.py`:
    46 s → 5 s.
- **Post-processing.** `tools/hyak/reion_post.slurm` now holds the
  post-processing that `reion_production.slurm` ran after every
  successful leg.
  - `POST=1` (default) runs it in the same job, as before.
  - `POST=0` skips it.
  - `POST=sbatch` submits it as a separate 16-core job (`POST_CPUS`,
    `POST_MEM`, `POST_OUT`), so the big node is released at once.
  - `POST=sbatch` has not been tried in a production job yet.
- **Tests.**
  - New: `movie_slab` (smoke full loop: slabs against the cubes, frames.txt,
    reion.py, identical make_movie and stills images from cubes and slabs,
    the error messages, the planar style) and `restart_smoke_det_movie_slab`.
  - Defaults: the main build (195426c) and this branch give byte-identical
    outputs on the smoke problem with movie frames, 337 files including 44
    frame files, for both legacy settings.
  - Fast suite on klone: 101/111 pass. The 10 failures are the known
    platform ones (`ref_*lognormal_16src*` ×6,
    `ray_dump_*_lognormal_16src_smallbuf` ×2, `ref_c01_*` ×2).

## Done 2026-10-10, branch `redshift-per-distance` (merged into `main`)

**Redshift per distance travelled with a reduced speed of light
(`docs/production_run.md` §1, "Redshifting per time or per distance").**
- Problem: `ray_redshifting` and `moment_redshift` scale frequencies by
  fscale = (a − da)/a once per hydro step, per unit time. At ĉ < c the
  radiation then loses c/ĉ times more energy per distance than real
  photons. At ĉ = 0.1c and z ≈ 5.5 the 15.7 eV group falls below 13.6 eV
  after ≈ 6–7 pMpc (≈ 70 pMpc at c), comparable to λ_mfp after overlap;
  at 0.01c after ≈ 0.7 pMpc. Γ_HI is biased low for a given emissivity.
- New option `redshift_per_distance` (params.txt / `--set`, `radhydro`
  and `rt_driver`; default 0 = the original code path):
  fscale_eff = fscale^s, s = min(nstep_rt·dl_rt/dl_c, 1), the distance the
  radiation moved in the step (actual `nstep_rt`, after nint,
  `rt_substep_carry` and the ≥ 1 guard) over the light distance of
  `calc_timestep`. No effect with `c_reduced ≥ 1`. Rays and M1/OTVET.
  Independent of `legacy_physics` (the running production runs use
  `legacy_physics = 0` and are unchanged). Not in the checkpoint: can be
  changed at a restart like `c_reduced`. Recorded in the `rt_driver` log,
  the `radhydro` summary line (when set) and the `ray_dump` json;
  `make_production.py` writes it; **production default `= 1` since
  2026-10-11 (user decision, `docs/PLAN.md`)**, `--redshift-per-distance 0`
  for the original.
- Tests: `redshift_distance` (unit) and `redshift_distance_smoke`
  (`rt_driver`, rays and m1; physics label), both pass. Fast suite on
  klone (gcc 13): 99/109 pass; the 10 failures are the same ones as on
  `main` without the change (`ref_*lognormal_16src*` ×6,
  `ray_dump_*_lognormal_16src_smallbuf` ×2, `ref_c01_*` ×2: mismatches with
  the stored references on this platform, e.g. last-digit rt_stats values).
- **Open (user):** whether production runs should use it. Not turned on
  anywhere.

## Done 2026-10-10, branch `igm-diagnostics` (merged into `main`)

**IGM diagnostics for the f_esc calibration (`docs/igm_diagnostics.md`).**
- `tools/analysis/igm_stats.py` measures a run's snapshots (z ≤ 7.1) to
  JSON:
  - Γ_HI from the `gamma` field (s⁻¹, last RT substep): volume average,
    ionized-cell averages, transmission-weighted;
  - the 912 Å mean free path from random, ionized and source starts
    (periodic random-direction sightlines, 1/e and exponential fit);
  - the Lyα mean flux in redshift space (u_gas velocities, Doppler,
    top-hat cells) and the τ rescaling A to Bosman+22;
  - ṅ_ion at every hydro step from `stats.txt` `Luv`:
    ṅ = 0.2·10^25.2·1.0001·Luv. **No C++ output is needed.**
- `tools/analysis/igm_plot.py` makes the 4-panel comparison with
  `data/igm/igm_z5_6.txt` (brought in from `fesc-z`) and a table.
  `tools/hyak/igm_diag.slurm` runs both on a node.
- 256³/10 Mpc/h runs: figure `docs/figures/igm_256.png`, outputs in
  `runs/igm_diag/`.
  - ⟨F⟩ and Γ both point to f_esc ≈ 0.23–0.25 at z = 5.4–6.
  - ṅ(f_esc = 0.2) = 10^51.0 at z = 5.42, within the Becker+21 errors.
  - At fixed Γ the simulated λ_mfp is ~2.7× longer than observed
    (unresolved absorbers, box 2.2 pMpc < λ), so λ cannot calibrate these
    runs.
  - The RSLA (ĉ = 0.1c: delay λ/ĉ, enhanced ray redshifting) probably
    lowers Γ/ṅ.
- Open:
  - run the tools on the 512³/40 Mpc/h run when it has z < 6.5 outputs
    (the memory and time were checked on a synthetic 512³ snapshot);
  - check the Gaikwad+23 Γ definition against the paper.

## Done, merged into `main` 2026-10-11 (branch `fesc-z`)

**Redshift-dependent escape fraction and the emissivity of the runs.**
- `f_esc_beta`, `f_esc_z0`, `f_esc_max` (abundance-matched sources):
  f_esc(z) = min(f_esc_max, f_esc ((1+z)/(1+f_esc_z0))^β), the power law of
  D'Aloisio et al. (declining emissivity near the end of reionization).
  β = 0 (default) is the constant f_esc exactly. Unit test in
  `test_abundance_matching`; `make_production.py --f-esc-beta`.
- `tools/analysis/emissivity.py` (on `main`): escaped ṅ_ion(z) of the runs
  from FoF catalogs (hmf256: same ICs), vs the Lopez et al. median read by
  eye from D'Aloisio et al. Fig. 1. `reion256_prod` (f_esc = 0.28) matches
  that median at z ≳ 9 but is 2–2.5× above it at z < 7.5. Matching it with
  these halos needs f_esc ≈ 0.123 ((1+z)/7)^1.31 (0.12 at z = 6, 0.22 at
  z = 10), a milder decline than their β = 2.6–3, because the
  abundance-matched intrinsic emissivity already flattens at z < 7.
- The z = 5–6 data (`data/igm/igm_z5_6.txt`, `emissivity.py --obs`) need
  more than that law at z ≈ 5: a constant f_esc ≈ 0.15–0.2 with the hmf256
  halos. The run with that law (`reion256_fz131`) was therefore cancelled,
  and calibration moved to the 40 Mpc/h 512³ runs (constant f_esc; see
  the IGM diagnostics entry above). No run uses `f_esc_beta` yet.
- Production runs now end at z = 5 (user, 2026-10-10: Lyα forest
  comparison); `make_production.py` default.

## Done, merged into `main` 2026-10-10 (branch `recomb-moments`)

**Recombination radiation with the moment RT, the photon energy, and the
helium lines (`docs/recomb_radiation.md` §2.1–2.3, §3, §4.3).**
- **Moment RT:** `recomb_rad = 1` now works with `rt_method = m1 | otvet`.
  - Every cell adds its recombination photons to its own photon density
    (isotropic, zero flux), with no blocks or pool, and case A everywhere.
  - Strömgren: R/R_B = 1.007 (M1) and 1.006 (OTVET).
  - 256³ at 0.1c (`reion256_c01_{m1,otvet}_caseA_rr`):
    - history within Δz ≤ 0.03 of the same method's case B, τ_es
      unchanged;
    - ionised gas 3–5% hotter before overlap and 18% hotter at z = 6 (rays:
      14%);
    - about 20% of the recombination photons are lost to redshift;
    - cost +4% RT time (rays: +43%).
  - Movie: `docs/movies/recomb_256_c01_moments.mp4`.
- **Photon energy:** the mean excess over threshold now comes from the
  Milne relation, integrated with the code's cross-sections, replacing kT.
  - H: 0.96 kT at 10⁴ K and 0.88 kT at 4×10⁴ K.
  - This matches the extra case-A cooling of the Hui & Gnedin fits to 2–3%
    of kT, so the energy is consistent with `legacy_physics = 0`.
- **`recomb_rad_lines = 1`:** H-ionising helium recombination photons:
  - HeI 19.8 and 21.2 eV lines and two-photon decay: 0.97 per HeII
    recombination;
  - HeII Lyα, two-photon decay (1.42 photons per decay) and Balmer
    continuum.
  - Strömgren with M1: +3.2% in R, as expected for ~8% extra photons.
  - **The branching fractions are from memory (Osterbrock & Ferland ch. 4;
    Nussbaumer & Schmutz 1984 two-photon fit) and need checking against
    the sources.**
- Fast suite: 94/104, the same 10 known failures. New tests:
  - fast: `thread_rt_recomb_m1_lognormal_16src` and
    `restart_smoke_recomb_otvet`, both with the lines;
  - nightly: `recomb_stromgren_{m1,otvet}`.
- Open: no 256³ run with the lines (user, 2026-10-10: helium is unimportant enough that one is not needed); helium branching fractions and the two-photon fit still to check against Osterbrock & Ferland (2006) and Nussbaumer & Schmutz (1984).

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
  With `ray_rotate = 3` (2026-10-11): the B16 holes in case B are gone,
  and the block-2 excess of A + RR falls from 5.0% to 3.7%; the rest is
  the wider edge (`docs/recomb_radiation.md` §3).
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
  Convergence checks done 2026-10-10 (`docs/hydro_unsplit.md` §6,
  `tests/hydro_suite/convergence.py`): the Sedov bias is a first-order
  0.6-cell offset, +0.63 % at 256³ (legacy −1.1 %, its energy leak grows
  with N); the cold-gas heating converges with the step, and in the
  production box (10 Mpc/h, 128³ and 256³, z ≥ 3) the original step
  already keeps the underdense gas on its adiabat. Fit for 256³ production.
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
- Movie, test movie `docs/movies/reion_test_128.mp4` (the thread's own
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
stills and movie calls); `tools/movie/make_movie.py` (sic) of `production-run`
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
  references are platform-dependent (libm). Decided 2026-10-12: skip them
  on klone (label `platform_ref`, left out by `tools/hyak/build_hyak.sh`);
  CI keeps running them.
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

## Threads to launch (2026-10-09, revised 2026-10-10)

*Assessment of 2026-10-10 (session `lucid-heisenberg`): PR #19 merged.
The production-run thread (item 2) was started on Hyak in a session that
did not finish; nothing from it reached the repository (no branch, no 256³
mass function, no 0.1c calibration). Every 256³ run on klone so far used
the global SFR model and no frame rotation, so none is a science run. The
B16 fix `ray_rotate` (item 5) sits on the unmerged branch
`claude/exciting-clarke-zn62ih` and now conflicts with `main` in five
files; production needs it.*

On Hyak (`claude remote-control` in tmux, or batch jobs):

1. ~~Merge and validate~~ Done 2026-10-09 (PRs #16, #17, #12, #7). Still
   to do on klone: the fast and nightly suites on the merged `main`
   (expect the 10 known RT stored-reference failures of `docs/cluster.md` §3
   and `docs/rt_moments.md` §6.1).
2. **Production run at 256³, 0.1c, abundance-matched sources.** Opening
   prompt for the Hyak session (it also covers item 4):

   > Read CLAUDE.md, docs/PLAN.md, docs/STATUS.md and docs/ROADMAP.md.
   > You are on UW Hyak klone (`docs/cluster.md`): account `astro`,
   > partition `cpu-g2`, `source tools/hyak/env.sh`, build with
   > `tools/hyak/build_hyak.sh` on a compute node, runs under
   > `RUNS=/gscratch/astro/$USER/radhydro/runs` (the older runs of
   > `docs/cluster.md` §7 may be under `/gscratch/astro/$USER/runs`).
   > Work on branch `production-256` off `main`; commit often and push.
   >
   > Decisions already taken: 256³, 10 Mpc/h, `c_reduced = 0.1`,
   > `nray_max_factor = 16`, `legacy_rt = legacy_physics = 0`,
   > `deterministic_parallel`, `rt_substep_carry`, `ray_rotate = 3`,
   > `sfr_model = abundance`, log ξ_ion = 25.38, f_esc calibrated to
   > τ_es = 0.054. `hydro_scheme`: HYDRO_SCHEME_DECISION (`unsplit`, which
   > fixes the cold-gas overheating of `docs/hydro_tests.md`, or `legacy`,
   > which every 256³ run so far used); `make_production.py` does not set
   > it, so add `hydro_scheme = ...` to params.txt or `--hydro-scheme` to
   > the script.
   >
   > 0. **Recover the earlier attempt.** `sacct -S 2026-10-09 -X
   >    --format=JobID,JobName%30,State,Elapsed,NodeList` and the run
   >    directories under both roots (`DONE`, `run/reion.txt`,
   >    `run/restart/`). Reuse any finished no-source 256³ run,
   >    calibration run or production run that matches the settings
   >    above; list what you found in docs/STATUS.md.
   > 1. **`ray_rotate` must be on `main`.** If `git grep ray_rotate
   >    include/radray/config.hpp` on `main` finds nothing, merge `main`
   >    into `claude/exciting-clarke-zn62ih` (conflicts in
   >    apps/rt_driver.cpp, include/radray/io.hpp, tests/CMakeLists.txt,
   >    docs/STATUS.md, docs/PLAN.md: keep both sides, the branch adds
   >    the `ray_rotate*` options, a `rays` unit test and a `ray_rotate`
   >    physics test), rebuild, run the fast suite (only the known
   >    stored-reference failures allowed), push, open a PR and merge it.
   >    Then branch `production-256` from the new `main`.
   > 2. **Build and test** `main` on a node: `sbatch
   >    tools/hyak/build_hyak.sh`, then the nightly label
   >    (`CMAKE_EXTRA=-DRADRAY_NIGHTLY_TESTS=ON CTEST_LABELS=nightly`).
   >    Record the results in docs/cluster.md §3.
   > 3. **256³ mass function** (`docs/hmf_validation.md`): a no-source
   >    256³ run (`tvir_min = 1e30`, `write_catalogs`, 96 threads,
   >    minutes), `hmf_compare.py` at z = 5.5, 6.25, 8, 10, 12.6 against
   >    Sheth–Tormen, then `fit_fof_hmf.py` → `data/hmf/fof_256_10mpc.txt`
   >    and `build/am_fit --lf data/uv_lf/lf_params.txt --hmf st --sim-hmf
   >    data/hmf/fof_256_10mpc.txt` → `data/uv_lf/am_params_256_10mpc.txt`
   >    (`docs/abundance_matching.md` §3c). Add the 256³ figure and table
   >    to docs/hmf_validation.md.
   > 4. **f_esc at 0.1c.** `tools/production/calibrate_fesc.py $RUNS/cal128
   >    --N 128 --am-params data/uv_lf/am_params_128_10mpc.txt --fesc 0.2
   >    0.35 0.5 0.8 --threads 96 --extra c_reduced=0.1 data_files=dm
   >    fields_at_data1=0 hydro_scheme=<decision>` (later lines of
   >    params.txt override earlier ones, so `--extra` wins over the
   >    script's 0.01c). About 16 min per run on 96 threads.
   >    Then one 256³ run at the interpolated f_esc with the 256³ L(M) to
   >    confirm τ_es (1.5 h on a whole node). Write the table into
   >    docs/production_run.md §8 and the chosen f_esc into docs/PLAN.md.
   > 5. **The production run:** `make_production.py $RUNS/reion256_prod
   >    --N 256 --c-reduced 0.1 --am-params
   >    data/uv_lf/am_params_256_10mpc.txt --f-esc <calibrated> --xi-ion
   >    25.38 --movie-dz 0.025 --movie-fields fHI,temp,rhog
   >    --restart-every 4`, the hydro-scheme line, then `sbatch
   >    --exclusive -c 192 --mem=256G --time=1-00:00:00
   >    --export=ALL,PROB=$RUNS/reion256_prod,THREADS=192
   >    tools/hyak/reion_production.slurm`. The job writes `run/reion.txt`,
   >    stills and `run/movie/cosmo.mp4`. Put the history table, τ_es,
   >    stills and a compressed movie (`docs/movies/`, a few MB) into a
   >    new docs/production_256.md, and compare the history with the
   >    global-model `reion256_c01` run of `docs/cluster.md` §7.
   > 6. **c convergence:** one 128³ run at `c_reduced = 0.3` with the same
   >    sources as the 128³ calibration run at the chosen f_esc (about 50
   >    min on 96 threads), and if time allows one at 1.0 (hours). Compare
   >    x_HII(z) and τ_es with the 0.1c run (`tools/analysis/reion_history.py`).
   > 7. **Ray-precision check** (docs/ray_compression.md §9, the one unrun
   >    item): `make_production.py $RUNS/p128 --N 128 --movie-dz 0` then
   >    `sbatch --export=ALL,ROOT=$RUNS/p128
   >    tools/production/hyak_precision.sbatch` (a 6-task array, 24
   >    threads each, so it runs alongside the jobs above), then
   >    `precision_metrics.py` and `slices.py` as the script's header says.
   >    Record the result in docs/ray_compression.md §6/§9.
   > 8. At the end update docs/STATUS.md (this list), docs/PLAN.md (f_esc,
   >    hydro_scheme as used) and push; open a PR for `production-256`.

3. ~~Framework-spike benchmarks~~ Done 2026-10-10 (PR #19, branch
   `spike-hyak-bench`; `docs/framework_decision.md` §3.2–4). AMReX
   approved and no ACCESS proposal for now (user, 2026-10-10).
4. ~~Ray-compression 128³ precision check~~ Folded into the item 2
   prompt (step 7); `tools/production/hyak_precision.sbatch` is a klone
   job.

In the cloud (branch + PR each):

5. **B16 fix behind a flag**: done on branch `claude/exciting-clarke-zn62ih`
   (`ray_rotate = 1|2|3`, production uses 3, user decision 2026-10-10),
   **not merged**; merging now conflicts with `main` in five files. Merge
   it before the production run (step 1 of the item 2 prompt, or in the
   cloud).
6. **Stage D: sorted-ray-array RT in AMReX** on a uniform grid, with the
   R-struct encoding (8 dithered photon bits), validated against
   `rt_driver --set legacy_rt=0` as the spike was. AMReX is confirmed
   (2026-10-10), so this can start. Opening prompt (2026-10-10):

   > Read CLAUDE.md, docs/PLAN.md, docs/STATUS.md, docs/ROADMAP.md, then
   > docs/framework_decision.md (§1, §2, §6 Stage D, §8),
   > docs/design_parallel_and_encoding.md (§0.1, §2, §3, §4.2) and
   > docs/ray_compression.md (the R-struct recommendation). This is
   > Stage D of the roadmap: the ray-tracing RT on AMReX, uniform grid,
   > no refinement, CPU first. Branch `stage-d-amrex` off `main`; commit
   > often and push; one PR at the end (or one per milestone below).
   >
   > Starting point: `spike/` (`spike/README.md`, `spike/amrex/README.md`).
   > `spike/amrex/rt_amrex.cpp` already does the whole RT substep on
   > AMReX (sorted rays, two-pass emission/split/merge, canonical merge
   > order, int64 deposits with SumBoundary, Redistribute) and is bitwise
   > identical to `spike/serial/rt_serial` for 1–8 ranks and to round-off
   > against `rt_driver --set legacy_rt=0` on the five RT problems. Build
   > AMReX 26.10 with `-DSPIKE_FETCH_AMREX=ON` (MPI, OpenMP, 3D, particles
   > in single precision) as the spike did; cache the install.
   >
   > Ground rules: the serial port stays the bit-for-bit legacy reference
   > and every stored-reference, live-Fortran, restart and thread test
   > must keep passing unchanged. The AMReX path is a new option
   > (`rt_backend = amrex`, default the current code), implemented only
   > for `legacy_rt = 0` (sorted arrays cannot reproduce the zombie-ray
   > thresholds, design doc §0.1). AMReX is an optional dependency:
   > `-DRADRAY_AMREX=ON`, found with `find_package` or fetched; the
   > default build and the fast suite must not need it.
   >
   > Milestones, each validated before the next:
   > 1. **One physics implementation.** Move the ray and RT-cell kernels of
   >    `spike/common` (`ray_kernels.hpp`, `rt_cell.hpp`, `healpix_hd.hpp`,
   >    `consts_hd.hpp`) to `include/radray/kernels/`, as PR #15 did for
   >    `hydro_mh.hpp`, and make `src/ray_tracing.cpp` and
   >    `src/rt_tools.cpp` call them, so cast, advance, split, merge,
   >    emission, redshifting and `photo_cell` exist once. Bit-identical
   >    results in every flag combination (the whole fast suite, then the
   >    nightly label). Carry the serial code's newer physics into the
   >    kernels where it belongs: `ray_rotate` (if merged by then; it is
   >    on `claude/exciting-clarke-zn62ih`), the recombination-radiation
   >    emission (`recomb_rad`), and the `rt_substep_carry` / `ra_max`
   >    bookkeeping. Leave the moment solvers (`rt_method = m1|otvet`) on
   >    the serial path.
   > 2. **Library component.** Turn `rt_amrex` into
   >    `src/amrex/rt_level.cpp` (+ header) working on one level: a
   >    `ParticleContainerPureSoA` of rays, cell `FabArray`s, int64
   >    deposits with `SumBoundary`, 2 ghost layers with the assertion on
   >    the ray reach, `Redistribute` after every substep. Add what the
   >    prototype lacks: the full source model (`halo_sources`, including
   >    `sfr_model = abundance` and the local-deposit branch for thick
   >    source cells), redshifting and `post_rt` inside the hydro loop,
   >    the ray pool guard, movie frames, and checkpoints (a new section
   >    of the checkpoint format, version 5). Wire it into `rt_driver` and
   >    `radhydro` behind `rt_backend = amrex`.
   > 3. **Encoding.** Store the rays as the R-struct of
   >    docs/ray_compression.md (direction as HEALPix index + quantised
   >    offset, position as an in-cell offset, float photon total + 8-bit
   >    dithered log-quantised spectral shape, never-merged flag), decoded
   >    to float/double in registers. Keep a `ray_encoding = float|rstruct`
   >    switch so the float layout remains for bitwise validation; measure
   >    bytes per ray and the x_HII(z), τ_es and T errors of the encoded
   >    run against the float one on `lognormal_16src` and the 64³
   >    production configuration (the acceptance metrics of
   >    docs/ROADMAP.md: x_HII within 1–2 %, τ_es within 1 %).
   > 4. **Validation and tests** (CTest label `amrex`, skipped when AMReX
   >    is off): bitwise against `rt_serial` for 1, 2, 4, 8 ranks × 1, 2
   >    threads and box sizes 16³/32³ (the spike's `validate.sh` matrix);
   >    round-off against `rt_driver --set legacy_rt=0` on the five RT
   >    problems with `tools/compare/compare.py`; rank- and thread-count
   >    independence; bit-identical restarts; a full 64³ `radhydro` run
   >    compared with the serial one (history, τ_es, per-cell fields).
   >    Timing: `rt_backend = amrex` against the serial OpenMP code on the
   >    4-core VM at 64³ and 128³ (`tools/production/timing_table.py`),
   >    and a klone job script in `tools/hyak/` for 24 ranks × 8 threads
   >    at 256³ to compare with `docs/cluster.md` §6 (the user runs it).
   > 5. **GPU, last and optional:** add the CUDA target
   >    (`-DAMReX_GPU_BACKEND=CUDA`, nvcc ≥ 12.2, `--fmad=false`), rate
   >    tables in device memory, a team-level step-1 for source cells, and
   >    a klone `ckpt` job script (`spike/hyak/build_gpu.sh` is the start).
   >    Compile-only is acceptable here; the A100/H200 runs are the user's.
   >    Record the ratio against one node in docs/framework_decision.md §4
   >    (PLAN.md: ≥ 5× on H200 reopens the ACCESS question).
   >
   > Known traps: AMReX's `ParallelCopy`/`FillBoundary` need `operator+=`,
   > so struct FABs do not work (use component FABs); `ParallelFor` is
   > serial on CPU without MFIter tiling + OpenMP; particle precision is
   > build-wide (SINGLE, the original's `real(4)` for DM); the merge must
   > use the canonical order (cone, payload hash, bytes) or results depend
   > on the rank layout. Don't fix the AMReX multi-node hydro slowdown of
   > `docs/framework_decision.md` §3.3 here (Stage E), but note anything
   > you learn about it.
   >
   > Write `docs/stage_d.md` (design as built, validation tables,
   > timings, encoding errors, what is left for Stage E), update
   > docs/STATUS.md, docs/PLAN.md (decisions taken) and CLAUDE.md's layout
   > table, and open the PR.
7. **Plane-parallel / RT-domain mode** (`docs/plane_parallel.md`).
8. Smaller: ~~the stored-RT-reference platform policy~~ decided
   2026-10-12 (skip on klone, label `platform_ref`); the 2 % growth
   deficit; ~~the unsplit open items (Sedov radius convergence, cold-gas
   time-step sensitivity at 128³)~~ done (`docs/hydro_unsplit.md` §6,
   both converge); the AMReX multi-node hydro anomaly
   (`docs/framework_decision.md` §3.3).

**Decisions needed from the user:** none open from the merge thread.
Decided 2026-10-11/12 (`docs/PLAN.md`): `hydro_scheme = unsplit` and
`redshift_per_distance = 1` for production; M1/OTVET for comparison only,
with our explicit OTVET kept as it is; the stored RT references that
differ on klone are skipped there (label `platform_ref`).

## Next

1. ~~Halo mass function at production resolution~~ Done at 128³
   (`docs/hmf_validation.md`); redo at 256³ before the production run.
2. ~~Abundance-matched sources~~ Done (`sfr_model = abundance`,
   `docs/abundance_matching.md`). **Open: recalibrate f_esc at 0.1c and
   256³** (`docs/production_run.md` §8, `calibrate_fesc.py`); the 0.01c
   value is 0.5.
3. ~~Production reionisation run and movie~~ **Done 2026-10-10**
   (`docs/production_256.md`: 256³, 0.1c, abundance-matched sources with
   f_esc = 0.2, `hydro_scheme = unsplit`, `ray_rotate = 3`; reionization
   ends at z = 5.6–5.7, τ_es = 0.057, 1.5 h on a whole node). The table
   below was the plan.

   | Setting | Value |
   |---|---|
   | Box | 10 Mpc/h |
   | c_reduced | 0.1 (0.01 delays reionisation by Δz ≈ 2–2.5, `docs/cluster.md` §7); convergence check with one 0.3c or 1c run at 128³ |
   | nray_max_factor | 16 |
   | Ngas = Ndm | 256³ (1.5 h at 0.1c with a movie on a whole node) |
   | Flags | legacy_rt = legacy_physics = 0, deterministic_parallel, rt_substep_carry, ray_rotate = 3 |
   | Sources | sfr_model = abundance, log ξ_ion = 25.38, f_esc = 0.2 (end of reionization at z ≈ 5.5, user decision) |
   | hydro_scheme | `unsplit` (production-256 thread; `make_production.py --hydro-scheme`) |
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
   papers). ~~Open: convergence check of the Sedov radius bias and the
   cold-gas time-step sensitivity at 128³~~ done 2026-10-10
   (`docs/hydro_unsplit.md` §6): Sedov +0.63 % at 256³ (first order);
   cold gas converged with the original step in the production box to
   z = 3 at 128³ and 256³. No obstacle to `unsplit` at 256³.
7. **Plane-parallel / RT-domain mode** (`docs/plane_parallel.md`).
8. After the production run, follow `docs/ROADMAP.md`:
   - ray compression study: done (`docs/ray_compression.md`); 64³
     production-configuration check done, 128³ job script ready for Hyak;
   - framework spike (AMReX vs Parthenon): done, AMReX recommended
     (`docs/framework_decision.md`), **approved 2026-10-10**; Hyak CPU and GPU
     benchmarks done (H200 ≈ 5 nodes on the RT; no ACCESS proposal yet);
   - ray data-structure rewrite: done (Stage D, `docs/stage_d.md`);
   - MPI: done through the AMReX backends (Stage E, `docs/stage_e_*.md`);
   - AMR for gas, then RT: built (Stage F1–F3, `docs/stage_f.md`); what
     is left is item 9;
   - multigrid/tree gravity (Stage G);
   - adaptive RT.
9. **Before AMR is ready for science** (2026-10-15, after PR #41):
   - **Gravity on the refined levels (Stage G).** Since `pm_gas_levels`
     (below) every level gets the force of the uniform PM mesh over its
     own cells; the force resolution is still the mesh's everywhere.
     Multigrid Poisson on the hierarchy (adaptive resolution, no huge
     uniform mesh) remains Stage G, TreePM optional.
   - **Hydro validation with refinement (F2):** the pancake and KH cases,
     the full tier, and the composite of the finest data for comparisons
     at the matched finest resolution.
   - **RT validation with refinement (F3):** `uniform_8src` with a refined
     region, the Iliev tests with refinement, a 128³
     production-configuration run with two levels against uniform 128³
     and 256³ (x_HII(z), τ_es, T_M, front positions). Only the
     one-source I-front is done (`amrex_rt_amr_ifront`). Scripts for the
     128³ runs and a 512³ cost probe (0 vs 2 levels on 5 nodes) are
     ready, not yet run (`docs/stage_f.md` §3, "Next").
   - **pm_refine and the refined gas: done 2026-10-15** (`pm_gas_levels
     = 1`, default; `docs/stage_f.md` §4): the PM deposit takes the gas
     of the finest level covering each mesh cell and every level is
     kicked on its own cells, so with `pm_refine = 4` both levels of a
     512³ run get gravity at the 2048³ mesh. A 16³ run refined everywhere
     now matches the uniform 32³ run on the same mesh to 0.3 % (gas) and
     1e-4 cells (particles), against 9 % and 8e-3 with the level-0
     coupling. `make_production.py --ndm` allows 1024³ particles on the
     2048³ mesh. To measure on klone: the 2048³ mesh's cost and memory.
   - **Production:** the sink's production threshold (`sf_overdensity`,
     "a few thousand"), a first refined production run, profiling the
     flux-register and FillPatch overheads (`reflux`, `fill` timers), the
     distributed FoF and P(k) (rank 0 today, the bottleneck at 1024³
     particles).
   - **Limits:** only the AMReX backends refine; M1/OTVET exist only in
     the serial RT (`rt_backend = amrex` needs `rt_method = rays`), so
     with refinement they could act on level 0 at best, untested and not
     refused; the coarse substep lags the fine-level
     fronts (1.5 % with the whole box on level 1) and fails light-limited
     fronts at k = 4 (two levels at 0.1c, `docs/rt_multicell.md` §7–8).

## Open questions for the user

- **The unified RT-methods suite** (`docs/PLAN.md`, "Unified RT-methods
  suite", 2026-10-17): (a) ~~ĉ~~ decided 2026-10-17 and revised
  2026-10-18 (PR #46): 0.1c / 0.2c / 0.3c / 0.5c switched at x_HII,V =
  0.25, 0.5 and 0.75; (b) case B for the whole suite, or recombination radiation
  ported to the backend first; (c) the f_esc law once the 512³ run
  reports; (d) the uniform-grid suite before the moment methods have AMR;
  (e) FLD proper alone, or P1 as well; (f) the paper for the FLD limiter.

- ~~RT on AMR: mixture, pure time sub-cycling or a pure coarse substep;
  whether `photon_conserve = 1` (B17) becomes part of `legacy_physics = 0`.~~
  Decided 2026-10-14 (`docs/PLAN.md`, "Stage F decisions"): a pure coarse
  substep; B17 is part of `legacy_physics = 0`.

- Original `Halo/` catalogs, `sfr.txt`, or the Trac & Cen halo finder. These
  would pin down the original source model; it is reconstructed, see
  `docs/halos_and_sources.md`.
- Whether to keep the original's default `nray_max_factor = 64` in legacy
  mode. Production uses 16.

## Known limits

- **MPI only through the AMReX backends** (`rt_backend = amrex`, `pm_backend = amrex`; the hydro, halo finder and sources still run redundantly on every rank). GPU is pending. On one Genoa node the RT reaches 94× on 192 threads at 256³ and about 50× at 128³ (`docs/cluster.md` §6).
- **RT dominates the cost.** Hydro+PM take under 1%. At full c, z=12→6 at
  128³ would take days on one core, hence `c_reduced`.
- **Network access from cloud sessions.** arXiv, ADS and journal sites are
  unreachable.
