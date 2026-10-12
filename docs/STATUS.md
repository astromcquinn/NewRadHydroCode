# Project status

*Last updated: 2026-10-12: the nightly workflow fixed (branch `claude/new-session-knktw3`, below). Earlier that day: rays crossing k cells per substep, `rt_cells_per_substep` with `photon_conserve` (B17) and `ray_deposit = exposure` (branch `rt-multicell`, PR #27, `docs/rt_multicell.md`, merged into `main`). Before that, 2026-10-11: production runs use `hydro_scheme = unsplit` and `redshift_per_distance = 1` (user decisions, now the `make_production.py` defaults); Stage D, the ray-tracing RT on AMReX (branch `stage-d-amrex`, PR #24, `docs/stage_d.md`, merged into `main`); the 256³ production run (branch `production-256`, PR #23, below) and slab-only movie frames `movie_slab`, the movie-reader speed fix and `POST=sbatch` (branch `movie-slab`, merged into `main`). Before that, 2026-10-10: `redshift_per_distance` (branch `redshift-per-distance`, merged); IGM diagnostics (branch `igm-diagnostics`, merged); recombination radiation with the moment RT and helium lines (branch `recomb-moments`, merged into `main`); status assessment and the Hyak production and Stage D prompts ("Threads to launch"); B16 fix `ray_rotate` merged into `main` (branch `claude/exciting-clarke-zn62ih`; production uses `ray_rotate = 3`). The same day: the framework-spike Hyak CPU and GPU benchmarks (PR #19, branch `spike-hyak-bench`: AMReX approved, no ACCESS proposal for now), case A recombination and recombination radiation (branch `recomb-radiation`), and moment RT (M1, OTVET) with its 256³ klone comparison (branch `claude/elegant-bell-rqd3l1`). Before that, 2026-10-09: the survey (PR #16), `hyak-setup` (PR #17), `movies` (PR #12) and `production-run` (PR #7).*

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
8. Smaller: the stored-RT-reference platform policy; the 2 % growth
   deficit; ~~the unsplit open items (Sedov radius convergence, cold-gas
   time-step sensitivity at 128³)~~ done (`docs/hydro_unsplit.md` §6,
   both converge); the AMReX multi-node hydro anomaly
   (`docs/framework_decision.md` §3.3).

**Decisions needed from the user:** the stored-reference policy on klone;
the OTVET formulation question below. (`hydro_scheme = unsplit` and
`redshift_per_distance = 1` for production: decided 2026-10-11,
`docs/PLAN.md`.)

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
   - ray data-structure rewrite;
   - MPI;
   - AMR for gas, then RT;
   - multigrid/tree gravity;
   - adaptive RT.

## Open questions for the user

- RT on AMR (`docs/rt_multicell.md` §9): whether to adopt the proposed
  mixture (a coarse substep with 2^l cells per substep on level l, plus
  photo-chemistry sub-cycling and transport sub-cycling only where a front
  outruns the substep), pure time sub-cycling, or a pure coarse substep.
  Also whether `photon_conserve = 1` (B17) should become part of
  `legacy_physics = 0`; it is a separate flag for now.

- Moment RT: M1 and OTVET are for comparison only (user, 2026-10-11;
  `docs/PLAN.md`). Still open: whether our OTVET (the explicit
  time-dependent two-moment system with the OTVET tensor) should follow a
  specific published formulation (Gnedin & Abel 2001, Petkova & Springel
  2009, ...; `docs/rt_moments.md` §5); that needs the paper.

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
