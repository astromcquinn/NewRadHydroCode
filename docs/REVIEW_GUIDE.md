# Reviewing the repository, with the RT methods in mind

*Written 2026-11-09 for the RT-closure paper. It says where things are,
in which order to read them, what to check in each RT method, and how to
prove that a change kept the results. `docs/STATUS.md` is the log of what
was done; `docs/PLAN.md` the decisions; this file is the map.*

## 1. The shape of the repository

Three layers, bottom up:

| Layer | Where | What it is |
|---|---|---|
| kernels | `include/radray/kernels/` | the arithmetic, header-only, written once: `ray_kernels.hpp` (cast, split, merge, emission, redshift), `rt_cell.hpp` (the per-cell physics: ionisation balance, heating, cooling, `fgam`), `moment_cell.hpp` (M1/OTVET fluxes, the GLF update), `hydro_mh.hpp`, `pm_kernels.hpp` |
| serial port | `include/radray/*.hpp`, `src/*.cpp` | the faithful port of the Fortran plus every option behind a flag; `rt_tools.cpp` (the substep loop), `ray_tracing.cpp` (the ray chains), `rt_moments.cpp` (M1, OTVET), `otvet_closure.cpp` (the alternative tensors), `recomb_rad.cpp`, `plane_parallel.cpp` |
| backends | `src/amrex/` | the same kernels on AMReX levels, MPI + OpenMP: `rt_level.cpp` (rays, sorted arrays, refinement), `moment_level.cpp` (M1/OTVET on z-slabs; refined levels in PR #31), `hydro_level.cpp`, `pm_level.cpp` |

The headers carry the design: every `include/radray/*.hpp` starts with a
description of the component, its equations and its flags, and the
`.cpp` files start with a reading guide of their functions. Read a header
before its implementation.

The documents are one report per workstream (`docs/*.md`, 50 files).
The ones that matter for the RT methods, in reading order:

1. `docs/port_notes.md` sections 1, 5, 5a, 5b: the routine map of the
   original, its quirks (B*, H*, F*) and the fixes behind `legacy_rt`
   and `legacy_physics` (P1 to P11, B17).
2. `docs/rt_moments.md`: what the moment solver does, its tests, the
   first comparison with the rays, the caveats (section 5).
3. `docs/otvet_closure.md`: the rays' Eddington tensor as ground truth
   (section 1), where each closure goes wrong (3), the candidates (4),
   the implementation and its flags (5), the 256 cubed and 512 cubed runs (6, 7).
4. `docs/rt_methods_512.md` and `docs/rt_suite_512.md`: the two 512
   cubed comparisons and their conclusions.
5. `docs/absorbers_closures.md` and `docs/taueff_21cmfast.md`: the
   paper's sections on the photon sink, the single absorber, the tau_eff
   distribution and 21cmFAST.
6. `docs/rt_multicell.md`, `docs/recomb_radiation.md`,
   `docs/stage_e_moments.md`: the k-cell paths, the recombination
   radiation, the moment RT across ranks.

`docs/STATUS.md` has a topic index near the top; `docs/PLAN.md` records
every user decision, so when a choice looks odd, search PLAN.md for the
flag's name before changing it.

## 2. Where each RT method lives

| Method | Selected by | Code | Documents | Tests |
|---|---|---|---|---|
| adaptive rays (the original) | `rt_method = rays` (default) | `src/rt_tools.cpp`, `src/ray_tracing.cpp`, `include/radray/rays.hpp`, `kernels/ray_kernels.hpp`, `kernels/rt_cell.hpp` | `port_notes.md`, `validation_spec.md`, `parallel.md`, `rt_multicell.md` | `ref_*`, `ref_legacy0_*`, `ref_fixed_*` (byte comparison with the Fortran outputs), `rays`, `end_to_end*`, `thread_rt_*`, `photon_conserve`, `ifront_rsla*` |
| M1 | `rt_method = m1` | `src/rt_moments.cpp`, `kernels/moment_cell.hpp` | `rt_moments.md` | `moments`, `thread_rt_m1_lognormal_16src`, `restart_smoke_m1_inject0`, `gamma_avg_m1` |
| OTVET, thin tensor | `rt_method = otvet` | `compute_otvet_tensor` in `rt_moments.cpp`, `src/rt_moments_fft.cpp` | `rt_moments.md` | `otvet_fft_threads_*`, `thread_rt_otvet_*`, `tensor_dump_otvet` |
| OTVET closures: split, screen, cone | `otvet_split = 1`, `otvet_closure = screen \| cone` | `src/otvet_closure.cpp`, `include/radray/otvet_closure.hpp` | `otvet_closure.md` sections 4, 5 | `restart_smoke_otvet_split`, `restart_smoke_otvet_cone`, `tensor_dump_otvet_cone`, `tensor_dump_otvet_split`, `thread_rt_otvet_{split,screen,cone}_*` |
| OTVET with an isotropic background | `otvet_background = ga01` | PR #32 (branch `otvet-background`) | `otvet_background.md` on that branch | on that branch |
| case A and recombination radiation | `recomb_case = A`, `recomb_rad = 1` | `src/recomb_rad.cpp` | `recomb_radiation.md` | `thread_rt_recomb_*`, `restart_smoke_recomb_otvet` |
| the rays' Eddington tensor (diagnostic) | `tensor_dump = 1` | `src/rt_tensor_diag.cpp` | `otvet_closure.md` section 1 | `tensor_dump_rays` |
| moments on the AMReX backend | `rt_backend = amrex` with `rt_method = m1 \| otvet` | `src/amrex/moment_level.cpp` | `stage_e_moments.md` | `amrex_moment_*` (label `amrex`) |

The physics is shared: every method fills `rt_grid%ngam`, `nu` and `dx`
and then runs the same `radiative_physics` (`rrk::photo_gas_cell`), which
returns the transmitted fraction `fgam` that the rays and the moments
apply in their step 3. The injection is the same number of photons per
source cell and substep in both (`rrk::source_photons` for the rays,
`mk::inject_cell` for the moments: `fResc * pop2_init * Luv * dt / vol_rt / nH_avg`
per group), and the redshift of a spectrum is one function
(`rrk::redshift_spectrum`). So a difference between methods is a
difference in transport, not in chemistry or sources.

## 3. What to check, method by method

Read with the paper's description in one hand. For each item the place
to look is given; "verified" means the code was read for this guide and
does what is stated.

### 3.1 Rays

- **Emission.** 192 rays per source cell per substep (`nside_init = 4`,
  `rrk::source_spec`), from the cell centre, equal photons per ray, the
  Pop II spectrum `pop2_init` in 5 groups. A cell with
  `exp(-min_g tau_g) < 0.01` deposits its photons directly instead
  (`rrk::source_is_thick`). Verified. With `ray_rotate = 3` (B16, used in
  every production run) the HEALPix frame of each emission is rotated by a
  hash of the seed, substep and cell (`rrk::frame_rotation`); without it
  the fixed 192 directions imprint streaks and a 12-fold ring.
- **Splitting.** `rmax(nside) = nside sqrt(3 / (pi ray_per_cell))` with
  `ray_per_cell = 1`, so a ray splits when it covers one cell face
  (`rrk::split`: four NESTED children, a quarter of the photons each).
  Verified; described at the top of `src/ray_tracing.cpp`.
- **Merging and the pool.** A ray fainter than `ngam_ray(16)` is merged
  with the faint rays of its cell in the same one of 48 cones
  (`nside_merge = 2`), one brighter than `ngam_ray(32)` may split. These
  thresholds come from the photon histogram of all live rays and the
  pool size `Nray_max = nray_max_factor N^3` (`calc_ngam_ray`), so the
  angular resolution of the method is set by the pool: the paper should
  say which `nray_max_factor` the suite used (15, `rt_suite_512.md`
  section 1) and that `legacy_rt = 0` adds the F3 pool guard
  (`pool_guard`), which raises the split threshold when the pool would
  overflow.
- **Deposit along the path.** `rrk::cast_1`: the segments of the ray's
  path (`ray_segment`; 3 segments in the original, B8, 4 with
  `legacy_rt = 0`), segment n receives `max(dl_n / remaining path, 1 - e^{-tau dl_n})`
  of the photons not yet deposited (the share deposit). `rrk::cast_3` then
  removes what the cells absorbed, `dng (1 - fgam)`, and moves the ray.
  Verified. The photoionisation rate in the cell is
  `Gamma = n_gamma sigma c_hat` with `n_gamma = ngam nH_avg` (the comment
  at "Maximum photoionization rate" in `kernels/rt_cell.hpp`), the same
  expression the moments use.
- **Photon conservation in the chemistry.** The original iterates the
  implicit ionisation balance at most 5 times and can ionise a cell with
  more photons than it was given (B17); `photon_conserve = 1` (part of
  `legacy_physics = 0`) iterates until every group is within 1e-3 of its
  photon limit. `photon_consume_dense = 1` (P11, also part of
  `legacy_physics = 0` since 2026-10-28) makes every photoionisation
  consume its photon; the original keeps a fraction
  `min(e^{-(delta - 50)/50}, 1)` of the photons in gas denser than 50
  times the mean. **The 512 cubed suite was built before P11** and ran
  with the original rule (`absorbers_closures.md` section 3.2), and
  `absorbers_closures.md` section 6 measures the difference (10 to 25 per
  cent in the median Gamma for every method, 2 to 3 points between M1
  and the rays): state this in the paper.
- **Redshifting.** Once per hydro step (`ray_redshifting`), each group
  loses the photons that cross its lower edge to the group below, with a
  power-law spectrum fitted from the ratio of neighbouring groups
  (`rrk::redshift_spectrum`). With `redshift_per_distance = 1` the shift
  follows the distance travelled at `c_hat` rather than the time.

### 3.2 M1

- **Equations and scheme** (`include/radray/rt_moments.hpp`): the two
  moment equations for N and F per group, closed by Levermore's
  `chi(f)` (`mk::m1_pressure`), transported by the first-order global
  Lax-Friedrichs scheme with wave speed `c_hat`, `moment_nsub` subcycles
  (>= 3) per RT substep so that the Courant number is <= 1/3. Verified.
  Consequences to state: fronts are 2 to 3 cells thick, small fronts run
  ahead (`rt_moments.md` sections 4, 5); there is no second-order
  variant.
- **Realizability.** After every subcycle N < 0 is zeroed and |F| is
  clamped to N (`mk::glf_cell`). This is a projection, not a limiter:
  photons are conserved, momentum is not where the clamp acts. The
  photon budget (injected, absorbed, lost) is printed to `ray.log` every
  hydro step and closes to round-off (`moment_photons_total`).
- **Injection order.** `moment_inject = 1` (the suite's choice) injects
  after the transport, so `photo_cell` absorbs the source photons in
  the source cell before they move; `= 0` is the rays' order
  (`rt_moments.md` section 5a explains why neutral source cells needed
  this).
- **Absorption.** `moment_step_3` scales N, F and E of each group by the
  cell's `fgam` (clamped to [0, 1]), the same factor the rays apply to
  their photons.
- **Energy.** E = nu N is advected with N and gives the group's effective
  frequency for the cross sections; the redshift is applied cell by cell
  once per hydro step with the rays' table (`mk::redshift_cell`, the
  flux following the photons).

### 3.3 OTVET

- **Formulation.** This is the time-dependent two-moment system of M1
  with the Eddington tensor replaced by the optically thin tensor
  `h = P / tr P`, `P = sum_s L_s d d^T / |d|^4` over minimum-image
  separations (`otvet_kernel`, `compute_otvet_tensor`; direct sum or an
  FFT convolution, `rt_moments_fft.cpp`), solved explicitly with the same
  GLF machinery. It is **not** the formulation of Gnedin & Abel 2001
  (their isotropic background, PR #32, nor their combination of the
  equations) nor the implicit diffusion form of Petkova & Springel 2009
  (`rt_moments.md` section 5; the user's decision to keep it,
  `PLAN.md`). The paper must say so.
- **The source cell's own weight.** `kSelfWeight = 4` (an isotropic
  tensor in the source cell, 1/r^2 at half a cell); off-diagonal kernel
  components vanish at separation N/2 (two equidistant images). Both are
  choices, visible only where sources compete within a few cells.
- **The closures** (`otvet_closure.hpp` has the definitions, the top of
  `otvet_closure.cpp` the reading guide):
  - split: `P_g = N_g I/3 + min(N_g, N_dir,g) (h - I/3)` with
    `N_dir,g = q_g tr P / (4 pi)` (`mk::otvet_split_pressure`);
  - screen: `e^{-r/lambda}` over the 27 nearest images plus an isotropic
    continuum, lambda from the photon budget unless `otvet_screen_lambda`;
  - cone: the thin near field inside `2 r_c` plus HEALPix cones marched
    from receivers on a grid `otvet_cone_coarse` times coarser through
    mip-maps of the emissivity and the capped HI column, a cutoff
    `e^{-s/lambda_x}` from the budget, and a Gaussian blend with the
    split thin tensor (`otvet_cone_blend`). The defaults (nside 4,
    coarse 4, r_c 8, tau cap 1, 8 boxes, budget 1, blend 4 cells) are in
    `OtvetClosureParams`; the suite used them.
  - "The photon budget" means: the Gamma-weighted photons of the moment
    solution, `sum_g sigma_g N_g` per cell, fix the normalisation of the
    direct light (`OtvetClosureIn::nw`, `OtvetCone::budget`).
- **Cost and diagnostics.** `ray.log` prints lambda, lambda_x, N_dir/N
  and the tensor's wall time every hydro step (`otvet_update_tensor`).

### 3.4 Checks that cut across the methods

- `tests/CMakeLists.txt`: `thread_rt_*_lognormal_16src` (1 vs 4 threads
  byte-identical for every method and closure), `restart_smoke_*`
  (bit-identical restarts), `tensor_dump_*`.
- The reference tests `ref_*` are the Fortran agreement of the rays;
  there is no Fortran for the moments, so their guard is the unit test
  `moments` (`tests/test_moments.cpp`: positivity, conservation,
  isotropy, thread independence) and the stored-output tests above.
- A cheap end-to-end consistency run for the paper: `rt_driver` on
  `tests/problems/lognormal_16src` with `rt_method = rays | m1 | otvet`
  (plus the closure keys), comparing the photon budgets in `ray.log`
  (`Photons`, `Moments : in flight ... injected ... absorbed ... lost`)
  and the I-front positions (`tools/analysis/ifront_multi.py`).

## 4. Proving a change keeps the results

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j
```

```bash
ctest --test-dir build -j4 -LE "long|nightly"
```

For a change that must be invisible, compare the outputs of a run before
and after with no tolerance at all:

```bash
python3 tools/compare/compare.py out_before out_after --rtol 0 --atol-frac 0
```

`rt_driver` runs the RT alone on `tests/problems/<name>`; `radhydro`
runs the full loop on `tests/problems_full/<name>`; `--set key=value`
overrides `params.txt`. The `ref_*` tests do this against the stored
Fortran outputs (`reference/*/outputs*`) on every build. A comment-only
edit can be proved so without running anything:

```bash
diff <(git show HEAD:F | gcc -x c++ -fpreprocessed -dD -E -P - | tr -d ' \t\n') <(gcc -x c++ -fpreprocessed -dD -E -P F | tr -d ' \t\n')
```

Nine `platform_ref` tests differ on klone by design (`docs/cluster.md`);
a few analysis self-tests skip (return 77) without SciPy or matplotlib.

## 5. Repository hygiene

- One branch and one PR per workstream; merge with `main` and run the
  fast suite before pushing (`CLAUDE.md`).
- Branches whose PR is merged can be deleted; `git branch -r --merged origin/main`
  lists them, and `git worktree list` / `git worktree prune` the stale
  worktrees of parallel agents (`.claude/worktrees/`).
- `docs/STATUS.md` is appended at the end of every session; its first
  paragraph is the newest-first summary and the index below it points at
  the sections by topic.

## 6. Stale comments and small things found on 2026-11-09 (fixed 2026-11-10)

Found while commenting; none affected a result. Fixed on 2026-11-10 (branch
`claude/trusting-feynman-2w7i8q`), except the last open item:

- Parameter lists brought up to date: `include/radray/rt_backend.hpp`
  (`amrex_moment_stats`, `amrex_moment_fft_threads`, `amrex_split_on_entry`,
  `amrex_ghost_rule`; the check's comment now lists what
  `check_rt_backend_param` refuses), `hydro_backend.hpp` (the Stage F
  `amrex_*` refinement keys, `sf_overdensity`, `sf_efficiency`),
  `pm_backend.hpp` (`amrex_pk`, `pm_gas_levels`, `pm_level_solve`,
  `pm_mlmg_*`, `pm_level_subcycle`), `hydro_unsplit.hpp` (the `dual_*`
  keys). `docs/stage_d.md` §2.1 now says its list was that of Milestone 2.
- `amrex_state.hpp` and `config.hpp` state `check_amrex_state_param`'s
  actual rules (the moment methods on slabs or the hierarchy, `recomb_rad`
  with them only).
- `xray_background.hpp`: the orphan comment now points to
  `rrk::xray_secondary`; `plane_parallel.hpp` names `check_plane_param`.
- `docs/stage_g.md`: `pm_mlmg_bottom_sweeps = 64` was the default *then*
  (8 now). `docs/pm_refine.md`: checkpoints of versions 1–6 (not 1–4) are
  read with an empty FoF cache.
- `src/halos.cpp` `make_halo_units`: the unused `Npart` removed.
- `src/io.cpp`, `apps/radhydro.cpp`: the `photon_conserve` and
  `photon_consume_dense` errors say -1 is accepted.
- `src/abundance_matching.cpp`: `rms_dex` is the weighted rms of the
  residuals without the monotonicity penalty (the minimised cost is
  unchanged, so the fits are too).
- Compiler warnings of the Release build removed: `checkpoint.cpp`
  (the `memset` that zeroes `DriverState`'s padding, which a
  value-initialisation would not guarantee, now goes through `void*`),
  `tests/test_rays.cpp` (`1e50f`, which is +inf in float, is now `1e30f`),
  `tools/ray_compression/ray_codec.cpp` (the unused `ctx_rec` parameter
  of `code_ray` removed).
- **Open:** coefficients written from memory and still marked so in the
  code: `rrk::xray_secondary` (`kernels/rt_cell.hpp`), the Tinker 2008 mass
  function (`mass_function.hpp`), the HS03 constants (`halo_sources.cpp`).
  None of them enters the RT-methods comparison. Checking them needs the
  papers (Shull & van Steenberg 1985 or Ricotti, Gnedin & Shull 2002;
  Tinker et al. 2008; Hernquist & Springel 2003), which the network here
  cannot reach.
