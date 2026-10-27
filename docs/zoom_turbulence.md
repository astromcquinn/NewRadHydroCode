# Zoom-in runs: reionization-driven turbulence (and later magnetic fields) in one region of a cosmological box

*Plan written 2026-10-22 from the user's request and Cain, McQuinn,
Scannapieco, D'Aloisio & Trac 2026 (PRD 113, 023050, "C26"; the user
uploaded the PDF). Nothing here is implemented yet. The decisions the plan
needs from the user are collected in §9; the request itself is recorded in
`docs/PLAN.md`, "Zoom-in turbulence programme".*

## 0. Summary

C26 resolved kpc-scale turbulence in the post-reionization IGM with
uniform 1024³ and 2048³ boxes of 0.25 h⁻¹ Mpc, ionized by plane-parallel
fronts at a chosen z_re. The goal now is the same physics inside a real
reionization simulation: one region of a 40 (or 20) h⁻¹ Mpc box, refined
to ≥ 1024³ cells at ≤ 0.3–0.6 kpc/h, that an I-front from the box's own
sources sweeps through, with the ICs seeded at the region's resolution and
re-usable for repeated, deeper zooms. Later: magnetic fields.

What the code already has (Stages D–F): an AMReX hierarchy for gas, RT
and (through a uniform mesh) gravity, validated to two refinement levels,
running 512³ over 5 klone nodes; a fixed refinement region per level
(`amrex_refine_box`); nested-phase ICs (`--phases-from-ndm`); particles
with a mass component (the star sink); checkpoints and restarts.

What a zoom needs on top, in order of size:

| Work package | What | Size |
|---|---|---|
| WP0 C26 reproduction | the uniform small boxes of C26 in this code: the plane-parallel RT-domain mode (`docs/plane_parallel.md`), the 0.25 h⁻¹ Mpc resolution ladder, E(k) and the turbulence diagnostics; the code-to-code check and the first science (user request 2026-10-22) | 3–4 weeks, in parallel with WP1–WP5 |
| WP1 zoom ICs | hierarchical white noise, nested displacement fields, multi-mass particles, per-level gas files, separate baryon/DM transfer functions, runtime cosmology | 2–3 weeks |
| WP2 static zoom levels | levels tagged by particle mass class (they follow the Lagrangian region), per-level IC read, multi-mass DM in the PM container and the FoF | 1–2 weeks |
| WP3 gravity on the hierarchy (Stage G) | multigrid Poisson per level with particles on their levels; the uniform mesh cannot reach 0.3 kpc/h in a 40 Mpc/h box | 3–4 weeks |
| WP4 RT into the zoom | raise the HEALPix cap (nside 8192 today) so rays can split into 0.3 kpc/h cells 10 Mpc/h from their source; measure; a plane-parallel "beam" boundary only if needed; the coarse-substep rule at 8 levels | 1–3 weeks |
| WP5 physics knobs | preheating floor T_min(z), z_init, Planck cosmology | 1 week |
| WP6 diagnostics | finest-level cube of the region, E(k), eddy turnover, ionization-time field, region selection | 1 week |
| WP7 production and cost | the ladder Z1 → Z3 on klone, scaling data, the proposal | continuous |
| WP8 MHD | kinematic induction with constrained transport first, full MHD second, AMR coupling, seeds | 6–10 weeks, after Z2 |

Rough cost of the target run (1024³ at 0.3 kpc/h, §7): about 10⁵
core-hours and under 1 TB, i.e. days on the group's 1000 klone cores. The
suite the proposal should ask for (several regions, preheating variants,
MHD, one 2048³ convergence run) is a few 10⁶ core-hours.

## 1. What C26 established, and what the zoom adds

C26 (§II–III, Appendix A): a modified RadHydro, plane-parallel rays cast
from the faces of cubic RT domains, 0.25 h⁻¹ Mpc, 1024³ gas, RT and DM
(Δx = 244 pc/h, 2048³ for convergence), from z = 300 with separate CAMB
transfer functions, hydro + gravity to z_re ≈ 7, then I-fronts with
Γ_HI = 3 × 10⁻¹³ s⁻¹ and J_ν ∝ ν⁻¹·⁵ in five bins, reduced c after the
fronts cross, to z = 4. Findings that set the requirements:

* Turbulence appears only for Δx ≲ 0.5 kpc/h (their Fig. 3 bottom, Appendix
  A: Re_num ≈ 2 Δx_turb/Δx_cell must reach ≈ 50 on the driving scale
  Δx_turb ≈ 10 kpc/h × (Δt/100 Myr)). Eq. (A4): the maximum cell size for
  any turbulence before pressure smoothing finishes is ≈ 1.2 kpc/h.
* The driving is the expansion of filaments and minihalos after impulsive
  heating; E(k) follows k⁻⁵ᐟ³ from ≈ 200 Myr after z_re, between the
  sound-crossing scale and the numerical dissipation scale.
* Preheating suppresses it: T_min = 100 K leaves patchy turbulence around
  the massive filaments, T_min = 1000 K removes it (eq. A6).
* The box size is converged at 0.25 vs 0.5 h⁻¹ Mpc; the fronts are
  plane-parallel and instantaneous for the whole box at z_re.

What a zoom in a full reionization run adds, and the proposal's science
case: the front's real arrival time, direction, speed and spectrum
(hardening, the ĉ-limited phase), the large-scale environment of the
region (its density, the bulk flow, the nearby filaments; C26's tiny
boxes are void-like), the self-consistent relaxation during an extended
reionization rather than a flash, and a volume-filling answer by zooming
on several regions chosen from the base run's ionization history. For
magnetic fields it adds the one thing a uniform box cannot: the Biermann
battery at the front itself, with the front's true geometry.

## 2. What exists (inventory, with the gaps that matter)

| Need | Exists | Gap |
|---|---|---|
| Refined gas levels | F1/F2: `amrex_max_level`, tagging on overdensity, sub-cycled hydro with refluxing, `cool_suppress_overdensity`, the sink | levels are density-tagged; a fixed region only through `amrex_refine_box` (one box per level, tests) |
| RT across levels | F3: the coarse substep (`rt_cells_per_substep` = 2^l), split on entry, ghost rule `conserve`, validated with one and two levels (8 sources) | never run beyond two levels; HEALPix `lside_max = 13` (nside 8192, `constants.hpp`): a ray cannot resolve cells smaller than r/8192 (1.2 kpc/h at r = 10 Mpc/h); two levels are rank-dependent at 1e-10 |
| Gravity with refinement | F4: `pm_gas_levels = 1`, every level kicked from the uniform `pm_refine` mesh (2048³ at 512³ gas) | force resolution is the mesh's everywhere: 20 kpc/h at best; Stage G (multigrid on the hierarchy) not started |
| Particles | `pm_level.cpp` container with a type and a mass component (stars), int64 fixed-point deposit, load balancing with the boxes | dark matter is single-mass (`mass` from `N_dm_m_c`), one `Ndm` per run, the FoF counts particles |
| ICs | `tools/ics/make_ics.py`: Zel'dovich at z = 100, one EH no-wiggle P(k), `--phases-from-ndm` (the coarse modes of a smaller run reused) | one uniform grid; one transfer function; the hierarchy only up to the whole-box Ndm; the cosmology is compile-time (`constants.hpp`, the original's Ω_m = 0.27, σ₈ = 0.80) |
| Plane-parallel RT, T_min floor, DC modes | specified in `docs/plane_parallel.md` | not implemented |
| Outputs | `fields_levels = 1` (a plane at the finest resolution), snapshots of the level-0 composite, movies | no cube of a region at the finest level |
| Production | 512³ + 2 levels on 5 nodes (`amr512_runs.sh`), checkpoints every 10 steps (~240 GB), `fesc512_b2.5_uniform` done (13.4 h), `fesc512_b2.5_amr` resumed after the crash fix | the AMR cost probe and the refined run's timing are not in the docs yet |

## 3. The zoom ladder

Base: the production configuration of the f_esc thread, 40 h⁻¹ Mpc, 512³
gas (Δx₀ = 78.1 kpc/h), 1024³ particles, the ĉ schedule, f_esc(z), the
sink. (A 20 h⁻¹ Mpc base saves one level; the 40 h⁻¹ Mpc box keeps the
calibrated sources and a more representative front geometry. User's
call, §9.) Refinement ratio 2; cell sizes and the region that holds 1024³
cells at each level:

| level l | Δx_l [kpc/h] | 1024³ region [Mpc/h] | particle mass, 1 per cell [M☉/h] | role |
|---|---|---|---|---|
| 0 | 78.1 | — | 4.5 × 10⁶ (the 1024³ base particles) | the box |
| 2 | 19.5 | 20 | 5.6 × 10⁵ | buffer |
| 4 | 4.88 | 5.0 | 8.7 × 10³ | **Z1**: 512³ over 2.5 Mpc/h, the machinery test |
| 6 | 1.22 | 1.25 | 136 | **Z2**: 1024³ over 1.25 Mpc/h; C26's "minimum for any turbulence" (eq. A4) |
| 7 | 0.610 | 0.625 | 17 | (optional rung) |
| 8 | 0.305 | 0.3125 | 2.1 | **Z3**: 1024³ over 0.31 Mpc/h; C26's fiducial class (244 pc/h, 1.2 M☉/h) |
| 9 | 0.153 | 0.156 | 0.27 | **Z4**: 2048³ over 0.31 Mpc/h, convergence (C26's 122 pc/h), proposal scale |

"Repeatedly zoom in" is this ladder: every rung is nested in the previous
one and shares every coarser mode of the ICs (WP1), so Z1, Z2 and Z3 are
the same region at increasing resolution, and a deeper rung (Z4, or a
sub-region of Z3 at level 10) is added without regenerating anything above
it. The intermediate levels are thin shells (a pad of 4–8 of their own
cells around the next finer region, §4.2), so the cell count is dominated
by the finest level: Z3 holds about 1.1 × 10⁹ refined cells plus the 1.3 ×
10⁸ base cells, and about 1.2 × 10⁹ particles.

Region selection (WP6): from the base run, by the ionization-time field
(the redshift at which each base cell passed x_HII = 0.5, a new cheap
output) and the density: a mean-density IGM patch far from the massive
halos, ionized at a chosen time by a front from a known direction; later
a void and a filament-crossing patch for the volume-filling question.
**User (2026-10-22): the regions are chosen on both axes, density and
ionization time** (early and late in the base run's history), since the
relaxation and the front's speed and spectrum differ between them.
The region is Lagrangian: its particles at z_re are traced back to the
ICs (the checkpoints hold the particles with their ids), and the zoom's
particle set is their Lagrangian bounding box plus a pad for the
differential displacement (a few 100 kpc/h at z = 4–6; the bulk flow is
followed by the levels, §4.2).

## 3b. WP0: the uniform small boxes of C26 first (user request 2026-10-22)

Before (and alongside) the zooms, reproduce C26's set-up in this code:
0.25 h⁻¹ Mpc uniform boxes, plane-parallel fronts from the faces of RT
domains at z_re, the resolution ladder of their Fig. 3 (bottom), and the
turbulence they found. It needs none of the zoom machinery, only:

* **The plane-parallel RT-domain mode** of `docs/plane_parallel.md`
  (specified from DM20 and C26, not built): rays from chosen faces of
  cubic RT domains, one per face cell per substep, no splitting or
  merging, deleted at the domain's far face; J_ν ∝ ν⁻¹·⁵ in five
  equal-photon bins; the amplitude set by Γ_HI in the source cells
  (C26: 3 × 10⁻¹³ s⁻¹); turn-on at z_re; ĉ reduced after the fronts have
  crossed. On the uniform grid, serial and backend. The open questions
  of that document (cooling in dense gas, one or both faces per domain,
  the Supplemental §C 2 procedure) go to the user; C26's text says the
  fronts come from the domain boundaries and the gas before z_re evolves
  with hydro and gravity only.
* **ICs:** MUSIC in uniform mode (no zoom): 0.25 h⁻¹ Mpc, 64³ to 1024³,
  z_init = 300, separate CAMB baryon and dark-matter transfer functions,
  Planck cosmology, 1LPT: the WP1 converter's uniform case, and the WP5
  runtime cosmology.
* **Physics:** cooling on everywhere (no sink, no 200× suppression: the
  `cool_suppress_overdensity = 0` flag exists), Compton heating (in the
  original), the T_min floor of WP5 for the preheated variants later.
* **Diagnostics:** WP6's E(k), N_τ, PDFs and slices, applied here first.

Rungs, each against C26's published curves (their Figs. 1–3):

| Rung | Runs | Where | Checks |
|---|---|---|---|
| C1 | 0.25 h⁻¹ Mpc at 64³, 128³, 256³, 512³ (Δx = 3.9 → 0.49 kpc/h), z_re = 7, no preheating | one klone node each (hours to a day) | the onset of turbulence with resolution (their Fig. 3 bottom: no k⁻⁵ᐟ³ range above 0.5 kpc/h, one at 0.47), the temperature maps, the sound-crossing bump and its drift |
| C2 | 1024³ (Δx = 0.24 kpc/h), z_re = 7; then T_min = 100 K and 1000 K | the backends on 2–5 nodes (≈ the cost of Z3's finest level, §7) | their fiducial run: E(k) at z = 7.2 … 4, the k⁻⁵ᐟ³ range, N_τ ≈ 1 at z = 4, the suppression by preheating (their Fig. 2) |
| C3 (optional) | 2048³, flash-ionized as in their convergence run | proposal scale | convergence of the dissipation scale |

What it buys: a direct code-to-code check against the modified RadHydro
of C26 (this code is a port of the same original), the diagnostics and
the coupling-interval test (T4, §4.4) in exactly the setting where the
answer is known, a reference for the zooms (a Z3 region in a C26-like
environment should reproduce C2's statistics), and the first science
(preheating and z_re variants) at a cost klone carries alone. C2 is also
the first item of the ACCESS supplement's plan (§8), with Z2 and Z3
following.

Thread prompt: "Read CLAUDE.md, docs/PLAN.md, docs/STATUS.md,
docs/zoom_turbulence.md §3b and docs/plane_parallel.md. Implement the
plane-parallel RT-domain mode behind parameters (defaults unchanged, the
bitwise suite untouched) with the tests listed in docs/plane_parallel.md,
then the C1 ladder on klone and the E(k) comparison with C26 Fig. 3."

## 4. Design

### 4.1 WP1: zoom initial conditions

**Route (user decision 2026-10-22): MUSIC (Hahn & Abel 2011)
generates, a converter in `tools/ics` writes the code's files.** MUSIC already does the hard parts: hierarchical white
noise with a seed per level (the ladder property for free), the nested
displacement solves with padding, CAMB transfer functions with separate
baryon and dark-matter fields, 1LPT (the user's choice) or 2LPT, and
years of zoom validation; rewriting that in numpy would cost two weeks
and reintroduce the failure modes it solved. The user installs MUSIC on
klone (cloud sessions cannot fetch it). The converter maps MUSIC's generic
per-level output (grids of density and velocity per level, particles per
level) onto the code's conventions below; `make_ics.py` keeps a small
built-in zoom (the steps below, 1–2 levels) for the fast-suite problems,
which cannot depend on MUSIC. To check when MUSIC is set up: how it sets
the gas velocity with separate transfer functions, and that its growth
rates match the code's runtime Planck cosmology (item 5).

The steps, whichever tool produces them:

1. **Hierarchical white noise.** Level 0: the base run's draws (the
   `--seed 7 --phases-from-ndm 512` stream of the production ICs, so the
   zoom's box is the base run's box exactly). Level l + 1 over its region:
   draw unit white noise on the 2× grid, subtract each block's mean of 8
   and add the parent's value scaled by 1/√8, so that every block average
   reproduces the parent's noise and the small-scale modes are new. The
   draws of each level are keyed by (seed, level, region origin), so a
   deeper rung later reuses every existing level's draws unchanged.
   (Hahn & Abel's refinement rule; check the coefficient against the paper
   before relying on it.)
2. **Nested displacement and density fields.** Per level: the long modes
   come from the parent's fields interpolated onto the level's grid (band
   limited at the parent's Nyquist), the short modes from an FFT of the
   level's noise on its padded region with the parent's band removed.
   Validate: a zoom that covers the whole box at level 1 must reproduce
   the uniform 2N run's ICs (as `--phases-from-ndm` is validated today, to
   float round-off); the density power spectrum inside the region matches
   the input P(k) to the level's Nyquist; the displacement is continuous
   across the region's boundary; a 2-level nest matches the equivalent
   1-level nest at 4× (the ladder property).
3. **Particles per level:** one particle per cell of the finest level
   covering each Lagrangian cell (mass 8⁻ˡ of a base-cell particle, or
   8⁻ˡ⁻¹ with 2× particles as the base run has), written as one `xv_dm`
   stream per mass class with the class's mass in a header file; ids
   encode the class.
4. **Gas per level:** `u_gas.init.L<l>` on the level's region grid
   (density, momentum, kinetic energy as now), read at the first regrid
   (WP2).
5. **Transfer functions and cosmology.** Separate baryon and dark-matter
   transfer functions from a CAMB table supplied by the user (the cloud
   cannot fetch them; `docs/plane_parallel.md`), 1LPT (user decision),
   z_init a parameter (C26 used 300; the gas Jeans-scale structure of a
   cold IGM starts earlier than our z = 100). Planck parameters (user
   decision) need the cosmology
   to become a runtime parameter set in the code (`constants.hpp` holds
   the original's float literals as `constexpr`; the defaults stay the
   original for the bitwise tests, `omegam`, `omegab`, `hubble0`, `sigma8`,
   `ns` as `params.txt` keys): a small but wide change (cosmology.hpp, the
   rates' units, the halo units, the IC tool).

Sizes: a level-8 region of 1024³ with a 128-cell pad is 1280³ complex
arrays, about 35 GB per field in double; the generator runs on one klone
node in an hour, like the 1024³ ICs today (293 GB, 10 min).

**Status (2026-10-23, thread D, branch `zoom-ics`):** built as
`tools/ics/music_convert.py` (MUSIC's generic HDF5 → `xv_dm.init`,
`u_gas.init`, `u_gas.init.L<l>`, `xv_dm.L<l>`, `ics_levels.txt`) and
`tools/ics/music_ladder.py` (the ladder's configurations), with T0, the
ladder check and a fast-suite zoom (`tools/ics/music/README.md`; decisions
in `docs/PLAN.md`). Two corrections to the text above: MUSIC2 constrains a
level's noise to its parent by Fourier splicing over the level's patch,
not by the block-mean rule of step 1, so the ladder property needs every
rung to put each shared level on the same grid (the generator does); and
the fast suite uses a committed MUSIC zoom (895 KB) instead of a built-in
nest in `make_ics.py`. The amplitude and velocities follow the code's
growth (step 5): with CAMB's velocity columns the CDM moves at 0.92× the
code's growing mode at z = 300 and the gas lags (0.74× at k = 50 h/Mpc,
0.06× at 2500), CAMB's physics rather than the code's single-fluid growth.

### 4.2 WP2: static zoom levels

* **Tagging by particle mass class** instead of fixed boxes: a level-l
  cell is tagged when it holds a particle of class > l (the Enzo and
  RAMSES "must-refine particles" rule). The levels then follow the
  Lagrangian region through its bulk flow, the pad of each level is the
  particles' own spread, and a coarse particle inside the finest region
  is the standard contamination diagnostic (count and report it every
  step). `amrex_refine_box` stays for the tests; density tagging stays
  available in addition (off by default inside the zoom: the finest level
  must be uniform for the turbulence measurement, as in C26).
* **Proper nesting:** `amrex_n_error_buf` and AMReX's `n_proper` give the
  shells; the IC pad (§4.1) must exceed the buffer of every level.
* **Per-level ICs:** at the first regrid the new levels are filled from
  `u_gas.init.L<l>` where the file covers them and from the parent
  elsewhere (today: always from the parent).
* **Multi-mass dark matter:** the container's mass component (today 1
  for DM, the stars' mass otherwise) carries the class mass; the deposit
  already weights by mass; `summarize_dm`, the host mirror, the
  checkpoint (version bump) and the FoF (mass-weighted, with the
  minimum halo mass a mass, not a count) follow. The serial PM
  (`src/pm.cpp`) is left single-mass: zooms run on the backends only.
* **Sources:** the halo finder sees halos of a few 10⁴ M☉ in the zoom;
  the source model keeps its global mass floor (`tvir_min`, the
  abundance-matching range), so the zoom adds no sources below the base
  run's floor. Sources inside the region that the base run already has
  must be the same objects (check the catalogue across the ladder).

### 4.3 WP3: gravity on the hierarchy (Stage G)

The uniform `pm_refine` mesh stops at 2048³ (20 kpc/h in 40 Mpc/h); the
zoom needs the force at 0.3 kpc/h. The route in `docs/ROADMAP.md` §6:

* Level 0 keeps the FFT (periodic, bitwise today).
* Level l ≥ 1: particles live on the finest level covering them (AMReX's
  `AmrParticleContainer` redistribution by level), CIC deposit on the
  level, plus the gas of the level; solve lap φ_l = (2a/3) δ_l with the
  Dirichlet boundary from φ_{l−1} interpolated (the level solve of Enzo
  and RAMSES) using AMReX's `MLMG` with `MLPoisson` (cell-centred); the
  force by the 4-point gradient as the PM does; the kick per level on
  its own cells (the F4 structure stays). Nyx's composite multigrid
  solve is the alternative; the level solve is simpler to validate
  bitwise across ranks and matches the sub-cycled hydro.
* Determinism: the multigrid residual norm is a reduction; use AMReX's
  fixed-iteration mode or a deterministic reduction so the result is
  rank-independent (the project's rule 5).
* Validation: the Ewald single-particle test of `docs/pm_refine.md` on a
  level; the refined-everywhere run against the uniform 2N run (the
  `pm_gas_levels` equivalence test, now with the level solve); a small
  zoom (32³ base, 2 levels) against the uniform 128³ run: halo masses
  and positions inside the region, gas density profiles; momentum
  conservation.

### 4.4 WP4: the RT into the zoom

* **The HEALPix cap.** `lside_max = 13` (nside 8192, int32 pixel indices
  from the original's ray type). A ray from a source at r must have
  nside ≳ r/Δx to deposit at the cell scale: at r = 10 Mpc/h that is
  3 × 10⁴ for Δx = 0.3 kpc/h (level 8) and 8 × 10³ for 1.2 kpc/h (level
  6). Z2 sits at the cap, Z3 beyond it. First step: `lside_max` 16–17 with
  64-bit pixel indices in the ray structs (both encodings), the NESTED
  routines checked at that depth (`healpix_hd.hpp` ports HEALPix 3.11,
  which handles int64), the merge bins untouched. The split-on-entry
  rule then carries the validated F3 machinery into every level.
* **Count and cost.** The number of rays needed to cover the region at
  the fine cell scale is (L/Δx)² per occupied direction, about 10⁶ ×
  (number of bright directions) ≈ 10⁸–10⁹ in flight at Z3 for as long as
  the region is being ionized, at 140 B (or less with `rstruct`) and 2 µs
  per ray-substep: 10²–10³ core-s per substep, tens of GB. Measured at
  Z1 and Z2 before Z3. After overlap the merge criterion thins them.
* **If that is too many: beams at the zoom boundary.** Rays entering a
  chosen level are binned per face cell and direction and re-emitted as
  non-splitting parallel rays, which is exactly C26's RT-domain scheme
  with the incident flux supplied by the coarse rays instead of a fixed
  Γ_HI (photon conserving by construction; the geometric 1/r² dilution
  across the region, 5–20 %, is then approximated per beam from its
  origin distance). This is a new flag and a new validation against the
  splitting rays on Z1; build it only if the measurement asks for it.
* **The coarse-substep rule at 8 levels.** Rays cross 2⁸ = 256 fine
  cells per base substep. The Stage F decision (one coarse substep, no
  transport sub-cycling) was taken for two levels and dense gas; in the
  zoom the front's passage is the measurement. C26's fronts crossed
  their domains at full c and the paper's Appendix A ties the turbulence
  to the heating impulse, not to the front's internal structure, so a
  front smeared over a fraction of a base cell per substep may be
  acceptable, but it must be shown: run the uniform plane-parallel
  problem at `rt_cells_per_substep` = 16, 64, 256 (the infrastructure of
  `docs/rt_multicell.md`) and compare the front's temperature and
  ionization profiles and the post-front E(k) at the level of C26's
  Fig. 3. If it fails, the options are a larger ĉ while the region is
  crossed (the schedule by ionization state exists; a zoom-aware
  trigger is a few lines) or transport sub-cycling inside the zoom
  levels only, which reverses the Stage F decision there (§9).
* **Spectrum:** the base run's 5 Pop II groups carry over.
* **Speed of light (user, 2026-10-22): a slow ĉ such as 0.1c is
  acceptable, under the condition that light crosses the smallest
  subregion more slowly than sound crosses one of its cells,
  L/ĉ > Δx/c_s, so that the gas can respond while the front sweeps the
  region.** With L/Δx = 1024 on every 1024³ rung and c_s = 20 km/s in
  the heated gas, that is ĉ ≲ 1024 c_s ≈ 0.07c: 0.1c misses it by 1.5×,
  0.05c meets it with margin, and a 2048³ rung (Z4) allows 0.14c. The
  ĉ schedule by ionization state can carry a zoom-aware entry (a lower
  ĉ while the region is being crossed), and a lower ĉ is cheaper (the
  substep count scales with ĉ).
* **The RT–hydro coupling interval is the real constraint, and today it
  is the hydro step.** The RT runs all its substeps inside one hydro
  step and the hydro of every level follows; the base step is set by
  `ra_max` (da/a ≤ 0.1) and the base-cell signal, 50–100 Myr at z ≈ 6–8
  (the 512³ runs take ~55 steps from z = 100 to 5.5), in which light at
  0.1c travels 16–24 Mpc/h and a front 1–3 Mpc/h: the whole zoom region
  is ionized and heated within one step and relaxes only in the next.
  C26's uniform boxes couple every cell-crossing of light. So the
  condition above only buys something together with a short coupling
  interval while the front is in the zoom. Two cheap parameters, both
  default-off: a cap on the hydro step in Myr (`dt_max_myr`, next to
  `ra_max`), and a cap in RT substeps per hydro step, down to one
  (hydro step = the light crossing of a base cell at ĉ: 0.46 Myr at
  0.1c and z = 7, about 17 finest-level CFL steps at Z3, 200 hydro
  steps during a 100-Myr passage). Coupling at the finest CFL step
  itself would need RT substeps shorter than a base-cell crossing (rays
  crossing ~15 fine cells), i.e. transport sub-cycling in the zoom, the
  Stage F question again. The T4 test decides how short the coupling
  must be: the uniform plane-parallel problem with the heating applied
  in one lump per 100, 10, 1 and 0.1 Myr against the per-substep
  coupling, compared on the post-front E(k) and the filament expansion.
  The base-level overheads per hydro step (FoF by time already, the
  summaries, `post_rt`) multiply with the step count and are to be
  measured at Z1.

### 4.5 WP5: physics knobs

* **Preheating floor** T_min from z_floor (C26: 100 K and 1000 K for
  z < 15; `docs/plane_parallel.md`), applied in the RT cell physics; the
  default 0 keeps every reference. **User decision 2026-10-22: the first
  runs are unheated; preheating is added later.**
* **Runtime cosmology and z_init** (§4.1 item 5).
* Cooling stays on in the dense gas (the Stage F decision) with the sink
  above a few thousand × mean; C26 had no sink and no star formation in
  its tiny boxes, so the sink threshold inside the zoom is a choice to
  record (the minihalos that drive the turbulence are far below it).

### 4.6 WP6: diagnostics and region selection

* **Ionization-time field** z_ion per base cell (first x_HII > 0.5),
  written with the snapshots and used to pick regions.
* **Finest-level cube** of the zoom region (`fields_cube`: density,
  temperature, velocity, x_HI, Γ, later B) in float at the snapshots,
  4 GB per field at 1024³.
* **Analysis** (`tools/analysis/zoom_turbulence.py`): E(k) per component
  (C26 eq. B2), the turnover time and N_τ (B1), the k⁻⁵ᐟ³ range, the
  density and temperature PDFs, the sound-crossing bump, the resolution
  ladder overlay (C26 Fig. 3), the Reynolds number estimate (B3), and
  for MHD the magnetic energy spectrum and the growth rate.
* **Contamination and coverage:** coarse particles in the finest region,
  the fraction of the Lagrangian region still inside the finest level.

### 4.7 WP7: production, cost and the proposal

See §7 and §8.

### 4.8 WP8: magnetic fields (after Z2 validates)

Staged so that the first science result (the kinematic growth rate) comes
from the smallest change:

1. **Kinematic induction.** Face-centred B on each level (three face
   MultiFabs), constrained transport with the upwind edge EMFs of
   Gardiner & Stone 2005/2008 computed from the hydro's face states and
   HLLC velocities; no back-reaction. ∇·B = 0 to round-off on each level.
   This measures the dynamo's kinematic growth from a seed, which is the
   quantity C26 §IV argues about, with no change to the gas dynamics
   (and bitwise identical gas when the field is off).
2. **Full MHD.** HLLD (Miyoshi & Kusano 2005) in `hydro_mh.hpp`'s
   predictor–corrector (the MUSCL-Hancock + CT scheme of Fromang,
   Hennebelle & Teyssier 2006 is the closest published form), the
   magnetic pressure and tension in the update, the dual-energy rule
   extended with the magnetic energy (a seed of 10⁻²⁰ G is 10⁻²² of the
   thermal energy: B must stay a separate variable, in double, and never
   be recovered from the total energy).
3. **AMR coupling.** Divergence-free prolongation at coarse-fine faces
   (Balsara 2001; Tóth & Roe 2002), restriction by face averaging, and
   EMF refluxing along the coarse-fine edges (AMReX's `FluxRegister` is
   for cell fluxes; the edge register is written by hand as Castro does).
   Validate with the Orszag–Tang vortex and a circularly polarized Alfvén
   wave on the hydro suite, refined and unrefined.
4. **Seeds and the cosmological form.** A uniform comoving seed B₀
   (∝ a⁻²) set at the ICs; the Biermann battery ∝ ∇p_e × ∇n_e as a source
   term evaluated from the RT's electron density and temperature in the
   front (new physics only this code can do self-consistently); the
   induction equation in the code's comoving variables to be derived and
   documented (the original has no MHD; nothing to keep bitwise).
   **User (2026-10-22): the Biermann battery and one more non-ideal MHD
   term at the I-fronts are wanted; which term, and the details, are
   deferred until the hydrodynamic zooms exist.**
5. **Not planned:** explicit resistivity beyond that; the effective
   magnetic Reynolds number is set by the grid and reported from the
   resolution ladder, as C26 does for the kinetic cascade.

Cost: HLLD + CT roughly 1.5–2× the hydro; three extra face fields and
the EMFs about +50 % memory on the fine levels.

## 5. Validation ladder

Every work package gets tests in the fast suite at tiny sizes and one
klone run at the next size up, in the repo's existing style
(`docs/stage_f.md` §3):

| Rung | Problem | Reference |
|---|---|---|
| T0 | zoom ICs over the whole box at level 1 (WP1) | the uniform 2N ICs, float round-off |
| T1 | 32³ base, 2 zoom levels over the central eighth, gravity only (WP2, WP3) | uniform 128³: halo masses, positions, density PDF in the region |
| T2 | the same with hydro and the sink | uniform 128³ gas fields in the region |
| T3 | 8 sources (`uniform_8src`) with a 3-level zoom; rays split to level 3 (WP4) | uniform 256³: x_HII,V, photons in flight, front radii (the F3 metrics) |
| C1, C2 | C26's uniform 0.25 h⁻¹ Mpc boxes with the plane-parallel mode (§3b) | C26's published E(k), maps and preheating dependence |
| T4 | plane-parallel front on the uniform grid at `rt_cells_per_substep` 16–256 (in the C1 setting) | k = 1: front profiles, post-front E(k) |
| T5 | 128³ base, 10 Mpc/h, Z1-like zoom (4 levels) on 2 nodes | the uniform 512³ run of the same box (hours on 5 nodes): everything inside the region |
| Z1 | 512³ base, level 4 over 2.5 Mpc/h | the base run's level-0 composite outside the region; cost and ray counts |
| Z2 | level 6, 1024³ over 1.25 Mpc/h | C26's convergence sequence (their Fig. 3 bottom); first science |
| Z3 | level 8, 1024³ over 0.31 Mpc/h | C26's fiducial; the proposal's benchmark |

Determinism: every rung byte-identical across rank counts, through a
restart, and with `amrex_max_level = 0` bitwise the Stage E level; the
two-level 1e-10 rank dependence noted in `docs/STATUS.md` is to be found
before the ladder goes beyond two levels.

## 6. Order and threads

```
WP0 C26 boxes: plane-parallel mode ── C1 ── C2 ── (T4 in this setting)   [parallel track, from now]
WP1 ICs ──────────┐
WP3 gravity ──────┼── T1/T2 ── WP2 levels ── T3 (WP4 cap) ── T5 ── Z1 ── Z2 ── Z3
WP5 knobs ────────┘                 │
WP6 diagnostics ────────────────────┘ (needed from Z1 on)
WP4 substep test T4 ── any time (uniform grid)
WP8 MHD ── after Z2 (kinematic first), on the hydro suite, then Z2-MHD
Proposal ── after Z1's cost and scaling numbers (§8), in parallel with Z2
```

WP0, WP3 and WP5 are independent and start now, as parallel threads
(worktrees, one branch and PR each); WP1 starts once MUSIC and the CAMB
tables are on klone; T4 belongs to the WP0 thread, since it needs the
plane-parallel front. The memory work (below, thread E) follows WP3
because both touch the PM level. Opening prompts, ready to paste
(revised 2026-10-22 after the decisions of §9):

**A. WP0, the plane-parallel RT-domain mode and the C1 ladder (branch
`plane-parallel`).** "Read CLAUDE.md, docs/PLAN.md, docs/STATUS.md,
docs/zoom_turbulence.md (§3b and §4.4) and docs/plane_parallel.md.
Implement the plane-parallel RT-domain mode behind parameters, defaults
unchanged and the bitwise suite untouched: rays launched from chosen
faces of cubic RT domains, one per face cell per substep, no splitting or
merging, deleted at the domain's far face; J_ν ∝ ν^-1.5 in five
equal-photon bins; the amplitude set by Γ_HI in the source cells; turn-on
at z_re; on the uniform grid, serial RT and the AMReX backend. Before
building the domain geometry, ask the user the open questions of
docs/plane_parallel.md (one or both faces per domain, cooling in dense
gas, the Supplemental §C 2 procedure). Add the tests listed there
(uniform-gas front against the analytic position and speed, slab
attenuation, photon conservation, the Γ normalisation, domains). Then
the T4 test of §4.4: two default-off caps, the hydro step in Myr and the
RT substeps per hydro step, and the uniform front with the heating
coupled per 100, 10, 1 and 0.1 Myr against per-substep coupling, compared
on the front profiles and the post-front E(k). Finally the C1 ladder on
klone (0.25 h⁻¹ Mpc at 64³ to 512³, z_re = 7, ICs from the existing
generator until WP1 lands) against C26 Fig. 3. Record decisions in
PLAN.md and results in STATUS.md."

**B. WP3, gravity on the hierarchy (branch `stage-g-multigrid`).** "Read
CLAUDE.md, docs/PLAN.md, docs/STATUS.md, docs/zoom_turbulence.md §4.3,
docs/stage_f.md §4 and docs/stage_e_pm.md. Stage G: particles on the
finest level covering them, the CIC deposit per level, and the level
Poisson solve lap φ_l = (2a/3) δ_l with AMReX MLMG (MLPoisson) and
Dirichlet boundaries from the parent's potential, the 4-point gradient
and the per-level kick; level 0 keeps the FFT, bitwise. Behind a flag
(`pm_level_solve`, default off), with `pm_gas_levels` kept as the
comparison. Deterministic across rank counts (fixed iteration counts or
a deterministic reduction). Tests: the Ewald single-particle force on a
refined level, the refined-everywhere run against the uniform 2N run
(the `pm_gas_levels` equivalence test with the level solve), momentum
conservation, rank independence and restart; then T1 of §5 (32³ base,
two levels over the central eighth, gravity only, against uniform 128³:
halo masses and positions, the density PDF in the region). Profile the
solve at 128³ with two levels on klone."

**C. WP5, runtime cosmology and the physics knobs (branch
`runtime-cosmology`).** "Read CLAUDE.md, docs/PLAN.md, docs/STATUS.md,
docs/zoom_turbulence.md §4.5 and §4.1 item 5. Make the cosmology a
runtime parameter set (omegam, omegab, hubble0, sigma8, ns as params.txt
keys) with the original's float literals as the defaults, so every
stored reference and Fortran comparison stays bitwise; thread it through
cosmology.hpp, the rate and halo units and tools/ics (the IC tool reads
the same keys). Add z_init as a parameter and a temperature floor
T_min from a redshift z_floor (C26: 100 K and 1000 K for z < 15) in the
RT cell physics, default off. Tests: the defaults bitwise (the fast
suite), a Planck run differing only where it should, the floor's
identity at 0 and its effect on a cooling test. Document in
docs/production_run.md."

**D. WP1, zoom ICs with MUSIC (branch `zoom-ics`; after MUSIC and the
CAMB tables are on klone).** "Read CLAUDE.md, docs/PLAN.md,
docs/STATUS.md, docs/zoom_turbulence.md §4.1 and §4.2, and
docs/pm_refine.md §2. Write the converter from MUSIC's generic per-level
output to the code's files: one particle stream per mass class with its
mass, `u_gas.init.L<l>` per level on the level's region grid, the base
files as today; with the Planck cosmology and 1LPT; a MUSIC
configuration generator for the ladder (seed per level, the region from
a Lagrangian bounding box plus pad). Tests: T0 (a MUSIC run with no zoom
converted and run against the existing generator's ICs at the same seed
statistics: P(k) and the growth), the ladder property (a level-l region
regenerated with a deeper rung keeps every coarser mode), and a tiny
2-level zoom problem for the fast suite built with the in-repo
generator's 1–2-level nest (§4.1 steps 1–2) so the suite does not depend
on MUSIC."

**E. Memory for the ACCESS nodes (branch `rank-memory`; after B
merges).** "Read CLAUDE.md, docs/STATUS.md ('Known limits'),
docs/stage_e_state.md and docs/stage_f.md §4 (the refined run's 39 GB
per rank and 245 GB on rank 0). Cut the redundant per-rank host state
and rank 0's mirror: the halo finder and the sources distributed (or run
on a gathered subset by mass class), the outputs and checkpoints written
from the levels without a full rank-0 mirror, and the IC read scattered.
Target: under 4 GB per core at 512³ with two levels and at 1024³
particles, rank 0 no larger than the others by more than the catalogue.
Bitwise outputs and checkpoints against the current code; measure on
klone at 512³."

The user's own items alongside: MUSIC and the CAMB transfer-function
tables on klone (for D); deleting the stale branch
`claude/kind-wozniak-4fiw74`; the open questions of thread A when it
asks.

## 7. Cost estimates (to be replaced by measurements at Z1 and Z2)

Per-cell and per-ray costs from `docs/production_run.md` §5 and
`docs/stage_d.md` §4.4; the finest level dominates.

| Item | Z2 (level 6) | Z3 (level 8) | Z4 (level 9, 2048³) |
|---|---|---|---|
| finest-level cells | 1.1 × 10⁹ | 1.1 × 10⁹ | 8.6 × 10⁹ |
| finest-level hydro steps, z_init → 4 (CFL at 1–30 km/s) | ≈ 2.5 × 10⁴ | ≈ 10⁵ | ≈ 2 × 10⁵ |
| hydro, at ≈ 3 µs per cell-update incl. AMR overhead | 2 × 10⁴ core-h | 8 × 10⁴ core-h | 1.4 × 10⁶ core-h |
| RT cell physics, ≈ 8000 base substeps × 1.7 µs per cell | 4 × 10³ | 5 × 10³ | 4 × 10⁴ |
| rays, 10⁸–10⁹ in flight while the region is crossed | ≤ 5 × 10³ | ≤ 5 × 10³ | ≤ 4 × 10⁴ |
| **total** | **≈ 3 × 10⁴ core-h** | **≈ 10⁵ core-h** | **≈ 1.5 × 10⁶ core-h** |
| memory (cells × ≈ 0.5 kB, particles × 64 B, rays) | 0.7 TB | 0.7 TB | 5 TB |
| checkpoint | ≈ 0.4 TB | ≈ 0.4 TB | ≈ 3 TB |
| on the owned 992 klone cores | ≈ 1.5 days | ≈ 5 days | not on klone |

MHD multiplies the hydro line by 1.5–2. The large uncertainties are the
AMR overhead with 8 levels of FillPatch and refluxing (never profiled;
the Sedov two-level run paid 60× for 8× the cells and 2× the substeps,
part of which is the whole-box refinement), the step count (the cold
pre-reionization gas with 10–30 km/s flows), and the ray population. The
ladder exists to pin these before the proposal is written.

**A first measurement of the AMR overhead** (2026-10-22: `fesc512_b2.5_amr`
against its unrefined control `fesc512_b2.5_uniform`, `docs/stage_f.md`
§4; 40 Mpc/h, 512³ base grid + 2 levels above 8 ρ̄ and 64 ρ̄, 1024³ particles
on the 2048³ mesh, 5 klone cpu-g2 nodes, 120 ranks × 8 threads, z = 90 → 5):

| | uniform 512³ | + 2 levels | ratio |
|---|---|---|---|
| refined cells at z = 5 (level 1 / level 2, lower bounds from the density thresholds) | – | ≥ 6.4 M / ≥ 2.0 M (≈ 6 % of the 134 M base cells) | |
| step time z = 90 → 5 (the per-step timers; without checkpoints) | 4.3 h | 11.3 h | 2.6 (3–4 per step after z = 6) |
| of which RT | 2.15 h | 7.96 h | 3.7 |
| of which hydro | 0.03 h | 1.15 h | 46 |
| rays in flight at z = 5 | 1.5 × 10⁹ | 4.4 × 10⁹ | 2.9 |
| hydro cell-substeps (the kernel's counter, all levels; sub-cycling 1482 substeps in 248 calls, up to 56 in one, against 1 per call) | 1.5 × 10¹⁰ | 4.9 × 10¹⁰ | 3.4 |
| hydro core-seconds per counted cell-substep | 6 µs | 78 µs | 13 |
| memory, typical rank (peak RSS, min over ranks) | 17 GB | 39 GB | 2.3 |
| one checkpoint (`restart_every = 1`, 4 × 10⁹ rays) | 290 s | 760 s | 2.6 |
| core-hours, step time / job wall-clock | 4 100 / 12 900 | 10 900 / 36 000 | 2.6 / 2.8 |

(The refined run's wall-clock includes two crashes, a 9.5 h hang before
the abort-on-error fix and the checkpoint of every step.) Two lessons for
the estimates above. (i) The hydro's AMR overhead is per box, not per cell:
with the refined cells scattered over the halos, FillPatch (2150 s),
refluxing (1010 s) and the CFL reduction (660 s) are 95 % of the 4020 s of
hydro, and a counted cell-substep costs 13× the uniform run's, far from
the 3 µs line, which came from the whole-box Sedov test. The zoom's one region
is the favourable case (few large boxes), so the truth lies between; Z1
measures it. (ii) The sources sit in the refined gas, so most rays are
born at the finest level and split there: 6 % more cells cost 2.9× the
rays and 3.7× the RT. The "rays" line above assumes 10⁸–10⁹ in flight;
with the zoom region holding the sources of interest, size the ray pool
from the finest level's cell count, not the base grid's, and expect the RT
to dominate the hydro until the refined cell count is large.

## 8. The allocation proposal

klone builds and validates the machinery (T0–T5, Z1, Z2 and probably one
Z3); the proposal asks for the suite that klone cannot hold at once:
several regions × preheating (none, 100 K, 1000 K) × MHD seeds, and Z4.

Venue: **NSF ACCESS**, through the group's existing Accelerate project
PHY250179 ("High-resolution cosmological simulations of intergalactic
gas", the same NASA grant as this work; it holds the 3 million credits of
the Accelerate ceiling and ends in late 2026). What the ACCESS portal
offers for it (checked by the user, 2026-10-22), and the route (user
decision):

1. **A supplement of 1.5 million credits, now.** The project is eligible
   for it; it is a short justification on a rolling review. It funds the
   near-term plan in full: the C26 reproduction (C2 and its variants,
   §3b), Z2 and the first Z3 (≈ 10⁵ core-hours each, §7), with margin for
   the measurements that the next request needs. Draft:
   `docs/access_supplement/access_supplement_2026.tex` (Overleaf-ready, with figure placeholders).
2. **An extension** of the project's end date, up to six months past the
   grant's end, when the end date approaches, so the credits are not lost
   while the plane-parallel mode and the zoom machinery are built.
3. **A renewal** (available from 30 days before the end date) as the next
   Accelerate, up to the 3 million ceiling, for the suite: several regions
   × preheating (none, 100 K, 1000 K), the MHD runs, written with the Z2
   and Z3 results.
4. **Maximize** (semi-annual committee review, scaling data on the target
   machine) only if the suite plus Z4 outgrows the Accelerate ceiling.

A Discover request is not needed: it would be a second project on the
same grant for less than the supplement. The policy point behind this:
ACCESS allows a PI several projects but ties them to merit-reviewed
supporting grants (one project per grant, or related work combined), so
everything under the NASA grant lives in PHY250179. NASA HEC remains a
fallback. A GPU request only makes sense once the AMReX backends run on
GPUs (Stage D §5: prepared, not built), where the H200 was 5× a
`cpu-g2` node on the RT; the hydro at the finest level is the part that
would gain most.

Skeleton, with the numbers each section needs:

1. **Science justification.** C26's result and its open questions (§1):
   volume filling, the preheating dependence, the front's real geometry,
   the B-field growth with a self-consistent battery. The ladder as the
   convergence argument.
2. **Code and readiness.** The faithful port validated bit for bit; MPI
   through AMReX, 512³ on 5 nodes with the measured scaling
   (`docs/stage_e_state.md`, `docs/stage_d.md` §4.4); AMR to two levels
   validated; Z1/Z2 results and timings; restarts and checkpoint
   cadence; the determinism guarantee.
3. **Computational plan.** The table of runs (Z3 × regions × variants, Z4,
   the MHD set), core-hours per run from the measured Z1/Z2 cost scaled
   as §7, storage (checkpoints, snapshots, the finest-level cubes) and
   the data retention plan.
4. **Scaling.** Strong scaling of a Z2 step on 1–5 klone nodes (the
   finest level's boxes over the ranks, the FillPatch share), the memory
   per node, and weak scaling Z1 → Z2 → Z3; a benchmark on the target
   machine if the venue requires it (its test allocation).
5. **Timeline and personnel**, data management, the publication plan.

Before writing: Z1 and Z2 done with their timings and scaling runs
(`tools/hyak/amrex_scaling.slurm` pattern), the Z3 cost extrapolated,
and the venue's deadline. A realistic date is after WP1–WP6 and Z2,
about three months of threads from now.

## 9. Decisions (user, 2026-10-22)

1. **Base box:** 40 h⁻¹ Mpc / 512³, the calibrated production
   configuration; 8 levels to 0.3 kpc/h.
2. **Rungs:** Z1 → Z2 → Z3 → Z4 as in §3.
3. **Regions:** one mean-density IGM patch first, then a void and a
   filament-crossing patch, and on the second axis regions ionized
   early and late; selection as in §3.
4. **ICs:** MUSIC + a converter (confirmed); 1LPT; Planck cosmology,
   which needs runtime cosmology in the code.
5. **RT:** raise the HEALPix cap first and measure the ray counts;
   beams only if needed. Whether the coarse-substep rule may be relaxed
   inside the zoom waits for the T4 test.
6. **Preheating:** none in the first runs; added later.
7. **MHD:** later; the Biermann battery and one more non-ideal term at
   the I-fronts are wanted (§4.8).
8. **Proposal:** a 1.5 million-credit supplement to the existing
   Accelerate project PHY250179 now, its extension, then its renewal for
   the suite; Maximize only beyond that (§8; revised 2026-10-22 from
   "Discover, then Accelerate").
10. **C26 reproduction first (user, 2026-10-22):** the uniform small
   boxes with the plane-parallel mode as WP0 (§3b), a parallel track
   from now, and the first content of the supplement's plan.
9. **Speed of light:** a slow ĉ (0.1c) is fine provided light crosses
   the smallest subregion more slowly than sound crosses a cell (§4.4:
   ĉ ≲ 0.07c on a 1024³ rung at 20 km/s).

Still open: the sink threshold inside the zoom (§4.5), the non-ideal
MHD term (7), and how short the RT–hydro coupling interval must be
while the front crosses the zoom (§4.4, the T4 test).
