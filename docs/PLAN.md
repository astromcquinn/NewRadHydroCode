# NewRadHydroCode — development plan

Goal: a cosmological radiative-transfer code in C++ whose RT is based on the
adaptive ray-tracing algorithm in `reference/Radhydro` (Trac & Cen 2007;
Trac, Cen & Loeb 2008), then MPI + GPU parallel, then compact ray encodings
and adaptive calculations.

## Reference algorithm (summary)

* RT grid of `Nrt^3` cells; gas state comes from the Eulerian grid.
* Each hydro step is subdivided into `nstep_rt` RT substeps so that rays move
  ~1 RT cell per substep (finite speed of light, `calc_timestep`).
* Rays carry `Nfreq=5` photon groups (`ngam`, effective frequency `nu`),
  position `x`, unit direction `dx`, radius from (effective) origin `r`,
  distance-to-next-split `dr`, HEALPix `nside`/`ipix` (NESTED scheme).
* Per RT substep (`rt_tools.f90: radiative_transfer`):
  1. `rt_step_1`: recompute cell optical depths; build the photon-count PDF of
     rays (`ngam_ray` thresholds); `ray_tracing_1` over every cell: sources emit
     `12*nside_init^2` rays (or deposit locally when optically thick); each ray
     is cast (`ray_cast_1`), split 1->4 (`ray_cast_split`, when `dr<0` and the
     ray is bright), removed, or binned into one of 48 `nside=2` cones and
     merged with the other rays in that cone and cell (`ray_cast_merge`).
     Casting deposits photons onto the RT grid along up to 3 segments.
  2. `rt_step_2`: `radiative_physics` — photoionisation/heating/chemistry per
     cell, which sets the transmitted fraction `fgam`.
  3. `rt_step_3`: `ray_tracing_3` — advance rays one cell, attenuate by `fgam`,
     delete faint rays, relink the per-cell ray chains.
* `post_rt`/`ray_redshifting` handle cosmological redshifting of rays.

## Phases

1. **Faithful port + validation (now).**
   * `reference/harness/`: build the original Fortran RT with gfortran as an
     RT-only driver (no hydro/PM; fixed gas field + source list), single thread.
   * `src/`, `include/radray/`: serial C++17 port of the RT module
     (ray tracing, RT physics/rates, cosmology/units needed by RT).
   * `tests/problems/`: shared test problems; `tools/compare/`: comparison.
   * Acceptance: C++ agrees with the Fortran harness on every test problem
     (see `docs/validation_spec.md`).
2. **Parallelism.** Kokkos for CPU/GPU portability; MPI domain decomposition
   with ray exchange across rank boundaries; replace linked-list chains with
   sort-by-cell arrays (GPU friendly).
3. **Ray encoding & adaptivity.** Compact ray representation (HEALPix
   hierarchy index + quantised offsets + packed photon groups), error-driven
   split/merge criteria, adaptive RT substepping.
4. **Coupling** to hydro/gravity (port or couple to an external code).

## Phase 1b (decided 2026-10-08)

* **RT fixes behind a flag.** `legacy_rt = 1` reproduces the original
  bit-for-bit, in both the Fortran harness and the C++ code. `legacy_rt = 0`
  turns on these bookkeeping fixes (see `docs/port_notes.md` §5, §5a):
  * B1: freed rays are zeroed, so they no longer count toward the thresholds;
  * B2: `ngam_ray` is recomputed every substep, with no stale values;
  * B3: `rt_grid%dx` is reset every substep;
  * B4: `tau` is updated in all cells;
  * B8: casts and advances may use up to 4 segments.

  The Fortran harness and the C++ code must agree bit-for-bit in both modes.
  The B12 physics oddities are untouched for now, pending a decision. B5–B7 are
  fixed unconditionally, since they are memory-safety bugs. B9, the int32
  overflow that only matters for N ≥ 323, is fixed only when `legacy_rt = 0`,
  so legacy mode stays exactly the original.
* **N-body + hydro port.** `pm_tools.f90` (PM gravity for dark matter and gas),
  `fft_tools.f90` and `hydro_tools.f90` (with `-DMF`, as the original
  Makefile sets). Both the Fortran harness and the C++ code use FFTW3
  single precision in place of MKL DFTI, so they can match bit-for-bit.
  Each part gets its own Fortran reference driver and its own validation
  problems, and then the full `radhydro` loop is integrated (Phase 1c).
* **Physics fixes behind a flag** (decided 2026-10-08). `legacy_physics = 1`
  keeps the original physics. `legacy_physics = 0` fixes:
  * the case-B recombination cooling coefficients, which currently use the
    case-A formulas;
  * the missing HeI collisional-excitation cooling;
  * the unset `coef_cool(15)`.

  Recombination radiation stays off for now. Both Fortran and C++ get the
  flag. This is scheduled after the `legacy_rt` work, because both touch the
  RT rate code. **Done** (`docs/port_notes.md` §5b, harness patch 9). The
  flag also fixes three more errors found in the rate code: the HeIII
  case B coefficient typo, collisional-ionisation cooling using kT
  instead of the ionisation potential, and the missing HeII dielectronic
  recombination. `coef_cool(15)` is now the dielectronic recombination
  cooling. The recommended configuration is `legacy_rt = 0`,
  `legacy_physics = 0`.
* **On-the-fly halo finder** (decided 2026-10-08): FoF on the DM particles, a
  source model reconstructed from Trac & Cen 2007 and Trac, Cen & Loeb 2008,
  and a reader for the original catalog format to keep legacy mode working.

## Phase 1d (requested 2026-10-08, after 1c and the physics fixes)

1. **Validate the halo mass function** from full runs against Sheth–Tormen
   (and Tinker et al.), as a gate before any science run.
2. **Hydro test suite:** Sod, Sedov–Taylor blast, Zel'dovich pancake with
   gas (Bryan et al. 1995), and a Santa Barbara-like cluster or
   adiabatic cosmological box if feasible. Compare against analytic
   solutions and the literature.
3. **Abundance-matched sources:** calibrate the halo UV luminosity
   (L_UV–M_halo) so that the simulated UV luminosity functions match the
   observed high-z LFs (z ≈ 6–10 HST, e.g. Bouwens et al. 2021; z > 10
   JWST). Ionising emissivity follows from L_UV through ξ_ion·f_esc;
   f_esc is calibrated to τ_es (Planck) and the end of reionisation.
4. **Movie of a full reionisation run:** the largest run this machine can do
   efficiently (4 cores, 15 GB RAM, no GPU). This first needs OpenMP in the
   hot loops and a ray-memory budget (`nray_max_factor`). The size gets
   decided from measured timings.

## Phase 1c (done 2026-10-08): the full radhydro loop

* One FFT layer, one `calc_timestep`, one `summarize_gas`, one set of gas
  domains, in C++ and in the Fortran harnesses (`reference/harness_common`).
* `apps/radhydro`: the complete `radhydro.f90` loop with run-time N, Ndm,
  box, zdata, source mode (legacy catalogs or on-the-fly FoF), `legacy_rt`,
  `legacy_physics`, `nray_max_factor`.
* `reference/harness_full`: gfortran build of the original program; the
  on-the-fly mode calls the C++ halo finder through a C interface.
* Bit-identical Fortran vs C++ on all validation cases
  (`tests/problems_full`), both source modes, both flag settings.
* Cost: RT dominates (rays × light-crossing substeps); see
  `docs/full_run.md` §5 for the measured cost model and the extrapolation
  that sizes the production run.

### Production-run decisions (2026-10-08)

* **Box:** 10 Mpc/h for the first reionisation movie run (changed from
  20 Mpc/h on 2026-10-08).
  * RT validation problems so far use 5 Mpc/h.
  * Full-loop tests and the demo use 4 Mpc/h.
* **Reduced speed of light:** `c_reduced` parameter, mapped onto the
  original's `cscale`.
  * Use 0.01c. The user accepts that light is too slow for now: at 0.01c it
    covers only about 15–20 comoving Mpc/h between z = 12 and 6, so overlap
    and the post-overlap background will be biased.
* **Ray buffer:** `nray_max_factor = 16`, giving about 5–10 live rays per
  cell (live rays saturate at 0.3–0.6 × Nray_max).
* **Resolution:** Ngas = Ndm.
  * 128³ is the default.
  * 192³ fits in memory if the timings allow.
* **Sources:** abundance-match the resolved halos to the observed UV
  luminosity functions. No sub-grid sources for unresolved galaxies for now
  (user decision, 2026-10-08). Only halos above the resolution limit emit,
  about 7×10⁸ M☉/h for 128³ particles in 10 Mpc/h.

### Abundance-matched sources (decided 2026-10-08; docs/abundance_matching.md)

* **Model:** a new `sfr_model = abundance`. The existing models stay.
* **LFs:** use recent JWST LFs in preference to Bouwens+21. At z ≥ 9 that
  is Donnan et al. 2024; the user will supply z ≈ 6–8 papers.
* **Shape:** the shape of L_UV(M_halo) is fitted at z ≤ 9, with the faint
  end extrapolated when matching (down to M_UV = −12).
* **High z:** at z > 9.5 only the normalisation is fitted, to the luminous
  galaxies observed.
* **Beyond the data:** L(M) is held fixed above the highest LF redshift.
* **Mass function:** match to the simulated FoF mass function (decided
  2026-10-09; this replaces Tinker), so the box's LF matches the
  observations.
* **f_esc:** start the production run with f_esc = 0.5 (2026-10-09).
* **Dust:** no correction.
* **ξ_ion:** a JWST-era value, log ξ_ion = 25.38, with a single f_esc
  calibrated to τ_es.

### Long-term decisions (2026-10-08; details in docs/ROADMAP.md)

* Follow the stage order in `docs/ROADMAP.md` §7. MPI comes after the
  framework spike. Single-node OpenMP runs on Hyak cover testing until then.
* **Cluster:** UW Hyak `klone`, about 15 group nodes, high memory, no queue
  time limits. GPU count still unknown.
* **Dense gas under AMR:** keep cold dense gas by default; the
  star-formation sink is an optional flag (revised 2026-10-08). Whether
  cooling stays on in dense gas is still open.
* **Plane-parallel RT mode** (rays cast from box faces): requested, see
  `docs/ROADMAP.md`.
* **RT accuracy target:** a few percent in astrophysical outputs. Ray-level
  precision is set empirically against acceptance metrics.

## Phase 1d progress: production readiness (2026-10-08)

Done (`docs/production_run.md`, `docs/parallel.md`): reduced speed of light
`c_reduced` (C++ and both Fortran harnesses, bit-identical, I-front test);
OpenMP in RT, hydro, PM, FoF and summaries, bit-identical for any thread
count; `deterministic_parallel` ray tracing (thread-count independent);
ray-pool guard for F3 (`legacy_rt = 0`, both codes); 88-byte rays; F1 fix;
bit-identical restart; movie frames independent of the hydro step;
nightly test label; measured cost and timing table for the production run.


## Fortran harness scope (decided 2026-10-08)

* The Fortran harnesses only have to reproduce the **original** code, so the
  C++ can be compared with it at the default settings.
* New C++ options that change results (fixes, research instruments such as
  the `ray_quant_*` lossy-precision options of the ray-compression study)
  go behind a flag whose default is the original behaviour. They need **no**
  mirrored implementation in Fortran. The already mirrored `legacy_rt` and
  `legacy_physics` stay as they are.

## Hydro solver thread (2026-10-09): the unsplit scheme behind `hydro_scheme`

Decided with the hydro-suite finding that the original moving-frame TVD
scheme spuriously heats cold supersonic gas (`docs/hydro_tests.md` §6, §7,
§9); details and results in `docs/hydro_unsplit.md`.

* **`hydro_scheme = legacy | unsplit`** (default `legacy`, which stays
  bit-identical to the Fortran in every test). `unsplit` is the unsplit
  MUSCL-Hancock + HLLC solver of the framework spike
  (`include/radray/kernels/hydro_mh.hpp`, header-only and
  framework-independent so the AMReX thread compiles the same file), with
  the **dual-energy formulation** of Bryan et al. 1995: a separately
  advected internal-energy field (`GasGrid::e`), updated with the −p∇·v
  work from the Riemann solver's interface velocities, and the η₁/η₂
  selection (`dual_eta1`, `dual_eta2`; `dual_energy = 0` turns it off).
  The entropy-conserving variant (Ryu et al. 1993) was not built: the
  dual-energy version passed the tests.
* **Loop placement:** the six sweeps around `particle_mesh` become one
  sub-cycled unsplit step over `dt_sim/2`, the PM kick, and a second one;
  same call sites, same `calc_timestep` (so the same `nstep_rt`). The PM
  kick of the gas acts on the conserved lab-frame variables with the
  thermal energy unchanged (`pm::particle_mesh_unsplit`).
* **RT coupling rule:** the RT's `photo_cell` takes the thermal energy from
  the internal-energy field when it exists and writes its heating/cooling
  back to both `u(5)` and the field (the physics coupling point). A second,
  diagnostic RT edit reads the field for the RT grid temperature
  (`average_gas`/`gas_state_from_avg`), as does `summarize_gas`, because
  the float `u` cannot hold a cold supersonic cell's thermal energy; both
  are no-ops (same arithmetic) with the legacy scheme. **Approved by the
  user (2026-10-09).**
* **Constants:** η₁ = 10⁻³ and η₂ = 0.1 confirmed against the Enzo method
  paper (Bryan et al. 2014, §4.1.1, eqs. 44–45) and the total-energy reset
  against the Nyx paper (Almgren et al. 2013, §3.2), both supplied by the
  user on 2026-10-09 (Bryan et al. 1995 itself is not accessible to the
  user either), and their exact forms tested on the suite problems
  (`docs/hydro_unsplit.md` §2). Deliberate deviations, measured there: the
  η₁ test and the η₂ maximum use the ±2 neighbourhood along all axes
  (Enzo: cell-local η₁, ±1 along the sweep for η₂), and the work term uses
  the time-centred pressure (Enzo: old time). All are parameters.
* **Validation:** the hydro suite with `--scheme unsplit` and its own
  tolerances, the restart and thread-count tests, and one coupled 64³ run
  per scheme. No Fortran harness copy (new physics beyond the original
  needs none, see the Fortran harness scope above).

## Framework and ACCESS decisions (2026-10-10)

* **AMReX is the framework for Stages D–G** (approved by the user after the
  Hyak CPU and GPU benchmarks, `docs/framework_decision.md`). The migration
  plan is §6 of that document.
* **No NSF ACCESS proposal for now.** The 128³ RT benchmark put one H200
  at 5.0–5.9 `cpu-g2` nodes and one A100 at 3.3–3.8
  (`docs/framework_decision.md` §4.4). GPU portability stays a requirement
  (`docs/ROADMAP.md`, NSF ACCESS); the question can be reopened with the
  AMReX CUDA numbers of Stage D.

## Merge thread (2026-10-09): no decisions

The survey (PR #16) and the merges of `hyak-setup` (PR #17), `movies`
(PR #12) and `production-run` (PR #7) into `main` followed the survey's
plan in `docs/STATUS.md`. **No decisions were taken in this thread.** The
open ones are listed in `docs/STATUS.md` ("Decisions needed from the
user"): AMReX for Stage D, the production `c_reduced` (0.1 recommended)
and N (256), `hydro_scheme` for production, the stored-RT-reference
policy on klone, and deleting the merged branches.

## Moment RT thread (2026-10-10)

* **Request (user):** "implement OTVET and M1 closure methods for the RT as
  an alternative compared to the ray tracing".
* **Done as:** a new runtime option `rt_method = rays | m1 | otvet`
  (default `rays`, bit-identical to the original), following ground rule 3
  (new behaviour behind a flag whose default is the original). No Fortran
  mirror, per the 2026-10-08 decision. The physics coupling
  (`radiative_physics`) is shared with the rays, unchanged.
* **Choices made in this thread** (open to revision; `docs/rt_moments.md`):
  - one explicit two-moment GLF solver for both closures, as in RAMSES-RT;
  - OTVET in the time-dependent form, not the implicit diffusion form of
    Petkova & Springel 2009;
  - first order in space and time;
  - `moment_nsub = 4` (C = 1/4; 1/3 decouples odd and even cells);
  - the OTVET source-cell self-weight 1/(0.5 cell)²;
  - `auto` switching between the direct sum and the FFT for the tensor.

## B16 fix thread (2026-10-10): `ray_rotate`

* **The user asked for the per-substep random rotation** of the source
  ray frame (after the question of its cost). It is behind the new flag
  `ray_rotate` (default 0 = the original). It is not tied to `legacy_rt`:
  the original stays the default, and the Fortran harnesses need no mirror
  (the 2026-10-08 rule for new flags).
* Implementation choice: rotate at emission only, and keep `ipix` in the
  grid frame (the pixel containing the rotated direction). This avoids a
  per-ray frame id and any re-indexing of live rays.
* **Production uses `ray_rotate = 3`** (user decision 2026-10-10), written
  by `make_production.py`. The user's caveat: the single-source tests
  exaggerate the effect for a full cosmological run, where many
  overlapping sources fill each other's coverage gaps, so the 3–10 %
  single-thread ray-tracing cost buys less there than the uniform tests
  suggest.

## Emissivity calibration thread (2026-10-10): 40 Mpc/h, 512³

* **Production runs end at z = 5** (Lyα forest comparison); the
  `make_production.py` default since 2026-10-10.
* **Calibrate the sources on bigger boxes against the observed IGM at
  z = 5–6**, not on the source models from the literature ("I don't think
  they are based on much"). Box: **40 Mpc/h, 512³** (m_p = 5.1e7 M☉, the
  same mass and cell resolution as the 10 Mpc/h 128³ box, so the 128³
  FoF-corrected abundance matching `data/uv_lf/am_params_128_10mpc.txt`
  applies).
* **Data** (papers supplied by the user; `data/igm/igm_z5_6.txt` on branch
  `fesc-z`): Becker+21 (λ_mfp, ṅ), Zhu+23 (λ_mfp), Bosman+22 (⟨F_Lyα⟩) and
  Gaikwad+23 **Γ_HI only**: the user is skeptical of Gaikwad's mean free path
  inferences (and so of their ṅ), not of their Γ.
* **Speed of light:** c_reduced = 0.2 is acceptable, but it must stay at
  that value through the end of reionization (z ≈ 5.5) before it is
  changed (raised for z = 5.5 → 5). Which reduced c is optimal is worth
  measuring. **Revised the same day: the first 512³ run uses
  c_reduced = 0.1 to z = 6.5, then 0.3 through overlap and on to z = 5**.
  The user first suggested returning to 0.1c after overlap. Then two effects
  of a low ĉ after overlap were pointed out: the slow relaxation of Γ
  (λ/ĉ ≈ 100–300 Myr) and redshifting per unit time, i.e. c/ĉ too much per
  distance (`ray_redshifting`). The user agreed to stay at 0.3c, and to a new
  flag that redshifts per distance travelled (default off).
* **IGM diagnostics** (`tools/analysis/igm_stats.py`, `docs/igm_diagnostics.md`):
  f_esc is calibrated on ⟨F⟩, Γ_HI and ṅ. λ_mfp is a check only: it comes out
  too long at fixed Γ in these boxes (unresolved absorbers, box ≲ λ).
* **f_esc:** with the bigger box err on the side of a smaller escape
  fraction, perhaps 0.1; the first 512³ run uses f_esc = 0.1 (constant).
