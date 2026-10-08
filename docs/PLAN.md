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

## Phase 1d progress: production readiness (2026-10-08)

Done (`docs/production_run.md`, `docs/parallel.md`): reduced speed of light
`c_reduced` (C++ and both Fortran harnesses, bit-identical, I-front test);
OpenMP in RT, hydro, PM, FoF and summaries, bit-identical for any thread
count; `deterministic_parallel` ray tracing (thread-count independent);
ray-pool guard for F3 (`legacy_rt = 0`, both codes); 88-byte rays; F1 fix;
bit-identical restart; movie frames independent of the hydro step;
nightly test label; measured cost and timing table for the production run.

