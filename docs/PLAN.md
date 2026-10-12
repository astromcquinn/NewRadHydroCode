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
  Superseded 2026-10-10: **f_esc = 0.28**, the 128³ calibration at 0.1c
  with the unsplit hydro (`docs/production_run.md` §8); the 256³
  production run uses it with the 256³ L(M)
  (`data/uv_lf/am_params_256_10mpc.txt`).
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

## Stage D thread (2026-10-10): the ray-tracing RT on AMReX

* **Request (user):** Stage D of the roadmap on branch `stage-d-amrex`:
  one physics implementation (`include/radray/kernels/`), the AMReX RT
  level as a library component behind `rt_backend = amrex` (`legacy_rt = 0`
  only, AMReX optional), the R-struct encoding behind `ray_encoding`,
  the `amrex` test label, timings, and the CUDA target last and optional.
  Report: `docs/stage_d.md`.
* **Choices made in this thread** (open to revision):
  - `radray::Ray` is `rrk::RayP` plus the chain links; the links moved to
    the end of the struct, so the checkpoint format is version 5 and
    versions ≤ 4 are converted on reading.
  - The kernels take pointer views (`rrk::CellRef`) of the RT cell, so
    the serial SoA grid and the backend's AoS cell share one function; the
    rt_stat sums are returned per gas cell (`PhotoStats`) and added by
    the host in its own order, which keeps the serial sums bit-identical.
  - Stage D's MPI model: the hydro, PM, halo finder and sources run
    redundantly on every rank on the full serial state; only the RT is
    distributed, and the level's fields are gathered on rank 0 and
    broadcast once per hydro step (plus before a movie frame). Stage E
    distributes the rest.
  - The RT statistics (the τ_es history among them) are summed in 128-bit
    fixed point and reduced as 32-bit limbs, so they are rank- and
    thread-independent like the fields.
  - `spike/serial rt_serial` (now built against the library kernels, with
    the pool guard) is the bitwise reference of the backend; the spike's
    AMReX and Parthenon prototypes are frozen and no longer compile.
  - R-struct: 28 B of ray data (7 int32 components), the photon shape as
    5 × 8-bit log fractions at 1/16 octave with deterministic dithering,
    the direction offset only for merged rays (on-centre flag), r at 1/8
    cell with dr dropped, nu at 8 bits per band.
  - Round-off tests against the serial port use rtol 1e-3 on the fields
    (single I-front cells amplify the summation-order round-off to 1e-4)
    with the source-cell `dx` excepted, and exact integer ray statistics.
  - `recomb_rad`, the ray dumps and `ray_quant` stay serial-only for now.
* **Not done here:** the GPU build (no nvcc ≥ 12.2 in the cloud
  container; the CUDA target, device rate tables and the `ckpt` scripts
  are in place), the klone 256³ timing (`tools/hyak/amrex_scaling.slurm`,
  the user runs it), `ray_rotate` (its branch is not on the remote).


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

## Production-256 thread (2026-10-10)

* Settings from the thread brief: 256³, 10 Mpc/h, `c_reduced = 0.1`,
  `nray_max_factor = 16`, `legacy_rt = legacy_physics = 0`,
  `deterministic_parallel`, `rt_substep_carry`, `ray_rotate = 3`,
  `sfr_model = abundance`, log ξ_ion = 25.38, f_esc calibrated to
  τ_es = 0.054.
* **`hydro_scheme = unsplit`** for the production run and the 0.1c
  calibration runs. The brief left the choice open (`unsplit`, which fixes
  the cold-gas overheating, or `legacy`, which every earlier 256³ run
  used); the thread took the documented recommendation of
  `docs/hydro_unsplit.md` §10 and `docs/STATUS.md`. `ra_max` was left at
  the original 0.1: the time-step sensitivity of the cold gas only
  appears at z ≲ 2 (`docs/hydro_unsplit.md` §6). `make_production.py
  --hydro-scheme` writes the line.
* **f_esc = 0.28** (interpolated at 0.1c; `docs/production_run.md` §8).
* **User decision (2026-10-10, during the thread): f_esc should ideally
  make reionization end at z = 5.5**, with τ_es as the check against
  Planck. At 128³ and 0.1c that is f_esc ≈ 0.19–0.20 (x_V = 0.99 at
  z = 5.57 and τ_es = 0.050 for 0.2, 0.6σ below Planck). The 256³ run at
  0.28 is kept as the 256³ calibration point; the final production run
  uses the f_esc that puts the 256³ end of reionization at z = 5.5.
* **User (same day): use round numbers for f_esc (0.3, 0.2, …) from now
  on, and do not try to match the Planck τ_es**; the end of reionization
  is the only calibration target, τ_es is just reported. The 256³ box at
  f_esc = 0.28 reionizes much earlier than the 128³ one (x_V = 0.80
  against 0.58 at z = 7: the resolved 10⁸–10⁹ M☉ halos), so the
  production candidates are 256³ runs at f_esc = 0.2 and 0.1; the one
  ending nearest z = 5.5 is the production run (0.15 if neither is).
* **Result: the production run is f_esc = 0.2** (`docs/production_256.md`):
  x_V = 0.99 at z = 5.74 and 0.999 at z = 5.60, τ_es = 0.057; f_esc = 0.1
  does not finish by z = 5.4, and 0.15 reaches only x_V = 0.87 at z = 5.5.
* The 256³ run at that f_esc is made in the production configuration
  (movie every Δz = 0.025, checkpoint every 4 steps), so it is both the
  256³ τ_es check of the brief's step 4 and the production run of step 5
  if τ_es comes out at Planck; otherwise it is the calibration point for
  a corrected f_esc and a second run.
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
* **Calibrate on ṅ_ion (2026-10-11).** With f_esc = 0.23 the 40 Mpc/h 512³
  run matches the Becker+21 / Gaikwad+Zhu ṅ, ends reionization at z ≈ 5.6
  and gives τ_es = 0.052. But Γ_HI is ~3× and ⟨F_Lyα⟩ ~4× the data, and
  λ_mfp is far too long. The user: the Lyα forest measures will be missed
  until AMR resolves the sinks that set the mean free path, and ṅ is the
  most important constraint. No sub-grid absorber model for now.
* **Speed-of-light study (2026-10-11):** judge c_reduced on x_H alone (the
  x_HII history and maps against 1c), in a bigger box than 10 Mpc/h. The
  goal is the c that is acceptable as a function of x_H. Runs: 20 Mpc/h
  256³ at c = 0.05–1, plus a 10 Mpc/h 128³ box check (`runs/sol/`).
* **Movies always (2026-10-11):** every run saves frequent movie frames,
  even when full snapshots are rare; slab-only frames to keep the disk
  manageable (branch `movie-slab`).
* **PM refinement** (`pm_refine = 2`, a PM mesh and particles 2× the gas
  grid, for halo mass resolution) goes to a separate thread started by the
  user, on Hyak.

## Production defaults, moment RT, stored references on klone (2026-10-11/12)

* **Production runs use `hydro_scheme = unsplit`** (user decision,
  confirming the choice the production-256 thread made from the
  recommendation of `docs/hydro_unsplit.md` §10). `make_production.py`
  writes `hydro_scheme = unsplit` by default (`--hydro-scheme legacy` for
  the original). The code default stays `legacy`, which is what the
  bit-for-bit validation against the Fortran uses (ground rule 2).
* **Production runs use `redshift_per_distance = 1`** (user decision): the
  radiation is redshifted per distance travelled, not per time
  (`docs/production_run.md` §1). With ĉ < c, per time redshifts the photons
  that set Γ_HI c/ĉ times too much over their path to absorption (10× at
  0.1c, 3.3× at 0.3c), which matters after overlap, when λ_mfp is large.
  `make_production.py` writes `redshift_per_distance = 1` by default
  (`--redshift-per-distance 0` for the original). The code default stays 0
  (the original, ground rule 3); the flag is not tied to `legacy_physics`.
* Expected effect: a higher Γ_HI and ⟨F_Lyα⟩ after overlap at fixed f_esc;
  the ṅ_ion calibration (and so f_esc) is unchanged. Runs already started
  keep the setting they started with unless changed at a restart (the flag
  is not in the checkpoint).
* **M1 and OTVET are only for comparison** with the ray tracing (user
  decision): no second-order moment transport, and no moment hand-over in
  transparent regions (ROADMAP §4.5) for now.
* **OTVET stays as implemented (2026-10-12):** the explicit, time-dependent
  two-moment system with the optically thin Eddington tensor (minimum-image
  1/r² sum over the source cells), solved with the same GLF scheme as M1.
  No need to follow a specific published formulation.
* **Stored RT references on klone (2026-10-12): skip them there.** The 10
  tests that differ on klone's libm (last-place double sums, threshold
  rays; `docs/cluster.md` §3) carry the label `platform_ref`;
  `tools/hyak/build_hyak.sh` leaves it out. They are not deleted, not
  regenerated and not given a tolerance: CI (the platform the references
  were made on) keeps running them bit for bit.

## RT multicell thread (2026-10-11): rays crossing k cells per substep

* **Request (user):** test a coarse RT substep for the RT on AMR (Stage F)
  instead of sub-cycling fine levels in time: on the uniform grid, rays
  travel k RT cells per substep with a k times longer substep at the same
  ĉ. Branch `rt-multicell`; report `docs/rt_multicell.md`.
* **The flag:** `rt_cells_per_substep = k` (default 1 = the original,
  byte-identical; `params.txt` and `--set`; `radhydro` and `rt_driver`;
  in the log, the summary line and the checkpoint, version 6). `nstep_rt`
  counts the light crossings of k cells (carry and ≥ 1 guard on that
  count), `ray_segment` covers k cells in up to 3k+1 segments, in the shared
  kernels (`RTP::cells`), so the serial code and the AMReX backend agree.
* **Restrictions:** `legacy_rt = 0` only; `rt_method = rays` only (the
  moment solvers refuse k > 1); 1 ≤ k ≤ 16; light must cross at least k
  cells per hydro step (otherwise the ≥ 1 guard makes it faster than ĉ).
  `deterministic_parallel` domains d ≥ 2k+2. The AMReX backend supports
  k > 1: its ghost layers and local `Redistribute` grow to k+1 (round-off
  tests at k = 2, 4).
* **Generalisations k forced (this thread's choices, in the k > 1 scheme,
  no extra flag):** the photo-ionisation rate uses the crossing speed
  k·len_rt/dtsec_rt (without it Γ in ionised gas is k times low); the split
  radius clamp is [rmin, rmin+k]; for k > 4 sources emit at the nside whose
  split radius covers the first path (`rrk::emit_nside`), since a ray splits
  only once per substep; `redshift_per_distance` uses the k-cell distance.
* **Two new options found necessary (defaults = the original):**
  `photon_conserve = 1`, the B17 fix (`photo_cell`'s photon limit iterated to
  convergence; the original lets photon-starved front cells ionise without
  the photons, by an amount that grows as k²), and `ray_deposit = exposure`
  (a ray's photons are kept for the cells further along its path unless
  absorbed; the share deposit lets them jump past I-fronts by up to k − 1
  cells). Neither is tied to `legacy_physics`/`legacy_rt` (no Fortran
  mirror, per the 2026-10-08 rule).
* **Defaults stay the original;** the recommended set for k > 1 is
  `photon_conserve = 1` with `ray_deposit = exposure` (it passes the I-front
  and production tests where the share deposit does not, `docs/rt_multicell.md`
  §4, §7).
* **Finding:** the transport over k cells is accurate; the error is the k
  times longer substep (v_front > ĉ/k fails first: light-limited fronts at
  0.1c for k ≥ 4, dense gas at k = 8). On the uniform grid k > 1 is
  equivalent to fewer substeps and is not offered as a speed option.
* **Recommendation for the AMR RT (proposed in `docs/rt_multicell.md` §9;
  the user decides):** a mixture. Use the coarse substep with 2^l cells per
  substep on level l (with the exposure deposit, `photon_conserve` and
  per-level angular resolution). Add photo-chemistry sub-cycling in fine
  cells, and time sub-cycling of the transport only in patches where a
  front outruns the substep (flagged from cells deposited as thick that
  transmit most of their deposit).
## RT-methods 512³ thread (2026-10-11/12): Γ_HI from rays, M1 and OTVET

* **Request (user):** three 40 Mpc/h, 512³ runs differing only in
  `rt_method`, at one constant `c_reduced = 0.2`, f_esc = 0.23,
  `redshift_per_distance = 1`, small output (slab movie frames every
  Δz = 0.025, snapshots at z = 7–5), whole-box Γ_HI statistics in
  `frames.txt`, Γ tools, and the comparison of Γ_HI after overlap
  (`docs/rt_methods_512.md`).
* **User decisions during the thread:**
  - **store only the HI photoionization rate** (new field `gammaHI`, one
    float per cell, instead of the 3-component `gamma`);
  - **case A with transported recombination radiation for all three
    methods, and make it the default** — done as the `make_production.py`
    default (`recomb_case = A`, `recomb_rad = 1`, `recomb_rad_block = 4`);
    the C++ default stays the original (ground rule 3);
    `--recomb-case B --recomb-rad 0` gives case B;
  - **one more run with case B, no recombination radiation** (taken with
    the rays, the reference);
  - **the 3D power spectrum of Γ after overlap** (analytic work predicts a
    difference between the moment methods and the rays): extra Γ cubes at
    every hydro step after z ≈ 5.4 from bit-identical continuations, and
    `tools/analysis/gamma_power.py`.
* **Choices made in this thread** (open to revision): the RT starts at the
  first `zdata1`, so the snapshot list starts at z = 15 (first sources at
  z ≈ 13); `data_files = gas` (only `Gas/u_gas` at the snapshots, enough for
  `igm_stats.py` without the source-start λ_mfp) to save disk;
  `restart_keep = n` for "the last two checkpoints"; the ICs of
  `reion512_40mpc_f0.1` reused (identical parameters).
* **Block-emission clumps (2026-10-12).** The rays' block emission of
  recombination radiation made spurious self-shielded clumps with case A +
  RR (`docs/rt_methods_512.md` §6c). Fixed behind `recomb_rad_cell_onspot = τ`
  (C++ default 0; `make_production.py` default 1). **User: do not rerun the
  512³ case A + RR ray run**; its post-overlap residual HI, ρ_HI/x_HII power
  and λ_mfp stay flagged as affected in the document.

## M1 neutral source cells (2026-10-12, branch `m1-clumps`)

* **User asked** whether M1's many neutral absorbers after overlap are
  surprising, and asked for a check (512³ snapshots, then a 128³ test with
  the source map written).
* **Found:** an injection-order artifact, not the closure. The moment RT
  injected a substep's source and recombination photons before the
  transport, so most left the source cell before `photo_cell` absorbed
  anything; bright, dense source cells stayed neutral
  (`docs/rt_moments.md` §5a).
* **Choices made in this thread** (open to revision): `moment_inject = 1`
  (inject after the transport) behind a flag, C++ default 0 (ground
  rule 3); `make_production.py --moment-inject` default 1, as for
  `recomb_rad_cell_onspot`; `fields = ...,Luv` writes the source
  luminosity per cell (only when listed). The 512³ M1 and OTVET runs are
  not rerun; `docs/rt_methods_512.md` flags the affected results.

