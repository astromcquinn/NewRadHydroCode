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


## Stage E2 thread (2026-10-12): the PM gravity on AMReX

* **Request (user):** Stage E2 on branch `stage-e-pm`: the PM gravity on
  AMReX behind `pm_backend = amrex` (default `serial`, bit-for-bit): the
  dark matter in an AMReX `ParticleContainer` (single-precision build,
  AMReX 26.10 as Stage D), the CIC deposit onto a MultiFab of the PM mesh
  (`pm_refine`: a BoxArray of its own), the gas density from the host
  grid, the Poisson solve with `amrex::FFT`, the force interpolation, the
  DM kick and the gas kick returned to the host (`calc_kick` and
  `calc_kick_unsplit`); the level copies the host state in at the start of
  the PM step and gathers the results back (the hydro and the RT still
  hold the full state on every rank); `src/amrex/pm_level.cpp` with a stub;
  minimal, additive edits to the shared files; the hydro files untouched
  (the E1 thread's). Validation under the label `amrex_pm`; timings at
  128³ and 256³; `docs/stage_e_pm.md`.
* **Choices made in this thread** (open to revision; `docs/stage_e_pm.md`):
  - the two drifts of `particle_mesh` stay on the host (per-particle
    arithmetic, identical on every rank), so the particles never move while
    on the level and no `Redistribute` is needed; the level holds no state
    between steps, so checkpoints and restarts are unchanged;
  - the gather back is a bitwise-OR MPI all-reduce of disjoint pieces
    (exact; a sum would turn −0.0 into +0.0);
  - the CIC deposit is int64 fixed point at 2⁻³⁸ (range 3.4e7 per mesh
    cell), the kick statistics 128-bit at 2⁻⁵⁶;
  - **`amrex_pm_fft = pencil` is the default**: per-line 1-D FFTW
    transforms on slab layouts, rank- and thread-independent by
    construction; `amrex_pm_fft = amrex` (`amrex::FFT::R2C`) is available
    and was bitwise rank/thread-independent in every test here, but its
    determinism rests on FFTW's planner choices for `plan_many` transforms,
    which is not guaranteed;
  - the PM kernels (CIC, force stencil, Green's function) moved into
    `include/radray/kernels/pm_kernels.hpp`, called by the serial code
    (byte-identical) and the level, as Stage D did for the RT.

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

## Refined PM mesh (2026-10-11): `pm_refine`

* **Goal: resolve halos of 1.3×10⁸ M☉ in the 40 Mpc/h, 512³ runs**, i.e.
  20 particles of 6.4×10⁶ M☉: Ndm = 1024 under the 512³ gas grid.
* **`pm_refine = 2`** gives a PM mesh of 2·Ngas with Ndm = 2·Ngas
  particles; **the default 1 stays byte-identical**. The gas density is
  spread onto the fine mesh and the accelerations are averaged back,
  mirroring each other (adjoint), in both `calc_kick` and
  `calc_kick_unsplit`. **Threaded FFTW.** `make_ics.py` and the checkpoint
  must handle Ndm > Ngas.
* **Validation asked for:** a single-particle force test, momentum
  conservation, a Zel'dovich test, and the halo mass function at
  40 Mpc/h with 1024³ DM against the 10 Mpc/h 256³ catalogs
  (`runs/hmf256`, the same particle mass and PM mesh spacing).
* Thread choices (not user decisions; `docs/pm_refine.md`): the gas
  spread is piecewise constant (each gas cell's mass evenly over its r³
  mesh cells), the transpose of the force average; the 1024³ ICs keep the
  large-scale modes of the 512³ production ICs (`--phases-from-ndm 512`,
  seed 7); the validation run uses the production hydro (`unsplit`); FFT
  threads are capped at 32 (faster, results identical).
* **FoF less often than every step** (user, 2026-10-11, after the cost
  of FoF on 10⁹ particles: 32 s per hydro step). Implemented as
  `fof_dt_myr` (default 0 = every hydro step; `docs/pm_refine.md` §5).
  (Superseded the same day: FoF stays every hydro step, below.)
* The pm_refine = 2 abundance matching uses a FoF correction fitted to
  the 40 Mpc/h 1024³ catalogs (the documented rule: refit at each
  resolution); `data/uv_lf/am_params_1024_40mpc_r2.txt`. (Superseded the
  same day by the Sheth–Tormen refit, below.)
* **FoF every hydro step** (user, 2026-10-11, revising the above): the
  hydro steps are 20–30 Myr apart until overlap, so `fof_dt_myr` would
  save only minutes; production keeps the default 0 (every hydro step).
  The option stays.
* **Abundance matching to an analytic mass function** (user, 2026-10-11,
  replacing the 2026-10-09 decision to match to the box's own FoF mass
  function): "calibrate to it so that the calibration doesn't change
  depending on the box". Plain ST fits the 40 Mpc/h FoF halos badly
  (FoF/ST 1.1–1.9 above 500 particles), so the user chose ST with its
  coefficients refit to the well-resolved FoF halos (option 2):
  A = 0.3458, a = 0.6164, p = 0.3, `data/uv_lf/am_params_st_fof.txt`
  (`docs/pm_refine.md` §6). f_esc has to be recalibrated with it.

## Stage E prerequisites on klone (2026-10-12)

* **What's next (assessment for the user):** Stage E of `docs/ROADMAP.md`
  §7, MPI for hydro, PM, halos/sources and I/O on AMReX (the RT is
  already distributed by Stage D), then Stage F, AMR. Prerequisites: the
  klone 256³ timing of the Stage D backend (`tools/hyak/amrex_scaling.slurm`,
  never run), the multi-node hydro slowdown of `docs/framework_decision.md`
  §3.3, and a production configuration the backend accepts.
* **User decision: case B is the production default again**
  (`make_production.py --recomb-case B --recomb-rad 0` are the defaults),
  so that the production configuration runs on `rt_backend = amrex` and
  Stage E does not wait for `recomb_rad` on the backend. Case A with
  transported recombination radiation stays available with
  `--recomb-case A --recomb-rad 1` (its `recomb_rad_block` and
  `recomb_rad_cell_onspot` defaults are unchanged).
* The Hyak scripts' module loading was broken under `set -u` (Lmod reads
  unset variables); `radray_module_load` in `tools/hyak/env.sh` fixes it.
* Done on klone the same day (`docs/stage_d.md` §4.4, `docs/STATUS.md`):
  the backend built and `ctest -L amrex` green; 256³ on one node and
  512³ on 2 and 4 nodes timed against the serial code, fields to
  round-off. Choices made without the user: `reion256_c01` has no
  checkpoint, so the 256³ timing uses `reion256_f0.2`; the 512³ run is
  the case B checkpoint `reion512_40mpc_f0.23`; 512³ at 24 ranks × 8
  threads per node runs out of memory (the redundant serial state is
  about 100 GB per rank), so the field broadcast and the restart scatter
  were made piecewise (bitwise identical results) and the 512³ layout is
  12 ranks × 16 threads per node, on 2 nodes or more. The next backend
  optimisation target is the deposit `SumBoundary` (25 % of the RT at
  512³); Stage E removes the redundant state.
* Later the same day (klone, after the Stage E1 merge): the build of
  `5c0c974` passes `ctest -L amrex` 81/81 and the fast suite 125/125;
  256³ re-timed on 1, 2 and 4 nodes (RT 156.8 / 100.2 / 76.8 s against
  the serial 144.4 s), the 512³ 2- and 4-node results of the morning
  kept (the RT level is unchanged). The multi-node hydro slowdown
  measured (`docs/stage_e_hydro.md` §5): it is the UCX rendezvous path
  of the klone OpenMPI, not the ghost volume. Choices made without the
  user: `tools/hyak/env.sh` exports `UCX_RNDV_THRESH=8m
  UCX_ZCOPY_THRESH=inf` for every AMReX job (`radray_ucx_eager`,
  `RADRAY_UCX_RNDV_THRESH` to change it, `RADRAY_UCX_EAGER=0` to turn
  off: eager for everything fixed the hydro but slowed the 512³ RT by
  8 %, eager below 8 MB fixes the hydro and speeds the RT up by 10–15 %); the production box size stays 32³;
  the ob1/TCP path (as fast, but a bus error at exit) was not pursued;
  512³ on one node with 8 fat ranks was not tried.

## Stage E4 (2026-10-13): load balancing

* Done after PR #34 (Stage E3) was merged; branch `stage-e-balance`.
* **Choices made in this thread** (open to revision; `docs/stage_e_state.md` §5):
  - The cost is cells + 10 × live rays per box (the measured ratio of a
    ray-substep to a cell-substep at 512³); `amrex_balance_ray_weight`.
  - The space-filling curve mapping by default (locality for the
    exchanges), knapsack as an option; a remap only when the efficiency
    gains 0.05 (`amrex_balance_min_gain`; negative remaps every step,
    which the tests use).
  - The mapping is rebuilt once per RT step before the step, on the ray
    counts at that moment; not inside the step, and the hydro and PM
    work is not weighed.
  - The checkpoint's ray order became canonical (box order, cell, bytes)
    so that a rebalanced run's checkpoints equal the fixed mapping's; the
    on-disk format is unchanged.
  - Production runs on the backend: `make_production.py --backend amrex`
    puts the backend configuration into `params.txt` (so every restart
    runs it), and `reion_production.slurm` launches MPI ranks when it
    sees `rt_backend = amrex` there; the layout is the job's `--ntasks`
    × `--cpus-per-task` (24 × 8 per node by the examples). Validated at
    128³ against the serial code and through a time-limit kill and
    resubmission (`docs/stage_e_state.md` §6).
  - One 256³ run on 4 nodes hung inside the RT step with identical data
    up to the hang and did not reproduce (five further runs of the case
    on other nodes, all bitwise identical); the backend now writes every
    rank's call stack on SIGTERM (`Backtrace.<rank>`, resolved by
    `tools/hyak/backtrace_lines.py`) so that the next hang can be read,
    and the job scripts run in the run directory. Not pursued further
    without a second occurrence.

## Stage E2 on klone (2026-10-13): the PM backend's validation and scaling

* **Request (user):** on branch `stage-e-pm-klone`, build `main` with
  AMReX 26.10 + FFT and run the `amrex_pm`, `amrex`, `amrex_state`, `pm`
  labels and the fast suite; the 64³ production-configuration comparison
  of the three backends with `amrex_state = 1` against the serial code
  (`docs/stage_e_pm.md` §3), as a nightly test if it fits; the PM backend
  at `pm_refine = 2` (`hmf512_40mpc_r2` restarted for 3 steps) on 1, 2, 4
  cpu-g2 nodes against the serial code, both FFTs, box sizes 32 and 64,
  bitwise across all of them; if the amrex FFT ever differs, say so and
  keep `pencil` the default. Nothing may change results; the serial PM
  stays byte-identical.
* **Choices made in this thread** (open to revision; `docs/stage_e_pm.md`
  §3.1, §4.1):
  - The 64³ comparison is the nightly test `amrex_state_full64_vs_serial`
    (95 s), compared as `amrex_full64_vs_serial` does for the RT, with
    tolerances set from the klone measurement (histories and τ_es 1e-3,
    per-cell fHI rms 0.15, T 1e-2, ρ 5e-4). The 5e-2 per-cell fHI
    difference is the rays' Γ noise reshuffled by the PM's round-off
    (fHI difference = Γ difference, correlation −0.9998), not an error.
  - `--set step_timers=1` prints the backends' cumulative phase timers
    after every step (output only); the PM level got a separate `drift`
    timer.
  - `pm::init_pm(..., host_arrays)`: with `amrex_state = 1` only the I/O
    rank allocates the host PM arrays and runs the host initial FFT
    (memory only; the 1024³ runs went OOM without it).
  - The amrex FFT's planner reset is guarded by
    `fftwf_planner_nthreads() > 1` (klone's AMReX never initialises
    FFTW's threads; the unguarded reset crashed at exit). Results
    unchanged.
  - The scaling script compares each FFT against its own 1-node
    reference; `pencil` stays the default (the amrex FFT was also bitwise
    across nodes and boxes, but is slower on klone and has no guarantee
    by construction).
  - No new user decisions. Open for the user: the FoF and P(k) on rank 0
    make the full step 3× slower than the serial code at this size
    (§4.1); a distributed halo finder (E3's open item) or a less
    frequent FoF is the remedy.


## Stage E3 (2026-10-13): the state on the levels

* **Request (user):** merge the branches, understand the hydro threshold
  oddity, try the 8 × 24 layout, then continue through the E stages.
* Done: PR #32 (Stage E1 and the klone runs) and PR #33 (Stage E2,
  opened from this thread since E3 needs both) merged into `main`; the
  threshold oddity is the gather's 16–32 MB broadcast blocks going eager
  above 16 MB (`docs/stage_e_hydro.md` §5), 8 MB stays; 8 × 24 fits one
  node at 512³ but is 1.26× slower than the serial code.
* **Choices made in this thread** (open to revision; `docs/stage_e_state.md`):
  - `amrex_state = 1` is a switch on top of the three backends (default
    0); the host grids become rank 0's mirror rather than disappearing,
    so that every file format and every statistic stays bitwise the
    host-state run's, which is the test (`ctest -L amrex_state`).
  - The halo finder stays serial on rank 0's particle mirror (gathered
    by id in chunks); a distributed FoF is left for later.
  - The per-step statistics gather the gas (40 B per cell) and five RT
    fields (20 B per cell) to rank 0 every step; the on-level 128-bit
    alternative is deferred until round-off there is accepted.
  - P(k) stays the host code on rank 0's mirror.
  - The checkpoint format is unchanged (version 7); the ray gather is
    chunked.
  - One decomposition: `amrex_hydro_max_grid_size` 0 (= `amrex_max_grid_size`),
    `amrex_pm_max_grid_size` 0 (= `pm_refine` × that).

## Stage E1 (2026-10-12): the unsplit hydro on AMReX

* **Request (user):** start Stage E with the hydro (E1) in this thread;
  E2 (PM on AMReX) in a parallel thread with the brief given in the
  session; the Hyak thread builds the backend and runs the Stage D
  scaling job.
* **Choices made in this thread** (open to revision; `docs/stage_e_hydro.md`):
  - `hydro_backend = serial | amrex` is a switch independent of
    `rt_backend`; `amrex` needs `hydro_scheme = unsplit` (the legacy
    split scheme never goes on the backend).
  - The level is bitwise identical to the serial unsplit step, not merely
    round-off close: the serial step is stencil-local with exact
    reductions, so the test is `compare.py --rtol 0 --atol-frac 0`.
  - The face and update code moved into the kernel header
    (`hydro_mh.hpp`) rather than being duplicated.
  - Stage E1 keeps the full host state on every rank (copy in per box,
    `MPI_Allgatherv` out), like Stage D's RT; the state moves onto the
    level for good when E2–E4 are in.
  - Flux temporaries are per box (never communicated); `Q` is computed on
    the grown box rather than exchanged.
  - `amrex_hydro_max_grid_size` defaults to 32 (the spike's 16³ boxes are
    the suspect of the multi-node slowdown; §5 of the doc).
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

## OTVET-closure thread (2026-10-12): a better Eddington tensor for `rt_method = otvet`

* **Request (user):** branch `otvet-closure`; improve OTVET's Eddington tensor,
  whose post-overlap P_Γ is 2–9× the rays' at k ≲ 6 h/Mpc (r ≥ 0.95), with a Γ
  floor at 0.62⟨Γ⟩ and low peaks (`docs/rt_methods_512.md`), while keeping its
  good pre-overlap morphology. An attenuated kernel exp(−r/λ_mfp) is **the
  baseline to beat**, not the answer. Steps: (1) ground truth — a ray-tracer
  diagnostic of the true second-moment tensor per cell, compared on 256³ runs
  with OTVET's h and M1; (2) brainstorm candidate closures with reasoning,
  cost and memory, **show the list and ranking to the user before
  implementing**; (3) implement the best two or three behind flags, measure
  against the rays' tensor, then full runs (`gamma_power.py`,
  `reion_power.py`, `rt_methods_runs.py`). Report: `docs/otvet_closure.md`.
* **Choices made in this thread** (open to revision):
  - the diagnostic is `tensor_dump = 1` (radhydro only, serial RT backend):
    on the zdata1 snapshot steps the Γ_HI-weighted moments S0, S1, S2 of the
    radiation, summed over the substeps (also over the odd substeps alone, for
    the ray noise), with the thin OTVET tensor of the step's sources and the
    source list. For the moment methods the same moments of their own N, F
    and closure D N. Off by default (byte-identical, a null pointer);
  - the ground-truth runs are new 256³ runs (`$RUNS/otvcl256{,_m1,_otvet}`)
    at **20 Mpc/h**, the 512³ run's cell size (78 kpc/h), mass resolution and
    configuration (case A + RR with `recomb_rad_cell_onspot = 1`, ĉ = 0.2,
    f_esc = 0.23, `redshift_per_distance`, unsplit), snapshots and tensor
    dumps at z = 15, 10, 9, 8.5, …, 5. The existing 256³ runs (10 Mpc/h,
    0.1c) kept only their final checkpoints, so any 256³ ground truth needed
    reruns;
  - fixed on the way (no flag: it only enables a restart that used to
    throw): a checkpoint written before the RT starts holds empty
    recombination buffers, and restarting it with `recomb_rad = 1` failed.
  - a priori choices for the cone (numpy prototypes against the rays'
    tensor, `docs/otvet_closure.md` §3):
    - nside 4 rather than 2: it fixes the cone's small-scale force;
    - receivers 4× coarser: as good as 2×, at ⅛ of the cost;
    - near field r_c = 8 cells and τ cap 1 per cell;
    - an extra cutoff exp(−s/λ_x), with λ_x from the photon budget
      (ΣN_dir = ΣN, applied per distance bin);
    - blended with the split's small scales at σ = 4 cells;
  - status: the ranking (1 cone + budget + blend, 2 split, 3 screened kernel)
    and the implementation order (split → screen → cone) were presented to
    the user on 2026-10-12. No solver change before their go-ahead.
* **User decision (2026-10-12):** implement all three closures (split,
  screened kernel, cone with the budget cutoff and the blend), then run them
  at 512³ like the existing rays, OTVET and M1 runs. The user's remark that
  the split alone does not capture the mean free path after overlap led to
  making the split a modifier on any direct tensor (thin or screened), with
  the cone's blend carrying it for the cone.
* **Implementation choices (this thread):**
  - flags `otvet_closure = thin | screen | cone`, `otvet_split`,
    `otvet_screen_lambda`, `otvet_cone_{nside,coarse,rc,taucap,smax,budget,blend}`
    (`otvet_closure.hpp`); the defaults run the original code path
    unchanged;
  - the split is applied per group and per GLF subcycle with the cell's
    current N_g. The blend's split uses the Γ-weighted N at the start of
    the hydro step, since the Gaussian cannot be redone every subcycle;
  - the split with the cone is refused; the blend replaces it;
  - λ (screen) and λ_x (cone) come from the solution's own Γ-weighted
    photons at the start of each hydro step. The cone keeps its light per
    half-octave distance bin, so the budget needs no second march and no
    state across steps (restart-safe);
  - the screened and cone tensors are isotropic below a noise floor
    (10⁻⁶ of the summed luminosity, normalised, the single-precision FFT's
    floor), and their eigenvalues are clamped to ≥ 0;
  - 512³ runs without `tensor_dump` (14 GB per snapshot); the 256³ runs
    keep it for the a posteriori tensor comparison.
* **User decisions (2026-10-13):** merge PR #35 (done, `fd9d425`); speed up
  the tensor FFTs (branch `otvet-fft-speed`: cached kernel spectra and
  component-parallel transforms, byte-identical); movies of the 512³ runs
  and the after-overlap P_Γ figure as for OTVET and M1 (`docs/otvet_closure.md`
  §7, `$RUNS/otvet_closure_512/compare.mp4`).

## Stage F decisions (2026-10-14): AMR for gas and RT

Taken by the user after Stage E closed (PR #37), on the questions listed
in `docs/STATUS.md` and `docs/rt_multicell.md` §9.

* **Rays across refinement levels: one coarse substep.** No time
  sub-cycling of the transport (neither pure nor the §9 mixture's
  patch-wise sub-cycling): in a level-l patch rays cross 2^l fine cells
  per coarse substep (`rt_cells_per_substep` = 2^l locally), with the
  exposure deposit (`ray_deposit = exposure`), `photon_conserve` and the
  angular resolution set per level (emission/entry nside resolving the
  fine cells, splitting continuing into fine levels). The user's reason:
  the I-fronts move much slower than light in the refined (dense) gas,
  except around quasars, which are not modelled; the uniform-grid study
  passed every target at one level (k = 2) and failed only light-limited
  fronts at 0.1c with two levels (k = 4, `docs/rt_multicell.md` §7–8).
  Per-cell photo-chemistry sub-cycling (no rays) stays available as a
  later option for dense cells with short time-scales.
* **B17 is part of `legacy_physics = 0`.** `photon_conserve` (the
  converged photon limit of `photo_cell`) defaults to `!legacy_physics`;
  an explicit 0 or 1 overrides in either mode. Mirrored in the Fortran
  harnesses (RT harness PATCHES.md 14) so that the C++ stays
  byte-identical to them in every combination of the mirrored flags
  (CLAUDE.md rule 2); the stored `legacy_physics = 0` references were
  regenerated on the CI platform (`.github/workflows/regen_references.yml`).
  Production (`make_production.py`, `legacy_physics = 0`) therefore gets
  it without a flag; `--photon-conserve 0` turns it off. The exposure
  deposit stays a separate flag.
* **Cooling in dense gas: on, with an aggressive star-formation sink.**
  The original's suppression above 200× the mean density (cooling ×
  (200/x)², SN heating × (200/x)^1.5, `rt_cell.hpp`) becomes a parameter
  with the original as its default (`cool_suppress_overdensity = 200`;
  0 = no suppression), and a sink converts gas above an overdensity
  threshold to collisionless star particles (`sf_overdensity`, default 0
  = off; "a few thousand for most runs", so the production default will
  be of that order and set when the sink exists; the efficiency per
  free-fall time is a second parameter). Dense gas keeps cooling and gets
  cold, as the user wants for self-shielded minihalos, and the sink keeps
  the hydro time step and the cooling update sane in collapsed gas. The
  cooling update itself (first order, positive, monotone;
  `docs/port_notes.md` §5b) is kept; sub-cycling it in cells with y ≫ 1
  is a follow-up if the tests show it is needed.
* **The rest is the implementer's choice** (user: "you decide"): hydro
  sub-cycling in time with AMReX `FluxRegister` refluxing; RT deposits per
  level with a coarse/fine rule at the boundary; validation on the Iliev
  tests with refinement and against uniform runs at the matched finest
  resolution, and bitwise against the uniform backend with refinement off.
* **Order (agreed):** (1) gas AMR with refinement off, bitwise identical
  to Stage E; (2) refined hydro on the hydro test suite, with the cooling
  flag and the sink; (3) the RT coarse substep across levels.

## Gravity on the refined levels (2026-10-15)

* **User:** do the cheaper intermediate step before multigrid gravity:
  feed the refined gas levels to the uniform PM mesh and kick every level
  from the mesh, with `pm_refine = 4` as the target (512³ gas, two
  levels, the 2048³ mesh).
* **Implementer's choices:** the parameter `pm_gas_levels`, default 1
  (the level-0 coupling of F1 was documented as interim, and refinement
  has no original behaviour to keep; without refinement the flag does
  nothing, so every existing result is unchanged), 0 for comparison;
  `pm_refine` a power of 2 with it; the particle count apart from the
  mesh (`make_production.py --ndm`), so `pm_refine = 4` runs can use
  1024³ particles rather than 2048³ (whose FoF and P(k) on rank 0 would
  not fit). Multigrid Poisson on the hierarchy stays Stage G.

## Checkpoint cadence (2026-10-18)

* **Decision (user):** future production runs checkpoint every 10 hydro
  steps (`make_production.py --restart-every`, default now 10; it was 1).
  A 512³ refined run's checkpoint is ~240 GB and took ~9 min per write, more
  than the step itself (`fesc512_b2.5_amr`). Runs already set up keep the
  value in their `params.txt`.

## No particle snapshots (2026-10-18)

* **Decision (user):** the dark-matter particles are never needed except for
  restarts, so production snapshots no longer write them. New `data_files =
  nodm` (`norays` without `DM/xv_dm`: the gas files `u, v, w, e` and
  `RT/rt_grid` are kept; the particles stay in the checkpoints), the
  default of `make_production.py` and `amr512_runs.sh`. Saves ~3 GB per
  snapshot at 512³ particles, ~26 GB at 1024³. `w_gas` (ρ_HI, ρ_HeI,
  ρ_HeII, ρ_Z) is also what a refined run's mean free path should use: the
  composite's ρ_HI is the correct mean, unlike rhog × fHI.

## f_esc evolution thread (2026-10-14): a falling escape fraction

* **Request (user):** lower the late-time emissivity with
  f_esc(z) = f6 ((1+z)/7)^β (`f_esc`, `f_esc_beta`, `f_esc_z0 = 6`,
  `f_esc_max`; no code change) so that Γ_HI and the Lyα forest come closer
  to the data. Undershoot ε on purpose and let reionization end later
  (z ≈ 5.3–5.5). Calibrate at 256³ on the box, ICs and parameters of the
  256³ rays run with only the sources changed; then one 512³ rays run to
  z = 5 with the parameters of `reion512_rays`; report how much of λ_mfp's
  excess ε alone removes. Options beyond ε (resolution, a subgrid
  Lyman-limit opacity) are to be discussed before anything is implemented.
  Report: `docs/fesc_evolution.md`.
* **256³ base: `otvcl256`** (20 Mpc/h, the 512³ cell size, case A + RR with
  `recomb_rad_cell_onspot = 1`, f_esc = 0.23). It is the β = 0 point; the
  grid runs differ from it only in the source lines (`tensor_dump`, an
  output-only diagnostic, dropped). Script `tools/hyak/fesc_runs.sh`.
* **Grid (thread choice): f6 set so that z50 stays near 7**, from a one-zone
  Q(z) model fitted to `otvcl256` (C ≈ 4.3): β = 3, 5, 7 as asked
  (f6 = 0.12, 0.075, 0.045) plus β = 2 (f6 = 0.15), added at once because
  the model predicted that β ≥ 4 does not finish by z = 5 at z50 = 7. The
  model underestimates the early photons' efficiency: z50 came out 7.16,
  7.27, 7.63, 8.04. β = 5 and 7 were not rerun at a lower f6: they already
  do not finish reionizing by z = 5 (x_V = 0.98, 0.77), and a lower f6 only
  delays them further.
* **User decision (2026-10-14): the 512³ run uses
  `recomb_rad_cell_onspot = 1`** (the block-emission clump fix, as in the
  256³ base and the production default), although `reion512_rays` was run
  without it; the comparison with `reion512_rays` therefore includes that
  fix (Γ −2 %, history unchanged at 128³; λ_mfp possibly longer).
* **512³ choice (thread): β = 2.5, f6 = 0.135**, between the β = 2 and 3
  points of the grid: A ≈ 1 (the forest) and Γ ≈ Gaikwad+23 by log
  interpolation, and z99 ≈ 5.7 at 256³, i.e. ≈ 5.5 at 512³ (that box ends
  0.18 later at fixed f_esc). A 256³ run at the same values checks the
  interpolation. Snapshots at the closure runs' list
  (15 7 6.5 6.25 6 5.75 5.5 5.25 5); otherwise `reion512_rays`.

## Reduced speed of light by ionisation state (2026-10-17)

* **User decision: production runs use ĉ = 0.1c for x_HII,V < 0.3, 0.2c
  for 0.3 ≤ x < 0.8 and 0.4c for x ≥ 0.8.** Motivation: the runs of
  different threads had used different constant speeds (0.1c, 0.2c,
  0.1c → 0.3c), and the end of reionization and the post-overlap Γ_HI
  depend on ĉ (`docs/production_run.md` §1).
* Implemented as `c_reduced_schedule` (apps/radhydro; the C++ default is no
  schedule, the original) and the `make_production.py` default
  (`--c-schedule`). Choices of the implementing thread: x is the
  volume-weighted x_HII of the box at the end of the previous hydro step;
  x = 0 until the RT has written its grid; no hysteresis (stateless, so
  restarts are exact without a checkpoint change).

## Reduced speed of light, revised schedule (2026-10-18)

* **User decision: ĉ = 0.1c for x_HII,V < 0.25, 0.2c to 0.5, 0.3c to 0.75,
  0.5c above** (`c_reduced_schedule = 0.1 0.25 0.2 0.5 0.3 0.75 0.5`, the
  `make_production.py` default; replaces 0.1c / 0.2c from 0.3 / 0.4c from
  0.8). Evidence discussed: the 20 Mpc/h 256³ constant-ĉ runs (`runs/sol`)
  delay x_HII,V = 0.99 by Δz = 0.09 (0.5c), 0.19 (0.3c), 0.30 (0.2c),
  0.56 (0.1c) against 1c, about half of the old schedule's delay coming
  from the 0.1c phase; and the post-overlap Γ_HI of the fesc512 pair
  (0.1c) at a third of the local-equilibrium value σ ṅ λ, the background's
  relaxation λ/ĉ (200–300 Myr at 0.1c) being longer than the time over
  which λ_mfp doubles (~70 Myr at z = 5.5–5.2): an intermediate step and a
  higher maximum. The RT substeps per hydro step scale with ĉ (the late
  phase at 0.5c costs ~5× the RT time of 0.1c).

## Rays at the refinement boundaries (2026-10-20/21)

* **Decisions (user):** `amrex_split_on_entry = 1` (rays split on arrival
  in a finer level, no radius clamp) and `amrex_ghost_rule = conserve` (a
  fine ray's deposit outside its level goes to the parent with weight
  1/(1 + √fgam): the parent absorbs what the ray loses) are the defaults.
  Evidence: `docs/stage_f.md` §3, "Rays at the level faces".

## Zoom-in turbulence programme (2026-10-22)

* **Request (user):** a plan to zoom repeatedly on one region of the box
  and resolve the turbulence an I-front drives as it passes, the physics
  of Cain, McQuinn, Scannapieco, D'Aloisio & Trac 2026 (PRD 113, 023050;
  uniform 0.25 h⁻¹ Mpc boxes with plane-parallel fronts) inside a full
  reionization run. The ICs are seeded at the region's resolution; the
  region holds at least 1024³ cells; klone builds it up and a proposal
  for allocation time funds the suite; magnetic fields come later.
* **Plan:** `docs/zoom_turbulence.md` (work packages WP1–WP8, the
  ladder Z1–Z4, the validation rungs T0–T5, cost estimates, the proposal
  skeleton). No decision taken yet: the open choices are its §9 (base
  box, rungs, regions, IC tool and cosmology, the RT's HEALPix cap or
  beams and the coarse-substep rule in the zoom, preheating variants,
  MHD scope, proposal venue).
* **Decisions (user, 2026-10-22):** base 40 h⁻¹ Mpc / 512³; the rungs
  Z1–Z4 of `docs/zoom_turbulence.md` §3; one mean-density region first,
  then a void and a filament patch; 1LPT ICs with Planck cosmology (so
  the cosmology becomes a runtime parameter set, the original as the
  default); raise the HEALPix cap first, beams only if the ray counts
  demand it; no preheating in the first runs; MHD later, with the
  Biermann battery and one more non-ideal term at the I-fronts; the
  proposal to NSF ACCESS as Discover (< ~10⁶ core-hours, a paragraph,
  quick) then Accelerate (< ~3 × 10⁶, a few pages). On the IC tool the
  user asked why extending `make_ics.py` was recommended over MUSIC; the
  answer (`docs/zoom_turbulence.md` §4.1) proposes MUSIC with a
  converter, to be confirmed.
