# Project status

*Last updated: 2026-10-09, production-run thread (branch `production-run`, PR #7).*

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

## Production run, sources, calibration, movie (PR #7, branch `production-run`)

**Abundance-matched sources** (`docs/abundance_matching.md`):
- `sfr_model = abundance`, a new source model; the existing models are
  unchanged and bit-identical.
- L_UV(M_halo, z) comes from abundance matching to observed UV LFs with
  the Tinker 2008 HMF. The shape is fitted at z ≤ 9, with the faint end
  extrapolated to M_UV = −12. At z = 10–12.5 only the normalisation is
  fitted, to the JWST LFs of Donnan+24. L(M) is held fixed beyond the
  data.
- Ionizing rate = f_esc ξ_ion L_UV, with log ξ_ion = 25.38. The RT code
  is untouched.
- Tools: `apps/am_fit`, `apps/hmf_table`, `tools/analysis/plot_am.py`,
  `tools/analysis/box_lf.py`.
- The LFs at z ≥ 9 are read from the paper (Donnan+24). The z ≈ 6–8 rows
  (Bouwens+21) are still unconfirmed; the user will supply JWST-era LFs.

**Halo mass function** at 128³ / 10 Mpc/h (`docs/hmf_validation.md`):
- FoF agrees with Sheth–Tormen above about 100 particles; it is 1.3–2.2×
  above Tinker.
- There is a 2–4× deficit at 20–60 particles, growing with z (PM force
  resolution and Zel'dovich ICs). The first halo forms at z ≈ 11.3.
- So the box's UV LF is 1.3–2× above the observed one brighter than
  M_UV ≈ −16, and below it at the faint end.

**f_esc calibration** at 128³ (`docs/reionization_production.md` §1):
- **f_esc = 0.50 gives τ_es = 0.0540** (Planck 2018: 0.054 ± 0.007), with
  the midpoint at z = 7.0. x_M = 0.985 at z ≈ 5.3, so the end is late.
- Each run takes 13–16 min on 4 threads.
- τ is a weak constraint: f_esc 0.2–0.8 gives τ = 0.046–0.058. Only
  f_esc ≳ 0.8 ends reionization by z ≈ 5.7.
- The biases are documented: 0.01c, no sub-grid sources and late first
  halos bias f_esc high; the box's excess of bright halos biases it low.

**Movie** (`tools/movie/make_movie.py`):
- Each frame has slices of x_HII, T and gas density, plus x_HII(z) and
  τ_es panels; ffmpeg assembles the MP4.
- Test movies: 64³ (`runs/`, not committed) and 128³
  (`docs/figures/reion_test_128.mp4`).
- Movie frames written before the first RT substep have an uninitialised
  rt_grid. The tools skip them.

**Hyak** (`tools/hyak/`, `docs/reionization_production.md` §3):
- Scripts: build, thread-scaling test, and a production job that resumes
  automatically from the latest checkpoint.
- Cost estimate: 128³ takes minutes on a node; 256³ about 0.5–1 h on
  32–64 cores (about 20 GB); 512³ about 8–16 h (about 160 GB). Thread
  scaling on Genoa is still to be measured.

**Update (user decisions, 2026-10-09):**
- L(M) is now matched to the simulated FoF mass function
  (`data/hmf/fof_128_10mpc.txt`, `am_fit --hmf st --sim-hmf ...`).
- The production run starts at f_esc = 0.5.
- **All LF parameters are now read from the papers** the user supplied:
  Bouwens+22 (z ≈ 6, 8), GLIMPSE/Atek+26 (z ≈ 7), Donnan+24 (z = 9–12.5),
  with cross-checks. The Bouwens+21 φ* values recalled from memory had
  been wrong.
- With these LFs and f_esc = 0.5: τ_es = 0.0516, midpoint z = 6.72,
  x_M = 0.967 at z ≈ 5.3.
- The fit is in `data/uv_lf/am_params_128_10mpc.txt`.
- ξ_ion = 10^25.38 is confirmed from Pahl et al. 2025.

**Not run yet:** the production run itself (the commands are in the doc),
and the scaling test on Hyak.

## Measured production costs (PR #10, `docs/production_run.md`)

Box 10 Mpc/h, Ngas = Ndm, factor 16, z = 100 → 5.5, fixed flags,
`deterministic_parallel`:

| N | c_reduced | Measured / projected time | Memory |
|---|---|---|---|
| 128 | 0.01 | 69 min measured on about 2 cores; about 30 min projected on 4 | 2.8 GB |
| 192 | 0.01 | about 2.5–3 h on 4 cores | 8–11 GB; checkpoints of 6–9 GB each |

At 0.01c the light horizon is only about 16 Mpc/h over z = 15 → 5.5.

## Next

1. **User decisions pending for the production run**
   (`docs/reionization_production.md`):
      - the production N on Hyak (128 or 256);
   - a 0.1c reference run.
2. **Production run on Hyak** with `tools/hyak/reion_production.slurm`,
   then the movie (`make_movie.py`). First run the scaling test.
3. Recalibrate f_esc when the LFs or the matching change: four runs of
   about 15 min at 128³ (`calibrate_fesc.py`).
4. **User's cluster: UW Hyak `klone`, account `astro`.**
   - Owned: about 1000 `cpu-g2` cores (about 5 nodes of AMD Genoa, 192 cores
     and 1.5 TB per node) with 7.8 TB of memory, a 64-core,
     1 TB `cpu-g2-mem2x` slice, and 40 `compute-bigmem` cores.
   - No owned GPUs: GPUs only via preemptible `ckpt` partitions. No time
     limits. See `docs/ROADMAP.md` (answers section).
   - Access: Claude Code installed in user space, `claude remote-control`
     inside tmux, batch jobs with restarts.
5. **Hydro fix for cold supersonic gas** (see the hydro-suite finding above),
   behind a flag, before the plane-parallel science runs.
6. **Plane-parallel / RT-domain mode** (`docs/plane_parallel.md`).
7. After the production run, follow `docs/ROADMAP.md`:
   - ray compression study, which can start now;
   - framework spike (AMReX vs Parthenon);
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

- **OpenMP only.** MPI/GPU is Phase 2. Thread scaling on Genoa is still to be measured on Hyak.
- **RT dominates the cost.** Hydro+PM take under 1%. At full c, z=12→6 at
  128³ would take days on one core, hence `c_reduced`.
- **Network access from cloud sessions.** arXiv, ADS and journal sites are
  unreachable.
