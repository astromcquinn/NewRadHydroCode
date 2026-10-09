# Project status

*Last updated: 2026-10-09 (after PR #10 and the ray compression study).*

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

## Measured production costs (PR #10, `docs/production_run.md`)

Box 10 Mpc/h, Ngas = Ndm, factor 16, z = 100 → 5.5, fixed flags,
`deterministic_parallel`:

| N | c_reduced | Measured / projected time | Memory |
|---|---|---|---|
| 128 | 0.01 | 69 min measured on about 2 cores; about 30 min projected on 4 | 2.8 GB |
| 192 | 0.01 | about 2.5–3 h on 4 cores | 8–11 GB; checkpoints of 6–9 GB each |

At 0.01c the light horizon is only about 16 Mpc/h over z = 15 → 5.5.

## Next

1. Check the halo mass function at production resolution against Sheth–Tormen/Tinker.
2. **Abundance-matched sources:**
   - Fit L_UV(M_halo) to the observed high-z UV LFs (HST z≈6–10, JWST z>10).
   - Calibrate f_esc to Planck τ_es.
   - No sub-grid sources for unresolved halos (user decision).
3. **Production reionisation run and movie.** The box is decided; the
   resolution waits on the timing table.

   | Setting | Value |
   |---|---|
   | Box | 10 Mpc/h |
   | c_reduced | 0.01 (the user accepts slow light) |
   | nray_max_factor | 16 |
   | Ngas = Ndm | 128³, or 192³ if the timing table allows |
   | Flags | legacy_rt = legacy_physics = 0 |
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
6. **Hydro fix for cold supersonic gas** (see the hydro-suite finding above),
   behind a flag, before the plane-parallel science runs.
7. **Plane-parallel / RT-domain mode** (`docs/plane_parallel.md`).
8. After the production run, follow `docs/ROADMAP.md`:
   - ray compression study: done (`docs/ray_compression.md`); 64³
     production-configuration check done, 128³ job script ready for Hyak;
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
