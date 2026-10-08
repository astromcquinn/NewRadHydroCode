# Project status

*Last updated: 2026-10-08, end of the first working thread.*

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

## In progress (when this file was written)

These are background agents in the first thread. Their work lands in a
follow-up PR.

1. **Production readiness:**
   - `c_reduced` (reduced speed of light via the original's `cscale`);
   - deterministic OpenMP (thread-count independent), with speedups;
   - the ray pool forcing merges instead of overflowing (legacy_rt=0);
   - the F1 restart fix;
   - long tests moved to a `nightly` label;
   - a timing/memory table for the production config.
2. **Hydro test suite** (`tests/hydro_suite`, `docs/hydro_tests.md`):
   - Sod (including strong-shock and 123 cases) with convergence;
   - Sedov–Taylor;
   - Zel'dovich pancake with gas (Bryan et al. 1995);
   - adiabatic cosmological box;
   - sound waves;
   - Kelvin–Helmholtz if feasible.

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
5. After the production run, follow `docs/ROADMAP.md`:
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

- **Serial code.** OpenMP is in progress. MPI/GPU is Phase 2.
- **RT dominates the cost.** Hydro+PM take under 1%. At full c, z=12→6 at
  128³ would take days on one core, hence `c_reduced`.
- **Network access from cloud sessions.** arXiv, ADS and journal sites are
  unreachable.
