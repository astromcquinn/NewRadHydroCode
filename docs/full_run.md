# The full radhydro loop (Phase 1c)

`apps/radhydro` runs the complete time loop of `reference/Radhydro/radhydro.f90`:
hydro (`hydro_tools.f90 -DMF`), PM gravity (`pm_tools.f90`), sources (legacy
halo catalogs or the on-the-fly friends-of-friends finder) and ray-tracing
RT. `reference/harness_full` is a gfortran build of the original program
(the unmodified modules plus the harness patches) that serves as its
reference. **The two are bit-identical** on every validation case, with
legacy catalogs and with on-the-fly sources, with the original code
(`legacy_rt = legacy_physics = 1`) and with all fixes on
(`legacy_rt = legacy_physics = 0`).

Contents: 1. unification of the sub-ports; 2. the driver; 3. the Fortran
harness; 4. validation problems and results; 5. demo run and cost model;
6. quirks; 7. open questions.

## 1. One implementation of the shared infrastructure

The hydro and PM ports were made in parallel and duplicated some pieces.
Each now exists once, in C++ and in the Fortran harnesses:

| piece | before | now (C++) | now (Fortran harnesses) |
|---|---|---|---|
| FFT | `Fft1d` in hydro.cpp (1-D, aligned buffer) and `fft.cpp` (3-D) | `include/radray/fft.hpp`: `fft::Fft1d` (MKL 1-D descriptor) and `fft::fft_3d`; `fft_1d` uses `Fft1d` | `reference/harness_common/mkl_dfti_shim.f90`: one `MKL_DFTI` module (1-D and 3-D) under the original `hydro_tools.f90` and `fft_tools.f90` (the PM harness no longer patches `fft_tools.f90`) |
| time step | `calc_timestep` (RT, fixed da), `calc_timestep_hydro`, `pm::calc_timestep_full` | `radray::calc_timestep(cfg, units, TimestepMode::{Original,Courant,Fixed}, da, clk, out)` in `cosmology.hpp` (the RT driver's fixed-da form is a thin wrapper) | `harness_common/harness_timestep.f90` (`calc_timestep_mode`, `init_variables_z`) for the hydro and PM harnesses; the full harness uses the original routine of `radhydro.f90` |
| `summarize_gas` | full routine (hydro) + `summarize_gas_speeds` (PM) | the full `summarize_gas` everywhere | the original `summarize_gas` everywhere |
| gas domains | `Domains` + `init_domain_gas` (hydro), `GasDomains` (PM) | `Domains`/`init_domain_gas` (`domain_gas.cpp`), passed to `init_pm` | the original `init_domain_gas` |
| speed globals | `HydroState::cmax_gas…`, `PMState::amax_dm…`, `GasSpeeds` | `SimClock` fields (the "Sim" variables of `radhydro_global.f90`) | the original globals |
| Units | no particle units | `Units::mp_unit, Msun_unit` (calc_sim_var), used by the halo finder | original |
| sizes | `box` compile-time, `Ndm = Ngas/2` | `Config::box`, `Config::Ndm`, `Npart`, `x_dm_box` at run time | `-DBOX`, `-DNDM` |
| gfortran MAX/MIN | local to hydro.cpp | `f::gmax`, `f::gmin` (fortran.hpp) | — |

Every pre-existing CTest stayed bit-identical through the unification, and
the stored outputs of `reference/harness_hydro` and `reference/harness_pm`
are reproduced bit for bit by the unified harnesses (-O2 and
-O0 -fcheck=all). Libraries: `radray_fft` (fft.cpp) ← `radray_hydro`
(hydro, domains, summarize_gas) ← `radray_pm`.

## 2. apps/radhydro

```
radhydro <problem_dir> <outdir> [--set key=value ...] [--quiet]
```

The sequence is `radhydro.f90`'s, statement for statement:

```
init_variables(z_init); init_domain_gas; init_domain_rt; init_hydro; init_pm;
init_rt (rates, legacy_physics); init_ray_tracing; init_io (zdata1/2); init_fft;
calc_sim_var; read_init_conditions; calc_sim_var; summarize_dm/gas/rt
do
  calc_sim_var; calc_timestep (original: istep 0 -> da/a = 1/9; then
     min(Courant, PM), ra_max); data1/data2 flags
  transpose_ij, transpose_ik, hydro forward, transpose_ik, hydro forward,
  transpose_ij, hydro forward; half clock step; calc_sim_var
  particle_mesh
  hydro reverse, transpose_ij, hydro reverse, transpose_ik, hydro reverse,
  transpose_ik, transpose_ij; half clock step; calc_sim_var
  if (idata1 > 1): read_halo | on-the-fly FoF; if (Nsource > 0) calc_halo,
     radiative_transfer
  summarize_dm, summarize_gas, summarize_rt
  write_data_2; calc_power_spectrum + write_data_1
  exit if redshift <= zdata1(Ndata1)
```

### Parameters (`params.txt`, `--set key=value`)

| key | default | meaning |
|---|---|---|
| `N` | 32 | Ngas = Nrt (multiple of 8) |
| `Ndm` | N/2 | DM particles per dimension (N multiple of Ndm) |
| `box` | 5 | Mpc/h |
| `nray_max_factor` | 64 | Nray_max = factor·N³ (see §6, quirk F3) |
| `z_init` | 300 | initial redshift (the original's `zinit`) |
| `zdata1`, `zdata2` | files `zdata1.txt`, `zdata2.txt` of the problem | output redshifts (list); RT starts after the first zdata1 output, the run ends at the last |
| `ics` | problem dir | directory of `xv_dm.init`, `u_gas.init` (`tools/ics/make_ics.py`, now with `--box`) |
| `legacy_rt`, `legacy_physics` | 1, 1 | the RT-bookkeeping and atomic-physics fixes (0 = fixed; recommended) |
| `source_mode` | onthefly | `legacy` (read_halo: `dir_halo`/sfr.txt + halo.z=XX.XXXX) or `onthefly` |
| `dir_halo` | `<problem>/Halo` | legacy catalogs |
| `fof_b`, `fof_min_npart`, `fof_so` | 0.2, 20, 0 | finder (docs/halos_and_sources.md) |
| `sfr_model`, `tvir_min`, `m_min`, `l_per_sfr`, `sfr_history` | global, 1e4, 0, 8e27, — | source model |
| `write_catalogs` | — | write the on-the-fly catalogs in legacy format |
| `add_metals` | 1 | the add_halo_to_rt metal injection |
| `fields_every`, `fields_at_data1` | 0, 1 | field output schedule |
| `max_steps` | 1000 | the original's `timing(1000)` |

`N`, `Ndm`, `box` and `nray_max_factor` are compile-time parameters in the
original; the C++ takes them at run time and the Fortran harness checks that
the problem matches its build.

### Outputs

* `step_%04d/stats.txt` every step (step 0 = initial state): `istep, idata1,
  idata2, nstep_rt, Nhalo, Nsource, nray_active, ncst1, nmrg1, nspl1, ndel1,
  redshift, scalefactor, year, da_sim, dt_sim, dt_yr, cmax_gas, vmax_gas,
  vrms_gas, c_limit, amax_dm, arms_dm, vmax_dm, vrms_dm, ngam_rays_total,
  Luv, sfr, xHII_V, xHII_M, xHeIII_V, T_V, T_M` (the speeds in grid units,
  as after the summaries; the last five are sums over the RT grid);
* `step_%04d/<field>.bin` at the zdata1 outputs and the last step:
  `fHI fHeI fHeII temp rhog ngam nu fgam tau gamma lambda dx` (rt_grid),
  `u_gas w_gas v_gas xv_dm` (float32, Fortran order);
* `Log/hydro.log` (the summarize_gas unit-100 records, byte-identical to the
  Fortran), `Log/ray.log`, `Pk/pk.<ztag>`, `DM/xv_dm.<ztag>`,
  `timing.txt` (per step: hydro, PM, halos, RT, summaries, IO seconds).

## 3. reference/harness_full

See `reference/harness_full/README.md` and `PATCHES.md`. In short: the
original `radhydro.f90` with HARNESS_FULL patches for the command line,
`z_init`, the flags, the source mode and the compact output; `box` by the
preprocessor; everything else from the existing harnesses.

**On-the-fly sources in Fortran (task 2b): done.** `src/halo_capi.cpp`
exports `radray_onthefly_sources`, a C function that runs the same
`onthefly_catalog` (find_halos + assign_sources) the C++ driver calls. The
Fortran harness calls it through `bind(C)` from `read_halo_onthefly` in
place of `read_halo` and links the C++ libraries. Both codes therefore use
identical sources whenever their particles are identical, and the on-the-fly
mode is validated Fortran-vs-C++ like the legacy mode.

**The `idata1 > 1` gate.** RT (and the source search) only starts after the
first zdata1 output: `idata1` is incremented by `write_data_1`. The
validation problems put the first zdata1 redshift where halos already exist.

## 4. Validation

### Problems (`tests/problems_full`, made by `make_problems.py`)

| problem | grid | sources | z range | steps | RT |
|---|---|---|---|---|---|
| `legacy` | Ngas 32, Ndm 16, box 4 Mpc/h | legacy catalogs at z = 10, 9, 8, 7, 6 made by the C++ FoF finder and legacy writer from a DM run of the same ICs; `sfr.txt` = Hernquist & Springel fit | 100 → 7.06 (zdata1 = 8.5, 7.5) | 24 | last step, 1669 substeps, 1 source |
| `onthefly` | same | C++ FoF (b = 0.2, ≥ 8 particles) every RT step | same | 24 | same |
| `smoke` | Ngas 16, Ndm 8, box 4 | one hand-made legacy source (z = 40 and 20 catalogs) | 100 → 19.8 (zdata1 = 30, 20) | 15 | 3 steps, 1479 substeps |

Why so small and so short: with the original time step (da/a ≤ 1/9 until
the Courant limit takes over, which never happens at this resolution) every
hydro step spans ~1500–3300 light-crossing RT substeps, and the RT cost is
∝ rays × substeps (§5). One RT hydro step at N = 32 takes 4–6 minutes, so
the main problems do RT in their last step only. Ndm = 16 resolves halos
only above ~1e10 Msun, so the box is 4 Mpc/h and the FoF minimum is 8
particles; the first halo appears at z ≈ 8 (seed 7, chosen from a scan of 8
seeds for an early halo).

### Results

`compare.py --rtol 0 --atol-frac 0` and `tests/compare_full_reference.py`
(stats of every step, SHA-256 of every field file, `hydro.log`, P(k) and
particle files byte for byte):

| case | legacy_rt / legacy_physics | Fortran vs C++ |
|---|---|---|
| smoke (legacy) | 1/1 and 0/0 | identical, all 15 steps |
| smoke (on the fly, no halos form) | 1/1 and 0/0 | identical |
| legacy | 1/1 and 0/0 | identical, all 24 steps |
| onthefly | 1/1 and 0/0 | identical, all 24 steps |

Also identical: the Fortran `-O0 -fcheck=all` build (smoke both modes,
onthefly 1/1; no runtime errors), and the ASan/UBSan C++ build (smoke, both
modes; clean). A 26-step hydro+PM-only run of the N = 32, box 2 ICs was
identical before any RT was added (hydro.log, P(k), particles).

CTests (label `full`):

| test | what | time |
|---|---|---|
| `full_ref_smoke_rt1/rt0` | C++ vs stored Fortran outputs | 2 s |
| `full_ref_{legacy,onthefly}_{rt1,rt0}` | same, label `long` | 5–6 min each |
| `full_vs_fortran_live_smoke_{legacy,onthefly}_{rt1,rt0}` | fresh gfortran build + run vs C++ | 2–7 s each (after the first build) |
| `full_sanitize_smoke_rt1/rt0` | ASan/UBSan radhydro vs stored outputs | ~10 s |

Stored references: `reference/harness_full/outputs/<case>` (stats of all
steps, field hashes of the output steps, the full last-step fields for
`legacy_rt1` and the smoke cases; 9 MB in total).

## 5. Demo run (C++ only) and cost model

### Set-up

Ngas = 64, Ndm = 32, box = 4 Mpc/h, Zel'dovich ICs at z = 100 (seed 7),
on-the-fly sources (FoF b = 0.2, ≥ 10 particles = 2.1e9 Msun; global SFR
normalised to the Hernquist & Springel fit), `legacy_rt = legacy_physics = 0`,
`nray_max_factor = 4` (Nray_max = 1.05M), zdata1 = 8.5, 8, …, 5.5. One core;
during part of the run three other jobs shared the machine (4 cores), so the
RT time is an upper bound by maybe 10–20 %.

**It was stopped after the first RT step**, as agreed (> 30 min): the
pre-RT part (23 steps, z = 100 → 7.95) took 33 s, the first RT step
(z = 7.95 → 7.06, 3338 RT substeps) 4945 s = 82 min.

### Results

* **Ionized fraction.** RT starts after the zdata1 output at z = 7.95, when 5
  halos exist. Over the next hydro step (z = 7.95 → 7.06, 17 sources at its
  end, comoving SFR density 0.12 Msun/yr/Mpc³) the box is completely
  ionized: x_HII = 0.999995 (volume), 0.999992 (mass), x_HeIII = 0.038,
  mean T = 1.7e4 K. A box of 4 Mpc/h holds ~1e69 H atoms; the HS03 SFR
  density gives ~2 ionizing photons per atom in one ~110 Myr hydro step
  (f_esc = 0.2), so a single step suffices. With the original time step
  (da/a = 1/9, §6 F4) the reionization history is therefore not resolved by
  the hydro-step outputs (the RT substeps resolve it internally). A
  meaningful x_HII(z) needs the abundance-matched source model and/or
  outputs inside the RT step.
* **Halo mass function** (`docs/figures/full64_hmf.png`): 17 FoF halos at
  z = 7.06 and 26 at z = 6.25, 2e9–1.2e10 Msun. Against Sheth–Tormen the
  counts are 0.3–0.5 of ST in the lowest bin (10–15 particles) and
  0.25–0.85 above (log M/(Msun/h) ≥ 9.4, 1–11 halos per bin, so Poisson
  errors of 30–100 %). The deficit, largest at the low-mass end, is what PM
  force resolution (one gas cell = 62 kpc/h) and 10–40-particle FoF groups
  give; ST is not expected to hold there.
* **Power spectrum** (`full64_pk.png`): the DM Δ²(k) (CIC-deconvolved) is
  0.72 and 0.93 of linear theory in the two lowest bins (k = 1.6, 3.1 h/Mpc;
  4 and 13 modes, cosmic variance ~50 % and ~28 %), 1.3–1.7 at
  k = 4.7–6.3 h/Mpc (non-linear growth at z ≈ 7), and shot noise
  (V/Npart = 0.002 (Mpc/h)³) dominates above k ≈ 15 h/Mpc. The gas follows
  the DM at low k and is smoothed by pressure at high k.
* **T–ρ** (`full64_trho.png`, z = 7.06, just after reionization): T0 =
  1.7e4 K at the mean density, γ − 1 = −0.08 (fit over |log δ| < 0.5): a
  nearly isothermal, slightly inverted relation, as expected right after
  photo-heating (denser gas cools faster); a tail to 3e4–5e4 K at δ > 10
  from shock heating.
* **Timing** (`full64_timing.png`): see the table below.

### Cost per component (one core, -O2)

| N (Ngas) | Ndm | hydro / step | PM / step | summaries / step | FoF / call | memory without rays |
|---|---|---|---|---|---|---|
| 32 | 16 | 0.19 s | 0.007 s | 0.009 s | ~5 ms | ~20 MB |
| 64 | 32 | 1.27 s | 0.048 s | 0.048 s | 34 ms | 99 MB |
| 128 | 64 | 9.3 s | 0.36 s | 0.40 s | — | 491 MB |

Hydro is ~0.6 µs per cell per sweep (6 sweeps per step) and scales as N³;
PM is 4 % of the hydro; the halo finder is negligible. **RT is everything
else**:

| run | rays at the end of the RT step | substeps | RT time | per substep | peak RSS |
|---|---|---|---|---|---|
| N = 32, nray_max_factor 4 (Nray_max 131k) | 66k | 1669 | 290–340 s | 0.2 s | 26 MB |
| N = 32, factor 16 (524k) | 207k | 1669 | 585 s | 0.35 s | 48 MB |
| N = 32, factor 64 (2.1M, the original) | 675k | 1669 | 1517 s | 0.91 s | 120 MB |
| N = 32, factor 1 (32k) | — | — | fails: the ray pool overflows (§6 F3) | — | — |
| N = 64, factor 4 (1.05M) | 629k (835k peak) | 3338 | 4945 s | 1.48 s | 168 MB |

* **Cost model.** RT time ≈ (number of RT substeps) × (≈1.5–3 µs per ray
  per substep) × (rays). The ray population grows during the step and
  saturates at roughly 0.3–0.6 × Nray_max with `legacy_rt = 0` (merging
  thresholds are fractions of Nray_max), so cost ∝ nray_max_factor·N³ ×
  substeps. Rays move one cell per substep, so the number of substeps over
  an interval is the comoving light travel distance in cells:
  S = 2c/(H0√Ωm)·(√a_end − √a_start)/Δx, Δx = box/(h N). For box = 4
  Mpc/h, z = 12 → 6: S ≈ 4.7e3·(N/16), i.e. 9.3k substeps at N = 32,
  18.6k at N = 64, 37k at N = 128.
* **Memory.** ≈ 20 MB + 230 B per cell (hydro, RT grid, PM, ray heads;
  99 MB at N = 64, 491 MB at N = 128) plus ≈ 140 B per ray actually used
  (ray, index, free lists). Nray_max only
  reserves address space; the pages are touched as rays are created.
* **Extrapolation to a run through reionization (z = 12 → 6, box 4 Mpc/h,
  factor 4):** N = 32: ~9k substeps × 0.2 s ≈ 0.5 h; N = 64: ~19k × 1.5 s
  ≈ 8 h; N = 128: ~37k × ~12 s (8× rays) ≈ 5 days, and ~1 GB of rays. Hydro
  and PM are < 1 % of that at every N. A production run at N ≥ 64 on this
  machine therefore needs one of: a reduced speed of light (the original
  has `cscale = 1` hard-wired in `calc_timestep`), fewer rays
  (nray_max_factor ≈ 2–4, but see F3), or the Phase 2 parallelisation.

## 6. Quirks found in this phase

* **F1. `write_data_1/2` with `Ncpu = 1` write only the particles** (the
  `select case (icpu)` loop, like `read_init_conditions`); not patched, no
  effect on the state (reference/harness_full/PATCHES.md).
* **F2. `Nhalo_max = Npart/1024`** is 4 for Ndm = 16 and `read_halo` has no
  bound check; the harness uses `NHALO_MAX = 16384`.
* **F3. The ray pool can overflow.** With `nray_max_factor = 1` the rays of
  one bright source exceed Nray_max = N³ within one RT step, in both
  `legacy_rt` modes: merging is triggered by thresholds on the photon
  content, not by the pool size. The C++ throws ("ray storage index out of
  range"); the original would write past `ray(:)`. Factor 4 was enough in
  every run here; the original default is 64.
* **F4. The Courant limit never binds at these resolutions.** Every step is
  the ra_max step da/a = 1/9 (z = 100 → 7 in 24 steps). An RT step then
  spans ~1500–3300 light-crossing substeps (N = 32–64, box 4) and
  ~100 Myr, which is why the RT validation is short and why reionization
  completes inside one hydro step in the demo.
* **F5. Before the first RT call the rt_grid is all zeros** (`fHI = 0`,
  `temp = 0`), so the stats' x_HII is 1 and T is 0 until RT starts
  (diagnostic only; both codes).
* **F6. `fft_tools.f90` uses `c_loc(a)` on a dummy without `target`**:
  ifort accepts it, gfortran does not compile it (harness_common fix).
* **F7. `init_data` does not initialise `data_2d%luv`** (D2D slices only).
* Pre-existing quirks that matter for the full loop: Q1 (one-thread
  `read_init_conditions`, patched as before), H6/Q13 (`summarize_gas`
  clamps the state), Q3 (speed globals converted to km/s in place by
  `calc_timestep`).

## 7. Open questions

1. **Production size.** Is a reduced speed of light (a `cscale < 1` in
   `calc_timestep`, behind a flag) acceptable for the production run, or
   should the RT substep stay light-crossing and the run wait for Phase 2?
2. **nray_max_factor.** 4 is ~4.5× cheaper than the original 64 at N = 32,
   and the global state after the RT step agrees to 1e-5 (x_HII
   0.999991 vs 0.999997, mass-weighted T within 0.1 %); should the
   default change, and should the pool grow (or force merging) instead of
   overflowing (F3)?
3. **Sources.** The global HS03 normalisation puts all star formation in the
   few resolved halos and ionizes a 4 Mpc/h box within one step; the
   abundance-matched model planned next should replace it for the
   reionization run.
4. **F1** (`write_data_*` with one thread): patch it like
   `read_init_conditions`, so restarts work with Ncpu = 1?
5. The four long CTests (~5 min each) could be moved to a nightly label if
   that is too slow for routine use.
