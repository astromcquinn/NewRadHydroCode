# Fortran reference harness (RT only)

This directory builds Hy Trac's original RadHydro radiative transfer with
gfortran. It runs one thread, with no hydro and no gravity, and implements
`docs/validation_spec.md`. Its outputs are the reference that the C++ port is
compared against with `tools/compare/compare.py`.

```
reference/harness/
  src/                 reference modules (minimally patched, see PATCHES.md) + driver
    rt_harness.f90     the driver (replaces radhydro.f90)
    ifport_shim.f90    module IFPORT: time(), qsort()
  healpix/             HEALPix 3.11 Fortran subset (pix_tools and its dependencies)
  Makefile
  make_test_problem.py fallback generator for tests/problems/*
  run_problem.sh       run one problem (sets the stack limit)
  run_all.sh           build what is needed and run every problem
  tests/check_healpix.py, hpx_check.f90   HEALPix-vs-healpy check
  tests/analyze_outputs.py                sanity checks (I-front, photon budget)
  outputs/             stored reference outputs (subset, see below), legacy_rt = 1
  outputs_legacy_rt0/  the same with legacy_rt = 0 (RT bookkeeping fixes on)
  outputs_legacy_rt0_phys0/  legacy_rt = 0 and legacy_physics = 0 (also the
                       atomic-physics fixes and B17; the recommended configuration)
  outputs_c0.1_fixed/  the same two flags with c_reduced = 0.1 (two problems)
  store_outputs.py     copy the stored subset of a run_all.sh output root
```

## Build

You need gfortran (tested with 13.3) and gcc. Nothing else is required: no
MKL, cfitsio or external HEALPix library.

```
make                    # build/N32_R64/rt_harness.x  (-O2 -g)
make N=16               # build/N16_R64/...
make N=64               # build/N64_R64/...
make NRF=4              # Nray_max = 4*N^3: build/N32_R4/...
make debug              # build/debug_N32_R64/ (-O0 -g -fcheck=all -fbacktrace -Wall)
make debug FPE=overflow # also trap FP exceptions (see the ne/ne_prev caveat below)
make test-healpix       # compare pix2vec_nest/vec2pix_nest with healpy (pip install healpy)
```

`N` (grid size, `Ngas = Nrt = N`) and `NRF` (`nray_max_factor`) are
compile-time parameters, because they size static arrays in the original code.
The executable checks that they match the problem's `params.txt` and stops
with a "rebuild with make N=.. NRF=.." message if not.

## Run

```
./run_problem.sh <problem_dir> <outdir> [nsteps] [exe] [--set key=value ...]
./run_all.sh [outroot=reference/harness/runs] [opt|debug]
RT_SET=legacy_rt=0 ./run_all.sh <outroot>      # RT bookkeeping fixes on
RT_SET="legacy_rt=0 legacy_physics=0" ./run_all.sh <outroot>   # + physics fixes, B17
python3 store_outputs.py <outroot> outputs_legacy_rt0_phys0      # keep the stored subset
```

`photon_conserve` (B17, PATCHES.md 14) follows `legacy_physics` unless set:
`--set photon_conserve=0` with `legacy_physics=0` gives the fixed physics
with the original photon limit. The stored sets are regenerated on the CI
platform by `.github/workflows/regen_references.yml` (docs/cluster.md
section 3: klone's libm differs in the last place on some problems).

The executable takes the same `--set key=value` arguments:
`rt_harness.x <problem_dir> <outdir> [nsteps] [--set key=value ...]`. Each one
overrides or adds a `params.txt` entry. The C++ `rt_driver` accepts the same
`--set`.

`run_all.sh` does the following:

1. Calls `make_test_problem.py`, which only writes problems that do not
   exist yet.
2. For each problem it reads `N` and `nray_max_factor`, builds the matching
   executable, runs it, and deletes the `D2D/` slice dumps.

Running the executable directly (`build/.../rt_harness.x <problem_dir>
<outdir> [nsteps]`) needs `ulimit -s unlimited`: `ray_tracing_1` keeps a 15 MB
automatic array on the stack. Without it, the run segfaults at the first RT
step.

## Inputs

These follow the spec. The files in `tests/problems/<name>/` are:
* `params.txt`: `N`, `z_init`, `da_sim`, `nsteps`, and optionally
  `nray_max_factor` and `legacy_rt`.
  * `legacy_rt = 1` (the default) is the original code.
  * `legacy_rt = 0` turns on the fixes for bugs B1–B4, B8 and B9
    (PATCHES.md 8, `docs/port_notes.md` §5a).
  * `legacy_physics` (default 1, the original physics) is independent of
    `legacy_rt`. `legacy_physics = 0` turns on the atomic-physics fixes
    P1–P10 (PATCHES.md 9–11, `docs/port_notes.md` §5b): case B
    recombination cooling, HeI collisional excitation (2³S and, P10, the
    1¹S ground state), HeII dielectronic recombination, the HeIII case B
    rate, collisional-ionisation cooling, HeII case A recombination, the
    free-free Gaunt factor and interpolated rate tables.
* `gas.bin`: float32 `(9,N,N,N)` holding `u_gas(1:5)` and `w_gas(1:4)`.
* `sources.txt`: one source per line as `i j k Luv sfr`.

`make_test_problem.py` (numpy only) writes these problems:

| problem | gas | sources | nray_max_factor |
|---|---|---|---|
| uniform_1src | uniform, `rho = Omega_b/Omega_m` | (16,16,16), Luv=2e28, sfr=1 | 64 |
| uniform_8src | uniform | 8 fixed cells, Luv=1e28, sfr=0.5 | 64 |
| lognormal_16src | log-normal, seed 20121, sigma_ln=1.3, delta 0.0023..140 | 16 densest cells, Luv=1e28(1+0.1 rank), sfr=Luv/2e28 | 64 |
| uniform_8src_smallbuf | as uniform_8src | as uniform_8src | 4 (spec problem 4) |
| lognormal_16src_smallbuf | as lognormal_16src | as lognormal_16src | 1 (extra: actually merges) |

In every problem the gas starts neutral: `w = (XHy rho, YHe rho, 0, 0)` with
zero velocity. The energy is `E = rho*Tgas_min`, as in
`read_init_conditions`, so `T = temp_min`. At z=10, calc_sim_var gives
`temp_min = 1.513 K`. It falls to 1.483 K at the last step, because
`temp_min ∝ a^-2`.

## What the driver does (vs. `radhydro.f90`)

The driver runs these steps in order:

1. **Initialise.** `init_variables` is copied, but the clock starts at
   `z_init` from params (the original uses the compile-time `zinit=300`). The
   driver then calls `init_domain_gas`, `init_domain_rt`, `init_data`,
   `init_halo`, `init_rt` and `init_ray_tracing`, followed by `calc_sim_var`.
   * `init_hydro` and `init_pm` are replaced by plain allocations. `hoc_dm` is
     zero, so `make_data_2d` finds no particles.
   * `init_io` is replaced by opening `rt.log` (unit 200) and `ray.log`
     (unit 201) in `<outdir>`.
2. **Load the initial state.** `gas.bin` is copied into `u_gas` and `w_gas`,
   then `calc_sim_var` and `summarize_rt` are called.
3. **Loop over each hydro step.**
   * `calc_sim_var`.
   * The RT part of `calc_timestep`, copied verbatim, with `da_sim` fixed
     from params. `dt_sim` and `dt_yr` come from `calc_dt_from_da`, as in the
     `istep == 0` branch.
   * The two half-step clock advances, each followed by `calc_sim_var` (as in
     `radhydro`; hydro and PM are skipped).
   * Sources go into `rt_grid%Luv` and `rt_grid%sfr` as `calc_halo` and
     `add_halo_to_rt` would do it (real(4) per source, real(8) sums, no metal
     injection).
   * `radiative_transfer`, which rewinds the clock by one step and
     re-advances it per RT substep, then `summarize_rt`, then the output
     files.

`make_data_2d` and `write_data_3` still run unmodified on every RT substep.
`make_data_2d` clamps `u_gas(1) >= Dgas_min` as a side effect.

## Outputs (`<outdir>/`)

* `step_%04d/`: one directory per hydro step, containing:
  * `<field>.bin`: float32 little-endian, Fortran order:
    * `fHI fHeI fHeII temp rhog`: (N^3)
    * `ngam nu fgam tau`: (5,N^3)
    * `gamma lambda dx`: (3,N^3)
    * `u_gas`: (5,N^3)
    * `w_gas`: (4,N^3)
  * `stats.txt`: the spec keys.
    * `nray_active = Nray - ray_save%nfree`.
    * `ncst1`, `nmrg1`, `nspl1` and `ndel1` are `ray_save` values. These are
      the counts of the **last RT substep** of the hydro step, because
      `rt_step_3` resets them every substep. Per-substep values are in
      `ray.log`.
    * `ngam_rays_total` is the real(8) sum, in storage order, of
      `ray%ngam(1:5)` over all slots with `nside > 0` after `post_rt`. At
      that point every inactive slot has been zeroed.
  * `rt_stats.txt`: extra, non-spec values. These are the cumulative
    `rt_save` counters (`ngam`, `nion`, `ncol`, `nrec`, `gamma`, `lambda`,
    `tau`), the `ngam_ray(1:64)` thresholds, `Nray` (which includes the
    free-buffer allocation), `dtsec_rt`, `temp_min` and similar.
* `ray.log` and `rt.log`: the original logs (units 201 and 200).
* `stdout.log`: written by the run scripts.
* `Log/`: `nside_of_r.dat` and `r_of_nside.dat`.
* `D2D/`: per-substep slice dumps. The scripts delete this directory.

### Stored reference outputs (`outputs/`)

A full run writes 5.6 MB per step per problem, 280 MB in total, which is too
much to keep in git. `outputs/<problem>/` therefore keeps:
* `stats.txt` and `rt_stats.txt` for **all** 10 steps.
* `ray.log`, `rt.log`, `stdout.log` and `Log/`.
* The full binary fields of **step_0010 only**: 5.6 MB per problem, 24 MB in
  total. `uniform_8src_smallbuf` has no binaries, because its spec outputs are
  bit-identical to `uniform_8src` (see below).

`outputs_legacy_rt0/` holds the same subset from
`RT_SET=legacy_rt=0 ./run_all.sh`. `outputs_legacy_rt0_phys0/` holds it from
`RT_SET="legacy_rt=0 legacy_physics=0" ./run_all.sh`. Again,
`uniform_8src_smallbuf` has no binaries in either; with the fixes its fields
are still identical to `uniform_8src`.

To regenerate everything (about 15 s), run `./run_all.sh outputs-full`.
Output is deterministic: repeated runs are bit-identical, apart from the
timing lines in the logs.

## Sanity checks (N=32, 10 steps, `tests/analyze_outputs.py`)

* **I-front (uniform_1src).** The volume-weighted xHII rises linearly:
  6.8e-4, 1.38e-3, ... up to 7.16e-3 at step 10. The I-front radius (sphere
  with the volume of the cells where xHII > 0.5) grows 1.66, 2.27, 2.39,
  2.83, 2.83, 3.21, 3.33, 3.55, 3.60, 3.65 RT cells. That is roughly
  `t^{1/3}`, as expected for a photon-limited front where recombinations are
  negligible (t_rec ≈ 1.3e16 s, much longer than the run's ≈ 2.6e14 s). The
  expected radius from photons emitted, 263 cell-equivalents of H, is about
  3.9 cells. Tmax is about 2.7e4 K.
* **Photon budget.** "absorbed" counts photoionisations of H, HeI and HeII
  from `rt_save%nion`; "in rays" is `ngam_rays_total/N^3`.
  * Uniform problems: (absorbed + in rays)/emitted = 0.9996 to 0.9999 at
    every step.
  * lognormal: the ratio is 1.005 to 1.022. This is expected: in cells with
    `rho > rho_source = 50` the code ionises without consuming photons
    (`x = exp(-(d-50)/50)` "self-absorption" factor), and the densest cells
    reach delta = 140.
  * Ionisation budget: mass-weighted xHII / (nion + ncol - nrec)(H) = 1.03
    at step 1 and 1.004 at step 10.
* **Ray counts** (`nray_active` at step 10, and `nmrg1` / `nspl1` in the last
  substep):

  | problem | nray_active (step 10) | merges | splits |
  |---|---|---|---|
  | uniform_1src | 816 | 0 | 640 |
  | uniform_8src | 5056 | 0 | 1120 |
  | lognormal_16src | 18187 | 0 | 6804 |
  | uniform_8src_smallbuf | 5056 | **0** | 1120 |
  | lognormal_16src_smallbuf | 12831 | 1891 | 4500 |

  In lognormal_16src_smallbuf, merges start at step 5. Over all substeps of
  the run the total is 50900 merges, from ray.log.
* **Merging coverage.** `ngam_ray(j)` stays 0 until the number of counted rays
  exceeds `j*int(Nray_max/64.)`. Merging needs `ngam_ray(16) > 0`, which
  means more than `Nray_max/4` counted rays. That is 524k rays for factor 64,
  and 32768 for factor 4, so `uniform_8src_smallbuf`, with about 5000 rays,
  never merges.

  Its spec outputs are bit-identical to `uniform_8src`; only `Nray` and the
  buffer counts in `rt_stats.txt` differ. A factor scan gives:
  * Factor 2 makes lognormal_16src merge from step 8, with 17056 merges.
  * Factor 1 makes it merge from step 5, with 50900 merges.
  * The uniform problems never merge at factors 1, 2 or 4.

  Hence the extra problem `lognormal_16src_smallbuf` (factor 1).
* **Runtime checks.** The `-fcheck=all` debug build runs all five problems
  with no runtime errors, after patch 5 in PATCHES.md. `FPE=overflow` also
  passes. The `-fcheck=all` build is also clean with `legacy_rt = 0`, and in
  both modes its outputs are bit-identical to the `-O2` build.
* **`legacy_physics = 0`.** The `-fcheck=all` build is clean in all four
  (`legacy_rt`, `legacy_physics`) combinations, and its outputs are
  bit-identical to the `-O2` build. With `legacy_physics = 1` the outputs
  equal `outputs/` and `outputs_legacy_rt0/` exactly. For the physics
  impact (HeIII +21% in lognormal_16src, ionised-gas temperatures −0.5% to
  +1.8%, I-front unchanged), see `docs/port_notes.md` §5b.
* **`legacy_rt = 0`.** Merging in lognormal_16src_smallbuf now starts at
  step 7 instead of step 5, with 41775 merges in total instead of 50900. The
  uniform problems and the factor-64 problems still never merge. For the
  full impact, see `docs/port_notes.md` §5a.
* **Compiler sensitivity.**
  * `-O0` (debug) and `-O2`: every output is bit-identical (`compare.py`
    gives max_rel = 0) for all five problems.
  * `-O3 -march=native` (FMA, vectorisation): ray stats are identical and the
    fields agree to at most 1.3e-5 relative (`fgam`, lognormal_16src_smallbuf).
    Most fields agree to 1e-10 or better.
* **Runtime.** Each problem takes 1.3 to 4.3 s at `-O2` on one core, and 2.6
  to 7.7 s for the debug build. About 0.3 s of each run is the
  initialisation of the 2M-entry ray arrays.
* **Other grid sizes.** N=16 (lognormal_16src, 10 steps, `nstep_rt = 6`)
  runs cleanly. So does N=64 (uniform_1src, 3 steps, about 6 s per step,
  1.2 GB for 16.7M ray slots). At N=64 the I-front radius after step 1 is
  3.44 cells, about 1.72 cells in N=32 units, consistent with the 1.66 cells
  at N=32.

## Parameter values found in the reference code

* `cosmo_tools.f90: calc_sim_var` sets these on every call:
  * `sfr_min = 1e-6`
  * `ngam_min = 1e-6`
  * `fgam_min = 1e-4`
  * `seg_min = 1e-4`
* `ray_tracing_tools.f90`:
  * `nside_init = 4`, so `nray_init = 192`.
  * `nside_min = 4`, `nside_max = 8192` (`lside_max = 13`).
  * `ray_per_cell = 1`.
  * Merge cones use `nside_merge = 2`, giving 48 cones.
* Time step at z=10, da=1e-4, N=32:
  * `nstep_rt = 12` in every step.
  * `dt_sim ≈ 5.45e-3` code units, `dt_yr ≈ 8.09e5` yr.
  * `dtsec_rt ≈ 2.126e12 s` at step 1 and 2.137e12 s at step 10.
* Domain decomposition with `Ncpu = 1`: n1 = 2, so there are 8 domains of
  16^3 cells, processed in order 1..8.

## HEALPix

The original links HEALPix 3.11 (`pix_tools`). The harness compiles the same
version's Fortran sources (see PATCHES.md).

`make test-healpix` compares the Fortran routines with healpy 1.20.1 (C++
HEALPix) for nside = 2^0 .. 2^13. The test inputs are:
* All pixels for small nside; for large nside, 128 edge pixels plus 20000
  random pixels.
* For each nside, 27009 vectors: random, special axes, float32-rounded, and
  unnormalised.

Results:
* `vec2pix_nest`: exact agreement everywhere.
* `pix2vec_nest`: the vectors differ from healpy by up to 1.1e-15, a few
  ulp. Round-trips are exact.

The C++ port should replicate the Fortran 3.11 algorithm rather than call
healpy or C++ HEALPix if it wants bit-identical ray directions. The
differences are usually invisible once directions are rounded to real(4), but
not always.

## Caveats and observations relevant to the C++ port

* **ifort vs gfortran.** The original was built with
  `ifort -O3 -ipo -ftz -fpp -openmp` and many threads. The ray-tracing order
  is nondeterministic with more than one thread: `set_domain_rt` hands out
  domains dynamically, and chain order and float sums depend on it. The
  harness is single-threaded and deterministic. Expect round-off-level
  differences from ifort builds because of FMA/x87, `-ftz` and math-library
  `exp`/`log`/`**`.
* **Deleted rays count towards `ngam_ray`.** `ray_remove` returns a slot to
  the free buffer without zeroing it. `calc_ray_1` scans all `Nray_max` slots
  for `nside > 0`, so rays deleted earlier in the same hydro step are still
  counted in `ngam_avg` and `ngam_pdf`. Also, `ngam_avg` is divided by `Nray`,
  which includes the free buffer. Slots are zeroed only in
  `ray_buffer_deallocate` (in `post_rt`). This behaviour is kept, and the C++
  port must reproduce it when the thresholds are active.
* **`ngam_ray` persistence.** When there are not enough rays, `ngam_ray(j)`
  keeps its previous value, which starts at 0. Before any threshold is set,
  every ray is cast or split, never merged.
* **IEEE semantics needed.** In `photo_cell`, a cell with exactly zero
  electrons (fully neutral, `w` rounding up to `nH`) evaluates
  `ne/ne_prev` with `ne_prev = 0`, which gives Inf or NaN. The `do while`
  test then continues (Inf > 1e-3) or exits (NaN). The C++ port must use the
  same comparisons, with no `-ffast-math`.
* **Histogram overrun.** The `post_rt` histogram overrun (patch 5) and the
  uninitialised `jrt` in `calc_gas_3` (patch 6) are genuine bugs in the
  original. Both affect diagnostics only.
* **Stats are last-substep values.** `ncst1`, `nmrg1`, `nspl1` and `ndel1`
  in stats.txt describe the last substep only.
* **Data race.** `end_domain_rt` updates `domain(1,:)` outside the critical
  section. This is a race under OpenMP and harmless with one thread.
* **Photons are not conserved by design.** This happens in four places:
  * The consumption factor `x < 1` above `rho_source` (see the photon
    budget).
  * Photons dumped into an optically thick source cell by `ray_source`
    (`exp(-min tau) < 1e-2`) are absorbed there, and the rest is not carried
    by any ray.
  * `ray_redshifting` moves photons below 13.6 eV out of the lowest group.
  * Faint rays (below `fgam_min`/`ngam_min`) are deleted.
* **`make_data_2d` writes to `u_gas`.** It clamps `u_gas(1)` to `Dgas_min`
  as a side effect. The generator keeps `delta >= 1e-3`, so the clamp is
  never active.
