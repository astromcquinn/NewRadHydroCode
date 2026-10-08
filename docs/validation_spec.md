# Validation spec (Phase 1)

Both the Fortran harness (`reference/harness`) and the C++ port must implement
exactly this, so their outputs can be compared file-by-file.

## Configuration

* Cosmology and RT parameters: exactly `reference/Radhydro/cosmo_parameters.f90`
  (box = 5 Mpc/h, Omega_m = 0.27, Omega_b = 0.045, h = 0.7, Nfreq = 5, ...).
* `Ngas = Nrt = N` (default `N = 32`), so `N_rt_gas = 1`.
* `Nray_max = nray_max_factor*N^3`, default factor 64 as in the original (it
  enters the `ngam_ray` thresholds and the ray-buffer logic — keep both).
  `params.txt` may set `nray_max_factor`.
* `legacy_rt` (default 1): 1 runs the original RT bookkeeping. 0 switches on
  the fixes for bugs B1–B4, B8 and B9 (`docs/port_notes.md` §5a). Both modes
  are part of the validation.
* Ray parameters as in `ray_tracing_tools.f90` (`nside_init=4`,
  `lside_max=13`, `ray_per_cell=1`).
* `RAYFAST` undefined (the default code path). `METAL`, `FREQ1` undefined.
* Fortran reference runs with **one OpenMP thread** (`Ncpu = 1`).
* No hydro and no gravity: the gas changes only through RT physics.
* Time stepping: start at `z_init = 10`; every "hydro step" has fixed
  `da_sim = 1e-4`; `dt_sim`, `dt_yr`, `dt_sec`, `nstep_rt`, `dt_rt`, `dtsec_rt`
  computed exactly as in `radhydro.f90: calc_timestep` (RT part) and
  `cosmo_tools.f90`. Default `nsteps = 10` hydro steps.
* Each hydro step then calls `radiative_transfer` as in `radhydro.f90`
  (with `idata1 > 1` treated as true and sources always present).

## Inputs (`tests/problems/<name>/`)

* `gas.bin`: raw little-endian float32, Fortran order, shape `(9, N, N, N)`:
  components `u_gas(1:5)` then `w_gas(1:4)` in code units at `z_init`
  (see the header comment of `radhydro.f90`).
* `sources.txt`: one source per line: `i j k Luv sfr` with 1-based RT cell
  indices; these are assigned to `rt_grid%Luv`/`%sfr` exactly as
  `halo_tools.f90: add_halo_to_rt` would (metal injection skipped).
* `params.txt`: `key = value` lines (`N`, `z_init`, `da_sim`, `nsteps`,
  optional `nray_max_factor` and `legacy_rt`). Both drivers accept
  `--set key=value` on the command line, which overrides or adds an entry:
  `rt_driver <prob> <out> --set legacy_rt=0`, or
  `rt_harness.x <prob> <out> [nsteps] --set legacy_rt=0`. The problem
  directories themselves are not duplicated per mode.

Problems:
1. `uniform_1src`: uniform gas at mean baryon density, neutral, T = temp_min,
   one source in the centre cell.
2. `uniform_8src`: as 1 with 8 sources at fixed positions.
3. `lognormal_16src`: deterministic log-normal density field (seeded), neutral,
   16 sources in the densest cells.
4. `uniform_8src_smallbuf`: as 2 but `nray_max_factor = 4` (in practice never
   merges; kept as a buffer-logic test).
5. `lognormal_16src_smallbuf`: as 3 but `nray_max_factor = 1`; merges from
   step 5 on — the merge-path test.

The canonical inputs are the files committed in `tests/problems/` (generated
by `reference/harness/make_test_problem.py`, seed 20121). Reference outputs
are in `reference/harness/outputs/` (`legacy_rt = 1`) and
`reference/harness/outputs_legacy_rt0/` (`legacy_rt = 0`, stored the same
way).

## Outputs (`<outdir>/step_%04d/`, after each hydro step)

Raw little-endian float32, Fortran order (component index fastest), one file
per field named `<field>.bin` (e.g. `fHI.bin`, `ngam.bin`):
`fHI, fHeI, fHeII, temp, rhog` (N^3), `ngam, nu, fgam, tau` (Nfreq,N^3),
`gamma, lambda, dx` (3,N^3) — i.e. the `rt_grid` fields — plus
`u_gas` (5,N^3) and `w_gas` (4,N^3).

`stats.txt`: `key = value` lines: `istep, redshift, scalefactor, nstep_rt,
nray_active, ncst1, nmrg1, nspl1, ndel1, ngam_rays_total` (sum of
`ngam` over all live rays, all groups, double precision).

## Acceptance

`tools/compare/compare.py ref_dir new_dir` reports per-field max/mean
absolute and relative differences per step and the ray-stat differences.
Target: ray stats identical; float fields agree to ~1e-5 relative (round-off
level) on all problems. Any larger discrepancy must be explained.

Achieved, and required from now on: bit-identical results in both modes
(`compare.py --rtol 0`, `tests/compare_reference.py`), checked by the
CTests `ref_<problem>` and `ref_legacy0_<problem>`. With `legacy_rt = 1` the
outputs must also stay identical to `reference/harness/outputs`.
