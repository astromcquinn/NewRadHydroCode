# tools/compare

Compares two RT output directories written in the phase-1 format of
`docs/validation_spec.md`. Typically the reference directory comes from the
Fortran harness and the new one from the C++ port.

Requires python3 and numpy. matplotlib is only needed for `--plots`.

## Usage

```
python3 tools/compare/compare.py REF_DIR NEW_DIR [options]

  --N N                grid size (default: inferred from the size of fHI)
  --nfreq NF           photon groups (default: inferred from the size of ngam)
  --rtol R             relative tolerance on float fields (default 1e-5)
  --atol-frac F        absolute floor per field = F * max|ref field| (default 1e-6)
  --field-rtol NAME=R  per-field rtol override, repeatable (e.g. dx=1e-4)
  --stats-int-tol R    relative tolerance on integer ray stats (default 0 = exact)
  --stats-rtol R       relative tolerance on float stats (default: --rtol)
  --fields a,b,...     compare only these fields
  --ignore a,b,...     skip these fields
  --steps 1,2,...      compare only these steps
  --plots OUTDIR       write PNGs: mid-plane (k=N/2) slices of fHI and temp
                       (ref / new / relative diff) and history curves
  --plot-steps S       'last' (default), 'all', or a comma-separated list
  --json FILE          write the full machine-readable report
  -q, --quiet          print only failures, histories and the summary
```

Exit status: `0` means everything is within tolerance. `1` means at least one
field or stat exceeded its tolerance. `2` means a structural problem: missing
steps or files, inconsistent sizes, or bad arguments.

Example:

```
python3 tools/compare/compare.py runs/fortran/uniform_1src runs/cpp/uniform_1src \
        --plots plots/uniform_1src --json uniform_1src.json
```

## Input format

* Step directories are named `step_NNNN`. Only steps present in both runs are
  compared, and steps present in just one run are reported (exit status 2).
* Fields are raw little-endian float32 in Fortran order, with the component
  index fastest: `fHI fHeI fHeII temp rhog` (N^3), `ngam nu fgam tau`
  (Nfreq,N^3), `gamma lambda dx` (3,N^3), `u_gas` (5,N^3) and `w_gas` (4,N^3).
  A file may be called `<name>` or `<name>.{bin,dat,raw,f32}`. Any other file
  whose size is a multiple of 4N^3 bytes is also compared, as an extra field.
* PM outputs (`apps/pm_driver`, `reference/harness_pm`) are also understood:
  `v_gas` (3,N^3) and `phi` (N^3) are known fields, N is inferred from
  `u_gas` when there is no scalar RT field, and `xv_dm` is read as a
  particle array `(6,Npart)` (component fastest, Npart from the file size).
  The worst particle value is reported as `(c=..,ip=..)` (1-based). For the
  bit-for-bit check use `--rtol 0 --atol-frac 0`.
* `stats.txt` holds `key = value` lines. Fortran `D` exponents are accepted.
  The integer keys `istep nstep_rt nray_active ncst1 nmrg1 nspl1 ndel1` are
  compared exactly by default. `redshift scalefactor ngam_rays_total` and any
  unknown numeric keys are compared with `--stats-rtol`.

## What is reported

For each step and each field:

| column | meaning |
|---|---|
| `max_abs` | max \|new − ref\| |
| `max_rel` | max \|new − ref\| / max(\|ref\|, \|new\|, atol) |
| `mean_abs` | mean \|new − ref\| |
| `n_diff` | number of values that differ at all (bitwise, as float32) |
| `n_fail` | number of values with \|new−ref\| > atol + rtol·max(\|ref\|,\|new\|), plus NaN/Inf mismatches |
| worst | 1-based Fortran location `(c,i,j,k)` of the worst value, with both values (a NaN mismatch takes precedence) |

`atol = atol_frac * max|ref field|`. Values that are tiny compared with the
field's dynamic range count as zero, so for example a direction component of
1e-9 that comes out of cancellation does not count as a 100% error. Lower
`--atol-frac` to make the comparison stricter.

For each step and each run, the global histories are:

* `xHII_V/M`, `xHeII_V/M`, `xHeIII_V/M`: volume- (`_V`) and mass-weighted
  (`_M`, weight `rhog`) ionised fractions. Here xHII = 1−fHI, xHeII = fHeII
  and xHeIII = 1−fHeI−fHeII.
* `T_V`, `T_M`: volume- and mass-weighted mean temperature.
* `ngam_grid`: Σ ngam over cells and groups. This is the grid photon density
  after `radiative_physics`, so it counts the photons that remain after
  absorption in the last substep.
* `ngam_rays_total`, `nray_active`: copied from `stats.txt`.

The summary prints the worst field overall, the **first step at which any
field or stat goes beyond tolerance**, and PASS or FAIL.

## Self-test

```
python3 tools/compare/test_compare.py
```

The test builds synthetic N=8, Nfreq=5 runs and injects these differences:
a scalar perturbation, a multi-component perturbation, a raw-offset poke that
checks the Fortran index convention, a NaN, integer and float stat changes,
missing steps and files, and sub-tolerance round-off. It checks detection,
worst-cell location, first-divergence step and exit codes.

## Notes on tolerances

Round-off between gfortran and C++ shows up at the 1e-7 to 1e-6 relative level
in float32 outputs. It can grow through the iterative chemistry solver and
through threshold decisions. For example, a ray whose photon count sits
exactly on `ngam_ray(16)` can switch between cast and merge. A single such flip
changes the integer ray stats and then propagates. When the ray stats differ,
the first divergent step and the location of the worst cell are the place to
start debugging.
