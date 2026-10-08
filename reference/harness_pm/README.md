# PM gravity reference harness (gfortran)

This harness builds the original particle-mesh code of RadHydro
(`pm_tools.f90`, with `-DMF` as the original Makefile sets) with gfortran,
single-threaded, using FFTW3 single precision instead of MKL. It runs the
gravity part of the `radhydro` step with no hydro update and no RT. The
C++ port (`src/pm.cpp`, `apps/pm_driver`) must reproduce its output bit for
bit.

Source changes are listed in [PATCHES.md](PATCHES.md): FFTW in `fft_tools`,
and the `Ncpu = 1` read bug in `read_init_conditions`.

## Build and run

Requirements: gfortran, and libfftw3-dev (`fftw3.f03` and `libfftw3f`).

```
make                       # NGAS=32 NDM=16 MF=1 -> build/G32_D16_MF1/pm_harness.x
make NGAS=64 NDM=32 MF=0   # other sizes / the non-MF branch
make debug                 # -O0 -fcheck=all
./run_all.sh [outroot] [debug]     # all tests/problems_pm -> runs/<p>
python3 store_outputs.py runs outputs   # refresh the stored reference
```

`Ngas`, `Ndm` and the MF branch are compile-time settings in the original.
The harness checks them against the problem's `params.txt` and stops with a
hint if they do not match. `run_all.sh` builds the right executable for each
problem.

```
pm_harness.x <problem_dir> <outdir> [nsteps]
```

The harness creates `<outdir>/Init`, a link to the problem directory, and
changes into `<outdir>`. The original routines then find their files where
the original expects them: `Init/xv_dm.init`, `Init/u_gas.init`, and `Pk/`
for `calc_power_spectrum`.

### params.txt

| key | meaning | default |
|---|---|---|
| `Ngas`, `Ndm` | grid and particle lattice sizes (must match the build) | 32, Ngas/2 |
| `mf` | 1 = `-DMF` build, 0 = non-MF build (must match the build) | 1 |
| `z_init` | initial redshift (the original uses the compile-time `zinit = 300`) | 300 |
| `nsteps` | number of steps | 10 |
| `da_sim` | 0 = the code's own `calc_timestep`; > 0 = fixed expansion step | 0 |
| `pk_every` | call `calc_power_spectrum` every n steps (0 = never) | 0 |

### Outputs

* `step_NNNN/` for NNNN = 0..nsteps. Step 0 is the state after
  `read_init_conditions`.
  * `xv_dm.bin`: float32 `xv_dm(6,Npart)`.
  * `u_gas.bin`: `(5,Ngas^3)`.
  * `v_gas.bin`: `(3,Ngas^3)`.
  * `phi.bin`: `density(1:Ngas,:,:)` after `particle_mesh`, which is the
    potential.
  * `w_gas.bin`: step 0 only.
  * `stats.txt`: `key = value` with 17 significant digits. Keys: `istep`,
    `nstep_rt`, `redshift`, `scalefactor`, `year`, `da_sim`, `dt_sim`,
    `dt_yr`, `amax_dm`, `arms_dm`, `vmax_dm`, `vrms_dm` (PM grid units, from
    `calc_kick`), `cmax_gas`, `vmax_gas`, `vrms_gas`, and the `summarize_dm`
    reductions `dm_{x,y,z}_{avg,min,max}`, `dm_v{x,y,z}_{avg,min,max}`
    (velocities in km/s).
* `Pk/pk.<ztag>`: written by the original `calc_power_spectrum`.
* `stdout.log`: the original's stdout, from `run_all.sh`.

## What the driver does (src/pm_harness.f90)

It mirrors `radhydro.f90`:

1. `init_variables` (shared copy, `harness_common/harness_timestep.f90`; it
   starts at `z_init`), then the **original** `init_domain_gas`,
   `init_hydro`, `init_pm`, `init_fft`, `calc_sim_var` and
   `read_init_conditions`, then `calc_sim_var`, `summarize_dm` and
   `summarize_gas` (original).
2. For each step:
   * `calc_sim_var`;
   * `calc_timestep` (shared verbatim copy `calc_timestep_mode`): on the first step the expansion step
     is limited by `ra_max` (da = a/9); after that it is
     `dt = min(1/cmax_gas, 0.5/sqrt(amax_dm))`, limited by `ra_max`. This
     step also computes the RT substep count, and converts the
     `amax/vmax/vrms` globals to physical units in place. Two changes:
     `da_sim > 0` fixes the step, and the zdata output flags are dropped;
   * half clock step + `calc_sim_var`;
   * **`particle_mesh`** (original);
   * half clock step + `calc_sim_var`;
   * `summarize_dm` and `summarize_gas` (both original);
   * optionally `calc_power_spectrum` (original).

   The hydro sweeps and the transposes are skipped, so the gas stays in x,y,z
   order, as the PM step expects. RT is skipped.

`summarize_gas` keeps its side effects on the state: `u_gas(1) >= Dgas_min`
and `w_gas >= 0`. (Until Phase 1c the driver had a copy of its time-step
part, `summarize_gas_speeds`; the results are identical.)

## Validation status

The C++ `pm_driver` output is byte-identical to this harness on every field
of every step, every stats value and every Pk file. This holds for:

* the four problems in `tests/problems_pm`;
* extra runs at Ngas = 48 (FFT size not a power of two) and Ngas = 128 /
  Ndm = 64;
* the `-O0 -fcheck=all` build of the harness, which is also bounds-check
  clean.

`reference/harness_pm/outputs/` stores the compact reference that CTest
uses (`pm_ref_*`): stats and the SHA-256 of every field for every step, the
last step's `xv_dm`/`phi`, and the Pk files. The hashes assume the same FFTW
build and CPU family. On another machine, rebuild and rerun the harness, and
compare with `tests/run_pm_all.sh`.
