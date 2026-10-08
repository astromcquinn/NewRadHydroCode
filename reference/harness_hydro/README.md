# Fortran reference harness (hydro only)

This directory builds Hy Trac's original RadHydro hydro solver
(`hydro_tools.f90` with `-DMF`, the moving-frame scheme) with gfortran and
runs the `radhydro.f90` step sequence **without PM and without RT**, on one
thread. Its outputs are the reference for the C++ port (`src/hydro.cpp`,
`src/domain_gas.cpp`, `apps/hydro_driver.cpp`); see
`docs/port_notes_hydro.md`.

```
reference/harness_hydro/
  src/hydro_harness.f90    the driver (replaces radhydro.f90)
  src/hydro_tools.f90      copy of ../Radhydro/hydro_tools.f90, one format fix
  src/mkl_dfti_shim.f90    module MKL_DFTI on FFTW3 (single precision)
  Makefile                 also compiles ../Radhydro/{physical_constants,
                           cosmo_parameters,general_tools,cosmo_tools,
                           domain_tools}.f90 unmodified and reuses
                           ../harness/src/{radhydro_global,ifport_shim}.f90
  make_hydro_problems.py   writes tests/problems_hydro/* (numpy only)
  run_problem.sh           run one problem
  run_all.sh               build and run every problem
  store_reference.py       copy the stored subset into outputs/
  outputs/                 stored reference outputs (see below)
  PATCHES.md               every source change
```

## Build and run

Needs gfortran (tested with 13.3) and FFTW3 single precision
(`libfftw3f`, headers incl. `fftw3.f03`). No MKL.

```
make                 # build/N32/hydro_harness.x (-O2 -g)
make N=64            # build/N64/...   (Ngas is a preprocessor define, -DNGRID)
make debug           # build/debug_N32/ (-O0 -g -fcheck=all -fbacktrace -Wall)
./run_problem.sh ../../tests/problems_hydro/sod_x runs/sod_x [nsteps]
./run_all.sh [outroot=runs] [opt|debug]
```

`run_problem.sh` builds the executable for the problem's `N` if needed. The
executable refuses to run if the problem's `N` differs from the build, or if
it was built with `NCPU /= 1`.

## What the driver does (vs. `radhydro.f90`)

1. **Initialise:** `init_variables` (copied; the clock starts at `z_init`
   from `params.txt` instead of `zinit = 300`), then the original
   `init_domain_gas` and `init_hydro`, then `calc_sim_var`. PM, RT, halo,
   data, IO and FFT initialisation are not called.
   `domain_gas.txt` records the domain table and the order in which
   `set_domain_gas`/`end_domain_gas` hand out domains over three passes.
2. **Initial state:** `u_gas.bin`, `w_gas.bin` and (if present) `v_gas.bin`
   are read as given (replaces `read_init_conditions`); without `v_gas.bin`
   the grid velocity is 0, as `init_hydro` leaves it. Then `calc_sim_var`
   and `summarize_gas` (written to `step_0000/`).
3. **Each step:**
   * `calc_sim_var`;
   * `calc_timestep`, hydro part, verbatim except for the step size choice:
     * `timestep = fixed`: `da_sim` from `params.txt`, `dt_sim` and `dt_yr`
       from `calc_dt_from_da` (the form of the `istep == 0` branch);
     * `timestep = courant`: the code's own step, the `istep > 0` branch,
       used from the first step on: `dt_vg = 2*dt_gas_factor/cmax_gas`,
       `dt_ap = dt_pm_factor/sqrt(amax_dm) = +Infinity` (no particles),
       `dt_sim = min(dt_vg,dt_ap)`, `calc_da_from_dt`, `ra_max` limit,
       `calc_dt_from_da`. (The original `istep == 0` branch takes a 10 %
       expansion step `da = a*ra_max/(1-ra_max)` meant for the start of a
       cosmological run; it is not used.)
     * `vmax_gas`, `vrms_gas` are converted to km/s in place, as in the
       original (they are recomputed by `summarize_gas`).
   * `radhydro`, verbatim minus `particle_mesh` and the RT block:
     `transpose_ij`, `transpose_ik`, `hydro_step('forward')`, `transpose_ik`,
     `hydro_step('forward')`, `transpose_ij`, `hydro_step('forward')`, half
     clock step + `calc_sim_var`, `hydro_step('reverse')`, `transpose_ij`,
     `hydro_step('reverse')`, `transpose_ik`, `hydro_step('reverse')`,
     `transpose_ik`, `transpose_ij`, half clock step + `calc_sim_var`,
     `summarize_gas`.
   * output of `step_%04d/`.

Nothing in the solver is changed to make the problems "non-cosmological":
`hydro_tools.f90` works in supercomoving variables (Trac & Pen 2004:
`x_unit ∝ a`, `time_unit ∝ a^2`, `vel_unit ∝ 1/a`), in which the adiabatic
Euler equations for gamma = 5/3 have no expansion terms. The expansion enters
only through the time step (`calc_dt_from_da`), the temperature floor of
`check_temp` (`1/temp1_unit`, i.e. 1 K) and `c_limit = max(2 cmax_gas,
2E8/vel_unit)`. The test problems are set up in code units.

## Inputs (`tests/problems_hydro/<name>/`)

Raw little-endian float32, Fortran order (component index fastest):
`u_gas.bin` (5,N^3): rho, rho*(v - vg) (3), energy density in the frame vg;
`w_gas.bin` (4,N^3): rho_HI, rho_HeI, rho_HeII, rho_Z; optional `v_gas.bin`
(3,N^3): grid velocity vg. `params.txt`: `N`, `z_init`, `nsteps`,
`timestep` (`fixed` | `courant`), `da_sim`. Units: lengths in cells,
density in units of the mean matter density (mean gas density
`Dgas_avg = 1/6`), velocity in cells per code time. At z = 10,
`vel_unit = 134 km/s` and e/rho = 1 corresponds to 1.77e6 K.

| problem | setup | steps |
|---|---|---|
| sod_x, sod_y, sod_z | Sod-like tube (rho 1 : 0.125, p 1 : 0.1, c_L = 2) along x, y, z; periodic, two interfaces; bulk velocity (0.25,-0.15,0.1) rotated with the tube | 5, Courant |
| sound_x | travelling sound wave, A = 1e-3, one wavelength, c0 = 2, bulk (0.2,0.3,-0.25) | 24, da = 0.004 |
| random3d | log-normal density (sigma_ln 1.5), compressive velocity (rms 2.5) + bulk (10,-8,6), T = 300 K (rho/<rho>)^(2/3), near-empty and fast low-density cells, random fHI/fHeI/fHeII/Z, smoothed initial v_gas | 10, da = 0.002 |
| random3d_freeze | same field, bulk (0.3,-0.2,0.1), fast cells at (100,-80) | 10, da = 0.002 |

Every problem needs a non-zero velocity in all three directions: an exactly
zero grid velocity along a whole pencil makes `advect` loop forever (bug H1
in `docs/port_notes_hydro.md`).

Branch coverage (counted by the C++ driver over the run):

| problem | euler sub-cycled pencils (max nt) | advect sub-cycled (max nt) | smooth_cells_1/2 | c_limit branches | freeze | check_temp fixes |
|---|---|---|---|---|---|---|
| sod_x (sod_y, sod_z) | 9344 (9280, 4096) of 30720, nt 2 | 0 (1) | 0/0 | 0 | 0 | 0 |
| sound_x | 0 (1) | 0 (1) | 0/0 | 0 | 0 | 0 |
| random3d | 205 (3) | 8161 (2) | 0/0 | 0 | 0 | 28695 |
| random3d_freeze | 280 (2) | 0 (1) | 19/12 | 12 | 3 | 51723 |

## Outputs (`<outdir>/`)

* `step_%04d/` for the initial state (0000) and every step:
  * `u_gas.bin`, `w_gas.bin`, `v_gas.bin`: the full state (float32);
  * `stats.txt`: `istep, redshift, scalefactor, year, da_sim, dt_sim,
    dt_yr, dt_hy, c_limit, cmax_gas, vmax_gas, vrms_gas, vel_unit,
    temp1_unit`, and totals summed in double in k,j,i order: `mass`
    (sum rho), `momentum` (lab frame, sum u(2:4) + rho vg), `energy` (lab
    frame: thermal + rho |dv+vg|^2/2), `energy_grid` (sum u(5)), `thermal`,
    `w_total` (4), `v_total` (3);
  * `gas_summary.txt`: the record `summarize_gas` writes to unit 100
    (format `50es13.5`), prefixed with `summary = `.
* `domain_gas.txt`, `stdout.log`.

### Stored reference outputs (`outputs/`, 7 MB)

All `stats.txt` and `gas_summary.txt`, `domain_gas.txt` and `stdout.log`
for every problem, plus the binary fields of the last step for sod_x,
sound_x, random3d and random3d_freeze. `./run_all.sh` followed by
`./store_reference.py` regenerates them (about 6 s). Runs are deterministic.

## Checks done

* `make debug` (-O0, `-fcheck=all`): all problems run without runtime
  errors, and every output is bit-identical to the -O2 build.
* The C++ port is bit-identical on all problems, all steps, all fields,
  stats and summary records (CTest `hydro_ref_*`, `hydro_vs_fortran_live`).
* Runtime (-O2, one core): sod_* 0.4 s, sound_x 1.9 s, random3d* 1.0 s
  (hydro part; C++: 0.31-0.36, 1.53, 0.94 s).
