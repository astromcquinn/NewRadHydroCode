# Port notes: Fortran RadHydro hydro → C++ `radray`

Phase 1b, hydro part. This covers `hydro_tools.f90` compiled with `-DMF` (the
moving-frame scheme, as the original Makefile sets it), the gas domain
decomposition of `domain_tools.f90`, the transposes, and the hydro parts of
`radhydro.f90` (`calc_timestep`, `radhydro`) and of the summaries. The
general method (precision rules, constant folding, index conventions) is in
`docs/port_notes.md`; only what is specific to hydro is repeated here.

**Validation result:** the C++ port is **bit-for-bit identical** to the
gfortran harness (`reference/harness_hydro`, one thread) on all six hydro
problems in `tests/problems_hydro`:

* all `u_gas`, `w_gas`, `v_gas` files of all steps (step 0 included):
  195 field files, `compare.py --rtol 0 --atol-frac 0` reports 0 differing
  values;
* every number of every `stats.txt` (time step, `c_limit`, `cmax_gas`,
  `vmax_gas`, `vrms_gas`, mass/momentum/energy/`w`/`v` totals);
* the `summarize_gas` unit-100 records (`50es13.5`), byte for byte;
* the gas domain table and the `set_domain_gas` order.

The gfortran `-O0 -fcheck=all` harness build gives the same bits as `-O2`
and has no runtime errors. The C++ ASan/UBSan build is clean and also
bit-identical. CTests: `hydro_ref_*` (against the stored outputs in
`reference/harness_hydro/outputs`), `hydro_vs_fortran_live` (fresh gfortran
build, all steps and fields), `hydro_physics`, `hydro_sanitize_*`.

## 1. Routine map

| Fortran | C++ (`include/radray/hydro.hpp`, `src/hydro.cpp`, `src/domain_gas.cpp`) |
|---|---|
| hydro_tools.f90 `init_hydro` | `init_hydro` (allocates and zeroes `GasGrid::u,w`, `HydroState::v`; the first-touch bookkeeping and its log are not ported) |
| `hydro_step(flag)` | `hydro_step(..., SweepOrder::Forward/Reverse, ...)` |
| `smooth` (`#ifdef MF` branch) | `smooth` (internal) |
| `euler` | `euler` (internal; the dead `if (.false.)` first-order half update is not ported) |
| `advect` (`#ifdef MF` branch) | `advect` (internal) |
| `freeze` | `freeze` (internal) |
| `smooth_cells_1` / `smooth_cells_2` | `smooth_cells<float>` / `smooth_cells<double>` (the unused `v` argument is dropped) |
| `check_temp` | `check_temp` (internal) |
| `transpose_ij`, `transpose_ik` | same names (same blocked loop order, nb = 8) |
| `summarize_gas` | `summarize_gas` → `GasSummary` (record + printed table); `format_gas_record` writes the `50es13.5` record |
| MKL DFTI descriptors in `hydro_step`/`smooth` | class `Fft1d` (FFTW3 single precision, see §3) |
| domain_tools.f90 `init_domain_gas`, `set_domain_gas`, `end_domain_gas` | same names, using the `Domains` struct of `simulation.hpp` (`box` = `indx_dom_gas`, `domain` = `domain_gas`) |
| radhydro.f90 `calc_timestep` (hydro part) | `calc_timestep_hydro` (`TimestepMode::Fixed` / `Courant`) |
| radhydro.f90 `radhydro` (hydro part) | `radhydro_hydro_step` (no `particle_mesh`, no RT, no summaries) |
| radhydro_global.f90 `cmax_gas, vmax_gas, vrms_gas, c_limit, dt_hy`, `v_gas` | `HydroState` |
| radhydro_global.f90 `ra_max, dt_gas_factor, dt_pm_factor`; hydro_tools `rmf, cfl` | `hydrop::` constants |
| cosmo_tools.f90 `calc_sim_var`: `Tmf_min = Tmf_max = 1` | added to `Units` / `calc_sim_var` (the only edit to shared files besides CMake) |
| radhydro.f90 main loop (hydro only) | `apps/hydro_driver.cpp` |

Not ported: timing/`omp_get_wtime` logging, `init_hydro`'s first-touch
table, the standard-output formatting of `hydro_step` and `calc_timestep`
(the C++ prints the same quantities in its own format). `summarize_dm`,
`particle_mesh` and the data/IO routines belong to the PM/IO ports.

## 2. Layout

* `u_gas(5,N^3)`, `w_gas(4,N^3)` stay in `GasGrid` (unchanged struct);
  `v_gas(3,N^3)` is `HydroState::v`. Component index fastest, cells in
  Fortran order, so a pencil `u_gas(:,:,j,k)` is a contiguous `5*N` float
  block, exactly as the Fortran passes it to `smooth/euler/advect/freeze`.
* Pencil kernels use 1-based cell indices `i` with the original periodic
  index arithmetic (`im = 1+mod(i-2+N,N)`, ...), so every statement maps one
  to one. Per-pencil scratch arrays (`a,c,p,uh,vh,wh,fr,fl,...`) live in a
  `Sweep` struct allocated once per `hydro_step`.
* Requirements inherited from the original and checked with an exception:
  `Ngas` even (euler's odd/even cell pairs) and a multiple of 8 (the
  transposes' 8-blocks run out of bounds otherwise).

## 3. Precision and Fortran semantics (hydro specifics)

* real(4) arrays `u`, `v`, `w` against real(8) locals: every promotion point
  is kept. Notably:
  * `sum(u(2:4,i)**2)` on a real(4) array is a **float** sum of float
    squares; the same on `uh` (real(8)) is a double sum;
  * `(v(:,ip)-v(:,im))/2` (euler and advect first half) and
    `(u(:,im)-u(:,i))/4` (`smooth_cells_1`) are float operations widened
    afterwards; the same expressions on `vh`, `uh` are double;
  * `uv(:,i) = u(:,i)*v(1,i)` in advect is a **float** product stored in a
    real(8) array; `v(1,i)*v(1,ip) > 0` is a float product;
  * `(/vx(i),vy(i),vz(i)/) - v(:,i)` in smooth is a float subtraction;
  * every assignment to `u`, `v`, `w` rounds to float at that statement,
    e.g. the full update stores `u(:,i)` and then subtracts the expansion
    term from the stored float `u(5,i)`, and the boundary update adds
    `fu(:,1)*dth` to the already-rounded `u(:,1)`.
* Single-precision literals: `0.1` in `min(max(dvg,-0.1/dt),0.1/dt)`,
  `2E8` in `c_limit`, `1E5` (km/s), `2.740`, `0.407`, `2.242`, `1.5` in the
  recombination fit of `summarize_gas`, `Zsolar = 0.02`, `ra_max = 0.1`.
  `bscale = 1./Ngas` is a real(4) division.
* **MAX/MIN.** gfortran 13 (-O0 and -O2) propagates NaN in either argument,
  and on ties returns the **first** argument of MAX and the **second** of
  MIN (so `max(-0.,0.) = -0.`, `min(-0.,0.) = +0.`); 3-argument MAX is
  `max(max(a,b),c)`. `gmax`/`gmin` in `hydro.cpp` reproduce this; it
  matters for signed zeros written back to `w_gas` by `summarize_gas`
  (`w = max(w,0.)`). (`f::max` in `fortran.hpp` differs only for NaN in the
  second argument and is not used here.)
* `wg = exp(-((pi*rmf*(i-1)/ngasr)**2)/2)`: `pi*rmf` is folded at compile
  time, which is an exact IEEE product, then `*(i-1)`, `/ngasr`, square.
* **FFT.** MKL DFTI is replaced by FFTW3 single precision in both the
  harness (`mkl_dfti_shim.f90`) and the C++ (`Fft1d`), identically:
  `FFTW_ESTIMATE` r2c/c2r plans on an `fftwf_alloc_real(n+2)` buffer, data
  copied in/out, backward result times the real(4) scale. MKL's default
  1-D real storage (CCS) is FFTW's r2c layout, which is what `smooth`'s
  kernel loop assumes. MKL would differ at round-off level, like ifort.
* libm `exp`, `sqrt`, `pow` are the same glibc calls as gfortran's.

## 4. Interpretations

1. **Driver step sequence** (`hydro_harness.f90` = `hydro_driver.cpp`):
   `calc_sim_var`; `calc_timestep` (hydro part); `radhydro` without
   `particle_mesh` and RT; `summarize_gas`. The initial state is read from
   `u_gas.bin`, `w_gas.bin` and optionally `v_gas.bin` (0 otherwise, as
   `init_hydro` leaves it), followed by `calc_sim_var` and `summarize_gas`.
2. **Time step.** `fixed`: `da_sim` given, `calc_dt_from_da` (the
   `istep == 0` form). `courant`: the original `istep > 0` branch from the
   first step on, with `amax_dm = 0`, i.e. `dt_ap = dt_pm_factor/sqrt(0) =
   +Infinity` (computed as such in both codes). The original `istep == 0`
   branch (`da = a*ra_max/(1-ra_max)`, a 10 % expansion step) is a
   cosmological start-up choice and is not used by the test problems.
   `calc_timestep` also converts `vmax_gas`, `vrms_gas` to km/s in place;
   ported (they are overwritten by the next `summarize_gas`).
3. **Stats.** The totals in `stats.txt` are diagnostics of the drivers,
   computed identically in both: `momentum = Σ (u(2:4) + rho vg)`,
   `energy = Σ (u5 - |u(2:4)|²/2rho + |u(2:4)+rho vg|²/2rho)` (lab frame).
4. **Test problems** are dimensionless set-ups in code units, starting at
   z = 10; nothing in the solver is altered (no expansion terms exist in the
   supercomoving formulation for gamma = 5/3). See
   `reference/harness_hydro/README.md` for the list and branch coverage.

## 5. Bugs and quirks in the original (all preserved unless noted)

* **H1. `advect`/`euler` never terminate when the signal speed of a whole
  pencil is exactly zero.** The sub-step count is
  `nt = ceiling(amax*dt/cfl)` and the loop exits on `it == nt`; with
  `amax = 0` (advect: the grid velocity `v(1,:)` of the pencil is exactly
  zero; euler: no pressure and no velocity) `nt = 0` and the loop runs
  forever. Any problem with a direction of exactly zero velocity triggers
  it in advect (e.g. a 1-D tube at rest: the transverse pencils, and the
  first x sweep, have `vg = 0` exactly). Cosmological runs never hit it.
  The C++ throws an exception where the Fortran would hang; all test
  problems carry a uniform bulk velocity with three non-zero components.
* **H2. The moving frame is smoothed along the current sweep only.**
  `smooth` replaces all three components of `v_gas` by the 1-D Gaussian
  (rmf = 8 cells) smoothing of the total velocity *along the pencil*. After
  a transverse sweep the grid velocity equals the local fluid velocity in
  the other directions. Results therefore depend on the orientation of a
  problem relative to the sweep order (z,y,x / x,y,z): sod_x and sod_y agree
  to 1e-7, sod_z does not. `cmax_gas` (`max |vel - vg| + cs`) and hence
  the Courant step also depend on which sweep came last: sod_z takes
  shorter steps (t = 1.76 vs 2.24 after 5 steps).
* **H3. MF temperature weighting is disabled.** `calc_sim_var` sets
  `Tmf_min = Tmf_max = 1`, so `wt = 1/sqrt(T)` is clamped to 1 everywhere
  and `smooth` is a plain Gaussian smoothing of the velocity.
* **H4. `vmax_gas` ignores vz.** `summarize_gas` uses
  `max(maxval(x4(2:3)), maxval(abs(x3(2:4))))`: the maximum of vz is
  missing (the minimum is included). Only used in the log of the next
  `calc_timestep`.
* **H5. Swapped normalisations in the printed velocity table.**
  `vv = x1(2:4)*.../x1(1)` divides the volume sum by the mass sum and
  `vm = x2(2:4)*.../Ncell` the mass-weighted sum by the cell count (the
  density, sound speed, T, ... lines do it the other way round). Also
  `x3(2:4) = 0`, `x4 = 0` initialise min/max velocity at 0. Print only.
* **H6. `summarize_gas` writes to the state:** `u_gas(1) = max(u_gas(1),
  Dgas_min)` and `w_gas = max(w_gas,0.)`. So a cell set below `Dgas_min`
  in the initial conditions is raised before the first sweep (and the
  `u(1) < Dgas_min` branch of the sweeps is reached only by cells emptied
  during a sweep). Clumping factors are NaN (0/0) for fully neutral gas.
* **H7. `euler`'s `rho > Dgas_frz` test after the `c >= c_limit` branch
  uses the density of the right neighbour** (`rho` is last set for
  `j = ip` in the k-loop), and the cell's freezing speed is the one
  recomputed after `smooth_cells_1`. `cmax` (which sets `nt`) is tracked in
  the first-order pass only. Preserved.
* **H8. `smooth_cells` does not conserve the passive scalars or mass
  exactly.** `w` is rescaled by `r = 1 ± du/rho` per cell, which conserves
  Σw only if the three cells have the same abundances, and `u(1)` is
  clamped to `Dgas_min` afterwards (adds mass). In random3d_freeze the `w`
  totals drift by 1.5e-3 in 10 steps; mass stays conserved to 5e-8 there.
* **H9. Lab-frame momentum and energy are not exactly conserved** by the MF
  scheme: the frame-velocity source terms (`v - a*dvg`, the
  `p/(g-1)(1-exp(-(g-1) dvg))` term) are not in conservation form, and
  `check_temp` imposes a 1 K floor (energy redistributed between three
  cells). Measured drifts over the runs: 2e-4–2e-3 (momentum), 4e-5–4e-3
  (energy), 6e-8 for the sound wave.
* **H10. `gamma_ideal = 5./3` is a real(4) division**
  (1.6666666269...). The supercomoving equations are free of expansion
  terms only for gamma = 5/3 exactly; the tiny inconsistency (2.4e-8) is
  ignored by the original and here.
* **H11. `hydro_step` log format `2f`** has no width (ifort extension,
  gfortran compile error); fixed in the harness copy only (PATCHES.md 1).
* **H12. Grid-size assumptions without checks:** `Ngas` must be even
  (euler's odd/even pairs) and a multiple of 8 (transposes); `init_domain_gas`
  leaves `n1` undefined for odd `Ngas`. The C++ throws instead.
* **H13. `c_limit = max(2*cmax_gas, 2E8/vel_unit)`** uses `cmax_gas` from
  the previous `summarize_gas`, which excludes cells below `Dgas_frz`;
  fast low-density cells therefore exceed it and are smoothed with their
  neighbours (`random3d_freeze`).
* **H14. Domain scheduling with threads:** `set_domain_gas` spins until a
  domain is free and `end_domain_gas` updates the table outside the
  critical section (same race as B14 for RT). Irrelevant with one thread;
  the C++ throws if no domain is eligible (the Fortran would spin forever).
  The gas domains are not used by the hydro sweeps themselves (only by
  `data_tools.f90` and `pm_tools.f90`).
* **H15. ifort vs gfortran:** as for RT (B13), the production build
  (`-O3 -ipo -ftz`, MKL) is not bitwise reproducible with IEEE gfortran and
  FFTW; the reference is the gfortran harness. With more than one OpenMP
  thread the hydro sweeps are still deterministic (pencils are
  independent, `vmax` is only logged), apart from `summarize_gas`'s
  reduction order.

## 6. Physics sanity checks (`tests/hydro_physics_check.py`)

* **Sod tubes** (N = 32, 5 Courant steps, t = 2.24 for x/y and 1.76 for z):
  L1 errors against the exact Riemann solution (gamma = real(4) 5/3, both
  periodic interfaces, advected by the bulk velocity) are 6–7 % in density
  and pressure and 5.7 % of c_L in velocity, dominated by the smeared
  contact; the shock is within 0.3 cells of the exact position
  (23.22 vs 22.97 for x/y, 21.75 vs 21.46 for z). Profiles are uniform
  across the tube to 1.5e-8.
* **Sound wave** (A = 1e-3, k = 2π/32): phase speed 2.2108 against
  c0 + vbx = 2.2 (+0.5 %, dispersion of the 2nd-order scheme at 32 cells
  per wavelength), amplitude damped by 0.6 % over t = 3.
* **Conservation:** mass to 1e-8 (float round-off) in every problem; see
  H8/H9 for `w`, momentum and energy.

## 7. Timings (one core, -O2, N = 32)

| problem | Fortran harness | C++ |
|---|---|---|
| sod_x / sod_y / sod_z (5 steps) | 0.40 / 0.40 / 0.39 s | 0.36 / 0.35 / 0.31 s |
| sound_x (24 steps) | 1.89 s | 1.53 s |
| random3d / random3d_freeze (10 steps) | 0.93 / 0.99 s | 0.94 / 0.95 s |

About 0.3 µs per cell per sweep (six sweeps per step), including the FFT
smoothing and `summarize_gas`.

## 8. Not faithfully portable / open points

* Only `Ncpu = 1` and the `-DMF` code paths are validated. Without `-DMF`
  `smooth` sets `v = 0` and `advect` is empty; not ported.
* The MKL FFT is replaced by FFTW (both sides); bit-identity with an actual
  MKL build is not possible.
* `hydro_step`/`calc_timestep` standard output is not reproduced
  character for character; the unit-100 `summarize_gas` record is.
* The harness's first Courant step uses the `istep > 0` formula (see §4.2).
