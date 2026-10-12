# Port notes: RadHydro N-body / PM gravity → C++ `radray::pm`

This is Phase 1b of `docs/PLAN.md`. It follows the method and conventions of
`docs/port_notes.md`: Fortran precision and operation order are kept, and
the C++ is compared with a gfortran reference build. This document records:

* the routine map;
* the interpretations made where the Fortran was ambiguous;
* the quirks and bugs found in the original. They are preserved and flagged
  unless stated otherwise.

**Validation result:** `apps/pm_driver` is **bit-for-bit identical** to the
gfortran harness `reference/harness_pm`, single-threaded, with both using
FFTW3f in place of MKL. The comparison covers every step of the four problems
in `tests/problems_pm`:

* the fields `xv_dm`, `u_gas`, `v_gas`, the potential `phi`, and the initial
  `w_gas`;
* every stats value: clock, `amax/arms/vmax/vrms_dm`,
  `cmax/vmax/vrms_gas`, and the `summarize_dm` reductions;
* every power-spectrum file.

Extra runs at Ngas = 48 (FFT size not a power of two) and Ngas = 128 /
Ndm = 64 (262k particles) are also identical, as is the `-O0 -fcheck=all`
harness. The C++ is ASan/UBSan clean. `compare.py --rtol 0 --atol-frac 0`
passes.

## 1. Routine map

| Fortran | C++ (`namespace radray::pm` unless noted) |
|---|---|
| radhydro_global.f90 sizes: `Ngas, Hgas, Ndm, Npart, ngasr, ndmr, N_dm_gas, N_dm_gas_c, x_dm_gas` | `PMConfig` (`PMConfig::make(Ngas, Ndm, mf)`; `mf` = the `-DMF` switch) |
| radhydro_global.f90 arrays `xv_dm, ll_dm, hoc_dm, hoc_cpu, density, dm_cpu` | `PMState` (SoA vectors; `xv[(ip-1)*6+c]`, 1-based `ip`, `ll[0]` dummy) |
| the "Sim" globals `amax_dm, arms_dm, vmax_dm, vrms_dm, cmax_gas, vmax_gas, vrms_gas` | `SimClock` fields (Phase 1c; they were in `PMState` and `GasSpeeds`) |
| `u_gas`, `w_gas` | plain `std::vector<float>` arguments (the `GasGrid::u/w` layout). `v_gas(3,Ngas^3)` is a separate vector, so `grid.hpp` stays untouched for the hydro port |
| domain_tools.f90 `init_domain_gas` | the hydro port's `init_domain_gas` (`Domains`, `src/domain_gas.cpp`), passed to `init_pm` (Phase 1c; the PM had its own `GasDomains` copy). `set_domain_gas`/`end_domain_gas` with one thread hand out the domains in order 1..Ndomain on every pass, so the C++ loops over them in that order |
| pm_tools.f90 `init_pm` | `init_pm` |
| `particle_mesh` | `particle_mesh` |
| `calc_drift` | `calc_drift` |
| `link_dm` | `link_dm`. Both phases are reproduced: per-cpu (j,k) lists, then relinking into the (i,j,k) chains. The chains end up in increasing `ip` |
| `calc_density_field`, `cic_mass_assignment` | same names |
| `calc_gravity_field` | `calc_gravity_field` |
| `calc_kick` (`accelerate_dm`, `accelerate_gas`, `force`) | `calc_kick`. `force` is the file-local `force()`. Both `accelerate_gas` branches (MF and non-MF) are ported |
| `calc_power_spectrum` | `calc_power_spectrum` (computation) + `write_power_spectrum` (`(7es14.6)` → `%14.6E`) |
| `summarize_dm` | `summarize_dm`, which returns the reductions and prints the summary |
| fft_tools.f90 `fft_3d`, `fft_1d`, `init_fft` | `radray::fft::fft_3d`, `fft_1d` (`src/fft.cpp`, the one FFT layer shared with the hydro's `Fft1d`). `init_fft` is a forward transform of the zero `density`, done inline in the driver |
| cosmo_tools.f90 `calc_pk_3d`; general_tools `sinc` | `calc_pk_3d`; file-local `sinc` |
| cosmo_tools.f90 `ztag` (in `calc_sim_var`) | `ztag(z)` |
| io_tools.f90 `read_init_conditions` | `read_init_conditions` (see Q1) |
| radhydro.f90 `calc_timestep`, complete | `radray::calc_timestep` (`cosmology.hpp`, modes Original/Courant/Fixed), shared by every driver since Phase 1c (it replaced `pm::calc_timestep_full`) |
| radhydro.f90 `radhydro`, gravity part | `apps/pm_driver.cpp` main loop |
| hydro_tools.f90 `summarize_gas` | the hydro port's full `summarize_gas` (Phase 1c; it replaced `summarize_gas_speeds`, see §4.2) |
| hydro_tools.f90 `init_hydro` allocations | `pm_driver` |

Not ported: `write_data_1/2` (snapshot output), the `timing` array and the
`(f8.2,2a10,a)` timing lines, and the multi-thread paths (§5).

## 2. FFT (MKL DFTI → FFTW3f)

`reference/harness_pm/PATCHES.md` §1 has the full mapping (since Phase 1c
the FFTW calls are made by the shared shim
`reference/harness_common/mkl_dfti_shim.f90` under the original
`fft_tools.f90`). In summary:

* **Layout.** `a(n1+2,n2,n3)` in place. The real data is in `a(1:n1,:,:)`.
  The conjugate-even half spectrum `c(n1/2+1,n2,n3)` (MKL
  `DFTI_COMPLEX_COMPLEX` with strides `(1,n1/2+1,(n1/2+1)n2)`) is exactly
  FFTW's in-place r2c layout for the row-major dimensions `(n3,n2,n1)`.
* **Scaling.** Forward is unnormalised, sign −1 (MKL = FFTW). Backward is
  MKL's `DFTI_BACKWARD_SCALE = 1E0/(n1 n2 n3)` (real(4)): FFTW's
  unnormalised c2r is followed by a real(4) multiplication of the n1 real
  outputs of each row.
* **Determinism.** Both codes call `fftwf_plan_dft_{r2c,c2r}_3d` with
  `FFTW_ESTIMATE | FFTW_UNALIGNED` on every call, then
  `fftwf_execute_dft_*`, then `fftwf_destroy_plan`. The plans, and so the
  arithmetic, are identical in Fortran and C++ regardless of the array
  alignment.
* **Versus MKL.** A build of the original with MKL differs at round-off
  level (FFT algorithms), on top of the ifort fast-math differences (B13 in
  `port_notes.md`).
* **Build.** FFTW3f is an optional CMake dependency (`RADRAY_WITH_PM`, ON by
  default). Without `fftw3.h`/`libfftw3f`, only the RT library and driver are
  built.

## 3. Precision and evaluation order (selected)

Everything follows `port_notes.md` §3. The PM-specific points:

* **The `force()` stencil is single precision.** `-2./3*(d(ip1)-d(im1)) +
  1./12*(d(ip2)-d(im2))`: the literals are real(4) and `density` is real(4),
  so the whole expression is computed in float. Only the result is stored in
  the real(8) function value.
* **Density updates.** In `cic_mass_assignment` and the kernel multiply,
  `density = density + w*mass` and `density*wphi` are evaluated in double
  and stored as float. Accumulation order: domains 1..8, cells k,j,i inside
  each domain, chain order inside each cell.
* **Kick reductions.** `vrms += sum(xv_dm(4:6,ip)**2)` sums in real(4) and
  then adds to the real(8) accumulator. `arms += sum(f**2)` is all real(8).
* **`Aphi = -4*pi*a/(6*pi)`.** The unary minus applies to the quotient
  (exact either way).
* **Real `MOD`.** gfortran inlines it with x87 `fprem`, which is exact; it
  matches `std::fmod`. Integer `INT` is `cvttsd2si` (`f::int4`).
* **Drift.** `mod(x + v*dt/2 + ndmr, ndmr)` is evaluated in double and stored
  as float. The stored coordinate can round up to exactly `Ndm`. `link_dm`
  (`mod(int(x*x_dm_gas),Ngas)`) and the CIC (`mod(x*x_dm_gas-0.5+ngasr,ngasr)`)
  both map that back into the box.
* **OpenMP reductions** with one thread equal the sequential loop: the
  private copy starts at the identity and is combined with the initial value
  0 (or ±huge).

## 4. Interpretations

1. **Gravity-only step.** The driver runs `radhydro` with the hydro sweeps,
   the transposes and RT skipped: half clock step, `calc_sim_var`,
   `particle_mesh`, half clock step, `calc_sim_var`, `summarize_dm`, and the
   time-step part of `summarize_gas`. The DM and gas are synchronised at the
   mid-step, as in the original. The kick uses `dt_sim` and the mid-step
   scale factor.
2. **`summarize_gas`.** Until Phase 1c the harness and the C++ used a copy
   of its time-step part (`summarize_gas_speeds`). Both now call the full
   original/ported `summarize_gas`, with its side effects `u_gas(1) =
   max(u_gas(1),Dgas_min)` and `w_gas = max(w_gas,0.)`; the stored outputs
   did not change.
3. **Time step.** `calc_timestep` is ported completely: the istep == 0
   `ra_max` step, `dt = min(2*dt_gas_factor/cmax_gas,
   dt_pm_factor/sqrt(amax_dm))` with the `ra_max` limit, the RT substep
   count, and the in-place unit conversions. The zdata output flags are
   dropped. A problem may fix `da_sim` (`pancake_z_nomf` does).
4. **`calc_power_spectrum`** is called every `pk_every` steps, after the
   step's summaries. That is the same place as in the main loop, where it
   was triggered by the zdata1 output redshifts. It changes no state that is
   used later: `density` and the linked lists are rebuilt by the next
   `particle_mesh`.
5. **Initial conditions** are read from a problem directory instead of
   `Init/` (the harness links `Init/` to it).

## 5. Quirks and bugs in the original (preserved unless noted)

* **Q1. `read_init_conditions` ignores the gas when `Ncpu = 1`.** The files
  are read in `do icpu=1,Ncpu` with `select case (icpu)`: case 1 reads DM,
  case 2 reads gas. With fewer than 2 threads `u_gas.init` is never read, and
  the gas starts with zero density. **Fixed** in the harness
  (`max(Ncpu,2)`, PATCHES.md §2) and in C++, which reads both files.
* **Q2. `vmax_gas` ignores the maximum of vz.** It is
  `max(maxval(x4(2:3)), maxval(abs(x3(2:4))))`, where `x4(2:3)` should
  presumably be `x4(2:4)`. It only enters the log line of `calc_timestep`.
  Preserved.
* **Q3. Globals change units in place.** `calc_timestep` converts `amax_dm`,
  `arms_dm`, `vmax_dm`, `vrms_dm`, `vmax_gas` and `vrms_gas` to km/s and
  km/s/Gyr in place. That works only because `calc_kick` and `summarize_gas`
  recompute them before the next use. If a step skipped the PM call or the
  summary, the next time step would use physical units as if they were grid
  units. Preserved. The stored stats are written after the step, in grid
  units.
* **Q4. `summarize_dm` initialises `vmax = 0` and `xmax = 0`, not
  `-huge`.** If all the velocities of a component are negative, its
  reported maximum is 0 instead of the true maximum. This affects
  diagnostics only. Preserved.
* **Q5. The first step is always `da = a/9`** (`ra_max = 0.1`, real(4)
  literal), whatever the forces. From z = 300 that is a jump to z ≈ 270.
  Preserved; the pancake checks show the leapfrog copes (§6).
* **Q6. `dt_vg = 1/cmax_gas` and `dt_ap = 0.5/sqrt(amax_dm)`** are
  infinite when the gas is static and cold, or when there is no force. If
  both are zero, `calc_da_from_dt` gets an infinite `dt` and produces NaNs.
  This cannot happen with real ICs: `cs > 0` from `Tgas_min`. Preserved.
* **Q7. Gas in a gravity-only run.** In the MF build the kick goes to the
  grid velocity `v_gas`, and `u_gas(2:4)` (momentum relative to the grid) is
  untouched. In the non-MF build the kick goes to `u_gas(2:5)`. In either
  case the gas density moves only through the hydro, so without hydro the
  gas is a static source of gravity. The DM growing mode then obeys
  `D'' = (2/3)(1-f_b) a D`, not `(2/3) a D`. The physics checks use this
  (`tools/ics/pmcosmo.py`, `--growth frozen-gas`). This is the expected
  behaviour of the configuration, not a bug.
* **Q8. `ztag` is not updated for z ≥ 1000** (no branch). The power-spectrum
  file name would then be stale or blank. The C++ returns "" there.
* **Q9. `#endif MF`** has a stray token, which gfortran accepts with a
  warning. **`fft_1d`** is dead code: the hydro calls to it are commented out.
* **Q10. Integer `Ngas/Ndm`.** `N_dm_gas = Ngas/Ndm` is integer division, and
  the particle mass uses `N_dm_gas**3`. If Ngas is not a multiple of Ndm, the
  mass is silently wrong. `PMConfig::make` rejects it, and odd `Ngas`, for
  which `init_domain_gas` finds no subdivision and the FFT padding is wrong.
* **Q11. Thread-count dependence.** With more than one thread:
  * `calc_density_field` deposits from domains handed out dynamically
    (`schedule(dynamic,1)` plus the neighbour exclusion), so the float
    accumulation order in `density` depends on timing;
  * the `calc_kick` reductions are split per thread;
  * `link_dm` merges per-thread lists in thread order.

  Multi-threaded runs of the original are therefore not reproducible bit for
  bit, not even run to run. The port reproduces `Ncpu = 1`. Phase 2 will need
  a deterministic deposit, such as sorted or owner-computes.
* **Q12. `calc_power_spectrum`** re-links the DM and overwrites `density`.
  This is harmless where the main loop calls it. The P(k) bins use
  `kk = 1+nint(|k|/k_f)` and are written for kk = 2..Hgas+1. Corner modes
  with |k| > Hgas k_f are binned but never written.
* **Q13. `summarize_gas` clamps `u_gas(1)` and `w_gas`.** A diagnostic
  routine modifies the state (Interpretation 2). Preserved.

* **Q14. DM speed and acceleration diagnostics in the wrong units.**
  `calc_timestep` converts `vmax_dm` and `vrms_dm` (DM-grid units) with the
  gas-grid `vel_unit`, and `amax_dm` and `arms_dm` (mesh units) with the
  gas-grid `acc_unit`. The printed km/s values are therefore off by
  `x_dm_gas` (a factor 2 with the original's default Ndm = Ngas/2). With
  `pm_refine = r` (`docs/pm_refine.md`) the accelerations are off by r as
  well. Only the log line is affected: the time step uses the grid values,
  which is correct. Preserved.

**`pm_refine` (not in the original).** A PM mesh r times finer than the gas
grid, for Ndm up to r·Ngas; `docs/pm_refine.md`. With r = 1 (the default)
every routine here is the original.

## 6. Initial conditions and physics validation

`tools/ics/make_ics.py` writes `xv_dm.init` and `u_gas.init` in the format
`read_init_conditions` expects (float32 stream):

* DM positions are in DM cells in [0, Ndm), on the lattice q = n + ½.
  Velocities are dx/dτ in DM cells per code time unit, where
  `time_unit = 2/(3 H0 h sqrt(Ωm)) a²`, so dτ = dt/time_unit.
* `u_gas`:
  * `u(1)` = ρ/ρ̄_m (mean Ω_b/Ω_m);
  * `u(2:4)` = ρ v, with v in gas cells per time unit;
  * `u(5)` = ½ρv². The code adds ρ·Tgas_min as thermal energy.
* `w_gas` is not read: it is set neutral. `v_gas` is not read: it is zeroed,
  so the full gas velocity goes in `u(2:4)`.
* Gas cell i is centred at gas coordinate i−½. This matches the PM's
  `x*x_dm_gas − 0.5`.

The linear P(k) uses an Eisenstein–Hu no-wiggle transfer function with n_s,
σ8, Ω's and h from `cosmo_parameters.f90` (the real(4)-rounded values).
Velocities use the growing mode `v = (dD/dτ)/D ψ`. D is integrated in the
code's own units, with radiation in E(a) (from the Meszaros constant mode),
so f(z=300) = 0.94 rather than 1.

Results (C++ = Fortran, so they hold for both):

* **Pancake** (`pancake_x`; shell crossing at a = 0.05; 25 code-chosen
  steps up to D/D_c = 0.93). Errors relative to the analytic amplitude:
  * mode amplitude: x 1.1%, v 1.3%;
  * single particles: 3.2% in x and 5.8% in v. The pattern is a grid-scale
    one, the usual CIC/finite-difference anisotropy. It is the same with
    Ndm = Ngas, so it is not aliasing from the 2:1 particle lattice;
  * time integration alone accounts for ≤ 0.2% in x and ≤ 0.5% in v. This
    was measured by integrating the continuum force law with the code's own
    DKD steps;
  * transverse motion: exactly 0.

  `pancake_z_nomf` (non-MF, fixed da = 0.002, so the first step has
  da/a = 0.6): mode error 1.7% in x and 1.9% in v, of which 1.3% is time
  integration.
* **Linear growth** (`growth_3d`, uniform static gas). The DM P(k) in the
  two lowest bins (k = 1.26 and 2.51 h/Mpc) follows (D/D₁)² to 1.9% up to
  a = 0.02. It then rises above linear theory: 15% by z ≈ 12, through mode
  coupling in a 5 Mpc/h box.
* **Momentum.** The total DM + gas momentum is conserved to |ΔP|/Σ|p| ≤
  3e-8 in all four problems. In `zeldovich_3d` the DM and the static
  gas exchange momentum, equal and opposite.

## 7. Tests

| CTest | what |
|---|---|
| `pm_unit` | FFT layout and scaling against a naive DFT (3D and 1D, non-power-of-2 sizes), round trip; domains; `link_dm` chains; CIC mass conservation; the potential of a single Fourier mode against the FD Green's function |
| `pm_run_<p>` | runs `pm_driver` (fixture) |
| `pm_ref_<p>` | bit-for-bit against `reference/harness_pm/outputs/<p>`: stats every step, SHA-256 of every field every step, last-step `xv_dm`/`phi`, Pk files |
| `pm_momentum_<p>` | total momentum |
| `pm_pancake_x`, `pm_pancake_z_nomf`, `pm_growth_3d` | §6 |

All the PM tests carry the label `pm` (`ctest -L pm`). For the sanitizer
run, configure with `-DRADRAY_SANITIZE=ON` and run `ctest -L pm`: all pass,
with no ASan/UBSan reports. Against a live harness run:

```
reference/harness_pm/run_all.sh
tests/run_pm_all.sh build/pm_driver reference/harness_pm/runs build/pm_runs
python3 tools/compare/compare.py reference/harness_pm/runs/<p> build/pm_runs/<p> --rtol 0 --atol-frac 0
```

## 8. Timings (single thread, -O2)

| run | Fortran `particle_mesh` | C++ `particle_mesh` |
|---|---|---|
| Ngas 32 / Ndm 16, 30 steps (`zeldovich_3d`) | 0.080 s | 0.092 s |
| Ngas 128 / Ndm 64, 10 steps | 1.79 s | 2.25 s |

That is about 3 ms per step at N = 32. A complete N = 32 problem, including
writing every field every step, takes 0.1–0.4 s in C++. The C++ is a
faithful transcription and has not been optimised; for example, the stencil
is recomputed 8 times per particle, as in the original.

## 9. Not faithfully portable / open items

* MKL itself: FFTW replaces it in both codes (§2). The original MKL/ifort
  binary agrees only to round-off.
* Multi-threaded behaviour (Q11): not reproducible in the original either.
* The stored SHA-256 references assume this FFTW build (3.3.10) and an
  x86-64 CPU with the same SIMD dispatch. Elsewhere, regenerate them from the
  harness (`run_all.sh` + `store_outputs.py`).
* `summarize_gas` (full), `write_data_*` and the zdata-driven output schedule
  are left for the hydro/integration phase (1c).
