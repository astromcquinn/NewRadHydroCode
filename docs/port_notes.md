# Port notes: Fortran RadHydro RT → C++ `radray`

This document covers phase 1, the faithful port. It records:

* which C++ function implements which Fortran routine;
* every place where the Fortran was ambiguous and had to be interpreted;
* everything in the original that looks like a bug. Each such behaviour is
  preserved unless stated otherwise, and is flagged here.

**Validation result:** the C++ port agrees **bit for bit** with the gfortran
harness (`reference/harness`, one thread) on all five problems in
`tests/problems`.

* 700/700 field files (all steps) are byte-identical.
* `stats.txt` and `rt_stats.txt` agree exactly for all 10 steps, including
  all 64 `ngam_ray` thresholds.
* The merge path is exercised: `lognormal_16src_smallbuf` merges from step 5.

The CTests `ref_*` check this continuously against
`reference/harness/outputs`.

## 1. Routine map

| Fortran | C++ |
|---|---|
| **physical_constants.f90, cosmo_parameters.f90** | `include/radray/constants.hpp` (`radray::pc`, `radray::cp`, `radray::folded`) |
| **radhydro_global.f90** sizes | `Config` (`config.hpp`). `Ngas=Nrt=N` is set at run time. `Nray_max = nray_max_factor*N^3` |
| radhydro_global.f90 types | `RTGrid` (rt_grid_type, SoA) and `GasGrid` (u_gas, w_gas) in `grid.hpp`. `Ray`, `RayStat` and `RayState` in `rays.hpp`. `RTStat` in `simulation.hpp` |
| radhydro_global.f90 module variables | `Simulation` (`simulation.hpp`), with `SimClock`, `Units`, `RayState`, `Domains`, `RTStat` |
| **radhydro.f90** `init_variables` | `init_variables` (cosmology.cpp) |
| radhydro.f90 `calc_timestep` (RT part) | `calc_timestep` (cosmology.cpp), with a fixed `da_sim` |
| radhydro.f90 `radhydro` (clock part) | `advance_clock_hydro_step` (cosmology.cpp) |
| radhydro.f90 main loop | `apps/rt_driver.cpp` |
| **cosmo_tools.f90** `calc_sim_var` | `calc_sim_var` |
| cosmo_tools.f90 `calc_dt_from_da`, `calc_da_from_dt` | same names |
| **domain_tools.f90** `init_domain_rt` / `set_domain_rt` / `end_domain_rt` | same names (rt_tools.cpp) |
| **halo_tools.f90** `calc_halo` / `add_halo_to_rt` (and the `link_halo` order) | `calc_halo` (rt_tools.cpp) |
| **general_tools.f90** `interpolate` (`'lin'` branch) | `interpolate_lin` (ray_tracing.cpp) |
| general_tools.f90 `index_sort_r4_i4` | `std::stable_sort` in `make_index_gas_rt` (only reached when `N_rt_gas > 1`) |
| **rt_rates.f90** `calc_atomic_rates`, `photoion_sigma` | same names (rt_rates.cpp). `METAL` is off |
| **ray_tracing_tools.f90** `init_ray_tracing` | `init_ray_tracing` + `RayState::init` |
| `ray_tracing_1`, `ray_source`, `ray_cast_split`, `ray_cast_merge` | same names |
| `ray_tracing_3` | same name |
| `ray_cast_1`, `ray_cast_3` (the `#else RAYFAST` versions) | same names |
| `ray_segment` | same name |
| `ray_create`, `ray_remove`, `ray_insert_head`, `ray_insert_tail` | same names |
| `lside_from_radius`, `nside_from_radius`, `nray_from_nside`, `rmin_from_nside`, `rmax_from_nside` | same names |
| `ray_redshifting` / `redshift_ray` (the `#else FREQ1` version) | `ray_redshifting` |
| `ray_buffer_allocate` / `_supplement` / `_deallocate` | same names |
| HEALPix 3.11 `pix_tools`: `pix2vec_nest`, `vec2pix_nest`, `mk_pix2xy`, `mk_xy2pix1` | `healpix::pix2vec_nest`, `healpix::vec2pix_nest` (healpix.cpp) |
| **rt_tools.f90** `init_rt` | `init_rt`. `init_simulation` bundles the init sequence |
| `radiative_transfer`, `prep_rt` (`calc_rt_i`), `post_rt` (`calc_rt_f`, `calc_ray_f`) | same names. The shared gas-averaging code is in `average_gas`, `gas_state_from_avg` and `store_gas_state` |
| `rt_step_1` (`calc_rt_1`, `calc_ray_1`), `rt_step_2`, `rt_step_3` | same names |
| `radiative_physics` (`photo_source`, `photo_cell`, `make_index_gas_rt`) | `radiative_physics`. `photo_source` and `photo_cell` are inlined and marked with comments |
| `summarize_rt` | `summarize_rt` (to stdout) |
| unit 201 `ray.log` | `<outdir>/ray.log` (same content, not byte-identical formatting) |

### Not ported

These do not affect the RT state or the spec outputs:

* `make_data_2d` and `write_data_3` (2-D slices);
* `calc_gas_3` (the `zreion` reionisation redshift);
* the `Log/nside_of_r.dat` and `r_of_nside.dat` dumps;
* `random_seed` (the random branch is `if (.false.)`);
* `mp_unit` and `Msun_unit` (particle units);
* the `timing` bookkeeping.

## 2. Index conventions

* RT cells are addressed with Fortran 1-based `(irt,jrt,krt)` everywhere.
  `RTGrid::cell(i,j,k)` maps them to the 0-based linear index
  `(i-1)+N*((j-1)+N*(k-1))`, with i fastest as in Fortran.
* Multi-component fields keep the component index fastest:
  `ngam(f,cell) -> ngam[cell*5+f-1]`. This matches the Fortran derived type
  and the output format, so the output is a plain memory dump.
* Ray storage indices `ir` are 1-based. Slot 0 is a dummy, and `0` means
  "no ray" in the chains, as in Fortran. `htoc_ray(m,cell)` with m = 1..4 is
  `RTGrid::htoc_ray(m, c)`.
* Loop counters over photon groups, components and segments are 0-based in
  C++. Wherever Fortran used a literal index (for example `ngam_ray(16)`,
  `alpha_rec(4,iT)`), the C++ keeps the Fortran number and converts it
  explicitly (`ngam_ray[16-1]`, accessors that take 1-based arguments).

## 3. Precision and Fortran semantics (how faithfulness was achieved)

* `real(4)` maps to `float` and `real(8)` to `double`. Every Fortran
  expression is reproduced with the same operation order and the same
  promotion points. Examples:
  * `ray%r**2` (merge) and `ray%ngam*ray%nu` (merge) are float products.
  * `sum(u(2:4)**2)/u(1)/2` in `photo_cell` is computed entirely in float.
  * `1-rt_grid%fgam` in `ray_cast_3` is a float subtraction.
  * `(/irt,jrt,krt/)-0.5` is float.
  * `sum()` of a `real(4)` array accumulates in float, left to right, from 0.
* **Single-precision literals.** Literals such as `0.99`, `1E-6`, `1E-3`,
  `0.98` and `1E-9` are single precision in Fortran. They are written as
  `0.99f` and so on and promoted exactly where Fortran promotes them.
  Several `real(8), parameter`s are initialised from single literals:
  * `omegam = 0.27` is really `(double)0.27f`, and so are `omegab`,
    `hubble0`, `XHy`, `YHe`, `freq_init`, `fResc`, `eSF`, `eSN` and
    `lgT_del`;
  * `HI_eV = 13.60` (and the HeI and HeII equivalents);
  * `sfr_min = ngam_min = 1E-6`, `fgam_min = seg_min = 1E-4`;
  * `gamma_ideal = 5./3` is a float division.

  Exponent literals in `rt_rates.f90` (`x**1.503`) are float literals
  promoted to double.
* **Constants folded at compile time.** gfortran folds some expressions with
  correctly rounded MPFR:
  * `10**25.2` in `pop2_init` (the whole expression is real(4));
  * `Tcmb0**4` in `omegar`;
  * `log10(1.01)` (real(4)) in `rt_step_1`;
  * `log(4.)` (real(4)) in `lside_from_radius`;
  * `1E3*(ob*h0**2)**0.4` in `calc_sim_var`;
  * `1D1**lgD_min`.

  Correct rounding can differ from libm by one ulp, so these are hardcoded as
  hex-float literals taken from a gfortran dump
  (`tests/fortran_consts/dump_consts.f90`, output in `expected_consts.txt`).
  All the other derived constants are computed in C++ and checked against
  the same dump by the `fortran_semantics` test.
* **Integer powers** `x**2`, `x**3`, `x**4` are expanded the way gfortran
  does it: `x*x`, `(x*x)*x`, `(x*x)*(x*x)`. Real powers use `pow`.
* **`INT`, `NINT`, `FLOOR`, `CEILING`, `MOD`** are implemented in
  `fortran.hpp`:
  * `INT` truncates toward zero. NaN and out-of-range values give
    `INT32_MIN`, as gfortran's `cvttsd2si` does. This matters for `log10(0)`
    in `calc_ray_1`.
  * `NINT` rounds half away from zero.
  * `MOD` for integers is C++ `%`. For reals it is `fmod`.
* **libm.** The C++ calls the same glibc `exp`, `expf`, `log`, `log10`,
  `pow`, `sqrt`, `atan2`, `sin` and `cos` that the gfortran build calls. The
  build uses `-ffp-contract=off -fno-fast-math`, so no FMA and no
  reassociation.
* **Linked lists and processing order are reproduced exactly:**
  * head insertion of source rays;
  * tail insertion of rays that move cells;
  * the `htoc_ray(1:2)` snapshot versus the live `htoc_ray(3:4)` (rays that
    are appended during a pass are not processed again in that pass);
  * the per-cone `merge1` lists in chain order;
  * the domain order. With `Ncpu=1`, `init_domain_rt` gives `n1 = 2`
    (8 domains), and `set_domain_rt` returns them in the order 1..8 every
    pass. Inside each domain the loops run k, j, i.
* **Ray storage pool.** It is emulated exactly: the `ray_index` permutation,
  `Nray`, `Nrbuf`, `ifree` with LIFO reuse, the `nrbuf/4` refill in
  `ray_create` and spill in `ray_remove`, and the
  allocate/supplement/deallocate logic with `Ncpu=1`. This is required, not
  cosmetic; see bug B1. It also makes the storage indices identical to the
  Fortran ones. Only the touched part of the pool is materialised: capacity
  for `Nray_max` rays is reserved but untouched. Positions above the
  high-water mark are implicitly the identity and zero, which is equivalent.
* **HEALPix.** It is a transcription of the HEALPix **3.11 Fortran**
  `pix2vec_nest`/`vec2pix_nest` (32-bit path, `nside <= 8192`), which is the
  version the original Makefile links. It is not a transcription of
  healpix_cxx/healpy, whose floating-point expressions differ:
  * healpy's `pix2vec` differs by up to 1.1e-15;
  * healpy's `vec2pix` can disagree right at pixel boundaries.

  Against healpy, 1e6 random vectors (nside 1..8192, plus special
  directions) give identical pixels.

## 4. Interpretations of ambiguous Fortran

1. **Time stepping in the driver.** For each step the driver does what the
   `radhydro.f90` main loop does without hydro or PM:
   * `calc_sim_var`;
   * `calc_timestep`. Both branches reduce to `da_sim` fixed, then
     `calc_dt_from_da` and `dt_sec = dt_yr*yr2sec`, then the RT step part;
   * the two half-step clock updates, each followed by `calc_sim_var`;
   * `calc_halo`, then `radiative_transfer`, which rewinds the clock by
     `da_sim` and re-advances it by `da_sim/nstep_rt` per substep;
   * `summarize_rt`.

   The initial state is `z = z_init` and `a = 1/(1+z_init)`, with `year`
   from the `init_variables` formula. This agrees with the harness.
2. **Sources.** Each line of `sources.txt` is one `source_type`, with `L` and
   `sfr` rounded to `real(4)`. `link_halo` reverses the per-cell order twice,
   so each cell's sources are summed in input order, in double, and stored
   as float.
3. **`gfortran MAX/MIN` with mixed kinds**, for example
   `max(u(1),Dgas_min)`. Both arguments are promoted to double, and the
   result is stored back to float where the left-hand side is real(4).
4. **`nray_active`** in `stats.txt` is `Nray - ray_save%nfree` after
   `post_rt`, which is `Nray`, because the free buffers have been returned.
   **`ngam_rays_total`** is summed in storage order over slots with
   `nside > 0` after `post_rt`. At that point dead slots have been zeroed,
   so these are the live rays. **`ncst1`, `nmrg1`, `nspl1` and `ndel1`** are
   the values of the last RT substep (`ray_save`).
5. **`int(Nray_max/64.)` and `j*int(...)`** in `rt_step_1` are int32 in
   Fortran. The C++ uses int64 for the product. This is identical unless
   `64*N^3 > 2^31`, roughly N ≥ 323, where the Fortran overflows (bug B9).

## 5. Bugs and quirks in the original (preserved unless noted)

* **B1. Dead rays steer the split/merge thresholds.** `ray_remove` only
  unlinks a ray and pushes its index onto the free list. `calc_ray_1` scans
  every storage slot with `nside > 0`. Rays that were freed during the
  current hydro step are therefore counted in `ngam_pdf`, and so in
  `ngam_ray`, until `ray_buffer_deallocate` zeroes them in `post_rt`.
  Preserved by emulating the pool exactly. `ngam_avg/Nray` divides by the
  pool count, but it is only logged.
* **B2. `ngam_ray(j)` is never reset** while the ray count is too small to
  reach `j*Nray_max/64`. It keeps its stale value, or its initial 0. The
  consequences:
  * with the default `Nray_max = 64 N^3`, merging is effectively disabled;
  * with `ngam_ray(32) = 0`, every ray with `dr < 0` is split.

  Preserved.
* **B3. `rt_grid%dx` is never reset.** `calc_rt_1` zeroes `nu`, `ngam`,
  `fgam`, `gamma` and `lambda`, but not `dx`. The new photon-weighted
  directions are added on top of the previous substep's unit vector.
  `radiative_physics` then normalises the sum, and that direction is used
  by `make_index_gas_rt` when `N_rt_gas > 1`. Preserved.
* **B4. `tau` is zero in neutral cells during the first RT call.**
  `init_rt` zeroes `tau`. `calc_rt_1` updates it only in cells with
  `fHI < 0.99` or `fHeI < 0.99`, and `calc_rt_f` (all cells) runs only in
  `post_rt`. During the first hydro step:
  * rays are attenuated only through `fgam`;
  * `ray_cast_1` spreads photons over segments in proportion to path length;
  * a source in a neutral cell always emits rays, instead of depositing
    photons locally.

  Preserved.
* **B5. Out-of-bounds histogram in `post_rt`/`calc_ray_f`.** The arrays are
  `ngam_pdf, nrad_pdf(-9:10)`, but `max(ngam,1E-9)` with a real(4) `1E-9`
  (9.99999972e-10) gives bin −10 for empty cells. This is a stack write out
  of bounds on every call, and the C++ port's UBSan/ASan build caught it.
  **Fixed** as in the harness (PATCHES.md 5): the arrays run −10..10. It
  affects diagnostics only.
* **B6. `calc_gas_3` bug.** The j-loop sets `irt` instead of `jrt`, so
  `jrt` is uninitialised in `zreion(2:3)`. It is not ported (zreion only);
  the harness fixed it.
* **B7. `merge1` overflow.** Each of the 192 cones holds at most 10000 ray
  indices, with no bound check. That is also about 15 MB of automatic stack
  per thread. The C++ uses dynamic vectors and throws if a cone exceeds
  10000, which would have been memory corruption in Fortran. Similarly,
  `ray_create` and `ray_buffer_supplement` throw where the Fortran would
  overflow `ifree`, or index it negatively when `nrbuf/4 < nr`.
* **B8. At most 3 segments per cast and advance.** A unit step along a
  near-diagonal can cross 3 faces and need 4 segments. The remainder is
  dropped, and the ray advances less than 1 cell. Preserved; the `rays` test
  confirms it happens.
* **B9. int32 overflow** of `j*int(Nray_max/64.)` for N ≥ ~323, and loss of
  float precision in `int(Nray_max/64.)` for very large `Nray_max`. These are
  irrelevant at N = 32. The C++ uses int64 here (see §4.5).
* **B10. The `1E-6` nudge is lost.** In `ray_cast_3`, `x0+(l+1E-6)*dx` is
  stored as float. Near x = 32 the float ulp is 3.8e-6, so rays usually end
  exactly on a cell face. The home cell is still `int(x)`. A stored `x` can
  also round up to exactly `N`, which maps to cell 1 through
  `mod(int(x),N)`. Preserved by storing float32 and computing in double.
* **B11. `dr` is redundant but stored.** `dr` and `rmax-r` round
  differently in float. It is kept as a field.
* **B12. Physics oddities.** These are preserved; they matter for later
  physics validation, not for the port:
  * the case-B recombination *cooling* coefficients (`coef_cool(4:6)`) use
    the case-A formulas;
  * HeI collisional excitation cooling (`coef_cool(11)`) is zero;
  * `coef_cool(15)` is never set;
  * recombination radiation is disabled (`if (.false.)`);
  * both branches of the `nside_merge` choice give 2;
  * `photo_cell`'s final gas update does not clamp `us(1)` to `Dgas_min`.
* **B13. ifort versus gfortran.** The original was built with
  `ifort -O3 -ipo -ftz`. ifort's default `-fp-model fast=1` allows
  reassociation and FMA, and `-ftz` flushes denormals to zero. Production
  RadHydro results are therefore not bitwise reproducible with IEEE
  gfortran or C++. Our reference is the IEEE gfortran harness.
* **B14. `end_domain_rt`** updates the shared domain table outside the
  critical section. This is a race with more than one thread, and
  irrelevant here.
* **B15. `ray_buffer_supplement`** contains `!$omp barrier` inside a
  `parallel do`, which is non-conforming OpenMP. The harness removed it, and
  the C++ is serial.

## 6. Open issues / limitations

* Only `N_rt_gas = 1` (Ngas = Nrt) is validated. The general path in
  `make_index_gas_rt` uses a stable sort, whereas the original uses ifort's
  `qsort`. Ties could be ordered differently.
* `Ncpu = 1` only. The buffer logic is written for one thread; the
  multi-thread `nr(:)` exchange in `ray_buffer_supplement` is reduced to the
  one-thread case.
* `ray.log` has the same content as the Fortran unit-201 log, with
  different number formatting. It is not compared.
* Bitwise agreement assumes x86-64, glibc libm and the compiler flags
  above. A different libm (or a GPU) changes `exp`, `log`, `pow` and so on
  by an ulp, and the results will then only agree to round-off.
