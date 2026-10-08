# Halos and ionising sources

The original RadHydro does not find halos itself. Each hydro step,
`read_halo` (io_tools.f90) reads an external catalog `Halo/halo.z=XX.XXXX`
made by a separate halo finder, which we do not have. It rescales the
catalog to a global star-formation history (`Halo/sfr.txt`), and
`calc_halo` / `link_halo` / `add_halo_to_rt` (halo_tools.f90) put the
sources on the RT grid. This module adds:

* an **on-the-fly friends-of-friends halo finder** for the DM particles,
  with optional spherical-overdensity (SO) refinement
  (`include/radray/halos.hpp`, `src/halos.cpp`);
* a **source model** that assigns SFR, UV luminosity, stellar mass and gas
  fraction to halos (`include/radray/halo_sources.hpp`,
  `src/halo_sources.cpp`);
* **legacy compatibility**: a reader and writer for the catalog format, and
  `read_halo` with its `sfr.txt` scan, redshift tag and rescaling,
  bit-identical to the Fortran;
* ports of `link_halo` and of the `add_halo_to_rt` metal injection, built on
  the existing `calc_halo` in `src/rt_tools.cpp`;
* analytic Press–Schechter and Sheth–Tormen mass functions, used for
  validation (`include/radray/mass_function.hpp`);
* `apps/find_halos`: runs the finder and the source model on an `xv_dm`
  snapshot and writes a legacy catalog.

The code is built as the `radray_halos` library. OpenMP is optional
(`-DRADRAY_HALOS_OPENMP=ON`), and the results are bit-identical for any
number of threads.

## 1. Where it plugs into the step

radhydro.f90 `radhydro` does, per hydro step:

```
forward hydro (3 sweeps)
year, a += dt/2 ; calc_sim_var
particle_mesh                       ! DM + gas gravity, drift/kick
reverse hydro (3 sweeps)
year, a += dt/2 ; calc_sim_var
if (idata1 > 1) then                ! after the first output time
   call read_halo                   ! -> source(:), Nhalo, Nsource
   if (nsource > 0) then
      call calc_halo                ! link_halo + add_halo_to_rt
      call radiative_transfer
   endif
endif
```

The C++ equivalent of the block inside `if (idata1 > 1)` is

```cpp
SourceCatalog cat = update_sources(sim, xv_dm, Npart, params);  // read_halo|FoF, calc_halo
if (cat.Nsource > 0) radiative_transfer(sim);
```

`update_sources` runs once per hydro step, at the end-of-step redshift
`sim.clk.redshift`, which is the redshift `read_halo` uses. It works in one
of two modes:

* `Mode::Legacy`: runs `read_halo(dir_halo, redshift)`;
* `Mode::OnTheFly`: runs `find_halos` on `xv_dm`, then `assign_sources`. If
  `write_catalog_dir` is set, it also writes each catalog in legacy format,
  so that a run can be replayed later in legacy mode.

Both modes then call `calc_halo_records`, which runs `link_halo`, the
existing `calc_halo` and, if `add_metals` is set (the original always does
this), the `w_gas(4)` metal injection.

In on-the-fly mode the finder sees `xv_dm` as it is when called: after the
full PM step and at the end-of-step redshift. A legacy catalog was made from
a snapshot at the catalog redshift `z2 >= z`. The `idata1 > 1` gate belongs
to the driver; `update_sources` does not apply it.

`Simulation` has no particle array yet, because the PM port is in progress.
`update_sources` therefore takes the particles explicitly.

### Particle units (correction to the task description)

`xv_dm(1:3,:)` is in **DM-grid units [0, Ndm)**, not gas-grid units. There
are three pieces of evidence: pm_tools.f90 `calc_drift` wraps with
`mod(..., ndmr)`, the mass assignment uses `xv_dm*x_dm_gas`, and
`link_dm` uses `int(xv_dm*x_dm_gas)`. Velocities are in DM-grid units per
`time_unit`. `summarize_dm` converts them to km/s with
`x_dm_gas*vel_unit/1E5`. `make_halo_units(cfg, clk, units)` encodes these
conversions:

* `L = Ndm`;
* `pos_to_mpch = box/Ndm`;
* `vel_to_kms = x_dm_gas*vel_unit/1e5`;
* `mpart_msun = Msun_unit`.

The FoF routine itself only needs the box length in the units of the
positions.

## 2. Friends-of-friends

* **Linking length:** `ll = b * L / Npart^(1/3)`, with `b = 0.2` by default.
* **Friends:** two particles are friends if their minimum-image separation
  satisfies `d^2 <= ll^2`. The separation is computed in double precision,
  per axis, with `dx -= L*nearbyint(dx/L)`.
* **Groups:** the connected components of the friendship graph. Groups with
  at least `min_npart` members (default 20) are halos.
* **Cell grid:** the grid has `nc = min(floor(L/ll), 2*Npart^(1/3), 1024)`
  cells per side, so each cell is at least `ll` wide. Particles are put in
  the cells with a counting sort, which is stable in particle index. Each
  cell is compared with itself and its 26 periodic neighbours. Only cell
  pairs with `c2 >= c` are visited, so each pair of cells is visited once.
  Neighbour lists are deduplicated, which keeps very small grids
  (`nc < 3`) correct.
* **Union-find:** the union-find is lock-free. It uses CAS on a parent
  array with path halving, and always links the larger root under the
  smaller one. Parent pointers only decrease, so the final root of every
  component is its smallest particle index, whatever the order of the
  unions. This is what makes the result independent of the thread count.
  The parallel loop over cells (`schedule(dynamic)`) only performs unions.
* **Output:** each group is labelled by its smallest particle index. Groups
  are listed in increasing label order, and members in increasing particle
  index. All per-halo sums run over members in that order.
* **Cost:** O(N) for a fixed clustering. Pair work grows like
  N × (density inside halos) × ll^3. For 128³ clustered particles the
  finder takes 2.0 s serial and 1.4 s with 4 threads. Grid building and
  group collection are still serial. If that matters, they are the next
  thing to parallelise with per-thread prefix sums.

### Halo properties

| quantity | definition |
|---|---|
| `N` | FoF member count, or particles inside `R_SO` when SO is on |
| `M` | `N * Msun_unit` [Msun]. `Msun_unit = rhoc*Ωm*(box Mpc)^3/h/Npart/Msun` (cosmo_tools.f90), with `Npart = Ndm^3`. The DM particles carry **all** the matter (Ωm, not Ωm−Ωb), as in the original PM. |
| `x` | comoving Mpc/h in `[0, box)`. **Default: centre of mass**, periodic-safe: members are unwrapped by minimum image around the lowest-index member. This is valid while the group extent is below L/2, which is always true for real halos. Alternative `HaloCentre::PotentialMinimum`: the member with the lowest softened potential `−Σ 1/sqrt(r²+ε²)`, with `ε = 0.1 ll`. The sum is exact up to `potmin_nmax = 5000` members; above that it runs over every k-th member. |
| `v` | mean member velocity, peculiar km/s |
| `R` | virial radius, comoving Mpc/h: `M = 4π/3 Δc ρcrit(z) R_phys³`, with **Bryan & Norman (1998)** `Δc = 18π² + 82x − 39x²`, `x = Ωm(z) − 1`, relative to the critical density. With SO, R is the SO radius instead. |
| `T` | virial temperature `T = μ m_H G M / (2 k R_phys)`, with `μ = 0.59` (ionised primordial gas), in K |

### SO refinement (`so_refine = true`)

1. Start at the FoF centre.
2. Gather all particles, not only the FoF members, within a search radius.
   The search radius is the larger of the FoF extent and `2 R_vir(M_FoF)`,
   and it grows until the crossing lies inside it, up to 0.45 L.
3. Sort the particles by (distance, index). `R_SO` is the **outermost**
   radius at which the enclosed mean density is still at least
   `Δ ρcrit(z)`. Δ is Bryan & Norman by default, or `so_delta_crit`.
4. Recentre on the centre of mass inside `R_SO`, and repeat
   `so_iterations` times (default 3).
5. Halos with fewer than `min_npart` particles inside `R_SO` are dropped.

SO spheres of neighbouring FoF groups may overlap, so mass may be
double-counted. There is no unbinding.

### Sorting

Halos are sorted by descending `M`, then descending `N`, then ascending
`id` (the FoF label). This makes the order fully deterministic.

### MPI (not implemented; plan)

* **Partition and ghosts:** slab or block partition, as in
  `domain_tools.f90`. Each rank imports ghost particles within `ll` of its
  faces. Each rank then runs the same cell-grid FoF on owned plus ghost
  particles, with global particle IDs as union-find keys.
* **Merging groups across boundaries:**
  1. For every ghost particle, each rank emits a `(global id, local root
     global id)` pair.
  2. The pairs are exchanged with the owning ranks.
  3. A distributed union-find (or iterated min-label propagation) over the
     cross-boundary edges merges groups until no label changes.

  The label is the smallest global ID, so the result equals the serial
  partition.
* **Properties:** member properties are reduced to the rank that owns the
  group's minimum-ID particle. All reductions use the member order sorted
  by global ID, which keeps the output deterministic. Groups larger than a
  subdomain need the iterative merge, but no special casing.

## 3. Source model

### What the original code tells us

These facts are certain, because they come from the code:

* `source_type`: `x(3), v(3)` real(4), `N` int(4), and `M, R, T, L, sfr,
  mstar, fgas` real(4). This is 56 bytes. Only `x`, `L` and `sfr` are used
  (in `link_halo` and `add_halo_to_rt`).
* `x` is in comoving Mpc/h, because `link_halo` uses `x/x_rt_box` with
  `x_rt_box = box/Nrt`.
* `sfr` is in Msun/yr. It drives the metal injection
  `dmstar = sfr*dt_yr*Msun`, and the comment in rt_tools.f90 says so.
* `L` is a UV luminosity in **erg/s/Hz**. The photons emitted per RT substep
  are `ngam = fResc*pop2_init*L*dtsec_rt`, with
  `pop2_init = 10**25.2 × (0.4185, …, 0.0054)`. These weights sum to 1.0001,
  so `pop2_init` is about 1.6e25 ionising photons per erg/Hz of UV energy,
  split over the 5 groups. The `ergs/s/Hz` comment in halo_tools.f90 fits
  this.
* The catalog is sorted by descending mass, and the sources are the
  **leading** records with `L > 0`: `read_halo` stops at the first `L <= 0`.
  So the source property is monotonic in mass. In practice, every halo above
  a mass or temperature threshold is a source.
* `read_halo` multiplies `L` and `sfr` by `sfr(z)/sfr(z_cat)`, interpolating
  column 10 of `sfr.txt` linearly in z. So the catalogs' `L` and `sfr` follow
  a **global SFR history**: halos are fixed between catalog redshifts, and
  their total tracks `sfr.txt`. `calc_halo` divides `ΣL` and `Σsfr` by
  `(box/h)^3` and calls the result "Comoving SFR in Msun/yr/Mpc^3". This
  makes column 10 of `sfr.txt` a comoving SFR density in those units.
* `eSF = 0.05` and `rho_source = 50` are used only in the SN heating of
  rt_tools.f90: `sn_rate = x*eSN*eSF*rhog/tdyn*c²`, with
  `tdyn = sqrt(3π/(32 G ρ))`. They are not used in the halo sources.

### What the papers say (from memory and abstracts; full texts were unreachable)

The full texts were not reachable from this environment. arXiv, IOP, ADS and
Harvard DASH were blocked by the network proxy. What follows comes from
abstracts, secondary citations and recollection.

* **Trac & Cen (2007, ApJ 671, 1; astro-ph/0612406):**
  * Hybrid code: N-body for the dark matter, "physically motivated
    prescriptions for baryons and star formation", and adaptive ray tracing.
  * Halos are resolved down to virial temperatures of 10^4 K, the atomic
    cooling limit.
  * The acknowledgements thank L. Hernquist for the code that computes the
    **semi-analytical SFR**. This is the Hernquist & Springel (2003, MNRAS
    341, 1253) model. The ten-column `sfr.txt`, with the SFR in the last
    column, is consistent with being the output of that code.
  * The stellar IMF depends on gas metallicity (Pop III/Pop II). Pop II only
    is what `pop2_init` represents in this code version.
  * My recollection, which I **could not verify**, is that the halo SFR was
    taken as proportional to halo mass above the cooling threshold, and
    normalised to the Hernquist & Springel global SFR history. That matches
    the code's `sfr.txt` rescaling exactly.
* **Trac, Cen & Loeb (2008, ApJ 689, L81):** the IGM temperature after
  inhomogeneous reionization, with early (z≈9) and late (z≈6) models. The
  abstract does not describe the source prescription. The two reionization
  histories were presumably obtained by changing the source efficiency or
  the SFR normalisation. In this code that would mean a different `sfr.txt`
  or `fResc`.

### Implemented models

**Default: `SfrModel::GlobalMassProportional`**, a reconstruction of TC07:

```
active halos: T_vir >= Tvir_min (1e4 K) and M >= M_min (0)
sfr_i   = M_i / Σ_active M_j  ×  ρ̇*(z) × (box/h)³          [Msun/yr]
L_i     = L_per_sfr × sfr_i,  L_per_sfr = 8.0e27 erg/s/Hz per Msun/yr
mstar_i = min(M_i / Σ M_j × ρ*(z) × (box/h)³, fb M_i)         [Msun]
fgas_i  = max(fb − mstar_i/M_i, 0)                             (our definition)
```

* **Global SFR density ρ̇*(z):** if `sfr_history_path` is set, column 10 of
  that `sfr.txt`, linear in z, and 0 outside the table. Otherwise the
  Hernquist & Springel (2003) fit
  `ρ̇* = ρ̇m χ²/(1+α(χ−1)³ exp(βχ^{7/4}))`, `χ = (H/H0)^{2/3}`, with
  `ρ̇m = 0.013 Msun/yr/Mpc³`, `α = 0.012`, `β = 0.041`. These constants are
  quoted from memory. With them the fit peaks at 0.14 Msun/yr/Mpc³ at
  z ≈ 5.6, which is the expected shape.
* **Stellar density ρ*(z):** `∫ ρ̇* dt` from `z_max = 50`.
* **Consistency with the original:** the total SFR of the box matches what
  `calc_halo` reports as `rt_save%sfr`.
* **Source of `L_per_sfr`:** the 1500 Å Salpeter conversion of Madau,
  Pozzetti & Dickinson (1998), `L_ν = 8.0e27 SFR`. With `pop2_init` it gives
  1.3e53 ionising photons/s per Msun/yr before `fResc`. That is a standard
  Pop II (Starburst99, Salpeter, low Z) value, which strongly suggests this
  is the convention the catalogs used.

**Alternative: `SfrModel::Dynamical`**

```
sfr = eSF × fb × M / t_dyn,   t_dyn = sqrt(3π / (32 G ρ_vir)),  ρ_vir = Δc ρcrit(z)
mstar = min(sfr × t(z), fb M)      (t = cosmic time since z_max; crude)
```

This is the same form as the star formation rate in the original SN-heating
term, with `eSF = 0.05` from cosmo_parameters.f90. At z = 8, `t_dyn` is about
1.2e8 yr, and a 1e9 Msun halo has an SFR of about 0.07 Msun/yr.

**Record order:** star-forming halos come first, in mass order, followed by
the rest. Records whose `L` underflows real(4) count as non-sources. This
keeps `read_halo`'s leading-sources contract even when SO radii make `T_vir`
non-monotonic in M.

### Confidence and uncertainties

| item | status |
|---|---|
| units of `x`, `sfr`, `L`; leading-sources contract; `sfr.txt` rescaling | **certain** (from the code) |
| `sfr.txt` column 10 = comoving SFR density from the Hernquist & Springel semi-analytic code | likely (acknowledgements + calc_halo comment) |
| SFR ∝ M_halo above T_vir = 1e4 K, normalised to the global SFR history | **reconstruction**; consistent with every code-level constraint, but not verified against the paper text |
| `L = 8e27 × sfr` | **inferred** from the `pop2_init` normalisation; the catalogs could have used a different per-halo conversion (e.g. age-dependent) |
| Hernquist & Springel fit constants (0.013, 0.012, 0.041) | from memory; check against HS03 eq. 46–47 |
| `v`, `R`, `T`, `mstar`, `fgas` definitions | ours; the original never reads them |
| mass definition of the original external finder (FoF b=0.2? SO?) | unknown |

## 4. Legacy compatibility

### Files

* **`Halo/halo.<tag>`:**
  * `read(un) Nhalo` with `Nhalo` integer(8), then `read(un) source(ih)`
    records.
  * The file is opened with Intel `form='binary'`, which means **no record
    markers**: it is a raw stream, the same as gfortran
    `access='stream', form='unformatted'`. This is `CatalogFormat::Stream`,
    the default for writing.
  * `CatalogFormat::Sequential` (4-byte markers, as gfortran
    `form='unformatted'` writes) is also supported.
  * The reader auto-detects the format from the file size.
  * Native little-endian.
* **`Halo/sfr.txt`:** `read(un,*) z,dummy(8),sfr`. The reader reproduces
  list-directed semantics:
  * blanks and commas separate values;
  * a row may continue on the next lines;
  * values left on a row's last line are discarded;
  * `d` exponents are accepted;
  * a partial row at end of file is dropped (`end=1`).

### `read_halo` behaviour, reproduced exactly (including quirks)

* **Scan:** scan rows with `z2 ← z1, z1 ← z`, and stop at the first row with
  `z <= redshift`. The catalog is at `z2`, the row **before** it. Even an
  exact match `z1 == redshift` uses the previous, higher-z catalog, with
  `sfr_scale = sfr1/sfr2`.
* **No catalog:** if `z2 <= 0` or `sfr2 <= 0` (for example, `redshift` is
  above the first row), then `Nhalo = Nsource = 0` and there is no RT that
  step.
* **End of file:** if the file ends before any row has `z <= redshift`, the
  last two rows are used and the SFR is linearly **extrapolated**.
* **`sfr_scale`:** `sfr_scale = (sfr2 + (sfr1−sfr2)(redshift−z2)/(z1−z2))/sfr2`
  in real(8). `L` and `sfr` of the leading `L > 0` records are multiplied by
  it and rounded to real(4).
* **Tag:** `z < 10` gives `'z=0'` followed by `f6.4`, so `z=08.0000`.
  `10 <= z < 100` gives `'z='` followed by `f7.4`.
  * gfortran rounds exact decimal ties to even, like printf. This is
    verified against gfortran: for example, 8.03125 gives `z=08.0312`.
  * A value of 9.99995 or above but below 10 overflows `f6.4` and gives
    `z=0******`.
  * For `z >= 100` the original leaves the tag undefined, and we throw.
  * Intel's tie rounding is assumed to be the same.
* **Path length:** the original's `character(100) :: fn` truncates long
  paths. Our port does not.
* **Array bound:** the original holds at most `Nhalo_max = Npart/1024`
  records. Neither the on-the-fly finder nor our reader applies that bound.

### `link_halo` / `calc_halo` / `add_halo_to_rt`

* **Cell of a source:** `link_halo(records, Nsource, cfg)` returns the
  1-based RT cell `1 + mod(int(x/x_rt_box), Nrt)` of each source.
* **Order within a cell:** `link_halo`'s two passes reverse each chain
  twice, so a cell's chain visits its sources in record order. This is the
  order the existing `calc_halo` (rt_tools.cpp) sums in. `calc_halo_records`
  passes the list straight to it; nothing is duplicated.
* **Metals:** `add_halo_metals` ports the `w_gas(4)` injection that
  `calc_halo` in rt_tools.cpp skips:
  * `w4 += fZesc2*Zyield2*sfr*dt_yr*Msun/vol_gas/rho_unit`, with
    `fZesc2 = 0.10` and `Zyield2 = 0.02` as real(4) literals;
  * at gas cell `1 + mod(int(x/x_gas_box), Ngas)`;
  * in the original visiting order: RT cells with k slowest, then record
    order within a cell;
  * the sum is done in real(8) and stored to real(4).

## 5. Parameters

| parameter | default | meaning |
|---|---|---|
| `HaloFinderParams::b` | 0.2 | linking length / mean interparticle spacing |
| `min_npart` | 20 | minimum particles per halo |
| `centre` | `CentreOfMass` | or `PotentialMinimum` |
| `potmin_nmax` | 5000 | exact potential up to this many members |
| `so_refine` | false | SO refinement |
| `so_delta_crit` | 0 | SO threshold relative to ρcrit (≤0: Bryan & Norman) |
| `so_iterations` | 3 | SO re-centring iterations |
| `mu` | 0.59 | mean molecular weight for T_vir |
| `ncell_max` | 0 (auto) | cap on FoF cells per dimension |
| `SourceModelParams::model` | `GlobalMassProportional` | or `Dynamical` |
| `Tvir_min` | 1e4 K | star formation threshold |
| `M_min` | 0 | extra mass threshold |
| `sfr_history_path` | "" | `sfr.txt` for the global SFR (else HS03 fit) |
| `hs_rho_m, hs_alpha, hs_beta` | 0.013, 0.012, 0.041 | HS03 fit |
| `L_per_sfr` | 8.0e27 | erg/s/Hz per Msun/yr |
| `eSF` | 0.05 | Dynamical model efficiency |
| `z_max` | 50 | start of the stellar-mass integration |
| `HaloSourceParams::mode` | `OnTheFly` | or `Legacy` |
| `dir_halo` | "Halo/" | legacy input directory (with trailing `/`) |
| `add_metals` | true | metal injection as in `add_halo_to_rt` |
| `write_catalog_dir` | "" | dump on-the-fly catalogs in legacy format |

## 6. Tests

These run in CTest as `halos` and `halo_catalog`.

* **FoF against brute force:** the partition is identical (labels, groups
  and members) to an O(N²) FoF on three particle sets:
  * uniform random particles, near percolation;
  * Gaussian blobs, including blobs centred on faces, edges and corners of
    the periodic box;
  * Soneira–Peebles hierarchical sets.

  This holds for `min_npart` = 1, 5 and 20, and for cell grids
  `nc` = auto, 1, 2, 3 and 7.
* **Determinism:** byte-identical catalogs on repeat runs. With OpenMP, they
  are also byte-identical for 1, 2, 3, 4 and 8 threads, for the FoF, SO and
  potential-minimum variants. `find_halos` on 128³ particles writes
  byte-identical catalogs serially and with 4 threads.
* **Properties:** for a halo across the box corner, the centre of mass and
  the mean velocity are recovered. The virial density and `T` match their
  definitions. SO mean density is at least Δρcrit. The potential-minimum
  centre is correct.
* **Mass function:** a synthetic halo population drawn from Sheth–Tormen at
  z = 3 is placed on 64³ particles as truncated isothermal spheres on a
  Poisson background. The FoF counts per mass bin agree with ST within
  Poisson noise plus 30%. Typical output: 160/171, 55/45, 7/11.
* **Source models:** the SFR sums to `ρ̇*(z) V`, `L/sfr = 8e27`, the
  threshold and ordering hold, and the dynamical formula is correct.
* **Legacy against gfortran** (`tests/halo_catalog_ref.f90`, built when
  gfortran is found): the C++ `read_halo` is bit-identical to the verbatim
  Fortran `read_halo` (stream access) at 11 redshifts. The tag, `Nhalo`,
  `Nsource`, `sfr_scale` bits, and the rescaled `L` and `sfr` bits all
  match. The redshifts cover:
  * above the table;
  * an exact row match;
  * an exact-tie tag;
  * extrapolation past the end of the file.

  The `sfr.txt` used includes continuation lines, commas, `d` exponents and
  extra values. Further checks:
  * gfortran's sequential and stream files read identically;
  * catalogs written by C++ are read back by gfortran bit for bit;
  * 18 redshift tags are identical to gfortran F editing.
* **`link_halo` / `calc_halo` / metals:** cell assignment, per-cell sums and
  `w_gas(4)` increments.

### Follow-up, needs the PM port

* **Mass function from a real snapshot:** run on an evolved N-body snapshot
  (PM, or 2LPT/Zel'dovich ICs evolved with PM). Compare with ST/PS: for
  example, 128³ particles in 5 Mpc/h at z = 6–10, expecting FoF b=0.2
  within ~20% of ST above about 100 particles. `apps/find_halos` already
  reads `DM/xv_dm.<ztag>` files.
* **Full-run comparison:** if the original catalogs ever become available,
  compare them with the on-the-fly catalogs for the same snapshot. That
  would settle the open questions below.

## 7. Open questions

1. Can the original `Halo/` catalogs and `sfr.txt`, or the TC07 halo finder
   and source code, be obtained? They would settle:
   * the mass definition;
   * the SFR(M) relation;
   * the L/sfr conversion;
   * what `mstar`, `fgas`, `R` and `T` meant.
2. Should the default normalisation be `sfr.txt`, if the user has one, or
   the HS03 analytic fit?
3. Should on-the-fly catalogs be held fixed between "catalog redshifts" and
   rescaled like the legacy ones, or recomputed every hydro step? The
   current code recomputes every step. That is more expensive but smoother.
4. Should the original's `Nhalo_max = Npart/1024` bound be kept?
5. Halo mass: total matter per particle (as in the original PM), or
   `(1 − fb)` × that, if a DM-only mass is wanted?
