# Refined PM mesh: `pm_refine`

*Branch `pm-refine`, 2026-10-11. Goal: resolve halos of 1.3×10⁸ M☉
(20 particles of 6.4×10⁶ M☉) in the 40 Mpc/h, 512³ runs. That needs
Ndm = 1024 and a PM mesh fine enough for 1024³ particles.*

## 1. What it does

`pm_refine = r` (params.txt or `--set`, in `radhydro`, `pm_driver` and the
hydro-suite driver; default 1) solves gravity on a mesh of **Nm = r·Ngas**
cells per dimension instead of the gas grid. The DM uses the mesh exactly as
the original uses the gas grid. The only changes are the mesh size and the
particle-to-mesh ratio:

| | original (r = 1) | r > 1 |
|---|---|---|
| mesh | Ngas | Nm = r·Ngas |
| particles per mesh spacing | `N_dm_gas = Ngas/Ndm` (integer) | `N_dm_m = Nm/Ndm` (integer ≥ 1) |
| CIC, `link_dm`, `accelerate_dm` | on the gas grid | on the mesh (`x_dm_m`, `nmr`) |
| particle mass (mesh units) | (1−f_b)·N_dm_gas³ | (1−f_b)·N_dm_m³ |
| gas density on the mesh | the gas cell | every mesh cell gets the density of the gas cell it lies in (the cell's mass spread evenly over its r³ mesh cells) |
| force on the gas | `force(i,j,k)` | mean of `force` over the r³ mesh cells of the gas cell, divided by r (mesh cells → gas cells) |
| FFTs | serial FFTW plans | threaded FFTW (`libfftw3f_omp`), all threads |

* **Adjoint spread and average.** The gas spread and the force average are
  transposes of each other, and the DM CIC deposit and interpolation are as
  well. With the symmetric Green's function and the antisymmetric 4-point
  gradient, every mesh element's force on every other is equal and
  opposite. Total DM + gas momentum is therefore conserved to round-off, as
  in the original (§3).
* **Units.** Poisson's equation has the same constant in mesh units as in
  gas units (lap φ = (2a/3)δ). The mesh force is then r times the gas-unit
  force, hence the 1/r on the gas. DM velocities stay in DM-grid units
  (`/x_dm_m` in the kick). `amax_dm` is in mesh units, so the PM time step
  `dt_ap = 0.5/sqrt(amax_dm)` keeps particles to the same fraction of a
  mesh cell per step. That is √r shorter in time for the same physical
  acceleration, which is the right criterion for the finer mesh.
* **Ndm > Ngas.** `x_dm_gas = Ngas/Ndm` is now a real number (0.5 for
  Ndm = 2·Ngas), and is used only where DM velocities are converted to
  gas-grid or physical units: `summarize_dm`, and the halo velocities in
  `make_halo_units` (`src/halos.cpp`). It was integer division there, which
  gives 0 for Ndm > Ngas. `PMConfig::make` requires Nm to be a multiple of
  Ndm. With r = 1 this is the original check, Ngas a multiple of Ndm (Q10).
* **Domains.** `init_pm` scales the gas domains of `init_domain_gas` by r
  and keeps their order, so the deposit `cic_mass_assignment_all` stays
  deterministic, in the same domain order.
* **The power spectrum** (`Pk/pk.*`): DM on the mesh, gas on the gas grid
  (separate buffer); rows k = 2..Hgas+1 as before.
* **`pm_driver`** writes `phi.bin` on the mesh (Nm³).
* **Checkpoint.** The PM mesh holds no state between steps (only the
  particles do, and they are checkpointed with their count checked), so the
  checkpoint format is unchanged. `pm_refine` is not stored: like
  `c_reduced`, it can be changed at a restart. The checkpoint for Ndm > N is
  tested (`restart_smoke_pm_refine`).

**pm_refine = 1 is byte-identical to before.** Every r > 1 branch is
separate (`gas_force`, `calc_density_field`), r = 1 uses the original
serial FFT plans and never initialises FFTW's threads, and all other
expressions reduce to the original ones exactly (`nmr = ngasr`,
`x_dm_m = N_dm_gas`, `Nmcell = Ncell`). Checked: the PM, hydro and
full-loop stored references, the live Fortran comparisons, and the rest
of the fast and nightly suites (§3).

## 2. Initial conditions

`tools/ics/make_ics.py` accepts Ndm > Ngas (a multiple or a divisor of
Ngas). The gas fields keep only the modes |m| < Ngas/2 of the particles'
field, since the gas grid cannot hold the rest.

`--phases-from-ndm N` (zeldovich) takes the white-noise modes |m| < N/2 from
an Ndm = N run with the same seed (the same draws, rescaled by
(Ndm/N)^{3/2}). The rest comes from an independent stream. The large scales
then match the N³ run: the 1024³-particle ICs of the 40 Mpc/h box
(`--seed 7 --phases-from-ndm 512`) have the gas field of the 512³ production
ICs to 4.5×10⁻⁸ of its rms (float round-off of the larger FFT).

`tools/production/make_production.py --pm-refine 2` writes `pm_refine = 2`,
`Ndm = 2N` and those ICs. `--data-files dm` writes only the particles at
the zdata1 outputs (26 GB per output at 1024³).

1024³ ICs: 9 min 52 s and 293 GB on one node (numpy, single-threaded FFTs).

## 3. Tests

RESULTS

## 4. Cost

COST

## 5. Notes

* **Q14 (preserved).** The DM speed and acceleration diagnostics of the
  `calc_timestep` log (`vmax_dm`, `vrms_dm`, `amax_dm`, `arms_dm`) are
  converted to km/s with the gas-grid units, but they are in DM-grid and
  mesh units. They are off by x_dm_gas (2 in the original's default
  Ndm = Ngas/2) and by r for the accelerations. They are diagnostics only;
  the time step uses the grid values.
* `env.sh` now selects `fftw-3.3.10-float-omp` (the same SSE2 build with
  `--enable-openmp`) unless `RADRAY_FFTW` is already set. A shell that
  sourced an older `env.sh` still has the serial prefix in `RADRAY_FFTW`,
  and CMake caches the library path. Build in a fresh directory, or pass
  `RADRAY_FFTW` explicitly. CMake reports "Threaded FFTW: …", and
  `radhydro` prints "(PM mesh …, threaded FFTW)" or "serial FFTW". Without
  the threaded library, pm_refine > 1 still works, with serial FFTs.
