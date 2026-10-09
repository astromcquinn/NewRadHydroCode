# The unsplit hydro scheme with dual energy (`hydro_scheme = unsplit`)

*Hydro solver thread, 2026-10-09. Fixes the spurious heating of cold
supersonic gas found by the hydro test suite (`docs/hydro_tests.md` §6, §7,
§9) by bringing in the unsplit solver the AMR plan needs anyway, behind a
flag, validated with the existing suite. The legacy scheme is untouched and
stays bit-identical to the Fortran harnesses.*

Contents: 1. the scheme; 2. the dual-energy formulation (every constant and
its provenance); 3. how the step sits in the loop; 4. the RT coupling and
the state; 5. parameters; 6. results of the hydro suite, both schemes side
by side; 7. quirks of the new scheme (U1–U5); 8. the coupled run; 9. cost;
10. what the new scheme does worse, and when to use which.

## 1. The scheme

`include/radray/kernels/hydro_mh.hpp` is the unsplit MUSCL-Hancock + HLLC
kernel of the framework spike (`spike/common/hydro_mh.hpp` on branch
`framework-spike`), moved into the library, extended, and maintained here
from now on. It is header-only and framework-independent (host/device
macros from `kernels/rr_hd.hpp`, fixed stencils, no library state) so that
the AMReX thread compiles the same header for its patches. The host driver
on the uniform periodic grid is `src/hydro_unsplit.cpp`.

One substep of length dt (dx = 1 cell, double precision throughout):

1. primitives W = (ρ, v, p) and monotonised-central limited slopes in each
   direction (cells i−1, i, i+1). The pressure comes from the thermal energy
   selected by the dual-energy rule (§2);
2. Hancock predictor W_half = W − dt/(2dx) Σ_d A_d(W) ∂_d W (the full 3-D
   sum, no transverse flux corrections: the RAMSES-type unsplit MUSCL);
3. face states W_half ± ∂_d W/2 and an HLLC flux (Toro 2009 §10.4, Davis
   wave-speed estimates) at every face. The face also returns the
   interface velocity v* (the contact speed S_M of the HLLC fan, or the
   velocity of the supersonic side) and the upwind side;
4. U_new = U − dt/dx Σ_d (F_{d,+} − F_{d,−}). The internal energy density
   and the four species densities (`w`) are advected with the HLLC mass
   flux times the upwind face value (the specific internal energy of the
   upwind face state; the species reconstructed along the face direction
   with a 1-D Hancock step), and the internal energy gets the −p∇·v work
   (§2);
5. per cell: the density floor `Dgas_min`, the 1 K temperature floor of
   `check_temp` applied to both energies (per cell, no redistribution to
   the neighbours), the η₂ synchronisation and the total-energy reset (§2).

**Passes.** Predictor on the interior plus one ghost layer, fluxes and
update slab by slab (each k slab computes its own x, y and both z faces,
so no data is shared between slabs), then the per-cell finalisation. Every
value is computed from a fixed stencil of the previous pass by the same
function, the only reductions are a max (the CFL speed) and integer
counters, so the result is **the same bits for any thread count**
(`thread_full_unsplit_hydro_pm` checks it on the full loop, see §8).

**Time step.** The outer step is the original `calc_timestep` (unchanged,
§3). Inside, each call advances by its interval in equal substeps of at
most dt_cfl = `unsplit_cfl` / max_cells Σ_d(|v_d| + c) (RAMSES's
definition of the signal speed; the unsplit MUSCL-Hancock scheme without
transverse corrections is stable for dt Σ_d(|v_d| + c)/dx < 1), with the
default `unsplit_cfl = 0.8`. With the Courant step of `calc_timestep`
(dt_sim = 1/max(|v| + c) in cell units) a half step dt_sim/2 needs 1–2
substeps for flows along one axis and up to 2 for isotropic ones, i.e.
2–4 substeps per hydro step (the legacy scheme does 6 sweeps); measured
counts are in §6 and §9.

**Precision.** The state stays in `GasGrid::u` (float: ρ, ρv, E in the
lab frame), `GasGrid::w` and the new `GasGrid::e` (float). A call converts
to double, sub-cycles, and rounds to float once at the end, so there are
two float round trips per hydro step (one per half step), as for the
legacy scheme. What this costs is measured in §6 (the sound waves at
A = 10⁻³ sit on the same float floor as the legacy scheme; the cold-gas
temperature is unaffected because it is carried by `e` with its own float
mantissa, see §2).

## 2. The dual-energy formulation

Any scheme that evolves the total energy E and recovers the thermal energy
as e_T = E − ρv²/2 loses it when ρv²/2 ≫ e: the truncation error of E,
of order the local total energy, is thermalised. In cold supersonic flows
(Mach 10–1000 in the pre-reionisation IGM) that is the dominant error of
the legacy scheme and would be the same for the unsplit kernel on its own.
The dual-energy formulation of Bryan et al. (1995, Comput. Phys. Commun.
89, 149), as used by Enzo and Nyx, carries the internal energy density
ρe as a second variable:

* **Evolution.** ∂_t(ρe) + ∇·(ρe v) = −p ∇·v. Discretised as a passive
  scalar advected with the HLLC mass flux (F_e = F_mass × e_upwind, where
  e_upwind is the specific internal energy of the upwind *face state* of
  the main reconstruction, so the advection is consistent with the
  pressure the Riemann solver saw), minus dt p_half Σ_d (v*_{d,+} −
  v*_{d,−})/dx with the cell's Hancock half-step pressure p_half and the
  Riemann interface velocities v* (time-centred, so the work term is
  second order). The species `w` are advected by the same mass flux.
* **Selection (η₁).** The pressure of a cell is (γ−1)e_T if
  e_T / E > η₁, else (γ−1)ρe. This is what every reconstruction uses
  (`thermal_energy` in the kernel), and what `pm::particle_mesh_unsplit`
  keeps fixed during the gravity kick.
* **Synchronisation (η₂).** After the update, ρe := e_T when
  e_T / max_nb(E) > η₂, the maximum being over the cell and its six face
  neighbours. Shocks, where the internal-energy equation misses the
  entropy production, thereby re-synchronise the field from the
  conservative total energy; cold cells next to energetic ones keep
  their own value.
* **Total-energy reset.** In cells whose pressure comes from ρe
  (e_T / E ≤ η₁) the stored total energy is reset to E := ρv²/2 + ρe
  (as Nyx does). So `u(5)` always carries the dual-energy temperature and
  every diagnostic that reads `u` (`summarize_gas`, the RT grid
  temperature, the step fields, the stats) is consistent with the hydro,
  at the cost of exact energy conservation in those cells only: the
  change is |e_T − ρe| ≤ η₁E per cell and step, i.e. the thermalised
  truncation error that a total-energy scheme would have put into heat.
  Cells with e_T / E > η₁ are never touched, so hot or subsonic flows
  (the Sedov blast, the tubes, the sound waves) conserve energy to
  round-off (§6).
* **Floors.** ρ ≥ `Dgas_min`; e/ρ ≥ 1 K (`check_temp`'s floor,
  `1/temp1_unit`) applied to both e_T and ρe. Floor hits are counted
  (`UnsplitCounters::floor_T`, `floor_rho`) and printed at the end of a
  run.
* **Initialisation.** `init_unsplit_state`: ρe = max(E − ρv²/2, ρ·1 K)
  from the initial conditions or a legacy checkpoint; a non-zero grid
  velocity (moving-frame initial state, `v_gas.bin`) is folded into `u`
  and zeroed. With `dual_energy = 0` the field is evolved but never used
  (e_T always), for comparison runs.

**Constants and their provenance.**

| Constant | Value | Where | Provenance |
|---|---|---|---|
| η₁ | 10⁻³ | `dual_eta1` | **Unconfirmed.** The value commonly quoted for Bryan et al. 1995 (and Enzo's default), written from memory: the cloud session cannot reach arXiv/ADS and the user has not yet supplied the paper. To be checked against the paper. |
| η₂ | 0.1 | `dual_eta2` | **Unconfirmed**, same source. The neighbour set (the six face neighbours) is this implementation's choice; Bryan et al. may use the sweep neighbours only. |
| CFL | 0.8 against Σ_d(\|v_d\| + c) | `unsplit_cfl` | RAMSES's `courant_factor` default for its unsplit MUSCL scheme (from memory; the stability limit itself, Σ_d(\|v_d\| + c) dt/dx < 1, is standard). |
| γ | 5/3 (double) | kernel | The supercomoving equations have no expansion terms for γ = 5/3 (`docs/hydro_tests.md` §1); the legacy scheme uses the real(4) 5./3 (H10), the summaries and the RT still do. |
| ρ_min, p_min | 10⁻¹⁰, 10⁻¹² (code units) | kernel | Positivity guards of the spike kernel, never reached in the tests (the `Dgas_min` floor is far above). |

The entropy-conserving variant (Ryu et al. 1993) was not built: the
dual-energy version met the acceptance criteria (§6).

## 3. How the step sits in the loop

The original (`radhydro.f90`) does three forward sweeps (each dt_sim/2 in
one direction), the PM drift–kick–drift with the full dt_sim, three
reverse sweeps: Strang splitting of hydro and gravity with the kick in the
middle. The unsplit scheme keeps exactly that structure and the same call
sites (`hydro_half_step` in `apps/radhydro.cpp`, the suite driver and
`radhydro_hydro_step`):

```
calc_timestep (unchanged: dt_sim = min(Courant, PM), ra_max; nstep_rt)
hydro_half_step First   : unsplit step over dt_sim/2 (sub-cycled on its CFL)
clock += dt_sim/2; calc_sim_var
particle_mesh_unsplit   : DM drift, link, density (gas + DM), gravity,
                          DM kick (unchanged), gas kick on (ρv, E) with the
                          thermal energy fixed, DM drift
hydro_half_step Second  : unsplit step over dt_sim/2
clock += dt_sim/2; calc_sim_var
radiative_transfer, summaries, output  (unchanged)
```

The gravity source is therefore a kick of the conserved variables at the
half step: ρv += ρ g dt_sim, E := ρv_new²/2 + e_sel, with e_sel the thermal
energy selected by the η₁ rule, so gravity changes the kinetic energy by
exactly the work done and never the thermal energy (the legacy scheme
kicks the frame velocity, which has the same property). This is second
order (Strang) like the original, needs no change of `calc_timestep` and
hence gives the same `nstep_rt` for the same `cmax_gas`; what does change
is `cmax_gas` itself (U1 below): the legacy scheme's `max(|v − v_g| + c)`
is `max(|v| + c)` in the lab frame, so the Courant step is shorter where
the bulk flow is fast, and the hydro step count of a run differs between
the schemes (§8 reports both).

The supercomoving variables need no expansion terms for γ = 5/3, so the
kernel is the Newtonian one; the temperature floor, `Dgas_min` and the
`summarize_gas` clamps (H6) see the scale factor through the units as
before. `HydroState::v` (the grid velocity) is identically zero with the
unsplit scheme, so every consumer of `u` that adds `v` (the summaries, the
stats, the RT, the data files) reads lab-frame values without change.

## 4. The RT coupling and the state

* **State.** `GasGrid::e` (float, N³, the internal energy density in the
  units of `u(5)`) is allocated only with `hydro_scheme = unsplit`; the
  legacy code never reads it. It goes through the checkpoint (format
  version 2; version-1 files are still read, with `e` rebuilt from `u` on
  the first step), the step field output (`e_gas.bin`, written when it
  exists), `write_data_1/2` (`Gas/e_gas.<ztag>`) and the restart test
  `restart_smoke_unsplit` (bit-identical restart with RT on, §8).
* **RT.** `radiative_physics` → `photo_cell` heats and cools a cell by
  changing `u(5)`. With the field present it takes the thermal energy from
  `e` instead of `u(5) − KE` and writes the result to both. That is the one
  RT edit (`src/rt_tools.cpp`, two lines); with the legacy scheme the
  pointer is null and the arithmetic is the original's, bit for bit
  (checked by every stored-reference and live-Fortran test, §8). The RT
  thus sees the dual-energy temperature exactly, and the hydro sees the RT's
  heating through both variables.
* **Everything else** (`summarize_gas`, `calc_timestep`, the RT grid
  temperature `rt_grid%temp`, the movie frames, the P(k)) reads `u`, which
  carries the dual-energy temperature through the total-energy reset of
  §2 to float precision: E − KE in float has a relative error 6×10⁻⁸·E/e,
  i.e. 0.03 % at e/E = 2×10⁻⁴ (Mach 100) and 6 % at e/E = 10⁻⁶ (Mach
  1300, below the 1 K floor in practice). The analysis of the suite reads
  `e_gas.bin` where it exists for that reason (`hs_common.read_state`).

## 5. Parameters

| Key | Default | Meaning |
|---|---|---|
| `hydro_scheme` | `legacy` | `legacy` (the original, bit-identical) or `unsplit` |
| `unsplit_cfl` | 0.8 | Courant number of the substeps against Σ_d(\|v_d\| + c) |
| `dual_energy` | 1 | 0: total energy only (the kernel without the formulation; for comparison) |
| `dual_eta1` | 10⁻³ | pressure from E − KE when (E − KE)/E > η₁ (unconfirmed, §2) |
| `dual_eta2` | 0.1 | ρe := E − KE when (E − KE)/max_nb E > η₂ (unconfirmed, §2) |

All drivers (`apps/radhydro`, `apps/hydro_driver`, the suite driver) take
them in `params.txt` and with `--set`. The run summary prints the
`UnsplitCounters` (calls, substeps, the largest number of substeps of one
call, the smallest CFL step, cell-substeps whose pressure came from `e`,
cell-substeps synchronised, floor hits); `apps/radhydro` also writes them
to `timing.txt`.
