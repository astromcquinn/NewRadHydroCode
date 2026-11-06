# The moment RT on the refined levels (Thread C, stage 2)

*Branch `moment-amr`, 2026-11-06. The user's decision after
`docs/absorbers_closures.md`: an AMR run with M1 or the OTVET cone where
the absorbers are resolved. The moment RT on the AMReX backend
(`src/amrex/moment_level.cpp`, `docs/stage_e_moments.md`) runs on level 0
only and `level_of` refuses `amrex_max_level > 0` with it; this document
is the design of the refined levels, mirroring the rays' hierarchy
(`docs/stage_f.md` §3, `src/amrex/rt_level.cpp` `Hierarchy`).*

## 0. Summary

A `MomentHierarchy` holds one moment level per refinement level of the gas
state: level 0 the existing `MomentLevel` on z-slabs, every refined level
a `MomentFine` on the hydro's grids of that level. Per coarse RT substep
(one base cell of light, the user's decision of 2026-10-14 for the rays)
a level-l transport runs 2^l × `moment_nsub` GLF subcycles of Courant
number 1/`moment_nsub` on its own cells, with its ghost cells outside the
level's grids holding the parent's moments; the physics (`photo_cell`)
runs on every level's own cells with the level's units, the parent's
cells under a finer level have neither sources nor physics (the finer
level has them), and after the substep the fine moments, cells and gas
are averaged down onto the parents, so level 0 is always the composite.
Milestone 1 is M1 (a local closure) without refluxing at the coarse-fine
faces, its photon non-conservation measured by the budget; milestone 2
adds the OTVET tensors (prolonged from level 0, the cone's attenuation
at the base resolution) and the refluxing if the budget asks for it; the
production run is the 512³ box with two levels against the rays' refined
run.

## 1. What exists

| piece | state |
|---|---|
| `MomentLevel` (`moment_level.cpp`) | M1, OTVET thin, split, screen, cone on z-slabs, MPI + OpenMP, bitwise the serial solver for any rank count; the sums in the serial order (chains over the ranks in z order) or in 128-bit fixed point (`amrex_moment_stats = fixed`); the checkpoint in the serial format; `recomb_rad` |
| the rays' `Hierarchy` (`rt_level.cpp`) | the pattern: `sync` (the levels follow the gas state's grids, rebuilt at a regrid), `set_covered` (`makeFineMask`), `in_patch`, the parent's TF in the ghosts outside a level (`fill_ghost_tf`), the deposits across the boundary (`exchange_deposits`), `average_down` (cells and gas, 8 children in double), `place_sources` (a covered base cell's L on the finest cell holding each source), per-level units (`units_of`: len and vol ÷ 2^l, 8^l), the checkpoint's level directory (`write_levels`, `restore_levels`), the statistics summed over the levels in fixed point (`publish_stats`) |
| the gas state (`gas_state.hpp`) | `fine[l-1]` = {ba, dm, G} per refined level, `finest`; `amrex_hydro_level0_is_current()` |
| the hydro level | `FluxRegister` refluxing (`hydro_level.cpp`: `CrseInit`, `FineAdd`, `Reflux`) as the template for milestone 2 |

## 2. Design

### 2.1 Objects

`MomentFine : CellLevel`-like (its own constructor on a given level,
BoxArray and mapping): `N = Nrt << lev`; `C` (the RT cells) on the
level's grids; `M` (the 25 moment components, one ghost layer), `Wn`;
`in_patch` (1 on the level's cells and the ghosts covered by its other
boxes, 0 outside the grids), `covered` (the cells under the next level);
`ul` the level's units (`units_of`); `GA` for `gammaHI_avg`; `RB` for
`recomb_rad` (milestone 2); the fixed-point statistics only (the serial
order needs slabs, which the hydro's grids are not: a hierarchy requires
`amrex_moment_stats = fixed`, or sets it). No `H`/`Nd` of its own in
milestone 1 (M1); milestone 2 prolongs the parent's.

`MomentHierarchy : LevelBase`: `levels[0]` the `MomentLevel`, `levels[l]`
the `MomentFine`s; `finest`; the budget (injected, absorbed, lost) in base
units summed over the levels.

### 2.2 The substep (as built)

```
calc_rt_1 on every level's cells (its own gas, its units)
inject (moment_inject = 0): level 0 on its uncovered source cells, level l
        on its source cells (place_sources: a covered base cell's L on the
        finest cell holding the source)
for l = 1 .. finest:  fill_ghost_moments(l): the parent's N, F, E in the
        ghosts outside the level (per volume, as they are), kept as Mg
level 0 transport: nsub GLF subcycles of C = 1/nsub; its flux across each
        face into a covered neighbour is recorded in RF (minus sign)
for l = 1 .. finest (coarse first, Berger & Colella):
        fill_ghost_moments(l) again (the parent after its transport), kept
        as Mg1; the level's 2^l x nsub subcycles with its outside ghosts at
        Mg + (it + 1/2)/(2^l nsub) (Mg1 - Mg); the photons its valid cells
        send into the outside ghosts go, per coarse volume, into RFup
for l = finest .. 1:  reflux: RF of level l-1 += RFup of level l; the
        uncovered cells of level l-1 next to the level get RF (the finer
        level's flux instead of their own; realizability enforced);
        average_down_moments(l): the level's moments onto the covered
        parent cells
inject (moment_inject = 1)
deposit N, E, sum F into every level's cells; physics on every level's
        uncovered cells (the level's units, X.cells = 2^l so Gamma = n sigma c_hat);
        fgam scaling of each level's M (absorbed, uncovered cells only)
average_down_moments again; the statistics in fixed point over the levels
```

**Units.** A level's N, F, E are photons per ⟨n_H⟩ and the level's cell
volume (`inject_cell` with the level's `vol_rt`), so a value is a
density: a child takes the parent's value as it is, a parent takes the
mean of its 8 children. Every budget (injected, absorbed, lost, in
flight) is kept in base units: a level-l total ÷ 8^l.

### 2.3 Conservation at the interface

The reflux registers (`RF`, `RFup`, 2^32 fixed point per moment
component, so the sums are layout independent) replace the parent's GLF
flux across a coarse-fine face by the sum of the finer level's fluxes
over its subcycles (`glf_face`: the GLF face flux from the same
quantities the update uses). With the fine level run *before* the parent
on the parent's start-of-step ghosts, the fine level drained ghost cells
whose photons the parent spent elsewhere in the same substep, and the
reflux drove parent cells negative: the realizability clamp then created
photons at 0.1 % of the in-flight total per substep, 6 % of the
injection over the I-front run. The coarse-first order with the ghosts
interpolated in time between the parent's old and new values (§2.2)
brings the transport's non-conservation to 10⁻¹¹ per substep and the
clamp to 10⁻⁶ over the run. `RADRAY_MOMENT_AMR_DEBUG=1` prints, per
substep, the photon change of the transport, the injection and the
physics against the counters, and per step the budget residual
J − A − L − I.

### 2.4 Regrid, restart, outputs

- **Regrid** (`sync` at every prep): a level whose grids changed is
  rebuilt; its cells and moments where the new grids hold them are copied
  from the old level (`ParallelCopy`), the rest prolonged from the parent;
  the covered masks follow.
- **Checkpoint**: level 0 in the serial format as now (the composite,
  which a single-level restart can read); the refined levels' `M` as
  VisMF files in the checkpoint's level directory (`write_levels`,
  `RT_moments_l`), restored by `restore_levels` when present, else
  prolonged from the parent.
- **Outputs**: level 0 is the composite after `average_down`, so `fields`,
  the movie and `gammaHI_avg` read it as now; `fields_levels = 1`
  (finest-level slices) is not implemented for the moments.
- `amrex_balance = 1` (ray-weighted) stays refused for the moments.

### 2.5 The OTVET tensors on the levels (milestone 2)

The tensor is computed once per hydro step on level 0 from the full
source list (every source is on the base grid's list); a refined level's
`H` (and `Nd`) take the parent's values per child (piecewise constant).
For the thin, split and screen kernels this loses only the sub-cell
structure of a smooth field; for the cone the HI attenuation of the
near-field march is the base grid's, so an absorber resolved only on a
fine level does not shadow the tensor. That is the single-absorber
test's warning (`docs/absorbers_closures.md` §5: the cone's advantage is
its attenuation) and the reason for a later fine-grid march; the
production comparison will carry M1 and the cone with this caveat.

## 3. Tests

1. `amrex_moment_amr_ranks_smoke`: the smoke problem with
   `amrex_max_level = 2`, `amrex_refine_box = 0,0,0,7,7,7`, M1, 1/2/3
   ranks byte-identical (`rank_check.py`), and the restart
   (`restart_check.py`) byte-identical.
2. `amrex_moment_amr_whole`: the I-front problem (`make_ifront.py`) with
   one level over the whole box against the uniform run at 2N: the front
   radius within the tolerance of the rays' `amrex_rt_amr_ifront`.
3. `amrex_moment_amr_ifront`: the refined region over the centre, the
   front crossing the boundary: no jump, the radius against the 2N run.
4. The photon budget residual printed per step (`RADRAY_MOMENT_AMR_DEBUG`).
5. `amrex_max_level = 0`: every existing moment test bitwise (the
   hierarchy is not built).

## 4. The production run

The rays' refined run exists: `fesc512_b2.5_amr` (`docs/stage_f.md` §3:
40 Mpc/h, 512³ + two levels above 8 ρ̄ and 64 ρ̄, 1024³ particles on the
2048³ mesh, 0.1c, case B, f6 = 0.135, the exposure deposit, the sink at
3000 ρ̄, cooling suppression 0, FoF every 30 Myr, before P11; 37.5 h on 5
nodes, 36 000 core-hours). The matched cone (and M1) runs take its
`params.txt` with `rt_method` and the closure keys and its ICs, so the
comparison differs in the RT only, and `photon_consume_dense = 0` is
set explicitly to match it (the comparison is of methods, not of P11).
Cost: the moment RT was 2.4× cheaper than the rays on the uniform grid;
on the hierarchy the fine levels' transport is 2^l subcycles on 6 % of
the base cell count, so ≈ 8–12 node-days per run; disk ≈ 0.3 TB per run
with the particles. Decide the pair (cone, M1) when milestone 2 passes.

## 5. Status

- 2026-11-06: design. Milestone 1 built in `src/amrex/moment_level.cpp`
  (`MomentHierarchy`; `CellLevel` on a given level and grids;
  `MomentLevel` on a refined level with fixed-point sums, the covered
  mask, the budgets in base units, the substep in phases): the smoke
  problem with two levels is byte-identical over 1, 2 and 3 ranks and
  through a restart (`amrex_moment_amr_ranks_smoke`,
  `amrex_moment_amr_restart_smoke`); the budget residual J − A − L − I is
  10⁻¹⁰ per step (`RADRAY_MOMENT_AMR_DEBUG=1`). Found on the way: the
  hydro's refined levels refuse z-slabs, so a hierarchy's base level runs
  on boxes with the fixed-point sums (`exact` only on slabs); the parent's
  covered cells must be left out of the absorbed, lost and in-flight sums
  (they counted the children's absorption again, 3× with three levels);
  a refined level injects with its own cell volume (the base units gave
  8× too few photons: the I-front test caught it).
- Milestone 2 built: the OTVET tensor on a slab shadow of level 0
  (`tensor_on_shadow`), H and N_dir prolonged onto the levels; the cone
  smoke tests pass (1/2/3 ranks, restart). The refluxing and the
  coarse-first subcycling (§2.3) built after the M1 I-front test
  (`amrex_moment_amr_ifront`: the refined front against the uniform M1
  runs at N and 2N) showed the front running ahead once it left the
  refined region; it passes with them. `tools/hyak/amr_moment_runs.sh`
  makes and submits the siblings of `fesc512_b2.5_amr`.
