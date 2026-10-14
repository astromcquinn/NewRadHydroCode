# Stage E3: the state on the levels (`amrex_state = 1`)

*Started 2026-10-13 on branch `stage-e-state`, after Stage D (the RT on
AMReX, `docs/stage_d.md`), E1 (the hydro, `docs/stage_e_hydro.md`) and E2
(the PM, `docs/stage_e_pm.md`). Those three keep the full serial state on
every MPI rank and copy it to and from their level every step. This stage
makes the levels the owners of the state, so that a rank holds only its
boxes and its particles, and the host arrays exist only transiently on
rank 0 for the file formats. It is what lets 512³ run on one node (today
every rank holds about 100 GB of redundant state at 512³,
`docs/stage_d.md` §4.4) and removes the per-step gathers (2–3 s of the
5 s hydro step at 256³ on 2 nodes, `docs/stage_e_hydro.md` §4.1).*

## 0. Summary

| Item | State |
|---|---|
| The map of every producer and consumer of the host state (the survey of 2026-10-13, §1) | done |
| `amrex_state = 1`: one gas BoxArray and DistributionMapping shared by the RT, hydro and PM levels; the gas state as a persistent MultiFab; the hydro and RT levels read and write it locally; the per-rank host gas arrays dropped | **done (phase 1, 2026-10-13)**, §4 |
| The DM particles persistent in the PM container (drift and `Redistribute` on the level, no rebuild from the host); the halo finder on a gather of the particles to rank 0 with the catalogue broadcast; the host PM arrays only on rank 0 | **done (phase 2)**, §4 |
| Rank-0 I/O through slab gathers: fields, movie frames and slabs, data files, checkpoints in the existing formats; chunked ray gather (the 2 GB `MPI_Gatherv` limit); restart by scatter | **done** (the mirror on rank 0 is refreshed before every output and checkpoint; a restart scatters rank 0's grids and particles onto the levels) |
| Statistics: maxima by exact reduction (the timestep uses only `cmax_gas` and `amax_dm`), the order-dependent serial sums of `summarize_gas`, `write_stats` and `write_frame` computed on rank 0 from the gathered fields (bitwise the serial order) | **done**: every sum stays the serial one on rank 0's mirror, so every file is bitwise that of the host-state run; the per-step cost is the gathers of §4 |
| Validation: the full loop with `amrex_state = 1` byte-identical to the same run with `amrex_state = 0` (the three backends with the host state), for 1/2/3 ranks and 1/2 threads; bit-identical restarts | **done**: `ctest -L amrex_state` (5 tests), §3; the whole `amrex` label 159/159 on klone |
| klone: 256³ on one node bitwise identical to the host-state run, the other ranks at 4.5 GB; 512³ on one node, 1–4 nodes scaling | 256³ done (§4), 512³ running |

Ground rules kept: nothing changes without `amrex_state = 1` (default 0);
the serial code and the Stage D/E1/E2 paths are untouched in behaviour.

## 1. What holds the state today (survey of 2026-10-13)

Per rank, with all three backends on, steady state:

| Host array | Bytes per gas cell (or per particle) | Written by | Read by |
|---|---|---|---|
| `gas.u`, `gas.w`, `gas.e` | 20 + 16 + 4 | hydro `copy_out` (every rank broadcasts its block), PM `gather` (BOR all-reduce), RT `sync_fields` (slab broadcast), `add_halo_metals`, `summarize_gas` clamps, `fold_frame`, `init_unsplit_state` | the three levels' copy-in, `prep_rt`'s host loop, summaries, I/O, checkpoint |
| `hs.v` | 12 | allocated by `init_hydro`; identically zero with `unsplit` | `fold_frame` scans it every half step |
| `RTGrid s.rt` (36 floats + `htoc`) | 144 + 32 | `sync_fields` (all but `Luv`, `sfr`, `htoc`), `calc_halo` (`Luv`, `sfr`) | the RT level's `prep` (diagnostic fields only), `write_stats`, `write_frame`, `write_fields`, `write_data_files`, checkpoint |
| `pm.xv`, `pm.ll`, `pm.hoc`, `pm.density`, `pm.fmesh` | 32 per particle; 8 + 4 (+ 12) per mesh cell | PM `gather`, `calc_drift` (host, twice per step) | PM `prep` (every rank walks all particles), `update_sources` (FoF), `summarize_dm`, P(k) (host `link_dm` + CIC + FFTW on every rank), `write_fields`, `write_data_files`, checkpoint |
| serial ray pool | virtual only under `rt_backend = amrex` | | |

About 244 B per cell per rank plus 32 B per particle (33 GB at 512³,
before the transients), against 40 B per local cell for the gas on the
level. The hard consumers, i.e. the ones that need all the data in one
place or in a layout of their own: the FoF and SO halo finder
(`src/halos.cpp`), the P(k) and Poisson FFTs, the particle order by
index (`xv_dm` files, the checkpoint, the BOR gather by id), the
monolithic rank-0 checkpoint with its 2 GB ray gather, the serial-order
double sums of the statistics, and the single-file outputs in Fortran
order.

## 2. Design

### 2.1 One decomposition

`amrex_state = 1` requires `rt_backend = amrex`, `hydro_backend = amrex`,
`pm_backend = amrex`, `hydro_scheme = unsplit`, `rt_method = rays`,
`recomb_rad = 0`, `N_rt_gas = 1`. The three levels then build their
BoxArray from the gas domain with the one `amrex_max_grid_size`
(`amrex_hydro_max_grid_size` and `amrex_pm_max_grid_size` must equal it
or be 0; with `pm_refine = r` the PM mesh boxes are `r` times the gas
boxes, as E2 already does), and share one `DistributionMapping` (built
once by the state and passed to the levels), so the same gas cell is on
the same rank in every level and the copies between them are local. The
hydro's `RealBox` is irrelevant to the copies (index space only).

### 2.2 The gas state

A persistent `GasState` (`src/amrex/gas_state.cpp`): one `FabArray` of
floats with 10 components (u 0–4, w 0–3, e), no ghosts, on the shared
decomposition. It is filled once from the initial conditions (rank 0
reads `u_gas.init` in slabs and `ParallelCopy`s them) or from the
checkpoint, and then:

- the hydro level's `copy_in` reads its own boxes from it (float →
  double, local) and `copy_out` writes them back (local): the per-rank
  broadcasts go away; `fold_frame` is dropped (v is zero with `unsplit`;
  the legacy checkpoint conversion keeps it on rank 0 at restart) and
  `init_unsplit_state`'s e initialisation becomes a per-cell map on the
  level;
- the RT level's `prep` reads u, w, e from it (local) and `post` writes
  w, u(5), e back (local, the same cells `sync_fields` wrote), so
  `sync_fields` is no longer called per step;
- the PM level reads the gas density and energy from it and `kick_gas`
  writes u (and e) back locally: the BOR all-reduce of the gas goes away;
- `summarize_gas`'s clamps are a per-cell map on it; its maxima are
  `ReduceRealMax` (exact; `cmax_gas` and `vmax_gas` feed the timestep
  and `c_limit`), its sums 128-bit fixed point (`vrms_gas`, the masses,
  the clumping: diagnostic, in `stats.txt` and `hydrolog`);
- `add_halo_metals` is a sparse scatter: every rank applies the sources
  whose gas cell it owns, in the same stable order.

The host `gas.u/w/e` and `s.rt` arrays are not allocated with
`amrex_state = 1`; `prep_rt`'s host loop, `hs.v` and the serial ray pool
are skipped. The RT level's diagnostic fields (dx, gamma, lambda, nu,
ngam, fgam, tau) stay in its `CellData` and are gathered only for output.

### 2.3 The particles

The PM container keeps its particles between steps (the particle id is
the host index `ip`, so the files and the checkpoint keep their order):
the two drifts become per-particle maps on the tiles followed by
`Redistribute`; the kicks write velocities in place; the BOR gather of
`xv` goes away; `summarize_dm` is a reduction (maxima exact, sums 128-bit).
P(k): the DM density is the level's CIC deposit (already a MultiFab), the
gas density its own component; the FFT is the level's pencil FFT (or
`amrex::FFT`) and `calc_pk_3d`'s shell histogram a reduction (double
sums per bin, order-independent by 128-bit accumulation; P(k) is a
diagnostic output and is validated to round-off against the host
version).

The halo finder stays the serial FoF, on rank 0: the positions (3 floats
per particle, 1.6 GB at Ndm = 512³) are gathered by id in chunks, rank 0
runs `find_halos`, the halo catalogue (small) is broadcast, and
`assign_sources`, `link_halo` and `calc_halo` run redundantly on the
catalogue on every rank; `Luv` and `sfr` are applied to the RT level's
cells as a sparse scatter in `prep` (each rank its own source cells) and
`rt_stat.Luv/sfr` summed on rank 0 in the catalogue order (bitwise the
serial sum). A distributed FoF is Stage E's later item; this gather costs
rank 0 a transient of 12 B per particle and the FoF time (10–17 s at
512³, as today). The FoF cache (`cached_catalog`) keeps working on rank
0's catalogue.

### 2.4 I/O and checkpoints

All files keep their formats. Rank 0 assembles each field from the level
with the slab `ParallelCopy` of `sync_fields` (8 k-planes at a time, one
field at a time, so rank 0's transient is one plane set, not the grid)
and writes it; the particle files gather `xv` by id in chunks. The movie
frames gather only the movie fields (or the `movie_slab` planes: three
`ParallelCopy`s of one plane each). `write_stats` and `write_frame`
compute their sums on rank 0 from the gathered fields in the serial k,j,i
order, so `stats.txt` and `frames.txt` are bitwise those of today; the
Gamma median is taken on rank 0 from the gathered field as today.

The checkpoint (version 7, unchanged on disk) is written by rank 0 field
by field through the same gathers: the RT grid block from the RT level's
cells, the gas block from the gas state, `pm.xv` from the container by
id, the backend rays by a chunked gather (replacing the single
`MPI_Gatherv` with its 2 GB limit, 25 M rays: 512³ holds 1.1 G). The
reader is the mirror: rank 0 reads a field in slabs and `ParallelCopy`s
it onto the level, reads the particles in chunks and broadcasts them so
that every rank creates its own (as E2's `prep` does today with the
locator), scatters the rays in chunks (as today). No rank holds the whole
grid at any time; rank 0's transient is one slab or chunk. A legacy
(serial) checkpoint is converted on rank 0 as today before the scatter.

### 2.5 What stays redundant on every rank

The small things: the clock, the units, the halo catalogue and the
source list, the FoF cache, the rate tables, the per-step scalars.

## 2.6 What was built (phases 1 and 2, 2026-10-13)

`include/radray/amrex_state.hpp` (AMReX-free: the parameter, the checks
and the collective calls), `src/amrex/gas_state.{hpp,cpp}` (the gas
MultiFab, the slab transfers `gather_slabs`/`scatter_slabs`, the clamps
and the sparse metals), and the three levels:

- **Decomposition.** `amrex_state_init` builds the gas `BoxArray`
  (`amrex_max_grid_size`) and its `DistributionMapping` once, before any
  level. The RT level takes them as they are, the hydro level too
  (`amrex_hydro_max_grid_size = 0` means that box size; the check refuses
  another), the PM level refines the boxes by `pm_refine` for the mesh
  and keeps the gas boxes for its gas FabArray. All three index the gas
  MultiFab through their own `MFIter`, which is valid because the
  FabArrays share the BoxArray and the mapping.
- **Gas.** 10 floats per cell (u, e, w: the hydro level's component
  order). The hydro level copies its boxes in (float → double) and out
  (rounded once, as the serial step) locally. The RT level's `prep` reads
  u, w, e from it and the sources from the host `Luv`/`sfr` arrays (every
  rank computes them from the catalogue), and keeps the radiation fields
  of the previous step in its own cells, which now persist (zeroed at
  construction, like the host arrays); its `post` writes u, w, e back
  locally instead of the all-rank `sync_fields`. The PM level fills its
  gas boxes from it and writes the kicked momenta and energy back
  locally; the velocity all-reduce is gone.
- **Particles.** `amrex_pm_load_particles` broadcasts rank 0's `xv` in
  chunks of 4 M particles and every rank keeps those whose home cell lies
  in its boxes (id = the host index `ip`, as before). `amrex_pm_drift`
  is `calc_drift` on the tiles (the host arithmetic), then the home cells
  and a `Redistribute`; `particle_mesh_unsplit` calls it around the level
  step. `amrex_pm_gather_particles` sends (id, x, v) records to rank 0 in
  chunks, which writes them into its `xv` by id: the mirror, refreshed
  when stale (after a PM step) and needed (the halo finder, the
  summaries, an output, a checkpoint).
- **Halos and sources.** `update_sources` is split: `source_catalog`
  (the FoF, the SO spheres, the source model and the FoF cache) runs on
  rank 0 on its particle mirror, `bcast_catalog` sends the records, and
  `apply_sources` (`Luv`/`sfr` per cell, the metals) runs on every rank;
  the metals go onto the gas MultiFab as a sparse scatter in the serial
  order (`amrex_state_add_metals`).
- **Mirror and outputs.** Rank 0 keeps the host `GasGrid`, `RTGrid`,
  `hs.v` and the PM arrays; the other ranks free them at start-up
  (`GasGrid::has_e()` reports the dual-energy field present without an
  array; `RTGrid::allocate_sources_only`). The driver refreshes the
  mirror collectively: `summaries` clamps the gas on the level, gathers
  it and runs `summarize_gas` on rank 0 (then broadcasts `cmax_gas`,
  `vmax_gas`, `vrms_gas`, which every rank's timestep and `c_limit`
  need); `write_stats` gathers the five statistics fields of the RT
  cells; `write_frame`, `write_fields`, the data files, P(k) and the
  checkpoint gather everything (`sync_mirror`). The checkpoint reader
  skips the grid and particle blocks on the other ranks; the driver then
  scatters rank 0's grids and particles onto the levels. The ray gather
  of the checkpoint is now chunked (8 M rays per message).
- **What every rank still holds:** `Luv` and `sfr` (8 B per cell) and the
  transient `per_cell` lists of `calc_halo` (24 B per cell) during the
  source assignment; the levels' own arrays (RT cells 184 B per local
  cell plus the deposits and the rays, hydro 330 B, PM about 30 B plus
  the particles); the catalogue and the small state. Rank 0 adds the
  mirror (about 230 B per cell plus 24 B per particle) and, at a restart,
  the serial ray pool and ray list of the checkpoint.

## 4. klone (2026-10-13)

`tools/hyak/amrex_scaling.slurm` with `EXTRA="--set hydro_scheme=unsplit
--set hydro_backend=amrex --set pm_backend=amrex --set
amrex_hydro_max_grid_size=0 --set amrex_state=1"` (the build of commit
`8ef34cc`'s phase 1, without the threaded slab loops), one hydro step
restarted from the production checkpoint, 24 ranks × 8 threads on one
`cpu-g2` node, the UCX setting of `docs/stage_e_hydro.md` §5.

**256³** (`reion256_f0.2`, step 49, jobs 41822519 and 41823171): the
fields and `stats.txt` of `amrex_state = 1` are **bitwise identical** to
the run with the same three backends and the host state on every rank
(`amrex_state = 0`); against Stage D's run with the serial hydro and PM
the differences are the PM backend's round-off (its fixed-point CIC and
pencil FFT, `docs/stage_e_pm.md`), 13 rays of 144 M.

| 256³, step 49, 24 × 8 on one node | RT [s] | hydro [s] | PM [s] | io [s] | step [s] | peak RSS rank 0 / min / mean [GB] |
|---|---|---|---|---|---|---|
| Stage D (serial hydro and PM on every rank) | 156.8 | 9.9 | 1.4 | 2.1 | 172.8 | 37.3 / 12.8 / 14.7 |
| the three backends, host state (`amrex_state = 0`) | 156.4 | 1.1 | 0.7 | 2.0 | 163.1 | 37.3 / 8.0 / |
| the state on the levels, phase 1 (`amrex_state = 1`, job 41822519) | 151.3 | 1.3 | 0.5 | 21.4 | 177.7 | 37.4 / 4.5 / |
| **the state on the levels, phase 2 + threaded transfers** (job 41823773) | 150.3 | 0.4 | 0.1 | 5.3 | **159.0** | 37.4 / **3.1** / |

**512³ on one node** (`reion512_40mpc_f0.23`, step 53, 1.12 G rays, job
41822520, the same phase-1 build, 24 ranks × 8 threads, `--mem=1400G`):
**it runs**, where Stage D ran out of memory at 24 × 8 and at 12 × 16 and
needed 8 × 24 to fit (3258 s). RT 1756.1 s, step 1848.9 s, against the
serial code's 2593.4 and 2878.7 s on the same node (1.48× faster) and
Stage D's 1460.7 s on 2 nodes and 782.1 s on 4. Peak RSS 293.0 GB on
rank 0 (the restart's serial ray pool and ray list, as before) and
35.3 GB at least on the others (the level's cells, deposits and rays,
plus the phase-1 host particles and PM arrays that phase 2 removes),
about 1.1 TB in all. Phases: `step1` 247.2, `cast1` 534.5 (`SumBoundary`
105.2), `phys` 279.4, `cast3` 386.4, `sort` 163.8, `migrate` 71.0,
`sync` 43.1 s; hydro 2.7 s, PM 3.6 s, halos 14.6 s (the FoF on rank 0),
P(k) 30.2 s (the host P(k) on rank 0's mirror), io 64.0 s (the mirror
gathers of the step and the fields of the last step, single-threaded in
this build). With phase 2 and the threaded transfers (job 41824737): RT 1783.2 s,
step 1837.5 s, io 29.4 s, P(k) 14.0 s, PM 0.7 s, the ranks other than 0
at 23.8 GB (rank 0 293 GB, the restart), the fields and statistics
bitwise identical to the phase-1 run, and **bitwise identical to the
2-node run with the same three backends and the host state on every
rank** (`amrex_state = 0`, 12 ranks × 16 threads, job 41823172: RT
1123.1 s, step 1164.6 s, the ranks other than 0 at 64.7 GB): every
field and `stats.txt`.

| 512³, step 53, the state on the levels (phase 2), 24 × 8 per node | RT [s] | step [s] | RT speed-up | peak RSS rank 0 / min over ranks [GB] |
|---|---|---|---|---|
| serial, 192 threads (1 node) | 2593.4 | 2878.7 | | 203.8 |
| Stage D, 2 nodes (12 × 16; one node does not fit) | 1460.7 | 1540.9 | | 292.9 / 103.0 |
| 1 node (job 41824737) | 1783.2 | 1837.5 | 1 | 293.0 / 23.8 |
| 2 nodes (job 41824738) | 974.4 | 1034.5 | 1.83 | 293.0 / 13.1 |
| 4 nodes (job 41824739) | 538.3 | 598.4 | 3.31 | 143.8 / 13.0 |
| 5 nodes, 120 ranks, the whole allocation (job 41851042, the E4 build) | 444.4 | 506.0 | 4.01 | 143.4 / 13.0 |

The fields of every node count are bitwise identical (the 5-node run
against the 4-node runs too). From 4 to 5 nodes the RT speeds up 1.21×
for 1.25× the ranks (`SumBoundary` 22.7 → 40.1 s is the part that does
not: the deposit exchange grows with the rank count); 4.0× the one-node
RT on 5 nodes, 5.7× the serial step. Phases, 1 → 2 → 4
nodes: `step1` 275.6 → 165.5 → 94.3, `cast1` 532.1 → 268.2 → 154.3
(`SumBoundary` 99.9 → 80.5 → 22.7), `phys` 277.3 → 140.1 → 70.3, `cast3`
390.2 → 221.9 → 123.7, `sort` 163.9 → 84.4 → 46.2, `migrate` 71.0 → 43.6
→ 26.2 s. So the RT scales 1.83× per doubling of the nodes (1.81× from 2
to 4), against Stage D's 1.87× from 2 to 4 nodes at 4.8× the memory per
rank; on 4 nodes the step is 598 s, 4.8× faster than the serial code on
one node. Rank 0's peak on 4 nodes is 144 GB (the restart's chunked
scatter leaves less in flight with 96 ranks), the others 13 GB.

Reading: the ranks other than 0 hold 3.1 GB with the particles on the
level (4.5 GB in phase 1: the level's share of the cells and rays,
`Luv`/`sfr`, the catalogue) against 8 GB with the host state and
12.8 GB in Stage D; rank 0's peak is the restart (the serial
ray pool and ray list read from the checkpoint, as before). The RT is 3 %
faster (no `sync_fields` broadcast after `post`), the hydro and PM are
the level steps alone. The io of 21 s is the mirror: the five-field
gather of `write_stats` and the gas gather of `summaries` every step,
and the full gather for the fields of the last step, unpacked by rank 0
in a single thread at about 200 MB/s in phase 1; with the slab loops
threaded (phase 2) the io is 5.3 s and the step 159.0 s, faster than
every host-state layout, with the hydro at 0.4 s and the PM at 0.1 s
(no particle rebuild, no gather). The per-step gathers (60 B per cell)
are the price of keeping every statistic bitwise the serial one; a
128-bit on-level path can replace them where round-off is acceptable.

## 5. Load balancing (`amrex_balance = 1`, Stage E4; branch `stage-e-balance`)

The shared decomposition of §2.1 is the RT level's, so a ray-weighted
mapping balances the RT and moves the gas and the particles with it.
Before every RT step (after the sources are assigned, whose work does not
depend on the mapping) `amrex_balance_step` prices each box at its cells
plus `amrex_balance_ray_weight` (default 10: at 512³ a ray-substep costs
about ten times a cell-substep, `docs/stage_d.md` §4.4) times its live
rays, builds a new mapping of the same BoxArray by AMReX's space-filling
curve (`amrex_balance_method = sfc`, the default, which keeps the boxes of
a rank contiguous for the exchanges) or knapsack, and compares the balance
efficiency (the mean rank cost over the largest) of the two mappings. If
the proposed one gains more than `amrex_balance_min_gain` (0.05) the state
moves: the gas MultiFab by `ParallelCopy`, the RT cells by `ParallelCopy`
and the rays by a `Redistribute` under the new mapping (their positions
are unchanged), the particles likewise, the work arrays rebuilt, the hydro
level dropped and rebuilt at its next step. The RT, hydro and PM levels
are bitwise independent of the mapping (the same argument as their rank
independence: every cross-box sum is an integer or a maximum, every merge
canonical), so a rebalanced run is byte-identical to the fixed mapping:
`ctest -L amrex_state` includes `amrex_balance_full_pair_r2` and `_r3`
(the smoke problem with 8³ boxes, remapped at every RT step with
`amrex_balance_min_gain = 0`, against the fixed mapping). The driver
prints the efficiency of every RT step and whether it remapped, and the
count at the end.

**klone, 512³** (jobs 41826700–1, the same checkpoint and layouts as
§4, `amrex_balance = 1` with the defaults; one RT step in the run, so
one remap):

| 512³, step 53, 24 × 8 per node | RT, fixed mapping [s] | RT, rebalanced [s] | gain |
|---|---|---|---|
| 2 nodes | 974.4 | 936.1 | 4 % |
| 4 nodes | 538.3 | 470.0 | 13 % |
| 5 nodes, 120 ranks (jobs 41851042–3, the whole allocation) | 444.4 | 398.1 | 10 % |

**klone, 256³** (`reion256_f0.2`, step 49, 92 RT substeps, 24 × 8 per
node, 512 boxes of 32³; jobs 41826697/41826699 rebalanced, 41840184
fixed; one remap in the one RT step):

| 256³, step 49 | RT, fixed mapping [s] | RT, rebalanced [s] | gain | step [s] |
|---|---|---|---|---|
| 2 nodes | 83.7 | 74.1 | 11 % | 93.3 → 84.1 |
| 4 nodes | 48.4 | 42.5 | 12 % | 56.7 → 50.7 |

The fields and statistics of the rebalanced runs are bitwise those of
the fixed mapping at both sizes (at 256³ the five runs, 2 and 4 nodes
with and without the balancing and the previous build, give one
checksum of the fields and the particles). At 512³ the gain grows with
the rank count, as the imbalance does (at 2 nodes 85 boxes of 32³ per
rank average out the ray clustering, at 4 nodes 43 do); at 256³ the
remap costs about the same fraction at 5 and 10 boxes per rank. The
whole `amrex` label passes on this build (161/161 on klone).

**A hang, once.** The first 256³ fixed-mapping run on 4 nodes (job
41826698, 96 ranks) hung inside the RT step, at about substep 80 of 92,
with every rank spinning until it was cancelled three hours later; its
ray log is an exact prefix of the good runs' logs, so the data were
identical up to the hang. The same binary, checkpoint and parameters on
four other nodes finished in 3 minutes, and the previous build did too,
so it is not reproducible from the data (a transport-level or hardware
event is the likeliest reading; UCX reported no error). To catch the
next one: the backend installs a SIGTERM handler, so `scancel` or the
time limit makes every rank write its call stack to `Backtrace.<rank>`
(unbuffered, `backtrace_symbols_fd`, then the signal's default: AMReX's
own handler, tried first, writes through stdio and aborts the other
ranks from the first one done, which left their files empty);
`amrex_scaling.slurm` and `reion_production.slurm` run the backend in
the run directory, where those files land; the klone build keeps line
tables (`-g1`); `tools/hyak/backtrace_lines.py <exe> Backtrace.*`
resolves the addresses and groups the ranks by stack. Tested on a 24 × 8
production leg killed by its time limit (job 41851905): 21 of the 24
ranks wrote their stack (the substep's kernels, as expected), three
files stayed empty (the signal caught those ranks where the unwinder
could not run before Slurm's kill, 10 s later). Give short runs a short
`--time`. Four more repeats of the case (jobs 41850991–41850994, two
of them on the node set of the 5-node runs) finished in 3 minutes each,
bitwise identical: nine runs of the case, one hang, nothing to read yet.

What this does not do yet: balance inside a step (rays migrate during the
substeps, so the cost of a box changes while the step runs), weigh the
hydro and PM work (small against the RT), or split boxes (the granularity
is the box; at 256³ on 4 nodes, 512 boxes of 32³ over 96 ranks, the
balance is limited by 5 boxes per rank and finer boxes would be the next
lever).

## 6. Production runs on the backend

`tools/production/make_production.py --backend amrex` writes the backend
configuration into `params.txt` (`rt_backend`, `hydro_backend`,
`pm_backend = amrex`, `amrex_max_grid_size = 32`,
`amrex_hydro_max_grid_size = 0`, `amrex_state = 1`, `amrex_balance = 1`,
`deterministic_parallel = 0`; it needs `--hydro-scheme unsplit` and no
recombination-radiation transport), so that every restart runs the same
configuration. `tools/hyak/reion_production.slurm` sees `rt_backend =
amrex` there and launches the ranks with `srun` (24 × 8 per node, the
module and UCX settings of `amrex_scaling.slurm`, the run directory as
the working directory); its restart, `DONE` and post-processing logic is
unchanged, since rank 0 writes every file the serial code writes
(`docs/cluster.md` §5).

**Validation (klone, 2026-10-13):** a 128³ production problem at 0.1c
(`make_production.py --N 128 --c-reduced 0.1 --f-esc 0.5`, the
abundance-matched sources, movie frames every Δz = 0.05, a checkpoint
every 4 steps) run to z = 5 from the same ICs by the serial code and by
the backend:

| 128³, 31 steps to z = 5, one node | wall [s] | RT [s] | hydro [s] | io [s] | peak RSS [GB] | z(x_HII = 0.5) | τ_es |
|---|---|---|---|---|---|---|---|
| serial, 192 threads (`prod128-serial`, job 41851096) | 605 | 534 | 35 | 22 | 5.6 | 7.831 | 0.06104 |
| backend, 24 × 8 (`prod128-amrex`, job 41851097) | 484 | 442 | 1.4 | 7.7 | 3.2 (rank 0) | 7.831 | 0.06104 |

The histories (`reion.txt`: x_HII by volume and mass, τ_es at every
frame) agree to 1 in the fifth digit; the fields differ cell by cell
from the first sourced step on, as the two ray orders do
(`docs/stage_d.md` §4: the backend is not the serial continuation), with
the global statistics (v_rms, c_max, the DM moments) to 2 × 10⁻⁵. The
job wrote the checkpoints, `DONE`, the stills and the movie through the
unchanged post-processing. **Restart:** the same problem under a
4-minute time limit (job 41851312, killed in step 25 by the limit's
SIGTERM) and resubmitted (job 41851314, resumed from the step-24
checkpoint, 5 min): every field of the last step, its statistics and
the final checkpoint are bitwise those of the uninterrupted run.

## 3. Phases and validation

| Phase | Work | Test |
|---|---|---|
| 1 | shared decomposition, `GasState`, the hydro and RT levels on it, the PM level reading and writing the gas on it, `summarize_gas` on the level, the host gas arrays gone, fields and stats output through the slab gather | `radhydro --set amrex_state=1` byte-identical to `amrex_state=0` on the 64³ smoke problems (fields, stats, frames) for 1/2/3 ranks × 1/2 threads; the fast suite and `ctest -L amrex` unchanged |
| 2 | particles persistent, drifts and `Redistribute` on the level, the halo gather, sparse `Luv`/`sfr` and metals, P(k) on the level, particle output by id | the same byte-identity (P(k) to round-off); the pm_refine = 2 loop |
| 3 | checkpoint and restart through the gathers, chunked rays | bit-identical restart (`restart_smoke` with `amrex_state = 1`); a checkpoint written with `amrex_state = 1` restarts the `amrex_state = 0` loop and the serial loop byte-identically |
| 4 | klone: 512³ on one node (`reion512_40mpc_f0.23`), 1–4 nodes; memory per rank | `docs/stage_d.md` §4.4 extended |

Why byte-identity against the current backends is the right test: every
per-cell and per-particle operation is the same function on the same
values (the levels already are bitwise identical to the serial step in
their arithmetic), the maxima are exact reductions, and the only
order-dependent sums that reach a file are computed on rank 0 in the
serial order. The sums that become 128-bit fixed point (`summarize_gas`,
`summarize_dm`, P(k)) are diagnostic; they are listed and compared to
round-off.

## 4. Open questions (for the user)

- Whether a distributed FoF is wanted before Stage F, or the rank-0
  gather is enough (it is at 512³ with Ndm = 512³: 1.6 GB and the FoF's
  own time).
- Whether the on-disk checkpoint should move to AMReX's native
  checkpoint (parallel, one file per rank) once the format no longer has
  to be shared with the serial code; version 7 is kept here so that a
  run can move between the serial code and the backend.
