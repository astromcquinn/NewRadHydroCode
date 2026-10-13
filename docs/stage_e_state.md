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
| `amrex_state = 1`: one gas BoxArray and DistributionMapping shared by the RT, hydro and PM levels; the gas state as a persistent MultiFab; the hydro and RT levels read and write it locally; the per-rank host gas arrays dropped | planned, phase 1 |
| The DM particles persistent in the PM container (drift and `Redistribute` on the level, no rebuild from the host); the halo finder on a gather of the particle positions to rank 0; `Luv`/`sfr` and the halo metals as sparse per-source scatters | planned, phase 2 |
| Rank-0 I/O through slab gathers: fields, movie frames and slabs, data files, checkpoints in the existing formats; chunked ray gather (the 2 GB `MPI_Gatherv` limit); restart by scatter | planned, phase 3 |
| Statistics: maxima by exact reduction (the timestep uses only `cmax_gas` and `amax_dm`), the order-dependent serial sums of `write_stats` and `write_frame` computed on rank 0 from the gathered fields (bitwise the serial order), the per-step `summarize_gas` sums as 128-bit fixed point (diagnostic only) | planned, phase 1 and 3 |
| Validation: the full loop with `amrex_state = 1` byte-identical to the same run with `amrex_state = 0` (the three backends with the host state), for 1/2/3 ranks and 1/2 threads; bit-identical restarts; the existing round-off comparison against the serial code carries over | planned |
| klone: 512³ on one node, 1–4 nodes scaling | phase 4 |

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
