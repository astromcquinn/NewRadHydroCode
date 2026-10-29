# Rank 0 without the full grid copy (`amrex_state = 1`)

*Branch `rank0-no-mirror`, 2026-10-29. User decisions in `docs/PLAN.md`,
"Rank 0 without the full grid copy".*

## 1. The problem

With the state on the levels (`amrex_state = 1`, `docs/stage_e_state.md`),
rank 0 keeps a full copy of the grid, the "mirror", for every output that
needs the whole box: the per-step stats, the gas summary, the field files,
the movie, the snapshots, the checkpoint and the restart. Per cell it holds:

| Item | Bytes per cell |
|---|---|
| RT fields (rhog … tau, 36 floats) | 144 |
| `htoc` (ray chain heads, 4 × int64; zero for the moment RT) | 32 |
| gas u (5), w (4), e (1) | 40 |
| grid velocity `hs.v` (3; zero with the unsplit hydro) | 12 |
| **total** | **228** |

That is 245 GB at 1024³ and about 2 TB at 2048³. The 1024³ cone run
(`docs/stage_e_moments.md` §10) peaked at 458 GB on rank 0 against about
57 GB on the other ranks. The rest of rank 0's excess is the host PM arrays
of the power spectrum (density and linked list on the 2048³ PM mesh) and
the particle copy of the halo finder.

## 2. The approach: stream, don't mirror

Every consumer makes one pass over the cells in the serial order
(c = i + N j + N² k), either summing in that order or writing a file in
that order. `gs::gather_slabs` already brings the level to rank 0 in
slabs of 8 k-planes, in k order. Each consumer becomes a function of a run
of consecutive planes, and it is called:

- in the serial code: once, over the host arrays (the same arithmetic as
  now);
- on the level: once per slab, as the slabs arrive on rank 0.

The sums are sequential in the same order, and the files are written in
the same order, so every output stays byte-identical. Rank 0's transient
is one slab: 8 N² cells of the components being streamed.

Each field is streamed separately, in the order the file layouts need:
one array per field in the checkpoint and the step fields, interleaved per
cell for `Gas/u_gas` (5) and `RT/rt_grid` (36). The checkpoint's zero
blocks (`htoc` and `hs.v` with the moment RT and the unsplit hydro) are
written as zeros without being stored.

## 3. The work, in order

| Step | Consumer | Change |
|---|---|---|
| 1 | streaming primitives | stream components of the RT cells and of the gas state to rank 0 (`amrex_state.hpp`) |
| 2 | `write_stats`, `update_sched_x` | sums over streamed planes |
| 3 | `summarize_gas` | the per-plane body shared by the serial and streamed paths |
| 4 | `write_fields` | one streamed pass per field file |
| 5 | movie frames | slabs or full cubes written per streamed slab. The ionized-cell Γ_HI values for the median are still collected on rank 0 (at most 4 B per cell) |
| 6 | `write_data_files` | `Gas/*` and `RT/rt_grid` streamed; `DM/xv_dm` from the particle copy (kept) |
| 7 | power spectrum | the gas density streamed into its FFT buffer. The DM part keeps the host PM mesh until step 11 |
| 8 | checkpoint write | the grid blocks streamed; zero blocks written as zeros |
| 9 | checkpoint read, ICs | rank 0 reads slab by slab and scatters (`gs::scatter_slabs`) |
| 10 | rank 0 allocates no grid | `RTGrid::allocate_none`, empty gas arrays; then find any remaining use |
| 11 | (stage 2) power spectrum on the level | the DM deposit and FFT on the PM level; then rank 0 frees the host PM mesh, keeping only the particles |

## 4. Validation

- Every output byte-identical:
  - the existing pairs of host state (`amrex_state = 0`) against level
    state (`amrex_state = 1`): fields, stats, movie frames, snapshots and
    checkpoints;
  - the restart tests;
  - the moment-RT pairs against the serial RT.
- The serial path changes only in structure. The full suite, including the
  stored Fortran references, must stay byte-identical.
- Rank 0's peak memory measured before and after, at 256³ on klone.
