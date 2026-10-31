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

## 5. Step 11 (2026-10-30, branch `pk-on-level`): the power spectrum on the levels

`amrex_pk = level` (`config.hpp`; default `host`, the previous behaviour;
`make_production.py --backend amrex` writes `level`) computes the `Pk/`
files of the data1 outputs on the levels (`amrex_pm_power_spectrum`,
`pm_level.cpp`):

- **Dark matter:** level 0's CIC deposit of the DM particles alone (no
  stars), exact in fixed point as the gravity's, then δ = ρ/(1 − f_b) − 1
  on the PM mesh.
- **Gas:** ρ/f_b − 1 from the gas state on the gas grid.
- **Transform:** each field goes through a pencil FFT of its own size
  (`PencilFFT::forward_visit`: the same 1-D FFTW line transforms whatever
  the layout).
- **Binning:** `pk_bins_level` uses the per-mode values and binning of
  `calc_pk_3d`. The bin sums are exact: 128-bit fixed point, with a scale
  per bin and column taken from a first pass (the bin's count and largest
  value), so two forward transforms per field. The normalisation
  `pm::calc_pk_finish` is shared with `calc_pk_3d`.
- The result is rank- and thread-independent. Rank 0 then allocates only
  the particle array of the host PM arrays (no density mesh, linked list,
  head-of-chain arrays or FFT plan): about 70–100 GB less on rank 0 at
  1024³ with `pm_refine` 2.
- **Not with refined PM levels:** with `pm_level_solve`, particles live on
  refined levels, and `amrex_pk = level` refuses to run.
- **Against the host spectrum:** the 16³ smoke problem, `pm_refine` 1 and
  2. The DM columns agree within 1.1e-6 relative; the gas columns are
  identical. Every other file is byte-identical. `Pk/` is byte-identical on
  1 and on 2–3 ranks (`tests/pk_level_check.py`, tests `amrex_pk_level_*`).
- Peak memory: 256³ gas, `pm_refine` 2 (a 512³ PM mesh), OTVET cone,
  4 ranks × 8 threads, 12 steps with two power spectra.

  | Peak RSS | `amrex_pk = host` | `amrex_pk = level` |
  |---|---|---|
  | rank 0 | 10.53 GB | 9.13 GB |
  | ranks 1–3 | 7.46–7.57 GB | 7.66–7.77 GB |

  Rank 0 saves the host density mesh and linked list of the 512³ mesh
  (1.4 GB; at 1024³ with `pm_refine` 2, the 2048³ mesh, about 80–110 GB).
  The other ranks gain about 0.2 GB: the pencil FFT's buffers, freed after
  each spectrum (about 20 B per mesh cell over the ranks, 2.7 GB per rank
  at a 2048³ mesh on 64 ranks).

## 6. Status (2026-10-29): steps 1–10 done

- Primitives (`amrex_state.hpp`):
  - `amrex_state_stream_rt` (a component range or any list of the
    `rrk::CellData` floats) and `amrex_state_stream_gas` bring the level to
    rank 0 slab by slab as a `StreamSlab`;
  - `amrex_state_fill_gas` / `amrex_state_fill_rt` are the reverse;
  - `amrex_gamma_avg_stream` streams `gammaHI_avg`.
- Consumers (`apps/radhydro.cpp`, `src/checkpoint.cpp`):
  - `GridSums` (write_stats and the ĉ schedule);
  - `GasSummer` (summarize_gas, split into begin / `add_planes` / finish;
    the serial path calls it over the host arrays);
  - `write_fields_streamed`, `write_frame_streamed` (one pass over the
    components a frame needs, so a default frame moves 4 floats per cell,
    against 6 for the old gather);
  - `write_data_files_streamed`;
  - the power spectrum's gas density;
  - the checkpoint's grid blocks (`stream_grids_a` / `_b`).
- Restarts and initial conditions:
  - `read_checkpoint` skips the grid blocks on every rank, and
    `read_checkpoint_grids` reads them on rank 0 slab by slab onto the
    levels;
  - `u_gas.init` goes onto the gas state the same way, with the per-cell
    `pm::ics_gas_cell` and `unsplit_init_e` shared with the serial path.
- Rank 0 allocates neither the RT nor the gas grid. The zero blocks
  (`htoc`, `v_gas`) are written as zeros, never stored.
- Tests on klone:
  - `amrex_state` 26/26 byte-identical: the pairs against the host state
    (fields, stats, movie frames, snapshots, checkpoints), the restarts, and
    the moment RT against the serial RT;
  - the `amrex` label 321/321 (not nightly);
  - the serial fast suite 163/163 (without `platform_ref`).
- Peak memory: 256³, OTVET cone, 4 ranks × 8 threads, 12 steps with a
  snapshot and checkpoints. The old executable is the 1024³ cone run's.

  | Peak RSS | before | after |
  |---|---|---|
  | rank 0 | 11.45 GB | 7.79 GB |
  | ranks 1–3 | 5.10–5.33 GB | 4.54–4.71 GB |

  Rank 0 saves 3.7 GB, the 228 B per cell of §1 (3.8 GB at 256³); at 1024³
  that is about 245 GB. Rank 0's remaining excess over the other ranks
  (about 3 GB here) is the particle copy of the halo finder and the host
  PM arrays of the power spectrum (step 11).
