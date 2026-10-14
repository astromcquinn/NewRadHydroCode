# The moment RT on AMReX (`rt_backend = amrex` with `rt_method = m1 | otvet`)

*Design, 2026-10-13, branch `moments-amrex` (from `main` at `e0b8591`,
after PR #36). Not implemented yet: this document is the plan for the
user's review. Goal: run the moment RT (M1, thin OTVET, `otvet_split`,
`otvet_closure = screen | cone`) across nodes, bitwise identical to the
serial moment RT, then a 1024³, 40 Mpc/h cone run on four `cpu-g2` nodes.*

## 0. Summary

- **Where it plugs in.** `rt_backend = amrex` now refuses anything but
  `rt_method = rays` (`src/rt_backend.cpp:79`, `src/amrex_state.cpp:38`).
  The plan lifts that: with `rt_method = m1 | otvet` the backend runs a
  new `MomentLevel` (`src/amrex/moment_level.cpp`), beside the ray level,
  not inside it. It reuses the per-cell kernels of `rrk::` (`prep_cell`,
  `calc_rt_1_cell`, `physics_cell`, `post_cell`) and the gas state of
  Stage E3. The serial path is untouched, so its defaults stay
  byte-identical.
- **The decomposition is z-slabs.** Each box is whole k-planes (N × N ×
  nz) and each rank holds a contiguous run of planes. That one choice
  makes every order-dependent sum of the serial code reproducible
  exactly (§2). With `amrex_state = 1` the hydro and PM levels share the
  slabs. The serial code sums plane by plane, and only whole planes on
  one rank reproduce that order without moving the addends.
- **Communication per RT substep:** one ghost layer per GLF subcycle,
  two faces per slab, and almost all of it inside a node. N doubles per
  plane-ordered sum.
- **Communication per hydro step:** the OTVET tensor. Its FFTs run on
  "component ranks" as the same serial FFTW plans as today, so the
  result is bitwise the serial one. No gather of the grid is needed:
  every rank already knows the sources.
- **The one sum that cannot be both exact and fast** is the
  photochemistry statistics `rt_stat`. In the serial code it is one
  sequential sum over every cell of the grid; it is diagnostic only.
  Proposal: `amrex_moment_stats = exact | fixed`.
  - `exact`, the default: a recorded, chained replay, bitwise the serial
    file. Its cost is O(N³) sequential per substep, which is fine up to
    256³.
  - `fixed`: 128-bit fixed point, as Stage D does, for production.
  - The fields never depend on this flag.

## 1. What the serial moment RT does, and what each piece needs

| Piece | Serial code | Data dependence | Per |
|---|---|---|---|
| `calc_rt_1` (tau, reset of the RT fields) | `moment_step_1` | per cell | substep |
| injection (sources; `recomb_rad` photons) | `moment_step_1`, `inject` | per cell (source cells) | substep |
| GLF transport, `moment_nsub` (4) subcycles per group | `transport_group` | Jacobi, 7-point stencil: 1 ghost layer of (N, F, E) of the group; P = D N on the grown box needs h, ndir there too | subcycle × group |
| deposit ngam, nu, dx | `moment_step_1` | per cell | substep |
| photochemistry | `rt_step_2` → `radiative_physics_plane` → `physics_cell` | per cell (N_rt_gas = 1) | substep |
| absorption N, F, E *= fgam | `moment_step_3` | per cell | substep |
| redshift of the group spectra | `moment_redshift` | per cell | hydro step |
| OTVET thin tensor | `compute_otvet_tensor` (direct sum or FFT) | all sources → every cell | hydro step |
| split, screen, cone | `otvet_closure_tensor` | FFT convolutions; cone: pyramids, march, budget, blend | hydro step |

Memory of the serial moment state: 25 floats per cell (N, F, E for 5
groups, 100 B), the transport scratch (56 B), h (24 B), ndir (4 B).

## 2. The sums: the reason for slabs

Every sum below was read in the code. "Feeds the fields" means a
different rounding would change the radiation or the gas, not only a
diagnostic.

| Sum | Serial order | Feeds the fields? | On the slabs |
|---|---|---|---|
| `injected` (sources) | sequential over source cells in cell order | no (log, checkpoint) | each rank's source addends (5 per source cell) gathered to rank 0 in slab order = cell order; exact |
| `absorbed` (`moment_step_3`), `lost` (redshift), `recomb` injection | `plane_sum`: sequential within a plane, then over planes | no | the plane partials are local; `MPI_Allgatherv` of N doubles, summed in plane order; exact |
| `ngam_rays_total` (stats.txt "in flight") | sequential over all cells | no | a chain: rank r adds its contiguous cells to the running sum of rank r−1; exact. ~1 s per call at 1024³, once per hydro step |
| `rt_stat` (photochemistry: ngam, nion, ncol, nrec, nrab, Γ, Λ, dtau) | sequential over all cells; per cell a fixed sequence of addends (`RPSink`) | no (stats.txt, rt_stats.txt, log, checkpoint) | `exact`: per-cell `CellPhys` records, replayed in a chain (the next rank waits); `fixed`: 128-bit fixed point (the Stage D helpers) |
| any photons in group g (skip) | OR | yes | exact OR reduce |
| OTVET `Lmax`, `Lsum`, `qw` | max / sequential over the source list | yes | every rank holds the source list in cell order |
| `trmin` = 1e-6 Σ rho | `plane_sum` over a sparse array | yes (screen, cone) | from the source list, plane by plane; adding the zeros is exact, so skipping them is too |
| screen λ from Σ nw | `plane_sum` | **yes** | plane partials; exact |
| cone `target_ion`, `near_tr` | `plane_sum` | **yes** (λ_x) | plane partials; exact |
| cone `tb[b]` (per-bin light) | sequential over receivers r | **yes** (λ_x) | a chain over the receiver planes (slab thickness a multiple of `otvet_cone_coarse`); exact, ~0.3 s at 1024³ |
| split and screen `ratio` | `plane_sum` | no | plane partials |

**Why not the 32³ boxes of Stage E3:** a plane split over several ranks
cannot reproduce a sequential sum without moving all its addends to one
place. `absorbed` alone is 5 doubles per cell per substep, and the screen
and cone budgets are 8 B per cell per hydro step. Transposing these to
slabs every substep would cost more than the transport's own halo
exchange. Slabs avoid it, and they suit the rest:

- the transport exchanges 2 large faces per slab, and only the two
  slabs at a node boundary cross the network;
- the cone's receivers and pyramid levels are slab-local;
- the photochemistry's per-cell cost varies, but whole planes average it
  well.

The cost: at most N / nz ranks. 1024 planes with nz ≥ 8 allow up to 128
ranks, for example 4 nodes × 16 ranks × 12 threads with nz = 16.

**Parameters:**

- `amrex_slab = 1` (new): z-slabs of `amrex_slab_planes` planes (0 = N /
  ranks), a multiple of 2, and of `otvet_cone_coarse` for the cone.
- The slabs are mapped to ranks in z order by a hand-built
  `DistributionMapping`; the chains need that order.
- `rt_method = m1 | otvet` with `rt_backend = amrex` requires
  `amrex_slab = 1`.
- With `amrex_state = 1` the shared decomposition of E3 becomes the
  slabs: `amrex_state_init` builds them, the hydro and PM levels take
  them, and `pm_refine` refines them as it refines cubes. The hydro
  needs nz ≥ its ghost width (checked).
- Stage E4's `amrex_balance` (PR #37, in `main`) is refused with
  moments for now and stays off for moments. A plane-contiguous rebalance (variable nz from
  measured cost) is possible later without breaking exactness.

## 3. The per-substep pieces on the level

`MomentLevel` holds:

- `C`: CellData, 46 floats, no ghosts, as Stage D. The deposits `A` and
  `TF` are not needed: there are no rays.
- `M`: 25 floats per cell, the moment state (component g·5 + {N, Fx,
  Fy, Fz, E}), no ghosts.
- `W`: the scratch of one group, 5 components with 1 ghost, plus P (9
  components) on the grown box and the new state (5).
- `H`: h (6) and ndir (1), 1 ghost, filled once per hydro step.

Step by step:

1. **`prep`** (per hydro step): `prep_cell` from the gas state (E3) and
   the sources, as the ray level does. With `amrex_state = 0` (the tests
   and `rt_driver`), copy in from the host grids as Stage D does.
2. **`calc_rt_1`**: `calc_rt_1_cell` per cell (tau, reset).
3. **Injection**: per source cell. The addends go to rank 0 for
   `injected` (§2). Both `moment_inject = 0 | 1` are kept.
4. **Transport**:
   - Per group with photons (the exact OR): copy the group's 5 fields
     into `W`, then `FillBoundary` for 1 ghost (periodic; in x and y it
     is a copy within the slab).
   - Compute P on the grown box with the serial's
     `m1_pressure`/`otvet_split_pressure` (moved to a kernel header), and
     the 7-point update with the same expression order, including the
     zero-neighbour shortcut and the realizability clamp. Copy back.
   - Each output cell is the same function of the same 7 inputs as in
     `transport_group`, so the result is bitwise the same for any slab
     thickness and thread count.
   - Messages per subcycle and group: 2 faces × N² × 20 B. At 1024³
     with nsub = 4 and 5 groups, that is 0.84 GB per rank per substep,
     inside the node except at the two node faces (0.42 GB per node per
     substep across the network).
   - Possible later: the groups are independent, so all 5 could share
     one exchange.
5. **Deposit**: ngam, nu and dx into `C`, per cell, with the serial
   expressions.
6. **Photochemistry**: the ray level's step 2 loop (`physics_cell` per
   cell, fgam out) without the ray decoding. `rt_stat` is accumulated as
   §2 sets.
7. **`step_3`**: N, F, E *= fgam per cell; `absorbed` by plane partials.
8. **`post`** (per hydro step):
   - The redshift per cell; `lost` by plane partials.
   - `ngam_rays_total` by the chain.
   - `post_cell`, and the gas written back (E3) or `sync_fields`
     (`amrex_state = 0`).

`tensor_dump` is refused at first. `recomb_rad = 1` (the user wants it
eventually) is phase 3b, §9. It needs:
- the recombination buffers `ngam`, `egam` per cell and group on the
  level, filled by the photochemistry;
- the injection into the moments per cell, with `launched` in plane
  partials;
- the buffers in the checkpoint through the slabs.

## 4. The OTVET tensor (once per hydro step)

The tensor needs the sources (every rank already has the catalogue and
the source list) and, for the closures, three per-cell fields of the
solution: nw (Γ-weighted photons), col (HI column) and x_HII. All three
stay on the slabs.

### 4.1 The FFT convolutions: component ranks, no grid gather

Today's `otvet_fft_convolve_kernel` already works as follows:

- one serial forward FFT of rho;
- six cached kernel spectra;
- six independent products and serial inverse FFTs, one per component.

Each inverse FFT is the same single-threaded FFTW plan on `fftwf_malloc`
arrays; that is what makes it thread-independent. The plan:

1. **Split the function into its three pieces**, and have the serial
   code call them. Same arithmetic, byte-identical; checked by the
   existing tensor tests.
2. **Assign each component (m = 0..5) of each convolution to a
   "component rank".** The thin tensor and the screen use one
   convolution (6 tasks). The cone uses two, the blend's thin and the
   near field (12 tasks). The tasks are spread over the nodes for memory
   and bandwidth.
3. **On each component rank:**
   - build rho from the source list, the serial `source_density`;
   - run the forward FFT itself (redundant, but no 4 GB broadcast);
   - take the kernel spectrum of its component from its own cache
     (filled on its first use: the thin and near kernels are fixed);
   - form the product and the inverse FFT.
   - Its output, `out[m]` (N³ floats), is scattered onto the slabs (one
     component of `P`) by `ParallelCopy` from a single-box FabArray on
     that rank.
   - The division by N³ and the conversion to double happen on the slab,
     in the same expression.
   - Memory per task at 1024³: r, rk, its spectrum, the scratch and out,
     about 22 GB.
4. **Normalisation on the slabs**, per cell: `normalise`,
   `normalise_safe`, ndir. Exact.
5. **Screen with the photon budget:** λ changes every step, so the
   kernel is not cached. Its fill (27 images × exp per cell) is too slow
   on one thread at 1024³ (about 9 min per component). Instead the
   kernel is filled on the slabs, a pure function of the separation,
   and `ParallelCopy`'d to the component ranks. That is 4 GB per
   component.
6. **Measured on klone** (2026-10-13, one `cpu-g2` node, FFTW 3.3.10
   float, a standalone benchmark in
   `/gscratch/astro/mcquinn/radhydro/scratch/fftbench`): one 3-D
   transform of the serial code's plan, and the alternatives.

   | N | serial plan, 1 thread (r2c / c2r) | the same 3-D plan with FFTW threads | bitwise = serial? | line-by-line (as the PM pencil FFT), one node, no transposes | its difference from the serial plan |
   |---|---|---|---|---|---|
   | 256 | 0.24 / 0.23 s | 0.03 / 0.03 s (8 threads) | yes, both | 0.02 s | 2.9e-9 of the max |
   | 512 | 2.69 / 2.62 s | 0.09 / 0.09 s (32 threads) | yes, both | 0.06 s | 1.1e-9 |
   | 1024 | 21.8 / 21.6 s | 0.70 / 0.72 s (32 threads) | yes, both | 0.48 s | 3.6e-10 |

   - **Decision: the component ranks run the serial 3-D plans with
     FFTW's threads** (the rank's `SLURM_CPUS_PER_TASK`). At every size
     and thread count tried they are bitwise the single-threaded plan.
     FFTW does not promise that, so the backend checks it: once per run,
     at the first tensor, one task compares the threaded and the
     single-threaded forward transform of its rho. It throws on a
     mismatch, with a key to fall back to one thread. A ctest does the
     same at 64³ and 128³.
   - **Cost at 1024³:** about 0.7 s per transform on 32 threads, about
     2 s on 12. That is a few seconds per convolution plus the transfers
     (4.3 GB per task output, 26–52 GB per step in all), so roughly
     5–10 s per hydro step for the cone. The single-threaded plans would
     take about 45 s per convolution.
7. **Not taken: a distributed FFT** (the PM's pencil transform, or
   `amrex::FFT`).
   - The line work alone is 0.48 s per node at 1024³, so over 4 nodes
     about 0.1 s. Its two all-to-all transposes (4.3 GB of complex data
     each, three quarters of it across the network) take an estimated
     1–2 s per transform.
   - So it would be no faster than the threaded component ranks, and it
     differs from the serial plan at about 1e-9 of the maximum.
   - It stays an option if the tensor's transfers ever dominate.
8. Side result: the serial code could thread its FFTs the same way,
   byte-identically. This is not part of this thread.

The thin direct sum (`otvet_tensor = direct`, or auto when
N_src × N³ ≤ 2.7·10⁸) is per cell over the source list on the slabs:
exact, no FFT.

### 4.2 The cone

| Step | On the slabs |
|---|---|
| h_s (thin, split by β = min(1, q_w ndir / N_w)) | the thin convolution of §4.1, then per cell |
| near field P_n | the near-kernel convolution of §4.1 (6 more components) |
| pyramids of rho and of the capped column | level 1 is built locally from 2×2×2 blocks (nz even). `MPI_Allgatherv` of level 1 (N³/8 floats per pyramid: 0.54 GB each at 1024³), then levels ≥ 2 redundantly on every rank, with the same arithmetic. The march samples only levels ≥ l_min = ⌊log₂ max(2 r_c tanθ, 1)⌋; with the defaults (r_c = 8, nside = 4) l_min = 1, so level 0 (the full grid) is never replicated. If l_min = 0, level 0 is allgathered too (with a memory warning). |
| the march (12 nside² pixels per receiver, out to 8 boxes) | each rank marches the receivers whose coarse block lies in its slab (nz a multiple of `otvet_cone_coarse`). Each receiver's result is independent of the others: exact. |
| budget: `target_ion`, `near_tr` | plane partials (§2) |
| budget: `tb[b]` over receivers in order | a chain in z order (§2), then the bisection for λ_x redundantly on every rank |
| `Pr` = Σ_b fac_b P_b and the trilinear prolongation | per receiver, then per cell with one ghost plane of receivers (ncr² × 6 doubles) |
| blend: Gaussian of σ = 4 cells (R = 16) per component, axes x, y, z | x and y are local to the slab. For z, each component is `ParallelCopy`'d to y-slabs (N × N/P × N), blurred along z, and copied back: 2 × 4 B per cell per component (52 GB per step at 1024³). A 16-plane ghost would move 4× the slab. |
| clamp and normalise | per cell |

The cone's march was "a few seconds" at 512³ serially. It splits over
all ranks and stays small.

## 5. Checkpoints

The format stays version 7. The moment block is 5 vectors of
kNfreq · N³ floats in group-major order, then the budgets.

- **Write:** rank 0 writes the vector length, then streams group by
  group and slab by slab, through the slab `ParallelCopy` of E3
  (`gather_slabs`, 8 planes at a time). The budgets are on rank 0. No
  rank ever holds the full moment state.
- **Read:** today every rank reads the moment vectors in full
  (`src/checkpoint.cpp:374–392`). The plan: rank 0 reads them in slabs
  and scatters them onto `M` (`scatter_slabs`); the other ranks skip the
  block.
- h and ndir are not checkpointed, as today: the first substep after a
  restart recomputes them (`M.h` empty).
- So a checkpoint moves freely between the serial code and the backend,
  as for E3.

## 6. What every rank holds at 1024³, and what has to change first

Per local cell, with `amrex_state = 1`:

| Item | Size per cell |
|---|---|
| moment state | 100 B |
| transport scratch | ~80 B |
| h, ndir (1 ghost) | 28 B |
| CellData | 184 B |
| gas | 40 B |
| hydro level | ~330 B |
| PM | ~30 B, plus 32 B per particle |
| **total** | about 820 B, so 0.9 TB over the job at 1024³ |

On top of that:

- **Rank 0's mirror** (E3): about 230 B per cell plus 24 B per particle,
  270 GB.
- **The tensor tasks:** about 22 GB each.

**Two E3 leftovers that are N³ on every rank.** Both must go before
1024³:

- the full `Luv` and `sfr` arrays (8 B per cell, 8.6 GB per rank);
- the transient `per_cell` lists of `calc_halo` (24 B per cell, 26 GB
  per rank).

With 16 ranks per node that would be 550 GB per node. Node 0 would then
hold the mirror (270 GB), its share of the levels (220 GB) and the cone
tasks: more than its 1.5 TB at 32 ranks, and too close at 16.

**Done (phase 0), simpler than first planned.** Every rank already
receives the halo catalogue (`bcast_catalog`), so no rank-0 step is
needed:

- `calc_halo` builds a sparse list of the source cells in ascending cell
  order (`RTGrid::src_cell`, `src_Luv`, `src_sfr`, 16 B per source
  cell) with a stable sort of the sources by cell. The per-cell sums keep
  the original input order. `rt_stat.Luv` and `rt_stat.sfr` add the
  source cells in cell order; the empty cells added 0, which changes no
  sum. So the serial results are byte-identical, and the N³ `per_cell`
  lists are gone on every backend.
- The dense `Luv` and `sfr` are filled from the list where they exist: on
  every rank of the serial backends and on rank 0's mirror. With
  `amrex_state = 1` the other ranks allocate no grid field at all
  (`RTGrid::allocate_none`), and a restart skips their `Luv` and `sfr`
  blocks.
- The RT level's prep sets each local cell's Luv and sfr from the list
  (binary search on the box's k-planes).
- The same list is the global source list the OTVET tensor needs.

## 7. Validation

The references and the tests:

- **RT only, the cleanest test:** `rt_driver` with `rt_backend = amrex,
  amrex_slab = 1`, against the serial `rt_driver` with the same
  `rt_method`.
  - Every file (fields, `stats`, `rt_stats`) is compared with
    `compare.py --rtol 0 --atol-frac 0`, with `amrex_moment_stats =
    exact`.
  - Cases: M1, thin OTVET (direct and FFT), split, screen with the
    budget, cone with the defaults.
  - Each case runs on 1, 2 and 4 ranks × 1 and 2 threads, and with two
    slab thicknesses.
  - The moment problems of `tests/problems*`, with the existing
    `run_amrex.cmake` pattern. Label `amrex;moments`.
- **Full loop:** `radhydro` with `rt_backend = amrex, amrex_state = 1`
  against `rt_backend = serial` with `hydro_backend = amrex, pm_backend =
  amrex, amrex_state = 0`.
  - The PM backend is round-off different from the serial PM, so the
    reference has to use it too. The serial moment RT runs on rank 0's
    host grids.
  - Smoke problem, 1/2/3 ranks, `backend_pair_check.py`.
- **Restart:** `restart_check.py` with `amrex_state = 1` and moments,
  for M1 and the cone. Also, a serial checkpoint restarted on the
  backend, and the reverse.
- **Fixed mode:** `amrex_moment_stats = fixed` gives bitwise the same
  fields as `exact` (a test). `stats.txt` agrees to round-off in the
  `rt_stat` columns only.
- **klone, 1, 2 and 4 nodes:**
  - 256³ (`otvcl256_*` checkpoints): one hydro step per closure from a
    restart. Every file bitwise against the serial restart step, in
    `exact` mode.
  - 512³, in `fixed` mode (§8): the fields bitwise against a serial
    512³ restart step of `reion512_rays_otvet_cone`.

## 8. Measurements, then the 1024³ run

- **512³ timing** (`reion512_40mpc_f0.23` ICs or its checkpoint, cone):
  one hydro step from a restart on 1, 2 and 4 nodes.
  - Timers: transport, halo exchange, photochemistry, tensor (FFT
    tasks, transfers, march, blend), sums.
  - Peak RSS per rank.
  - The rank × thread layout (e.g. 16 × 12, 8 × 24) and nz chosen from
    the measurements.
- **1024³, 40 Mpc/h, cone defaults, 4 nodes, to z = 5:**
  - movie slices (`movie_slab`, `movie_dz` ≤ 0.05);
  - ICs from `tools/ics` at 1024³;
  - the same f_esc as 512³ unless the user says otherwise;
  - on the preemptible `ckpt-g2` partition (user decision 2026-10-13):
    `PreemptMode = REQUEUE` with no grace time, so a preemption kills the
    job at once. The job script therefore restarts from the newest
    complete checkpoint. With `restart_every` = 2–3 hydro steps (about
    1 h of wall time at 1024³) and the two newest checkpoints kept. At 1024³ a checkpoint is about 0.4 TB,
    written by rank 0 through the slab gathers; its write time is
    measured at 512³ first. If it is too slow, a per-rank checkpoint
    is a separate decision.
  - Time estimate: 16 × (2.0 + 1.7) h of transport and photochemistry
    on one node ≈ 59 h, so about 16–20 h on 4 nodes at 75–90 %
    efficiency, plus the tensor (1–2 h) and the hydro, PM, FoF and I/O
    (~1–2 h).
  - Nodes: every `ckpt-g2` node has at least 1.5 TB (sinfo
    2026-10-14: 68 × 192-core cpu-g2, the rest GPU nodes with 128
    cores); request `--constraint` on the 192-core nodes.
  - Disk, scaled ×8 from the 512³ cone run (Gas 23 GB, movie 38 GB,
    restart 2 × 50 GB): about 190 GB of snapshots, 150 GB of movie
    slices (×4: N² per frame) and 800 GB of checkpoints, about 1.2 TB.
    `/gscratch/astro` had 2.0 TB free (95 % full) on 2026-10-14 and
    `/gscratch/scrubbed` is over its quota.
  - Then the comparison with the 512³ cone: `gamma_power.py`,
    `reion_power.py`, `rt_methods_runs.py`, `igm_stats.py`.

## 9. Phases

| Phase | Work | Test |
|---|---|---|
| 0 | `amrex_slab` (BoxArray, z-ordered mapping, the E3 levels on slabs); the 128-bit, chain and plane-partial helpers in an internal header (today duplicated in `rt_level.cpp` and `pm_level.cpp`); the sparse per-rank sources (§6) | the `amrex` label unchanged; E3's `amrex_state` tests with `amrex_slab = 1` byte-identical |
| 1 | `MomentLevel`: M1 (prep, calc_rt_1, injection, transport, deposit, physics, step 3, redshift, budgets, `amrex_moment_stats`) in `rt_driver` and `radhydro` | §7 for M1 |
| 2 | OTVET: the FFT split into pieces (serial byte-identical), component ranks, thin direct and FFT, split, screen | §7 for thin, split, screen |
| 3 | the cone (pyramids, march, chains, blend with the z-transpose) | §7 for the cone |
| 3b | `recomb_rad = 1` with moments on the level (and case A, `recomb_rad_lines`) | §7 with `recomb_rad = 1` for M1 and the cone |
| 4 | checkpoint and restart through the slabs | restart tests |
| 5 | klone 256³ / 512³ validation and timing on 1/2/4 nodes | §7, §8 |
| 6 | 1024³ cone run and the comparison | §8 |

**Progress (2026-10-14):** phases 0-4 done, bitwise against the serial
moment RT (label `amrex_moments`). As built, beyond the plan:

- The sources needed no rank-0 `calc_halo` (§6).
- The serial per-cell arithmetic of the moments, the FFT convolution and
  the cone are shared code (`moment_cell.hpp`, `OtvetConvolver`,
  `OtvetCone`, `otvet_mip_coarsen`, `otvet_cone_prolong`, the blend
  helpers), the serial path calling the same functions (byte-identical,
  checked against the previous build on 10 serial moment cases).
- The cone gathers the finest mip level a march samples (`lmin`: 1 with
  the defaults; 0 when r_c tan θ < 1, as with `otvet_cone_rc = 1`), not
  always level 1.
- `amrex_slab_planes` must be a multiple of `otvet_cone_coarse` and of
  2^lmin for the cone (checked).
- Checkpoints: the version-7 format unchanged; a serial moment checkpoint
  restarts on the level byte-identically.
