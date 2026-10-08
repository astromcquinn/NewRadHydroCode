# Design: parallel ray tracing (Phase 2) and compact rays and adaptivity (Phase 3)

Status: proposal. Grounded in `reference/Radhydro/ray_tracing_tools.f90` (RTT),
`reference/Radhydro/rt_tools.f90` (RTL), `reference/Radhydro/domain_tools.f90`
(DT), `reference/Radhydro/radhydro_global.f90` (RG) and `radhydro.f90` (RH).
`RTT:296-306` means lines 296–306 of `ray_tracing_tools.f90`. Statements I am
not certain of are marked **[unsure]**, and the findings that matter for the
Phase 1 port are collected in §0.

---

## 0. Findings that affect Phase 1 (faithful port / harness)

I found these while reading the code. Both the Fortran harness and the C++
port need to reproduce them exactly. They are also the places where a
"natural" re-implementation will silently differ from the original.

1. **Dead ("zombie") rays enter the photon-count PDF.** `ray_remove`
   (RTT:1273-1340) only unlinks a ray and pushes its index onto the free list.
   It does not zero `nside`/`ngam`. `calc_ray_1` (RTL:1085-1116) scans the
   *whole* storage range `ray_cpu(1:2)` and counts every slot with
   `nside > 0`. As a result, rays deleted earlier in the same hydro step are
   still counted in `ngam_pdf`, and through it in `ngam_ray(:)`, during RT
   substeps 2..`nstep_rt`. That holds until a slot is reused by `ray_create`
   or zeroed by `ray_buffer_deallocate` in `post_rt` (RTT:2011-2031). The
   merge/cast/split thresholds therefore depend on deletion history and on
   free-list (LIFO) reuse order. A port that keeps a compact array of live
   rays will get different thresholds. Recommendation for Phase 1: emulate
   the storage pool and free lists exactly, with `Ncpu = 1`. Note also that
   `ngam_avg/Nray` (RTL:950) divides by the pool high-water mark `Nray`, not
   by the number of live rays. It is only printed, so it is harmless.
2. **`ngam_ray(j)` is not reset when the ray count drops.** In RTL:952-961,
   `ngam_ray(j)` is set only when the cumulative count exceeds
   `j*int(Nray_max/64.)`. Otherwise it keeps the value from the previous
   substep, or its initial 0 (RTT:88). Consequences:
   * While fewer than `Nray_max/4` slots are counted, `ngam_ray(16) = 0`.
     Every ray with any photons is then cast (RTT:300) and never merged, and
     with `ngam_ray(32) = 0` every ray with `dr < 0` is split (RTT:303).
     Merging is purely a **memory-pressure** mechanism.
   * The values are stale after the count falls again. They are only reset
     when `ngam_avg == 0` (RTL:947-948).
3. **Merge buffer overflow is unchecked.** `merge1_type` holds at most
   `ir(10000)` per cone (RTT:205-208), and RTT:333-335 writes without a bound
   check. The local array is 192 × 10000 × 8 B ≈ 15 MB of stack per thread
   (`ulimit -s` / `OMP_STACKSIZE` matter). Only cones 0..47 are used, because
   `nside_merge = 2` in both branches of RTT:259-266. The array is
   dimensioned 0:191 (the `nside_init = 4` pixel count).
4. **float32 positions and the `1E-6` nudge.** `ray%x` is `real(4)`
   (RG:138). `ray_cast_3` advances to `x0 + (l+1E-6)*dx` (RTT:1114) and
   stores the result as float32. The ulp of a float32 near x = 32 is 3.8e-6,
   near 512 it is 6.1e-5 and near 1024 it is 1.2e-4. The nudge is therefore
   usually lost, and rays end up exactly on cell faces. `ray_segment` handles
   rays on a face (RTT:1193-1200). The *home cell* used for the chain is still
   `int(x)` (RTT:1143-1145), so a ray sitting on face `x = m` and moving in
   −x is chained in cell `m+1` but deposits first into cell `m`. The
   behaviour is consistent but subtle. The port must store float32 positions
   and do the arithmetic in double, exactly as Fortran does: `real(4)`
   storage with `real(8)` locals.
5. **At most 3 segments per substep.** `ray_segment` (RTT:1164-1217) computes
   up to 3 segments with total length ≤ 1. A unit step along a near-body
   diagonal can cross 3 faces, which needs 4 segments. The remainder is
   dropped, and the ray then advances by less than 1 (`r = r0 + l`,
   RTT:1115).
6. **`dr` is redundant.** Every code path keeps `dr = rmax_from_nside(nside) - r`:
   source (RTT:440-441), split (RTT:541-542), merge (RTT:722-723) and advance
   (RTT:1115-1116). It is only redundant mathematically, though. In float32
   `dr` and `rmax-r` round differently, so the faithful port must keep `dr`.
   Compact encodings (§4) can drop it.
7. **`end_domain_rt` (DT:314-331) writes the shared `domain` array outside
   the critical section** that `set_domain_rt` uses. This is a benign-looking
   race that is irrelevant at `Ncpu = 1` **[unsure whether it can deadlock
   or double-schedule with many threads]**.
8. HEALPix `pix2vec_nest`/`vec2pix_nest` come from the external HEALPix
   library (`use pix_tools`, RTT:28). HEALPix is GPL-2. Anything we ship
   (and every GPU kernel) needs its own implementation of these two
   functions. Górski et al. 2005 describe the algorithm. The implementation
   must match the library bit-for-bit on the pixel index and to round-off on
   the vector.

---

## 1. The current algorithm: data flow, parallel structure, and why it is hard

### 1.1 Data structures

* `rt_grid(Nrt,Nrt,Nrt)` is an array of structures with 36 float32 per cell
  (144 B): `rhog fHI fHeI fHeII temp Luv sfr`, `dx(3) gamma(3) lambda(3)`
  and `nu ngam fgam tau` (each `Nfreq = 5`) (RG:117-122).
* `ray(Nray_max)` is an array of `ray_type` (RG:134-140): `nside, ipix`
  (int4), `next, prev` (int8), `r, dr`, `x(3)`, `dx(3)`, `nu(5)` and
  `ngam(5)` (real4). That is 8 + 16 + 8 + 24 + 40 = **96 B**. The pool is
  sized `Nray_max = 64 Nrt^3` (RG:71). On top of that come `ray_index`
  (8 B/ray) and per-thread free lists (`Nfree_max = Nray_max/(4 Ncpu)`
  int8 per thread, i.e. 2 B/ray in total).
* `htoc_ray(4,Nrt,Nrt,Nrt)` holds the head and tail of a doubly linked list
  of rays per cell. Slots 1:2 are a frozen copy taken at the start of the
  sweep. Slots 3:4 are the live list (RTT:269-272, RTL:997-1003).

### 1.2 Per hydro step

`radiative_transfer` (RTL:167-280) runs `prep_rt` (gas → `rt_grid`, buffer
allocation, RTL:286-493) and then `nstep_rt` substeps. After that, `post_rt`
runs (free buffers, `ray_redshifting` RTT:1524-1771, final tau and
diagnostics). `nstep_rt = nint(dl_c/dl_rt*cscale)` with `cscale = 1`
(RH:333-341), so a ray moves one RT cell per substep at the true speed of
light. `cscale < 1` would give a **reduced speed of light** for free, because
rays still move one cell per substep.

### 1.3 Per RT substep

| phase | code | what happens | access pattern |
|---|---|---|---|
| 1a | `calc_rt_1` RTL:1016-1082 | zero `nu, ngam, fgam, gamma, lambda`; recompute `tau(5)` where `fHI<0.99 or fHeI<0.99` | per cell, independent |
| 1b | `calc_ray_1` RTL:1085-1116, RTL:946-969 | global histogram (2083 bins of 0.0043 dex) of `Σ ngam` over all *slots with nside>0*; `ngam_ray(j)` = the value at which the bright tail holds `j·Nray_max/64` rays | global reduction then a serial scan |
| 1c | `ray_buffer_supplement` RTT:1840-1959 | rebalance per-thread free lists | bookkeeping |
| 1d | `ray_tracing_1` RTT:195-391, scheduled by the domain colouring of RTL:982-990 | for every cell with rays or a source: emit source rays; classify every ray in the frozen chain; cast, split, remove, or bin into one of 48 cones and merge per cone | per cell; **writes to neighbour cells** |
| 2 | `radiative_physics` RTL:1374-1916 | per cell: `nu = Σngam·nu/Σngam`, normalise `dx`, implicit ionisation and heating, photon consumption (RTL:1695-1705); `fgam = nrad_left/ngam` (RTL:1461-1467, 1809) | per cell, independent |
| 3 | `ray_tracing_3` RTT:772-851 → `ray_cast_3` RTT:1037-1158 | per ray: recompute the same ≤3 segments, multiply the photons deposited in each segment by that cell's `fgam`, advance by the length travelled (≤1), delete if `Σngam < max(ngam_min, fgam_min·Σngam_old)`, relink into a new cell | per ray; reads ≤3 neighbour cells; relinks |

Classification in step 1d (RTT:296-306), applied in order, with the last
match winning:

```
merge   default
cast    if Σngam > ngam_ray(16)            (bright: never merged)
cast    if r == 0                          (fresh source ray)
split   if dr < 0 and Σngam > ngam_ray(32) (passed rmax(nside) and bright enough)
remove  if Σngam < ngam_min (1e-6)
remove  if nside == nside_max (8192)       (r ≳ 4000 cells)
```

* **Source** (RTT:394-498). If the cell is optically thick
  (`exp(-min tau) < 1e-2`), the substep's photons are added directly to the
  cell. Otherwise the cell emits `nray_init = 12·4² = 192` rays at the
  `nside = 4` pixel centres, from the cell centre, with `r = 0`.
* **Cast** (`ray_cast_1`, RTT:970-1034). The ray does *not* change. Its
  photons are distributed over ≤3 path segments. Segment n receives the
  fraction `max(fdl, 1-exp(-tau·dl))` of what remains, where
  `fdl = dl_n/Σ_{m≥n} dl_m`, so the last segment receives everything left.
  Each segment adds `ngam, ngam·nu, dx·Σngam` (13 floats) to its cell. **All
  of the ray's photons are put on the grid.** Step 2 decides the absorbed
  fraction, and step 3 gives each ray back its share `fgam_cell · deposited`.
  Ray–cell coupling is therefore exactly photon-conserving.
* **Split** (`ray_cast_split`, RTT:501-621). nside doubles. r is clamped to
  `[rmin, rmin+1]` of the new level, and the effective origin is
  `x0 = x − r·dx`. Each of the 4 children gets direction = child pixel centre
  + parent's deviation from its own pixel centre (RTT:532-534, 569-572), and
  position `x0 + r·dx_child` (≤ ~0.5–1 cell from the parent). Photons are
  split ¼ each. The children are chained into whatever cell they land in, and
  each is cast once.
* **Merge** (`ray_cast_merge`, RTT:624-766). All merge-class rays in the same
  cell and the same `nside=2` cone (48 cones) are combined. The new ray gets
  photon-weighted position, direction, `<r²>^{1/2}` and `dr`, and the summed
  photons per group with energy-weighted `nu` (energy is conserved). Its
  nside is recomputed from r, and r is reset to `max(rmax−dr, rmin)`. The
  merged ray is cast once.

### 1.4 Current parallel structure (OpenMP, shared memory)

* Cell sweeps in steps 1d and 3 are parallel over **27-coloured domains**
  (DT:178-331). The RT grid is cut into `n1³` cubic domains, with `n1` even
  and `≥ 2·nint(Ncpu^{1/3})`. For Ncpu = 16 and Nrt = 512 that gives
  n1 = 8, so 512 domains of 64³. A thread may take a domain only when none of
  its 26 neighbours is in progress (`set_domain_rt`, a busy loop inside a
  critical section). This makes the neighbour writes (deposits, relinking
  into neighbour chains) race-free without atomics.
* Ray storage is a global pool with per-thread free lists. `ray_create` and
  `ray_remove` take a critical section only when the local list
  under- or overflows (RTT:1223-1340).
* Steps 1a, 1b and 2 are trivially parallel over k-slabs.
* Results are reproducible only up to floating-point summation order, both
  in the deposits and in the chain order that determines merge order. The
  integer ray stats can also differ between thread counts, because chain
  order changes which slots are reused and so the zombie content of §0.1
  **[unsure how often in practice]**.

### 1.5 What makes it hard to port to GPU and MPI

1. **Pointer-chasing per-cell linked lists.** The traversal is serial within
   a cell, and relinking across cells needs exclusive access to the
   neighbours. That is why the 27-colour scheme exists. Neither carries over
   to GPU, and on CPU the colouring limits parallelism to ≤ `n1³/27`
   concurrent domains with coarse granularity.
2. **Dynamic population.** Rays are created (192 per emitting source per
   substep; ×4 per split) and destroyed every substep. The free-list
   allocator is inherently serial-ish.
3. **Data-dependent grouping.** Merging needs all merge-class rays of a
   (cell, cone) together, which is a group-by.
4. **Extreme imbalance.** The rays-per-cell histogram (RTL:825-858 bins up to
   ">1024") spans 0 to thousands. Sources are clustered, and before overlap
   all rays live inside a few HII bubbles.
5. **Ray counts.** The design budget is 64 rays/cell. Splitting keeps about
   one ray per cell² of the sphere at radius r (`ray_per_cell = 1`,
   RTT:1420-1518), and a source emits a new generation every substep. A cell
   at distance r < c·t from a source therefore holds about 1–4 rays per
   source until merging starts. With N_vis sources in view, a cell holds
   about N_vis·(1–4) rays before merging. Merging caps the merge-class
   population at 48 per cell, but bright rays above `ngam_ray(16)` are
   unlimited.
6. **Finite c at 1 cell per substep.** The number of substeps per run is
   huge. From z = 20 to 6 the comoving light-travel distance is ≈ 2630 Mpc,
   or ≈ 370 crossings of a 5 Mpc/h box:

   | Nrt | substeps per hydro step (z=10, da=1e-4) | substeps z=20→6 |
   |---|---|---|
   | 32 | 12 | 1.2·10⁴ |
   | 512 | 196 | 1.9·10⁵ |
   | 1024 | 392 | 3.8·10⁵ |

   Per-substep overheads (global PDF, sort, MPI exchange, kernel launches)
   are paid 10⁵ times. They have to be O(Nray) with small constants and
   latency-hidden.

---

## 2. GPU design

### 2.1 Programming model

| option | pros | cons |
|---|---|---|
| **Kokkos (recommended)** | One source for CUDA, HIP, SYCL and OpenMP. `parallel_for/reduce/scan`, `TeamPolicy` with scratch, atomics on any type, `ScatterView` for CPU deposits, `Kokkos::sort`/`BinSort`, unmanaged views over MPI buffers. Mature on DOE machines. CPU builds via the OpenMP backend reuse the same kernels | Its sort is not always the fastest **[unsure of current `sort_by_key` performance; Kokkos 4.x has `Kokkos::Experimental::sort_by_key`]**. Build complexity |
| CUDA + HIP direct (CUB / rocPRIM) | Best-in-class radix sort, segmented reduce, `DeviceSelect`, warp intrinsics | Two code paths plus a CPU path. More code to validate |
| OpenMP target | Minimal syntax, Fortran-like | No portable sort, scan or segmented primitives. Compiler quality varies. Hard to control shared memory and warps |
| C++ `std::execution::par` (nvc++ stdpar, roc-stdpar) | Standard C++; sort, scan and reduce are available | Relies on unified memory. No team or scratch control and no ScatterView. AMD support is young **[unsure]**. Weak for the atomics-heavy deposit kernel |

**Recommendation: Kokkos**, with a thin `primitives/` layer (`sort_pairs`,
`segmented_reduce`, `exclusive_scan`, `compact`) that defaults to Kokkos and
can be switched to CUB or rocPRIM per backend if profiling says so. All
physics kernels (`ray_segment`, cast, split, merge, `radiative_physics`)
become `KOKKOS_INLINE_FUNCTION`s shared with the serial C++ port. Each piece
of physics then has a single implementation.

### 2.2 Replace chains with sorted ray arrays

The state between substeps is a structure of arrays of N live rays, **sorted
by home cell**, plus `cell_start[Ncell+1]`, which is compressed sparse row
(CSR) offsets. No links and no free lists are needed. Each substep:

```
 0. tau update (cell kernel)                                   [= calc_rt_1]
 1. thresholds: histogram of log10(Σngam) over live rays (int64 bins,
    device atomics or per-team histograms) -> cumulative scan -> ngam_ray(1:64)
    (exact integer arithmetic, deterministic)                  [= calc_ray_1]
 2. classify: per ray -> {cast, split, merge(cone), remove}  (pure function)
 3. merge grouping: key = (cell, cone) for merge-class rays. Rays are
    already sorted by cell, so only a *segmented sort by cone inside each
    cell* is needed (team-level sort in scratch: most cells hold <=64 rays;
    long cells go to a global fallback sort)
 4. counts:  remove 0, cast 1, split 4, merge 1 per (cell,cone) segment head;
    source cells: +192 (or 0 if optically thick, deposit directly)
 5. exclusive_scan(counts) -> output offsets; allocate new arrays (size known)
 6. build new rays (expand): copy cast rays, write 4 children per split,
    segmented reduction per merge segment (19 doubles: w, w*xd(3), w*r^2,
    w*dr, w*dx(3), ngam(5), ngam*nu(5)), emit source rays
 7. cast_1: per new ray, deposit 13 floats into <=3 cells   (atomics, §2.4)
 8. radiative_physics: cell kernel                          [= rt_step_2]
 9. cast_3: per ray, attenuate/advance; flag deletions; compute new home cell
10. re-sort by home cell: rays moved at most ~2 cells, so use a counting
    sort (histogram of new cells -> scan -> scatter) or a radix sort on a
    32/64-bit key; stream-compact the deleted rays in the same pass
```

Steps 4–6 make split and merge **stream expansion and compaction with prefix
sums**. Merging is a **segmented reduction** over contiguous (cell, cone)
runs. All of this is O(N) and allocation-free if the output buffers are
sized with headroom: the output is at most 4·N + 192·N_src.

Equivalence with the Fortran: within a cell, the order of operations only
changes floating-point summation order. Merging uses only the frozen
pre-sweep rays of that cell. Rays that relocate into a neighbour during
step 1d are not reprocessed in this substep. Fresh rays are cast exactly
once. The array design preserves all of this if the relocated rays are not
regrouped until step 10. §0.1 (zombie PDF) is the one *semantic* difference.
The Phase 2 code should **not** reproduce it. Instead the C++ serial port
needs a switch `--pdf-live-only` so that the two can be validated against
each other (§6).

### 2.3 Determinism

Validation (§6) needs *bitwise* reproducibility across thread counts,
backends and MPI decompositions. Recipe:

* **Canonical in-cell order.** Sort by `(cell, class, cone, payload-hash)`,
  where the tie-breaker is a 32-bit hash of the ray's bit pattern. Rays with
  identical payloads contribute identically, so their relative order does not
  matter. The hash makes the order independent of where the ray came from.
* **Merge** reductions run sequentially in canonical order inside each
  segment. Segments are short: at most 48 per cell, typically 1–10.
* **Deposits** use **fixed-point int64 atomics**: integer addition is
  associative. Use one scale per quantity, for example `ngam` in units of
  2⁻⁴⁴ with a range to ~10⁶. Values below `ngam_min = 1e-6` are discarded
  by `radiative_physics` anyway (RTL:1409-1414), so the absolute resolution
  of ≈ 6·10⁻¹⁴ is far below anything that matters. Signed `dx·Σngam` works
  in two's complement. Cost: 8 B accumulators (13 × 8 = 104 B/cell of
  scratch). This is on by default in a "reproducible" build and can be
  switched to float atomics for speed.
* **Global sums** (diagnostics, `ngam_rays_total`) use the same fixed point,
  or a binned reproducible sum. The integer histogram for `ngam_ray` is
  exact.
* **Math library.** `exp`, `log`, `sqrt` and `atan2` differ in the last ulp
  between glibc, CUDA and HIP. A reproducible build ships its own
  `rr_exp/rr_log` (correctly rounded or fixed-polynomial) and compiles with
  `-ffp-contract=off` (nvcc `--fmad=false`). Expect a ~10–30% slowdown
  **[unsure]**. The fast build uses the vendor functions.

### 2.4 Deposition: atomics versus sort-based scatter

Each cast deposits into ≤3 cells, with 13 values each. Rays are sorted by
home cell, so adjacent threads target the *same* few cells, which is
worst-case contention. Options:

1. **Warp- or team-aggregated atomics (recommended).** Inside a team that
   owns a run of rays from the same cell, reduce in scratch memory by target
   (home cell plus ≤26 neighbours, a 27-slot table × 13 values), then issue
   one atomic per (team, target, value). This cuts the global atomic count by
   about the mean rays per cell.
2. **Sort-based scatter.** Emit `(target_cell, 13 values)` records (≤3N),
   radix-sort by target, then segmented-reduce. This is deterministic even in
   float, but it costs ~3× the ray traffic. Keep it as the fallback and for
   checking.
3. **Gather (pull)** from each cell over rays in a 5³ neighbourhood. This is
   wasteful with dense rays. Not recommended.

### 2.5 Layout and precision

* SoA for rays: `x[3][N]`, `dx[3][N]`, `ngam[5][N]`, `nu[5][N]`, `r[N]`,
  `nside_ipix[N]`. Kernels read all fields, so SoA mainly helps coalescing
  and selective movement during the sort. Sort **an index permutation with
  64-bit keys**, then gather the payload once.
* SoA for the grid: `tau[5][Ncell]`, `fgam[5][Ncell]`, and so on.
  `radiative_physics` needs ~25 values per cell.
* Precision: keep **storage in float32** as in the Fortran and do
  **arithmetic in double** in kernels. Consumer GPUs have 1/32–1/64 FP64
  throughput. H100, MI250X and MI300 have full-rate FP64, and the kernels are
  memory-bound anyway. Positions are stored as cell index plus a 16-bit
  offset (§4), not as global float32, which also removes §0.4.

### 2.6 Memory estimate

Rays: payload P bytes plus sort scratch (8 B key, 4 B index) plus a second
payload buffer for the expand/gather, so about **2P + 12 B per ray**.

| | P | per ray |
|---|---|---|
| Fortran `ray_type` without links (float, SoA) | 80 B | 172 B |
| Compact R32 (§4) | 32 B | 76 B |

Grid: 36 float32 + 9 gas float32 + 13 int64 accumulators + 4 B CSR offset ≈
**288 B/cell**.

| Nrt | cells | grid | rays at 16/cell (= Nray_max/4, where merging starts) | at 64/cell (= Nray_max) |
|---|---|---|---|---|
| 512 | 1.34·10⁸ | 39 GB | 2.1·10⁹ rays: **370 GB** (float) / **163 GB** (R32) | 8.6·10⁹: 1.48 TB / 0.65 TB |
| 1024 | 1.07·10⁹ | 309 GB | 1.7·10¹⁰: **2.95 TB** / **1.30 TB** | 6.9·10¹⁰: 11.8 TB / 5.2 TB |

For comparison, the Fortran pool at Nrt = 512 is `8.6e9 × (96+10) B ≈
0.9 TB` plus 4.3 GB of `htoc_ray`. With 80 GB GPUs at about 60 GB usable,
Nrt = 512 needs **4 GPUs (R32) or 7 GPUs (float)** at 16 rays/cell, and
Nrt = 1024 needs **~27 or ~55 GPUs**. Rounded up: **8 GPUs for 512³ and
64 GPUs for 1024³**. Phase 3 should replace the 64·Nrt³ budget with an
error-driven criterion (§5). The average ray density is the main memory
knob.

### 2.7 Throughput estimate

Per ray per substep, with the R32 payload: read and write in classify/expand
(~64 B), sort (key and index read+write over 2–4 radix passes ≈ 50–100 B, or
~30 B for the counting sort), payload gather (64 B), cast_3 read+write
(64 B), and ≤3 aggregated deposits. That is about **250–300 B of HBM
traffic** per ray-substep. Compute per ray is about 3 segment searches, 15
`exp` and 1–2 `pix2vec`. That is ~300–500 flop, so the kernels are
bandwidth-bound.

At 3.35 TB/s (H100) the ideal is ~10¹⁰ ray-substeps/s. **Planning figure:
2–3·10⁹ ray-substeps/s per GPU** (25–30% of peak, accounting for atomics
and launch overheads), and roughly 4–5× lower with the float 80 B payload
**[estimate; to be measured in stage 2c]**.

* Nrt = 512: 2·10⁹ rays on 8 GPUs gives ~0.1 s per substep. Over 1.9·10⁵
  substeps that is about **5 h wall and 40 GPU-h**, an upper bound because
  the early phase has far fewer rays. With `cscale = 0.2`, ~1 h.
* Nrt = 1024: 1.7·10¹⁰ rays on 64 GPUs gives ~0.1 s per substep. Over
  3.8·10⁵ substeps that is **~11 h wall and 700 GPU-h**, before MPI
  overheads.

The per-substep fixed cost (≈10 kernel launches, 1 small allreduce, halo and
migration messages) must stay under ~5 ms to remain below 10%. Use CUDA
graphs or HIP graphs for the substep sequence, and overlap halo exchange
with interior work.

---

## 3. MPI design

### 3.1 Decomposition

* The RT grid gets its **own decomposition**, independent of the hydro or
  gravity one. Gas state (9 floats per cell) is remapped in and out once per
  hydro step. That costs ≤1/(nstep_rt ≈ 200) of the RT communication volume.
* **Unit: bricks** of 16³ (or 32³) cells, ordered along a **Hilbert
  space-filling curve**. Each rank owns a contiguous SFC range chosen so
  that the **weighted cost** is equal. The cost model is
  `w_brick = a·N_cells + b·N_rays + c·N_src_cells`, with `a, b, c` fitted
  from per-kernel timers. A first guess is b/a ≈ 3, since physics is ~3×
  the cost of a ray. Weights are measured each substep and **repartitioned
  every hydro step**, or when imbalance exceeds 20%. Bricks are migrated
  with their rays.
* Slabs are simplest but scale only to Nrt ranks and have a large halo per
  volume. Uniform blocks are good for uniform load but bad for clustered
  sources. Weighted SFC bricks handle clustering, and with Hilbert ordering
  the surface/volume ratio is within ~1.5× of blocks **[typical,
  unverified for this workload]**.
* Within a rank, the GPU works on its whole subdomain as one array. CPU-only
  runs use Kokkos OpenMP, so the colouring of DT is no longer needed.

### 3.2 Ghost zones

Reach of one substep from a ray's home cell at the start of step 1:

* split children relocate ≤1 cell (displacement ≈ `r_new·|dx_child−dx|`
  ≈ 0.5 + 1/nside cells; RTT:541-576 **[unsure about the worst case near
  HEALPix poles; add an assertion]**);
* cast_1 from the new position reaches ≤1 more cell (3 segments, total length ≤1);
* cast_3 then reads the same ≤3 cells.

So deposits and reads reach **Chebyshev distance 2**. Two ways to handle it:

* (A) a **2-cell halo**. For a 128³ subdomain that is ~9.6% extra volume.
* (B) a 1-cell halo plus a migration step right after step 1d.

**Recommend (A)**: one migration per substep instead of two.

Per substep exchanges:

| when | what | direction | size per ghost cell |
|---|---|---|---|
| after 1a | `tau(5)` | owner → ghosts | 20 B |
| after 7 (cast_1) | `ngam, ngam·nu, dx·Σngam` accumulators (13 × int64) | ghosts → owner, summed | 104 B |
| after 8 (physics) | `fgam(5)` | owner → ghosts | 20 B |
| after 9 (cast_3) | rays whose new home is off-rank | to ≤26 neighbour ranks | 32–80 B per ray |
| step 1b | 2083-bin int64 histogram | allreduce | 17 KB |

For a 128³ subdomain with a 2-cell halo there are 132³−128³ ≈ 2.0·10⁵ ghost
cells. At 144 B each over the three grid exchanges, the grid halo is ~30 MB
per substep. Ray migration: only rays in the outer cell layer that move
outward leave, which is ~1–5% of local rays per substep for isotropic rays
(¼–1 × 6/128). Use GPU-aware MPI with persistent neighbour
communicators (`MPI_Neighbor_alltoallv` or persistent sends) and overlap the
interior compute.

### 3.3 Merging across rank boundaries

Merging groups by (home cell, cone), and migration completes before step 1.
All rays of a cell are therefore always on the owner rank, so **no
cross-rank merge is ever needed**. Split children and merged rays that move
into a ghost cell during step 1d stay on the source rank for this substep,
which matches the Fortran rule that relocated rays are not reprocessed.
They migrate with everything else after step 9.

### 3.4 Decomposition-independent results

These rules make the result bitwise identical for any rank count and
partition. They are required for validation:

1. **Thresholds** come from an integer histogram, which is an exact
   allreduce.
2. Canonical in-cell order (§2.3), which depends only on ray content.
3. **Deposits** use fixed-point int64, including the ghost→owner sums,
   since integer addition is associative.
4. **Source emission** is per cell and deterministic. Any randomisation
   (currently disabled, RTT:446) would use a counter-based RNG keyed by
   (cell, substep), never by rank or thread.
5. **Ray IDs** are not needed. If one is wanted for debugging, derive it
   from (cell, substep, slot in canonical order).
6. **Diagnostics** use fixed-point or reproducible sums.

Test (stage 2d): 1, 8 (2³ blocks), 27 and a random weighted SFC partition
must give `compare.py` with `--rtol 0`, i.e. bitwise-identical outputs and
identical ray stats.

---

## 4. Ray encoding

### 4.1 Information actually used by the algorithm

* `ipix` is used only in split (RTT:521-572), as the parent pixel centre and
  the children `4·ipix+n`. `nside` is used in split, removal (`nside_max`)
  and `rmin/rmax`.
* `dx` is used everywhere. Two facts matter here:
  * Source rays start exactly at pixel centres (RTT:468), and children get
    `centre(child) + (dx_parent − centre(parent))`. For a ray that was never
    merged the deviation is therefore **exactly zero** up to round-off: its
    direction is the HEALPix centre `(nside, ipix)`.
  * Only merged rays and their descendants have a non-zero deviation, and
    that deviation is bounded by about one pixel radius, because `ipix` is
    recomputed with `vec2pix` (RTT:572, 727).
* `x`, `r`: `dr` is redundant (§0.6). For never-merged rays,
  `x = x_src + r·dx` holds up to the `1E-6` nudges.
* `ngam(5), nu(5)`: in the default (non-`FREQ1`) path, `nu` changes only by
  energy-weighted merges and by redshifting, which clamps it to its band
  (RTT:1749).

### 4.2 Proposed compact encoding "R32" (storage only; decode to float/double in registers)

| field | bits | encoding | error |
|---|---|---|---|
| home cell | 0 stored (CSR position); 32 in MPI messages | — | — |
| position in cell | 3 × 16 | fixed-point offset, 1/65536 cell | ≤ 7.6·10⁻⁶ cell (better than float32 at x ≳ 64) |
| level `lside` | 4 | 2..13 | exact |
| `ipix` (nested) | 30 | `12·8192² < 2³⁰` | exact |
| direction deviation | 2 × 8 (or 0 with flag `on_centre`) | offsets along the local (θ, φ) tangent basis in units of the pixel radius / 128 | ≤ 0.4% of a pixel ⇒ lateral error ≈ 0.004·(pixel spacing at r) ≲ 0.004 cell, because spacing ≈ 1 cell at the split radius |
| r | 16 | fixed-point, 1/8 cell (r < 8192) | ≤ 1/16 cell, which shifts the split radius by ≤ 1/16 cell; nside follows `lside` exactly |
| `Σngam` | 32 | float32 | exact (as now) |
| spectral shape | 4 × 12 | `log2(ngam_i/Σngam)` quantised to 1/128 over [−32, 0] | rel ≤ 2.7·10⁻³ per group; photon number exact via the total |
| `nu_i` | 5 × 8 | linear position inside band `[freq(i), freq(i+1)]` | ≤ band/512 ≈ 0.01–0.05 eV |
| flags | 2 | `on_centre`, spare | |

Total: 48 + 4 + 30 + 16 + 16 + 32 + 48 + 40 + 2 = **236 bits, packed into
32 B**. The original is 96 B (80 B without links), so this is a 3× saving.
A minimal **R24** variant drops `nu` (use band means, or store one global
spectral index) and uses 8-bit shapes. Error grows to ~2% per group, so
R24 should be tested, not assumed.

Error control:

* Re-quantising every substep accumulates error. Shape errors of ~10⁻³ per
  step random-walk to ~10⁻³·√(10³) ≈ 3% after 1000 substeps. Use
  **round-to-nearest in log space with deterministic dithering** (a hash of
  (ray, substep)) so the errors are unbiased, and only re-encode fields that
  changed. The total is kept exactly, so photon conservation is unaffected.
* Experiment (stage 3a): run uniform_8src and lognormal_16src with float
  versus R32 versus R24 storage. Report the field differences from
  `compare.py`, the global photon budget, and the I-front radius.
  Acceptance: R32 has max rel ≤ 10⁻³ on `fHI` where `fHI > 10⁻³`, ≤ 1% on T,
  and the volume-weighted xHII history within 10⁻³.

### 4.3 Origin-based encodings and merging

For a never-merged ray, (source id, `lside`, `ipix`, r) determines both
position and direction: `x = x_src + r·centre(ipix)`. That is about 30 + 4 +
30 + 16 bits, with no direction or position fields. It is attractive in
principle. However:

* **It does not survive merging.** A merged ray mixes origins. Its
  photon-weighted mean position is generally *not* on any ray from a single
  point, and its r is set to `max(rmax−dr, rmin)` (RTT:722), not to a
  geometric distance.
* **Split already breaks it.** Split builds an *effective* origin
  `x0 = x − r·dx` from the clamped r (RTT:541-547). After a merge-then-split,
  x0 is an arbitrary point. Storing x0 costs as much as storing x.
* The saving over R32 is ~46 bits (position, deviation) minus ~30 bits
  (source id), or ~2 B per ray, at the cost of two code paths. **Not
  recommended.** The `on_centre` flag (deviation = 0) gets the cheap part of
  the gain.

### 4.4 Alternative: per-cell angular maps (moments)

After a full merge, the merge-class rays of a cell are at most 48 "macro-rays",
one per `nside = 2` cone. Each carries its own sub-cell position, mean
direction, r and spectrum. That is already a **Lagrangian discrete-ordinates
scheme with 48 adaptive bins**. Storing it instead as a fixed per-cell
HEALPix map:

* Dense map at nside 2 with 5 groups only: 48 × 5 × 4 B = 960 B/cell.
  Adding the per-cone mean direction, offset and r, needed to avoid
  ray-effect artefacts and numerical diffusion, gives ~17 floats per cone,
  or 3.3 KB/cell. Sparse rays at 16/cell × 32 B = 512 B/cell. **The map only
  wins when most cones are occupied in most cells**, which happens in the
  fully ionised, post-overlap IGM.
* Transporting a map is upwind discrete ordinates: photons in cone c go to
  the neighbours along the cone axis with geometric weights. This has the
  known DOM ray effects (fixed directions) and angular diffusion.
  Near sources and I-fronts the ray representation is much better: sharp
  shadows, exact 1/r² dilution from the HEALPix geometry.
* A **hybrid** is the natural Phase 3+ direction. Use rays near sources and
  in cells with an I-front indicator (§5). Use per-cell maps, or further a
  per-cell M1-like moment (energy density + flux, ~9 floats per group) in
  the transparent ionised bulk where Γ depends almost only on the local
  photon density. Convert between them by merging rays into the map when a
  cell becomes "transparent", and by emitting 48 rays from a map cell when
  an absorber appears. Larger research item; only start it after §5
  experiments show that ray counts in ionised regions dominate the cost.

---

## 5. Adaptivity

### 5.1 What is fixed today

* Merging and splitting are driven by **quantiles of the ray-brightness
  distribution relative to the memory budget**. `ngam_ray(16)` is the cast
  threshold and `ngam_ray(32)` the split threshold (RTL:946-969, §0.2).
  Accuracy therefore depends on `Nray_max`, not on any error estimate.
* Cone size for merging is fixed at `nside_merge = 2` (48 cones) everywhere
  (RTT:259-266).
* Angular sampling is one ray per cell² at radius r (`ray_per_cell = 1`).
* Removal uses an absolute floor `ngam_min = 1e-6`, a per-substep relative
  floor `fgam_min = 1e-4` (RTT:1074, 1124), and an age limit `nside_max`
  (r ≈ 4000–8000 cells).
* Time: one cell per substep at c (`cscale = 1`, RH:334). `Nfreq = 5` fixed
  groups.

### 5.2 Proposals

1. **Error-driven merge.** Merge a cone's rays only if the merged ray would
   deposit nearly where the individual rays would. Let σ_θ be the
   photon-weighted angular spread of the rays in the (cell, cone), and
   `L = min(L_abs, L_max)` the distance over which their photons matter.
   `L_abs` is the distance to τ_HI = 1 along the cone, estimated from the
   local `tau(1)` (`L_abs ≈ 1/tau_cell`, capped). Merge if `σ_θ·L < ε_x`,
   with `ε_x` around 0.25–1 cell. In transparent ionised gas L is large but
   the deposited energy is small, so weight by the absorbed fraction
   `1−exp(−τ L)`.
2. **Adaptive cone resolution.** Use `nside_merge = 1` (12 cones) in
   transparent ionised cells (`fHI < 10⁻³` in the cell and its 26
   neighbours). Use 2 by default, and 4 (192 cones) or no merging in cells
   flagged as I-front (`10⁻³ < fHI < 0.9`, or `max_nbr τ_HI > 1` and
   `min_nbr τ_HI < 0.1`). The Fortran already has the `nside_merge` hook
   (RTT:259-266, `merge1` sized for 192 cones).
3. **Error-driven split.** Split if the ray's photons can change a cell
   noticeably: `Σngam·(1−e^{−τ_ahead}) > ε_s · nH_cell·vol / nH_avg·vol`,
   i.e. it can ionise more than a fraction ε_s of a cell. Otherwise let it
   continue coarse, as merge-class rays already do when `dr < 0`. The
   quantile thresholds can stay as a **memory safety valve**: they apply
   only if the ray count exceeds a hard budget.
4. **Reduced speed of light / substepping.** Expose `cscale`. Add
   **multi-cell casting**: m cells per substep, with ≤3m segments and
   `dt_rt = m·dx/c`. Also add an **infinite-c mode** for the time-independent
   limit: iterate substeps without advancing time until rays leave the
   domain or are absorbed, with physics applied once per hydro step
   **[needs care: the current physics consumes photons per substep;
   infinite-c requires accumulating Γ over the sweep]**. RSLA (reduced speed
   of light approximation) is known to bias the post-overlap photon density
   and Γ when the mean free path is long **[literature: Gnedin 2016 and
   follow-ups; check before relying on it]**.
5. **Adaptive frequency groups.** Compare Nfreq = 3, 5 and 8, and check the
   physically motivated shape encoding
   `ngam_i ∝ s_i·exp(−Σ_s σ_{s,i} N_s)` with 3 column densities. That is
   exact for a single unmerged ray in the absence of the
   photon-limiting clamp (RTL:1669-1676), and only approximate after merges.
6. **Photon budget diagnostics.** Track per substep: emitted, absorbed
   (`nion`), held in rays, lost to removal (`ngam_min`, `fgam_min`,
   `nside_max`), deposited directly at thick sources, and lost to the merge
   clamp `ngam < ngam_min → 0` (RTT:708-713). These are needed to judge any
   adaptive scheme.

### 5.3 Experiments and metrics

Test set: the RT comparison project problems of Iliev et al. 2006
(MNRAS 371, 1057), Test 1 (isothermal Strömgren sphere, pure H), Test 2
(with heating) and Test 4 (multiple sources in a cosmological density
field), adapted to our H+He chemistry. Add our `uniform_8src` and
`lognormal_16src`.

Reference solution: `ray_per_cell = 4`, merging disabled, `fgam_min = 10⁻⁶`.

| experiment | knob | metrics | pass criterion |
|---|---|---|---|
| E1 merge criterion | ε_x ∈ {0.25, 0.5, 1, 2} vs quantile | L1(xHII), L1(log T), I-front radius r_I(t) vs analytic (Test 1), cost = Σ ray-substeps, peak rays | ≥2× fewer ray-substeps at ≤1% r_I error |
| E2 adaptive cones | nside_merge 1/2/4 by indicator | as E1, plus shadow sharpness: an opaque clump test, ionised fraction behind the clump | shadow edge width ≤1 cell |
| E3 split criterion | ε_s ∈ {10⁻³, 10⁻², 10⁻¹} | as E1 | — |
| E4 RSLA | cscale ∈ {1, 0.3, 0.1} | xHII(z), Γ PDF in ionised cells, time of overlap | xHII(z) within 1% before overlap; report Γ bias after |
| E5 multi-cell cast | m ∈ {1, 2, 4} | as E1 | — |
| E6 groups/encoding | Nfreq 3/5/8; R32/R24 | T in ionised regions, xHeIII | ΔT ≤ 2% |

Each experiment produces an **accuracy-versus-cost Pareto plot** from
`compare.py --json` output against the reference run.

---

## 6. Roadmap with milestones and validation

| stage | deliverable | validation |
|---|---|---|
| **1** (current) | Faithful serial C++ port, Fortran harness, `tools/compare` | `compare.py` default tolerances (rtol 1e-5, ray stats exact) on the 3 problems |
| **2a** | Serial C++ with SoA ray arrays, sort-by-cell, expand/compact (no links). Switch `--pdf-live-only` in *both* the faithful port and 2a | 2a vs faithful port, both with `--pdf-live-only`: fields rtol 1e-5. Ray stats identical, or differences explained by threshold ties. Then faithful-zombie vs live-only to quantify the effect of §0.1 |
| **2b** | Reproducible mode: canonical order, int64 deposits, deterministic sums, portable math | bitwise identical (`--rtol 0`) for 1/2/8/32 OpenMP threads |
| **2c** | Kokkos port (OpenMP, then CUDA/HIP); primitives layer; benchmarks | GPU and CPU bitwise identical in reproducible mode; fast mode at rtol 1e-5 against 2b. Measure ray-substeps/s against the §2.7 estimate |
| **2d** | MPI: blocks, then weighted Hilbert bricks; 2-cell halos; ray migration | bitwise identical across 1, 8 and 27 ranks and random partitions; weak/strong scaling on uniform_8src and lognormal_16src scaled to 256³ |
| **2e** | Load balancing (cost model, repartition every hydro step) | imbalance ≤ 1.2 on clustered-source tests; time per substep vs rank count |
| **3a** | R32 encoding (storage only) + R24 variant | §4.2 criteria; memory per ray measured |
| **3b** | Adaptivity knobs (§5.2), each behind a flag, default off (= Phase 1 behaviour) | experiments E1–E6 with Pareto plots; pick defaults |
| **3c** | (optional) Hybrid map/moment representation for transparent regions | Iliev tests + cost |
| **4** | Coupling to hydro and gravity | — |

Every stage keeps the earlier behaviour available behind flags, so that
`compare.py` against the previous stage stays meaningful. Every PR runs
`tools/compare/test_compare.py` plus a 16³, 2-step smoke comparison.

---

## 7. Open questions and uncertainties

* How large §0.1 (zombie rays in the PDF) is in practice. It only matters
  once the ray count reaches `Nray_max/4`, which may never happen at
  N = 32 with `Nray_max = 64 N³`. In that case merging never triggers in the
  Phase 1 tests and the merge path goes **untested**. Suggest adding a test
  problem with a reduced `Nray_max` (e.g. `4·N³`) so that merge and the
  quantile logic are exercised.
* The worst-case relocation distance of split children (§3.2), to be
  checked with an assertion in 2a.
* The kernel throughput figures in §2.7 are estimates from memory traffic,
  not measurements.
* Which Kokkos version and sort API to use. To be decided at stage 2c start.
