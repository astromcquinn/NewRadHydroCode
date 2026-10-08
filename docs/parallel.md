# OpenMP parallelisation (Phase 1d)

The C++ code runs its hot loops with OpenMP (CMake option `RADRAY_OPENMP`,
default ON; the thread count is `OMP_NUM_THREADS` or the parameter
`threads`). The requirement was: results bit-identical for any thread count,
the serial path bit-identical to the Fortran harnesses, and the original
semantics preserved wherever a parallel scheme would change them.

* **Everything except the ray tracing proper** (`ray_tracing_1`/`3`) is
  parallel in every mode and gives the same bits as the serial code for any
  thread count. All Fortran-vs-C++ references pass unchanged with 1, 2, 3 and
  4 threads.
* **The ray tracing** runs in parallel only with the new flag
  `deterministic_parallel = 1` (which requires `legacy_rt = 0`). It uses its
  own, fixed sweep order (8³-cell domains, colour by colour), so its results
  are bit-identical for any thread count but differ from the legacy order at
  the round-off level (section 3). With `deterministic_parallel = 0` (the
  default) the ray tracing runs on one thread in exactly the original order.

## 1. Rules

1. A parallel loop writes only what its iteration owns (a pencil, a plane,
   a cell, a domain of cells) or per-thread buffers.
2. Floating-point sums are never OpenMP reductions. The addends are recorded
   per iteration and added afterwards **serially, in the serial loop order,
   with the original expressions**, so every accumulator sees the same
   addends in the same order as the serial loop. Integer sums and max/min
   are combined from per-thread partials (order independent).
3. Small grids (< 32³ cells) skip the parallel regions: the RT does
   thousands of substeps, and the fork/join cost would dominate (results
   are the same either way, `par::worth`).

## 2. What runs in parallel, and how

| component | scheme | identical to serial? |
|---|---|---|
| hydro sweeps (`hydro_step`) | pencils (j,k) in parallel; per-thread sweep workspace and FFT descriptor (FFTW planner called serially); counters and `vmax` combined as integer sums / maxima | yes |
| transposes | planes of the untouched index | yes |
| `summarize_gas` | per-cell values (the `pow`/`exp`/`sqrt` work and the `u_gas`/`w_gas` clamps) in parallel per batch of k-planes; the 60 sums then serially in k,j,i order with the original expressions; `cmax` a maximum | yes (hydro.log identical to the Fortran with 3 threads) |
| PM `calc_drift`, `link_dm` rows, density init, k-space Green's function, `accelerate_gas` | independent particles / rows / cells | yes |
| PM CIC deposit (`cic_mass_assignment_all`) | thread t owns the target planes [k_a,k_b]; it walks the source cells in exactly the serial order (domains 1..8, k,j,i, chains), restricted to the planes that can reach its slab (±2), and adds only into its own planes. Every density cell receives the same addends in the same order as the serial loop (no atomics, no reduction buffers). | yes |
| PM kicks (`accelerate_dm`) | particles independent; `arms`/`vrms` addends recorded per plane and summed in plane order; `amax`/`vmax` maxima | yes |
| FoF (`find_halos`) | the existing OpenMP union-find (CAS on the smallest index), independent of thread order | yes |
| RT `prep_rt`, `calc_rt_1`, `calc_rt_f`, `htoc` sync, `ray_redshifting` | per cell / per ray slot | yes |
| `calc_ray_f` (ray.log PDFs) | per-thread integer histograms | yes |
| `radiative_physics` (`rt_step_2`) | k-planes in parallel in batches of 4×threads; each plane records its `rt_stat` addends (`RPSink`: accumulator id, addend) in order; the batch is replayed in plane order; `ngfill`/`nsfill` integer sums | yes |
| `calc_ray_1` (threshold histogram) | legacy: serial (the `ngam_avg` sum is in storage order). `deterministic_parallel`: per-thread integer histograms, `ngam_avg` (ray.log only) in exact 2⁻⁴⁰ fixed point | yes / thread-independent |
| `ray_tracing_1`, `ray_tracing_3` | legacy: one thread, the original 8 domains in order 1..8. `deterministic_parallel`: section 3 | yes / thread-independent |

## 3. Deterministic parallel ray tracing (`deterministic_parallel = 1`)

The original parallelises the two ray sweeps over cubic RT domains with a
27-colour lock (`domain_tools.f90`): a thread may take a domain when none of
its 26 neighbours is in progress. The order in which domains are processed
then depends on timing, and with it the order of the deposits into shared
cells, the order of the chains (which decides merge order) and the ray
storage slots: the original is not reproducible with more than one thread.

### The scheme

* **Domains.** Cubes of d³ cells, d = 8 (the smallest d ≥ 8 with Nrt/d even;
  N = 64: 512 domains, N = 128: 4096, N = 192: 13824). The domain size does
  **not** depend on the thread count.
* **Colours.** Domain (i,j,k) has colour (i mod 2) + 2 (j mod 2) + 4 (k mod 2).
  Two domains of one colour are at least d = 8 cells apart in some
  dimension.
* **Reach.** Everything a cell's ray work touches lies within 3 cells of it:
  split children land ≤ ~1 cell away and are cast (≤ 1 more cell, ≤ 4
  segments), merged rays sit inside the cell (photon-weighted mean of
  positions inside it) and are cast, `ray_cast_3` advances ≤ 1 cell. Every
  chain insertion is checked at run time against the limit d/2 − 1 = 3
  (`det_check_reach`, an exception otherwise; it has never fired). The
  footprints of two domains of one colour (domain ± 4 cells) are therefore
  disjoint: no shared cell, chain, ray or `rt_grid` element.
* **Order.** For colour 0..7: all domains of the colour in parallel (dynamic
  schedule), each swept k,j,i as in the original. Because same-colour
  domains are independent, the result equals that of the serial
  colour-major sweep, whatever the thread count and schedule.
* **Ray storage.** Per-thread free lists refilled from / spilled to a shared
  pool under a mutex (`RayState::det_take`, `det_give`). Which slot a ray
  gets depends on timing, but nothing depends on slots with
  `legacy_rt = 0`: chains are linked lists ordered by insertion, freed slots
  are zeroed (B1), the threshold histogram is integer, and the two
  storage-order sums are replaced by canonical ones (`ngam_rays_total`
  summed over the chains in cell order; `ngam_avg` in fixed point). The
  `Ray/` restart files are written renumbered in canonical chain order.
* **Counters** (`ncst1`, `nspl1`, …) per thread, folded as integer sums and
  maxima; `Nray` is the live count.

### Proof of thread-count independence

1. Within a colour, the domains' footprints are disjoint (reach check), so
   the parallel execution of a colour is equivalent to any serial order of
   its domains; all state a domain reads that another domain of the colour
   could write is outside its footprint.
2. Colours are separated by barriers and executed in a fixed order.
3. The only shared mutable objects are the ray pool (slot identity only,
   which no result depends on — above) and the per-thread counters (integer
   sums/maxima).
4. All other RT phases are bit-identical to serial for any thread count
   (section 2).

Hence every output is a function of the input alone. Tested by the CTests
`thread_rt_det_lognormal_16src` and `thread_rt_det_uniform_8src_smallbuf`
(N = 32, 64 domains, 8 per colour: every field file and every stat of every
step byte-identical for 1, 2 and 4 threads), `sanitize_thread_rt_det` (the
same under ASan/UBSan), `restart_smoke_det_movie`, and the production runs
(N = 64 and 128: identical outputs with 1 and 4 threads restarted from the
same checkpoint, `tools/production/speedup.sh`).

### Difference from the legacy order

`tools/production/compare_order.py` runs the RT problems with
`legacy_rt = legacy_physics = 0` in the legacy order and with
`deterministic_parallel = 1` (N = 32, 10 steps). Relative differences of the
global quantities, and per-cell differences:

| problem | ⟨x_HII⟩_V | ⟨x_HII⟩_M | ⟨T⟩_M | nray_active | ngam in rays | max \|Δx_HII\| per cell | rms Δx_HII | max \|ΔT/T\| (ionised cells) |
|---|---|---|---|---|---|---|---|---|
| uniform_8src | ≤ 3.3e-9 | ≤ 3.3e-9 | ≤ 5.6e-9 | 0 | ≤ 6.1e-9 | 1.8e-7 | 3.0e-9 | 2.8e-7 |
| lognormal_16src | ≤ 9.0e-10 | ≤ 1.4e-9 | ≤ 2.9e-9 | 0 | ≤ 3.8e-9 | 1.8e-7 | 4.8e-9 | 2.8e-7 |
| lognormal_16src_smallbuf (merges) | ≤ 7.7e-10 | ≤ 1.4e-9 | ≤ 2.9e-9 | 0 | ≤ 1.4e-9 | 2.4e-7 | 5.4e-9 | 2.7e-7 |

The differences are float32 round-off of the deposit order (one ulp of a
cell's x_HII is 6e-8); the ray counts, splits and merges are identical in
these problems. (In general a threshold decision can flip at a round-off
tie and change the ray population locally; the statistics are unaffected.)

## 4. Ray memory

* `Ray` is 88 B (was 96 B): the chain links are int32 instead of int64
  (`Config::make` refuses `Nray_max ≥ 2³¹`, i.e. more than 2.1·10⁹ rays,
  which needs ~190 GB anyway). `ray_index` and the free buffer are int32 as
  well. No result changes (all references pass).
* The pool is one `calloc` of `Nray_max + 1` rays: zero pages are mapped
  lazily, so the resident memory is that of the slots actually touched
  (the high-water mark), as before.
* The Fortran-emulating free buffer `ifree` (`Nray_max/4` entries) is only
  allocated when used (not with `deterministic_parallel`).
* Per ray actually used: 88 B (ray) + 4 B (`ray_index`, legacy pool only).

## 5. Speedups

See section 6 of `docs/production_run.md` for the measurements (4 cores,
N = 64 and 128), which were taken while another job shared the machine.
