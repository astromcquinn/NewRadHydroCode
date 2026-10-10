# Moment RT: M1 and OTVET as alternatives to the ray tracing

*Written 2026-10-10 (branch `claude/elegant-bell-rqd3l1`, merged into `main`
the same day). New methods, not
in the original code. The default, `rt_method = rays`, is unchanged and stays
bit-identical to the Fortran: the full fast suite (93 tests, including the
live gfortran comparisons) passes unchanged.*

## 1. Usage

```
rt_driver  <problem> <out> --set rt_method=m1        # or otvet
radhydro   <problem> <out> --set rt_method=otvet
```

| Parameter | Values | Meaning |
|---|---|---|
| `rt_method` | `rays` (default), `m1`, `otvet` | radiative transfer method |
| `moment_nsub` | integer ≥ 3, default 4 | transport subcycles per RT substep (Courant number 1/`moment_nsub`) |
| `otvet_tensor` | `auto` (default), `direct`, `fft` | how the OTVET tensor is summed over the sources |

All other parameters keep their meaning. `c_reduced` sets the light speed of
the moment solver too. `legacy_physics` selects the atomic physics, as with
the rays. `legacy_rt` and `deterministic_parallel` only concern the ray
bookkeeping and have no effect on the moment methods, which always behave as
the fixed path: the direction field is reset, and τ is updated in every cell.
They need not be set. The ray-only options (`ray_dump_*`, `ray_quant_*`, the
ray pool) do nothing with the moment methods.

## 2. What it does

`src/rt_moments.cpp` (`include/radray/rt_moments.hpp`) replaces the ray
casting of each RT substep (`ray_tracing_1` and `ray_tracing_3`). It solves
the two-moment equations for the photon number density N and the flux F of
each of the 5 frequency groups, on the RT grid, at the reduced speed of light
ĉ:

```
dN/dt + ĉ div F        = S − absorption
dF/dt + ĉ² div (D N)   =   − absorption
```

**Closures.** The two methods differ only in the Eddington tensor D:

- **M1** (Levermore 1984): D is local. With f = |F|/(ĉN),
  χ = (3 + 4f²)/(5 + 2√(4 − 3f²)) and D = (1−χ)/2 I + (3χ−1)/2 u uᵀ, where
  u = F/|F|. This is the closure of RAMSES-RT (Rosdahl et al. 2013) and
  ATON (Aubert & Teyssier 2008).
- **OTVET** (Gnedin & Abel 2001): D = h, the optically thin variable
  Eddington tensor. It is h = P / tr P, with P(x) = Σₛ Lₛ d dᵀ / |d|⁴ summed
  over the source cells; d is the minimum-image separation in cells, and Lₛ
  is the cell's Luv. The source cell itself gets an isotropic Lₛ·4/3, the
  self-weight of a radius of half a cell. Off-diagonal terms with a
  component at exactly N/2 vanish, since the two periodic images average
  out. h is recomputed every hydro step, when the sources change.
  - `direct` sums in double precision.
  - `fft` is a periodic convolution with single-precision FFTW
    (`src/rt_moments_fft.cpp`). It agrees with `direct` to 3·10⁻⁷ in h
    (`test_moments`).
  - `auto` uses `direct` while N_src·N_cell ≤ 2.7·10⁸, and `fft` otherwise.

**Transport.** The scheme is the first-order Global Lax–Friedrichs (Rusanov)
flux with wave speed ĉ, unsplit, on the periodic box, as in RAMSES-RT. One RT
substep is one cell of light travel (calc_timestep), so `moment_nsub`
subcycles give a Courant number of C = 1/`moment_nsub`.

- Each cell's state splits into right- and left-going parts (N ± F_d)/2 ≥ 0,
  so the update is a positive combination for C ≤ 1/3. After every subcycle,
  N ≥ 0 and |F| ≤ N are enforced; this is the realizability of M1.
- **C = 1/3 exactly gives the centre cell zero weight.** Odd and even cells
  then decouple: a source cell is emptied every third subcycle, and
  `photo_source` skips it. The default is therefore 4 (C = 0.25, 0.75 of the
  limit, like RAMSES-RT's 0.8).
- The photon-weighted group energy E = ν N moves with the photons (flux
  E F_d / N), so each group's mean frequency is carried. That mean is what
  the rays carry per ray, and what `photo_cell` reads.

**Coupling to the physics: unchanged.**

- Sources inject fResc·pop2_init·Luv·dtsec_rt photons per substep into their
  cell, isotropically (F = 0). That is the number `ray_source` emits.
- After transport, N is deposited in `rt_grid%ngam`, E in `rt_grid%nu` and
  Σ_g F in `rt_grid%dx`.
- `radiative_physics` (`photo_cell`) is then called unchanged. N, in units
  of ⟨n_H⟩, is the photon density in the cell during the substep: exactly
  what the rays deposit, with Γ = ĉ n_γ σ, the photon limiter and the
  rho_source consumption.
- After `photo_cell`, N, F and E are scaled by the surviving photon fraction
  `rt_grid%fgam`. This is the absorption term above, with the same κ for
  both moments. As with the rays, photons below `ngam_min` in a cell are
  dropped (absorbed).

**Redshift.** In `post_rt`, each cell's 5-group spectrum gets exactly the
rays' `redshift_ray`. The code is factored out of `ray_redshifting` as
`redshift_spectrum`, with unchanged arithmetic: every stored ray reference
still passes. F is moved between groups with the photons, using the same
keep and move-down fractions.

**Photon budget.** `ray.log` gets one line per hydro step:
`Moments : in flight … injected … absorbed … lost …`, in photons per H atom
of the box. Photons are conserved to float precision: injected − absorbed −
lost = in flight to 6 digits in every run. "Lost" is the redshift out of the
lowest group plus groups below `ngam_min`. `stats.txt`
`ngam_rays_total` holds the moment photons in flight.

**Determinism and restarts.**
- Every loop writes only its own cells (the transport is a Jacobi stencil),
  and every sum is taken in a fixed order. The results are byte-identical
  for any thread count (`thread_rt_m1_*`, `thread_rt_otvet_*`, 1/2/4
  threads).
- Checkpoint version 3 holds N, F, E and the budget, and restarts are
  bit-identical (`restart_smoke_m1`, `restart_smoke_otvet`).
- Versions 1–2 still read as ray runs. A checkpoint can't be restarted with
  a different `rt_method`.

**Memory.** 20 floats per cell for the state (5 groups × N, F, E), plus 14
of scratch, plus 6 for the OTVET tensor: about 140–160 B per cell, against
88 B per ray for the rays, of which there are up to 64 per cell.

## 3. Tests

Run with `ctest -L moments`, about 40 s.

| Test | What |
|---|---|
| `moments` (`tests/test_moments.cpp`) | photon and energy conservation of the transport, realizability, light-cone spread of a pulse, an M1 plane beam moving exactly 1 cell per substep with no sideways spread, the OTVET tensor (single source: rank 1, trace 1, the N/2 images; random sources: direct vs FFT), thread independence of the transport, redshift (photons, F following, frequency range) |
| `ifront_rsla_m1`, `ifront_rsla_otvet` | the I-front against the finite-c solution (as `ifront_rsla`), checked from R = 5 cells (`--rmin 5`) |
| `thread_rt_m1_lognormal_16src`, `thread_rt_otvet_lognormal_16src` | byte-identical outputs for 1, 2 and 4 threads |
| `restart_smoke_m1`, `restart_smoke_otvet` | bit-identical restart of the full loop |

## 4. Comparison with the rays

Produced by `tools/analysis/rt_methods.py build/rt_driver WORKDIR --threads 4`,
which writes `docs/figures/rt_methods_ifront.png`, `rt_methods_maps.png` and
`rt_methods.txt`. All runs use N = 32, c_reduced = 0.1 and
legacy_rt = legacy_physics = 0.

![I-front](figures/rt_methods_ifront.png)

**I-front in a uniform medium** (one source of Luv = 8e28, z = 10):

| Method | Error at ĉt = 2 cells | at 6 cells | at 10 cells | at 30 cells |
|---|---|---|---|---|
| rays | −1% | +4.7% | +3.9% | +1.9% |
| M1 | +63% | +16% | +8.7% | +2.5% |
| OTVET | +63% | +15% | +7.0% | +1.6% |

The first-order scheme spreads the photons over a few cells: the numerical
diffusion of upwind at C = 1/4 is σ ≈ 2 cells after 6 substeps. The moment
fronts therefore run ahead of the finite-c solution while they are small,
and converge to it, and to the rays, as R grows. The photon budget is exact
(§2), so the ionised volume is right once the front is resolved.

![maps](figures/rt_methods_maps.png)

Maps of x_HII in the plane through the sources, after 15 steps:

- **Shadow.** A 3× brighter source, and a clump of radius 2 cells at 50× the
  mean density, 4 cells away.
  - All three methods cast a shadow. The mean x_HII on the axis behind the
    clump is 0.000 for the rays, 0.061 for M1 and 0.008 for OTVET.
  - At the same distance beside the clump it is 0.52, 0.39 and 0.32.
  - OTVET's radial tensor keeps the shadow sharpest. M1 fills it in a little:
    its flux direction relaxes towards isotropy behind the obstacle.
  - The rays' front shows neutral holes at large radius, where the HEALPix
    sampling gets sparse.
- **Two sources** 12 cells apart. The fronts agree: the width of the ionised
  region through the midpoint is 7 cells with every method, and through a
  source it is 13, 13 and 15 cells. The classic M1 failure, where crossing
  beams merge, needs beams that cross; two isotropic point sources do not
  provoke it at this resolution.
- **16 sources in a lognormal field** (`tests/problems/lognormal_16src`,
  60 steps):

| step | x_HII (volume): rays / M1 / OTVET | x_HII (mass) | partly ionised cells | cells off by > 0.1 vs rays |
|---|---|---|---|---|
| 10 | 0.021 / 0.058 / 0.056 | 0.084 / 0.093 / 0.094 | 2% / 7% / 7% | – / 7% / 7% |
| 20 | 0.101 / 0.149 / 0.142 | 0.152 / 0.169 / 0.169 | 6% / 13% / 13% | – / 12% / 10% |
| 40 | 0.310 / 0.359 / 0.344 | 0.293 / 0.315 / 0.310 | 12% / 21% / 20% | – / 18% / 16% |
| 60 | 0.525 / 0.569 / 0.552 | 0.434 / 0.456 / 0.450 | 15% / 24% / 24% | – / 22% / 19% |

  The moment methods start ahead, while the fronts are a few cells across,
  as in the I-front test. Later they are 5–8% (M1) and 4–5% (OTVET) above
  the rays in x_HII. Their fronts are thicker: about 1.6× as many partly
  ionised cells. The median temperatures agree at step 60 (16 000–16 300 K).
  OTVET is the closer of the two to the rays everywhere.

**Cost** (the box, 60 steps of one RT substep each, whole run including
output; 4-core container):

| | 1 thread | 4 threads |
|---|---|---|
| rays (`deterministic_parallel = 1`) | 25.5 s | 11.0 s |
| M1 | 3.9 s | 2.1 s |
| OTVET | 4.6 s | 2.1 s |

The moment transport costs about 0.7 µs per cell per substep with one
thread (5 groups × 4 subcycles). It does not depend on the number of
sources or rays. Most of the rays' time here is `post_rt` (redshifting every
ray, every hydro step) and the ray casting. In a production run with many
substeps per hydro step and up to 64 rays per cell the ratio will differ:
measure it on the cluster. Cells whose stencil holds no photons are skipped.

## 5. Caveats and open questions

- **First-order diffusion.** Fronts are thick (about 2–3 cells), and small
  fronts run ahead (§4). A second-order (MUSCL + RK2) transport would reduce
  this, but needs realizability-preserving limiting. Not done.
- **M1 with crossing beams** (Rosdahl et al. 2013, Fig. 3; Aubert & Teyssier
  2008) is not tested here. The point-source problems above do not
  provoke it.
- **OTVET formulation.** This is the time-dependent two-moment system with
  the OTVET tensor, solved explicitly with the same GLF machinery as M1.
  Published OTVET implementations differ:
  - Gnedin & Abel 2001 combine the moment equations differently;
  - Petkova & Springel 2009 drop ∂F/∂t and solve the resulting diffusion
    equation implicitly.

  The papers could not be consulted (no network access to journals), so
  these remarks are from memory and should be checked. The tensor itself
  (1/r² weighted, optically thin) is the defining OTVET ingredient.
  `kSelfWeight` (the source cell's own weight, 4) is a choice; it only
  matters where several sources compete within a few cells.
- **No recombination radiation** (as with the rays, where it is off).
- **Not validated against the Iliev et al. 2006/2009 tests** quantitatively.
  That would be the next step, with plane-parallel sources for the shadow
  test (Test 3).
- **AMR, MPI and GPU** (ROADMAP): the moment solver is a 7-point stencil on
  the grid, with a ghost layer of one cell per subcycle. It is a natural fit
  for every one of them, and for the "moment hand-over in transparent
  regions" of ROADMAP §4.5.

## 6. Production comparison on Hyak (256³)

*Run on klone 2026-10-09/10 (branch `claude/elegant-bell-rqd3l1`, one
whole cpu-g2 node, 192 threads per run).* Two sets, each with three runs that
share the ICs and every parameter except `rt_method`:
- `reion256_c01`: c_reduced = 0.1.
- `reion256`: c_reduced = 0.01.

The ray runs are the existing production runs (`docs/cluster.md` §4–§7).
Every run uses `legacy_rt = legacy_physics = 0`, 10 Mpc/h, 256³, on-the-fly
sources, and goes from z = 100 to z = 5.4.

### 6.1 Reproduce

```
# on klone, repository at this branch (or main after the merge)
sbatch tools/hyak/build_hyak.sh                     # build/ + fast suite
export RUNS=/gscratch/astro/$USER/radhydro/runs
tools/hyak/rt_methods_runs.sh                       # DRYRUN=1 to only print
#  -> $RUNS/reion256_c01_{m1,otvet}, $RUNS/reion256_{m1,otvet}: params.txt +
#     rt_method, ics -> the base ICs; one whole-node job each
#     (reion_production.slurm: auto-restart, DONE, reion.py, movie)
sbatch --dependency=afterany:<the four jobs> tools/hyak/rt_methods_post.slurm
#  -> $RUNS/rt_methods_{c01,c001}/: compare.txt, history.png, maps_z*.png,
#     compare.mp4 + compare_stills/ (tools/movies/compare_movie.py)
```

- The fast suite on klone: 90 of 100 pass. The 10 failures are the known
  1-ulp platform differences of the stored RT references (`docs/cluster.md`
  §3): its 8, plus `ray_dump_{legacy,fixed}_lognormal_16src_smallbuf`, which
  compare with the same `lognormal_16src_smallbuf` references and show the
  same mismatches. Every `moments`, restart and thread test passes.
- `rt_methods_runs.py` compares the fields at the **movie frames** nearest
  to the reference's zdata1 redshifts. The zdata1 outputs themselves drift
  apart once the Courant step feels the different temperatures: by up to
  Δz = 0.1 after z ≈ 11 at 0.01c.
  - The movie frames are written at the same `movie_dz` redshifts from
    inside the RT substeps. They agree to Δz ≤ 0.003 at 0.1c, and to within
    one RT substep (Δz ≤ 0.05) at 0.01c. The `z_run` column of
    `compare.txt` gives each run's frame.
  - ⟨x⟩_M comes from `frames.txt`, since the frames carry no `rhog`.

Full tables: [`figures/rt_methods_256_c01.txt`](figures/rt_methods_256_c01.txt),
[`figures/rt_methods_256_c001.txt`](figures/rt_methods_256_c001.txt).

### 6.2 Movies

Side by side (rays | M1 | OTVET), with slices of x_HII and T at the box
centre and x_HII,V(z) of all three runs:
- [`movies/rt_methods_256_c01.mp4`](movies/rt_methods_256_c01.mp4)
  (0.1c, 173 frames, 2.6 MB);
- [`movies/rt_methods_256_c001.mp4`](movies/rt_methods_256_c001.mp4)
  (0.01c, 171 frames, 3.1 MB).

Each run's own movie (x_HII, T, ρ and τ_es) is
`$RUNS/<run>/run/movie/cosmo.mp4`.

![0.1c, z = 10.35](figures/rt_methods_256_c01_z10.35.png)

*0.1c at z = 10.35, x_HII,V ≈ 0.51.*

![0.01c, z = 8.94](figures/rt_methods_256_c001_z8.94.png)

*0.01c at z = 8.94, x_HII,V ≈ 0.50.*

### 6.3 Reionization history

![history 0.1c](figures/rt_methods_256_c01_history.png)

*0.1c: x_HII (left) and the difference from the rays in x_HII,V (right).*

![history 0.01c](figures/rt_methods_256_c001_history.png)

*Same for 0.01c.*

| set | run | z(x_V = 0.1) | 0.5 | 0.9 | 0.99 | z(x_M = 0.5) | τ_es |
|---|---|---|---|---|---|---|---|
| 0.1c | rays | 12.639 | 10.379 | 9.248 | 8.826 | 10.565 | 0.0904 |
| | M1 | 12.645 | 10.373 | 9.262 | 8.909 | 10.573 | 0.0905 |
| | OTVET | 12.624 | 10.367 | 9.282 | 8.897 | 10.577 | 0.0906 |
| 0.01c | rays | 11.645 | 8.960 | 7.291 | 6.494 | 9.352 | 0.0761 |
| | M1 | 11.710 | 9.029 | 7.309 | 6.517 | 9.416 | 0.0768 |
| | OTVET | 11.679 | 8.934 | 7.217 | 6.373 | 9.373 | 0.0764 |

- **The global history is the same to |Δx_HII,V| ≤ 0.018 at every z.**
  - At 0.1c: up to x_V = 0.9, z agrees to 0.04, and τ_es to 0.2%. The
    largest differences, ±0.012, come in the last stage of overlap
    (z ≈ 9–10). The end of reionization (x_V = 0.99) is Δz = 0.07–0.08
    later with the moment methods.
  - At 0.01c: M1 runs ahead by up to +0.018 (z ≈ 8.4), and OTVET falls
    behind by up to −0.015 (z ≈ 7.8). τ_es agrees to 1%. OTVET finishes
    Δz = 0.12 earlier at x_V = 0.99.
- **Small early fronts (§4) show at 0.01c only.** At z = 13.12 the moment
  runs have 17–20% more ionised volume (x_V = 0.0199 and 0.0195 against
  0.0166). They also have 1.4–1.5× the x_HII power at k = 2, and the median T
  of their ionised gas is 13% lower (15 100–15 200 K against 17 400 K).
  - The difference shrinks to 2–5% by z = 11.5.
  - At 0.1c the fronts are already larger at that redshift (x_V = 0.057),
    and the three runs agree to 2%.

### 6.4 Morphology and temperature

| set | z | x_V (rays) | cells with \|Δx\| > 0.1: M1 / OTVET | corr(log T): M1 / OTVET | partly ionised: rays / M1 / OTVET |
|---|---|---|---|---|---|
| 0.1c | 12.30 | 0.137 | 5.1% / 2.5% | 0.902 / 0.977 | 1.5 / 2.5 / 2.3% |
| | 10.90 | 0.364 | 13.3% / 5.9% | 0.818 / 0.959 | 4.4 / 5.3 / 4.6% |
| | 9.80 | 0.700 | 19.9% / 7.9% | 0.597 / 0.912 | 5.4 / 5.0 / 4.8% |
| | 9.00 | 0.965 | 6.0% / 2.2% | 0.083 / 0.796 | 1.0 / 0.9 / 0.8% |
| | 6.80 | 1.000 | – | 0.934 / 0.989 | – |
| 0.01c | 11.54 | 0.108 | 2.9% / 2.6% | 0.962 / 0.976 | 2.2 / 3.1 / 2.9% |
| | 8.99 | 0.494 | 10.3% / 7.0% | 0.889 / 0.951 | 4.2 / 5.9 / 5.3% |
| | 7.44 | 0.870 | 9.9% / 5.8% | 0.674 / 0.893 | 3.5 / 3.5 / 3.5% |
| | 6.85 | 0.966 | 4.2% / 2.4% | 0.541 / 0.867 | 1.4 / 1.5 / 1.5% |

- **OTVET matches the rays' ionisation and temperature morphology much better
  than M1.**
  - It has 2–2.5× fewer cells off by more than 0.1 in x_HII during the
    middle of reionization.
  - Its log T correlates with the rays' at 0.91–0.99 at 0.1c and 0.87–0.98
    at 0.01c.
  - For M1 the correlation drops through overlap (0.60 at x_V = 0.7, at
    0.1c).
  - This is consistent with the known weakness of M1 where radiation from
    many sources crosses: the local closure gives a single flux direction
    per cell. That is an interpretation, not tested here.
- **The last neutral islands are in different places** (0.1c, z = 9.0,
  x_V = 0.965–0.977, below). In the M1 run they have moved so far that its
  log T correlates with the rays' at only 0.08. OTVET's are closer: 0.80.
  - Tests that depend on where the islands survive late, such as the
    Lyα forest after overlap, are therefore sensitive to the RT method at
    this level.
- **Fronts.** The moment fronts are smooth, while the rays' are ragged:
  - the rays' HEALPix sampling leaves holes at large radius, worst at 0.01c
    (§4);
  - the moment methods have up to 1.4× as many partly ionised cells
    (5.9% against 4.2%) while the bubbles grow;
  - at overlap the fraction of partly ionised cells is the same for all
    three methods.
- **The temperature of ionised gas agrees to 1–2%** from z ≈ 11 on. Medians
  for rays / M1 / OTVET:
  - 0.1c, z = 6.80: 6194 / 6290 / 6225 K;
  - 0.1c, z = 5.45: 5042 / 5123 / 5099 K;
  - 0.01c, z = 7.44: 9584 / 9473 / 9648 K.

  The moment runs are slightly hotter after overlap at 0.1c.
- **x_HII power after overlap.** At 0.1c, M1's residual neutral-fraction
  power is 0.82–1.09 of the rays' at k = 2–16, and OTVET's is 0.58–0.89.
  OTVET's residual x_HI fluctuations are smoother. With x_V > 0.999 (0.01c,
  z ≤ 5.9) the ratio only says whether a last patch survives in the box,
  and is meaningless.

![last islands, 0.1c](figures/rt_methods_256_c01_maps_z9.0.png)

*0.1c at z = 9.0: the last neutral regions (dark) in the middle plane.*

### 6.5 Cost

Whole node, 192 threads. The RT time is the sum of the RT timers in
`run.log`. Wall time is the whole run, including hydro, PM, FoF and output.

| set | run | substeps | wall | RT | transport | photo-chemistry | other RT | peak RSS |
|---|---|---|---|---|---|---|---|---|
| 0.1c | rays | 4064 | 1:32:42 | 5032 s | 3552 s | 1382 s | 98 s | 27.0 GB |
| | M1 | 4037 | 49:35 | 2823 s | 1206 s | 1549 s | 68 s | 13.5 GB |
| | OTVET | 4035 | 51:39 | 2952 s | 1229 s | 1553 s | 170 s | 14.3 GB |
| 0.01c | rays | 404 | 26:17 | 882 s | 641 s | 181 s | 60 s | 23.1 GB |
| | M1 | 403 | 6:50 | 265 s | 116 s | 136 s | 13 s | 14.2 GB |
| | OTVET | 401 | 8:28 | 366 s | 116 s | 132 s | 118 s | 14.4 GB |

What the columns contain:
- **transport:** for the rays, `calc_ray_1 + ray_tracing_1 + ray_tracing_3`.
  For the moment methods, `ray_tracing_1` (the GLF subcycles and the
  deposit) plus `ray_tracing_3` (scaling by the surviving fraction).
- **photo-chemistry:** `radiative_physics`.
- **other RT:** `prep`, `calc_rt_1`, `sync` and `post_rt`. OTVET's `prep`
  is the FFT tensor, 106 s in all for the hydro steps with sources.

Findings:
- **Speed-up over the rays: 1.8× in RT time at 0.1c** (1.9× wall for M1,
  1.8× for OTVET), and 2.4–3.3× at 0.01c (3.1–3.9× wall).
- **Memory:** about half the rays' at 0.1c.
- **The transport met the estimate; the total did not.** The moment
  transport costs 0.29–0.30 s per substep at both light speeds (0.26–0.27 s
  of it in `ray_tracing_1`). That is inside the 0.15–0.3 s expected from
  the 32³ timings, and 2.9× cheaper than the ray tracing per substep at
  0.1c.
  - It is 3.1 µs of thread time per cell per substep, against 0.7 µs on
    one thread at 32³: the 192 threads give a 44× speed-up, at the low end
    of the 40–80× guessed for memory-bound scaling.
  - The 10–20 min estimate counted only the transport, though.
    `radiative_physics` (`photo_cell`, unchanged) costs as much as with the
    rays (0.33–0.45 s per substep) and is now **51–55% of the moment RT
    time**.
  - A faster RT now needs faster photo-chemistry. Its cost does not depend
    on the RT method.
- **Time per substep follows the ionised volume, not the number of sources**
  (wall time per hydro step ÷ substeps, 0.1c):

  | z | 11.6 | 9.8 | 8.6 | 7.8 | 6.8 | 5.5 |
  |---|---|---|---|---|---|---|
  | x_V (rays) | 0.24 | 0.70 | 1.00 | 1.00 | 1.00 | 1.00 |
  | N_source | 278 | 1020 | 2140 | 3200 | 4500 | 5300 |
  | rays [s] | 1.15 | 1.42 | 1.56 | 1.56 | 1.45 | 1.55 |
  | M1 [s] | 0.39 | 0.70 | 0.86 | 0.88 | 0.88 | 0.88 |
  | OTVET [s] | 0.40 | 0.71 | 0.90 | 0.92 | 0.92 | 0.95 |

  - The moment cost per substep doubles between x_V = 0.24 and overlap,
    then stays flat while N_source grows 2.5×.
  - The transport skips cells whose stencil holds no photons (§4). The
    rise is consistent with more cells holding photons, for the transport
    and the photo-chemistry, until photons fill the box at overlap. Not
    checked with per-step timers.
  - OTVET's tensor grows only with the source count (`auto` → FFT here, a
    fixed cost per hydro step).
- **Photon budget** (M1, 0.1c, `ray.log` at the end): injected 11.38,
  absorbed 2.54, lost 2.47 and in flight 6.36 photons per H atom. That
  closes to the printed digits.

### 6.6 Conclusions

- At 256³ in a 10 Mpc/h box, both moment methods reproduce the rays'
  reionization history to |Δx_V| ≤ 0.02 and τ_es to 1%. They take
  roughly half the RT time and half the memory at 0.1c.
- **OTVET is the better match** for morphology and temperature: 2–3× fewer
  discrepant cells than M1, and corr(log T) ≥ 0.87 throughout.
  - M1's last neutral islands and its temperature structure near overlap
    differ substantially from the rays'.
  - OTVET costs 5% more than M1 at 0.1c (the FFT tensor). At 0.01c, with
    few substeps per hydro step, it costs 38% more.
- **Next for speed:** `radiative_physics`. Beyond that, a second-order
  transport would sharpen the fronts (§5).
