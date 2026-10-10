# Moment RT: M1 and OTVET as alternatives to the ray tracing

*Written 2026-10-10 (branch `claude/elegant-bell-rqd3l1`). New methods, not
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

Not yet run; nothing here was measured on klone. To compare the methods on
the production boxes, run M1 and OTVET siblings of the existing ray runs
(`reion256_c01`, 0.1c; `reion256`, 0.01c; `docs/cluster.md` §4–§7). They
share the ICs and every parameter except `rt_method`.

```
# on klone, repository at this branch (or main after the merge)
sbatch tools/hyak/build_hyak.sh                         # build/ + fast suite
tools/hyak/rt_methods_runs.sh                           # DRYRUN=1 to only print
#  -> $RUNS/reion256_c01_{m1,otvet}, $RUNS/reion256_{m1,otvet}: params.txt +
#     rt_method, ics -> the base ICs; one whole-node job each
#     (reion_production.slurm: auto-restart, DONE, reion.py, movie)
# when they are done:
python3 tools/analysis/rt_methods_runs.py \
    $RUNS/reion256_c01/run:rays $RUNS/reion256_c01_m1/run:M1 \
    $RUNS/reion256_c01_otvet/run:OTVET --outdir $RUNS/rt_methods_c01
python3 tools/analysis/rt_methods_runs.py \
    $RUNS/reion256/run:rays $RUNS/reion256_m1/run:M1 \
    $RUNS/reion256_otvet/run:OTVET --outdir $RUNS/rt_methods_c001
```

`rt_methods_runs.py` writes:
- `history.png` and the z(x_HII) table with τ_es (from `movie/frames.txt`);
- for each `zdata1` output common to the runs: x_HII and T slices, the
  volume- and mass-weighted x_HII, partly ionised cells, the median T of
  ionised gas, cells off by more than 0.1 from the rays, the correlation of
  log T, and P_x(k) relative to the rays;
- the RT time and peak memory from `<run>.log`.

Expected cost (an extrapolation from the 32³ timings, to be checked):
- The transport is about 0.7 µs per cell per substep on one thread, about
  12 s per substep at 256³. With memory-bound scaling of perhaps 40–80× on
  192 threads, that is 0.15–0.3 s per substep: 10–20 min for the 4064
  substeps of the 0.1c run (the rays took 1 h 32 min), and a few minutes at
  0.01c (404 substeps).
- The OTVET tensor is about 7 serial 256³ FFTs per hydro step.
- Memory is about 3 GB of moment state, and checkpoints are 1–2 GB (no
  rays).

Things to look at:
- **Small fronts.** At 256³ the early fronts around the first, faint sources
  are a few cells across, where the first-order diffusion runs ahead (§4).
  Expect earlier ionisation at the start and a smaller difference later.
- **The 0.01c runs.** There the light-travel resolution per hydro step is
  coarse for both methods.
