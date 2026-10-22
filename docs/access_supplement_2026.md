# Draft: ACCESS supplement request for PHY250179 (1.5 million credits)

*Drafted 2026-10-22 for the user to edit and submit through the ACCESS
portal ("Request a supplement" on the project page). The planned-work
part comes from `docs/zoom_turbulence.md` (§1, §3, §3b, §7); the
progress part and the resource choice are placeholders the user fills
in. Keep it short: a supplement is reviewed on a rolling basis against a
brief justification.*

## Project

PHY250179, "High-resolution cosmological simulations of intergalactic
gas" (Accelerate, NASA grant 80NSSC24K1220). PI: Matthew McQuinn,
University of Washington.

## Progress to date (to fill in)

[What the 3 million credits bought so far: the runs, the resources used,
the papers or preprints they produced, and the credits remaining. If
Cain et al. 2026 (PRD 113, 023050) or its follow-ups ran under this
project, say so and cite it.]

## Why more credits are needed

Cain et al. (2026) showed, with uniform 0.25 h⁻¹ Mpc boxes at 1024³ and
2048³, that the heating by ionization fronts drives kiloparsec-scale
turbulence in the intergalactic medium (IGM) once the cell size is below
about 0.5 kpc/h, and that this turbulence can grow volume-filling
magnetic fields to the strengths the TeV-blazar limits require. Those
boxes were ionized by plane-parallel fronts at a chosen redshift, with no
information about where the gas sits in the cosmic web or when and from
which direction the real fronts arrive. The requested credits fund the
next step on both fronts, with a code that did not exist when PHY250179
was written:

1. **The same boxes, in the new code (C2).** A C++ port of Trac's
   RadHydro, validated bit for bit against the original and now parallel
   through AMReX (MPI + OpenMP, 512³ runs over five 192-core nodes),
   reproduces the Cain et al. configuration with plane-parallel fronts:
   1024³ at 244 pc/h, then the preheating variants (T_min = 100 K and
   1000 K) and two reionization redshifts. This is the code-to-code
   check and the baseline for the zooms.
2. **Zoom-in runs of one region in a full reionization simulation (Z2,
   Z3).** A 40 h⁻¹ Mpc, 512³ box with the sources, escape fraction and
   speed-of-light schedule calibrated to the Lyα forest, in which one
   region of 1.25 h⁻¹ Mpc (1024³ cells at 1.2 kpc/h) and then 0.31 h⁻¹ Mpc
   (1024³ at 0.3 kpc/h) is refined by six and eight levels of mesh
   refinement, with the initial conditions seeded at the region's
   resolution (multi-mass particles, nested phases). The region then sees
   the real front: its arrival time, speed, direction and spectrum, and
   the surrounding filaments and bulk flow. The regions are chosen by
   density and by ionization time from the base run. The result is the
   turbulence driven by a real front, and the first self-consistent
   setting for the Biermann battery at the front.

## Resources requested

[Resource: e.g. the CPU system used so far; credits converted with the
ACCESS exchange calculator at submission.]

Estimated cost, from the measured per-cell and per-ray costs of this code
on 192-core AMD Genoa nodes (one hydro cell-update ≈ 3 µs with mesh
refinement, one ray-substep ≈ 2 µs) and from the first refined 512³
production run (2.6× the unrefined step time):

| Run | Cells at the finest level | Core-hours |
|---|---|---|
| C2: 1024³ uniform, 0.25 h⁻¹ Mpc, z = 300 → 4 | 1.1 × 10⁹ | ≈ 10⁵ |
| C2 variants: T_min = 100 K, 1000 K; a second z_re | | ≈ 3 × 10⁵ |
| Z2: 40 h⁻¹ Mpc base + 6 levels, 1024³ over 1.25 h⁻¹ Mpc | 1.1 × 10⁹ | ≈ 3 × 10⁴ |
| Z3: the same region, 8 levels, 1024³ over 0.31 h⁻¹ Mpc | 1.1 × 10⁹ | ≈ 10⁵ |
| a second Z3 region (a different ionization time) | 1.1 × 10⁹ | ≈ 10⁵ |
| the resolution ladders below each (64³ to 512³), analysis | | ≈ 10⁵ |
| **Total** | | **≈ 7 × 10⁵ core-hours**, requested as 1.5 million credits with the conversion and a margin for the mesh-refinement overhead measured at 2.6× |

Memory: under 1 TB per run (≈ 0.5 kB per cell, 64 B per particle, the
rays); checkpoints of ≈ 0.4 TB every ten steps; the finest-level cubes
of the region 4 GB per field per snapshot. Storage: [to fill in].

## Readiness

* The code reproduces the original RadHydro bit for bit on every
  validation problem, with the physics and ray-tracing fixes behind
  flags.
* Parallel through AMReX: 512³ reionization runs on 1 to 5 nodes, with
  the RT scaling measured (2.0× from 1 to 4 nodes at 256³, 1.9× from 2
  to 4 at 512³), bitwise identical across rank counts.
* Mesh refinement for gas, radiation and gravity validated to two levels
  against uniform runs at the matched resolution; the first refined 512³
  production run has finished and its cost is measured.
* The plane-parallel front mode, the zoom initial conditions (MUSIC) and
  the multigrid gravity on the refined levels are the three pieces being
  built now, on the group's own cluster; the supplement funds the runs,
  not the development.

## Timeline

[Months 1–2: C2 and its variants. Months 2–5: Z2, then Z3 and the second
region. Month 6: analysis and the renewal request for the full suite.]
