# ACCESS supplement request, PHY250179: progress report and plan (draft, under three pages)

*Drafted 2026-10-22 for the user to edit and submit ("Request a supplement"
on the project page; 1.5 million credits). Bracketed items are for the
user. Sources: the three arXiv listings, the Willard et al. draft
(2026-10-22), `docs/zoom_turbulence.md` §3, §3b, §7, `docs/stage_f.md` §4.*

**Project.** PHY250179, "High-resolution cosmological simulations of
intergalactic gas" (Accelerate; NASA grant 80NSSC24K1220). PI Matthew
McQuinn, University of Washington. Resources used: [SDSC Expanse; add
any other system], [credits spent / remaining].

## 1. Progress to date

The allocation supported high-resolution simulations of the
intergalactic medium (IGM) at z = 3–7, in two directions: the response
of IGM gas to reionization at the kiloparsec scale, and the structure of
cold gas in the cosmic web. Four papers resulted, three on arXiv and one
about to be submitted.

**(a) SAGUARO: a suite of IGM reionization simulations** (Cain, Das,
D'Aloisio, Foreman, Scannapieco, Moreno, Lugatiman, Cohon, Tohfa & Trac
2026, arXiv:2603.25788). Over two hundred coupled radiative-hydrodynamics
simulations of IGM gas during and after reionization, in boxes of 2 and
0.25 h⁻¹ Mpc at resolutions that resolve the Jeans scale of the cold,
neutral IGM, spanning a grid of photoionization rates, reionization
redshifts and box-scale densities, with X-ray preheating, recombination
radiation, baryon–dark matter streaming and alternative dark-matter
cosmologies. The paper sets out the suite and its first results on the
gas dynamics and thermal structure, the ionizing opacity and
self-shielding, the IGM's drain on the reionization photon budget, the
halo mass function and Lyα transmission. [State which of these runs used
this allocation and the credits they took.]

**(b) An emulator for the ionizing-photon mean free path** (Tohfa, Cain,
McQuinn & D'Aloisio 2026, arXiv:2602.03923; JCAP 06 (2026) 094). A
deep-learning emulator trained on the high-resolution radiative-transfer
simulations of ionization fronts in ≈ 2 h⁻¹ Mpc volumes, which captures
the kiloparsec-scale clumping that large reionization volumes cannot
resolve. Combined with the measured photoionization rate and mean free
path at z = 5–6, and without power-law assumptions, it implies a factor
2–3 decline in the ionizing emissivity between z = 6 and 4.8; it also
serves as a sub-grid opacity for large simulations.

**(c) A self-consistent analytic model of the photoionization rate and
the reionization history** (Cain, Croker, D'Aloisio, Georgiev, Tohfa, Zhu
& Windhorst 2026, arXiv:2606.19449). An analytic formalism derived from
the cosmological radiative-transfer equation that predicts the ionized
fraction and Γ_HI together, tested against the radiative-transfer
simulations: percent-level in the ionized fraction and 20–30 % in Γ_HI at
z ≲ 6, which makes Γ_HI from the Lyα forest usable in fast inference of
the reionization history.

**(d) Multiphase gas in high-z cosmic sheets** (Willard, Faerman,
McCarty, McQuinn, Mandelker & Quinn, to be submitted). A suite of
moving-mesh (AREPO) simulations of individual z = 3–6 cosmic sheets, 4
and 8 cMpc/h boxes at 2 × 128³ to 2 × 1024³ particles (cold-gas cells of
0.4–3 pkpc), with an SPH comparison suite in MP-Gadget; all runs at the
San Diego Supercomputer Center. The sheets' cold cores fragment through
the nonlinear thin-shell instability and then Kelvin–Helmholtz, not the
thermal instability proposed before; fragmentation grows with resolution
and is not converged even at the highest; the covering fractions of
partial Lyman-limit and Lyman-limit systems rise by factors of 2–4 across
the resolution range; and cosmic sheets plausibly supply of order 1 % of
LLSs and 10 % of pLLSs at z = 4, a large share of the metal-free
absorbers. [Credits used.]

**Software built with the allocation's science in view.** A C++ port of
Trac's RadHydro (the code behind (a)–(c)), validated bit for bit against
the original, now runs with MPI + OpenMP through AMReX on 192-core nodes
(512³ reionization runs on five nodes, deterministic across rank counts),
with adaptive mesh refinement for gas, radiation and gravity validated
to two levels, a star-formation sink, and a first refined 512³
production run whose cost is measured. This development ran on the
group's own cluster; the supplement funds the runs it enables.

## 2. Why more credits: turbulence driven by real ionization fronts

Cain, McQuinn, Scannapieco, D'Aloisio & Trac (2026, PRD 113, 023050)
showed with uniform 0.25 h⁻¹ Mpc boxes at 1024³ and 2048³ that the
heating by ionization fronts drives kiloparsec-scale turbulence once
cells are smaller than ≈ 0.5 kpc/h, and that this turbulence can grow
volume-filling magnetic fields to the strengths the TeV-blazar limits
require. Those boxes were ionized by plane-parallel fronts at a chosen
redshift. The new code makes two things possible that were not when the
project was written:

1. **The same boxes in the new code (C2):** 1024³ at 244 pc/h with
   plane-parallel fronts, then preheating (T_min = 100 K, 1000 K) and
   reionization-redshift variants; the code-to-code check and the
   baseline.
2. **Zoom-in runs of one region of a full reionization simulation (Z2,
   Z3):** a 40 h⁻¹ Mpc, 512³ box with sources and escape fraction
   calibrated to the Lyα forest, in which one region of 1.25 h⁻¹ Mpc
   (1024³ cells at 1.2 kpc/h) and then 0.31 h⁻¹ Mpc (1024³ at 0.3 kpc/h)
   is refined by six and eight levels, with initial conditions seeded at
   the region's resolution. The region then meets the real front: its
   arrival time, speed, direction and spectrum, and the surrounding
   filaments and bulk flow. Regions are chosen by density and by
   ionization time. This gives the turbulence driven by a real front
   and the first self-consistent setting for the Biermann battery at
   the front.

## 3. Resources requested

[Resource and conversion from the ACCESS exchange calculator.] From the
measured costs of this code on 192-core AMD Genoa nodes (one hydro
cell-update ≈ 3 µs with refinement, one ray-substep ≈ 2 µs) and the
refined 512³ production run (2.6× the unrefined step time):

| Run | Finest-level cells | Core-hours |
|---|---|---|
| C2: 1024³, 0.25 h⁻¹ Mpc, z = 300 → 4 | 1.1 × 10⁹ | ≈ 1 × 10⁵ |
| C2 variants: T_min = 100 K, 1000 K; a second z_re | | ≈ 3 × 10⁵ |
| Z2: 512³ base + 6 levels, 1024³ over 1.25 h⁻¹ Mpc | 1.1 × 10⁹ | ≈ 3 × 10⁴ |
| Z3: the same region, 8 levels, 1024³ over 0.31 h⁻¹ Mpc | 1.1 × 10⁹ | ≈ 1 × 10⁵ |
| a second Z3 region (a different ionization time) | 1.1 × 10⁹ | ≈ 1 × 10⁵ |
| resolution ladders (64³–512³) and analysis | | ≈ 1 × 10⁵ |
| **Total** | | **≈ 7 × 10⁵ core-hours**, requested as 1.5 million credits with the conversion and a margin for the refinement overhead |

Memory under 1 TB per run; checkpoints ≈ 0.4 TB every ten steps; the
finest-level cubes 4 GB per field per snapshot. Storage: [to fill in].

## 4. Timeline

[Months 1–2: C2 and variants. Months 2–5: Z2, Z3, the second region.
Month 6: analysis, and the renewal request for the full suite (several
regions × preheating, the MHD runs).]
