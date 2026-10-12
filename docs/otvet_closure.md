# A better Eddington tensor for OTVET

*Branch `otvet-closure`, started 2026-10-12. Work in progress: steps 1
(ground truth) and 2 (candidates and ranking) of the user's plan; nothing
is implemented in the solver yet.*

**The problem** (`docs/rt_methods_512.md` §5–6). After overlap, OTVET's
P_Γ is 2–9× the rays' at k ≲ 6 h/Mpc with r ≥ 0.95, its Γ has a floor at
0.62⟨Γ⟩ and its peaks are too low. Before overlap it matches the rays'
morphology and temperatures better than M1. The tensor is optically thin,
h = P/tr P with P = Σₛ Lₛ d dᵀ/|d|⁴ over minimum-image separations (FFT).

## 1. Ground truth: the rays' Eddington tensor

**`tensor_dump = 1`** (radhydro, serial backend; off by default, so a
null pointer and byte-identical) accumulates, on the hydro steps that write
a zdata1 snapshot, the Γ_HI-weighted angular moments of the radiation in
every cell, summed over the RT substeps of that step:

```
S0 = Σ w,   S1 = Σ w n,   S2 = Σ w n nᵀ,    w = Σ_g σ_HI(ν_g) n_g
```

- **Rays:** every deposit of a ray segment (`rrk::cast_1`), with the ray's
  direction n. The thick-source deposit counts as isotropic. Recombination
  rays are included.
- **Moment methods:** each group's N, F and D N after the transport, with D
  the closure actually used (OTVET h or M1's D(f)).
- h_true = S2/S0 (trace 1), f = S1/S0 (the reduced flux F/ĉN).
- Written per snapshot step: `tensor.bin` (10 floats per cell, means over
  the substeps), `tensor_odd.bin` (the odd substeps only: the difference of
  the odd and even halves estimates the rays' sampling noise),
  `otvet_h.bin` (the thin tensor of the step's sources, as the code
  computes it) and `sources.txt` (source cells and Luv; the header gives the
  photons per substep per unit Luv per group and σ_HI(ν_g), so the direct
  light of the sources can be put in the rays' units).
- Tests: the outputs other than the dump are byte-identical with and without
  it, the dump is identical for 1 and 4 threads, tr h = 1 and |f| ≤ 1 to
  float precision; fast suite 108/108.

**Runs.** `$RUNS/otvcl256{,_m1,_otvet}`: 256³ in 20 Mpc/h, the cell
(78 kpc/h), mass resolution and configuration of the 512³ comparison (case
A + RR with `recomb_rad_cell_onspot = 1`, ĉ = 0.2c, f_esc = 0.23,
`redshift_per_distance`, unsplit hydro), dumps at z = 10, 9, 8.5, 8, 7.5,
7, 6.75, 6.5, 6.25, 6, 5.75, 5.5, 5.

**Analysis.** `tools/analysis/tensor_compare.py` evaluates closures *a
priori*, on the rays' own radiation field (the tensor each closure would
give in the rays' state), against h_true: the S0-weighted Frobenius
distance |h − h_true| and the anisotropy A = |h − I/3|/√(2/3) (0 isotropic,
1 one beam), by regime (fronts, ionised gas by Γ/⟨Γ⟩, shadowed and diffuse
cells, distance to the nearest source) and on 8³ blocks (the large-scale
tensor behind P_Γ at low k, nearly free of ray noise).

*Results: to come (§3).*

## 2. Literature

Read through the web (arXiv is not reachable from klone); quotes and
details in [`otvet_closure_literature.md`](otvet_closure_literature.md)
(tagged READ / INFER / RECALL), checked only as far as the HTML versions
allowed. Items marked *(memory)* were not verified.

- **Wu, McQuinn & Eisenstein 2021** (JCAP 02, 042; arXiv:2009.07278):
  linear theory of the Γ fluctuations after overlap with moment closures.
  OTVET implementations sum the sources only out to half a box, so the
  tensor isotropises on the scale L/2 instead of λ_mfp: "For box sizes
  smaller than twice the mean free path, the OTVET Eddington tensor is more
  anisotropic than the exact Eddington tensor", and simulations "likely
  overpredict fluctuations in the ionizing background by factors of 2−4".
  M1 suppresses the power on scales below λ_mfp. OTVET and M1 over-ionise
  self-shielding absorbers and OTVET cannot cast shadows. This is our
  post-overlap OTVET excess (same pattern, too much power, growing as λ
  grows). Their remedy, a kernel screened by exp(−r/λ) and summed over all
  periodic images, is the baseline to beat here.
- **Gnedin & Abel 2001** (astro-ph/0106278): the OTVET tensor without e^(−τ);
  a separate mean background; a short-range Green function to suppress the
  periodic images distorted the tensor; fronts leak behind obstacles at
  ≈ 1/5 of the normal speed; suggests Monte Carlo or short characteristics
  for better tensors.
- **Gnedin 2014** (CROC; arXiv:1403.4245, App. C): J = J̄ f + L (g − ḡ f), an
  isotropic background with an empirical Lyman-limit absorption plus a local
  part with the thin tensor of the sources in the box: the closest published
  analogue of a direct + diffuse split. **Gnedin 2025** (arXiv:2504.10571):
  CROC (OTVET) and Thesan (M1) Γ–source transfer functions; CROC suppressed
  below 1 cMpc by its approximate tensor; no ray tracing.
- **Finlator, Özel & Davé 2009** (arXiv:0808.3578): tensors from long
  characteristics with e^(−τ) to every source, 26 periodic replicas, updated
  only where J changed by > 5 %; cost ∝ n_grid^3.4–6.1. Their objection to
  the thin tensor is before overlap (fronts elongate towards each other),
  where ours already matches the rays.
- **Petkova & Springel 2009** (arXiv:0812.1801): tree tensor, non-periodic;
  shadows fill in ("the radiation also diffuses around the obstacle").
- **Gnedin & Katz 2026** (arXiv:2603.22400): second-order generalisations of
  M1 and OTVET are unstable as equations; first-order OTVET is stable. So we
  stay with first-order (h-only) closures.
- **Davies & Furlanetto 2016; Davies et al. 2017** (arXiv:1509.07131,
  1611.02711): a mean free path that varies with Γ amplifies the Γ
  fluctuations (semi-analytic) — the reason a single global λ may not be
  enough.
- Tree/attenuated sums *(memory)*: TreeCol (Clark, Glover & Klessen 2012),
  TRAPHIC (Pawlik & Schaye 2008), TREVR (Grond et al. 2019), TreeRay (Wünsch
  et al. 2021), LEBRON (Hopkins et al. 2020: thin tree sum with extinction at
  source and absorber). None was found used as a reionization Eddington
  tensor. Voxel cone tracing (Crassin et al. 2011, graphics) *(memory)*.

*Candidates and ranking: §4 (after the ground truth).*
