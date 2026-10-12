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

## 4. Candidate closures

Costs are for the 512³ run, where the present tensor (13 single-precision
FFTs of 512³ per hydro step) costs ≈ 28 s per hydro step (0.43 h of 4.8 h)
and 6 floats per cell (3.2 GB of 146 GB). All are first-order closures
(a tensor h, or P(N) = D N): Gnedin & Katz 2026 show that second-order
generalisations are unstable.

**What sets the direction of the radiation.**
- *Before overlap:* inside an H II region the gas is thin; the radiation
  is the direct light of the region's own sources, at ≤ ĉ t from them. The
  thin tensor is right except for sources in other regions (blocked by the
  neutral gas in between) and near the fronts.
- *During overlap:* light from several regions meets; neutral islands and
  dense filaments cast shadows, so the islands should be lit from the open
  sides only.
- *After overlap:* λ_mfp (≈ 100 cMpc/h at z ≈ 5.5 in these runs, several
  boxes) ≫ the source spacing. The radiation in a void is mostly light from
  many distant sources and periodic images, nearly isotropic, plus the
  direct light of the nearby sources; self-shielded clumps cast small
  shadows. The minimum-image thin tensor keeps only the nearby sources up to
  L/2, so it is too anisotropic away from them (Wu+21).

**Candidates** (physics; cost; memory; risks):

1. **Screened periodic kernel (baseline, Wu+21).** K = Σ_images e^(−r/λ)
   d dᵀ/r⁴ (27 images plus an isotropic continuum for the rest), with λ the
   photons' mean path from the moment solution's own budget,
   λ = ĉ N_flight/Ṅ_inj (absorption, redshifting and the light horizon
   together). Same FFTs; a kernel rebuilt when λ changes; no memory. Gets
   the box-size normalisation right; no shadows, one λ everywhere.
2. **Direct + diffuse split, from the solution's own N** (new; CROC's
   background split made local and self-consistent).
   P_g = N_g I/3 + min(N_g, N_dir,g)(h_dir − I/3): the direct light of the
   sources (thin sum in absolute units, the trace the FFT already computes)
   carries h_dir, and whatever the solution holds beyond it — wrapped and
   old photons, light from beyond λ, recombination radiation — is
   isotropic. Before overlap N ≈ N_dir in the H II regions, so it stays
   OTVET; after overlap the voids become diffuse by the amount the solution
   itself says. Per group for free (harder groups are more diffuse). Cost:
   one float per cell and a min per cell; P is piecewise linear in N with
   realizable slopes (I/3 and h_dir), so the GLF stays stable. Risks: the
   diffuse part has no direction (a streaming diffuse field becomes P1-like);
   unattenuated N_dir behind absorbers falls back to the thin tensor.
3. **Cone-traced attenuated direct light** (new; TreeCol / voxel cone
   tracing). For each receiver (a coarse grid) and N_pix HEALPix cones,
   march through mip-maps of the emissivity and the HI column with steps
   ∝ the cone width, accumulating e^(−τ) with the group-weighted
   transmission, periodically out to the horizon; the near field from the
   fine thin FFT. Gives the shadows of islands and large absorbers, a
   spatially varying λ, the periodic images. Combined with 2 for the rest.
   Cost (C++ estimate): 128³ receivers × 48 cones × ~20 samples × 2 fetches
   ≈ 2·10⁹ fetches, seconds per update; 192 cones on 256³ receivers ≈ 30×
   more. Memory: the pyramids and a coarse P, ≲ 1 GB. Risks: e^(−⟨τ⟩) of a
   cone average is not ⟨e^(−τ)⟩ for clumpy absorbers (a per-cell τ cap as a
   covering-fraction proxy); shadows narrower than a cone are lost.
4. **Coarse long characteristics from the sources** (Finlator+09 on a
   coarse grid, fine near field). Exact e^(−τ) per source; cost ∝ N_src
   × (L/2Δ)³ ≈ 10⁹–10¹⁰ segments per update at 512³ (minutes), minimum image
   only (images through 2).
5. **Flux-consistent hybrid.** The direct part as in 1–3 with its thin flux
   F_dir; the remainder (N − N_dir, F − F_dir) closed by M1, so a streaming
   diffuse field keeps its direction. +3 FFT components and 3 floats per
   cell; M1's crossing artefacts confined to the remainder.
6. **Per-group screened tensors** (λ_g): 5× the FFTs and +13 GB; matters
   more for He II and the temperature than for Γ_HI.
7. **Receiver-dependent λ(x):** screened kernels at a few λ, interpolated
   per cell by a λ(x) from the solution's smoothed κ (the Γ–λ coupling of
   Davies & Furlanetto on the receiving side). n_λ × the FFT cost.
8. **Regime mixing OTVET/M1** by f or the local optical depth. Cheap, but
   M1 fails where the radiation of many sources crosses, i.e. in the regime
   we would switch to it.
9. **Retarded sources** L(t − r/ĉ): store the luminosities of the last few
   hydro steps. Delays across the box are tens of Myr at ĉ = 0.2c; a small
   effect except near new sources.
10. **Fitted closure** h = a h_dir + b I/3 + c f̂f̂ with a, b, c functions of
    N/N_dir, f and x_HII, fitted to the rays: a diagnostic for the family of
    2 and 5 (how much a local mix can gain), not a closure on its own.
11. **The ray tracer itself at low resolution** (steady state, coarse
    grid, every hydro step): the most faithful, at roughly the cost of 4.

*Ranking: after the ground truth (§3).*
