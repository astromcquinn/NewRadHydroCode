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

Results in §3.

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

## 3. Where each closure goes wrong (a priori, 256³)

The rays' snapshots from z = 9.6 (x_V = 0.05) to 5.0, each compared with
the tensors the closures give for the rays' own sources and HI field
(`tensor_compare.py`; tables across z from `tensor_summary.py`, in
[`figures/otvet_closure_apriori.txt`](figures/otvet_closure_apriori.txt)).
Closures tested:
- `otvet`: the code's thin minimum-image tensor.
- `m1`: Levermore's D(f) on the rays' true flux, an oracle M1 that is
  better than M1 in practice, since the real M1's flux is its own.
- `split`: candidate 2 on the thin tensor.
- `atten`: candidate 1 at λ_eff, ½λ_eff and 2λ_eff, with λ_eff the
  Γ-weighted photon budget path ΣS0/(4πΣρ).
- `cone`: candidate 3 in numpy (`cone_tensor.py`). nside 2 (48 cones) or
  nside 4 (192 cones); receivers on a grid C = 2 times coarser; near field
  r < r_c = 4 cells; per-cell group-1 τ capped at 1 (`t1`).
- `post_otvet`, `post_m1`: the moment runs' own tensors at the nearest
  snapshot (another radiation field, so for orientation only).

The rays' sampling noise, from the odd and even substeps, is
0.005–0.02 on cells, far below every closure's error.

![a priori errors](figures/otvet_closure_apriori.png)

*Errors against z. Top: the S0-weighted Frobenius error |h − h_true| on
ionised cells, on ionised cells with Γ < ½⟨Γ⟩, and on 8³ blocks (the
tensors behind P_Γ at k ≲ 3 h/Mpc). Bottom: block anisotropy (black: the
truth), and the error power of the radiation force div(S0 h) relative to
the true force, at k = 0.47 and 3.8 h/Mpc.*

**Before overlap (z = 9.6–8.1, x_V ≤ 0.17).**
- OTVET is good in the bulk: cell error 0.06–0.08 at Γ > ⟨Γ⟩/2, block
  error 0.06–0.09, and force error power ≤ 4 % at k ≤ 1 h/Mpc.
- Its errors are concentrated in the 10–20 % of ionised cells where
  Γ < ⟨Γ⟩/2 or the cell is shadowed. There the cell error is 0.27–0.40,
  and the truth is more beamed than OTVET (A = 0.75 against 0.54): the
  light comes only from the region's own sources, while OTVET also adds
  sources in other regions, behind neutral gas.
- Attenuation halves those errors: the cone gives 0.15–0.18 and the
  screened kernel 0.15–0.24.
- The oracle M1 is best in those cells (0.08–0.10; a single beam is its
  exact case) and worst everywhere else.
- The split equals OTVET here (N ≈ N_dir), except at the sources, where
  the recombination radiation makes the truth isotropic (A_true = 0.44
  against OTVET's 0.55 within 1.5 cells). The split lowers the force error
  at k ≥ 7.5 h/Mpc by a third to a half.

**During overlap (z = 7.7–6.4, x_V = 0.25–0.86).** OTVET's large-scale
tensor is **too isotropic**, not too anisotropic.
- Block anisotropy: 0.32 against the truth's 0.42 at z = 7.4.
- Far from sources it is worse: at 24–48 cells from the nearest source,
  A = 0.41 against 0.72 at z = 7.4, and 0.33 against 0.53 at z = 6.6.
- A cell at the edge of an ionised region is lit from the inside of its
  region. OTVET also adds the regions behind the neutral gas, from every
  side.
- The cone tensor follows the truth in every distance bin (0.54/0.72
  against 0.54/0.72 at z = 7.4) and in the block anisotropy (0.415 against
  0.418). Its block error is 0.052–0.062 against OTVET's 0.068–0.100.
- Its force error at k = 0.47 h/Mpc is 0.002–0.013 against OTVET's
  0.004–0.05, and at 0.94 h/Mpc 0.02–0.25 against 0.09–0.40.
- Its low-Γ and shadowed cell errors are 0.11–0.16 against OTVET's
  0.21–0.41.
- One global λ cannot do this. The screened kernel is worse than OTVET on
  blocks and in the force at every k ≥ 1 h/Mpc, because which directions
  should be cut is a property of the geometry, not of the distance.
- **The cone's weakness is small scales.** Its force error at
  k ≥ 3.8 h/Mpc is 2–5× OTVET's at nside 2. Finer angles help: at z = 9.6,
  nside 4 has 0.22 against 0.46 for nside 2 at k = 3.8 (OTVET 0.18). The
  coarse receivers and the near-field hand-over are being tested in
  `cmp_conevar2` (r_c = 8).

**After overlap (z ≤ 6.2): Wu+21's normalisation error, measured.** The
minimum-image thin sum per unit emissivity is
∫_cube d³r/(4πr²) = 0.610 L. Wu+21's L_O is 4π × 0.610 L = 7.66 L, which
confirms the literature subagent's inference. So OTVET's anisotropy is
right only while the photons' path λ_eff ≈ 0.61 L. Beyond that, the
anisotropy is too large by λ_eff/0.61 L:

| z | x_V | λ_eff [cells] | λ_eff/0.61 L | A_OTVET/A_true (8³ blocks) | median N/N_loc |
|---|---|---|---|---|---|
| 6.19 | 0.94 | 155 | 1.01 | 1.07 | 1.06 |
| 5.85 | 0.9997 | 297 | 1.93 | 1.78 | 2.04 |
| 5.67 | 1.0000 | 366 | 2.38 | 2.37 | 2.56 |
| 5.34 | 1.0000 | 504 | 3.23 | 2.86 | 3.61 |
| 4.92 | 1.0000 | 687 | 4.40 | 3.76 | 4.99 |

(At z ≤ 5.34 the ratio falls 10–15 % short of λ_eff/0.61 L: the truth
keeps some anisotropy from local structure. OTVET's block anisotropy stays
at 0.25 while the truth's falls to 0.066. Its block error grows to 0.158;
the split's falls to 0.026, the screened kernel's to 0.023.)

- With the 512³ run's box twice as large and λ_mfp of 100–200 cMpc/h at
  z = 5.4–5.0, the factor is several. That is the 2–9× in P_Γ, growing in
  time.
- The direction is also wrong. On large scales OTVET's radiation force is
  *anti-correlated* with the true force at z ≤ 5.85 (r = −0.6 to −0.9 at
  k ≤ 4 h/Mpc). Its error power is 25–300× the true force's, which is
  small after overlap (κ is small).
- **The split fixes the normalisation with no λ.** Block anisotropy is
  0.138/0.110 against the truth's 0.147/0.108. Block error is 0.051/0.036
  against OTVET's 0.116/0.133. Force error at k ≤ 4 h/Mpc is 5–30×
  smaller, though r(k) stays weak (−0.5 to 0.6): the direction of the
  remaining anisotropy is still the thin one.
- The screened kernel at λ_eff is as good here as the split (0.049/0.033 on
  blocks), but it is worse than OTVET during overlap.
- The cone tensor with periodic images out to 8 boxes counts 3× the true
  light (ΣN_dir/ΣS0 = 2.7–3.1) and is too isotropic (A = 0.03 against
  0.11). The HI column does not stop that light; the light horizon
  (ĉ t since overlap ≈ 2 boxes at ĉ = 0.2c) and redshifting do.
- **An extra screening exp(−s/λ_x) on the cones** (`cmp_lamx`; cone with
  t1, nside 2):

  | z | λ_x | ΣN_dir/ΣS0 | cell err | block err | block A (truth) | force err, k = 0.47 / 0.94 / 3.8 |
  |---|---|---|---|---|---|---|
  | 7.40 | ∞ | 1.01 | 0.095 | 0.062 | 0.415 (0.418) | 0.002 / 0.063 / 1.46 |
  | 7.40 | λ_eff | 0.55 | 0.136 | 0.099 | 0.397 | 0.017 / 0.30 / 1.89 |
  | 7.40 | 4λ_eff | 0.84 | 0.101 | 0.067 | 0.409 | 0.003 / 0.094 / 1.42 |
  | 6.19 | ∞ | 0.87 | 0.077 | 0.053 | 0.255 (0.258) | 0.039 / 0.53 / 4.98 |
  | 6.19 | λ_eff | 0.52 | 0.090 | 0.064 | 0.268 | 0.077 / 0.79 / 7.18 |
  | 5.85 | ∞ | 2.67 | 0.097 | 0.086 | 0.048 (0.147) | 57 / 27 / 56 |
  | 5.85 | λ_eff | 0.73 | 0.064 | 0.046 | 0.120 | 10 / 5.0 / 24 |
  | 5.85 | 2λ_eff | 1.16 | 0.072 | 0.058 | 0.088 | 22 / 9.7 / 25 |
  | 5.85 (OTVET) | | | 0.126 | 0.116 | 0.261 | 185 / 25 / 31 |
  | 5.85 (split) | | | 0.066 | 0.051 | 0.138 | 20 / 4.2 / 5.8 |

  - The cut that helps after overlap hurts during it, where the HI already
    gives the right total light (ΣN_dir/ΣS0 ≈ 1).
  - The self-consistent rule is to choose λ_x so that the cone's direct
    light does not exceed the photons the solution holds (ΣN_dir = ΣN):
    λ_x = ∞ at z ≥ 6.2 and ≈ 1.5–2 λ_eff at z = 5.85.
  - With λ_x ≈ λ_eff the cone matches or beats the split at k ≤ 1 h/Mpc
    and on cells and blocks. It is still 4× worse at k ≥ 3.8 h/Mpc.
  - Combining it with the split makes it worse (force 23 against 10 at
    k = 0.47). Once normalised, the cone's own h = P/tr P is better than
    adding an isotropic remainder.
- The oracle M1 is too isotropic (A = 0.015–0.03 against 0.11): the net
  flux of a nearly isotropic field says little about its anisotropy. Its
  force error is as large as OTVET's.

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

**Ranking** (from §3; for the user's decision before anything goes into
the solver):

1. **Split (2) on the thin tensor.**
   - It is the only candidate that is never worse than OTVET, from
     z = 9.6 to 5.7.
   - After overlap it removes Wu+21's normalisation error with the
     solution's own N and no λ: block error 2–4× lower, force error 5–30×
     lower.
   - Near sources it handles the recombination radiation (force error at
     k ≥ 7.5 h/Mpc ⅓–½ lower).
   - Cost: one float per cell and no extra FFT.
   - It does nothing for the overlap-era geometry.
2. **Cone-traced attenuated direct light (3), with a budget-matched
   cutoff.**
   - It gets the overlap-era geometry right where nothing else does:
     anisotropy against distance and on blocks, beaming far from sources
     and in shadow. Large-scale force error is 5–20× lower than OTVET's,
     and block error a third to a half lower.
   - After overlap the light the HI does not stop must be limited by an
     extra exp(−s/λ_x). λ_x is chosen so that ΣN_dir = ΣN of the solution:
     ∞ during overlap, ≈ λ_eff after it. This is the solution's own photon
     budget, adjusted each hydro step from the last one.
   - With that cutoff it matches the split after overlap at k ≤ 1 h/Mpc.
     Not combined with the split, which makes it worse.
   - Cost: a C++ cone marcher. At 512³, 128³–256³ receivers and 48–192
     cones take roughly 10–75 s per hydro step, against the present 28 s.
     Memory ≲ 1.5 GB.
   - Weakness: small scales (k ≥ 4 h/Mpc), at present worse than OTVET.
     Finer cones help; the near-field radius is being tested.
3. **Screened periodic kernel (1), the baseline.**
   - After overlap it matches the split.
   - During overlap it is worse than OTVET at every k ≥ 1 h/Mpc.
   - Implemented for the comparison, with λ from the solution's own budget.

**Proposed implementation (step 3, after the user's go-ahead).**

One new flag, `otvet_closure`, whose default `thin` is today's tensor,
byte-identical. All options are deterministic and independent of the
thread count.

- **`split`**
  - The FFT path of `compute_otvet_tensor` already forms P. Keep tr P in
    absolute units: N_dir,g = q_g tr P / (4π), with q_g the moment
    injection per unit Luv of group g. That is one more float per cell; the
    per-group factor is a scalar.
  - `transport_group` then forms
    P_g = N_g I/3 + min(N_g, N_dir,g)(h − I/3) in place of N_g h.
  - With `min(N, N_dir) ≤ N`, the eigenvalues of P stay in [0, N], so the
    GLF fluxes keep their wave speed ĉ.
- **`screen`** (the baseline)
  - The kernel is e^(−r/λ) d dᵀ/r⁴ summed over the 27 images, plus the
    isotropic continuum beyond them.
  - λ = ĉ ΣN / ΣṄ_inj, from the moment solution's own budget at the
    previous hydro step.
  - The kernel is rebuilt when λ changes (6 more FFTs).
- **`cone`**
  - A new `src/otvet_cone.cpp`, reusing `healpix` from the ray tracer.
  - Mip-maps of the source emissivity and of the HI column. Receivers on a
    grid `otvet_cone_coarse` = 2 times coarser, `otvet_cone_nside` = 2 or
    4, near field r < r_c from the FFT with the thin kernel.
  - Per-cell group-1 τ capped at `otvet_cone_taucap` = 1.
  - λ_x chosen so that ΣN_dir = ΣN of the moment solution. It is updated
    from the previous hydro step: the cone records its direct light in a
    few distance bins, so the matching λ_x is found without a second march.
  - OpenMP over receivers, each receiver independent.
  - If the blend test (§3) holds up, the cone's large scales and the
    split's small scales are combined (`cone_blend`).

**Checks.**
- Each C++ tensor is computed on the rays' z = 7.4, 6.2 and 5.67 snapshots
  and compared with the numpy prototypes and with h_true (`tensor_compare`
  reads them through `--moment`).
- Then full 256³ runs of each option with `tensor_dump`, about 45 min each
  for the moment methods, against `otvcl256` (rays). They are validated
  with `gamma_power.py` (P_Γ, r(k), PDF, spread), `reion_power.py` and
  `rt_methods_runs.py` (history, morphology, T).
- Then the winner at 512³ against the existing 512³ rays.

Not recommended, on the evidence of §3:
- **8, mixing with M1:** the oracle M1 is the best closure only in shadowed
  pre-overlap cells, which the cone also fixes, and it is the worst on
  large scales in every regime.
- **7, receiver λ(x):** the cone gives a spatially varying attenuation
  directly, along each direction.
- **4 and 11, long characteristics and a low-resolution ray tracer:** the
  same physics as the cone at 10–100× the cost.
- **6, per group:** second order for Γ_HI; revisit for He II.
- **9, retarded sources:** small.
- **10, fitted closure:** a diagnostic, not a closure.
