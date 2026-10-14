# Literature review: OTVET Eddington tensors, M1, and Γ_HI fluctuations after overlap

Compiled 2026-10-12 for the OTVET-closure investigation (40 h⁻¹ Mpc, 512³; OTVET P_Γ is 2–9× the rays' after overlap
with r ≥ 0.95, the Γ PDF has a floor near 0.62⟨Γ⟩, and proximity peaks are too low).

**Method and caveats.** The two local PDFs could not be read: the Read tool needs `pdftoppm` and the node has no
pdftotext, pypdf or fitz. All papers were read as HTML (ar5iv or arxiv.org/html) through WebFetch. WebFetch runs the page
through a summarising model, which caps verbatim quotes at about 125 characters. So the quotes below are short and are
as that tool returned them. **Spot-check them against the PDF before citing them in a paper.** Where a page was too long,
some sections were not read; each such gap is noted below. Each item is tagged:

- **[READ]**: taken from the fetched text.
- **[INFER]**: my own derivation or reading.
- **[RECALL]**: from memory, not verified this session.

---

## 1. Wu, McQuinn & Eisenstein 2021, JCAP 02 (2021) 042, arXiv:2009.07278

"On the accuracy of common moment-based radiative transfer methods for simulating reionization"

Authors are from the arXiv abs page. The first ar5iv pass misattributed it to Eisenstein alone; that was an extraction error.

**What they do [READ].** The paper is analytic or quasi-analytic. It does not post-process any specific code's run. There
are three problems:

- an expanding H II region (§3);
- a spherical self-shielding absorber in a uniform background, with n_H ∝ r^−m, Γ(∞) = 5×10⁻¹³ s⁻¹ (§4);
- linear perturbation theory for Γ/intensity fluctuations in a static post-overlap IGM at z = 5.5, with λ_mfp = 30 and
  60 cMpc and halo-model sources (§5).

Closures: OTVET is h = P/TrP with P_mn = ∫ j(x′) G_mn(x−x′) d³x′ and G_mn = x_m x_n/|x|⁴, that is τ = 0 (Eqs. 2.10,
2.13–2.14). M1 is treated generically, citing several codes. Ray tracing is discussed through its angular resolution.

**How the OTVET tensor is computed: periodic images [READ].**

- §2.1.1: "Fourier transforming G_mn leads to a divergent zeroth mode … implying that the Eddington tensor is isotropic
  everywhere in a periodic volume."
- Hence "implementations of OTVET only include image sources out to half a box size away". App. B.2 says the same for
  the CROC simulations of Gnedin 2014, ref. [3].
- The tensor isotropises beyond r ~ 1/√(n_* λ_mfp) (Eqs. 2.15–2.16). OTVET effectively uses λ_mfp → L_box/2, "an
  incorrect scale if λ_mfp ≠ L_box/2". For example, a 20 cMpc box "implies an effective λ_mfp = 10 Mpc".
- No attenuation, background term or smoothing appears in the tensor.

**Linear fluctuation results [READ].**

| | δ_E (Fourier) | Eq. |
|---|---|---|
| Exact | S(k) δ_j, with S(k) = arctan(kλ)/(kλ), aκ̄ = 1/λ | 5.6–5.7 |
| M1 (isotropic at linear order) | (δ_j − b_κ,δ δ)/(1 + b_κ,Γ + k²λ²/3) | 5.14 |
| OTVET | [(1 + 2π²kλ²/(3L_O)) δ_j − b_κ,δ δ]/(1 + b_κ,Γ + k²λ²/3) | 5.15 |

M1 is correct only for k ≲ 1/λ. For OTVET, L_O = ∫_V r⁻² d³x over the simulated volume, with periodic wrap. Eq. 5.15 is
valid only for k ≫ 2π/L.

Quantitative statements:

- "For box sizes smaller than twice the mean free path, the OTVET Eddington tensor is more anisotropic than the exact
  Eddington tensor." The result is an overestimate of sub-box fluctuations.
- With λ = 60 Mpc: L = 30 Mpc overestimates the 50 h⁻¹ Mpc-scale variance by 3–4×; L = 1000 Mpc underestimates it by
  20–30%.
- §5.3 and §7: simulations "likely overpredict fluctuations in the ionizing background by factors of 2−4". Box sizes near
  2λ_mfp "likely fair much better".
- M1: "tends to isotropic after reionization" and "leads to significantly suppressed power on scales smaller than the
  photon mean free path". The 50 h⁻¹ Mpc variance is low by 20–30%, and by 60–70% if λ is overestimated by 2×.
- Ray tracing: "a high number of angular directions need to be followed". Closing at order ℓ is roughly ℓ² directions,
  which converges only for k ≲ ℓ/λ. Capping at 64 rays per cell is roughly ℓ = 8 (§5.2.3).

**Absorbers and mean free path [READ].**

- OTVET and M1 "over-ionize self-shielding absorbers". The recombination rate is about 1.3× (OTVET) and 1.4× (M1) for
  r⁻².
- Matching Γ needs 30–40% more emissivity. At fixed emissivity Γ is about 0.4× exact, and the effective λ_mfp is biased
  by about 30% (fixed Γ) to 2.5× (fixed emissivity).
- OTVET/FLD cannot "cast shadows behind dense absorbers".
- A flux limiter "is likely required to ensure |F| ≤ E in the proximity zone".
- From the reduced-c dependence of x_HI, simulated λ_mfp is "several times larger" than observed (§6.1, abstract).

**Remedies named [READ].** Long-characteristics or exact tensors with 26 periodic replicas (Finlator et al.), which are
costly. Higher-ℓ closures (App. C, only partly read). Many-direction ray tracing. No concrete OTVET fix (screening,
background term) is proposed in the text read. App. B and §5.2.3–6 were only partly read because of truncation.

**[INFER] Implication for us.**

- At kλ ≫ 1, Eq. 5.15 gives δ_E^OTVET → 2π²δ_j/(kL_O), and the exact result is → πδ_j/(2kλ).
- So OTVET has the **correct 1/k shape and the wrong amplitude**: ratio = 4πλ/L_O, k-independent. That fits our
  "r ≥ 0.95 but 2–9× power".
- M1 instead goes as 1/k². That explains its suppression and decorrelation.
- For our minimum-image cube, L_O = ∫_cube r⁻² d³x ≈ 7.66 L (computed numerically), so
  P_OTVET/P_exact ≈ (1.64 λ/L)². This equals 1 at λ ≈ 0.61 L. A 2–9× excess implies λ_mfp(sim) ≈ 0.86–1.8 L, i.e.
  35–75 h⁻¹ cMpc.
- **Measure λ_mfp in the run to test this.**
- The natural fix is a screened kernel, G_mn = x_m x_n e^(−r/λ)/r⁴, summed over all periodic images. Screening makes
  the sum converge, so it can be done by FFT on a padded or periodic grid. At kλ ≫ 1 this replaces L_O by
  ∫e^(−r/λ)r⁻²d³x = 4πλ, so the small-scale amplitude is right by construction.
- Equivalently, keep minimum image and add an isotropic remainder (δ_mn/3)·j̄·4πλ·e^(−L/2λ) (approximately) to P. This is
  a "background" term in the Gnedin & Abel sense.
- A local-κ variant, e^(−τ) with τ along the ray (or a TreeCol-like column), would also give shadows. It costs more.

## 2. Gnedin 2025, arXiv:2504.10571, "Do we know how to model reionization?"

**[READ]**

- Compares CROC (ART; OTVET; three 40 h⁻¹ cMpc boxes C40A/C/F) with Thesan-1 (AREPO-RT; M1; about 65 h⁻¹ cMpc). Both use
  "Variable Eddington Tensor" moment methods.
- CROC's OTVET solves a single Newtonian-limit equation (Eq. 4), (a/c)∂E/∂t = ∂_i[(κ+ε)⁻¹ ∂_j(E h^ij)] − κE + L. The ε
  is there "to avoid occasional numerical overflow".
- CROC uses ε = 10⁻³/Δx, which can amplify noise in low-opacity gas. ε = 0.01/Δx removes the excess small-scale power
  but breaks photon conservation.
- The comparison statistic is T_ΓL(k) = (P_Γ/P_L)^½ at matched volume-weighted x_HI. The expected form is
  arctan(kλ)/(kλ) (∝ 1/k at small scales), with λ fitted on large scales.
- Large-scale slopes agree, but amplitudes differ.
- CROC: T_ΓL "exhibits a suppression at scales below about 1 cMpc", which is "due to the highly approximate nature of
  the optically thin Eddington tensor calculation in CROC". ALTAIR runs with an accurate OT tensor (same scheme) "do not
  show such a suppression". ALTAIR M1 and OTVET both match the analytic expectation.
- Thesan: "exhibits much larger power on small scales, and the origin of that power is not clear". A stochastic SFR and
  rebinning noise are considered and judged insufficient.
- CROC has a separate "background" radiation component, for radiation entering or leaving the box and for quasars. It was
  excluded from the comparison "so that the comparison … is as direct as possible".
- Conclusion: "the answer to the paper title is 'not sure'".
- No ray tracing. Wu+21 is not cited in the text read. There is no mean-free-path analysis beyond λ as a fit parameter,
  and no explicit RT recommendations.

**[INFER]**

- CROC's suppression comes from an approximate (coarse or softened) tensor evaluation. Our FFT tensor is exact for
  minimum-image sources, so we see the opposite problem, the Wu+21 box-size anisotropy excess.
- The ALTAIR OTVET runs being "accurate" suggests that a correctly evaluated OT tensor in a box near 2λ behaves.
- The CROC-vs-Thesan comparison was made with the background component off.

## 3. Gnedin & Abel 2001 (New Astron. 6, 437), astro-ph/0106278; Gnedin 2014 (ApJ 793, 29), arXiv:1403.4245

### Gnedin & Abel 2001 [READ]

**Tensor.**

- h = P/TrP (Eq. 22), with P^ij = ∫d³x₁ ρ_*(x₁)(x−x₁)^i(x−x₁)^j/|x−x₁|⁴ (Eq. 24). This is Eq. 23 with e^(−τ) dropped,
  so there is **no attenuation**.
- It can be computed "in O(N) operations by standard algorithms used to compute high-resolution gravitational forces such
  as P3M, Tree, or Adaptive Mesh Refinement".
- The tensor is built from the full ρ_*, with no mean subtracted.

**Background split.**

- The mean J̄_ν is evolved separately (Eq. 4), with k̄_ν = ⟨k_ν J_ν⟩/J̄_ν.
- The local field is J = J̄ E₁ + g D₂ (Eq. 14). D₂ is sourced by ρ_* − ρ̄_* Λ_ν (Eq. 17).
- Λ_ν = exp(−⟨κ_ν⟩L), with L = (1/4π)∫d³x/r². It goes to 0 when λ ≪ box and to 1 when λ ≫ box. "Exact choice matters
  little".
- The authors say a first splitting (Eq. 16) does not work in practice.

**Periodic images.**

- A second tensor eigenvector "arises due to periodic boundary conditions as a contribution from periodic images of the
  source".
- They tried a short-range Green function (falling faster than 1/r²) to suppress images. It "introduce[s] deformations
  into the Eddington tensor" larger than the original error, so there was no improvement.

**Shadows.** Behind a dense wedge (64³), the front moves at about 1/5 of the normal I-front speed. In a 2× box this
drops by about 0.6, which they attribute to images. With the exact tensor there is no propagation into the shadow.

**Improvements.**

- "One may envision to employ a Monte Carlo or short characteristic methods to supply more accurate tensors."
- Exact tensors at selected points seem unlikely "with acceptable computational cost".
- Their scheme is the reduced-c, Newtonian-limit diffusion form (Eqs. 13, 20–21), with no flux limiter mentioned.

### Gnedin 2014, CROC I [READ, App. C]

**Two-field ansatz.**

- J_ν = J̄_ν f_ν + L_ν(g_ν − ḡ_ν f_ν) (C2).
- g, the local sources, obeys Dg/dl = −kg + ρ_S (C4). Its formal solution includes e^(−τ) (C5). In the code it is
  solved through the angle-averaged OTVET moment equation (C10). The appendix shows no ray integration of C5, so the
  attenuation is formal only.
- The tensor for G is "the optically thin Eddington tensor from all sources inside a periodic box". How the sum is
  computed (images or not) is not stated in the text read. Wu+21 App. B.2 says images reach out to half a box.
- f, the background with local absorption, obeys C7, stabilised as C9 with q(f̄) = 2/f̄ − 1. Its tensor is isotropic:
  "the Eddington tensor for the background radiation is taken to be isotropic".
- Background Lyman-limit absorption uses the Songaila & Cowie 2010 fit and is not computed from the box, because
  "periodic boundary conditions will extend that box over the whole universe". Quasars enter only the background.
- No fluctuation or mean-free-path results; no "OTVET with attenuation" tensor.

**[INFER]** CROC's design is the closest published analogue of "minimum-image anisotropic local part plus isotropic
background". Our implementation has only the first piece.

## 4. Petkova & Springel 2009 (MNRAS 396, 1383), arXiv:0812.1801 [READ]

- OT tensor (h = P/TrP) from a gravity-like tree. Nodes store luminosity and a luminosity-weighted centre, with the
  gravity opening criterion ("typical accuracy of ∼1%"). It is computed during the gravity walk at "negligible" cost.
- It is "computed non-periodically". No softening of the 1/r² kernel and no attenuation are mentioned.
- Implicit backward-Euler diffusion form with an SPH discretisation and CG with a Jacobi preconditioner. Optional flux
  limiter. No reduced c.

**Problems.**

- Anti-diffusion is possible; the fix is to blend in an isotropic part (Eq. 53).
- Errors at points "equidistant from two sources".
- Two-source H II regions are elongated along the joining axis; the spherical average is within 20%.
- Shadows: "the radiation also diffuses around the obstacle". Setting the tensor to zero in the shadow gives a sharp
  shadow.
- "the long-term evolution of sharp geometric shadows is clearly not followed accurately".
- The optically thin tensor errors "should be quite moderate or negligible" when a bright nearby source dominates.

## 5. Finlator, Özel & Davé 2009 (MNRAS 393, 1090), arXiv:0808.3578 [READ, last ~15% not read]

**Tensor.**

- Long characteristics from every cell to every source cell, with attenuation: N(n̂) = Σ_i (η_i/χ_i)e^(−τ_i)δ(n̂−n̂_i)
  (Eq. 14).
- Rays stop at τ_max = 6.
- Periodicity is handled by brute-force replicas: "nd = 1 corresponds to positioning 26 replica volumes". A hybrid
  switches to nd = 2 when x̄_HI < 0.5.
- The tensor is smoothed with a 27-cell tophat. There is no source merging.

**Cost and update.**

- t_LC ∝ n_grid³ n_S n_l, empirically ∝ n_grid^3.4–6.1; the optically thin limit would be n_grid⁷.
- Recomputed per cell only when J changes by more than f_J = 0.05 (t ∝ f_J^−0.5).
- In the cosmological test, a full update happens when J changes more than 2× in at least 20% of the volume.

**Why the optically thin tensor is inadequate.**

- Regime: multiple sources in an optically thick medium **before overlap**. The test has 5×10⁴⁸ and 5×10⁴⁷ s⁻¹ sources
  2.2 kpc apart at n_H = 10⁻³.
- "the optically thin approximation causes both H II regions to elongate along the axis that connects the two sources",
  which gives early overlap.

**Other.**

- Time-dependent moments with the full c; reduced or infinite c is criticised.
- Shadows are 4–5 dex fainter at 1 Myr and fill to about 1% by 15 Myr.
- No diffuse or background component is discussed.

**[INFER]** Their criticism is pre-overlap, and pre-overlap is exactly where our OTVET already matches rays. Our failure
is post-overlap, which is Wu+21's regime.

## 6. Attenuated, tree or hybrid tensors and closures (brief)

**[READ]**

- Gnedin & Katz 2026, arXiv:2603.22400 (OJAp 9), "Second Order Closures … Some Are Unstable":
  - Direct second-order generalisations of M1 and OTVET are unstable at the level of "the closure equations themselves".
  - OTVET generalised to any order above first is unstable.
  - Any local second-order closure depending only on E and F (not P) is physically unstable.
  - Abstract only; no new closure proposed.
- Davis, Stone & Jiang 2012, arXiv:1201.2223: a VET from a short-characteristics formal solution in Athena captures
  shadows (search snippet only).
- Aubert & Teyssier 2008 (ATON), arXiv:0709.1544: M1 "local Eddington tensor" for reionization (title only).
- Cain & D'Aloisio 2024 (FlexRT), arXiv:2409.04521: ray tracing. The intro names M1 and OTVET as the two prominent
  closures (search snippet).
- Searches found no 2022–26 paper using an attenuated or tree-based source sum for a reionization Eddington tensor, and no
  direct+diffuse moment decomposition beyond CROC's background split.

**[RECALL, unverified]**

- Rosdahl et al. 2013, RAMSES-RT (arXiv:1304.7126): M1.
- Kannan et al. 2019, AREPO-RT (arXiv:1804.01987): M1.
- Hopkins et al. 2020 (FIRE RT methods, arXiv:1811.12462): compares LEBRON with M1. The abs page confirms LEBRON and M1.
  - LEBRON is an optically thin gravity-tree sum with local extinction applied at source and absorber, a cheap "TreeCol-
    lite" attenuation.
  - GIZMO also has OTVET with a gravity-tree tensor.
- TreeCol (Clark, Glover & Klessen 2012), TRAPHIC (Pawlik & Schaye 2008), TREVR (Grond et al. 2019) and TreeRay (Wünsch
  et al. 2021) give tree-based attenuated columns or reverse ray tracing. None is used as an Eddington-tensor source in
  reionization as far as I know.
- Hayes & Norman 2003: VTEF with a short-characteristics tensor.

## 7. M1/OTVET and Γ fluctuations or mean free path after overlap

**[READ]**

- Wu+21 (§1 above) is the main reference. M1 suppresses small-scale P_Γ. OTVET's error is set by L vs 2λ.
- Simulated λ_mfp is several times too long because of reduced c.
- Gnedin 2025 T_ΓL comparison: CROC is suppressed below about 1 cMpc; Thesan (M1) has an unexplained small-scale excess
  or plateau; neither follows 1/k.
- Davies & Furlanetto 2016 (arXiv:1509.07131) and Davies et al. 2017 (arXiv:1611.02711): fluctuating λ_mfp that couples
  to Γ amplifies Γ fluctuations. These are semi-analytic, not moment RT (search snippets).
- Thesan reportedly reproduces the rapid λ_mfp evolution (search snippet; Garaldi et al. 2022 [RECALL]).

**[RECALL]** Keating et al. 2020 and Kulkarni et al. 2019 (ATON/M1, Sherwood) and CoDa-III (Lewis et al. 2022,
RAMSES-CUDATON M1) model τ_eff fluctuations mostly through late neutral islands, not post-overlap Γ structure.

---

## Bottom line for the OTVET closure work [INFER]

1. Our post-overlap excess is exactly the Wu+21 OTVET box-size effect. Minimum-image OTVET normalises anisotropy by
   L_O ≈ 7.66 L instead of 4πλ_mfp. That gives the right 1/k shape, hence the high correlation, with power biased by
   (1.64 λ/L)².
2. **Cheapest principled fix:** an exponentially screened OT kernel e^(−r/λ_eff)x xᵀ/r⁴ with all periodic images. Equivalently,
   minimum image plus an analytic isotropic remainder (a CROC-like background with an isotropic tensor). Take λ_eff from
   the measured mean free path, or from ⟨κ⟩ in ionised gas.
   - Testable prediction: P_Γ should then match the rays at kλ ≫ 1. Check this by measuring λ_mfp in the run.
   - Before overlap, λ_eff ≫ bubble size, so behaviour should be largely unchanged.
3. **The Γ floor and the missing low tail** are consistent with OTVET's documented lack of shadows and over-ionisation of
   self-shielding absorbers (Gnedin & Abel; Petkova & Springel; Wu+21 §4).
   - Fixing that needs local attenuation in the tensor: e^(−τ) along the ray, as in Finlator et al. with long
     characteristics, or a TreeCol/LEBRON-like column.
   - A screened global λ will not fix it.
4. **Low proximity peaks** are not a tensor problem. The OT tensor is exact for an isolated source (Wu+21; Gnedin & Abel).
   Suspect explicit-GLF numerical diffusion, reduced-c, or source-cell injection or flux limiting instead (Wu+21 note that a
   flux limiter is needed in proximity zones).
5. **Avoid second-order closures:** Gnedin & Katz 2026 show the direct generalisations are unstable.
