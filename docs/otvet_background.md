# OTVET with an isotropic background (Gnedin & Abel 2001; CROC, Gnedin 2014)

*Branch `otvet-background`, started 2026-11-06. Step 1 (the design,
§1–7) and the user's decisions (§7); step 2, `otvet_background = ga01`, is
implemented (§8); step 3, the 256³ validation, follows (§9).*

**The question.** Our `rt_method = otvet` gives the whole radiation field
the optically thin tensor h. Both published OTVET schemes give h only to a
local part and carry the rest as an isotropic background. The paper should
be able to compare against OTVET as those papers use it. This note maps both
schemes onto our solver, lists the candidate formulations and their costs,
and recommends one.

## 1. Sources

Both papers were downloaded from arXiv on 2026-11-06, as the LaTeX source
and the PDF (arXiv was reachable from klone that day, unlike in October).
The equation numbers below were counted from the LaTeX source.

- **Gnedin & Abel 2001**, New Astron. 6, 437, astro-ph/0106278,
  "Multi-dimensional cosmological radiative transfer with a Variable
  Eddington Tensor formalism". Read: §2 (eqs. 1–10), §3 (eqs. 11–24).
- **Gnedin 2014**, ApJ 793, 29, arXiv:1403.4245, "Cosmic Reionization On
  Computers I". Read: §2 (the paragraphs on radiative transfer and the
  background), App. A–C (eqs. C1–C12).

Quotes are verbatim from the sources. The notes in
[`otvet_closure_literature.md`](otvet_closure_literature.md) §3 agree with
the full text. One correction to them: G&A's Λ uses the volume mean of the
opacity, ⟨κ_ν⟩, and not the radiation-weighted k̄ of their eq. 4 (§2.1).

## 2. The two schemes as published

### 2.1 Gnedin & Abel 2001

**Splitting off the mean.**
- Eq. 4 evolves the mean intensity J̄_ν (the "background") with the box
  mean of the sources S̄_ν, and with k̄_ν ≡ ⟨k_ν J_ν⟩/J̄_ν, the
  radiation-weighted absorption. They note "k̄_ν is not a space average of
  k_ν since it is weighted by the local value of the specific intensity".
- Eqs. 5–10 write the rest relative to the mean. Eq. 5 defines J_ν = f_ν J̄_ν,
  with ⟨f_ν⟩ = 1. Eq. 7 is the transfer equation of f in the Newtonian
  limit. Its absorption is κ̂ = κ − κ̄ + S̄/J̄ (eq. 8). Eqs. 9–10 are its two
  moments, with ⟨E_ν⟩ = 1.
- Eq. 13 combines the two moments into one RSLA diffusion equation,
  2(a/ĉ)∂E/∂t = ∂_j(κ̂⁻¹ ∂_i(E h^ij)) − κ̂E + ψ. Eqs. 18–19 solve it
  semi-explicitly, with the absorption implicit.

**The two-field ansatz** (eq. 14). In the optically thin limit, with every
source sharing one spectrum S_ν = g_ν ρ_*, the field is ⟨J_ν⟩_Ω = J̄_ν E₁ +
g_ν D₂. Here E₁ = 1, and D₂ is the 1/r² sum of ρ_* − ρ̄_* (eq. 15). The
split serves to factor the frequency dependence: "this procedure is only
required in order to reduce the number of frequency bins in the
calculation".

**The abandoned split (eq. 16).** E₁ gets the isotropic tensor and the full
mean source ψ̄. D₂ gets the thin tensor h and the source ρ_* − ρ̄_*. Why it
fails, in their words: "the scheme (16) does not work, because Ẽ₂ =
g_ν D̃₂/J̄_ν can then be negative and positive, and numerical truncation
errors will prevent the exact cancellation for E_ν = Ẽ₁ + Ẽ₂ outside the
ionization front, where E_ν should be exactly zero."

What goes wrong, in more detail:
- Before overlap, the mean source ψ̄ feeds E₁ everywhere, including the
  neutral gas beyond the fronts. The −ρ̄_* of D₂ has to remove it there
  again.
- The two fields are solved with different tensors, so the cancellation
  is not exact.
- The local field changes sign. Its source is negative wherever there is
  no star, which is almost every cell.

**The adopted split (eq. 17).**
- E₁: tensor δ^ij/3, absorption κ̂, source ψ̄Λ_ν.
- D₂: tensor h (eq. 22 with the thin P of eq. 24, built from the full ρ_*),
  absorption κ, source ρ_* − ρ̄_* Λ_ν.
- Λ_ν "is only required to go to 0 when the mean free path is much shorter
  than the size of the periodic computational box, and to 1 when the mean
  free path is much larger than the box size". They "adopt Λ_ν =
  exp(−⟨κ_ν⟩𝓛)", with 𝓛 = (1/4π)∫d³x₁/(x − x₁)² over the box ("apparently
  independent of x").
- On the choice: "other choices of Λ_ν give different solutions when the
  mean free path is close to the box size, but converge to the same
  solution when the mean free path is much shorter or much larger".
- Before overlap ⟨κ⟩ is large, so Λ = 0. E₁ then has no source, D₂ has the
  positive ρ_*, and there is nothing to cancel.
- D₂'s source is still negative off the sources once Λ > 0. The paper does
  not discuss that further.

**For a cube with minimum-image separations**, 𝓛 = ⟨1/(2 max_i |n_i|)⟩
over directions = **0.611 L**. This is a Monte Carlo over 4·10⁶ directions,
consistent with the 7.66 L/4π = 0.610 L of `otvet_closure_literature.md`
§1.

**Absolute variables.** Multiply E₁ by J̄. Using eq. 4 without the Hubble
term, (a/c)∂J̄/∂t = −κ̄J̄ + (a/c)S̄, and with eq. 8:

```
(a/c) ∂(J̄E₁)/∂t = J̄[−∂_iF₁^i − (κ − κ̄ + ψ̄)E₁ + ψ̄Λ] + E₁[−κ̄J̄ + (a/c)S̄]
                 = −J̄ ∂_iF₁^i − κ (J̄E₁) + (a/c) Λ S̄
```

- The background in absolute units, N_B = J̄E₁, obeys the ordinary transfer
  equation, with the local κ and a uniform source ΛS̄. The local part,
  N_L = g D₂, has the same κ and the source S − ΛS̄.
- κ̂ only normalises the field to J̄. J̄ itself is the box mean of
  N_B + N_L.
- **So G&A conserve photons in the box and need no outside input except
  Λ.** The photons are the box's own: emitted S, absorbed κ, the periodic
  images included.

### 2.2 CROC (Gnedin 2014, App. C)

**The ansatz** (eq. C2): J_ν = J̄_ν f_ν + L_ν(g_ν − ḡ_ν f_ν).
- J̄ obeys the cosmological equation that follows C2 (unnumbered), with
  k̄ ≡ ⟨kJ⟩/J̄.
- g carries the box's sources: Dg/dl = −kg + ρ_S (eq. C4). Its formal
  solution is the attenuated 1/r² sum (C5).
- f is the background's local modulation: Df/dl = −kf + f⟨kf⟩ (eq. C7).
- Eq. C7 is unstable in the mean (eq. C8). It is stabilised by
  q(f̄)f⟨kf⟩ with q(1) = 1 and q′(1) < −1 (eq. C9), "a fiducial choice"
  q(x) = 2/x − 1.
- The OTVET forms are eqs. C10 (G, with h_G "the optically thin Eddington
  tensor from all sources inside a periodic box") and C11 (F, with
  h_F = δ^ij/3).
- The paper's case for the form: it "ensures that both f_ν and g_ν remain
  always non-negative".
- f and g are Newtonian: the Hubble terms are dropped for them and kept
  for J̄.

**The background's absorption comes from observations** (§2): "the
post-reionization evolution of the IGM would not be captured correctly if
we computed the mean free path for the cosmic background from within one
simulation box … Hence, we use the fit from [Songaila & Cowie 2010] to
account for Lyman Limit absorptions in the cosmic background; the
background is then still subject to local absorptions inside shielded
regions". Quasars and helium recombination radiation enter the background
only.

**What the background does.** Rewrite C2 as J = Lg + (J̄ − Lḡ) f.
- g is all of the box's own light: every source, the thin tensor, the
  periodic box. In our solver that is the present OTVET field, unchanged.
- The background adds A = J̄ − Lḡ times f. A is the difference between the
  externally set mean intensity and the box's own mean light.
- **With the box's own absorption the background vanishes.** Averaging C4
  gives dḡ/dt = ρ̄_S − ⟨kg⟩. With k̄J̄ = ⟨kJ⟩ = ⟨kg⟩ + A⟨kf⟩ (in L = 1
  units), the J̄ equation then gives dA/dt = −A⟨kf⟩. So A ≡ 0 if it starts
  at 0. CROC's background is therefore not a change of closure: its whole
  content is the outside mean free path, plus sources placed only in the
  background.
- **A < 0 when the outside mean free path is shorter than the box's.**
  Then J = Lg − |A|f can go negative where g is small. C2 keeps f and g
  non-negative but not J.

**Our runs against the observations.** These are the 512³ suite's λ_mfp
(the rays, `figures/taueff_21cmfast_z_snap.txt`) and the measurements in
`data/igm/igm_z5_6.txt`, transcribed from the papers the user supplied on
2026-10-10:

| z | λ_mfp, rays 512³ [pMpc] | observed [pMpc] |
|---|---|---|
| 5.0–5.1 | 45 | 9.1 (Becker+21), 9.3 (Zhu+23) |
| 5.25–5.31 | 28 | 5.4 (Zhu+23) |
| 5.49 | 15 | – |
| 5.65–5.70 | 6.6 | 3.3 (Zhu+23) |
| 5.93–6.00 | 2.6 | 0.81 (Zhu+23), 0.75 (Becker+21) |

- The box's own λ is 2–5× the observed one.
- With the background tied to the observations, A ≈ (λ_obs/λ_box − 1) ḡ ≈
  −0.5 to −0.8 ḡ. This is a crude estimate: the mean intensity scales
  roughly with λ, and Γ feeds back on λ.
- **CROC's background would then lower the box's mean intensity** toward
  the observed one. With A clamped to ≥ 0, it would not act at all.
- Either way it changes ⟨Γ⟩. A comparison with the rays at the same
  sources would then mix the closure with a different mean intensity.

## 3. Our solver, and how each scheme maps onto it

**What we solve** (`src/rt_moments.cpp`, `include/radray/kernels/moment_cell.hpp`,
`docs/rt_moments.md` §2):
- Per frequency group g (5 groups) and RT cell, the state is N_g (photons
  per mean H atom, ngam/vol_rt/n̄_H), F_g and E_g = ν N_g.
- The equations are ∂N/∂t + ĉ∇·F = S − absorption and
  ∂F/∂t + ĉ²∇·(D N) = −absorption.
- Transport is GLF (Rusanov) with wave speed ĉ, `moment_nsub` subcycles of
  Courant number C = 1/`moment_nsub`, and one cell of light travel per
  substep. After each subcycle, N ≥ 0, |F| ≤ N and E ≥ 0 are enforced.
- **The pressure is P = D N.** M1 uses D(f). OTVET uses h (6 floats per
  cell, recomputed each hydro step from that step's sources). The split
  uses P = N I/3 + min(N, N_dir)(h − I/3) (`otvet_split_pressure`).
- **Sources:** `inject_cell` puts fResc·pop2_init_g·Luv·dt_rt/vol_rt/n̄_H
  into the source cell, isotropically. The recombination photons
  (`recomb_rad = 1`) go into every cell where they were made. Both happen
  before the transport, or after it with `moment_inject = 1`.
- **Absorption:** `photo_cell` reads the deposited N, E and ΣF unchanged.
  The Γ it computes is ĉσN, so the steady state Sλ/ĉ gives the right Γ for
  any ĉ. It returns the surviving fraction fgam per cell and group, and
  every moment of the cell is scaled by fgam (`absorb_group`).
  - fgam is linear in the moments, so **any split of N into fields that
    share the cell's fgam absorbs each field exactly**. photo_cell must see
    the sum.
- **Redshifting** is per cell and per hydro step (`redshift_cell`; per
  distance with `redshift_per_distance`). Photons move between groups and F
  goes with them.
- **The budget** is one ray.log line per hydro step: injected − absorbed −
  lost = in flight, to float precision.

**The symbols:**

| G&A / CROC | ours |
|---|---|
| J_ν (angle-averaged), E | N_g in photons per n̄_H (ĉ-scaled: Γ = ĉσN) |
| F^i | F_g (in units of ĉN) |
| h^ij, E h^ij | D, P = D N |
| S_ν, ψ, ρ_* | the injection per substep (`inject_cell`; plus recombination photons) |
| S̄_ν | s̄_g: the box mean of a substep's injection |
| κ (per comoving length), e^(−κΔx) | τ_cell,g (calc_rt_1's `tau`) and fgam from photo_cell |
| ⟨κ_ν⟩𝓛 | ⟨τ_cell,g⟩ · 0.611 N_rt |
| k̄ = ⟨kJ⟩/J̄ | (absorbed per substep)/(ĉ dt_rt ΣN): the budget |
| J̄ | ⟨N_g⟩ over the box |
| ξ = ĉt/2a (eq. 18), C = 0.9 of (3/2)Δx²κ | substeps of one cell of light travel, `moment_nsub` |
| ĉ: G&A eq. 13 in the fields; CROC's J̄ with c | everything at ĉ, the `c_reduced` schedule |
| the Hubble terms: kept for J̄, dropped for f, g | `redshift_cell` per cell for every field |

Our formulation stays the explicit two-moment system. The user kept it on
2026-10-12 (`docs/PLAN.md`). G&A eq. 13 and CROC's elliptic iteration are
not adopted; only the field split and the tensors carry over.

## 4. Candidate formulations

### (a) `ga01`: two transported fields, the Λ share of the mean emission moved to an isotropic field

**The fields** (per group):
- **L:** the present field. Tensor h, all sources: stars and recombination
  photons.
- **B:** a second state (N, F, E) with D = I/3.

**The transfer.** At each injection, every cell moves
T_g = min(N_L,g, Λ_g s̄_g) photons from L to B. s̄_g is the box mean of the
substep's injection. The photons leave L in proportion to its (N, F, E),
and arrive in B isotropically (F_B unchanged), with L's mean frequency.

**Its relation to eq. 17:**
- Wherever N_L ≥ Λs̄, this is G&A eq. 17 in absolute variables (§2.1): B
  gets the uniform source Λs̄, and L the source S − Λs̄.
- Where G&A's D₂ would go negative, the transfer stops at what L holds.
  That happens in neutral gas, in the cores of shadows, and in voids where
  the local light is below Λ times the mean.
- So B gets less there than eq. 17 gives it, and L is 0 instead of
  negative. This is the only deliberate departure from G&A.

**Limits:**
- Λ = 0: T = 0 and B stays empty. The code then skips B, so L is
  **bitwise thin OTVET**.
- **With G&A's definition, Λ is negligible until overlap.** At 256³ a
  neutral cell has τ_cell ≈ 10–16 at z = 6–8 (my estimate), so
  0.611·256·⟨τ_cell⟩ exceeds 745 and exp underflows to exactly 0 while
  more than about a third of the volume is neutral. Λ stays below 10⁻⁵
  until the last ≈ 1 % of the volume ionises.
- Λ → 1: the whole mean emission is isotropic. L holds the excess of the
  local light over the mean.

**Λ, per group** (G&A's Λ_ν; harder groups get more background). The
options, which the user chooses (§7):
1. **Box, volume mean:** Λ_g = exp(−0.611 N_rt ⟨τ_cell,g⟩). This is G&A
   literally, from calc_rt_1's τ, computed each hydro step.
   - It is 0 until the last neutral islands vanish.
   - After overlap it is dominated by the residual HI, mass-weighted.
   - Estimated from the 512³ rays' λ: Λ ≈ 0.03 / 0.13 / 0.44 / 0.69 /
     0.82 / 0.88 at z = 6.2 / 6.0 / 5.7 / 5.5 / 5.25 / 5.0 (40 Mpc/h). At
     20 Mpc/h it is about 0.2–0.9. These values assume ⟨κ⟩ ≈ 1/λ_mfp;
     the volume mean exceeds that wherever self-shielded or neutral cells
     remain.
2. **Box, radiation-weighted:** Λ_g = exp(−0.611 L/λ_g), with
   λ_g = ĉ N_g/Ṅ_abs,g from the previous step's budget. This is G&A's own
   k̄ (eq. 4). It is nonzero during overlap: inside the bubbles the photons
   see only the fronts.
3. **Fixed λ** (a parameter).
4. **Observed λ_mfp(z)** from `data/igm/igm_z5_6.txt` (Becker+21, Zhu+23;
   z = 5.1–6.0 only). Outside that range it needs a rule from the user, and
   the papers for it.

**Properties:**
- **Positivity:** both fields stay ≥ 0 (T ≤ N_L), and each GLF keeps its
  own realizability. photo_cell sees N_L + N_B ≥ 0.
- **Photon conservation:** the transfer is internal, so the budget closes
  as now. ray.log gets the photons moved to B, and B's in-flight count.
- **Absorption in dense or shielded gas:** B and L share fgam, so in
  absolute units this is G&A's κ. B enters a clump isotropically from all
  sides. G&A's own shadow test (their §4.2) leaks the same way.
- **Memory:** one more N, F, E per group, 25 floats per cell. That is
  +1.7 GB at 256³ and +13.4 GB at 512³. The GLF scratch is shared.
- **Cost:** the transport doubles, since B is a second GLF (with the I/3
  pressure, no tensor read). photo_cell is unchanged, and was 51–55 % of
  the moment RT at 256³ (`docs/rt_moments.md` §6.5), so the RT costs
  about +40 %. The box means are cheap.
- **Redshift:** `redshift_cell` on each field separately. `ngam_min`
  applies per field, a small difference from one field.

**Against what exists:**
- **`otvet_split`:** a single field. It makes max(0, N − N_dir) of the
  cell's photons isotropic, with N_dir the unattenuated thin light. There
  is no λ and no memory. It acts only where N exceeds the thin light, and
  it is identical to OTVET before overlap (`otvet_closure.md` §6).
- **ga01** carries the isotropic photons as their own field, with their
  own history and transport. Its share is set by a global Λ, and it acts
  everywhere, near the sources too, once Λ > 0.
- **`otvet_closure = screen`:** a single field with one screened tensor.

### (b) `croc`: the present field plus a modulated background with an outside mean free path

**The fields:**
- **L:** the present OTVET field. All box sources, the thin tensor:
  CROC's g.
- **f:** per group (N, F; E is not needed, since f takes J̄'s spectrum)
  with D = I/3 and no source. It is absorbed by fgam and restored by the
  C9 term. Per substep, f ← f·(1 + q(f̄)·a), where a is the box-mean
  fraction absorbed and q(x) = 2/x − 1.
- **J̄_g:** an ODE per group: dJ̄/dt = s̄_g − ĉJ̄/λ_g, redshifted with
  `redshift_spectrum` each hydro step.
- **The deposit** is N_L + A f, with A_g = J̄_g − ⟨N_L,g⟩.

**Properties:**
- **Limits:** with λ from the box, A ≡ 0 and this is bitwise thin OTVET
  (§2.2). Only an outside λ (Songaila & Cowie in CROC, or `data/igm`) gives
  a background.
- **Positivity:** f ≥ 0, but A < 0 with the observed λ (§2.2). It is
  either clamped (A ≥ 0: no background, i.e. thin OTVET) or allowed, with
  the deposit clamped at 0 (a loss).
- **Conservation:** not by design. The background adds or removes A f
  photons from outside the box; ray.log would carry them as a separate
  term.
- **Memory:** +4 floats per cell per group (+1.3 GB at 256³, +10.7 GB at
  512³). The cost is as in (a).
- Optionally, as in CROC, the recombination photons go into the background
  instead of L.

### (c) `lambda`: no new field, the pressure of (a) in the thin uniform limit

P_g = N_g h + min(N_g, Λ_g N̄_g)(I/3 − h), where N̄_g is the box mean.

- In a thin, uniformly absorbing box, eq. 17 gives B = Λ N̄ isotropic and
  L = N − Λ N̄ along h: this pressure. Where N < Λ N̄, G&A's L < 0; here
  the cell is fully isotropic instead.
- **Positivity and realizability:** it is a convex mix of N h and N I/3,
  as for the split.
- **Conservation:** exact, since only the tensor changes.
- **Limits:** Λ = 0 takes the present path and is bitwise. With Λ → 1, the
  share up to the mean is isotropic.
- **Cost:** none. One box mean per group, taken with the existing plane
  sums.
- **It differs from (a)** in that the background has no transport or
  memory of its own: in each cell its share is pinned to Λ N̄. It is the
  split with a global share in place of N − N_dir.
- **It can be evaluated a priori** on the rays' dumps with
  `apps/otvet_tensor`, using the dump's Γ-weighted S0 as N. (a) cannot:
  its B field depends on the history. (c) is therefore the a-priori proxy
  for (a) in step 3.

### (d) Not recommended

- **Literal eq. 17 with a signed L:**
  - The GLF of a signed field is stable, since the system is linear with
    eigenvalues ≤ ĉ. But the realizability clamps must go for L.
  - E/N, L's mean frequency, is undefined where N_L crosses 0.
  - photo_cell needs N_L + N_B ≥ 0, so the deposit needs a clamp, i.e. a
    loss.
  - (a) is the same scheme wherever L ≥ 0.
- **Three positive fields** (stars with h, an "anti" field with h and the
  uniform source, B with I/3; total S + B − A): the same negative total.
- **CROC's J̄ with c instead of ĉ:** our stored N is ĉ-scaled. A J̄ at c
  would relax λ/c instead of λ/ĉ, inconsistently with the field it is
  added to.

## 5. Which photons are isotropic in each closure

| closure | isotropic part of the field | conserves the box's photons | extra memory |
|---|---|---|---|
| thin (now) | none | yes | – |
| `otvet_split` | max(0, N − N_dir) of the cell's photons | yes | 1 float |
| `screen` | (kernel e^(−r/λ), 27 images + continuum) | yes | – |
| (a) ga01 | the B field: Λ s̄ per substep, transported with I/3 | yes | 25 floats |
| (b) croc | A f added: the outside mean minus the box's | no (by design) | 20 floats |
| (c) lambda | min(N, Λ N̄) of the cell's photons | yes | – |

## 6. Recommendation

**(a) `ga01`, with the clipped transfer and Λ per group.** It is the only
candidate that is the published scheme. It is G&A eq. 17 in absolute
variables wherever their local field is positive.

- Like the closures we compare against the rays, it changes only how the
  photons move, at the box's own photon budget. ⟨Γ⟩ stays the box's, so the
  comparison isolates the closure.
- It stays positive and conserves photons to float precision. With Λ = 0
  it reproduces thin OTVET bitwise.
- **Λ:** I recommend option 1 (G&A's volume-mean ⟨κ⟩) as the default,
  since that is the paper's definition, with option 2 (the budget's k̄) as
  a switch. The a-priori step can show how much the choice matters
  (step 3: (c) with each Λ).
- (c) goes into `apps/otvet_tensor` as the a-priori proxy. It is cheap
  enough to also be a runtime option if wanted.

**(b) `croc` only if the paper should also show OTVET with the mean
intensity tied to observations.**
- With the box's own absorption it is thin OTVET exactly.
- With the observed λ_mfp it lowers ⟨Γ⟩ (A < 0 in our runs). That is a
  different mean, not a different closure.
- It reuses (a)'s second field, plus the ODE and the q term. About a day
  more.

**Implementation scope (step 2, after the decisions):**
- A key `otvet_background = none | ga01 | croc` (none by default), with
  `otvet_background_lambda = box | budget | fixed | observed` and
  `otvet_background_mfp` (fixed λ in RT cells).
- Allowed with `otvet_closure = thin` and `otvet_split = 0` only.
- On the serial moment solver first. With `rt_backend = amrex` it is an
  error until the moment level gets it (`docs/stage_e_moments.md`).
- Checkpoint version 10: the B field and the background budget (empty
  otherwise).
- `tensor_dump` records both fields.
- Tests in `test_moments`:
  - positivity;
  - photon conservation of the two-field transfer;
  - Λ = 0 bitwise equal to thin OTVET;
  - a uniform emissivity in a uniform box: B uniform and isotropic, L
    zero away from the sources at Λ = 1;
  - the thread and restart ctests, as for `otvet_split`.
- Byte-identity with the key off: radhydro smoke and rt_driver
  lognormal_16src against `main`, with `compare.py --rtol 0 --atol-frac 0`.

## 7. What the user decides

1. **Which formulation:** (a) ga01 (recommended), (b) croc, both, or (c)
   as a runtime option too.
2. **Where ga01's Λ comes from:**
   - the box's volume-mean opacity (G&A, recommended);
   - the box's radiation-weighted absorption (the budget);
   - a fixed λ;
   - the observed λ_mfp(z) (`data/igm`, z = 5.1–6.0; a rule is needed
     outside that range).
3. **For croc, the background's mean free path:** the box (identical to
   thin OTVET) or the observations. If observed, how to treat A < 0:
   clamp it (no background) or let the background subtract (deposit
   clamped).
4. **Recombination photons:** into L and split by Λ like the stars (G&A),
   or straight into the background (CROC).
5. **Λ per group** (G&A's Λ_ν, recommended) or one Γ-weighted Λ.

**The user's decisions (2026-11-07):** (a) `ga01` as recommended; Λ from
the box's volume-mean opacity, as G&A; Λ per frequency group; case B
recombination as the default (`recomb_case = B`: no recombination photons
are transported, so only the stars feed s̄). (b) `croc` and (c) as a
runtime option are not built; (c) stays the a-priori proxy in
`apps/otvet_tensor` (step 3).

## 8. Implementation (step 2)

**Keys** (`set_otvet_closure_param`, `src/otvet_closure.cpp`):
- `otvet_background = none | ga01` (default none).
- `otvet_background_fix`: −1 (default) = Λ from the box; a value in
  [0, 1] fixes Λ for every group (for the tests).
- ga01 needs `rt_method = otvet`, `otvet_closure = thin` and
  `otvet_split = 0`; with `rt_backend = amrex` it is an error.

**Λ** (`ga01_lambda`, `src/rt_moments.cpp`): once per hydro step in
`moment_prep`, Λ_g = exp(−0.611 N_rt ⟨τ_cell,g⟩). τ_cell,g is calc_rt_1's
per-cell optical depth at the group frequency (HI, HeI and HeII, rounded to
float as calc_rt_1 stores it), averaged over the cells in plane order.
`kGa01Path = 0.611` (`rt_moments.hpp`). ray.log gets a line
`Background : ga01 Lambda … <tau_cell> …` per hydro step.

**The substep** (`moment_step_1`):
1. Injection as before (both `moment_inject` orders). The photons injected
   per group are summed over the cells: s̄_g = Σ/N_rt³. With `recomb_rad =
   1` the recombination photons count too (G&A: all sources).
2. The transfer: every cell moves T_g = min(N_L,g, Λ_g s̄_g) to B
   (`mk::ga01_transfer`, `moment_cell.hpp`). L loses (N, F, E) in
   proportion; B gains N and E, not F.
3. Transport: L as before; B with `Closure::Iso` (P = N I/3), the same ĉ and
   `moment_nsub`. A group of B with no photons is skipped.
4. Deposit of the sum (`mk::deposit_cell2`): photo_cell sees N_L + N_B,
   E_L + E_B and the summed F.

`moment_step_3` scales both fields by the cell's fgam; `moment_redshift`
redshifts each field separately (`ngam_min` per field); the photon total
and the ray.log budget include B (`Background : in flight … moved …`).
`tensor_dump`'s diagnostic adds B with D = I/3, Γ-weighted at B's own mean
frequency.

**Memory:** B is 25 floats per cell, allocated at the first `moment_prep`
with ga01 (none otherwise).

**Checkpoint version 10:** B's five arrays and the cumulative transfer
after the star particles (empty arrays when off). A ga01 checkpoint
restarted without ga01 is an error; a version ≤ 9 checkpoint, or one
written before the RT started, restarts with an empty B.

**Checks (2026-11-07):**
- With the key off, byte-identical to `main` (6dd6a96), `compare.py
  --rtol 0 --atol-frac 0`: rt_driver lognormal_16src with rays, otvet and
  m1; radhydro smoke with rays and otvet: all PASS.
- ga01 with Λ = 0 is byte-identical to thin OTVET (rt_driver
  lognormal_16src).
- `test_moments`: Λ = 0 bitwise thin OTVET with B empty; positivity and
  realizability of both fields, and the photons of both conserved, at
  Λ = 0.3 and 1; a uniform emissivity at Λ = 1 gives a uniform B with F_B
  = 0 and an empty L; the box Λ (1 ionised, exp(−0.611 N τ/N³) for one
  partly neutral cell, ≈ 0 neutral); 1 and 3 threads bitwise.
- ctest `thread_rt_otvet_ga01_lognormal_16src` and
  `restart_smoke_otvet_ga01` (Λ fixed at 0.5, since the box's Λ is ≈ 0 in
  these neutral problems).
- Fast suite: 171/180 pass. The 9 failures are the `platform_ref` stored
  references, which fail the same way on `main` on klone (another libm;
  CI runs them).
