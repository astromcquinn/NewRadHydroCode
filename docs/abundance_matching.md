# Abundance-matched UV sources

*Status (2026-10-08): **draft**. The code is in place and tested, but the
luminosity-function parameters are not yet confirmed against the papers, and
the z > 10 (JWST) entries are still missing. Do not use the fitted
parameters for science until `data/uv_lf/lf_params.txt` has only
`confirmed` rows.*

This is a third source model, `sfr_model = abundance`, beside the
reconstructed Trac & Cen model (`global`, the default) and `dynamical`
(`docs/halos_and_sources.md` §3). Selecting it changes nothing else: the
existing models are bit-identical to before.

## 1. Method

1. **Observed LFs** (`data/uv_lf/lf_params.txt`): one Schechter or
   double-power-law fit per (redshift, paper), with the magnitude range the
   data constrain and a weight. φ* is in the observers' Mpc⁻³ mag⁻¹ (h = 0.7).
2. **Matching** (`abundance_match`, `include/radray/abundance_matching.hpp`):
   at each LF redshift, M_UV(M_halo) from n_halo(>M) = n_gal(<M_UV).
   * n_halo comes from the analytic Sheth–Tormen (default) or Tinker et al.
     2008 (Δ = 200m) mass function, in the simulation cosmology
     (cosmo_parameters.f90: Ωm = 0.27, h = 0.7, σ8 = 0.8, n = 0.96), BBKS
     transfer function.
   * The small cosmology difference from the observers' (Ωm = 0.3, h = 0.7) is
     ignored. It is a few percent in volume.
   * No scatter in L at fixed M, and no duty cycle.
   * Optional dust correction to intrinsic L: Meurer et al. 1999
     IRX–β, A_UV = 4.43 + 1.99 β with a linear β(M_UV). It is off by default
     (§4).
   * L_ν = 4π(10 pc)² 10^(−0.4(M_AB + 48.60)), i.e. M_UV = 51.595 − 2.5 log L.
3. **Smooth fit** (`fit_luv_model`, Nelder–Mead in log L over the matched
   points inside each entry's fitted magnitude range):

   ```
   L_UV(M, z) = 2 L_p / ((M/M_p)^-a + (M/M_p)^-b),   x = clamp(z, z_lo, z_hi) - z0
   log10 L_p = l0 + l1 x + l2 x²,   a = a0 + a1 x,   log10 M_p = m0 + m1 x,   b
   ```

   M is in M☉ and L in erg/s/Hz. z is clamped to the redshift range of the LF
   entries, so the relation is **held fixed** outside the constrained
   redshifts, rather than extrapolated (§4).
4. **Ionising photons.** The RT emits `fResc * pop2_init * L` photons per
   second. That is a fixed 0.2 × 10^25.2 photons per erg/Hz, split over 5
   groups (cosmo_parameters.f90). The source model stores, for each halo:
   * `L = L_UV × (ξ_ion / Σ pop2_init) × (f_esc / fResc)`, so
     N_ion,esc = f_esc ξ_ion L_UV exactly, without touching the RT or its
     spectrum;
   * `sfr = L_UV / l_per_sfr`, with the default 8×10²⁷ (Madau et al. 1998,
     Salpeter). It only drives the metal injection and the
     `sfr > sfr_min` gate of `calc_halo`;
   * `mstar = min(sfr × t_age, f_b M)`, informational only.

   ξ_ion and f_esc are degenerate: only f_esc ξ_ion enters the RT, and f_esc
   is calibrated to τ_es. ξ_ion matters only for quoting f_esc.

## 2. Usage

```
build/am_fit --lf data/uv_lf/lf_params.txt --out am_out [--hmf st|tinker] \
             [--extend-faint -13] [--dust 1 --beta0 -2.0 --dbeta -0.1] \
             [--free l0,l1,...] [--set am_b=0.3] [--only-confirmed]
python3 tools/analysis/plot_am.py am_out          # am_out/am_fit.png
```

`am_out/am_params.txt` holds `params.txt` lines (`sfr_model = abundance`,
`am_l0 … am_b`, `am_z0`, `am_zlo`, `am_zhi`), to which `xi_ion` (log10)
and `f_esc` are added for a run. `am_out/lf_check.txt` lists the input LF
against the LF the fit implies.

| key | default | meaning |
|---|---|---|
| `sfr_model` | `global` | `abundance` selects this model |
| `am_l0 … am_b` | see `LuvModel` | fit parameters (from `am_fit`) |
| `am_z0`, `am_zlo`, `am_zhi` | 8, 0, ∞ | pivot and clamp of the z dependence |
| `xi_ion` | 25.2 | log10 ξ_ion, photons s⁻¹ / (erg s⁻¹ Hz⁻¹) |
| `f_esc` | 0.2 (`fResc`) | escape fraction (calibrated to τ_es, step 4) |
| `l_per_sfr` | 8e27 | L_UV per M☉/yr (metals only) |
| `tvir_min`, `m_min` | 1e4 K, 0 | halos below do not emit (as for the other models) |

## 3. Draft fit (HST rows only, unconfirmed)

![draft fit](figures/am_fit_draft_unconfirmed.png)

With Sheth–Tormen and all eight parameters free:
* the rms is 0.04 dex in log L;
* the implied LFs match the inputs to within 0.1–0.17 dex over most of the
  fitted range, and to within 0.35 dex at the steep bright end.

The resolution limit of 128³ in 10 Mpc/h is 20 particles, i.e. 7.1×10⁸
M☉/h = 1.0×10⁹ M☉. These halos map to **M_UV ≈ −12 to −14**, far fainter
than the HST blank-field limit (M_UV ≈ −17). So most sources in the box lie
on the extrapolated faint end. At logM = 9, the choices below shift M_UV by:

| choice | z = 6 | z = 8 | z = 10 |
|---|---|---|---|
| fit observed range only (default) | −12.1 | −13.2 | −13.6 |
| `--extend-faint -13` (fit the extrapolated α down to −13) | −12.4 | −13.8 | −14.0 |
| + Tinker instead of Sheth–Tormen | −12.7 | −14.3 | −14.7 |

That is up to about 1 mag (0.4 dex in L). Much of a uniform offset is
absorbed by the f_esc calibration. The redshift and mass dependence is not
absorbed, so it shapes the timing and morphology of reionisation.

## 4. Decisions needed (user)

1. **LF parameters:** confirm the rows of `data/uv_lf/lf_params.txt`, or
   supply the papers (§5).
2. **z > 10 data set:** Donnan+24, Finkelstein+24, McLeod+24, Harikane+24
   or Whitler+25, or a combination.
3. **Faint end:** fit only the observed range, or include the
   extrapolated faint-end slope (`--extend-faint`), and down to which
   magnitude (lensing LFs reach about −13 to −15 at z ≈ 6–8)?
4. **Halo mass function** for the matching: ST, Tinker, or the measured
   mass function of the step-2 run (not yet implemented; it would remove
   the analytic-HMF vs FoF mismatch at the resolution limit).
5. **Dust:** none (default), or the Meurer IRX–β with which β(M_UV, z)
   relation (e.g. Bouwens et al. 2014)?
6. **ξ_ion:** 10^25.2 (the original's `pop2_init`, the canonical
   Robertson et al. 2013 value), or a JWST-era value (log ξ_ion ≈ 25.3–25.8;
   from memory, unconfirmed)?
7. **z outside the data:** hold L(M) fixed above the highest LF redshift
   (current), or extrapolate the linear z dependence?

## 5. Parameters used, with provenance

The cloud sessions cannot reach arXiv, ADS or the journals. "Snippet" means
the number appeared in a quoted search-result excerpt of the paper; "memory"
means recalled, not verified. **All need confirmation.**

**UV LFs (`data/uv_lf/lf_params.txt`),** Bouwens et al. 2021, AJ 162, 47
(arXiv:2102.07775), Schechter fits:

| ⟨z⟩ | M* | φ* (10⁻³ Mpc⁻³) | α | provenance |
|---|---|---|---|---|
| 5.9 | −20.93 ± 0.09 | 0.19 (+0.05/−0.04) | −1.93 ± 0.04 | memory |
| 6.8 | −21.15 ± 0.13 | 0.12 | −2.06 ± 0.05 | memory |
| 7.9 | −20.93 ± 0.28 | 0.09 (+0.09/−0.05) | −2.23 ± 0.20 | snippet |
| 8.9 | −21.15 (fixed) | 0.021 (+0.014/−0.009) | −2.33 ± 0.19 | snippet |
| 10.2 | −21.19 (fixed) | 0.0042 (+0.0045/−0.0022) | −2.38 | M*, φ*: snippet (Oesch et al. 2018 row); α: memory |

The Bouwens+21 evolving fit was partly seen in snippets and is not used:
M* = −21.03 − 0.04(z−6), α = −1.94 − 0.11(z−6),
φ* = 0.40×10⁻³ × 10^(−0.33(z−6) − 0.024(z−6)²).

JWST z > 10: **none entered.** Candidates are listed in the data file. Seen
in snippets: Donnan+24 fixes β = −4.05 at z = 10, and McLeod+24 fixes
α = −2.35 at z ≈ 11.

**Other published constants (from memory, unconfirmed):**

| constant | value | source |
|---|---|---|
| Tinker HMF, Δ = 200m | A0 = 0.186, a0 = 1.47, b0 = 2.57, c = 1.19; A ∝ (1+z)^−0.14, a ∝ (1+z)^−0.06, b ∝ (1+z)^−α, log α = −(0.75/log(Δ/75))^1.2 | Tinker et al. 2008, ApJ 688, 709 |
| Sheth–Tormen | A = 0.3222, a = 0.707, p = 0.3 | Sheth & Tormen 1999 (already in the code) |
| IRX–β | A_UV = 4.43 + 1.99 β | Meurer, Heckman & Calzetti 1999 |
| L_UV / SFR | 8.0×10²⁷ erg/s/Hz per M☉/yr | Madau, Pozzetti & Dickinson 1998 (already in the code) |
| ξ_ion | 10^25.2 Hz/erg | Robertson et al. 2013; equals the original's `pop2_init` |

## 6. Tests

`tests/test_abundance_matching.cpp` (CTest `abundance_matching`) checks:
* the AB conversion;
* the LF integrals;
* the abundance identity at every matched point;
* that `model_lf` integrates back to the halo abundance;
* recovery of known fit parameters;
* the Tinker mass function against Sheth–Tormen and its z freeze;
* the source model: `L`, `sfr`, the f_esc and ξ_ion scalings, the z clamp
  and the `params.txt` keys.
