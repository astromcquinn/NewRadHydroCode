# Plane-parallel / RT-domain mode: specification

*Taken from the user's papers, which were uploaded to the first thread. The
PDFs are not in the repo: ask the user if you need them.*

- **D'Aloisio, McQuinn, Trac, Cain & Mesinger 2020** (arXiv:2002.02467),
  "Hydrodynamic Response of the IGM to Reionization", §3.1–3.2 and
  Appendices A–B. Referred to below as **DM20**.
- **Cain, McQuinn, Scannapieco, D'Aloisio & Trac 2026** (Phys. Rev. D 113,
  023050; arXiv:2504.21082), "Kiloparsec-scale turbulence driven by
  reionization may grow intergalactic magnetic fields", §II. Referred to
  below as **C26**. It follows the procedure of Cain et al. (ref. [17] there)
  and its Supplemental Material §C 2, which we do not have.

Both papers used a modified version of this RadHydro code. Nothing in that
modified version is in `reference/Radhydro`, so the features below are new.
Implement them behind parameters; the defaults must leave current behaviour
and all legacy bit-for-bit tests unchanged.

## What the papers did

**Geometry and sources**

- **Plane-parallel ray tracing on the uniform RT grid,** at the same
  resolution as the hydro grid (Nrt = Ngas = Ndm).
- **Rays can be sent from any of the 6 directions at once.** Both papers sent
  rays along two orthogonal directions. This gives planar I-front
  propagation without unrealistic single-direction shadowing.
- **RT domains:**
  - The box is divided into cubical RT domains: DM20 used 32³ domains of
    side 32 kpc/h in a 1.024 Mpc/h box with 1024³ cells, i.e. 32 cells per
    domain.
  - Ray sources sit on the domain boundaries, and rays are deleted when they
    reach the far end of their domain.
  - So all gas is ionised at nearly the same time, the intensity stays
    controlled despite self-shielding systems, and there is no strong
    gradient across the box.
  - A single domain (Ndom = 1, sources on the box faces) is also used, e.g.
    the Fig. 1 illustration in DM20 and one convergence run.
- **No periodic RT boundary is needed.** Rays never cross the box boundaries
  because they end at domain edges (C26 note [49]). Hydro and gravity stay
  periodic.
- **Source cells** are excluded from all analysis quantities (DM20).

**Radiation**

- **Spectrum:**
  - J_ν ∝ ν^−1.5 between 1 and 4 Ry, in 5 frequency bins.
  - The bin edges are chosen so that each bin holds the same photon number
    density.
  - This replaces the code's Pop II `freq_init`/`pop2_init` in this mode.
- **Normalisation:**
  - Set the intensity so that the HI photoionisation rate in the source
    cells is a chosen Γ_HI. The parameter is Γ_−12 in units of 10⁻¹² s⁻¹;
    DM20 used 0.1, 0.3 and 3. C26 used Γ_HI = 3×10⁻¹³ s⁻¹ "everywhere
    except within self-shielded systems".
  - Photon number densities do not dilute with cosmological expansion, so
    Γ_HI stays fixed in time.
- **Turn-on:**
  - Radiation switches on at redshift z_re (DM20: 12, 8, 6; C26: about 7).
  - Before z_re the run is hydro + gravity only, with no RT.

**Speed of light**

- **Adaptive reduced speed of light (DM20):**
  - c_sim = 0.1c until the volume-weighted mean x_HI of the box drops below
    1%, then c_sim = 0.01c.
  - DM20 Appendix A found this accurate to better than 10% against c_sim = c.
  - C26 uses a reduced c after the I-fronts have crossed the domains.
- **Implementation:** a parameterised schedule (c_before, threshold,
  c_after) built on the `c_reduced` parameter (first-thread follow-up PR).

**Gas physics and initial conditions**

- **Cooling:** full primordial heating and cooling (DM20: "tracks all of the
  relevant heating and cooling processes for gas of primordial composition",
  citing D'Aloisio et al. 2018). Cold, dense, self-shielded gas is the point,
  so no star-formation sink.
- **Pre-heating option:** a temperature floor T_min imposed from a set
  redshift, mimicking X-ray heating. DM20 used T_i = 1000 K from z = 20; C26
  used T_min of 100 and 1000 K for z < 15. Otherwise the gas evolves
  adiabatically after decoupling.
- **Initial conditions:**
  - Start at z = 300 with first-order perturbation theory.
  - Use **separate transfer functions for baryons and dark matter**
    (CAMB). `tools/ics/make_ics.py` currently uses one Eisenstein–Hu
    no-wiggle spectrum. Add separate baryon and DM transfer functions, either
    from a CAMB/CLASS table read from a file (CAMB may not be installable
    here) or from Eisenstein & Hu fits with baryon suppression.
  - Ndm = Ngas.
- **DC modes (DM20 §3.2):**
  - A uniform box-scale overdensity δ/σ (δ/σ = 0, ±√3 are the
    Gauss–Hermite nodes), with full non-linear evolution of the DC mode, as
    in the separate-universe approach of Gnedin et al. 2011.
  - DM20 Appendix B describes averaging over DC-mode runs.
  - Implement this as a later, separate step: it changes the cosmology and
    expansion of the box.

## Typical sizes

| Paper | Box | Cells (gas = RT = DM) | Δx | Notes |
|---|---|---|---|---|
| DM20 production | 1.024 Mpc/h | 1024³ | 1 kpc/h | 32³ domains of 32 kpc/h |
| DM20 convergence | 0.256 Mpc/h | 64³–1024³ | | c_sim = c, 0.1c, 0.01c tests |
| C26 fiducial | 0.25 Mpc/h | 1024³ | 244 pc/h | runs to z = 4 |
| C26 high-res | 0.25 Mpc/h | 2048³ | 122 pc/h | |

## Tests to add (physics, not Fortran comparison)

- **Uniform-gas I-front:** a plane-parallel front in uniform gas, compared
  against the analytic front position and speed (photon flux balance, with
  recombinations).
- **Slab attenuation:** photons through a uniform slab, compared with e^(−τ)
  in each frequency bin.
- **Photon conservation:** photons injected = absorbed + escaped at domain
  ends + still in flight.
- **Γ_HI normalisation:** the measured Γ_HI in the source cells equals the
  set Γ_−12.
- **Domains:** with Ndom > 1, all domains ionise at nearly the same time, and
  rays never cross domain edges.
- **Adaptive c:** the schedule switches at x_HI = 1%, plus a diagnostic for
  I-fronts limited by c_sim.
- **Physics sanity check:** a scaled-down DM20 problem (e.g. 0.256 Mpc/h,
  128³, z_re = 8, Γ_−12 = 0.3). Compare the post-I-front temperature
  (about 20,000 K) and the ionised-gas clumping factor qualitatively with
  DM20.

## Open questions for the user

- **Cooling in dense gas:** did the modified code keep the original's
  no-cooling-above-200× switch, or was cooling fully on?
- **Domain sources:** do rays launch from one face of each domain per
  direction, or from both faces of each domain?
- **Supplementary procedure:** do we need Cain et al. [17] / C26
  Supplemental §C 2, e.g. for a different turn-on or domain procedure?
