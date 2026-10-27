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
modified version is in `reference/fortran/Radhydro` (the private submodule), so the features below are new.
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

## Implementation (2026-10-26, branch `claude/eager-sagan-h6a0tu`)

Built as specified above, with the user's answers to the open questions
(2026-10-26): cooling fully on (`cool_suppress_overdensity = 0` in the
plane-parallel problems; the code default stays 200), one face per domain
and direction (the default `plane_faces = +x +y`), no Supplemental
procedure (hydro and gravity to z_re, then every source face switches on
at the full Γ_HI at once; the existing `c_reduced_schedule` gives the
reduced ĉ after the fronts have crossed), and fixed 8³ domains for the C1
ladder (31.25 kpc/h at every resolution). A single front crossing the
whole box, which the user also asked for, is `plane_ndom = 1` with one
face (`plane_faces = +x`).

### Parameters (`include/radray/config.hpp`, `set_plane_param` in `plane_parallel.hpp`; both drivers)

| key | default | meaning |
|---|---|---|
| `plane_parallel` | 0 | 1: the faces of cubic RT domains are the only sources (no halo finder, no halo sources); legacy_rt = 0, rays, `ray_encoding = float`, no `recomb_rad`, uniform grid only |
| `plane_ndom` | 1 | domains per dimension; N/ndom cells each (a multiple of `rt_cells_per_substep`) |
| `plane_faces` | `+x +y` | the emitting faces of every domain, any of `+x -x +y -y +z -z` |
| `plane_gamma12` | 0.3 | Γ_HI in the source cells, 10⁻¹² s⁻¹ (C26: 0.3) |
| `plane_z_re` | −1 | the turn-on redshift (required by `radhydro`; `rt_driver` runs the RT every step) |
| `plane_alpha`, `plane_nu_max_ry` | 1.5, 4 | J_ν ∝ ν^−α between 1 Ry and ν_max |
| `photon_consume_dense` | 0 | 1: every photoionisation consumes its photon (the original's `rho_source` factor, which lets cells denser than 50× the mean ionise without consuming photons, is off); the plane-parallel problems use 1 |
| `dt_max_myr`, `rt_substeps_max` | 0, 0 | the step caps of the T4 test (below); in plane mode they act from the turn-on |
| `land_outputs` | 0 | 1: the hydro step before a `zdata1` output ends exactly on it (so runs with different step caps are compared at the same times); the plane-parallel problems use 1 |

`tools/production/make_production.py --plane-parallel` writes a production
`params.txt` of the mode (`--plane-ndom`, `--plane-faces`,
`--plane-gamma12`, `--plane-z-re`, `--dt-max-myr`, `--rt-substeps-max`):
DM20's schedule `c_reduced_schedule = 0.1 0.99 0.01` unless given, cooling
on, `photon_consume_dense = 1`, no abundance matching.

### Geometry and rays

* Rays: one per emitting face and face cell per RT substep
  (`rrk::plane_rays`, `include/radray/kernels/ray_kernels.hpp`), at the
  centre of the cell's lower face along the axis (cell index ≡ 0 mod L,
  so the domains tile the periodic box), directed along +axis into the
  domain or along −axis into the domain below (whose high face it is);
  `r = 0`, `dr = L`, `nside = 1` and `ipix` = the face as markers. With
  `plane_parallel = 1` `rrk::classify` never splits or merges any ray; a
  ray is deleted once it has crossed L cells (`rrk::plane_done`, checked
  after `cast_3` by the serial `ray_cast_3` and the backend's
  `substep_advance`), or when it is fainter than `ngam_min`.
* Serial RT: `plane_emit` (`src/ray_tracing.cpp`) inserts the rays at the
  head of their cells' chains before the sweep of `rt_step_1` (in parallel
  over a plane with `deterministic_parallel = 1`); the sweep casts them.
  AMReX backend: `rrk::step1_count/step1_write` take the face rays of a
  cell as extra input rays (`src/amrex/rt_level.cpp`).
* The rays are not redshifted (they live one domain crossing; the
  redshift tables assume the Pop II bands), so Γ_HI stays fixed in time.
* The pool: up to (faces) N³ / k rays in flight; `check_plane_param`
  requires `nray_max_factor` ≥ 1.25 faces / k (+1).

### Spectrum and normalisation (`plane_spectrum`, `src/plane_parallel.cpp`)

J_ν ∝ ν^−α between 1 Ry (the code's `HI_eV` = 13.60 eV) and 4 Ry cut into
5 bins of equal photon number (α = 1.5: edges 13.60, 15.46, 18.12, 22.34,
30.35, 54.40 eV). Each group's effective frequency is its bin's
photon-weighted mean (14.48, 16.70, 20.05, 25.84, 39.50 eV); the cell
physics evaluates the Verner cross sections there, as it does for every
ray, and the optical depths of the grid (`calc_rt_1`, `post_rt`, the
level's `sig_init`) use these frequencies instead of `freq_init`
(`group_frequencies`). 4 Ry = 54.40 eV is below the code's HeII edge
(54.42 eV): no HeII ionisation, as in the papers. Photons per ray and
substep: N = Γ_HI dt_rt / (⟨n_H⟩ len_rt σ_eff), σ_eff = Σ_n f_n σ_HI(ν_n),
split by f_n = 1/5: the photons that cross one cell face per substep, so
that the photoionisation rate Γ = n_γ σ ĉ of a lit cell is Γ_HI
(independent of `rt_cells_per_substep`). `<outdir>/plane_spectrum.txt`
lists the bins, cross sections and fluxes.

### Bookkeeping

`stats.txt` gains `plane_emitted`, `plane_absorbed` (deposited minus
transmitted, summed over the cells) and `plane_escaped` (the photons of
the rays deleted at their far face) per hydro step, in the rays' units
(per ⟨n_H⟩ per RT-cell volume); the absorbed and escaped sums are exact
fixed-point sums (thread- and rank-independent). With the photons in
flight (`ngam_rays_total`): emitted = absorbed + escaped + Δ(in flight)
to a few 10⁻⁵ (the rays deleted fainter than `fgam_min` of their photons).

### Turn-on and the step caps (`timestep_core`, `src/cosmology.cpp`)

In `radhydro` the RT runs in every hydro step whose start is at or below
`plane_z_re`; the step before lands on z_re exactly. `dt_max_myr` caps the
hydro step in Myr and `rt_substeps_max` shortens it to hold at most that
many RT substeps (1: the hydro step is one light crossing of
`rt_cells_per_substep` cells at ĉ, i.e. the RT and the hydro are coupled
every substep); in plane mode both act from the turn-on, otherwise always.

### Tests (`tests/CMakeLists.txt`)

| test | what |
|---|---|
| `plane_rt` (`tests/plane_check.py`, RT-only driver, uniform gas built per check) | Γ_HI normalisation to 5e-8 in the source cells and to the attenuation elsewhere, for 1, 2 and 4 cells per substep; the bins and the optical depths; slab attenuation per bin = Π(1 − τ_j) to 4e-4 and e^−τ to 2 %; one front across the whole box against the photon budget (thin gas: 0.5 %; 100× the mean density with recombinations: the column = (photons − the code's recombinations)/n_H and α_eff against case B); the light-limited front at ĉ; photon conservation (3e-5); 2³ domains ionising identically (1e-6) with no ray crossing a domain edge; 1 vs 2 threads byte-identical |
| `plane_full` (`tests/plane_full_check.py`, `radhydro`, 16³ uniform box) | the landing on z_re, no photons before it, DM20's schedule switching once at x_HI,V = 1 % (substeps 1263 → 133), the caps (`dt_max_myr = 2`, `rt_substeps_max = 50`) with the T_M(z), x_HII,V(z) of the uncapped run to 0.6 % |
| `plane_dm20` (nightly) | the scaled-down DM20 problem: 0.256 Mpc/h, 32³, z_re = 8, Γ_−12 = 0.3, cooling on, to z = 7: x_HII,V > 0.99, T_M of the ionised gas in 1–3 × 10⁴ K, the clumping factor in 1–10 |
| `amrex_plane_rt`, `amrex_plane_full` | the same on the AMReX backend with 2 ranks; the full loop against the serial run (x_HII,V and T_M to 1e-4, fHI to 1e-3) |

The bitwise suite is untouched: every new flag defaults to the original.

### Tools

* `tools/analysis/plane_front.py`: folded front profiles (⟨x_HI⟩, ⟨T⟩
  from the source face), the front position per output and its speed in
  units of ĉ (the light-limited diagnostic).
* `tools/analysis/velocity_power.py`: the kinetic energy spectrum E(k) of
  the gas velocity (C26 Fig. 3), with `--ref` for digitised points.
* `tools/analysis/t4_coupling.py`: the T4 test (make / run / analyse).
* `tools/hyak/c1_ladder.sh`: the C1 ladder and the T4 variants on klone.

### Results (2026-10-27)

* **T4** (`docs/STATUS.md`, "Done 2026-10-26", the table): coupling the
  heating to the hydro every 1 Myr reproduces the per-substep coupling to
  8 % in E(k) at 25 Myr after z_re and to 2–5 % later; 10 Myr leaves
  40–120 % differences during the first 100 Myr; the uncapped hydro step
  (44 Myr) lags by factors 2–4 in E(k) for 200 Myr. The ionisation
  profiles do not depend on the coupling. Figure
  `docs/figures/plane_t4_N64.png`, table `docs/data/plane/t4_N64.txt`.
* **C1 at 64³** (local; `docs/figures/plane_c1_N64_ek.png`): the fronts
  cross in the first step, T_M = 17 900 K behind them, E(k) steepens from
  k^−3.6 to k^−5 with no k^−5/3 range, as C26 Fig. 3 (bottom) at Δx ≥ 0.5
  kpc/h. The 128³–512³ rungs run on klone (`tools/hyak/c1_ladder.sh`).
