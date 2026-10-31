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
| `plane_transport` | `rays` | `sweep`: the packets in flight are a field of 5 floats per cell and lit face, shifted along its axis by `rt_cells_per_substep` cells every substep through the same cast, physics and attenuation kernels (`src/plane_sweep.cpp`; the backend's `PF` fields). No ray pool (no `Nray_max` limit), a fifth of the memory, no sorting or migration; the same photons as the rays to float round-off (test `sweep`). The packets in flight are not checkpointed: a restart refills the domains within one crossing. `make_production.py --plane-parallel` writes `sweep` |
| `plane_ndom` | 1 | domains per dimension; N/ndom cells each (a multiple of `rt_cells_per_substep`) |
| `plane_faces` | `+x +y` | the emitting faces of every domain, any of `+x -x +y -y +z -z` |
| `plane_gamma12` | 0.3 | Γ_HI in the source cells, 10⁻¹² s⁻¹ (C26: 0.3) |
| `plane_z_re` | −1 | the turn-on redshift (required by `radhydro`; `rt_driver` runs the RT every step) |
| `plane_alpha`, `plane_nu_max_ry` | 1.5, 4 | J_ν ∝ ν^−α between 1 Ry and ν_max |
| `photon_consume_dense` | 0 | 1: every photoionisation consumes its photon (the original's `rho_source` factor, which lets cells denser than 50× the mean ionise without consuming photons, is off); the plane-parallel problems use 1 (user decision, 2026-10-29) |
| `dt_max_myr`, `rt_substeps_max` | 0, 0 | the step caps of the T4 test (below); in plane mode they act from the turn-on |
| `land_outputs` | 0 | 1: the hydro step before a `zdata1` output ends exactly on it (so runs with different step caps are compared at the same times), and the run ends once the last output is written (its clock can stop a few ulp above that redshift); the plane-parallel problems use 1 |
| `plane_turnon_cfl`, `plane_turnon_temp`, `plane_turnon_courant` | 1, 2e4, 0.05 | the step that starts at `plane_z_re` is at most `plane_turnon_courant` × the CFL step with the sound speed of ionised gas at `plane_turnon_temp` (K) added (below; 2026-10-28; the T4 test needs a Courant fraction ≈ 0.05 to remove the turn-on lag, the default by the user's decision); `plane_turnon_cfl = 0` = the uncapped turn-on step of before |
| `T_min`, `z_floor` | 0 (off), unset | the pre-heating option of §1: a temperature floor in K from z < `z_floor` on, in the RT cell physics and the hydro, a plain switch-on (user decision 2026-10-29; `docs/production_run.md` §10, `docs/thermal_history.md` §5). C26: 100 or 1000 K for z < 15; `make_production.py --T-min K` writes `z_floor = 15` |

`tools/production/make_production.py --plane-parallel` writes a production
`params.txt` of the mode (`--plane-ndom`, `--plane-faces`,
`--plane-gamma12`, `--plane-z-re`, `--dt-max-myr`, `--rt-substeps-max`,
`--T-min`, `--z-floor`, `--thermal-floor`): DM20's schedule
`c_reduced_schedule = 0.1 0.99 0.01` unless given, cooling on,
`photon_consume_dense = 1`, `thermal_floor = 0.1` (the hydro floor and the
clamp of the RT physics floor, instead of the original's 1 K: the
mean-density adiabat is below 1 K at z < 8; user decision 2026-10-29,
`docs/thermal_history.md` §4), no abundance matching.

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

### The sweep transport (`plane_transport = sweep`, 2026-10-27)

A packet about to enter cell e along a face's axis has travelled m(e)
cells from its source face (e mod L for a + face, L − 1 − e mod L for a −
face); packets sit at m = 0 mod k, the one at m = 0 is the source packet
and is never stored. In `rt_step_1` every packet is cast as the ray it
stands for (`rrk::plane_packet_ray`, at the face it is crossing, with
`rrk::cast_1`), so the deposits are the rays' to the last bit; in
`rt_step_3` it is attenuated with `rrk::cast_3`, moved k cells on (in
place along its column, against the direction of motion), counted as
escaped at m + k = L, or dropped when faint. One column of cells per
thread (a packet deposits only into its own column), fixed-point sums:
thread-independent. On the backend the fields have k ghost cells along
the axis (`FillBoundary` before the advance), the deposits go into the
level's fixed-point accumulator, and the fields follow their cells on a
rebalance. Memory: faces × 5 floats per cell (two faces at 1024³: 43 GB
over the ranks) against 80 B per ray in flight (170 GB).

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

**The turn-on step (`plane_turnon_cfl`, `plane_turnon_temp`; user
request 2026-10-28).** Every global step runs the hydro and the PM first
and then the RT over the same interval (`apps/radhydro.cpp`), so the RT
heating reaches the hydro only at the next step, and the step is sized at
its start: dt = 2 dt_gas_factor / cmax_gas (one cell crossing at the
fastest |v| + c_s), limited by `ra_max` and the PM. The step that starts
at z_re therefore sees only the cold gas (v_max ≈ 16 km/s at z = 7 in the
0.25 Mpc/h box: 44 Myr at 64³), and the freshly heated gas sits that long
without a hydro response (the T4 lag of ⟨v²⟩ below). After it cmax_gas
holds the hot gas, so only this step needs care. With `plane_turnon_cfl =
1` (the default) the step that starts at the turn-on (a = a_re to 1e-9,
the first step with the RT on) is at most

  dt = C_on 2 dt_gas_factor / (cmax_gas + c_hot),  c_hot² = γ k T_hot / (μ m_H),

C_on = `plane_turnon_courant` (default 0.05, user decision 2026-10-28; 1 is the plain CFL step with c_hot), T_hot = `plane_turnon_temp`
(2e4 K), μ = 1/(2X + 3Y/4) = 0.588 (H II and
He III; the code's γ = 5/3, X = 0.76, Y = 0.24, k_cgs, mH_cgs of
`constants.hpp`), i.e. the hydro's own c_s² = γ(γ − 1) e/ρ of fully
ionised gas at T = `temp3_unit` e/ρ: 21.6 km/s at 2e4 K. c_hot enters in
code units through `vel_unit` = x_unit / time_unit, the conversion
`calc_timestep` applies to vmax_gas (checked in `plane_full`: dt_sim
(cmax_gas + c_hot) = C_on to 1e-7 with c_hot converted independently in
the test). It shortens da like the other caps and comes after them, so where
it does not bind (`dt_max_myr` ≤ the CFL step, or a box where `ra_max`
is shorter) nothing changes bit for bit. As every CFL step of the code it
holds the code velocities fixed through the step: in physical time dt
exceeds Δx_phys/(v_max + c_hot) at z_re by the growth of time_unit ∝ a²
over the step (≈ da/a, 5 % for the 45 Myr test step). The rule is
`timestep_core` (`src/cosmology.cpp`), which the serial and the AMReX
paths share (the backend broadcasts cmax_gas before it): no separate
backend code. **Restart:** the trigger is a function of the step's start
a only and the clock (cmax_gas included) is checkpointed, so a restart
from the checkpoint written at z_re repeats the capped step bit for bit
(checked: stats and fields of the turn-on step identical to the
uninterrupted run), and a restart anywhere else cannot fire it (the only
step that starts at a_re is the turn-on step). `make_production.py
--plane-parallel` needs nothing (the default is on). With `plane_parallel
= 0`, or without `plane_z_re`, nothing changes.

### Tests (`tests/CMakeLists.txt`)

| test | what |
|---|---|
| `plane_rt` (`tests/plane_check.py`, RT-only driver, uniform gas built per check) | Γ_HI normalisation to 5e-8 in the source cells and to the attenuation elsewhere, for 1, 2 and 4 cells per substep; the bins and the optical depths; slab attenuation per bin = Π(1 − τ_j) to 4e-4 and e^−τ to 2 %; one front across the whole box against the photon budget (thin gas: 0.5 %; 100× the mean density with recombinations: the column = (photons − the code's recombinations)/n_H and α_eff against case B); the light-limited front at ĉ; photon conservation (3e-5); 2³ domains ionising identically (1e-6) with no ray crossing a domain edge; 1 vs 2 threads byte-identical |
| `plane_full` (`tests/plane_full_check.py`, `radhydro`, 16³ uniform box) | the landing on z_re, no photons before it, DM20's schedule switching once at x_HI,V = 1 % (substeps 1263 → 133), the caps (`dt_max_myr = 2`, `rt_substeps_max = 50`) with the T_M(z), x_HII,V(z) of the uncapped run to 0.6 %; the turn-on step (`plane_turnon_temp = 4e4`, `plane_turnon_courant = 0.5`, as the default does not bind in this cold uniform box: the `ra_max` step is shorter): dt_sim (cmax_gas + c_hot) = 0.5, dt ≤ 0.5 Δx_phys/(v_max + c_hot) up to the expansion over the step, the steps before z_re bit for bit those of the default run, and with `plane_turnon_cfl = 0` the old step (96.8 against 37.1 Myr) |
| `plane_dm20` (nightly) | the scaled-down DM20 problem: 0.256 Mpc/h, 32³, z_re = 8, Γ_−12 = 0.3, cooling on, to z = 7: x_HII,V > 0.99, T_M of the ionised gas in 1–3 × 10⁴ K, the clumping factor in 1–10 |
| `amrex_plane_rt`, `amrex_plane_full` | the same on the AMReX backend with 2 ranks; the full loop against the serial run (x_HII,V and T_M to 1e-4, fHI to 1e-3), and the capped turn-on step: the steps to it bit for bit those of the serial run (the later CFL steps to 4e-8) |

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
  8 % in E(k) at 25 Myr after z_re and to 2–5 % later (0.1 Myr: 1–2 %); 10 Myr leaves
  40–120 % differences during the first 100 Myr; the uncapped hydro step
  (44 Myr) lags by factors 2–4 in E(k) for 200 Myr. The ionisation
  profiles do not depend on the coupling. Figure
  `docs/figures/plane_t4_N64.png`, table `docs/data/plane/t4_N64.txt`.
* **The turn-on step** (2026-10-28; `docs/STATUS.md`, "Done 2026-10-28",
  the table; `docs/figures/plane_t4_turnon_N64.png`,
  `docs/data/plane/t4_turnon_N64.txt`): the default CFL step with c_hot
  (20.2 Myr at 64³) leaves the lag of the uncapped step (E(k) off by
  1–2 at 25–100 Myr); a turn-on step of ~1 Myr (`plane_turnon_courant =
  0.05`, the default since the user's decision of 2026-10-28) brings E(k)
  to 9–20 % of the 0.1-Myr coupling and ⟨v²⟩ at 25 Myr from 9 to 26
  (km/s)² (reference 27.1); the lag is set by the first step alone.
* **C1 at 64³** (local; `docs/figures/plane_c1_N64_ek.png`): the fronts
  cross in the first step, T_M = 17 900 K behind them, E(k) steepens from
  k^−3.6 to k^−5 with no k^−5/3 range, as C26 Fig. 3 (bottom) at Δx ≥ 0.5
  kpc/h. The 128³–512³ rungs run on klone (`tools/hyak/c1_ladder.sh`).
