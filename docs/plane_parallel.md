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
| `photon_consume_dense` | −1 | 1: every photoionisation consumes its photon (the original's `rho_source` factor, which lets cells denser than 50× the mean ionise without consuming photons, is off); 0: the original. −1 follows `legacy_physics` (on with `legacy_physics = 0`: P11, `docs/port_notes.md` §5b, user decision 2026-10-28), so the plane runs of `make_production.py` and the tests, all with `legacy_physics = 0`, have it on; the plane-parallel problems use 1 (user decision 2026-10-29) through it. Since 2026-10-28 a general key, parsed by the drivers rather than `set_plane_param` |
| `dt_max_myr`, `rt_substeps_max` | 0, 0 | the step caps of the T4 test (below); in plane mode they act from the turn-on |
| `plane_packet_rescale` | 0 | 1: the sweep's packets in flight are rescaled to the current RT substep when dt_rt changes (below, "Packet rescaling"; `plane_transport = sweep` only) |
| `land_outputs` | 0 | 1: the hydro step before a `zdata1` output ends exactly on it (so runs with different step caps are compared at the same times), and the run ends once the last output is written (its clock can stop a few ulp above that redshift); the plane-parallel problems use 1 |
| `plane_turnon_cfl`, `plane_turnon_temp`, `plane_turnon_courant` | 1, 2e4, 0.05 | the step that starts at `plane_z_re` is at most `plane_turnon_courant` × the CFL step with the sound speed of ionised gas at `plane_turnon_temp` (K) added (below; 2026-10-28; the T4 test needs a Courant fraction ≈ 0.05 to remove the turn-on lag, the default by the user's decision); `plane_turnon_cfl = 0` = the uncapped turn-on step of before |
| `T_min`, `z_floor` | 0 (off), unset | the pre-heating option of §1: a temperature floor in K from z < `z_floor` on, in the RT cell physics and the hydro, a plain switch-on (user decision 2026-10-29; `docs/production_run.md` §10, `docs/thermal_history.md` §5). C26: 100 or 1000 K for z < 15; `make_production.py --T-min K` writes `z_floor = 15` |

`tools/production/make_production.py --plane-parallel` writes a production
`params.txt` of the mode (`--plane-ndom`, `--plane-faces`,
`--plane-gamma12`, `--plane-z-re`, `--dt-max-myr`, `--rt-substeps-max`,
`--T-min`, `--z-floor`, `--thermal-floor`): DM20's schedule
`c_reduced_schedule = 0.1 0.99 0.01` unless given, cooling on,
`thermal_floor = 0.1` (the hydro floor and the
clamp of the RT physics floor, instead of the original's 1 K: the
mean-density adiabat is below 1 K at z < 8; user decision 2026-10-29,
`docs/thermal_history.md` §4), no abundance matching. Dense gas consumes its photons through
`legacy_physics = 0` (P11; the explicit `photon_consume_dense = 1` line the
script wrote before 2026-10-28 is gone, `--photon-consume-dense 0|1`
overrides).

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

### Packet rescaling (`plane_packet_rescale`, 2026-11-05)

A packet holds the photons its face emitted in one RT substep
(N ∝ Γ_HI dt_rt at launch) and moves k cells every substep, whatever the
substep's length; a cell's rates are its absorbed photons over the current
dt_rt. When dt_rt changes, the packets launched at the old one give rates
off by dt_old/dt_new until they leave the domain, a wave that crosses every
domain at once and shows in the movies as the domain boxes sliding along
+x and +y. dt_rt changes at three places:
- **the `c_reduced_schedule` switch.** 0.1c → 0.0025c makes dt_rt about
  40× longer. In the C1 1024³ runs Γ_HI,V was 16× low at z = 6.99
  (3.5×10⁻¹⁴ against 5.7×10⁻¹³ s⁻¹) and recovered by z ≈ 6.92;
- **a hydro step cut short by a hot cell.** The substeps per step drop
  (5 → 1) and their lengths vary from step to step. In the X-ray run a
  ±11–15 % ripple in Γ crossed the domains at z ≈ 6.64;
- **the step that lands on an output.** It is short, so its Γ_HI is
  ~10× high (frames.txt at z = 6.4, 6.2), and its dim packets then cross
  the domains.

The total photons are conserved, so the time-averaged rate is right, but
the instantaneous rate is not.

`plane_packet_rescale = 1` (`plane_transport = sweep` only; serial and
backend) multiplies every packet in flight by dt_new/dt_old at the first
substep with a new dt_rt (`plane_packet_scale`, before the cast), so a
packet always holds the photons of one current substep. Lowering ĉ then
raises the photon density in flight, n_γ = F/ĉ. The change counts as
emitted in the bookkeeping (`plane_emitted`), which therefore still
closes. The default 0 is the behaviour before. The ray transport
(`plane_transport = rays`, and the halo sources' rays outside this mode)
has the same property and no rescaling.

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

**The coupling cap (`plane_couple_courant`; user decision 2026-10-30).**
After the turn-on cmax_gas holds the hot gas, so the CFL step is about
0.45 Δx/c_hot (16 Myr at 64³). The RT heating and cooling meet the hydro
only once per step, and in the C1 boxes that is long enough for a
clump to cool and contract while a neighbouring cell is drained and heated
numerically (Δ 47 → 0.4 with T 1.6 × 10⁴ → 7 × 10⁵ K, then 2 × 10⁶ K,
at 64³ without SN heating; cmax_gas 10× its value for ~30 Myr; below).
`plane_couple_courant` = C caps every step from the turn-on (a ≥ a_re) at

  dt_sim = C / c_hot  (code units: Δx = x_unit), i.e. dt = C Δx_phys / c_hot,

with c_hot of the turn-on rule (`plane_turnon_temp`). It scales with the
cell like the CFL step, so every rung of a ladder has the same coupling
as a fraction of the hot gas's cell crossing (C = 0.15: 5, 2.6, 1.3,
0.65 Myr at 64³, 128³, 256³, 512³ in the 0.25 Mpc/h box). Default 0 =
off; `make_production.py --plane-parallel` writes 0.15. It comes after
the turn-on rule (which is shorter at the turn-on step). A function of a
only, so restarts repeat it. The extra hydro, PM and summary steps are
cheap next to the RT, whose substeps are set by ĉ: at 64³ the capped run
costs 602 s against 700 s for the uncapped one, which the hot cells had
slowed down.

### Tests (`tests/CMakeLists.txt`)

| test | what |
|---|---|
| `plane_rt` (`tests/plane_check.py`, RT-only driver, uniform gas built per check) | Γ_HI normalisation to 5e-8 in the source cells and to the attenuation elsewhere, for 1, 2 and 4 cells per substep; the bins and the optical depths; slab attenuation per bin = Π(1 − τ_j) to 4e-4 and e^−τ to 2 %; one front across the whole box against the photon budget (thin gas: 0.5 %; 100× the mean density with recombinations: the column = (photons − the code's recombinations)/n_H and α_eff against case B); the light-limited front at ĉ; photon conservation (3e-5); 2³ domains ionising identically (1e-6) with no ray crossing a domain edge; 1 vs 2 threads byte-identical |
| `plane_full` (`tests/plane_full_check.py`, `radhydro`, 16³ uniform box) | the landing on z_re, no photons before it, DM20's schedule switching once at x_HI,V = 1 % (substeps 1263 → 133), the caps (`dt_max_myr = 2`, `rt_substeps_max = 50`) with the T_M(z), x_HII,V(z) of the uncapped run to 0.6 %; the turn-on step (`plane_turnon_temp = 4e4`, `plane_turnon_courant = 0.5`, as the default does not bind in this cold uniform box: the `ra_max` step is shorter): dt_sim (cmax_gas + c_hot) = 0.5, dt ≤ 0.5 Δx_phys/(v_max + c_hot) up to the expansion over the step, the steps before z_re bit for bit those of the default run, and with `plane_turnon_cfl = 0` the old step (96.8 against 37.1 Myr); `plane_couple_courant = 0.3`: dt_sim c_hot = 0.3 where it binds (most steps after z_re), at most 0.3 in all of them, the steps before z_re bit for bit, T_M(z) and x_HII,V(z) of the uncapped run to 0.5 % |
| `plane_diag` (`tests/plane_diag_check.py`, `radhydro`, 32³ Zel'dovich box, C26 cosmology) | `diag_dz = 0.2`: E(k) at the first step below each of z = 9, 8.8, …, 8 and at the zdata1 steps, none above z = 9; E(k), the mass-weighted E(k) and ⟨v²⟩ agree with `velocity_power.py` on the steps' `u_gas.bin` to 6 × 10⁻⁷; the T–Δ histogram holds N³ cells and mass N³ and its mass-weighted ⟨log T⟩ matches the fields; `movie_fields vmag` writes |v| slabs; every step's stats and the last step's fields byte-identical to the run without the diagnostics |
| `plane_xray` (`tests/plane_xray_check.py`, `radhydro`, 16³ Zel'dovich box) | `xray_heating = 1` before z_re = 9 at log L_X/SFR = 40.5: mass-weighted T 200 K at z_re; `plane_xray_substeps = 10`: 10 substeps per step before z_re, mass-weighted T and x_HII at z_re within 3 % of the light-crossing count (−1.6 %), T_V, T_M, x_HII,V at the end within 2 % |
| `plane_tcap` (`tests/plane_tcap_check.py`, `radhydro`, 16³ Zel'dovich box) | `T_max = 1e15` byte-identical to off (every step's stats and fields); `T_max = 1.2e4` binds: the hydro clamps 3664 cell-substeps, the RT temperature after z_re at most T_max (11 981 K against 84 420 K uncapped), still ionised (x_HII,V = 0.99999) |
| `plane_dm20` (nightly) | the scaled-down DM20 problem: 0.256 Mpc/h, 32³, z_re = 8, Γ_−12 = 0.3, cooling on, to z = 7: x_HII,V > 0.99, T_M of the ionised gas in 1–3 × 10⁴ K, the clumping factor in 1–10 |
| `plane_rescale` (`tests/plane_rescale_check.py`, the `plane_full` box with sweep, `dt_max_myr = 0.3`, a movie frame every dz = 5e-4) | the default dips at the 0.1c → 0.01c switch (Γ_HI,V 5.98 → 2.06 × 10⁻¹², recovering as the domains refill); `plane_packet_rescale = 1` holds Γ_HI,V within 0.23 % of its value before the switch to the end, the bookkeeping closes to 2e-7, 1 and 4 threads give identical stats and frames; `amrex_plane_rescale`: the backend's Γ_HI,V history matches the serial run to 1e-4 |
| `amrex_plane_rt`, `amrex_plane_full` | the same on the AMReX backend with 2 ranks; the full loop against the serial run (x_HII,V and T_M to 1e-4, fHI to 1e-3), and the capped turn-on step: the steps to it bit for bit those of the serial run (the later CFL steps to 4e-8) |

The bitwise suite is untouched: every new flag defaults to the original.

### On-the-fly diagnostics (`diag_dz`, 2026-11-01)

For runs without snapshots (user decision 2026-11-01: the 1024³ rung
saves movies and statistics only). `diag_dz > 0` writes, every `diag_dz`
in redshift from `diag_z0` (default the first zdata1 redshift) and at the
zdata1 outputs, at the end of the hydro step that reaches it:

* `Diag/ek.<ztag>`: E(k) of the gas velocity with `velocity_power.py`'s
  conventions (proper peculiar km/s, `diag_nbins` = 32 logarithmic bins
  from k_f to k_Nyquist in h/kpc, Σ E dk = ⟨v²⟩/2), volume-weighted and
  mass-weighted (√(ρ/⟨ρ⟩) v) columns, ⟨v²⟩ of both in the header.
* `Diag/tdelta.<ztag>`: cells and mass in bins of 0.05 in log Δ (from −2,
  120 bins) and log T (from 0, 160 bins), non-empty bins only; T from the
  RT grid, so not at the landing step of z_re.

Six FFTs per output on FFTW's threads (capped at 32, byte-identical for
every count) and per-plane sums, so the files do not depend on the thread
count; the run itself is unchanged (`plane_diag`). Needs the host grids,
`Nrt = Ngas` and `hydro_scheme = unsplit` (its lab-frame momenta). The
movie field `vmag` (|v|, km/s) has the same requirements. `data_files =
none` writes no snapshot files; with `fields_at_data1 = 0` the zdata1
outputs then hold only P(k) (`Pk/`). `make_production.py --diag-dz
--data-files none --fields-at-data1 0`; `c1_ladder.sh` SNAPSHOTS=0,
DIAG_DZ, MOVIE_DZ, CHAT_AFTER. Readers: `tools/analysis/plane_ek_panels.py`
for E(k) json; `snapshot_slices.py` turns saved snapshots into a slide show
(for the 512³ rung, which had `movie_dz = 0`).

### X-ray preheating before the turn-on (`plane_xray_substeps`, 2026-11-01)

`xray_heating = 1` (`docs/xray_heating.md`) works in the plane-parallel
mode: before `plane_z_re` the RT runs every hydro step for the uniform
X-ray background alone (no face sources), processing every cell. Its
substep count is still the light crossing at `c_reduced` (~7 700 per step
at 1024³), although there are no rays to move. `plane_xray_substeps = n`
(default 0 = the light crossings) gives those steps n substeps instead;
the steps from z_re on are unchanged. The error is first order in the
substep length, ~0.35 % in T per Myr of substep (`plane_xray`, 16³ at
log L_X/SFR = 40.5, ~200 K by z_re = 9; 30–100 Myr hydro steps):

| substeps per step | 1 | 3 | 10 | 30 | light crossings (50–1200) |
|---|---|---|---|---|---|
| mass-weighted T at z_re [K] | 167.0 | 189.0 | 197.1 | 199.4 | 200.2 |

The C1 X-ray run (1024³, hydro steps of 1–3 Myr) uses n = 2 and
log L_X/SFR = 39.9, tuned with `xray_bkg` in the C26 cosmology so the
mean-density gas reaches 100 K at z ≈ 7.1 (T_CMB crossed at z ≈ 9.5;
43 K at z = 8.6, 72 K at 7.6), with the production abundance matching
(`data/uv_lf/am_params_st_fof.txt`) for the star formation. Near-empty
cells (Δ < 10⁻³, from the 16³ Zel'dovich ICs here) reach 10⁵ K under the
X-rays and dominate T_V; use mass-weighted temperatures before z_re.

### Temperature cap (`T_max`, 2026-11-04)

The C1 1024³ main and T_min = 100 runs (no floor, or a 100 K one) grew a
few very hot cells right after z_re = 7: c_s up to 1.2–1.4 × 10⁴ km/s
(~10⁹ K), v up to ~5000 km/s, at z ≈ 6.87–6.63 on and off. The T–Δ
histograms put them at low density (Δ ≈ 0.01–5), 10–100 cells of 10⁹
(3 × 10⁻⁸ of the mass). They set the global step (0.02–0.06 Myr instead
of the 0.3 Myr coupling cap) and make the unsplit hydro take 15–30 CFL
substeps of ~33 s each inside one step (t_hydro 400–1100 s instead of
33–65 s). The T_min = 1000 run had none (c_s < 500 km/s); the X-ray run
a milder episode (c_s ~ 1400 km/s at z ≈ 6.7). At 512³ they only appeared
near z ≈ 5.1.

`T_max` (kelvin; 0 = off, default, bit-identical) clamps the thermal
energy at the end of every unsplit-hydro CFL substep (E − KE and the
dual-energy e, converted with temp2_unit, the ionised gas's mean
molecular weight, so the cap is T_max in ionised gas and ~2 T_max in
neutral gas) and the RT cell physics' temperature after the heating and
cooling. The hydro's count of clamped cell-substeps is printed at the end
(stdout and timing.txt `# cap_T`). Host unsplit hydro only (refused for
the split scheme and the AMReX hydro). The user's choice for the C1
reruns: 3 × 10⁵ K, far above the photoheated gas.

### Tools

* `tools/analysis/plane_front.py`: folded front profiles (⟨x_HI⟩, ⟨T⟩
  from the source face), the front position per output and its speed in
  units of ĉ (the light-limited diagnostic).
* `tools/analysis/velocity_power.py`: the kinetic energy spectrum E(k) of
  the gas velocity (C26 Fig. 3), with `--ref` for digitised points.
* `tools/analysis/t4_coupling.py`: the T4 test (make / run / analyse).
* `tools/hyak/c1_ladder.sh`: the C1 ladder and the T4 variants on klone.
  Since 2026-10-30 the ladder starts from C26's ICs (MUSIC + CAMB from
  z = 300, the C26 cosmology, `thermal_init = recfast` with its
  `recfast_table`; `ICS=zeldovich` for the earlier Zel'dovich ICs from
  z = 100; `docs/PLAN.md`, "Plane-parallel RT-domain mode").

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

### The C1 ladder from C26's ICs (2026-10-30; klone, `tools/hyak/c1_ladder.sh`)

Runs in `/gscratch/astro/mcquinn/radhydro/runs/c1_c26` (`c1_N<N>/run`;
the cancelled first launch with SN heating in `c1_N<N>/run_sn1`). MUSIC +
CAMB from z = 300, the C26 cosmology, `thermal_init = recfast` with its
table, z_re = 7, Γ_−12 = 0.3, 8³ domains lit from +x and +y, cooling on
everywhere, no preheating, the serial build (`sweep` transport).

**64³ checks** (16 threads):

* Before z_re the mean-density gas (median T at 0.9 < Δ < 1.1, from
  `e_gas`) follows the RECFAST table to ≤ 0.6 % from z = 50 to 15, then
  +1.7 % at z = 11 and +4–5 % at z = 8–7 (the non-linear warming of
  `docs/thermal_history.md` §4). The volume-weighted T of `Log/hydro.log`
  follows the table to 0.1–0.5 % down to z ≈ 20; below that shocked gas
  dominates it. `stats.txt` shows T_M = −nan before z_re (known).
* The turn-on step is 0.81 Myr (`plane_turnon_courant = 0.05`). The
  fronts cross the 8³ domains in that first step: x_HII,V = 0.986,
  T_V = 19 800 K, T_M = 18 500 K behind them; x_HII,V = 0.9998 at
  z = 6.9, when ĉ drops to 0.01c. T_V falls to 12 700 K at z = 6 and
  7 400 K at z = 4.2.
* E(k) against the Zel'dovich 64³ rung of 2026-10-27
  (`docs/data/plane/c1_N64_ek.json`): at z = 6.6 the new rung has 5× its
  ⟨v²⟩ (47 against 9.4 (km/s)²), because the old rung predates the short
  turn-on step and its ⟨v²⟩ only starts rising after z = 6.5. From
  z ≈ 6.2 on the two agree: ⟨v²⟩ 49 against 51.5 at z = 6, 44.5 against
  46.3 at 5.5, 40.7 against 40.9 at 5; E(k) band by band within 0.7–1.9×
  (different ICs and cosmology); the same slope between 0.25 and
  0.6 h/kpc (−4.3 at z = 6). The domain imprint is small: ≤ 1 % of ⟨v²⟩
  is on the domain harmonics (k a multiple of 2π/31.25 kpc/h), gone by
  z = 6.

**SN heating in source-free boxes (found in the first launch).** The
128³ and 256³ rungs slowed to 0.2-Myr steps from z ≈ 5.5: cmax_gas rose
from 60 to 2 500 and vmax_gas to 2 000 (code units; the DM stayed at
vmax ≈ 40). `photo_cell` adds the original's sub-grid SN heating,
eSN eSF ρ c²/t_dyn, to every cell above 10⁴ K and 50× the mean density,
with or without sources. Its (200/x)^1.5 damping above 200× is part of
`cool_suppress_overdensity`, which these problems set to 0, so the
photoheated minihalo gas exploded. Per unit mass it injects ~200× the
thermal energy of 10⁴ K gas per dynamical time. At 64³ it changes
⟨v²⟩ by only 1–7 % (the episode at z ≈ 5.4–5.0, vmax ~600, relaxes). Now off
everywhere (`sn_heating = 0`, user decision; `docs/PLAN.md`).

**Numerical hot cells, and the coupling cap.** Without SN heating one
cell next to a cooling clump (Δ ~ 1 000 at 64³) was evacuated and heated
to 2 × 10⁶ K (cmax_gas 10× for ~30 Myr). Variants at 64³, peak
cmax_gas / vmax_gas between z_re and z = 4 (baseline ≈ 28 / 13):

| variant | peak cmax / vmax |
|---|---|
| `sn_heating = 0`, the CFL steps (12–16 Myr) | 150–280 / 60 |
| the same with outputs every Δz = 0.01 (steps ≤ 2.8 Myr) | 27 / 13 |
| `dt_max_myr` = 1, 3 or 5 | 31 / 14 |
| `plane_couple_courant = 0.15` (5.4–6.4 Myr) | 31 / 14 |
| `hydro_scheme = legacy` (20 Myr steps) | 21 / 10 |
| `dual_e_self = 0` | 500 / 1 300 |
| `dual_energy = 0` | 164 / 38 |

The hot cell depends on the coupling interval of the unsplit hydro, not on
the dual-energy rules alone. E(k) at 0.05 < k < 0.4 h/kpc moves by
≤ 13 % between the uncapped run and the 1-Myr cap (≤ 8 % for 5 Myr).
Hence `plane_couple_courant = 0.15` (above).

**The ladder with `sn_heating = 0` and `plane_couple_courant = 0.15`**
(relaunched 2026-10-30 21:25; 512³ on 2026-10-31):

| N | Δx [kpc/h] | threads | wall | peak RSS | RT substeps | steps |
|---|---|---|---|---|---|---|
| 64 | 3.9 | 16 | 10.0 min | 0.29 GB | 2 844 | 179 |
| 128 | 1.95 | 32 | 2 h 39 | 2.1 GB | 5 424 | 366 |
| 256 | 0.98 | 96 | 5 h 40 | 15.7 GB | 10 625 | 938 |
| 512 (to z = 5) | 0.49 | 192 | 22 h 13 | 110 GB | 12 562 | 1 006 |

64³–256³ ran to z = 4, 512³ was stopped at z = 5. To z = 5 the RT took
0.20, 1.7, 1.0 and 4.5 s per substep (16, 32, 96, 192 threads) and 96 %
(128³) to 72 % (512³) of the time. At 256³ the hydro takes 40 % of the time
to z = 4 (the late hot cells below); 512³ spends 14 % in the hydro, 6 % in
the PM. After z_re the 512³ steps are 0.71–0.79 Myr (the coupling cap),
with 9 RT substeps each.

**1024³ (extrapolated, not run).** The cost per RT substep scales as N³
at a fixed thread count, and the cap halves dt, so 1024³ needs ~16× the
512³ run: ~350 h on one 192-core node. Memory ~0.9 TB fits one cpu-g2
node (1.5 TB), so the serial build with `thermal_init = recfast` works
(recfast is refused with `amrex_state = 1`). Disk ~0.6–0.8 TB (512³:
~90 GB). The RT after overlap dominates: ~21 substeps per step at
ĉ = 0.01c; a lower ĉ after overlap would cut it most.

**Thermal history** (`docs/data/plane/c1_c26_history.txt`, every step
from z_re). x_HII,V = 0.9998 at z = 6.9 and > 0.99995 from z = 6.2 in
every rung (`plane_front.png`: the fronts sit at the domain edge from the
first step). T_V and T_M converge:

| z | T_V 64 / 128 / 256 / 512 [K] | T_M 64 / 128 / 256 / 512 [K] |
|---|---|---|
| 6.9 | 19 135 / 19 081 / 19 036 / 18 997 | 18 837 / 17 971 / 16 932 / 16 270 |
| 6.6 | 16 661 / 16 533 / 16 521 / 16 623 | 15 933 / 15 448 / 15 257 / 15 180 |
| 6.2 | 13 989 / 14 031 / 14 195 / 14 370 | 14 004 / 14 038 / 14 122 / 14 170 |
| 6.0 | 12 927 / 13 050 / 13 243 / 13 415 | 13 191 / 13 270 / 13 474 / 13 513 |
| 5.5 | 10 822 / 11 031 / 11 232 / 11 373 | 11 427 / 11 470 / 11 618 / 11 730 |
| 5.0 | 9 243 / 9 454 / 9 631 / 9 739 | 9 966 / 10 028 / 10 163 / 10 257 |

From z = 6.2 T rises by 1–1.5 % per doubling, not yet converging (probably shock
heating of the better-resolved small-scale flows). At z = 6.9 T_M falls with N (−14 % from 64³ to 512³),
probably because more of the mass sits in resolved dense gas that cools fast.

**E(k) ladder** (`docs/figures/plane_c1_c26_ek.png`; by redshift with k^5/3 E(k): `docs/figures/plane_c1_c26_ek_panels.png`, `tools/analysis/plane_ek_panels.py`;
`docs/data/plane/c1_c26_ek.{json,txt}`). ⟨v²⟩ [(km/s)²]:

| z | 64 | 128 | 256 | 512 |
|---|---|---|---|---|
| 6.9 | 22.6 | 33.1 | 45.1 | 54.2 |
| 6.6 | 53.2 | 69.7 | 79.5 | 85.2 |
| 6.2 | 56.4 | 68.7 | 76.7 | 82.0 |
| 6.0 | 54.4 | 65.7 | 72.9 | 77.9 |
| 5.5 | 48.7 | 57.5 | 63.5 | 68.0 |
| 5.0 | 43.6 | 50.7 | 55.9 | 59.9 |

(Corrected 2026-11-01: `velocity_power.py` had converted the velocities
with the original cosmology, Ω_m = 0.27, instead of the run's 0.305; it now
reads `omegam` and `hubble0` from params.txt. The fix multiplies v by
√(0.305/0.27) = 1.063 and E(k) and ⟨v²⟩ by 1.13; slopes and k are
unchanged. The ⟨v²⟩ values quoted in "64³ checks" and in the sink tests
predate the fix and are 13 % low for these runs; the old Zel'dovich rung
used the original cosmology and is unaffected, so at 64³ the new rung's
⟨v²⟩ is ~13 % above it from z ≈ 6.2, not equal.)

The increments per doubling shrink by ~0.7 (11.3, 7.2, 5.0 at z = 6), so
the converged ⟨v²⟩ is ~10 (km/s)² above 512³. Local slopes d ln E / d ln k
(least squares in each band, z ≤ 6.2):

| N | 0.25–0.5 h/kpc | 0.5–1 | 1–2 | 2–4 |
|---|---|---|---|---|
| 64 | −3.7 … −4.3 | −4.7 … −5.0 | — | — |
| 128 | −2.4 … −2.75 | −3.8 … −4.0 | −5.1 | — |
| 256 | −1.7 … −2.0 | −2.5 … −2.7 | −3.8 … −4.1 | −5.1 |
| 512 | −1.5 … −1.7 | −1.85 … −1.95 | −2.6 … −2.8 | −3.8 … −4.15 |

**At 512³ a k^−5/3 range appears**: slopes −1.5 to −1.95 from 0.25 to
1 h/kpc (wavelengths 25–6 kpc/h, 51–13 cells), for z = 6.2 to 5. Each
doubling moves the steep (−4 to −5) dissipation range up by a factor 2 in
k, so the cutoff is numerical, at ~10 cells per wavelength. The inertial range is a factor
~4 in k at 512³, and its slope at 0.25–0.5 h/kpc is still changing
(−1.8 → −1.6 from 256³). Below 0.25 h/kpc E(k) is the domain-scale
forcing (flat or rising at z = 6.9, then ≈ −1 … −3 with the domain
harmonics, 2π/31.25 kpc/h).

512³ developed the hot cells before z = 5: a spike at z ≈ 5.5 (cmax_gas
290 → 560) and from z ≈ 5.1 an episode with cmax_gas up to 4 400,
vmax_gas 880 and 0.1–0.35 Myr steps (144 steps with cmax_gas > 500),
still running at z = 5. Its z = 5 E(k) shows no break at high k
(slope −4.15 at 2–4 h/kpc against −4.0 at z = 5.5), and ⟨v²⟩ follows the
other rungs, so the episode stays in a few cells; ending at z = 5 does
not avoid it at 512³, though.

**Late hot cells, and the serial sink.** With the cap the 64³ rung stays
clean to z = 4, but 128³ has an episode from z ≈ 4.15 (cmax_gas 52 →
2 200, vmax_gas 25 → 420, ~40 small steps) and 256³ from z ≈ 4.9 (vmax_gas
up to 9 400). Each time a near-empty cell (Δ ≈ 0.01, T ~ 10⁸ K) sits next
to a clump that cools without limit (Δ up to 10⁴ at 128³), the same
pathology as above at higher density. The ladder therefore ends at z = 5
for now (user decision, option 2), and the sink was tried alongside (option 3).

`sf_overdensity > 0` now works on the serial unsplit hydro too (before,
only on `hydro_backend = amrex`; `src/hydro_unsplit.cpp`, `serial_sink`).
After each hydro step every cell above `sf_overdensity` × the mean gas
density loses the fraction f = min(1, `sf_efficiency` dt / t_ff) of its
mass, momentum, internal and kinetic energy and species, with
t_ff = (3π/4) / √(a Δ) in code units. The removed mass is not turned into
particles (no gravity from it; the PM still sees the DM and the remaining
gas). Sums are per k-plane, so the result does not depend on the thread
count. The removed total goes to `timing.txt` (`# serial_sink events …
mass … gas_fraction`). The default (0) is byte-identical to before.

128³ to z = 4 (runs `sink_N128_sf<threshold>`, `sf_efficiency` = 1):

| threshold | peak cmax / vmax for z < 5.2 | gas removed by z = 4 | T_V(z = 4) | wall |
|---|---|---|---|---|
| none | 2 228 / 420 | — | 7 207 K | 2 h 39 |
| 3 000× | 74 / 33 | 2.0 % (285 events) | 7 197 K | 2 h 28 |
| 1 000× | 58 / 29 | 3.1 % (1 587 events) | 7 189 K | 2 h 33 |

Both thresholds remove the episode, as well as the smaller bumps at
z ≈ 5.1 and 4.7 in the run without a sink. x_HII,V and T_V change by
≤ 0.25 %. E(k) at k < 2 h/kpc changes by ≤ 0.1 % (3 000×) or 1.2 % (1 000×) at z = 6, ≤ 1.5 % (3 000×)
or 2.7 % (1 000×) at z = 5 and ≤ 5 % at z = 4.5. At z = 4 the 7 % difference
is the run without a sink's hot-cell episode
(`runs/c1_c26/post_sink`). 3 000× is the less invasive choice; whether the
threshold should scale with resolution (Δ_max grows with N) is open;
256³ and 512³ with the sink are not run yet.
