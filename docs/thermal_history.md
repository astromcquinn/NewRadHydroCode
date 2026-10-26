# The thermal history before heating: Compton coupling and the temperature floor

*2026-10-26.* The user asked whether the simulations reproduce the expected
temperature history before any preheating, Compton cooling off the CMB in
particular. **They do not.** Unheated gas follows the code's temperature
floor, an instantaneous decoupling at z ≈ 216. That is 0.76 of the RECFAST
temperature at the production start (z = 100) and 0.60 of it from z ≈ 30
on. The cause is two choices of the original (quirk **T1**,
`docs/port_notes.md` §5b):

* the initial temperature is the floor;
* neutral gas has no residual electrons, so it never couples to the CMB.

Nothing is changed yet; the fix proposed in §4 waits for the user's decision.

## 1. What the code does

* **Floor.** `calc_sim_var` sets `temp_min` = T_CMB for 1+z > x and
  T_CMB₀(1+z)²/x below, with x = 10³(Ω_b h²)^0.4 = 217.45. This is
  instantaneous decoupling at z ≈ 216, after which the gas is adiabatic.
  `Tgas_min = temp_min/temp1_unit` (neutral μ) is then used in three
  places:
  * **The ICs:** the thermal energy is set to it (`read_init_conditions`),
    so a run at z_init = 100 starts at 127.8 K.
  * **The legacy PM step:** it clamps every cell's thermal energy
    (`pm.cpp`, `particle_mesh`). The unsplit path (`particle_mesh_unsplit`)
    does not; the unsplit hydro has its own floor of **1 K**
    (`hydro_unsplit.cpp`, `Tmin = 1/temp1_unit`).
  * **The RT cell physics:** `temp_min` floors every cell it processes.
* **Compton** (`coef_cool(14)`, exact, `docs/xray_heating.md` §0) acts only
  in cells that `radiative_physics` processes (ionised or with photons). The
  ICs are fully neutral, so there are no residual electrons
  (x_e ≈ 2–3×10⁻⁴ after recombination) and no coupling.

## 2. The reference

The reference is `data/thermal/recfast_sim_cosmology.txt`: CAMB 1.5.9 with
RECFAST for the simulation cosmology (Ω_m = 0.27, Ω_b = 0.045, h = 0.70,
T_CMB = 2.725 K, Y_He = 0.24), giving x_e(z) and T_gas(z) with no heating.
HyRec is not compiled into the CAMB wheel; RECFAST is accurate to well
under a per cent for T_gas at these redshifts.

| z | T_CMB | RECFAST T_gas | x_e | code floor | floor/RECFAST |
|---|---|---|---|---|---|
| 300 | 820 | 771 | 4.1×10⁻⁴ | 820 | 1.06 |
| 216 | 591 | 516 | 3.4×10⁻⁴ | 591 | 1.14 |
| 150 | 411 | 314 | 3.0×10⁻⁴ | 286 | 0.91 |
| 100 | 275 | 168 | 2.7×10⁻⁴ | 128 | 0.76 |
| 50 | 139 | 50.9 | 2.3×10⁻⁴ | 32.6 | 0.64 |
| 30 | 84.5 | 19.9 | 2.2×10⁻⁴ | 12.1 | 0.61 |
| 20 | 57.2 | 9.35 | 2.1×10⁻⁴ | 5.53 | 0.59 |
| 10 | 30.0 | 2.61 | — | 1.52 | 0.58 |

* Below z ≈ 30, T_gas ≈ 0.021 (1+z)² K. The floor is 0.0126 (1+z)² K.
* The residual electrons keep the gas partly coupled to the CMB down to
  z ~ 100 and heat it noticeably to z ~ 40.

**One-zone check** (`build/xray_bkg --zstart 100 --T0 … --xe0 …`, X-rays
negligible): the code's own Compton coefficient, RECFAST's T and x_e at
z = 100 and adiabatic expansion reproduce RECFAST to 1% at z = 20 (9.28
against 9.35 K). The deficit splits about evenly between its two causes.
At z = 20:

* the code as it is: 5.53 K;
* the right initial temperature, no Compton: 7.21 K;
* Compton with the residual electrons, the code's initial temperature:
  7.88 K;
* both: 9.28 K.

## 3. The simulations

Setup:

* no-source runs, `tools/analysis/thermal_history.py`;
* 10 h⁻¹ Mpc, 64³ and 128³ (the 128³ cells are those of 40 h⁻¹ Mpc at
  512³);
* Zel'dovich ICs at z = 100 (`make_ics.py --seed 7`), `zdata1 = 6` so
  the RT never starts;
* both hydro schemes.

The temperature is the median over cells with 0.9 < Δ < 1.1, taken from
`e_gas` for unsplit.

![thermal history without heating](figures/thermal_history_noheat.png)

| z | RECFAST | unsplit 64³ | unsplit 128³ | legacy 64³ | legacy 128³ |
|---|---|---|---|---|---|
| 100 | 168 | 128 | 128 | 128 | 128 |
| 65 | 81.7 | 55.0 | 55.0 | 66.6 | 62.1 |
| 42 | 37.8 | 23.7 | 23.8 | 45.5 | 36.8 |
| 28 | 17.0 | 10.3 | 10.3 | 47.2 | 32.0 |
| 18 | 7.46 | 4.59 | 4.45 | 71.3 | 43.9 |
| 11 | 3.25 | 2.44 | 1.97 | 134 | 81.7 |
| 7 | 1.41 | 2.67 | 1.02 | 278 | 160 |

**Unsplit (the production scheme).** Mean-density gas follows the floor
curve to within 1–2% at both resolutions down to z ≈ 15: 0.60 of
RECFAST. Below that:

* at 64³ it warms numerically (2.7 K at z = 7);
* at 128³ it stays on the curve until the hydro's 1 K floor takes over
  at z ≈ 8. At z ≈ 6, 70% of the cells are at 1 K.

Voids cool adiabatically, about Δ^(2/3) below the mean, down to the
1 K floor.

**Legacy scheme.** The total-energy scheme heats the cold gas numerically,
the error `hydro_scheme = unsplit` was made to fix (`docs/hydro_unsplit.md`):

* mean-density gas reaches 30–50 K by z ≈ 30 and 160–280 K by z = 7, 100
  to 200× RECFAST;
* voids reach 80–480 K;
* the error falls with resolution, but slowly.

The legacy scheme cannot represent the unheated IGM at all.

**Consequences.**

* Before any heating, the production runs' cold IGM is 40% too cold at
  z ≲ 30. That changes the gas pressure, the Jeans and filtering scales
  of minihalo gas, and the gas clumping before reionization.
* For the C26-type turbulence runs (`docs/zoom_turbulence.md`), the
  unheated case is the point. Starting at z = 300 they would sit at T_CMB
  (6% too hot) until z ≈ 216, then follow the floor curve.
* Once X-ray heating (`docs/xray_heating.md`) or photo-heating takes over,
  the initial history is forgotten. In the one-zone model, with the
  default L_X/SFR, that happens from z ≈ 13.

## 4. Proposed fix (not implemented; for the user to decide)

Behind a flag whose default keeps the original:

1. **Initial temperature** T_gas(z_init) from the RECFAST table (the ICs'
   thermal energy) instead of `Tgas_min`.
2. **Compton coupling of the neutral gas through its residual electrons**,
   x_e,res(z) from the table, in every cell at every hydro step:
   T ← T_CMB + (T − T_CMB) e^{−Δt/t_C}, with
   t_C = 3 m_e c (1 + f_He + x_e)/(8 σ_T a T_CMB⁴ x_e).
   This is thermal only (the ionisation state is unchanged, so the RT is
   unaffected) and exact for any Δt. Where the RT physics already has
   electrons (ionised cells, the X-ray electrons) the existing Compton term
   applies and the residual one is the larger of the two x_e.
3. **The floor:** keep `temp_min` as the original but make the residual
   coupling the physical model. For the unsplit hydro, also lower its 1 K
   floor (with the right history the mean-density gas is 1.4 K at z = 7;
   voids are colder).

The one-zone check (§2) shows that 1 + 2 reproduce RECFAST at the mean
density. Open choices:

* the flag: `legacy_physics = 0`, which changes its stored references, or
  a new key such as `thermal_init = recfast`;
* whether the hydro floor should go below 1 K.
