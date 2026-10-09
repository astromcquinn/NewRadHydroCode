# Halo mass function at production resolution

*2026-10-09, step 2 of the production-run thread.*

## Run

* **Setup:** `tools/production/make_production.py runs/hmf128 --N 128`, the
  production configuration. That is box 10 Mpc/h, Ngas = Ndm = 128,
  Zel'dovich ICs at z = 100 with seed 7, and on-the-fly FoF with b = 0.2
  and ≥ 20 particles.
* **No sources:** `tvir_min = 1e30`, so the RT never runs and this is
  gravity plus adiabatic hydro.
* **Catalogs:** `write_catalogs` writes one at every hydro step from
  z = 15. There are 10 catalogs, z = 12.6 to 5.5.
* **Cost:** 2 min on 4 threads, 0.86 GB.
* **Particle mass:** 5.10×10⁷ M☉; the DM particles carry all the matter.
  20 particles are 1.02×10⁹ M☉.

```
python3 tools/production/make_production.py runs/hmf128 --N 128
cat >> runs/hmf128/params.txt <<X
tvir_min = 1e30
write_catalogs = runs/hmf128/run/halos
restart_every = 0
movie_dz = 0
X
OMP_NUM_THREADS=4 build/radhydro runs/hmf128 runs/hmf128/run --quiet
python3 tools/analysis/hmf_compare.py runs/hmf128/run/halos --box 10 --N 128 --z 5.5 6.25 8 10 12.6
```

`apps/hmf_table` gives the analytic mass functions, from the same code as
the abundance matching (`mass_function.hpp`). They use the simulation
cosmology, a BBKS transfer function, and Tinker et al. 2008 with Δ = 200m
and its z-evolution extrapolated.

## Result

![HMF](figures/hmf_128_10mpc.png)

The table is in `figures/hmf_128_10mpc.txt`.

* **Above about 100 particles (M ≳ 5×10⁹ M☉), FoF agrees with
  Sheth–Tormen** at z ≈ 5.5–8: FoF/ST = 1.0–1.5, within Poisson errors.
  * It is 1.3–2.2× above Tinker. This is the usual FoF (b = 0.2) vs SO
    (Δ = 200m) mass offset, plus Tinker's z-evolution extrapolated beyond
    its calibration (z ≲ 2.5).
* **Below about 100 particles the counts fall short of ST.**

  | particles | FoF/ST, z = 5.5 | FoF/ST, z = 6.2 | FoF/ST, z = 8 | FoF/ST, z = 10 |
  |---|---|---|---|---|
  | 20–35 | 0.50 | 0.40 | 0.28 | 0.23 |
  | 35–60 | 0.55 | 0.48 | 0.49 | (noisy) |
  | 60–110 | 0.63 | 0.86 | 0.56 | 0.6 |

  * The Warren et al. (2006) low-N FoF correction makes this worse, since
    it lowers the masses.
  * The deficit is the PM force resolution. The PM mesh has one cell per
    mean particle spacing (Ngas = Ndm), which softens small halos.
  * Zel'dovich ICs add to it: their transients suppress rare early halos
    most, which is why the deficit grows with z.
* **No halos at z = 12.6,** where ST predicts a few dozen above 10⁹ M☉.
  The first FoF halo, and so the first source, forms at z ≈ 11.3.

## Consequences

1. **The halo counts are reliable from about 100 particles (5×10⁹ M☉).**
   Below that the simulation has 2–4× fewer halos than the analytic HMFs,
   increasingly at high z.
2. **Abundance matching against an analytic HMF** (Tinker, the user's
   choice) assigns L(M) as if every halo existed. So the box's own UV LF:
   * falls below the observed LF at the faint end (M_UV ≳ −15 at z ≈ 8),
     because of the missing low-mass halos;
   * rises above it for M ≳ 5×10⁹ M☉ by the FoF/Tinker ratio, 1.3–2.2.

   Matching instead to the **simulated** FoF mass function would make the
   box's LF equal the observed LF by construction, down to the resolution
   limit. This is a user decision, and the options are being put to the
   user. Either way f_esc is calibrated to τ_es, which absorbs an overall
   offset.
3. **The late first halos delay the onset of the RT.** This adds to the
   biases from 0.01c and the missing sub-grid sources.
4. **Possible improvements (not done):**
   * 2LPT ICs, or a higher start redshift;
   * a finer PM mesh (Ngas = 2 Ndm, the original's ratio, at more cost in
     hydro and RT);
   * Ndm > N.
