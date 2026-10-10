# The 256³ production run

*Production-256 thread, 2026-10-10, on klone (`docs/cluster.md`). Branch
`production-256`. Runs under `$RUNS = /gscratch/astro/mcquinn/radhydro/runs`,
Slurm logs in `/gscratch/astro/mcquinn/radhydro/logs/prod256/`.*

## 1. Configuration

`tools/production/make_production.py $RUNS/<run> --N 256 --c-reduced 0.1
--am-params data/uv_lf/am_params_256_10mpc.txt --f-esc <f> --xi-ion 25.38
--movie-dz 0.025 --movie-fields fHI,temp,rhog --restart-every 4
--hydro-scheme unsplit`, submitted with
`sbatch --exclusive -c 192 --mem=256G --time=1-00:00:00
--export=ALL,PROB=$RUNS/<run>,THREADS=192 tools/hyak/reion_production.slurm`.

| Setting | Value |
|---|---|
| Box, resolution | 10 Mpc/h, Ngas = Ndm = 256³, Zel'dovich ICs at z = 100, seed 7 (the ICs of `reion256`) |
| Speed of light | `c_reduced = 0.1` (convergence: `docs/production_run.md` §9) |
| Rays | `nray_max_factor = 16`, `ray_rotate = 3` (B16 fix), `deterministic_parallel`, `rt_substep_carry` |
| Physics flags | `legacy_rt = legacy_physics = 0` |
| Hydro | **`hydro_scheme = unsplit`** (dual energy; `docs/hydro_unsplit.md`), `ra_max` at the original 0.1 |
| Sources | on-the-fly FoF (≥ 20 particles = 1.3×10⁸ M☉), `sfr_model = abundance` with the 256³ L(M) (`data/uv_lf/am_params_256_10mpc.txt`, `docs/hmf_validation.md`), log ξ_ion = 25.38, f_esc below |
| Outputs | fields at the 13 `zdata1` redshifts (z = 15 … 5.5), movie frames of x_HI, T and ρ_gas every Δz = 0.025, checkpoint every 4 hydro steps (`restart_keep_latest`) |
| Node | one whole `cpu-g2` node, 192 threads |

## 2. The f_esc = 0.28 run (`$RUNS/reion256_prod`, job 41706017)

f_esc = 0.28 is the 128³ calibration to τ_es = 0.054 at 0.1c
(`docs/production_run.md` §8). This run was made before the user's
decisions of the same day (`docs/PLAN.md`: round f_esc values, and
reionization ending at z = 5.5 as the only target), so it serves as the
256³ calibration point rather than the production run.

| | |
|---|---|
| Wall time | **1 h 31 min 54 s** on 192 threads (RT 4800 s, IO 452 s, hydro 165 s, FoF 24 s, PM 21 s) |
| Peak RSS | 30.8 GB |
| Hydro steps / RT substeps | 48 / 4042 (RT from z = 13.1) |
| Sources at the end | 5144 halos |
| Rays | peak 1.0×10⁸ live rays (0.37 Nray_max); the pool guard never acted |
| Disk | 139 GB (296 movie frames of 3 × 64 MB, the last checkpoint 16 GB) |

The history (`tools/analysis/reion_history.py`), against the global-model
run `reion256_c01` of `docs/cluster.md` §7 (same ICs and c, legacy hydro,
`ray_rotate = 0`) and the 128³ run with the same f_esc:

| run | weighting | z(10%) | z(25%) | z(50%) | z(75%) | z(90%) | z(99%) | z(99.9%) | τ_es |
|---|---|---|---|---|---|---|---|---|---|
| 256³, abundance-matched, f_esc = 0.28 | V | 10.34 | 8.58 | 7.63 | 7.10 | 6.78 | 6.44 | 6.30 | 0.0631 |
| | M | 10.89 | 9.11 | 7.84 | 7.19 | 6.81 | 6.44 | 6.28 | |
| 256³, global model (`reion256_c01`) | V | 12.64 | 11.48 | 10.38 | 9.66 | 9.25 | 8.83 | 8.65 | 0.0904 |
| | M | 12.89 | 11.76 | 10.56 | 9.75 | 9.29 | 8.82 | 8.64 | |
| 128³, abundance-matched, f_esc = 0.28 | V | 8.66 | 7.80 | 7.14 | 6.67 | 6.39 | 6.05 | 5.92 | 0.0541 |
| | M | 8.99 | 8.03 | 7.24 | 6.70 | 6.38 | 6.03 | 5.90 | |

| run | weighting | z = 10 | z = 9 | z = 8 | z = 7.5 | z = 7 | z = 6.5 | z = 6 |
|---|---|---|---|---|---|---|---|---|
| 256³, AM, f_esc = 0.28 | V | 0.123 | 0.200 | 0.374 | 0.556 | 0.801 | 0.982 | 1.000 |
| | M | 0.171 | 0.261 | 0.452 | 0.622 | 0.828 | 0.981 | 1.000 |
| 256³, global model | V | 0.622 | 0.965 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| 128³, AM, f_esc = 0.28 | V | 0.027 | 0.065 | 0.195 | 0.341 | 0.577 | 0.850 | 0.995 |

![history](figures/production_256_history_fesc0.28.png)

* **Against the global source model** (`reion256_c01`): the
  abundance-matched sources reionize the box 2.7 units of redshift later
  (midpoint z = 7.6 against 10.4, end z = 6.4 against 8.8) with τ_es =
  0.063 against 0.090. The global model's emissivity was never calibrated
  (`docs/cluster.md` §7 noted τ_es = 0.07–0.09 as "too early"); the
  abundance-matched one is tied to the observed UV LFs. The run also
  differs in `hydro_scheme` and `ray_rotate`, which `docs/hydro_unsplit.md`
  and `docs/port_notes.md` B16 found to leave the global history
  unchanged at 64³, so the source model is the whole difference.
* **Against 128³ at the same f_esc:** the 256³ box reionizes earlier by
  Δz = 0.5 at the midpoint and 0.4 at the end, with τ_es 0.063 against
  0.054, and starts much earlier (x_V = 0.12 against 0.03 at z = 10). Two
  causes, in order: the 256³ box resolves halos down to 1.3×10⁸ M☉
  (128³: 10⁹ M☉), and with L(M) matched to the observed LF down to the
  resolution limit those halos add their share of the emissivity, most
  of it at z ≥ 9 where they are the only halos; and the two N are
  different realisations (`docs/cluster.md` §7: the 256³ fundamental mode
  is 1.17× linear, the 128³ one 0.78×). A 128³ calibration therefore does
  not transfer to 256³; the end of reionization must be tuned at 256³.
* The cost is the same as the global-model run (1 h 32 min against
  1 h 33 min) although reionization is later: the RT cost per substep
  follows the ionised volume and the ray population, which saturates at
  0.4–0.5 Nray_max either way.
* Stills (`$RUNS/reion256_prod/run/stills/`): `still_z{12,9,7.5,6.5}.png`,
  `stills_grid.png`.

## 3. Tuning f_esc to end reionization at z = 5.5

By the 128³ slope of the end point, Δz(99 %) = 1.4 per e-fold in f_esc
(`docs/production_run.md` §8), the 256³ end at f_esc = 0.28 (z = 6.44)
moves to z = 5.5 at f_esc ≈ 0.14. Runs at the round values f_esc = 0.2
(`$RUNS/reion256_f0.2`, job 41707183) and 0.1 (`$RUNS/reion256_f0.1`, job
41707184) bracket it; see §4.
