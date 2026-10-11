# Rays, M1 and OTVET at 40 Mpc/h, 512³: the HI photoionization rate

*Run on klone 2026-10-11 (branch `rt-methods-512`; one whole `cpu-g2` node,
192 threads, per run).* Three reionization runs that differ only in
`rt_method` (rays, M1, OTVET), with the same ICs, sources, speed of light
and recombination treatment, for movies and for the comparison of the HI
photoionization rate Γ_HI after overlap. A fourth run, rays with case B
(on-the-spot) recombination and no recombination radiation, separates the
effect of the recombination radiation from that of the RT method. The 256³
comparison at 10 Mpc/h is `docs/rt_moments.md` §6.

## 1. Configuration

| | |
|---|---|
| Box, resolution | 40 Mpc/h, N = Ndm = 512 (cell 78 kpc/h, m_p = 5.1×10⁷ M☉: the mass and cell resolution of the 10 Mpc/h 128³ box) |
| ICs | Zel'dovich, z = 100, seed 7: the ICs of `reion512_40mpc_f0.1` (identical parameters, linked) |
| Flags | `legacy_rt = legacy_physics = 0`, `deterministic_parallel`, `rt_substep_carry`, `hydro_scheme = unsplit`, `ray_rotate = 3` (used by the rays only) |
| Sources | `sfr_model = abundance` with `data/uv_lf/am_params_128_10mpc.txt`, log ξ_ion = 25.38, constant f_esc = 0.23 (the ṅ_ion calibration of 2026-10-11, `docs/PLAN.md`) |
| Light | `c_reduced = 0.2`, constant for the whole run; `redshift_per_distance = 1` |
| Recombination | **case A with the recombination radiation transported** (`recomb_case = A`, `recomb_rad = 1`, `recomb_rad_block = 4` for the rays; `docs/recomb_radiation.md`), the production default since this thread (user, 2026-10-11). Fourth run: case B, no recombination radiation |
| Ray pool | `nray_max_factor = 15` (16 × 512³ overflows the int32 ray links) |
| Redshifts | z = 100 → 5; RT from z = 15 (the first `zdata1`; the first halos above 20 particles appear at z ≈ 13.2) |
| Movie | every Δz = 0.025 (≈ 400 frames), slabs of 8 cells normal to x, y, z (`movie_slab = 8`) of x_HI, T, ρ_gas and Γ_HI (`gammaHI`); `frames.txt` with the whole-box Γ_HI statistics (§2) |
| Snapshots | z = 15, 7, 6.5, 6, 5.5, 5: x_HI, x_HeI, x_HeII, T, ρ_gas, Γ_HI and `Gas/u_gas` (`data_files = gas`) |
| Checkpoints | every 4 hydro steps, the newest two kept (`restart_keep = 2`) |

```
R=/gscratch/astro/$USER/radhydro/runs
python3 tools/production/make_production.py $R/reion512_rays --N 512 --box 40 --c-reduced 0.2 \
    --nray-max-factor 15 --am-params data/uv_lf/am_params_128_10mpc.txt --f-esc 0.23 --xi-ion 25.38 \
    --movie-dz 0.025 --movie-slab 8 --movie-fields fHI,temp,rhog,gammaHI \
    --fields fHI,fHeI,fHeII,temp,rhog,gammaHI --zdata1 "15 7 6.5 6 5.5 5" --restart-every 4 \
    --restart-keep 2 --hydro-scheme unsplit --redshift-per-distance 1 --data-files gas
# the same with --recomb-case B --recomb-rad 0 -> $R/reion512_rays_caseB
SBATCH_ARGS="--exclusive --cpus-per-task=192 --mem=300G --time=20:00:00" \
    tools/hyak/rt_methods_runs.sh $R/reion512_rays          # -> reion512_rays_{m1,otvet}
sbatch --exclusive --cpus-per-task=192 --mem=400G --time=1-20:00:00 \
    --export=ALL,PROB=$R/reion512_rays,THREADS=192,POST=sbatch tools/hyak/reion_production.slurm
# (and the same for reion512_rays_caseB)
SETS=reion512_rays:rt_methods_512 EXTRA=reion512_rays_caseB:rays_B IGM=1 \
    sbatch --time=8:00:00 tools/hyak/rt_methods_post.slurm
```

## 2. What was added for this comparison

- **Whole-box Γ_HI in `movie/frames.txt`.** Every frame line ends with
  `GammaHI_V GammaHI_ion_med GammaHI_ion_mean fV_GammaHI_gt_1e-14`: the
  volume mean of Γ_HI, its median and mean over the cells with
  x_HII > 0.5, and the volume fraction with Γ_HI > 10⁻¹⁴ s⁻¹. They are
  taken from `rt_grid%gamma(1)` of the last RT substep before the frame
  (the same instantaneous value as the snapshots, `docs/igm_diagnostics.md`
  §2), so the Γ history needs no cube. Serial sums, so they do not depend
  on the thread count. Old six-column files, and files that change at a
  restart, still read (`reion.py` `frames_table`, `read_gamma`).
- **`gammaHI`** for `movie_fields` and `fields`: Γ_HI alone, one float per
  cell (the user asked to store only the HI rate: `gamma` has three
  components and 3× the disk).
- **`restart_keep = n`**, **`data_files = gas`**, and
  `make_production.py --fields --restart-keep --data-files --recomb-case
  --recomb-rad --recomb-rad-block` (`docs/production_run.md`).
- **Tools.** `compare_movie.py`: a log Γ_HI row (fixed range
  10⁻¹⁴–10⁻¹¹ s⁻¹) and a panel with the median Γ_HI of the ionised cells;
  `stills.py`: a Γ_HI column; `rt_methods_runs.py`: the Γ history (the four
  columns against z and their ratios to the rays) in `history.png`,
  `maps_gamma_zXX.png` and `gamma_pdf.png` (log Γ_HI in ionised cells) from
  the snapshots; `igm_stats.py` reads `gammaHI.bin`;
  `rt_methods_post.slurm` takes `SETS`, `EXTRA` and `IGM=1`.
- **Tests.** `movie_slab` now also writes `gammaHI` frames and checks the
  four columns of every frame against values recomputed from the full
  cubes. Fast suite on klone: 102/112, the 10 known stored-reference
  failures (`docs/cluster.md` §3).
- **Tried first on** the 256³ runs of `docs/rt_moments.md` §6 (snapshot
  Γ maps and PDFs from their `gamma.bin`) and on a 128³, 10 Mpc/h set with
  this exact configuration (`$RUNS/rtm128_test{,_m1,_otvet,_caseB}`,
  comparison in `$RUNS/rtm128_test_cmp/`).
