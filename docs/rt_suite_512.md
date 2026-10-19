# The uniform 512³ suite: rays, M1 and OTVET at one ĉ, f_esc and recombination treatment

*Run on klone from 2026-10-18 (branch `rt-suite-512`; one whole `cpu-g2`
node, 192 threads, per run; `tools/hyak/rt_suite512.sh`).* The rerun of
every RT method with the decisions of `docs/PLAN.md` ("Unified RT-methods
suite", 2026-10-17/19) for the paper on the moment methods against the
rays: the same ICs, sources, reduced speed of light (the schedule on the
ionised fraction), f_esc(z), recombination treatment (case B) and outputs
in every run, plus one ray run with case A and the recombination radiation
transported. The first 512³ set (constant 0.2c, constant f_esc = 0.23,
case A + RR, 512³ particles) is `docs/rt_methods_512.md`; the closure runs
of that set are `docs/otvet_closure.md`. Both documents and their figures
stay as they are; the runs of those names were overwritten by this set
(user decision, 2026-10-18).

## 1. Configuration

| | |
|---|---|
| Box, resolution | 40 Mpc/h, N = 512 gas cells (78 kpc/h), **Ndm = 1024** particles (m_p = 6.4×10⁶ M☉) on the **2048³ PM mesh** (`pm_refine = 4`, the uniform-grid refined mesh of `docs/pm_refine.md`; no AMR) |
| ICs | Zel'dovich, z = 100, seed 7 with the 512³ large-scale phases: `$RUNS/hmf512_40mpc_r2/ics` (linked) |
| Flags | `legacy_rt = legacy_physics = 0` (`photon_conserve` with it), `deterministic_parallel`, `rt_substep_carry`, `hydro_scheme = unsplit`, `ray_rotate = 3` (rays), `moment_inject = 1` (moments), `redshift_per_distance = 1` |
| Sources | `sfr_model = abundance` with **`data/uv_lf/am_params_st_fof.txt`** (the Sheth–Tormen refit to the well-resolved FoF halos of the 40 Mpc/h 1024³ run, `docs/pm_refine.md` §6; the sources of `fesc512_b2.5_{uniform,amr}`), `fof_min_npart = 20` (1.3×10⁸ M☉ with 1024³ particles), log ξ_ion = 25.38, **f_esc(z) = 0.135 ((1+z)/7)^2.5** (`f_esc = 0.135`, `f_esc_beta = 2.5`, `f_esc_z0 = 6`; no bump, user decision 2026-10-19) |
| Light | **the schedule** `c_reduced_schedule = 0.1 0.25 0.2 0.5 0.3 0.75 0.5`: 0.1c below x_HII,V = 0.25, 0.2c to 0.5, 0.3c to 0.75, 0.5c after (each hydro step's ĉ from the previous step's x_HII,V) |
| Recombination | **case B, no recombination radiation** (`recomb_case = B`, `recomb_rad = 0`) for the base rays and every moment method; the `caseA` ray run: `recomb_case = A`, `recomb_rad = 1`, `recomb_rad_block = 4`, `recomb_rad_cell_onspot = 1` |
| Cooling, sink | `cool_suppress_overdensity = 200` (the original), no star-formation sink |
| Ray pool | `nray_max_factor = 15` |
| Redshifts | z = 100 → 5; RT from z = 15 (the first `zdata1`) |
| Movie | every Δz = 0.025 (≈ 400 frames), slabs of 8 cells (`movie_slab = 8`) of x_HI, T, ρ_gas, Γ_HI (`gammaHI`, the last substep's) and **`gammaHI_avg`** (the running substep average of the step) |
| Snapshots | z = 15, 7, 6.5, 6.25, 6, 5.75, 5.5, 5.25, 5: x_HI, x_HeI, x_HeII, T, ρ_gas, `gammaHI`, `gammaHI_avg` and `Gas/u_gas`; **no dark matter** (`data_files = gas`; the particles are only in the checkpoints) |
| Checkpoints | every 10 hydro steps, one kept (`restart_keep = 1`; the writes are atomic) |
| Job | one `cpu-g2` node each, `--exclusive --cpus-per-task=192 --mem=1000G --time=4-00:00:00`, `reion_production.slurm` (resumes from the checkpoint by itself, writes `DONE`) |

The runs (`$RUNS/<name>`, the names of the 2026-10-11 set):

| run | method | `rt_method` and closure keys | recombination |
|---|---|---|---|
| `reion512_rays` | rays (the base) | `rays` | B |
| `reion512_rays_caseA` | rays | `rays` | A + RR, `recomb_rad_cell_onspot = 1` |
| `reion512_rays_m1` | M1 | `m1` | B |
| `reion512_rays_otvet` | OTVET, thin tensor | `otvet` | B |
| `reion512_rays_otvet_split` | OTVET, direct + diffuse split | `otvet`, `otvet_split = 1` | B |
| `reion512_rays_otvet_screen` | OTVET, screened kernel | `otvet`, `otvet_closure = screen` | B |
| `reion512_rays_otvet_cone` | OTVET, HI-attenuated cone tensor | `otvet`, `otvet_closure = cone` | B |

```
source tools/hyak/env.sh                 # RUNS=/gscratch/astro/$USER/radhydro/runs
sbatch tools/hyak/build_hyak.sh          # then ctest -R "gamma_avg|c_schedule|movie_slab"
export AM_PARAMS=data/uv_lf/am_params_st_fof.txt   # the script's default is the 128-box table
OVERWRITE=1 tools/hyak/rt_suite512.sh probe      # 3 steps at the 2048^3 mesh: memory, PM time
OVERWRITE=1 tools/hyak/rt_suite512.sh base       # reion512_rays
OVERWRITE=1 tools/hyak/rt_suite512.sh siblings   # caseA, m1, otvet, split, screen, cone
tools/hyak/rt_suite512.sh post                   # rt_methods_post.slurm -> $RUNS/rt_methods_512/
```

## 1b. The first launch, and why the table changed

The first launch (2026-10-18 20:55, the 128-box table
`am_params_128_10mpc.txt` of the 512³ f_esc calibration, decision (1) of
2026-10-18) ionised the box by z ≈ 9: x_V = 0.29 at z = 12.1, 0.83 at
z = 10.8, 0.99999 at z = 9.2. With 1024³ particles the 20-particle FoF
floor is 1.3×10⁸ M☉, 8× below the floor of that calibration, and the
fitted L(M) (L ∝ M^1.19 below the break, plus the +0.4–0.9 dex high-z
normalisation) applied to those halos gave 20–60× the calibration run's
source luminosity density at z = 11–13 (the table in `docs/STATUS.md`,
"STOPPED 2026-10-19"). The seven runs were cancelled after 9–11.5 h. The
user chose the Sheth–Tormen table (2026-10-19), the source population of
`fesc512_b2.5_uniform` and `fesc512_b2.5_amr`, so those two runs (0.1c,
the sink, cooling suppression 0, `ray_deposit = exposure`, the AMReX
backend) are the reference for the base's early history; f6 = 0.135 and
β = 2.5 were kept without the precheck.

## 2. The probe: memory and the PM at the 2048³ mesh

Three hydro steps of the base from z = 100 (no RT yet) on one node
(`$BASE/probe`, job 42013527, 13.7 min):

| | |
|---|---|
| peak resident memory | **311 GB** (`/usr/bin/time -v`: 326178764 kB) |
| PM per step | 235, 225, 220 s (`timing.txt`; the hydro 13, 7.5, 6.4 s; summaries 5.5 s) |

The estimate of the pre-launch checklist was ≈ 275 GB of PM host arrays
at the 2048³ mesh; the measured peak is that plus the gas fields. With the
RT's 150–250 GB the runs were expected near 500–560 GB, within the 1000 GB
requested. The PM adds ≈ 4 min per hydro step, ≈ 3.5 h per run over ≈ 55
steps, against ≈ 1 h at the 512³ mesh of the first set.

## 3. Cost

*(filled in when the runs are done)*

## 4. Reionization history

*(filled in when the runs are done)*

## 5. Γ_HI with the substep average

*(filled in when the runs are done)*

## 6. The IGM at z = 5–6

*(filled in when the runs are done)*

## 7. Movie

*(filled in when the runs are done)*

## 8. Conclusions and differences from the first set

*(filled in when the runs are done)*
