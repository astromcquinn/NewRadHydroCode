# MUSIC on klone: zoom and C26-box initial conditions

MUSIC (Hahn & Abel 2011; the maintained code is `cosmo-sims/MUSIC`,
"MUSIC2") generates the zoom ICs of `docs/zoom_turbulence.md` §4.1 and the
uniform C26 boxes of §3b (user decision 2026-10-22). This directory holds
the setup, the transfer-table script and the configurations; the converter
to the code's files is thread D of §6.

## Install (once, on a klone login node)

```
source tools/hyak/env.sh          # gcc 13 + the radray conda env
tools/hyak/setup_music.sh         # conda: gsl hdf5 fftw camb; clone; cmake; ~10 min
```

Result: `$RADRAY_SW/music/MUSIC` (sources in `$RADRAY_SW/src/MUSIC`). The
cmake recipe was checked in the cloud container on 2026-10-22 (MUSIC
builds and runs the configs below); on klone the libraries come from the
conda env instead of the system.

Smoke test (seconds):

```
mkdir -p $RUNS/ics_test && cd $RUNS/ics_test
cp $REPO/tools/ics/music/test32.conf .
OMP_NUM_THREADS=8 $RADRAY_SW/music/MUSIC test32.conf     # -> test32.hdf5
h5ls test32.hdf5 | head
```

## Transfer functions (CAMB, at zstart)

MUSIC's `transfer = camb_file` reads a 13-column CAMB table (k/h, CDM, baryon,
photon, massless ν, massive ν, total, no-ν, total-DE, Weyl, v_CDM, v_b,
v_b − v_CDM) and wants it at the starting redshift when `baryons = yes`:

```
python tools/ics/music/camb_transfer.py camb_z300.dat --zstart 300 \
    --omega-m 0.305 --omega-b 0.048 --h 0.68 --ns 0.9667 --sigma8 0.82
```

(the defaults are these, the cosmology of Cain et al. 2026). The script
rescales A_s so that σ₈(z = 0) matches; keep `sigma_8` in the `.conf` equal
to it. The table reaches k = 5000 h/Mpc (the 0.3 kpc/h cells of Z3 need
k ≈ 10⁴ h/Mpc in the input spectrum; MUSIC extrapolates beyond the table
with the last slope, so raise `--kmax` for Z4).

## Configurations

| File | What |
|---|---|
| `test32.conf` | 0.25 Mpc/h, 32³, Eisenstein–Hu: the build check |
| `c2_uniform.conf` | the C26 box: 0.25 Mpc/h, 1024³, z = 300, CAMB, 1LPT; the C1 ladder by changing `levelmin = levelmin_TF = levelmax` to 6–9 with the same `seed[10]` |
| `zoom_z1.conf`, `zoom_z2.conf`, `zoom_z3.conf` | the zoom ladder (`music_ladder.py`, below): 40 Mpc/h, base level 9 (512³), the 0.3125 Mpc/h Z3 region at the box centre given at level 17, `levelmax` 13, 15, 17; the grids of every shared level and the seeds are the same in all three |

Keys that matter (`[setup]`): `levelmin` is the base grid (2^levelmin
cells), `levelmax` the finest, `ref_center` with `ref_dims` (cells at
levelmax) or `ref_extent` (box units) the zoom region, `padding` the cells
of each level around the next finer one (must exceed the code's
`amrex_n_error_buf` and proper nesting; §4.2), `baryons = yes` gives
separate gas fields, `use_2LPT = no` is the user's 1LPT choice.
`[random] seed[l]` fixes the white noise of level l. The region for the
zoom comes from the base run: the chosen Eulerian region's particles traced
back to the ICs and their Lagrangian bounding box plus a pad (§3;
`music_ladder.py --ids`).

## The ladder generator: `tools/ics/music_ladder.py`

```
python tools/ics/music_ladder.py Z1.conf --box 40 --base-level 9 --levelmax 13 --ladder-max 17 \
    --center 0.5 0.5 0.5 --extent-mpch 0.3125 [--seed 2026] [--pad 0] [--padding 8]
python tools/ics/music_ladder.py Z1.conf ... --ids region_ids.txt --ndm 512 --pad 0.2
python tools/ics/music_ladder.py C1_64.conf --box 0.25 --base-level 6 --uniform --seed-level 10
```

* **The region is given once, at the deepest level of the ladder**
  (`--ladder-max`), and every level's grid follows from it by MUSIC's own
  rule (a level is the next finer one halved and grown by `padding` cells
  each side, rounded outward to even cells; `mesh.hh`). Each rung is
  written with its finest grid exact (`ref_center` + `ref_dims`, dyadic
  numbers that MUSIC's floating point reproduces), so every rung has the
  same grid at every level it shares with another; the `.conf` lists the
  whole chain in its header. This is what makes the ladder property hold
  in MUSIC2 (next section). A shallower rung's finest region is the
  deepest region plus the shells: for the 0.3125 Mpc/h Z3 region, Z1's
  level 13 is 96³ cells (0.47 Mpc/h), Z2's level 15 280³.
* **Seeds:** seed[l] = seed + 1000 (l − base level) for l = base … levelmax,
  the same in every rung. `--uniform` writes one seed at `--seed-level`
  (the C1 ladder: MUSIC restricts it to each rung's level, so the rungs
  share their modes).
* **Region from particle ids:** 1-based ids of a run started from an
  `xv_dm.init` with `--ndm` particles per dimension (the id is the IC order,
  lattice cell id − 1 with i fastest); the bounding box of their lattice
  cells on the periodic circle per axis (the shortest arc), plus `--pad`
  Mpc/h on every side.
* `--cosmology planck_c26` (default) or `original`, `--transfer camb_file`
  (default, `--transfer-file camb_z300.dat`) or `eisenstein`, `--zstart`,
  `--padding`, `--no-shift`.

## The generic output

`format = generic` writes one HDF5 file with the whole hierarchy, levels
0 … `levelmax` (the levels below `levelmin` are MUSIC's internal coarse
grids): per level `level_LLL_DM_rho`, `_DM_dx/dy/dz` (the displacements),
`_DM_vx/vy/vz`, `_DM_potential`, and with baryons `level_LLL_BA_rho`,
`_BA_vx/vy/vz`, all float64 indexed [x][y][z]; the arrays carry 4 ghost
cells per side (a level of n cells is stored as n + 8). `/header/grid_off_*`
and `grid_len_*` list levelmin … levelmax only: the offset of a level in
its parent's cells relative to the parent's grid, and its size in its own
cells (absolute offset: abs_l = 2 (abs_{l−1} + off_l)), in MUSIC's shifted
frame. Units: δ for the densities, box units for the displacements, km/s
per box length for the velocities (below). The values sit at the cell
centres. At z = 300 with the CAMB table the baryon δ has 0.35× the rms of
the CDM δ.

## The converter: `tools/ics/music_convert.py`

```
python tools/ics/music_convert.py RUN.hdf5 OUTDIR --conf RUN.conf [--cosmology planck_c26]
```

(MUSIC's log, `RUN.conf_log.txt`, must sit next to the `.conf`; `--log`
otherwise; `--levels 0,1` writes only those levels.) Output, float32,
Fortran order, no record markers (`make_ics.py`'s conventions):

| File | What |
|---|---|
| `xv_dm.init` | the base level as today: N³ particles (Ndm = Ngas = N), cell centre + displacement in base DM cells, v = dx/dτ in base cells per code time unit, in lattice order (particle id n + 1 sits at cell n, i fastest) |
| `u_gas.init` | the base gas: ρ (mean total matter density units, mean Ω_b/Ω_m), ρv (base gas cells per time unit), ρv²/2 (no thermal part) |
| `u_gas.init.L<l>` | refined level l (MUSIC level levelmin + l) on the level's whole grid (padding included), the same five fields, velocities in level-0 cell units (the AMReX hierarchy's convention) |
| `xv_dm.L<l>` | mass class l: one particle per level-l cell not covered by level l + 1, in the grid's order with the covered cells skipped; base units |
| `ics_levels.txt` | per level: the grid (offset in its own cells from the box origin, modulo the box, and size), the next level's grid in this level's cells, the particle count and the class mass (mean matter density × base gas-cell volume, (1 − f_b) 8^−l; and M☉/h). The base particles under level 1 are superseded by the finer classes |
| `ics_info.txt` | the cosmology (the five keys), z_init, N, box, the growth factors used, MUSIC's shift, the checks |

The conversion (1LPT):

* **Amplitude:** every field × r_D = D_code(a)/D_MUSIC(a). MUSIC normalises
  σ₈ today and scales back with its own growth, which has no radiation
  (Ω_r = 0 in its log): at z = 300 for Planck D_MUSIC = 0.0042486 (its log)
  against the code's 0.0044993 (`pmcosmo.Growth`, all matter clustering,
  radiation in H), r_D = 1.059. The converter computes MUSIC's D from the
  exact integral (it reproduces the log's D₊) and refuses a run whose log
  disagrees or shows Ω_r ≠ 0.
* **Velocities:** MUSIC writes v = vfact × ψ_θ with vfact = f a H/h in
  km/s per Mpc/h and ψ_θ the displacement of the velocity field, in box
  units, so **the stored numbers are km/s per box length** (× boxlength in
  Mpc/h for km/s; the 32³ test's "1.6" is 0.41 km/s). The code's
  velocities are v = F r_D ψ_θ N with F = (dD/dτ)/D (the code's growing mode)
  for transfer functions without velocities (Eisenstein–Hu, where ψ_θ = ψ),
  and F = d ln a/dτ (the code's expansion rate) for `camb_file`, whose
  v_CDM, v_b columns already carry f(k) at z_start (v_CDM/δ_CDM = 0.87 at
  k = 1.4 h/Mpc in `camb_z300.dat`).
* **Gas:** ρ = f_b (1 + r_D δ_b) from `BA_rho` (`DM_rho` without baryons),
  v from `BA_v`.
* **Shift:** MUSIC re-centres the box on the zoom region (`Domain will be
  shifted by` in its log, unless `no_shift = yes`); the converter undoes it
  (the base arrays rolled, the refined grids' offsets moved).
* **Cosmology:** until the code reads it at run time (thread C), the target
  is `--cosmology original | planck_c26` or the five keys `--omegam
  --omegab --hubble0 --sigma8 --ns`; a MUSIC run made for another
  cosmology is refused (relative 1e-6), and the target goes into
  `ics_info.txt` for the run to check (params.txt needs `N = Ndm = N`,
  `box`, `z_init`).
* **Guards:** on every level MUSIC's DM velocity must be vfact × its
  displacement (relative 1e-3) without velocity transfer functions, and
  correlate above 0.99 with it with them (next section).

## What MUSIC does that matters here

* **The noise of a level depends on its grid.** MUSIC2 (unlike MUSIC 1's
  "subtract the mean of 8 and add the parent") constrains a level's white
  noise to its parent by Fourier splicing over the level's whole patch
  with a Meyer blend (`random_music_wnoise_generator.cc`). A rung that puts
  level l on another grid redraws its noise (by ~3 % rms here). Hence the
  generator's fixed grid chain; with it the shared levels' noise is
  identical (≤ 1e-15) and the gas outside the deeper rung's new level is
  bit for bit the same.
* **What sits under a finer level.** On the intermediate levels the
  displacements and velocities under the next level are the mean of its 8
  children exactly; the densities are MUSIC's Fourier restriction (the
  same mean over the region, to 1e-17, but cell by cell only 0.7
  correlated with the block mean: the restriction is band-limited, the
  block mean aliases). The base level keeps its own fields there. The
  code's average-down at the first regrid replaces those cells.
* **Corrupt velocities at levelmin = 4.** MUSIC 2.0a (967651b) wrote
  garbage (uncorrelated, NaN or overflowing) or partly corrupted (0.97–0.999
  correlated) DM and gas velocities in about half of the runs with
  levelmin = 4 (a 16³ base; same config, different outcomes with the thread
  count), and correct ones in every run with levelmin ≥ 5 (≈ 15 runs). The
  converter's guards catch it; use levelmin ≥ 5.
* **`ref_dims` with spaces.** MUSIC strips spaces from the triples with
  `remove_if` without `erase`: `ref_dims = 26, 28, 28` parses as 26, 28,
  2828. Write the triples without spaces (the generator does).
* **The white noise sits at the cell centres**; `make_ics.py` puts it at the
  integer points and evaluates at n + ½, so the same noise gives fields half
  a cell apart (T0 accounts for it).
* MUSIC keeps the Nyquist planes; `make_ics.py` zeroes them.

## Tests

* `tests/music_zoom_check.py` (fast suite, `music_zoom`; numpy only): the
  converter and the generator on `tests/data/music/zoom2.npz`, a two-level
  MUSIC zoom (4 Mpc/h, base 32³, levels 6 and 7, original cosmology,
  z = 100, Eisenstein–Hu) packed with `music_convert.py --pack
  --mantissa-bits 14` (895 KB): the generator reproduces the `.conf` and
  MUSIC's grids; a cosmology mismatch is refused; the classes tile the box
  (the masses sum to the box's DM mass to 1e-15); every class's v/ψ is the
  code's (dD/dτ)/D; the restriction on the intermediate level; the region
  from ids across the periodic boundary; `radhydro` runs from the base
  files to z = 47.
* `tests/music/t0_check.py` (T0; needs MUSIC: `cmake -DRADRAY_MUSIC=…`,
  label `music`): `make_ics.py`'s own Eisenstein–Hu no-wiggle P(k) as a
  CAMB-format table, a uniform 32³ MUSIC box from it, converted, against
  `make_ics.py --Ndm 32` run on MUSIC's white noise (`wnoise_0005.bin`).
  Mode by mode (Nyquist planes excluded, the half-cell phase applied): the
  amplitude ratio 0.99999998 (MUSIC's σ₈ integral and growth against
  pmcosmo's, after r_D), the per-mode scatter 1.5 × 10⁻⁶ for the gas
  density and 3–5 × 10⁻⁸ for the particle and gas velocities (float32);
  the displacements 2 × 10⁻⁵ (float32 positions ~N); v/ψ = (dD/dτ)/D to
  1e-7 in both. `radhydro` runs from both to z = 30 and, from the
  converter's files moved half a cell onto `make_ics.py`'s sampling, agrees
  with the `make_ics.py` run to 2.4 × 10⁻⁵ (gas density, rms) and
  1.5 × 10⁻⁶ cells (particles).
* `tests/music/ladder_check.py` (label `music`): a uniform box U and the
  rungs A (levelmax base + 1) and B (base + 2) of one ladder: the base
  noise of U, A, B identical (shift undone; 3 × 10⁻¹⁷), level base + 1's
  noise identical in A and B (≤ 1.3 × 10⁻¹⁵, the threads' FFT round-off);
  the level-1 gas δ of A and B identical outside B's level 2 and 1.4 %
  (of the rms) different under it, where B carries its finer level's
  restriction; 0.4 % averaged onto level-0 cells; the level-1 displacements
  1.8 % (the multigrid Poisson solve sees B's finer density); A's and B's
  base 0.2 % from U's outside the region (MUSIC's hierarchical
  convolution; a wrong unshift would be O(1)).

## Runs on klone (2026-10-23)

One cpu-g2 node (AMD EPYC 9654), 16 threads, `$RUNS/zoom_ics/`
(`run_music.sh`, `z1wide/run.sbatch`); MUSIC 2.0a at 967651b, the CAMB
table `camb_z300.dat`, `--cosmology planck_c26`.

| Run | Grids | MUSIC | Converter | HDF5 | Code files |
|---|---|---|---|---|---|
| C2 at 256³ (`c2_uniform.conf`, levels set to 8; `seed[10]`) | 256³ | 72 s, 17.9 GB | 6 s, 1.6 GB | 2.0 GB | 369 MB |
| Z1 (`zoom_z1.conf`: base 9, levelmax 13, the Z3 region at level 17) | 512³ + 40³, 48³, 64³, 96³ | 102 s, 6.9 GB | 100 s, 10.6 GB | 15 GB | 5.5 GB |
| Z1 of §3's table (`z1wide`: 2.5 Mpc/h at level 13, `--ladder-max 13`) | 512³ + 92³, 152³, 272³, 512³ | 229 s, 12.1 GB | 163 s, 13.4 GB | 30 GB | 12 GB |

* **C2:** the seed at level 10 makes MUSIC draw the 1024³ noise and
  restrict it to 256³ (the C1 ladder's shared modes), which sets the 18 GB;
  at 1024³ itself, estimated (not run) from the same noise plus ~12 1024³ float64 grids, 100–150 GB: a bigmem or a whole cpu-g2 node. The
  converter's checks: the DM velocity–displacement correlation ≥ 0.9999999,
  the mean gas density f_b to 3 × 10⁻¹². CAMB's velocities: the CDM at
  0.92× the code's growing mode at every k (CAMB's f_cdm ≈ 0.87 against the
  code's total-matter f = 0.944 at z = 300), the gas at 0.74× (k = 50
  h/Mpc) down to 0.06× (k = 2500 h/Mpc): the baryons have not fallen into
  the small-scale CDM wells at z = 300.
* **Z1:** no shift (the region is at the box centre); MUSIC built exactly
  the generator's grids (offsets 246, 8, 8, 8 in parent cells); classes of
  50 176, 77 824, 151 552 and 884 736 particles (4.3 × 10⁶ … 8.3 × 10³
  M☉/h) besides the 512³ base particles (3.4 × 10⁷ M☉/h); the base dominates
  the files and the time. §3's Z1 (512³ finest over 2.5 Mpc/h, 1.3 × 10⁸
  particles of 8.3 × 10³ M☉/h, plus 0.34, 1.0, 3.3 M in the shells) costs
  only twice as much: four minutes and 12 GB for MUSIC on one node.
* **Settings that mattered:** `padding = 8` (the shells: 40³ … 96³ around
  a 0.3125 Mpc/h region at level 17); `overlap = 4` and MUSIC's default
  `convolution_margin = 4` unchanged; `levelmin ≥ 5` (above); `ref_dims`
  without spaces; MUSIC's memory scales with the largest grid (the base or
  the finest level) and with `seed[]` given above levelmin (the noise is
  drawn at that level over the whole box).
