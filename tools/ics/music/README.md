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
| `zoom_z3.conf` | the Z3 zoom: 40 Mpc/h base at level 9 (512³), `levelmax = 17` (8 levels, 0.305 kpc/h) over `ref_extent = 0.0078` of the box (0.3125 Mpc/h) around `ref_center`; Z1 and Z2 are `levelmax` 13 and 15 with the same seeds |

Keys that matter (`[setup]`): `levelmin` is the base grid (2^levelmin
cells), `levelmax` the finest, `ref_center`/`ref_extent` the zoom region in
box units, `padding` the cells of each level around the next finer one
(must exceed the code's `amrex_n_error_buf` and proper nesting; §4.2),
`baryons = yes` gives separate gas fields, `use_2LPT = no` is the user's
1LPT choice. `[random] seed[l]` fixes the white noise of level l: MUSIC's
hierarchical generator draws each level from its own seed and constrains
it to the coarser level, so a deeper rung with the same seeds keeps every
coarser mode (the ladder property of §3). The region for the zoom comes
from the base run: the chosen Eulerian region's particles traced back to
the ICs and their Lagrangian bounding box plus a pad (§3; `region = box`
is MUSIC's default, `ellipsoid` and `convex_hull` exist for a point file).

## The generic output (what thread D converts)

`format = generic` writes one HDF5 file with the whole hierarchy, levels
0 … `levelmax` (the levels below `levelmin` are MUSIC's internal coarse
grids): per level `level_LLL_DM_rho`, `_DM_dx/dy/dz` (the displacements),
`_DM_vx/vy/vz`, `_DM_potential`, and with baryons `level_LLL_BA_rho`,
`_BA_vx/vy/vz`; the arrays carry MUSIC's ghost zones (a level of n cells
is stored as n + 2·4), and `/header/grid_off_{x,y,z}`, `grid_len_{x,y,z}`
give each level's offset and size in its own cells. Units (from
`src/main.cc` and the Gadget plugin): the densities are overdensities δ,
the displacements are in units of the box length (the Gadget plugin
writes (cell + dx) × boxlength), the velocities are peculiar velocities
in km/s (displacement × vfact, vfact = f a H(a)/h). Checked here on the
32³ test with the CAMB table: at z = 300 the baryon δ has 0.35× the rms
of the CDM δ (the two transfer functions), the velocities are ~1.6 km/s
rms. The code's files want ρ in units of the mean total density (Ω_b/Ω_m
for the gas), velocities in grid cells per code time unit, and particles
on the lattice plus displacement (`tools/ics/make_ics.py`'s docstring
has the conventions): that mapping is the converter's job, with the T0
test (a uniform MUSIC box against `make_ics.py` at the same P(k): the
growth factor and the velocity factor must agree to the
float round-off of the FFTs).
