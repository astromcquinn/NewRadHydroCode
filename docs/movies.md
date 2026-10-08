# Reionisation movies

`tools/movies/make_movie.py` turns an `apps/radhydro` run directory into an
MP4 and a few still PNGs. It works on a run in progress, on any N, and on any
box size.

![test movie still](movies/reion64_test_z7.17.png)

Test movie (64³, 10 Mpc/h, 0.9 MB): [`docs/movies/reion64_test.mp4`](movies/reion64_test.mp4).

## Quick start

```
# set up the run with movie frames of x_HI, T and rho_gas every dz = 0.05
python3 tools/production/make_production.py RUN --N 128
OMP_NUM_THREADS=4 build/radhydro RUN RUN/run --quiet

# movie + stills (RUN/run/movie/cosmo.mp4, RUN/run/movie/cosmo_stills/)
python3 tools/movies/make_movie.py RUN/run --jobs 4
```

## What it shows

**Style `cosmo` (the default):**

- **Slices** through the box, or **thin projections** (`--thick n` averages n
  cells), of three fields:
  - x_HII (`cividis`, linear 0–1), or log₁₀ x_HI (−5 to 0) with `--xhi-log`.
    After overlap the linear x_HII panel is uniform; the log x_HI view still
    shows structure there.
  - log₁₀ T (`inferno`, 2.5–4.7).
  - log₁₀ ρ_gas/⟨ρ_gas⟩ (`viridis`, −1 to 2).
- The colour maps are perceptually uniform, and the scales are fixed for the
  whole movie. Override a scale with `--range-xHII`, `--range-xHI`,
  `--range-temp` or `--range-rhog LO HI` (log₁₀ for T and ρ).
- **Time panel:**
  - x_HII,V and x_HII,M against z, drawn up to the current frame, with a
    marker at the current z.
  - On the right axis, τ_es accumulated since the first frame, and the total
    τ_es if the IGM were fully ionised below the current z. This total is the
    accumulated part plus the integral from 0 to z, with HeIII below z = 3.
  - A band for Planck 2018, τ = 0.054 ± 0.007. **This value is quoted from
    memory**, not checked against the paper here. It is a visual guide only.
- **Labels:** redshift; cosmic time (flat ΛCDM with the original's Ω_m = 0.27,
  h = 0.7 and Ω_r); box size and N; the slice axis and position; the current
  x_HII,V, x_HII,M and both τ_es values.

**Style `planar`** (`--style planar --ray-axis x [--source-face low|high]`),
for plane-parallel runs (`docs/plane_parallel.md`, not implemented yet):

- Slices of x_HII and T in a plane that contains the ray direction. The ray
  direction runs horizontally.
- The I-front position against cosmic time. The front is where the x_HII
  averaged over planes perpendicular to the rays first drops below 0.5,
  counting from the source face. It is linearly interpolated between cells.
- The current proper front speed, in km/s and as a fraction of c.

The front finder is checked on a synthetic planar front (see Tests) but not
yet on a real plane-parallel run. On the cosmological test run the front
stays at 0, as it should.

## Inputs

| What | Where | Notes |
|---|---|---|
| frames | `movie/frames.txt`, `movie/frame_NNNNN_<field>.bin` | Written every `movie_dz` from inside the RT substeps (`docs/production_run.md` §4). Each field is float32 N³ in Fortran order. `rhog` is in units of the mean. |
| fields missing from the frames | `step_NNNN/<field>.bin` | If `movie_fields` lacks one of fHI, temp or rhog, the tool uses the regular field output nearest in redshift and labels the panel "step n, z = …". Field outputs come only at the zdata1 steps, so prefer putting all three fields in the movie frames. |
| no movie frames | `step_NNNN/{fHI,temp,rhog}.bin` | Each step with fields becomes one frame. |
| x_HII,V, x_HII,M | `frames.txt` (or `stats.txt`) | Sums over the whole box, so they don't depend on the slice. |
| box | `params.txt` in RUNDIR or RUNDIR/.. | Or `--box`. |

**Frames before the RT starts are skipped.** The C++ writes movie frames from
`movie_z0` (by default the first zdata1, z = 15). Before the first step with
a source, though, the RT grid is not set up yet: fHI = rhog = 0, and
frames.txt has x_HII,V = 1 and x_HII,M = nan. The tool drops these frames and
reports how many. In the 64³ test, 99 of 203 frames were dropped: the first
source forms in the step ending at z = 9.95. For a run whose sources form
late, set `movie_z0` to about the first-source redshift to save disk. A C++
change to skip these frames would be easy: they are not legacy output, and
`restart_smoke_det_movie` compares a restarted run with an uninterrupted
one, so it would still pass.

## Options

| Option | Default | Meaning |
|---|---|---|
| `-o OUT.mp4` | `RUNDIR/movie/<style>.mp4` | Stills go to `OUT_stills/` (or `--still-dir`). |
| `--axis x\|y\|z` | `z` | Slice normal. `z` reads contiguous slabs and is the fastest for big N. |
| `--pos P` | box centre | Slice position along `--axis`, in Mpc/h. |
| `--thick n` | 1 | Average over n cells along `--axis`: a thin projection. |
| `--zoom H0,H1,V0,V1` | whole box | Region of the slice plane to show, in Mpc/h. |
| `--zmax`, `--zmin`, `--every n` | all | Frame selection. τ_es still accumulates from the first frame of the run. |
| `--xhi-log` | off | Show log₁₀ x_HI instead of x_HII. |
| `--fps`, `--crf`, `--size W H`, `--dpi` | 12, 20, 1920×1080, 120 | Video settings: x264, yuv420p, faststart. |
| `--jobs n` | 1 | Number of processes that render frames, in order. |
| `--stills` | `auto` | `auto` takes the frames where x_HII,V first reaches 0.1, 0.5 and 0.9, plus the last frame. Also accepts a list of frame indices, or `none`. |
| `--no-mp4` | | Render only the stills. |

## Memory and speed (for 128³ and 192³ on Hyak)

- **Bounded memory.** The tool reads frames one at a time and memory-maps each
  field file, reading only the slab of the slice (`--thick` cells). Rendered
  frames are piped as raw RGB into ffmpeg; no frame directory is written.
- **Memory doesn't grow with N or the number of frames.** Each render process
  holds one figure, about 150–200 MB of Python and matplotlib. The planar
  style reads one whole x_HII cube per frame for the front: 28 MB at 192³,
  as float64 57 MB.
- **Speed.** The 64³ test took 16 s for 104 frames with `--jobs 4`, about
  0.15 s per frame per process. Rendering dominates, so at 128³ or 192³ a
  190-frame movie should take well under a minute on 4 cores. Slices along
  `x` or `y` read strided pages, which is slower on a network file system;
  copy the frames to local scratch or use `--axis z`.
- **ffmpeg.** The tool uses `ffmpeg` from `PATH`, or the static binary of the
  `imageio-ffmpeg` package. Without root, run
  `pip install --user imageio-ffmpeg` (on Hyak, or try `module load ffmpeg`).
  `--no-mp4` needs no ffmpeg.
- **Disk for the frames.** Each frame takes 4N³ bytes per field. Writing
  fHI, temp and rhog every Δz = 0.05 from z = 15 (190 frames) needs
  0.6 GB at 64³, 4.8 GB at 128³ and 16 GB at 192³. At 192³, use
  `movie_dz = 0.1`, or `movie_fields = fHI,temp` (rhog then comes from the
  zdata1 field outputs), or set `movie_z0` close to the first source.

## The test movie (64³)

How the test movie was made:

- **Setup:** `tools/production/make_production.py RUN --N 64`. That gives a
  10 Mpc/h box, Ngas = Ndm = 64, Zel'dovich ICs at z = 100 (seed 7),
  on-the-fly FoF sources, `c_reduced = 0.01`, `nray_max_factor = 16`,
  `legacy_rt = legacy_physics = 0`, `deterministic_parallel`,
  `rt_substep_carry`, and `movie_dz = 0.05` with fHI, temp and rhog.
- **Run:** 4 threads, z = 100 → 4.87, 27 hydro steps, 90 s, peak RSS 384 MB.
- **Movie:** `make_movie.py RUN/run --jobs 4` gave 104 frames from z = 9.95
  to 4.87 (after dropping the 99 pre-RT frames).

What the movie shows:

- x_HII,V passes 0.1, 0.5 and 0.9 at z ≈ 8.2, 7.2 and 6.5, and reaches
  0.9999 by z = 5.5.
- τ_es accumulated in the run is 0.023. Adding full ionisation below
  z = 4.87 gives a total of 0.052.
- **Bubbles.** Ionised regions grow and reach about 2×10⁴ K, and they have
  merged by z ≈ 6.5. The last neutral patches are dense filaments. Light
  travel at ĉ = 0.01c limits the growth: the bubbles stay small for several
  frames after a source appears.
- **Bubble shapes.** The bubbles are nearly spherical because at 64³ only a
  few halos are resolved (2 sources at z ≈ 8, 7 by z ≈ 7). The 128³
  production run with abundance-matched sources will be much clumpier.
- **Gas ahead of the I-front is too hot.** Filaments sit at 10³–10⁴ K and
  voids near 10²–10³ K, hotter than an adiabatic IGM. This is the hydro
  suite's finding that the scheme spuriously heats cold supersonic gas
  (`docs/hydro_tests.md`). Don't read pre-reionisation temperatures from
  these movies.

Artefacts, all in the data rather than the rendering:

- **Grid-aligned streaks (quirk B16, `docs/port_notes.md`).** Straight,
  one-cell-wide streaks of partly neutral gas run along x and y through
  sources, giving small "+" crosses on the sources. In the 64³ run one of
  them runs down the centre of the z = 32 slice (x = 31, y = 39–45, at
  x_HII ≈ 0.1 inside a bubble).
  - **Cause, as far as the uniform test shows:** the original HEALPix ray
    frame is fixed to the grid, and sources sit at cell centres. In that
    frame ±x and ±y are pixel vertices at every nside, so the rays around
    each of these axes are mirror-symmetric about it. They all leave the
    column of cells along the axis at the same distance, and those cells
    get too few photons. The rays around ±z (the poles) are diagonal, so
    the deficit there is much weaker.
  - **It is in the original code, not from our changes** (it is not B10,
    and not `deterministic_parallel` or `c_reduced`). In the uniform
    single-source test (`tests/problems/uniform_1src`, 30 steps;
    `tools/analysis/axis_deficit.py`), the axis cells near the front fall
    short of their shell mean by up to 0.25 in x_HII along x and y, and by
    0.003 along z. This holds with the legacy flags (bit-identical to
    gfortran), with the fixed flags, and with `deterministic_parallel`. With
    `c_reduced = 0.1` the shortfall reaches 0.45.

    ![axis deficit](figures/axis_deficit.png)

  - **Still open:** in uniform gas the front lags by about one cell on the
    axes. The cosmological streaks are deeper, about 7 cells. They need a
    dedicated look before any science run, and a fix behind a flag, e.g. a
    rotation of the HEALPix frame per source or per step.
- **Radial spokes.** Faint spokes inside the bubbles follow the ray
  directions.

## Tests

`tools/movies/test_make_movie.py` (CTest `movie_tool`; skipped without
matplotlib) builds a synthetic 16³ run directory and checks:

- that pre-RT frames are skipped;
- the fallback to the field outputs;
- the front positions of a planar front, to 0.03 Mpc/h (0.06 cells);
- the τ_es integral, against a fine-grid evaluation (1e-3);
- t(z = 6) = 963 Myr;
- that both styles render an MP4 and the stills with zoom, axis, position
  and thickness options.
