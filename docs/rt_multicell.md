# Rays crossing k RT cells per substep (`rt_cells_per_substep`)

*Written 2026-10-11 on branch `rt-multicell`. The test of a coarse RT
substep for the RT on AMR (Stage F, `docs/framework_decision.md` §6):
instead of sub-cycling a level-l region with 2^l RT substeps per coarse
substep, keep the coarse substep everywhere and let rays travel 2^l fine
cells per substep there. On the uniform grid that is "rays travel k cells
per substep, with a k times longer substep at the same ĉ". If k = 2 and 4
agree with k = 1 to the acceptance metrics of `docs/ROADMAP.md`, the AMR RT
needs no fine-level time sub-cycling.*

<!-- SUMMARY -->

## 1. The scheme

### 1.1 What the original substep does

`calc_timestep` sets `nstep_rt = nint(dl_c/dl_rt · ĉ/c)`, the light
crossings of one RT cell in the hydro step, and `dtsec_rt = dt/nstep_rt`.
In each substep every ray advances one cell length: `ray_segment` splits
the path at the cell faces into up to `nseg` segments (3, or 4 with
`legacy_rt = 0`); `cast_1` (step 1) deposits in each segment's cell
max(dl/sum(dl remaining), 1 − e^(−τ dl)) of the photons not yet deposited,
with τ from the previous substep; `photo_cell` (step 2) absorbs part of
each cell's deposit and sets the transmitted fraction `fgam`; `cast_3`
(step 3) attenuates the ray by (1 − fgam) of each segment, moves it to the
end of its path and relinks it. A ray already crosses up to 3 faces per
substep with the lagged optical depth; k only makes the path longer.

### 1.2 What changes with k (all in the shared kernels, so the serial code and the AMReX backend agree)

| Piece | k = 1 (original) | k > 1 |
|---|---|---|
| `calc_timestep` | nstep_rt = nint(dl_c/dl_rt · ĉ/c), carry and ≥ 1 guard on that count | nint(dl_c/(k dl_rt) · ĉ/c): k times longer substeps, the light distance per hydro step unchanged; carry and guard on the k-cell count |
| `ray_segment` | path of 1 cell, `nseg` = 4 (`legacy_rt = 0`) | path of k cells, `nseg` = 3k + 1 (at most k faces per axis), `rrk::kNsegMax` = 49 (k ≤ 16) |
| `cast_1`, `cast_3` | loop over the segments | the same loop; the `fdl` normalisation, `seg_min`, the `1e-6f` overshoot and the removal threshold unchanged |
| photo-ionisation rate in `photo_gas_cell` | Γ = n_dep σ len_rt/dtsec_rt | Γ = n_dep σ **k** len_rt/dtsec_rt (§1.3) |
| split radius clamp in `split` | [rmin, rmin + 1] | [rmin, rmin + **k**] (§1.3) |
| source emission | nside_init = 4 (192 rays) | the same for k ≤ 4; for k > 4 the nside whose split radius covers the first path (§1.3; 768 rays for k = 5–8) |
| `redshift_per_distance` | s = nstep_rt dl_rt/dl_c | s = nstep_rt **k** dl_rt/dl_c |
| `deterministic_parallel` | 8³ domains; chain insertions within d/2 − 1 cells | d ≥ 2k + 2 (8 for k ≤ 3, 16 at 64³ for k = 4); insertions in step 1 within d/2 − k, in step 3 within d/2 |
| AMReX backend | 2 ghost layers, local `Redistribute` over 3 cells | both at least k + 1 (the reach assertion is unchanged and passes, §2) |
| checkpoint | version 5 | version 6 records k; a restart with another k converts the substep carry |

Unchanged, because they are per substep: emission (ngam ∝ dtsec_rt, so a
source's photons per hydro step are the same), the thick-source local
deposit, the split and merge criteria, the ray-pool guard (it now counts
the rays a source emits), the merge cones, `recomb_rad` (its photons are
buffered per substep and emitted at the next one; the photon budget is
unchanged), `ray_rotate` (keyed on the substep number: fewer substeps,
fewer rotations, the same statistics), the movie frames (written at the
substep that crosses the frame redshift: k times coarser in time).

Restrictions: `legacy_rt = 0` (the legacy mode is the Fortran's, and its 3
segments cannot cover even one cell, B8), `rt_method = rays` (the moment
solvers keep k = 1 and refuse k > 1), 1 ≤ k ≤ 16. **Light must cross at
least k cells per hydro step**: otherwise the nstep_rt ≥ 1 guard (with the
substep carry it is applied silently) makes the radiation faster than ĉ.
The `ifront_rsla` problem at 0.1c has 2.3 cells of light per hydro step, so
k = 4 runs at 1.7 ĉ there (§3).

### 1.3 Three generalisations that k forced (found in this study)

**The photo-ionisation rate needs the crossing speed.** A ray carries the
photons emitted in its substep. With k cells per substep they occupy a
path of k cells, and a cell crossed in full receives the fraction 1/k of
them as its deposit. Physically every photon of the ray crosses that cell
during the substep, so the absorption probability per photon is τ_cell and
the rate is Γ = (photons crossing per unit area and time) σ = N_ray σ/(A
dtsec_rt) = k n_dep σ len_rt/dtsec_rt with n_dep = deposit/volume. The
original's len_rt/dtsec_rt is ĉ for one cell per substep; the general form
is k len_rt/dtsec_rt = ĉ. Without the factor, Γ in ionised gas is k times
too low and x_HI there k times too high: 4.4× at k = 4 in the uniform test
(§3, the `NORATE` diagnostic build), while the front radius hardly changes.
The `multicell` CTest checks the interior x_HI for this reason.

**The split radius clamp.** A ray past its split radius has moved at most
one substep, 1 cell in the original, beyond it, and `split` clamps the
daughters' radius to [rmin, rmin + 1] (it never binds for k = 1). With k
cells per substep it would bind and move the effective origin forward by
up to k − 1 cells; [rmin, rmin + k] keeps it a no-op. In the uniform tests
the two clamps give identical results up to k = 4, because rays emitted at
r = 0 and moving 4 cells per substep overshoot the split radii 3.9, 7.8,
15.6, 31.3 cells by less than one cell until nside = 64; at k = 8 the
overshoot is 4.1 cells at the first split (§5).

**The emission nside for k > 4.** A ray splits at most once per substep. A
source ray emitted at nside_init = 4 (split radius 3.9 cells) that travels
8 cells in its first substep stays at nside 4 over that path and the next,
about 4 cells² per ray at r = 8 against `ray_per_cell` = 1: the uniform
I-front at k = 8 had partly neutral shells inside the sphere (§3).
`rrk::emit_nside(k)` emits at the smallest nside ≥ nside_init with
rmax(nside) ≥ k − 1 (nside_init for k ≤ 4); after the first substep one
split per substep keeps up, since rmax doubles while r grows by k.

## 2. Two more problems, and two options

### 2.1 B17: `photo_cell` ionises photon-starved cells without the photons (`photon_conserve`)

`photo_gas_cell` iterates the implicit ionisation balance (ne to 1e-3, at
most 5 times). In each iteration the densities are computed with the
current rates, then the rates are rescaled so that the photo-ionisations
do not exceed the cell's photons. In a photon-starved cell with Γ dt ≫ 1 (a
neutral front cell that gets fewer photons than it has atoms), backward
Euler ionises the cell almost fully at every iteration (x_HI → x_HI/(1 +
Γdt)), so ne "converges" at once while the rates still exceed the photons;
the densities are kept and only the photons of the rescaled rates are
taken. The ionised hydrogen then exceeds the counted net ionisations. It
is a quirk of the original (now **B17** in `docs/port_notes.md`), small at
k = 1, but Γ_max dt in a thick front cell scales as k² (k times the photons
per ray and the factor k of the rate). On `lognormal_16src` (c = 1, 20
steps):

| k | x_HII,M | net ionisations | excess ("leak") | with `photon_conserve = 1` |
|---|---|---|---|---|
| 1 | 0.1745 | 0.1714 | 0.0031 (1.8 %) | 4e-5 |
| 2 | 0.1819 | 0.1700 | 0.0119 (6.5 %) | |
| 4 | 0.1962 | 0.1640 | 0.0323 (16.5 %) | 5e-5 |

**`photon_conserve = 1`** (default 0 = the original; independent of
`legacy_rt` and `legacy_physics`, so the Fortran-mirrored modes are
untouched; rays, moments and the AMReX backend) also iterates while any
frequency group is limited by more than 1e-3, at most 1000 times. The
budget then closes to 4e-5 per H atom; `radiative_physics` costs the same
(1.99 s against 1.97 s on that run). CTest `photon_conserve`.

### 2.2 The in-flight photons jump to the end of the path (`ray_deposit = exposure`)

The share deposit gives segment n the fraction max(dl_n/Σ_{m≥n} dl_m,
fabs_n) of the photons **not yet deposited**, so the deposits sum to the
ray's photons. Let a ray cross an ionised stretch and then reach a neutral
cell within one substep. The neutral cell receives only the share of the
remaining path; the photons deposited in the ionised stretch are mostly
transmitted there (fgam ≈ 1) and leave with the ray, whose new position is
the end of the path, beyond the front. In the next substep they are
deposited up to k − 1 cells beyond the I-front. With one cell per substep
the jump stays inside a cell; with k cells it smears every front over
about k cells, whatever the front speed. This is what makes the uniform
I-front at k = 8 partly neutral inside its volume-equivalent radius even
at c = 1 (§3).

Physically every photon of the ray crosses every cell of its path during
the substep, so a front reached in mid-path should get all photons that
arrive there. **`ray_deposit = exposure`** (default `share` = the original;
`legacy_rt = 0`; any k) deposits in segment n

  dep_n = A_n · max(dl_n/L, fabs_n),  A_{n+1} = A_n (1 − fabs_n),

with A_1 = the ray's photons, L the path length and fabs_n = 1 − e^(−τ_n
dl_n) from the previous substep's τ as before: dl_n/L is the time the ray's
photons spend in the cell (with the factor k of the rate this gives the
exact thin-cell absorption A_n τ_n dl_n), fabs_n the absorbed fraction in a
thick cell, and only the absorption estimate, not the deposit, is taken
from the photons passed on. `cast_3` removes what each cell actually
absorbed, dep_n (1 − fgam_n), from the photons that arrive and moves the
ray only while some are left. In thin cells it is the share deposit. The
deposits no longer sum exactly to the ray's photons, so photon
conservation holds only to the accuracy of the lagged absorption estimate:
in the uniform test the total ionisations stay within 0.2 % of k = 1 at
k = 4 and 8 (the share deposit: +1.5 % and +2.3 %).

## 3. Validation of k = 1

Every new option defaults to the original (`rt_cells_per_substep = 1`,
`photon_conserve = 0`, `ray_deposit = share`), and the generalisations are
written so that k = 1 is the same arithmetic (1 × dl_rt, `rmin + 1`, the
4-segment step, nside_init, `nray_init × nsrc`):

- **Fast suite:** `main` (07aa378) and this branch both pass 112/112 in the
  cloud container (gcc 13.3, gfortran 13.3; `healpix_vs_healpy` skipped in
  both), including every stored-reference and live-Fortran comparison in
  both legacy modes, the restart and the thread tests. <!-- FINALSUITE -->
- **Smoke diffs:** `rt_driver` outputs byte-identical to `main` on
  `uniform_8src` (`deterministic_parallel`, 2 threads) and, after each later
  commit, to the previous binary on `lognormal_16src` (fixed flags,
  `ray_rotate = 3`) and `uniform_8src` (legacy flags).
- **AMReX backend** (`-DRADRAY_AMREX=ON`, AMReX 26.10, OpenMPI 4.1.6):
  `ctest -L amrex` 50/50 before the exposure option (the bitwise tests
  against `rt_serial` at k = 1 included); the multicell tests below after
  it.
- **New tests** (all pass): `rays` (unit: k-cell segments cover the path in
  3k+1 segments and k = 1 is the 4-segment step exactly; the clamp),
  `multicell` (k = 2 vs 1 on the uniform I-front at c = 1: front within half
  a cell, finite-c solution, interior x_HI, front profile),
  `multicell_exposure` (the same at k = 4 with `ray_deposit = exposure` and
  `photon_conserve = 1`), `photon_conserve` (B17 resolved and fixed),
  `thread_rt_det_multicell_lognormal_16src` (k = 3, the largest k with 8³
  domains, byte-identical for 1, 2, 4 threads), `restart_smoke_det_multicell`
  (k = 2 with the substep carry, bit-identical restart), and in the AMReX
  build `amrex_roundoff_multicell{2,4,4x}_{uniform_8src,lognormal_16src}`
  (the backend with k = 2 and 4, and 4 with the exposure deposit, against
  the serial port with the same options to round-off; ghosts and local
  `Redistribute` grow to k + 1, so the backend supports k > 1 rather than
  refusing it).


## 4. Step 2: the uniform I-front (`uniform_1src`, finite-c solution)

`tools/analysis/rt_multicell_ifront.py` builds the problem of
`tests/ifront_rsla_check.py` (uniform neutral gas at the mean density,
z = 10, one source of Luv = 8e28 in the centre cell) and runs `rt_driver`
(`legacy_rt = legacy_physics = 0`, `rt_substep_carry = 1`) for every k and
c. It measures the front radius R = (3 V_ion/4π)^(1/3) after every hydro
step, against k = 1 and against the finite-c solution R = R_inf(t − R/ĉ)
of Shapiro et al. (2006) (R_inf from the c = 1, k = 1 run; ĉ = the light
front of the run, k cells per substep, over t), and the shell-mean x_HI and
T around the source. "share" is the scheme as specified (share deposit,
original photon limit); "pc" adds `photon_conserve = 1`; "exposure" uses
`ray_deposit = exposure`. Columns: the largest |R − R(k=1)| over the steps
with R ≥ 2 cells and at the last step; the relative error against the
finite-c solution (max / rms, R ≥ 2); the largest difference of the
shell-mean x_HI from k = 1 at the profile steps (the front profile); the
mean x_HI of the shells inside R/2 at the last step over k = 1 (the
photo-ionisation rate); the casts over the run; the wall time (one
thread; the runs shared the machine with other jobs).

**The original problem** (N = 32, da = 2e-4, 15 steps). At 0.1c it has
2.3 cells of light per hydro step, so k = 4 cannot keep ĉ (light 1.67× too
far: the ≥ 1 guard) and k = 2 has 1–2 substeps per step:

| scheme | c | k | max \|ΔR\| (R ≥ 2) | ΔR last | finite-c max / rms | shell x_HI max \|Δ\| | x_HI inside R/2 | casts | wall [s] |
|---|---|---|---|---|---|---|---|---|---|
| share | 0.1 | 1 | 0.00 | +0.00 | 0.107 / 0.040 | 0.000 | 1.00 | 1.35e+05 | 1 |
| share | 0.1 | 2 | 0.91 | +0.16 | 0.113 / 0.055 | 0.618 | 1.05 | 5.86e+04 | 0 |
| share | 0.1 | 4 | 0.58 | +0.57 | 0.052 / 0.045 | 0.999 | 1.00 | 6.03e+04 | 1 | light 1.67×
| share | 1 | 1 | 0.00 | +0.00 | 0.000 / 0.000 | 0.000 | 1.00 | 1.81e+06 | 5 |
| share | 1 | 2 | 0.06 | +0.02 | 0.014 / 0.007 | 0.072 | 1.04 | 7.38e+05 | 3 |
| share | 1 | 4 | 0.15 | +0.07 | 0.030 / 0.018 | 0.378 | 1.15 | 3.93e+05 | 2 |
| exposure+pc | 0.1 | 1 | 0.00 | +0.00 | 0.106 / 0.034 | 0.000 | 1.00 | 1.59e+05 | 1 |
| exposure+pc | 0.1 | 2 | 0.90 | -0.02 | 0.105 / 0.041 | 0.580 | 1.02 | 5.47e+04 | 1 |
| exposure+pc | 0.1 | 4 | 0.35 | +0.29 | 0.113 / 0.036 | 0.999 | 1.00 | 5.16e+04 | 1 | light 1.67×
| exposure+pc | 1 | 1 | 0.00 | +0.00 | 0.000 / 0.000 | 0.000 | 1.00 | 1.97e+06 | 5 |
| exposure+pc | 1 | 2 | 0.03 | -0.00 | 0.010 / 0.003 | 0.004 | 0.98 | 6.90e+05 | 3 |
| exposure+pc | 1 | 4 | 0.10 | -0.00 | 0.029 / 0.008 | 0.053 | 0.96 | 3.13e+05 | 2 |

**A 64³ version** (same box and source, cells half the size, da = 4e-4,
15 steps; 9.3 cells of light per hydro step at 0.1c, 93 at c = 1), with
k = 8 as the first step-5 case. The "share" k = 8 row predates the
emission rule of §1.3 (sources at nside 4: partly neutral shells inside
the sphere, x_HI inside R/2 at 36–41× that of k = 1); the other k = 8 rows
use it:

| scheme | c | k | max \|ΔR\| (R ≥ 2) | ΔR last | finite-c max / rms | shell x_HI max \|Δ\| | x_HI inside R/2 | casts | wall [s] |
|---|---|---|---|---|---|---|---|---|---|
| share | 0.1 | 1 | 0.00 | +0.00 | 0.034 / 0.011 | 0.000 | 1.00 | 7.52e+06 | 21 |
| share | 0.1 | 2 | 0.30 | +0.14 | 0.044 / 0.021 | 0.114 | 1.00 | 2.33e+06 | 9 |
| share | 0.1 | 4 | 2.28 | +0.59 | 0.254 / 0.075 | 0.496 | 1.00 | 9.28e+05 | 6 |
| share | 0.1 | 8 | 5.21 | +0.48 | 0.341 / 0.096 | 1.000 | 36.29 | 3.88e+05 | 4 |
| share | 1 | 1 | 0.00 | +0.00 | 0.000 / 0.000 | 0.000 | 1.00 | 9.06e+07 | 216 |
| share | 1 | 2 | 0.06 | -0.00 | 0.006 / 0.002 | 0.112 | 1.00 | 2.73e+07 | 83 |
| share | 1 | 4 | 0.22 | +0.03 | 0.026 / 0.009 | 0.388 | 1.00 | 1.03e+07 | 58 |
| share | 1 | 8 | 0.59 | +0.21 | 0.068 / 0.029 | 0.941 | 40.87 | 4.35e+06 | 30 |
| share+pc | 0.1 | 1 | 0.00 | +0.00 | 0.024 / 0.007 | 0.000 | 1.00 | 7.40e+06 | 27 |
| share+pc | 0.1 | 2 | 0.23 | -0.06 | 0.031 / 0.009 | 0.127 | 1.00 | 2.23e+06 | 14 |
| share+pc | 0.1 | 4 | 2.24 | +0.02 | 0.257 / 0.068 | 0.424 | 1.00 | 8.63e+05 | 8 |
| share+pc | 0.1 | 8 | 5.14 | -0.28 | 0.294 / 0.082 | 1.000 | 1.00 | 4.34e+05 | 6 |
| share+pc | 1 | 1 | 0.00 | +0.00 | 0.000 / 0.000 | 0.000 | 1.00 | 9.05e+07 | 236 |
| share+pc | 1 | 2 | 0.01 | -0.01 | 0.001 / 0.000 | 0.118 | 1.00 | 2.72e+07 | 114 |
| share+pc | 1 | 4 | 0.03 | -0.02 | 0.001 / 0.001 | 0.365 | 1.00 | 1.02e+07 | 68 |
| share+pc | 1 | 8 | 0.03 | -0.03 | 0.004 / 0.002 | 0.829 | 1.00 | 4.80e+06 | 47 |
| exposure+pc | 0.1 | 1 | 0.00 | +0.00 | 0.031 / 0.009 | 0.000 | 1.00 | 9.22e+06 | 28 |
| exposure+pc | 0.1 | 2 | 0.18 | -0.04 | 0.034 / 0.010 | 0.090 | 1.00 | 2.66e+06 | 13 |
| exposure+pc | 0.1 | 4 | 2.28 | +0.01 | 0.256 / 0.068 | 0.326 | 1.00 | 8.15e+05 | 8 |
| exposure+pc | 0.1 | 8 | 5.19 | -0.32 | 0.296 / 0.082 | 1.000 | 1.00 | 4.00e+05 | 7 |
| exposure+pc | 1 | 1 | 0.00 | +0.00 | 0.000 / 0.000 | 0.000 | 1.00 | 1.02e+08 | 260 |
| exposure+pc | 1 | 2 | 0.01 | -0.01 | 0.000 / 0.000 | 0.006 | 1.00 | 2.87e+07 | 118 |
| exposure+pc | 1 | 4 | 0.04 | -0.00 | 0.003 / 0.001 | 0.037 | 1.00 | 9.66e+06 | 68 |
| exposure+pc | 1 | 8 | 0.15 | +0.02 | 0.018 / 0.005 | 0.067 | 1.00 | 3.79e+06 | 40 |

![front radius against k = 1](figures/rt_multicell_ifront.png)

*R − R(k=1) for the three variants at 0.1c (top) and c = 1 (bottom); the
grey band is ±0.5 cell.*

![profiles c = 1](figures/rt_multicell_profiles_c1.png)

![profiles c = 0.1](figures/rt_multicell_profiles_c01.png)

*Shell-mean x_HI and T at steps 2, 6 and 15 (t = 6.5, 19.5, 49 Myr) for
k = 1, 2, 4, 8: share deposit (top two rows) and exposure deposit (bottom
two), both with `photon_conserve`; c = 1 and 0.1c. The spikes of x_HI at
r ≈ 7.5 and 14 cells are the B16 ring (`ray_rotate = 0` here), the same for
every k.*

What the tables say:

- **c = 1 (fronts photon-limited, v ≪ ĉ):** the front radius agrees to
  0.01–0.04 cells for k ≤ 4 and 0.03–0.15 cells for k = 8 with
  `photon_conserve`; with the original limiter it drifts by 0.22 (k = 4)
  and 0.59 cells (k = 8), which is B17. The front *profile* needs the
  exposure deposit: with the share deposit the shell x_HI is off by 0.37
  (k = 4) and 0.83 (k = 8) across the front, because photons jump past it
  (§2.2); with the exposure deposit by 0.004–0.04 (k ≤ 4) and 0.02–0.07
  (k = 8), the profiles lying on top of each other. The interior x_HI is
  that of k = 1 to 1 % for every variant with the rate factor of §1.3.
- **0.1c (the front is light-limited early):** k = 2 stays within 0.2–0.3
  cells (front profile 0.09–0.13 at step 2, then ≤ 0.04). k = 4 is 2.3
  cells behind at step 1 and k = 8 5.2 cells, whatever the deposit; both
  recover later (k = 4 to 0.01–0.02 cells at the end, k = 8 to 0.3 cells
  with a profile still off by 0.6). This is the lagged optical depth: in
  step 1 the front moves at 0.7ĉ, i.e. 2.8 cells per k = 4 substep, and a
  front cell that the previous substep's τ calls thick receives all of the
  ray's photons; if it ionises within the substep the transmitted photons
  leave with the ray to the end of its path, past the cells that should
  have been ionised next. The criterion is **v_front · dtsec_rt < 1 cell,
  i.e. v_front < ĉ/k**: from step 3 the front moves at 0.18ĉ (0.7 cells
  per k = 4 substep, 1.5 per k = 8 substep) and k = 4 is fine, k = 8 not
  yet.
- **Cost:** k = 2, 4, 8 take 2.2×, 3.8× and 6.5× less wall time at c = 1
  (exposure + pc), with 3.5×, 10.6× and 27× fewer casts: the saving is in
  the per-substep passes over the grid (calc_rt_1, the photo-chemistry,
  the emission) and in fewer, longer casts, not in fewer rays (§7).

## 5. Step 3: `lognormal_16src` and `uniform_8src`

`tools/analysis/rt_multicell_fields.py`: `rt_driver` on the two 32³
problems with the production ray settings of
`tools/ray_compression/run_precision.sh` (`legacy_rt = legacy_physics = 0`,
`nray_max_factor = 16`), `rt_substep_carry = 1`, c = 1 (the problems' own
setting: 11.6 cells of light per hydro step, enough for k = 8),
`lognormal_16src` for 100 steps (x_HII,V from 0 to 0.95) and `uniform_8src`
for 80 (to 0.23). Each k > 1 is compared with k = 1 of the same variant.
The per-cell columns use the relative difference |a − b|/max(|a|, |b|,
1e-6 max|a|) at the last step; `compare.py --steps last` gives the maxima
(1.00 everywhere: single cells always differ completely somewhere).

| problem, scheme | k | x_HII,V / x_HII,M at the end (rel. to k = 1) | max over steps: x_V, x_M, T_M | front shift max / last [cells] | rms rel. diff. last step: x_HI (ion.), T (ion.), ngam, Γ (ion.) | max rel. diff. (compare.py): fHI, temp, ngam |
|---|---|---|---|---|---|---|
| lognormal, share | 1 | 0.9495 / 0.8167 | | | | |
| lognormal, share | 2 | +0.48% / +0.92% | 109.74%, 12.15%, 17.1% | 0.25 / 0.12 | 32.0%, 4.9%, 18%, 31% | 1.00, 1.00, 1.00 |
| lognormal, share | 4 | +1.29% / +2.01% | 277.78%, 28.21%, 32.6% | 0.54 / 0.18 | 32.9%, 8.0%, 20%, 34% | 1.00, 1.00, 1.00 |
| lognormal, share | 8 | +2.97% / +6.06% | 968.10%, 38.71%, 51.0% | 1.13 / 0.31 | 34.8%, 15.2%, 34%, 33% | 1.00, 1.00, 1.00 |
| lognormal, share+pc | 1 | 0.9460 / 0.8100 | | | | |
| lognormal, share+pc | 2 | -0.41% / -1.18% | 46.25%, 2.94%, 7.3% | 0.17 / 0.11 | 34.4%, 5.1%, 18%, 32% | 1.00, 1.00, 1.00 |
| lognormal, share+pc | 4 | -0.34% / -2.12% | 71.74%, 7.99%, 11.2% | 0.33 / 0.16 | 33.4%, 7.6%, 19%, 34% | 1.00, 1.00, 1.00 |
| lognormal, share+pc | 8 | -4.44% / -13.37% | 91.75%, 21.48%, 26.0% | 0.77 / 0.51 | 65.7%, 15.1%, 49%, 55% | 1.00, 1.00, 1.00 |
| lognormal, exposure+pc | 1 | 0.9487 / 0.8177 | | | | |
| lognormal, exposure+pc | 2 | -0.83% / -1.79% | 5.22%, 2.77%, 2.9% | 0.14 / 0.13 | 35.2%, 5.1%, 19%, 32% | 1.00, 1.00, 1.00 |
| lognormal, exposure+pc | 4 | -1.26% / -2.87% | 34.71%, 8.39%, 7.3% | 0.16 / 0.16 | 33.1%, 6.2%, 21%, 32% | 1.00, 1.00, 1.00 |
| lognormal, exposure+pc | 8 | -7.04% / -14.18% | 103.38%, 22.19%, 21.1% | 0.66 / 0.66 | 67.9%, 18.1%, 54%, 57% | 1.00, 1.00, 1.00 |
| uniform_8src, share | 1 | 0.2310 / 0.2310 | | | | |
| uniform_8src, share | 2 | -0.09% / -0.09% | 2.50%, 2.50%, 25.3% | 0.24 / 0.17 | 45.1%, 4.9%, 36%, 28% | 1.00, 1.00, 1.00 |
| uniform_8src, share | 4 | +0.26% / +0.26% | 16.50%, 16.50%, 41.1% | 0.84 / 0.08 | 40.4%, 6.3%, 35%, 45% | 1.00, 1.00, 1.00 |
| uniform_8src, share | 8 | +2.70% / +2.70% | 17.21%, 17.21%, 44.9% | 1.13 / 1.13 | 67.2%, 48.2%, 55%, 58% | 1.00, 1.00, 1.00 |
| uniform_8src, share+pc | 1 | 0.2308 / 0.2308 | | | | |
| uniform_8src, share+pc | 2 | -0.19% / -0.19% | 1.51%, 1.51%, 24.9% | 0.24 / 0.17 | 45.2%, 4.9%, 36%, 29% | 1.00, 1.00, 1.00 |
| uniform_8src, share+pc | 4 | -0.30% / -0.30% | 2.78%, 2.78%, 36.9% | 0.84 / 0.08 | 41.0%, 5.8%, 36%, 45% | 1.00, 1.00, 1.00 |
| uniform_8src, share+pc | 8 | -0.08% / -0.08% | 33.09%, 33.09%, 45.3% | 0.94 / 0.94 | 69.5%, 41.9%, 56%, 57% | 1.00, 1.00, 1.00 |
| uniform_8src, exposure+pc | 1 | 0.2314 / 0.2314 | | | | |
| uniform_8src, exposure+pc | 2 | -0.04% / -0.04% | 3.16%, 3.16%, 14.1% | 0.02 / 0.00 | 4.1%, 0.9%, 27%, 6% | 0.40, 0.94, 1.00 |
| uniform_8src, exposure+pc | 4 | +0.06% / +0.06% | 3.91%, 3.91%, 22.4% | 0.02 / 0.00 | 4.2%, 1.3%, 29%, 9% | 0.35, 1.00, 1.00 |
| uniform_8src, exposure+pc | 8 | +0.64% / +0.64% | 27.81%, 27.81%, 37.0% | 0.16 / 0.06 | 33.3%, 4.2%, 40%, 19% | 1.00, 1.00, 1.00 |

**The per-cell differences are mostly ray noise.** Two equivalent k = 1
runs (share and exposure deposit, both with `photon_conserve`) differ by
20 % (lognormal) and 32 % (uniform_8src) rms in x_HI of ionised cells, 25 %
and 20 % in Γ, 14 % and 19 % in the photon counts and 2.7 % and 2.2 % in T:
the instantaneous Γ of a cell comes from the few rays that crossed it in
the last substep. The x_HI, ngam and Γ columns at k = 2–4 are at this
floor. T of ionised cells (5–8 % rms) is above it; §6 shows why.

What k does:

- **The front:** with the share deposit the mean I-front displacement
  reaches 0.25 (k = 2), 0.54 (k = 4) and 1.1 cells (k = 8) on lognormal and
  0.24, 0.84, 1.1 on uniform_8src; with the exposure deposit 0.14, 0.16,
  0.66 and 0.02, 0.02, 0.16. Uniform gas with the exposure deposit is the
  only case that meets half a cell at k = 8.
- **The global history:** on uniform_8src x_HII agrees to 0.3 % at the end
  for every variant up to k = 4 (0.06 % with exposure + pc; the maxima in
  the table are the first two steps, when the fronts are a few cells and
  light-limited). On lognormal the end state depends on the limiter: with
  the original one x_HII,M rises with k (+0.9, +2.0, +6.1 %), because the
  B17 leak grows with k (the excess over the counted net ionisations is
  0.005, 0.016, 0.057, 0.145 per H atom for k = 1, 2, 4, 8); with
  `photon_conserve` it falls with k (−1.2/−1.8 % at k = 2, −2.1/−2.9 % at
  k = 4, −13/−14 % at k = 8, share/exposure) with more recombinations
  (+4 %, +6 %, +11 % at k = 2, 4, 8). The deposit rule hardly matters for
  these global numbers; §6 shows that the k = 8 deficit is a substep-length
  effect.
- **The rays:** the live rays per substep are about the same for every k
  (k = 1: 2.0e5 mean, 3.4e5 at the end; k = 8: 1.2–2.3e5 and 3.8–4.3e5):
  a ray lives as long in time, emissions are k times rarer but k times
  brighter. Splits per substep rise with k (1.5e4 → 2.3e4, 3.0e4, 3.3e4 on
  lognormal), since a ray passes its split radius by up to k cells per
  substep; merges first rise and fall again at k = 8. The casts over the run
  fall by 1.8×, 3.3×, 6× (lognormal) and 2.3×, 3–6×, 9–15× (uniform_8src).
  At k = 8 the live rays oscillate from substep to substep (the burst of an
  emission, then the merges).

![ray statistics](figures/rt_multicell_rays_lognormal.png)

*Live rays, splits and merges per substep on `lognormal_16src` (exposure +
pc), on a common time axis (k = 1 substeps).*

![metrics vs k](figures/rt_multicell_fields.png)

*The mean front displacement (max over the run), the final x_HII,M and the
rms T of ionised cells against k for the six problem/variant pairs; dotted:
half a cell, 1 %, 3 %.*

![uniform_8src slices, share](figures/rt_multicell_slices_uniform8_share.png)
![uniform_8src slices, exposure](figures/rt_multicell_slices_uniform8_exposure.png)

*`uniform_8src` at step 80, slice through the box centre: share deposit
with the original limiter (top), exposure deposit with `photon_conserve`
(bottom). The share deposit's fronts at k = 8 are ragged, partly neutral
over several cells; the exposure deposit's are those of k = 1.*

![lognormal slices](figures/rt_multicell_slices_lognormal.png)

*`lognormal_16src` at step 100 (exposure + pc): the self-shielded
filament through the sources is larger at k = 8.*

## 6. Path length or substep length? (diagnostic runs)

A run with k cells per substep at speed ĉ differs from k = 1 in two ways:
each ray crosses k cells per substep (the transport), and every substep
is k times longer (everything else per substep: the optical depth lagged
by one substep, the photo-chemistry and thermal update over dtsec_rt, the
emission bursts). Two pairs of runs separate them: k = 8 at ĉ = 8c has the
substep length of k = 1 at c, and k = 1 at ĉ = c/8 has that of k = 8 at c
(the speed of light then differs, which matters little in these
photon-limited problems).

`lognormal_16src` (share + pc unless noted, 100 steps):

| run | substep length | x_HII,V | x_HII,M | T_M [K] |
|---|---|---|---|---|
| k = 1, c = 1 | 1 | 0.946 | 0.810 | 14254 |
| k = 8, c = 8 | 1 | 0.961 | 0.817 | 14085 |
| k = 8, c = 1 | 8 | 0.904 | 0.702 | 12631 |
| k = 1, c = 1/8 | 8 | 0.861 | 0.696 | 12601 |
| k = 8, c = 1, exposure + pc | 8 | 0.882 | 0.702 | 12306 |
| k = 8, c = 1, exposure + pc, photo-chemistry sub-cycled 8× | 8 | 0.911 | 0.745 | 12754 |

`uniform_8src` (exposure + pc), mean T of the ionised cells:

| run | substep length | step 3 | step 10 |
|---|---|---|---|
| k = 1, c = 1 | 1 | 22793 | 20524 |
| k = 4, c = 4 | 1 | 23159 | 20365 |
| k = 4, c = 1 | 4 | 25790 | 21122 |
| k = 1, c = 1/4 | 4 | 25735 | 20931 |
| k = 4, c = 1, photo-chemistry sub-cycled 4× | 4 | 25766 | 20795 |

**The transport across k cells per substep is accurate; the errors at
fixed ĉ come from the k times longer substep.** k = 8 at 8c reproduces
k = 1 to 1 % in x_HII,M and T_M, while k = 1 at c/8 shows the same 14 %
deficit in x_HII,M as k = 8 at c; on `uniform_8src` the early T excess of
k = 4 appears at k = 1 with 4× longer substeps and disappears at k = 4
with substeps of the k = 1 length. Sub-cycling the photo-chemistry inside
the substep (the same deposit and rate, the photons spread evenly over
the sub-steps: a diagnostic build, not committed) recovers about a third
of the x_HII,M deficit and none of the T excess. The rest is the optical
depth lagged by one, now longer, substep: a cell that the previous
substep's τ calls thick receives all arriving photons as its deposit, and
its rate, with the crossing-speed factor k, assumes the time-share deposit
of a thin cell; if the cell has become thin within the substep it absorbs
up to k times the physical amount, over-ionising and over-heating it,
and photons skip the cells behind a front that moved more than one cell
(§4). Both effects are strongest where τ changes fastest: near sources,
at fronts and in dense, quickly recombining gas.

What would remove the rest (not done here): a separate exposure field for
the rate (Γ from the time-share A dl/L of every cell, the photon limit
from the deposit), and a second, corrector pass of step 1 with the updated
τ in cells whose τ changed by more than a factor 2 within the substep.

<!-- RESULTS -->
