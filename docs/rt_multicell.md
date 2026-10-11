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

<!-- K1 -->

<!-- RESULTS -->
