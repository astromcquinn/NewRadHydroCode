#!/usr/bin/env python3
"""Write MUSIC configurations for the zoom ladder (docs/zoom_turbulence.md
section 3 and 4.1; tools/ics/music/README.md) and the uniform C26 ladder
(section 3b).

Zoom: the base level (2^base_level cells over the box) refined to levelmax.
The region is given once for the whole ladder, at the deepest level the
ladder will reach (--ladder-max, default levelmax): a centre and an extent
(box units or Mpc/h), or the Lagrangian bounding box of a list of particle
ids of a base run plus a pad. From it the grid of every level follows by
MUSIC's own rule (mesh.hh, refinement_hierarchy: a level is the next finer
one halved and grown by `padding` cells on each side, rounded outward to
even cells), and each rung is written with its finest grid exact
(ref_center + ref_dims, dyadic numbers, so MUSIC's floating point
reproduces the integers). Every rung of a ladder therefore has the same grid
at every level it shares with another, which is what the ladder property
needs: MUSIC2 constrains a level's white noise to its parent by Fourier
splicing over the level's whole patch (random_music_wnoise_generator.cc),
so the noise of a level depends on its grid, and a rung that put level l
elsewhere (as ref_extent at each rung's own levelmax does) would redraw it.
The seeds are a function of (seed, level) only,

    seed[l] = seed + 1000 (l - base_level),   l = base_level .. levelmax,

the same in every rung (tests/music/ladder_check.py checks the noise and
the fields). The levels between base and levelmax are thin shells of
`padding` cells each (MUSIC's default 8; it must exceed the code's
amrex_n_error_buf and proper nesting, section 4.2), so a shallower rung's
finest region is the deepest one's grown by those shells.

Uniform ladder (--uniform): levelmin = levelmin_TF = levelmax = base_level and
the seed given at --seed-level (>= every rung's level), from which MUSIC
restricts the noise down to the rung: every rung of the C1 ladder shares its
modes when --seed-level is the finest rung's.

Region from particle ids (--ids FILE): 1-based ids of a run started from
xv_dm.init with --ndm particles per dimension (the id is the IC order, the
lattice cell n = id - 1 with i fastest: tools/ics/music_convert.py,
tools/ics/make_ics.py); the bounding box of their lattice positions is taken
on the periodic circle per axis (the shortest arc holding all of them), then
grown by --pad Mpc/h on every side.

Usage:
  music_ladder.py Z1.conf --box 40 --base-level 9 --levelmax 13 --ladder-max 17 \\
      --center 0.5 0.5 0.5 --extent-mpch 0.3125 [--seed 2026] [--pad 0]
  music_ladder.py Z2.conf ... --levelmax 15 --ladder-max 17 (the same region)
  music_ladder.py Z1.conf --box 40 --base-level 9 --levelmax 13 --ladder-max 17 \\
      --ids region_ids.txt --ndm 512 --pad 0.2
  music_ladder.py OUT.conf --box 0.25 --base-level 8 --uniform --seed-level 10
Cosmology and transfer: --cosmology planck_c26 (default) with
--transfer camb_file --transfer-file camb_z300.dat (default), or
--transfer eisenstein.
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pmcosmo as pc  # noqa: E402


def level_seeds(seed, base_level, levelmax):
    return {lev: seed + 1000 * (lev - base_level) for lev in range(base_level, levelmax + 1)}


def periodic_bbox(x):
    """Smallest periodic interval [lo, lo + w] (box units, lo in [0, 1))
    holding every x in [0, 1): the complement of the largest gap."""
    s = np.sort(np.mod(x, 1.0))
    if len(s) == 1:
        return float(s[0]), 0.0
    gaps = np.diff(np.concatenate([s, [s[0] + 1.0]]))
    i = int(np.argmax(gaps))
    lo = s[(i + 1) % len(s)]
    return float(lo), float(1.0 - gaps[i])


def region_from_ids(ids, ndm):
    """Lagrangian cell centres (box units) of 1-based particle ids."""
    n = np.asarray(ids, dtype=np.int64) - 1
    if n.min() < 0 or n.max() >= ndm ** 3:
        sys.exit("ids out of range for --ndm")
    i, j, k = n % ndm, (n // ndm) % ndm, n // (ndm * ndm)
    return [(c + 0.5) / ndm for c in (i, j, k)]


def region(args):
    """Lower corner and extent of the region (box units, unshifted; the
    corner may be negative when the region wraps)."""
    if args.ids:
        ids = np.loadtxt(args.ids, dtype=np.int64, ndmin=1)
        q = region_from_ids(ids, args.ndm)
        lo, w = zip(*(periodic_bbox(c) for c in q))
        # cell centres: the cells themselves reach half a cell further
        lo = [li - 0.5 / args.ndm for li in lo]
        w = [wi + 1.0 / args.ndm for wi in w]
    else:
        if args.center is None or (args.extent is None and args.extent_mpch is None):
            sys.exit("give --center and --extent/--extent-mpch, or --ids and --ndm")
        ext = args.extent if args.extent is not None else [e / args.box for e in args.extent_mpch]
        w = ext * 3 if len(ext) == 1 else list(ext)
        lo = [c - 0.5 * e for c, e in zip(args.center, w)]
    pad = args.pad / args.box
    lo = [li - pad for li in lo]
    w = [wi + 2 * pad for wi in w]
    if max(w) >= 1:
        sys.exit(f"the region ({w}) covers the box")
    return lo, w


def even_floor(x):
    return x - x % 2


def even_ceil(x):
    return x + x % 2


def grid_chain(lo, w, base_level, top_level, padding):
    """The grid of every level base_level + 1 .. top_level, (lo, hi) in that
    level's cells (unshifted; Python ints), from the region at top_level and
    MUSIC's coarsening rule (gridding unit 2). MUSIC works in its shifted
    frame, where an even shift per level leaves the rounding unchanged."""
    n = 1 << top_level
    g = {top_level: ([even_floor(int(np.floor(x * n))) for x in lo],
                     [even_ceil(int(np.ceil((x + e) * n))) for x, e in zip(lo, w)])}
    for lev in range(top_level - 1, base_level, -1):
        l1, h1 = g[lev + 1]
        g[lev] = ([even_floor(a // 2 - padding) for a in l1], [even_ceil(-(-b // 2) + padding) for b in h1])
    for lev, (l1, h1) in g.items():
        if any(h - lo_ >= (1 << lev) for lo_, h in zip(l1, h1)):
            sys.exit(f"level {lev}'s grid ({l1}, {h1}) covers the box: fewer levels or a smaller region")
    return g


def write_conf(args, path):
    cosmo = pc.PRESETS[args.cosmology]
    lines = [f"# written by tools/ics/music_ladder.py: {' '.join(sys.argv[1:])}",
             "[setup]", f"boxlength   = {args.box!r}", f"zstart      = {args.zstart!r}"]
    if args.uniform:
        lv = args.base_level
        lines += [f"levelmin    = {lv}", f"levelmin_TF = {lv}", f"levelmax    = {lv}"]
    else:
        lo, w = region(args)
        top = args.ladder_max or args.levelmax
        if top < args.levelmax:
            sys.exit("--ladder-max must be >= --levelmax")
        g = grid_chain(lo, w, args.base_level, top, args.padding)
        l1, h1 = g[args.levelmax]
        n = 1 << args.levelmax
        centre = [float(np.mod((a + b) / 2 / n, 1.0)) for a, b in zip(l1, h1)]
        dims = [b - a for a, b in zip(l1, h1)]
        lines += [f"levelmin    = {args.base_level}", f"levelmin_TF = {args.base_level}",
                  f"levelmax    = {args.levelmax}",
                  # no spaces: MUSIC strips them with remove_if but not erase,
                  # so "26, 28, 28" would parse as 26, 28, 2828
                  f"ref_center  = {centre[0]!r},{centre[1]!r},{centre[2]!r}",
                  f"ref_dims    = {dims[0]},{dims[1]},{dims[2]}"]
        lines.insert(1, f"# the ladder's grids (region at level {top}; unshifted, own cells, [lo, hi)):")
        for lev in range(args.base_level + 1, top + 1):
            mark = "" if lev <= args.levelmax else "   (deeper rungs)"
            lines.insert(1 + lev - args.base_level, f"#   level {lev}: lo {g[lev][0]} hi {g[lev][1]} "
                         f"size {[b - a for a, b in zip(*g[lev])]}{mark}")
    lines += [f"padding     = {args.padding}", f"overlap     = {args.overlap}", "align_top   = no",
              f"baryons     = {'yes' if args.baryons else 'no'}", "use_2LPT    = no"]
    if args.no_shift:
        lines.append("no_shift    = yes")
    lines += ["", "[cosmology]", f"Omega_m       = {cosmo.omegam!r}", f"Omega_L       = {cosmo.omegal!r}",
              f"Omega_b       = {cosmo.omegab!r}", f"H0            = {100 * cosmo.hubble0!r}",
              f"sigma_8       = {cosmo.sigma8!r}", f"n_s           = {cosmo.ns!r}",
              f"transfer      = {args.transfer}"]
    if args.transfer == "camb_file":
        lines.append(f"transfer_file = {args.transfer_file}")
    lines += ["", "[random]"]
    if args.uniform:
        lines.append(f"seed[{args.seed_level or args.base_level}] = {args.seed}")
    else:
        for lev, s in level_seeds(args.seed, args.base_level, args.levelmax).items():
            lines.append(f"seed[{lev}] = {s}")
    lines += ["", "[output]", "format   = generic",
              f"filename = {args.output or os.path.splitext(os.path.basename(path))[0] + '.hdf5'}",
              "", "[poisson]", "fft_fine      = yes", "accuracy      = 1e-5",
              "grad_order    = 6", "laplace_order = 6", ""]
    with open(path, "w") as f:
        f.write("\n".join(lines))
    return lines


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out")
    ap.add_argument("--box", type=float, required=True, help="Mpc/h")
    ap.add_argument("--base-level", type=int, required=True)
    ap.add_argument("--levelmax", type=int)
    ap.add_argument("--ladder-max", type=int, default=None,
                    help="the deepest level of the ladder, where the region is given (default levelmax)")
    ap.add_argument("--uniform", action="store_true")
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--seed-level", type=int, default=None, help="--uniform: the level the seed is given at")
    ap.add_argument("--center", type=float, nargs=3)
    ap.add_argument("--extent", type=float, nargs="+", help="box units (1 or 3 values)")
    ap.add_argument("--extent-mpch", type=float, nargs="+", help="Mpc/h (1 or 3 values)")
    ap.add_argument("--ids", help="file of 1-based particle ids of the base run")
    ap.add_argument("--ndm", type=int, help="particles per dimension of that run")
    ap.add_argument("--pad", type=float, default=0.0, help="Mpc/h added on every side of the region")
    ap.add_argument("--padding", type=int, default=8)
    ap.add_argument("--overlap", type=int, default=4)
    ap.add_argument("--no-shift", action="store_true", help="keep MUSIC from re-centring the box")
    ap.add_argument("--zstart", type=float, default=300.0)
    ap.add_argument("--cosmology", choices=sorted(pc.PRESETS), default="planck_c26")
    ap.add_argument("--transfer", default="camb_file")
    ap.add_argument("--transfer-file", default="camb_z300.dat")
    ap.add_argument("--no-baryons", dest="baryons", action="store_false")
    ap.add_argument("--output", help="MUSIC's output file name (default OUT.hdf5)")
    args = ap.parse_args(argv)
    if not args.uniform and (args.levelmax is None or args.levelmax <= args.base_level):
        ap.error("--levelmax > --base-level is needed for a zoom (or --uniform)")
    if args.ids and not args.ndm:
        ap.error("--ids needs --ndm")
    write_conf(args, args.out)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
