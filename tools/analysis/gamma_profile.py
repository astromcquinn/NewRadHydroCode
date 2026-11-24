#!/usr/bin/env python3
"""Stacked radial profiles of the HI photoionization rate around the sources,
for runs that differ in the RT method (the proximity zones of
docs/rt_suite_512.md §5; rays, M1, OTVET and its closures).

usage: gamma_profile.py RUN:label [RUN:label ...] --outdir DIR
                        [--z 6.0,5.5,5.0] [--field gammaHI_avg]
                        [--sources auto|luv|rho[:T]|gamma[:T]] [--ngroups 3]
                        [--plot-group top] [--rmax 28] [--box 40] [--tol 0.1]
                        [--workers 16]

RUN is a run directory (<problem>/run) with step_NNNN/{stats.txt, fHI.bin,
<field>.bin} snapshots (gamma_power.py's layout). The first run is the
reference: for each target redshift its nearest snapshot is used (within
--tol), and every other run's snapshot nearest to that redshift.

Sources, taken once per redshift from the reference's snapshot and used for
every run (the source halos come from the dark matter and are the same in
every run up to the small gas differences):
  luv       Luv.bin (the source luminosity per RT cell; written when 'Luv'
            is in the run's `fields`); the cells with Luv > 0, grouped by Luv
  rho[:T]   the local maxima (3^3, periodic) of rhog.bin above T times its
            mean (default 100): a stand-in that does not depend on the RT;
            grouped by the peak density
  gamma[:T] the local maxima of the reference's Gamma above T times its mean
            over the ionised cells (default 10); grouped by the peak value.
            It selects on the reference's own field, so its central bins
            are biased up for the reference: use it as a check only
  auto      luv if Luv.bin exists, else rho
The sources are split into --ngroups groups of equal number by luminosity
(or its stand-in), from the faintest to the brightest; 'all' is every source.

The profile of a group is the mean over its sources of Gamma at separation
d (minimum image), averaged over spherical shells of |d| (edges 0, 0.5,
1.5, 2.5, 3.5, 5, 7, 10, 14, 20, 28, 40 cells, up to --rmax), and divided by
the run's mean Gamma over the ionised cells (x_HII > 0.5). It is computed
for all d at once as the cross-correlation of the source grid with Gamma
(single-precision FFTs): C(d) = sum_s Gamma(x_s + d) / N_s.

Writes DIR/gamma_profile.txt (every run, group, shell and redshift) and
DIR/gamma_profile.png: one column per redshift; top, Gamma(r)/<Gamma>_ion of
--plot-group with a 1/r^2 guide (streaming from the source: 1 + E (r0/r)^2,
anchored at the reference's 2-cell shell); bottom, the ratio to the
reference. The innermost point is the source cell itself (r < 0.5 cell),
drawn at r = 0.3 cell.
"""
import argparse
import os
import sys

import numpy as np
import scipy.fft as sfft
from scipy import ndimage

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gamma_power as gp  # noqa: E402

EDGES = np.array([0, 0.5, 1.5, 2.5, 3.5, 5, 7, 10, 14, 20, 28, 40], dtype=float)
# the colours and line styles of rt_suite_figures.py (the reference black)
COLORS = ["#000000", "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6"]
STYLES = ["-", "-", "-", "--", "-.", ":", "-", "--"]


def load_cube(path):
    a = np.fromfile(path, dtype=np.float32)
    N = round(a.size ** (1 / 3))
    if N ** 3 != a.size:
        raise ValueError(f"{path}: {a.size} floats is not a cube")
    return a.reshape(N, N, N)


def shell_index(N, rmax):
    """Shell index of every cell's minimum-image distance from the origin
    (-1 beyond rmax), and the shells' mean radii in cells."""
    d = np.minimum(np.arange(N), N - np.arange(N)).astype(np.float64)
    r = np.sqrt(d[:, None, None] ** 2 + d[None, :, None] ** 2 + d[None, None, :] ** 2)
    edges = EDGES[EDGES <= rmax + 1e-9]
    idx = np.digitize(r, edges) - 1
    idx[(r >= edges[-1])] = -1
    idx = idx.astype(np.int16)
    m = idx >= 0
    cnt = np.bincount(idx[m], minlength=len(edges) - 1).astype(np.float64)
    rmean = np.bincount(idx[m], weights=r[m], minlength=len(edges) - 1) / np.maximum(cnt, 1)
    return idx, cnt, rmean, edges


def local_maxima(f, thr):
    mx = ndimage.maximum_filter(f, size=3, mode="wrap")
    return np.argwhere((f == mx) & (f > thr))


def find_sources(mode, d_ref, gamma_ref, ion):
    """(cells (n, 3), luminosity or stand-in (n,), description)."""
    kind, _, t = mode.partition(":")
    if kind == "auto":
        kind = "luv" if os.path.exists(os.path.join(d_ref, "Luv.bin")) else "rho"
    if kind == "luv":
        L = load_cube(os.path.join(d_ref, "Luv.bin"))
        cells = np.argwhere(L > 0)
        return cells, L[tuple(cells.T)].astype(np.float64), "Luv > 0"
    if kind == "rho":
        rho = load_cube(os.path.join(d_ref, "rhog.bin"))
        T = float(t) if t else 100.0
        cells = local_maxima(rho, T * float(rho.mean()))
        return cells, rho[tuple(cells.T)].astype(np.float64), f"rhog peaks > {T:g} x mean"
    if kind == "gamma":
        T = float(t) if t else 10.0
        mean = float(gamma_ref[ion].mean()) if ion.any() else float(gamma_ref.mean())
        cells = local_maxima(gamma_ref, T * mean)
        return cells, gamma_ref[tuple(cells.T)].astype(np.float64), f"reference Gamma peaks > {T:g} x mean_ion"
    raise ValueError(f"unknown --sources {mode}")


def groups_of(cells, L, ngroups):
    order = np.argsort(L, kind="stable")
    out = [("all", cells)]
    if ngroups > 1 and len(cells) >= ngroups:
        parts = np.array_split(order, ngroups)
        names = ["bottom", "middle", "top"] if ngroups == 3 else [f"g{i + 1}" for i in range(ngroups)]
        out += [(n, cells[p]) for n, p in zip(names, parts)]
    return out


def stacked(cells, ghat, N, idx, cnt, workers):
    """Mean over the sources of the shell averages of the field with spectrum ghat."""
    S = np.zeros((N, N, N), dtype=np.float32)
    np.add.at(S, tuple(cells.T), 1.0)
    shat = sfft.rfftn(S, workers=workers)
    C = sfft.irfftn(np.conj(shat) * ghat, s=(N, N, N), workers=workers) / len(cells)
    m = idx >= 0
    return np.bincount(idx[m], weights=C[m].astype(np.float64), minlength=len(cnt)) / np.maximum(cnt, 1)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--z", default="6.0,5.5,5.0")
    ap.add_argument("--field", default="gammaHI_avg")
    ap.add_argument("--sources", default="auto")
    ap.add_argument("--ngroups", type=int, default=3)
    ap.add_argument("--plot-group", default="top")
    ap.add_argument("--rmax", type=float, default=28)
    ap.add_argument("--box", type=float, default=None)
    ap.add_argument("--tol", type=float, default=0.1)
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--dpi", type=int, default=200)
    a = ap.parse_args()

    gp.FIELD = a.field
    gp.S0_MODE = False
    os.makedirs(a.outdir, exist_ok=True)
    specs = [(r.split(":", 1)[0], r.split(":", 1)[1] if ":" in r else os.path.basename(r.rstrip("/")))
             for r in a.runs]
    box = a.box or gp.box_of(specs[0][0]) or 40.0
    zt = [float(x) for x in a.z.split(",")]
    snaps = [gp.snapshots(run, -1.0, 1e9) for run, _ in specs]
    if not snaps[0]:
        sys.exit("no snapshots of the reference")

    rows = []      # (z, run label, group, nsrc, shell, r_cells, value)
    panels = []    # (z, {label: profile}, rmean, nsrc of the plotted group, source description)
    idx_cache = {}
    for z0 in zt:
        zr = min(snaps[0], key=lambda z: abs(z - z0))
        if abs(zr - z0) > a.tol:
            print(f"z = {z0}: no reference snapshot within {a.tol}", file=sys.stderr)
            continue
        d_ref = snaps[0][zr]
        g_ref = gp.load_gamma(d_ref)
        N = g_ref.shape[0]
        if N not in idx_cache:
            idx_cache[N] = shell_index(N, a.rmax)
        idx, cnt, rmean, _ = idx_cache[N]
        ion_ref = load_cube(os.path.join(d_ref, "fHI.bin")) < 0.5
        cells, L, desc = find_sources(a.sources, d_ref, g_ref, ion_ref)
        if len(cells) == 0:
            print(f"z = {zr:.3f}: no sources ({desc})", file=sys.stderr)
            continue
        grp = groups_of(cells, L, a.ngroups)
        print(f"z = {zr:.3f}: {len(cells)} sources ({desc}), {d_ref}")
        prof = {}
        for i, ((run, label), sn) in enumerate(zip(specs, snaps)):
            zz = min(sn, key=lambda z: abs(z - zr)) if sn else None
            if zz is None or abs(zz - zr) > a.tol:
                print(f"  {label}: no snapshot near z = {zr:.3f}", file=sys.stderr)
                continue
            d = sn[zz]
            g = g_ref if i == 0 else gp.load_gamma(d)
            ion = ion_ref if i == 0 else load_cube(os.path.join(d, "fHI.bin")) < 0.5
            mean_ion = float(g[ion].mean()) if ion.any() else float(g.mean())
            ghat = sfft.rfftn(g.astype(np.float32), workers=a.workers)
            for gname, gc in grp:
                p = stacked(gc, ghat, N, idx, cnt, a.workers) / mean_ion
                for s in range(len(p)):
                    rows.append((zr, label, gname, len(gc), s, rmean[s], p[s]))
                if gname == a.plot_group:
                    prof[label] = p
            print(f"  {label}: z = {zz:.3f}, <Gamma>_ion = {mean_ion:.4g}")
        nplot = dict((n, len(c)) for n, c in grp).get(a.plot_group, 0)
        panels.append((zr, prof, rmean, nplot, desc))

    with open(os.path.join(a.outdir, "gamma_profile.txt"), "w") as f:
        f.write(f"# gamma_profile.py: field {a.field}, sources {a.sources}, box {box} Mpc/h\n")
        f.write("# Gamma(r)/<Gamma>_ion stacked over the sources of each group; r in cells (shell mean)\n")
        f.write("# z run group nsrc shell r_cells value\n")
        for zr, label, gname, ns, s, rc, v in rows:
            f.write(f"{zr:.4f} {label} {gname} {ns} {s} {rc:.3f} {v:.6g}\n")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if not panels:
        sys.exit("nothing to plot")
    labels = [l for _, l in specs]
    ref = labels[0]
    fig, axs = plt.subplots(2, len(panels), figsize=(3.3 * len(panels), 5.2), squeeze=False,
                            gridspec_kw={"height_ratios": [2, 1]}, sharex=True)
    for j, (zr, prof, rmean, nplot, desc) in enumerate(panels):
        dx = box / N
        x = rmean.copy()
        x[0] = 0.3
        x = x * dx
        top, bot = axs[0, j], axs[1, j]
        for i, label in enumerate(labels):
            if label not in prof:
                continue
            c, ls = COLORS[i % len(COLORS)], STYLES[i % len(STYLES)]
            lw = 2.0 if i == 0 else 1.4
            top.plot(x, prof[label], color=c, ls=ls, lw=lw, label=label)
            top.plot(x[:1], prof[label][:1], marker="o", ms=4, color=c, ls="none")
            if ref in prof:
                bot.plot(x, prof[label] / prof[ref], color=c, ls=ls, lw=lw)
        if ref in prof and len(prof[ref]) > 2:
            E = prof[ref][2] - 1.0
            r0 = rmean[2]
            rr = np.linspace(rmean[1], rmean[-1], 50)
            if E > 0:
                top.plot(rr * dx, 1 + E * (r0 / rr) ** 2, color="#9a9a96", lw=1.0, ls=(0, (4, 3)),
                         label=r"$1+E\,(r_0/r)^2$")
        top.set_xscale("log")
        top.set_yscale("log")
        top.set_title(f"z = {zr:.2f}, {a.plot_group} third ({nplot} sources)", fontsize=9)
        top.grid(alpha=0.3, lw=0.5)
        bot.set_xscale("log")
        bot.set_yscale("log")
        bot.set_ylim(0.2, 5)
        bot.axhline(1, color="#9a9a96", lw=0.8)
        bot.grid(alpha=0.3, lw=0.5)
        bot.set_xlabel(r"$r$ [cMpc/$h$] (first point: the source cell)")
        if j == 0:
            top.set_ylabel(r"$\Gamma(r)\,/\,\langle\Gamma\rangle_{\rm ion}$")
            bot.set_ylabel(f"ratio to {ref}")
            top.legend(fontsize=7, frameon=False)
    fig.suptitle(f"stacked {a.field} around the sources ({panels[0][4]})", fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(a.outdir, "gamma_profile.png"), dpi=a.dpi)
    fig.savefig(os.path.join(a.outdir, "gamma_profile.pdf"))
    print("wrote", os.path.join(a.outdir, "gamma_profile.{txt,png,pdf}"))


if __name__ == "__main__":
    main()
