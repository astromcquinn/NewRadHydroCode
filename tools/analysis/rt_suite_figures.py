#!/usr/bin/env python3
"""Paper figures of the uniform 512^3 RT-methods suite (docs/rt_suite_512.md):
side-by-side maps of x_HII and Gamma_HI of every method at a few ionised
fractions, and the 3D power spectra of both fields with the cross-correlation
coefficient against the reference (the case-B rays).

usage: rt_suite_figures.py RUN:label [RUN:label ...] --outdir DIR
                           [--xv 0.5,0.8,0.97] [--xv-gamma 0.5,0.8,0.97,1.0]
                           [--field gammaHI_avg]
                           [--box 40] [--nbins 24] [--slab 8] [--clip 0]

RUN is a run directory (<problem>/run) with step_NNNN/{fHI.bin,<field>.bin,
stats.txt} snapshots. The first run is the reference. --xv lists the target
ionised fractions of the x_HII figures, --xv-gamma those of the Gamma_HI
figures (x_V = 1 is empty for x_HII and the post-overlap point for Gamma).
For each target x_V
the reference snapshot with the nearest volume-weighted x_HII is chosen and
every other run's snapshot nearest in redshift is paired with it (the runs
write their snapshots at the same zdata1 redshifts); the panels are labelled
by the reference's x_V and z.

Writes to DIR:
  maps_xHII.png       rows: the chosen x_V; columns: the runs; the mean of a
                      --slab-cell slab through the box middle of x_HII
  maps_gamma.png      the same for log10 Gamma_HI (the --field cube, by
                      default the substep-averaged rate), fixed range
  pk_xHII.png         per x_V: Delta^2(k) of x_HII for every run, the ratio to
                      the reference, and r(k) = P_ab / sqrt(P_aa P_bb)
  pk_gamma.png        the same for delta_Gamma = Gamma/<Gamma>_V - 1 (with
                      --clip C, Gamma is capped at C <Gamma>_V first)
  pk.txt              the numbers behind the two spectra figures
  panels.txt          which snapshot of which run went into each panel
Uses the binning of gamma_power.py (log bins of |k| from the fundamental to
the Nyquist of the diagonal; Delta^2 = k^3 P / 2 pi^2).
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gamma_power as gp  # noqa: E402

# categorical colours in a fixed, colourblind-validated order (dataviz
# palette: blue, orange, aqua, yellow, magenta, green, violet); the
# reference is black
COLORS = ["#000000", "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6"]
STYLES = ["-", "-", "-", "--", "-.", ":", "-", "--"]


def read_stats(d):
    out = {}
    for line in open(os.path.join(d, "stats.txt")):
        if " = " in line:
            k, v = line.split(" = ", 1)
            try:
                out[k.strip()] = float(v)
            except ValueError:
                pass
    return out


def load_cube(path):
    a = np.fromfile(path, dtype=np.float32)
    N = round(a.size ** (1 / 3))
    return a.reshape(N, N, N)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--xv", default="0.5,0.8,0.97")
    ap.add_argument("--xv-gamma", default="0.5,0.8,0.97,1.0")
    ap.add_argument("--field", default="gammaHI_avg")
    ap.add_argument("--box", type=float, default=None)
    ap.add_argument("--nbins", type=int, default=24)
    ap.add_argument("--slab", type=int, default=8)
    ap.add_argument("--clip", type=float, default=0.0)
    ap.add_argument("--zmin", type=float, default=4.8)
    ap.add_argument("--zmax", type=float, default=7.2)
    ap.add_argument("--tol", type=float, default=0.1)
    ap.add_argument("--dpi", type=int, default=200)
    a = ap.parse_args()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    gp.FIELD = a.field
    os.makedirs(a.outdir, exist_ok=True)
    specs = [(r.split(":", 1)[0], r.split(":", 1)[1] if ":" in r else os.path.basename(r.rstrip("/")))
             for r in a.runs]
    box = a.box or gp.box_of(specs[0][0]) or 40.0
    targets_x = [float(x) for x in a.xv.split(",")]
    targets_g = [float(x) for x in a.xv_gamma.split(",")]
    targets = sorted(set(targets_x) | set(targets_g))

    # the reference's snapshots and their x_V; the panels
    snaps = [gp.snapshots(run, a.zmin, a.zmax) for run, _ in specs]
    ref = snaps[0]
    if not ref:
        sys.exit("no snapshots of the reference in the redshift range")
    ref_xv = {z: read_stats(d)["xHII_V"] for z, d in ref.items()}  # stats.txt of the full run
    panels = []
    for xv in targets:
        z = min(ref, key=lambda zz: abs(ref_xv[zz] - xv))
        dirs = []
        for i, (run, label) in enumerate(specs):
            if i == 0:
                dirs.append(ref[z])
                continue
            if not snaps[i]:
                dirs.append(None)
                continue
            zz = min(snaps[i], key=lambda q: abs(q - z))
            dirs.append(snaps[i][zz] if abs(zz - z) <= a.tol else None)
        panels.append((xv, z, ref_xv[z], dirs))
    with open(os.path.join(a.outdir, "panels.txt"), "w") as f:
        f.write("# target_xV ref_z ref_xV then one line per run: label step_dir z xHII_V\n")
        for xv, z, rx, dirs in panels:
            f.write(f"{xv} {z:.4f} {rx:.4f}\n")
            for (run, label), d in zip(specs, dirs):
                if d is None:
                    f.write(f"  {label} -\n")
                else:
                    st = read_stats(d)
                    f.write(f"  {label} {d} {st.get('redshift', float('nan')):.4f} {st.get('xHII_V', float('nan')):.4f}\n")

    nrun = len(specs)
    use = {"xHII": [i for i, p in enumerate(panels) if p[0] in targets_x],
           "gamma": [i for i, p in enumerate(panels) if p[0] in targets_g]}
    N = None
    binner = None
    # ---- maps and spectra in one pass over the cubes
    maps_x = [[None] * nrun for _ in panels]
    maps_g = [[None] * nrun for _ in panels]
    pk = {}  # (panel, field) -> list per run of (P_auto, P_cross_with_ref, mean)
    for ip, (xv, z, rx, dirs) in enumerate(panels):
        ref_fx = ref_fg = None
        for ir, d in enumerate(dirs):
            if d is None:
                continue
            fHI = load_cube(os.path.join(d, "fHI.bin"))
            if N is None:
                N = fHI.shape[0]
                binner = gp.Binner(N, box, a.nbins)
            xHII = 1.0 - fHI.astype(np.float64)
            del fHI
            gam = gp.load_gamma(d).astype(np.float64)
            lo, hi = N // 2 - a.slab // 2, N // 2 - a.slab // 2 + a.slab
            maps_x[ip][ir] = xHII[lo:hi].mean(axis=0)
            maps_g[ip][ir] = gam[lo:hi].mean(axis=0)
            # x_HII: Delta^2 of the field itself (the 21-cm convention)
            fx = np.fft.rfftn(xHII)
            # Gamma: the contrast of the volume mean, optionally capped
            gm = gam.mean()
            if a.clip > 0:
                gam = np.minimum(gam, a.clip * gm)
                gm = gam.mean()
            fg = np.fft.rfftn(gam / gm - 1.0)
            del xHII, gam
            if ir == 0:
                ref_fx, ref_fg = fx, fg
            for key, fa, fr in (("xHII", fx, ref_fx), ("gamma", fg, ref_fg)):
                pa = binner.power(fa)
                pc = binner.power(fa, fr) if fr is not None else pa
                pk.setdefault((ip, key), [None] * nrun)[ir] = (pa, pc)
            del fx, fg
        print(f"panel {ip}: target x_V {xv}, reference z = {z:.3f} x_V = {rx:.3f}", flush=True)

    k = binner.k
    pref = {"k": k}
    # ---- figures: maps
    plt.rcParams.update({"font.size": 9, "axes.titlesize": 9, "axes.labelsize": 9,
                         "xtick.labelsize": 8, "ytick.labelsize": 8, "figure.dpi": a.dpi})
    for name, key, maps, kw, cblabel, floor in (
            ("maps_xHII.png", "xHII", maps_x, dict(vmin=0, vmax=1, cmap="cividis"), r"$x_{\rm HII}$", None),
            ("maps_gamma.png", "gamma", maps_g, dict(norm=LogNorm(vmin=1e-14, vmax=1e-11), cmap="magma"),
             r"$\Gamma_{\rm HI}$ [s$^{-1}$]", 1e-15)):
        rows = use[key]
        npan = len(rows)
        w = 2.0
        fig, axes = plt.subplots(npan, nrun, figsize=(w * nrun + 0.6, w * npan + 0.3), squeeze=False)
        im = None
        for irow, ip in enumerate(rows):
            xv, z, rx, dirs = panels[ip]
            for ir, (run, label) in enumerate(specs):
                ax = axes[irow, ir]
                ax.set_xticks([])
                ax.set_yticks([])
                m = maps[ip][ir]
                if m is None:
                    ax.text(0.5, 0.5, "no snapshot", ha="center", va="center", transform=ax.transAxes)
                    continue
                if floor is not None:
                    m = np.maximum(m, floor)  # neutral cells at the bottom of the scale, not blank
                im = ax.imshow(m.T, origin="lower", extent=(0, box, 0, box), interpolation="nearest", **kw)
                if irow == 0:
                    ax.set_title(label)
                if ir == 0:
                    ax.set_ylabel(f"$x_V$ = {rx:.2f}, z = {z:.2f}")
        fig.subplots_adjust(left=0.04, right=0.92, bottom=0.02, top=0.95, wspace=0.03, hspace=0.03)
        cax = fig.add_axes([0.93, 0.1, 0.012, 0.8])
        fig.colorbar(im, cax=cax, label=cblabel)
        fig.savefig(os.path.join(a.outdir, name))
        plt.close(fig)

    # ---- figures: spectra
    with open(os.path.join(a.outdir, "pk.txt"), "w") as f:
        f.write(f"# box {box} Mpc/h, N {N}, field {a.field}, clip {a.clip}; Delta^2 = k^3 P / 2 pi^2\n")
        for key, ylab, fname in (("xHII", r"$\Delta^2_{x_{\rm HII}}(k)$", "pk_xHII.png"),
                                 ("gamma", r"$\Delta^2_{\delta_\Gamma}(k)$", "pk_gamma.png")):
            cols_use = use[key]
            npan = len(cols_use)
            fig, axes = plt.subplots(3, npan, figsize=(3.0 * npan + 0.5, 6.4), squeeze=False, sharex=True)
            for icol, ip in enumerate(cols_use):
                xv, z, rx, dirs = panels[ip]
                rows = pk[(ip, key)]
                pr = rows[0][0]
                f.write(f"# {key}: panel x_V = {rx:.3f} z = {z:.3f}; columns k [h/Mpc] then per run Delta2, ratio, r\n")
                f.write("# runs: " + " ".join(l for _, l in specs) + "\n")
                cols = [k]
                for ir, (run, label) in enumerate(specs):
                    if rows[ir] is None:
                        continue
                    pa, pc = rows[ir]
                    d2 = k ** 3 * pa / (2 * np.pi ** 2)
                    ratio = pa / np.maximum(pr, 1e-300)
                    r = pc / np.sqrt(np.maximum(pa * pr, 1e-300))
                    cols += [d2, ratio, r]
                    c, ls = COLORS[ir % len(COLORS)], STYLES[ir % len(STYLES)]
                    lw = 1.6 if ir == 0 else 1.2
                    axes[0, icol].loglog(k, d2, color=c, ls=ls, lw=lw, label=label)
                    if ir > 0:
                        axes[1, icol].semilogx(k, ratio, color=c, ls=ls, lw=lw, label=label)
                        axes[2, icol].semilogx(k, r, color=c, ls=ls, lw=lw, label=label)
                for row in np.array(cols).T:
                    f.write(" ".join(f"{v:.6g}" for v in row) + "\n")
                axes[0, icol].set_title(f"$x_V$ = {rx:.2f}, z = {z:.2f}")
                axes[1, icol].axhline(1, color="0.6", lw=0.8)
                axes[2, icol].axhline(1, color="0.6", lw=0.8)
                axes[1, icol].set_yscale("log")
                axes[1, icol].set_ylim(0.03, 30)
                axes[2, icol].set_ylim(0, 1.05)
                axes[2, icol].set_xlabel(r"$k$ [$h$ Mpc$^{-1}$]")
                for r_ in range(3):
                    axes[r_, icol].grid(True, which="major", color="0.9", lw=0.6)
                    for s in ("top", "right"):
                        axes[r_, icol].spines[s].set_visible(False)
            axes[0, 0].set_ylabel(ylab)
            axes[1, 0].set_ylabel(r"$P / P_{\rm ref}$")
            axes[2, 0].set_ylabel(r"$r(k)$ with ref")
            axes[0, 0].legend(frameon=False, fontsize=7, ncol=2)
            fig.tight_layout()
            fig.savefig(os.path.join(a.outdir, fname))
            plt.close(fig)
    print("wrote", a.outdir)


if __name__ == "__main__":
    main()
