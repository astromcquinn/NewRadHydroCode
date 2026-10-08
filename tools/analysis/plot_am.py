#!/usr/bin/env python3
"""Figure for an abundance-matching fit (apps/am_fit output directory).

Left: L_UV(M_halo) from abundance matching (points; filled = inside the
fitted magnitude range) and the parametric fit (lines), per LF redshift.
Right: input UV LFs (dashed) and the LFs the fit implies (solid).

usage: plot_am.py AM_OUT_DIR [out.png]
"""
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    d = sys.argv[1]
    out = sys.argv[2] if len(sys.argv) > 2 else d + "/am_fit.png"
    pts = np.loadtxt(d + "/am_points.txt")
    rows = [l.split() for l in open(d + "/lf_check.txt") if not l.startswith("#")]
    zs = sorted(set(pts[:, 0]))
    # One hue, light -> dark with redshift (ordered quantity)
    cmap = plt.get_cmap("Blues")
    col = {z: cmap(0.45 + 0.55 * i / max(1, len(zs) - 1)) for i, z in enumerate(zs)}

    fig, (a, b) = plt.subplots(1, 2, figsize=(11, 4.6), constrained_layout=True)
    for z in zs:
        s = pts[pts[:, 0] == z]
        inr = s[:, 7] == 1
        a.plot(s[:, 1], s[:, 6], "-", lw=2, color=col[z], label=f"z = {z:g}")
        a.plot(s[inr, 1], s[inr, 5], "o", ms=4, color=col[z])
        a.plot(s[~inr, 1], s[~inr, 5], "o", ms=4, mfc="none", color=col[z])
    a.axvline(9.0, color="0.5", lw=1, ls=":")
    a.text(9.03, a.get_ylim()[0] + 0.3, "128³ / 10 Mpc/h\nresolution limit", fontsize=8, color="0.35")
    a.set_xlabel(r"log$_{10}$ M$_{\rm halo}$ [M$_\odot$]")
    a.set_ylabel(r"log$_{10}$ L$_{\rm UV}$ [erg s$^{-1}$ Hz$^{-1}$]")
    a.set_title("Abundance matching (points) and fit (lines)", fontsize=10)
    a.legend(frameon=False, fontsize=8)

    for z in zs:
        r = [x for x in rows if float(x[0]) == z]
        M = np.array([float(x[2]) for x in r])
        pin = np.array([float(x[3]) for x in r])
        pm = np.array([float(x[4]) for x in r])
        rng = np.array([x[5] == "1" for x in r])
        b.plot(M, np.log10(pin), "--", lw=1, color=col[z])
        b.plot(M[rng], np.log10(pin[rng]), "-", lw=3, alpha=0.35, color=col[z])
        b.plot(M, np.log10(np.maximum(pm, 1e-30)), "-", lw=1.5, color=col[z], label=f"z = {z:g}")
    b.set_xlim(-23.5, -12)
    b.set_ylim(-8, 1)
    b.set_xlabel(r"M$_{\rm UV}$ (observed)")
    b.set_ylabel(r"log$_{10}$ $\phi$ [Mpc$^{-3}$ mag$^{-1}$]")
    b.set_title("Input LF (dashed; thick = fitted range) vs implied by fit", fontsize=10)
    b.legend(frameon=False, fontsize=8, loc="upper left")
    for ax in (a, b):
        ax.grid(alpha=0.2, lw=0.5)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    fig.savefig(out, dpi=130)
    print(out)


if __name__ == "__main__":
    main()
