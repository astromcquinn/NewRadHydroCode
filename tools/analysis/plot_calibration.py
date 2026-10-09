#!/usr/bin/env python3
"""Reionization histories and tau_es of an f_esc scan (calibrate_fesc.py).

usage: plot_calibration.py CAL_DIR [--out fig.png]
"""
import argparse
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import reion  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("caldir")
    ap.add_argument("--out")
    a = ap.parse_args()
    runs = []
    for d in glob.glob(os.path.join(a.caldir, "fesc_*")):
        if os.path.exists(os.path.join(d, "run", "movie", "frames.txt")) and \
                os.path.exists(os.path.join(d, "finished")):
            runs.append((float(os.path.basename(d)[5:]), os.path.join(d, "run")))
    runs.sort()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    cmap = plt.get_cmap("Blues")
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    taus = []
    for i, (f, r) in enumerate(runs):
        c = cmap(0.4 + 0.6 * i / max(1, len(runs) - 1))
        _, z, xV, xM = reion.read_frames(r)
        zz, tau = reion.tau_history(z, xM)
        taus.append(tau[-1])
        ax.plot(z, xV, "-", color=c, lw=2, label=f"f_esc = {f:g}  (τ = {tau[-1]:.4f})")
        ax.plot(z, xM, "--", color=c, lw=1)
    ax.set_xlabel("z")
    ax.set_ylabel(r"$\langle x_{\rm HII}\rangle$ (solid: volume, dashed: mass)")
    ax.set_xlim(12, 5)
    ax.legend(frameon=False, fontsize=8)
    fs = np.array([f for f, _ in runs])
    bx.axhspan(reion.PLANCK_TAU - reion.PLANCK_SIG, reion.PLANCK_TAU + reion.PLANCK_SIG,
               color="#ffb347", alpha=0.3, lw=0)
    bx.axhline(reion.PLANCK_TAU, color="#d08a1e", lw=1)
    bx.text(fs.min(), reion.PLANCK_TAU + 0.0005, "Planck 2018 (±1σ band)", color="#d08a1e")
    bx.plot(fs, taus, "o-", color="#2a6fbd", lw=2)
    bx.set_xscale("log")
    bx.set_xlabel(r"f$_{\rm esc}$")
    bx.set_ylabel(r"τ$_{\rm es}$")
    for x in (ax, bx):
        x.grid(alpha=0.2, lw=0.5)
        for sp in ("top", "right"):
            x.spines[sp].set_visible(False)
    out = a.out or os.path.join(a.caldir, "calibration.png")
    fig.savefig(out, dpi=130)
    print(out)


if __name__ == "__main__":
    main()
