#!/usr/bin/env python3
"""Reionization history and Thomson optical depth of apps/radhydro runs.

    reion_history.py RUN1[:label] RUN2[:label] ... [--out fig.png] [--table t.md]

Reads RUN/movie/frames.txt (frame, istep, istep_rt, z, x_HII,V, x_HII,M;
written every movie_dz from inside the RT substeps). Frames before the RT
starts (x_HII,M = nan) are skipped. Prints the redshifts at which x_HII,V and
x_HII,M reach 10, 25, 50, 75, 90, 99 and 99.9 per cent, the volume- and
mass-weighted fractions at fixed redshifts, and tau_es.

tau_es = sigma_T c n_H0 int x_e (1+z)^2 / H(z) dz, with x_e = x_HII,M (1 +
Y/4X) (HeII follows HII; HeIII and the residual ionisation before the RT are
ignored) above the last frame, x_e = 1 + Y/4X from there down to z = 3 and
1 + Y/2X below z = 3 (HeII reionization assumed instantaneous at z = 3). The
cosmology is the code's (include/radray/constants.hpp): Omega_m = 0.27,
Omega_b = 0.045, h = 0.70, Y = 0.24, flat.
"""
import argparse
import os

import numpy as np

trapezoid = np.trapezoid if hasattr(np, "trapezoid") else np.trapz  # numpy < 2.0

OM, OB, H, Y = 0.27, 0.045, 0.70, 0.24
SIGMA_T = 6.6524587e-25          # cm^2
C = 2.99792458e10                # cm/s
MP = 1.67262192e-24              # g
H0 = 100 * H * 1e5 / 3.0856776e24  # 1/s
RHOC = 3 * H0**2 / (8 * np.pi * 6.6743e-8)
NH0 = (1 - Y) * OB * RHOC / MP   # comoving H density [cm^-3]
FHE = Y / (4 * (1 - Y))          # He/H by number
LEVELS = [0.10, 0.25, 0.50, 0.75, 0.90, 0.99, 0.999]
ZFIX = [10, 9, 8, 7.5, 7, 6.5, 6, 5.5]


def read(run):
    d = np.genfromtxt(os.path.join(run, "movie", "frames.txt"), comments="#")
    d = d[np.isfinite(d[:, 5])]
    return d[:, 3], d[:, 4], d[:, 5]  # z, x_V, x_M (z decreasing)


def z_at(z, x, level):
    i = np.nonzero(x >= level)[0]
    if len(i) == 0:
        return float("nan")
    i = i[0]
    if i == 0:
        return z[0]
    # linear interpolation between the bracketing frames
    return z[i - 1] + (level - x[i - 1]) * (z[i] - z[i - 1]) / (x[i] - x[i - 1])


def tau_es(z, xm):
    hz = lambda zz: H0 * np.sqrt(OM * (1 + zz) ** 3 + 1 - OM)
    zz = np.concatenate([np.linspace(0, z[-1], 2000, endpoint=False), z[::-1]])
    xe = np.empty_like(zz)
    low = zz < z[-1]
    xe[low] = 1 + FHE + np.where(zz[low] < 3, FHE, 0)
    xe[~low] = xm[::-1] * (1 + FHE)
    f = xe * (1 + zz) ** 2 / hz(zz)
    return SIGMA_T * C * NH0 * trapezoid(f, zz)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--out", help="figure (png)")
    ap.add_argument("--table", help="write the Markdown tables to this file")
    a = ap.parse_args()
    runs = []
    for r in a.runs:
        path, _, label = r.partition(":")
        runs.append((path, label or os.path.basename(os.path.normpath(path))))
    data = {lab: read(p) for p, lab in runs}

    out = []
    out.append("| run | weighting | " + " | ".join(f"z({l * 100:g}%)" for l in LEVELS)
               + " | tau_es |")
    out.append("|---" * (len(LEVELS) + 3) + "|")
    for lab, (z, xv, xm) in data.items():
        t = tau_es(z, xm)
        for w, x in (("V", xv), ("M", xm)):
            out.append(f"| {lab} | {w} | " + " | ".join(f"{z_at(z, x, l):.2f}" for l in LEVELS)
                       + (f" | {t:.4f} |" if w == "V" else " | |"))
    out.append("")
    out.append("| run | weighting | " + " | ".join(f"z = {zf:g}" for zf in ZFIX) + " |")
    out.append("|---" * (len(ZFIX) + 2) + "|")
    for lab, (z, xv, xm) in data.items():
        for w, x in (("V", xv), ("M", xm)):
            vals = [np.interp(zf, z[::-1], x[::-1], left=np.nan, right=np.nan) for zf in ZFIX]
            out.append(f"| {lab} | {w} | " + " | ".join("—" if not np.isfinite(v) else f"{v:.3f}" for v in vals) + " |")
    text = "\n".join(out)
    print(text)
    if a.table:
        open(a.table, "w").write(text + "\n")

    if a.out:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        colors = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
        ink, muted, grid = "#2b2b2b", "#6b6b6b", "#e4e4e0"
        fig, axs = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
        for (lab, (z, xv, xm)), c in zip(data.items(), colors):
            axs[0].plot(z, xv, color=c, lw=2, label=f"{lab}, volume")
            axs[0].plot(z, xm, color=c, lw=2, ls="--", label=f"{lab}, mass")
            axs[1].semilogy(z, np.clip(1 - xv, 1e-6, 1), color=c, lw=2, label=lab)
        axs[0].set_ylabel(r"ionized fraction $x_{\rm HII}$")
        axs[1].set_ylabel(r"neutral fraction $1 - x_{\rm HII,V}$")
        for ax in axs:
            ax.set_xlabel("redshift z")
            ax.invert_xaxis()
            ax.grid(color=grid, lw=0.8)
            ax.tick_params(colors=muted)
            for s in ("top", "right"):
                ax.spines[s].set_visible(False)
            for s in ("left", "bottom"):
                ax.spines[s].set_color(muted)
            ax.yaxis.label.set_color(ink)
            ax.xaxis.label.set_color(ink)
        axs[0].legend(frameon=False, fontsize=8, labelcolor=ink)
        axs[1].legend(frameon=False, fontsize=8, labelcolor=ink)
        fig.savefig(a.out, dpi=150)
        print("wrote", a.out)


if __name__ == "__main__":
    main()
