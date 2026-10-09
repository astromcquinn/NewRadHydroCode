#!/usr/bin/env python3
"""Power spectra and linear growth of apps/radhydro runs.

    pk_growth.py RUN1[:label] RUN2[:label] ... [--out fig.png] [--kmax-growth 3]

Reads RUN/Pk/pk.z=ZZ.ZZZZ (calc_power_spectrum: per k bin, k [h/Mpc], the DM
Delta^2 raw and its scatter, the DM Delta^2 with the CIC window deconvolved
and its scatter, the gas Delta^2 and its scatter; Delta^2 = k^3 P / 2 pi^2)
and the box from RUN/../params.txt. Prints

  * growth: the DM Delta^2 of the lowest k bins at each output divided by
    the one at the first output, against the linear [D(z)/D(z0)]^2 (the same
    modes, so the cosmic variance of the small box cancels);
  * the ratio of the DM Delta^2 (CIC-deconvolved) to the linear theory of
    the ICs (tools/ics/pmcosmo.py: Eisenstein & Hu no-wiggle, sigma8, n_s,
    full growth) at the first and last outputs;
  * the gas/DM ratio of Delta^2 at k = 2, 5, 10, 20 h/Mpc.
"""
import argparse
import glob
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "ics"))
import pmcosmo as pc  # noqa: E402


def read_run(run):
    out = []
    for f in glob.glob(os.path.join(run, "Pk", "pk.z=*")):
        z = float(re.search(r"z=([0-9.]+)", f).group(1))
        d = np.loadtxt(f)
        out.append((z, d))
    out.sort(key=lambda t: -t[0])
    return out


def box_of(run):
    with open(os.path.join(os.path.dirname(run.rstrip("/")), "params.txt")) as f:
        for line in f:
            if line.split("=")[0].strip() == "box":
                return float(line.split("=")[1].split("#")[0])
    return 10.0


def linear_delta2(k, z, growth):
    return k ** 3 * pc.power_z0(k) * growth(1 / (1 + z)) ** 2 / (2 * np.pi ** 2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--out")
    ap.add_argument("--kmax-growth", type=float, default=3.0,
                    help="k [h/Mpc] up to which bins enter the growth check")
    a = ap.parse_args()
    growth = pc.Growth(1.0)
    runs = []
    for r in a.runs:
        path, _, label = r.partition(":")
        runs.append((path, label or path, read_run(path)))

    for path, label, pks in runs:
        z0, d0 = pks[0]
        k = d0[:, 0]
        sel = (k > 0) & (k <= a.kmax_growth)
        print(f"\n== {label} (box {box_of(path):g} Mpc/h); growth of the DM Delta^2 for "
              f"k <= {a.kmax_growth:g} h/Mpc ({sel.sum()} bins) relative to z = {z0:.2f}")
        print("|   z   | measured / linear [D(z)/D(z0)]^2, per k bin | mean |")
        print("|---|---|---|")
        for z, d in pks:
            lin = (growth(1 / (1 + z)) / growth(1 / (1 + z0))) ** 2
            r = d[sel, 3] / d0[sel, 3] / lin
            print(f"| {z:5.2f} | " + " ".join(f"{x:.3f}" for x in r) + f" | {r.mean():.3f} |")
        for z, d in (pks[0], pks[-1]):
            kk = d[:, 0]
            m = (kk > 0) & (d[:, 3] > 0)
            ratio = d[m, 3] / linear_delta2(kk[m], z, growth)
            idx = [int(np.argmin(np.abs(kk[m] - q))) for q in (0.63, 1.26, 2.5, 5, 10, 20)]
            print(f"z = {z:.2f}: DM Delta^2 / linear at k = " +
                  ", ".join(f"{kk[m][i]:.2f}: {ratio[i]:.2f}" for i in idx))
        for z, d in pks:
            kk = d[:, 0]
            idx = [int(np.argmin(np.abs(kk - q))) for q in (2, 5, 10, 20)]
            if z in (pks[0][0], pks[len(pks) // 2][0], pks[-1][0]):
                print(f"z = {z:.2f}: gas/DM Delta^2 at k = " +
                      ", ".join(f"{kk[i]:.1f}: {d[i, 5] / d[i, 1]:.3f}" for i in idx))

    if a.out:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        colors = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
        ink, muted, grid = "#2b2b2b", "#6b6b6b", "#e4e4e0"
        fig, axs = plt.subplots(1, 3, figsize=(15, 4.4), constrained_layout=True)
        kl = np.logspace(-0.5, 2.2, 200)
        for (path, label, pks), c in zip(runs, colors):
            for (z, d), ls in zip((pks[0], pks[-1]), ("-", "--")):
                m = d[:, 3] > 0
                axs[0].loglog(d[m, 0], d[m, 3], color=c, lw=2, ls=ls, label=f"{label}, z = {z:.1f}")
            z0, d0 = pks[0]
            sel = (d0[:, 0] > 0) & (d0[:, 0] <= a.kmax_growth)
            zz = np.array([z for z, _ in pks])
            meas = np.array([np.mean(d[sel, 3] / d0[sel, 3]) for _, d in pks])
            axs[1].plot(zz, meas, "o", color=c, ms=8, mec="white", mew=2, label=label, zorder=3)
            zl, dl = pks[-1]
            m = (dl[:, 1] > 0) & (dl[:, 0] > 0)
            axs[2].semilogx(dl[m, 0], dl[m, 5] / dl[m, 1], color=c, lw=2, label=f"{label}, z = {zl:.1f}")
        for z, ls in ((pks[0][0], "-"), (pks[-1][0], "--")):
            axs[0].loglog(kl, linear_delta2(kl, z, growth), color=muted, lw=1.2, ls=ls,
                          label=f"linear, z = {z:.1f}")
        z0 = runs[0][2][0][0]
        zs = np.linspace(5, z0, 100)
        axs[1].plot(zs, (np.array([growth(1 / (1 + z)) for z in zs]) / growth(1 / (1 + z0))) ** 2,
                    color=muted, lw=1.5, label=r"linear $[D(z)/D(z_0)]^2$")
        axs[0].set_xlabel(r"k [h/Mpc]")
        axs[0].set_ylabel(r"DM $\Delta^2(k)$ (CIC deconvolved)")
        axs[1].set_xlabel("redshift z")
        axs[1].set_ylabel(rf"$\Delta^2(z)/\Delta^2(z_0)$, k ≤ {a.kmax_growth:g} h/Mpc")
        axs[1].invert_xaxis()
        axs[2].set_xlabel(r"k [h/Mpc]")
        axs[2].set_ylabel(r"gas / DM $\Delta^2(k)$ (raw)")
        axs[2].axhline(1, color=muted, lw=0.8)
        for ax in axs:
            ax.grid(color=grid, lw=0.8)
            for s in ("top", "right"):
                ax.spines[s].set_visible(False)
            ax.legend(frameon=False, fontsize=7.5, labelcolor=ink)
        fig.savefig(a.out, dpi=130)
        print("wrote", a.out)


if __name__ == "__main__":
    main()
