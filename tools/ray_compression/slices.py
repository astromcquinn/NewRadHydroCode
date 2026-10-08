#!/usr/bin/env python3
"""Slice images of the lossy-precision runs (docs/ray_compression.md §6).

usage: slices.py PRECDIR PROBLEM [--step N] [--axis-index k] [--out PNG]

Rows: baseline, then each variant. Columns: x_HII, T, and the variant minus
the baseline for x_HII (and the ratio for T), on one z-slice (through the
brightest source by default: the slice with the most ionised cells).
"""
import argparse
from pathlib import Path

import numpy as np


def load(d, step, N):
    s = Path(d) / ("step_%04d" % step)
    f = lambda n: np.fromfile(s / (n + ".bin"), dtype=np.float32).reshape(N, N, N)
    return 1 - f("fHI"), f("temp")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("precdir")
    ap.add_argument("problem")
    ap.add_argument("--variants", nargs="+", default=["d1", "d4", "r3", "r6", "struct_o10_d4"])
    ap.add_argument("--step", type=int, default=0)
    ap.add_argument("--N", type=int, default=32)
    ap.add_argument("--k", type=int, default=-1)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    pdir = Path(args.precdir) / args.problem
    step = args.step or len(list((pdir / "base").glob("step_*")))
    N = args.N
    xb, Tb = load(pdir / "base", step, N)
    # array index order: file is Fortran (i fastest) -> numpy [k, j, i]
    k = args.k if args.k >= 0 else int(np.argmax(((xb > 0.5).sum(axis=(1, 2)))))
    names = {"base": "baseline", "d1": "1 bit, dithered", "d2": "2 bits, dithered", "d4": "4 bits, dithered",
             "d8": "8 bits, dithered", "r3": "3 bits, round to nearest", "r6": "6 bits, round to nearest",
             "struct": "R-struct round trip", "struct_o10_d4": "R-struct, 10-bit offsets, 4-bit dithered",
             "nrmf32": "nray_max_factor 32"}
    rows = ["base"] + args.variants
    fig, axes = plt.subplots(len(rows), 4, figsize=(10.5, 2.35 * len(rows)))
    for r, v in enumerate(rows):
        x, T = load(pdir / v, step, N)
        a = axes[r]
        im0 = a[0].imshow(x[k], origin="lower", cmap="Blues", vmin=0, vmax=1)
        im1 = a[1].imshow(np.log10(T[k]), origin="lower", cmap="magma", vmin=2, vmax=4.6)
        dx = x[k] - xb[k]
        im2 = a[2].imshow(dx, origin="lower", cmap="RdBu_r", vmin=-0.1, vmax=0.1)
        rt = T[k] / Tb[k] - 1
        im3 = a[3].imshow(rt, origin="lower", cmap="RdBu_r", vmin=-0.1, vmax=0.1)
        a[0].set_ylabel(names.get(v, v), fontsize=8)
        if v == "base":
            a[2].text(0.5, 0.5, "(reference)", ha="center", va="center", transform=a[2].transAxes, fontsize=8)
            a[3].text(0.5, 0.5, "(reference)", ha="center", va="center", transform=a[3].transAxes, fontsize=8)
        else:
            a[2].set_title("max |Δx| = %.3f, mean Δx = %+.1e" % (np.abs(dx).max(), (x - xb).mean()), fontsize=7)
            a[3].set_title("rms in ionised cells = %.3f" % np.sqrt(np.mean(rt[xb[k] > 0.5] ** 2)) if (xb[k] > 0.5).any()
                           else "", fontsize=7)
        for ax in a:
            ax.set_xticks([])
            ax.set_yticks([])
        if r == 0:
            a[0].set_title("x_HII", fontsize=9)
            a[1].set_title("log10 T [K]", fontsize=9)
        ims = (im0, im1, im2, im3)
    fig.suptitle("%s, hydro step %d, slice k = %d of %d" % (args.problem, step, k + 1, N), fontsize=10)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    for c, lab in enumerate(("x_HII", "log10 T", "x_HII − baseline", "T / baseline − 1")):
        pos = axes[-1, c].get_position()
        cax = fig.add_axes([pos.x0, 0.045, pos.width, 0.008])
        fig.colorbar(ims[c], cax=cax, orientation="horizontal").set_label(lab, fontsize=8)
    fig.savefig(args.out, dpi=130)


if __name__ == "__main__":
    main()
