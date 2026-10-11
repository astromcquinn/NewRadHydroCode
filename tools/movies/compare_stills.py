#!/usr/bin/env python3
"""Side-by-side neutral-fraction slabs of several runs at the same redshifts.

    compare_stills.py OUT.png RUN1[:label] RUN2[:label] ... [--z 11,10,9,8,7]
                      [--field fHI] [--axis 2] [--slab 2] [--pos K]

Rows are redshifts (the nearest movie frame of each run), columns are runs.
All runs must have the same N. The slab position defaults to the plane that
cuts the most ionised gas in the last run at the middle redshift (within the
stored slab for slab frames, movie_slab). Uses the frame reader and colour
scales of tools/movies/stills.py.
"""
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stills import (INK, MAPS, MUTED, check_slab, frame_slab, most_ionised, panel,  # noqa: E402
                    read_frames)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--z", default="11,10,9,8,7")
    ap.add_argument("--field", default="fHI")
    ap.add_argument("--axis", type=int, default=2)
    ap.add_argument("--slab", type=int, default=2)
    ap.add_argument("--pos", type=int, default=-1)
    ap.add_argument("--box", type=float, default=10.0)
    a = ap.parse_args()
    runs = []
    for r in a.runs:
        path, _, label = r.partition(":")
        runs.append((path, label or os.path.basename(os.path.dirname(path.rstrip("/")))))
    zs = [float(z) for z in a.z.split(",")]
    frames = {p: read_frames(p) for p, _ in runs}

    def pick(p, z):
        fr = frames[p]
        return fr[int(np.argmin(np.abs(fr["z"] - z)))]

    if a.pos >= 0:
        k0 = a.pos
    else:
        p = runs[-1][0]
        k0 = most_ionised(p, pick(p, zs[len(zs) // 2])["frame"], a.axis, a.slab)
    for z in zs:
        for p, _ in runs:
            check_slab(p, pick(p, z)["frame"], a.field, a.axis, a.slab, k0)

    plt.rcParams.update({"font.size": 10, "text.color": INK, "axes.labelcolor": INK})
    fig, axs = plt.subplots(len(zs), len(runs), figsize=(3.4 * len(runs) + 0.8, 3.3 * len(zs) + 0.8),
                            constrained_layout=True, squeeze=False)
    im = None
    for r, z in enumerate(zs):
        for c, (p, label) in enumerate(runs):
            f = pick(p, z)
            im = panel(axs[r, c], frame_slab(p, f["frame"], a.field, a.axis, a.slab, k0), a.field, a.box)
            axs[r, c].set_title(f"z = {f['z']:.2f}, " + rf"$x_{{\rm HII,V}}$ = {f['xv']:.2f}",
                                fontsize=9, color=MUTED)
            if r == 0:
                axs[r, c].set_title(f"{label}\n" + axs[r, c].get_title(), fontsize=10, color=INK)
        axs[r, 0].set_ylabel(f"z ≈ {z:g}", fontsize=11)
    fig.colorbar(im, ax=axs, location="bottom", shrink=0.6, pad=0.01).set_label(MAPS[a.field]["label"])
    fig.supxlabel(f"slab of {a.slab} cells at cell {k0} along {'xyz'[a.axis]}; "
                  f"box {a.box:g} Mpc/h (comoving)", fontsize=8, color=MUTED)
    fig.savefig(a.out, dpi=110)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
