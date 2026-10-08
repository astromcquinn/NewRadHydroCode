#!/usr/bin/env python3
"""Still images from the movie frames of an apps/radhydro run.

    stills.py RUNDIR OUTDIR [--z 12,9,7.5,6.5] [--axis 2] [--slab 2]

RUNDIR is the radhydro output directory with movie/frames.txt and
movie/frame_NNNNN_<field>.bin (float32 (N,N,N) rt_grid fields, Fortran
order). It needs the fields fHI, temp and rhog (set movie_fields =
fHI,temp,rhog). Writes:
  stills_grid.png      rows = redshifts; columns = gas density, neutral
                       fraction, temperature (thin slab through the box)
  stills_history.png   volume- and mass-weighted x_HII(z), frames marked
  still_zXX.png        one three-panel image per requested redshift
"""
import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# Sequential, perceptually uniform maps (no rainbow); fixed scales across
# frames so images can be compared and later strung into a movie.
MAPS = {
    "rhog": dict(cmap="bone", norm=LogNorm(0.1, 30), label=r"gas density  $\rho_b/\langle\rho_b\rangle$"),
    "fHI": dict(cmap="Blues_r", norm=LogNorm(1e-5, 1), label=r"neutral fraction  $x_{\rm HI}$"),
    "temp": dict(cmap="magma", norm=LogNorm(1e2, 1e5), label=r"temperature  $T$ [K]"),
}
INK = "#2b2b2b"
MUTED = "#6b6b6b"


def read_frames(run):
    rows = []
    with open(os.path.join(run, "movie", "frames.txt")) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = line.split()
            rows.append((int(p[0]), float(p[3]), float(p[4]), float(p[5])))
    return np.array(rows, dtype=[("frame", int), ("z", float), ("xv", float), ("xm", float)])


def load(run, frame, field):
    a = np.fromfile(os.path.join(run, "movie", f"frame_{frame:05d}_{field}.bin"), "<f4")
    n = round(a.size ** (1 / 3))
    return a.reshape((n, n, n), order="F")


def slab(a, axis, width):
    n = a.shape[axis]
    k0 = n // 2 - width // 2
    s = np.take(a, range(k0, k0 + width), axis=axis).mean(axis=axis)
    return s.T  # x right, y up


def params(run):
    p = {}
    with open(os.path.join(os.path.dirname(run.rstrip("/")), "params.txt")) as f:
        for line in f:
            if "=" in line and not line.lstrip().startswith("#"):
                k, v = line.split("=", 1)
                p[k.strip()] = v.split("#")[0].strip()
    return p


def panel(ax, img, field, box):
    m = MAPS[field]
    im = ax.imshow(np.clip(img, m["norm"].vmin, m["norm"].vmax), origin="lower", cmap=m["cmap"],
                   norm=m["norm"], extent=(0, box, 0, box), interpolation="nearest")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    return im


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("out")
    ap.add_argument("--z", default="12,9,7.5,6.5")
    ap.add_argument("--axis", type=int, default=2)
    ap.add_argument("--slab", type=int, default=2)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    fr = read_frames(a.run)
    try:
        box = float(params(a.run).get("box", 10))
    except OSError:
        box = 10.0
    want = [float(z) for z in a.z.split(",")]
    picks = []
    for z in want:
        i = int(np.argmin(np.abs(fr["z"] - z)))
        if fr[i]["frame"] not in [p["frame"] for p in picks]:
            picks.append(fr[i])
    fields = ["rhog", "fHI", "temp"]
    plt.rcParams.update({"font.size": 10, "text.color": INK, "axes.labelcolor": INK,
                         "xtick.color": MUTED, "ytick.color": MUTED})

    # grid: rows = redshift, columns = field
    fig, axs = plt.subplots(len(picks), 3, figsize=(10.5, 3.3 * len(picks) + 0.6),
                            constrained_layout=True, squeeze=False)
    ims = {}
    for r, p in enumerate(picks):
        for c, fld in enumerate(fields):
            img = slab(load(a.run, p["frame"], fld), a.axis, a.slab)
            ims[fld] = panel(axs[r, c], img, fld, box)
        axs[r, 0].set_ylabel(f"z = {p['z']:.2f}\n" + rf"$\langle x_{{\rm HII}}\rangle_V$ = {p['xv']:.2f}",
                             fontsize=11)
    for c, fld in enumerate(fields):
        cb = fig.colorbar(ims[fld], ax=axs[:, c], location="top", shrink=0.85, pad=0.01)
        cb.set_label(MAPS[fld]["label"])
    axs[-1, 0].text(0.02, -0.06, f"slab of {a.slab} cells through the box centre; "
                    f"box {box:g} Mpc/h (comoving)", transform=axs[-1, 0].transAxes,
                    fontsize=8, color=MUTED, va="top")
    fig.savefig(os.path.join(a.out, "stills_grid.png"), dpi=110)
    plt.close(fig)

    # one image per redshift
    for p in picks:
        fig, axs = plt.subplots(1, 3, figsize=(12, 4.4), constrained_layout=True)
        for c, fld in enumerate(fields):
            im = panel(axs[c], slab(load(a.run, p["frame"], fld), a.axis, a.slab), fld, box)
            fig.colorbar(im, ax=axs[c], location="bottom", shrink=0.85, pad=0.02).set_label(MAPS[fld]["label"])
        fig.suptitle(f"z = {p['z']:.2f}   volume-weighted ionised fraction {p['xv']:.2f}   "
                     f"({box:g} Mpc/h box)", color=INK)
        fig.savefig(os.path.join(a.out, f"still_z{p['z']:.2f}.png"), dpi=110)
        plt.close(fig)

    # ionisation history: one axis, two series, legend + direct labels
    fig, ax = plt.subplots(figsize=(6.4, 4.0), constrained_layout=True)
    c1, c2 = "#2a6fdb", "#e0832b"
    ax.plot(fr["z"], fr["xv"], color=c1, lw=2, label="volume-weighted")
    ax.plot(fr["z"], fr["xm"], color=c2, lw=2, label="mass-weighted")
    for p in picks:
        ax.plot(p["z"], p["xv"], "o", ms=8, color=c1, mec="white", mew=2, zorder=3)
    ax.set_xlabel("redshift z")
    ax.set_ylabel(r"ionised fraction $\langle x_{\rm HII}\rangle$")
    ax.set_ylim(0, 1.02)
    ax.invert_xaxis()
    ax.grid(color="#e6e6e6", lw=0.8)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(frameon=False, loc="upper left")
    ax.set_title("Ionisation history (markers: frames shown in the stills)", fontsize=10, color=INK)
    fig.savefig(os.path.join(a.out, "stills_history.png"), dpi=110)
    plt.close(fig)


if __name__ == "__main__":
    main()
