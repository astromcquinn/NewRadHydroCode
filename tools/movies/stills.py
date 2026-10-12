#!/usr/bin/env python3
"""Still images from the movie frames of an apps/radhydro run.

    stills.py RUNDIR OUTDIR [--z 12,9,7.5,6.5] [--axis 2] [--slab 2]

RUNDIR is the radhydro output directory with movie/frames.txt and
movie/frame_NNNNN_<field>.bin (float32 (N,N,N) rt_grid fields, Fortran
order), or slab frames frame_NNNNN_<field>.slab (movie_slab > 0, geometry
in movie/slab.txt; tools/movies/movie_frames.py). It shows fHI and temp,
and the gas density and Gamma_HI (10^-14 - 10^-11 s^-1) when the run wrote
them (movie_fields = fHI,temp,rhog,gammaHI).
With slab frames --axis must be one of movie_slab_axes, the slab
(--pos .. --pos + --slab - 1) must lie in the stored slab, and the default
position is the most ionised --slab cells within the stored slab. Writes:
  stills_grid.png      rows = redshifts; columns = gas density, neutral
                       fraction, temperature, Gamma_HI (thin slab through
                       the box)
  stills_history.png   volume- and mass-weighted x_HII(z), frames marked
  still_zXX.png        one three-panel image per requested redshift
"""
import argparse
import os
import sys

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
    "gammaHI": dict(cmap="plasma", norm=LogNorm(1e-14, 1e-11),
                    label=r"HI photoionization rate  $\Gamma_{\rm HI}$ [s$^{-1}$]"),
}
INK = "#2b2b2b"
MUTED = "#6b6b6b"

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import movie_frames as mf  # noqa: E402


def read_frames(run):
    rows = []
    with open(os.path.join(run, "movie", "frames.txt")) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = line.split()
            rows.append((int(p[0]), float(p[3]), float(p[4]), float(p[5])))
    fr = np.array(rows, dtype=[("frame", int), ("z", float), ("xv", float), ("xm", float)])
    # frames written before the first RT call hold an all-zero rt_grid
    # (quirk F5: x_HII reads 1, the mass weighting is NaN); skip them
    return fr[np.isfinite(fr["xm"])]


def load(run, frame, field):
    a = np.fromfile(os.path.join(run, "movie", f"frame_{frame:05d}_{field}.bin"), "<f4")
    n = round(a.size ** (1 / 3))
    return a.reshape((n, n, n), order="F")


def slab(a, axis, width, k0=None):
    n = a.shape[axis]
    if k0 is None:
        k0 = n // 2 - width // 2
    s = np.take(a, [(k0 + i) % n for i in range(width)], axis=axis).mean(axis=axis)
    return s.T  # x right, y up


def has_field(run, frame, field):
    return mf.frame_file(run, frame, field) is not None


def frame_slab(run, frame, field, axis, width, k0):
    """slab() of a frame's field, reading only the slab's cells (the same
    float32 values as slab(load(...))), from a full cube (.bin) or a slab
    frame (.slab, movie_slab)."""
    p = mf.frame_file(run, frame, field)
    if p is None:
        raise FileNotFoundError("no frame %d of %s in %s/movie" % (frame, field, run))
    return mf.read_cells(p, axis, k0, width).mean(axis=axis).T


def most_ionised(run, frame, axis, width):
    """First cell of the `width` planes along axis cutting the most ionised
    gas: over the whole box (periodic) for full cubes, within the stored
    slab for slab frames."""
    p = mf.frame_file(run, frame, "fHI")
    if p is not None and mf.is_slab(p):
        info = mf.info_of(p)
        if width > info["thick"]:
            sys.exit(f"--slab {width} is thicker than the stored slab (movie_slab = {info['thick']})")
        x = 1 - np.asarray(mf.slab_block(p, axis), dtype=np.float64)
        prof = x.sum(axis=tuple(i for i in range(3) if i != axis))
        return info["start"] + int(np.argmax(np.convolve(prof, np.ones(width), "valid")))
    x = 1 - load(run, frame, "fHI")
    other = tuple(i for i in range(3) if i != axis)
    prof = x.sum(axis=other)
    return int(np.argmax(np.convolve(np.r_[prof, prof[:width - 1]], np.ones(width), "valid")))


def check_slab(run, frame, field, axis, width, k0):
    p = mf.frame_file(run, frame, field)
    err = mf.plane_error(p, axis, k0, width) if p else None
    if err:
        sys.exit(err.replace("--pos/--thick", "--pos/--slab"))


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
    ap.add_argument("--pos", type=int, default=-1,
                    help="first cell of the slab (0-based); default: the plane cutting the "
                         "most ionised gas at the second requested redshift")
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
    # rhog and gammaHI only if the run wrote them (movie_fields)
    fields = [f for f in ("rhog", "fHI", "temp", "gammaHI") if has_field(a.run, picks[0]["frame"], f)]
    # slab frames thinner than --slab (movie_slab = 1, the default since
    # 2026-10-12): use the stored thickness
    p0 = mf.frame_file(a.run, picks[0]["frame"], "fHI")
    if p0 is not None and mf.is_slab(p0) and a.slab > mf.info_of(p0)["thick"]:
        print(f"note: --slab {a.slab} is thicker than the stored slab; using "
              f"{mf.info_of(p0)['thick']} cell(s)")
        a.slab = mf.info_of(p0)["thick"]
    if a.pos >= 0:
        k0 = a.pos
    else:
        ref = picks[min(1, len(picks) - 1)]
        k0 = most_ionised(a.run, ref["frame"], a.axis, a.slab)
    for p in picks:
        for fld in fields:
            check_slab(a.run, p["frame"], fld, a.axis, a.slab, k0)
    plt.rcParams.update({"font.size": 10, "text.color": INK, "axes.labelcolor": INK,
                         "xtick.color": MUTED, "ytick.color": MUTED})

    # grid: rows = redshift, columns = field
    fig, axs = plt.subplots(len(picks), len(fields), figsize=(3.5 * len(fields), 3.3 * len(picks) + 0.6),
                            constrained_layout=True, squeeze=False)
    ims = {}
    for r, p in enumerate(picks):
        for c, fld in enumerate(fields):
            img = frame_slab(a.run, p["frame"], fld, a.axis, a.slab, k0)
            ims[fld] = panel(axs[r, c], img, fld, box)
        axs[r, 0].set_ylabel(f"z = {p['z']:.2f}\n" + rf"$\langle x_{{\rm HII}}\rangle_V$ = {p['xv']:.2f}",
                             fontsize=11)
    for c, fld in enumerate(fields):
        cb = fig.colorbar(ims[fld], ax=axs[:, c], location="top", shrink=0.85, pad=0.01)
        cb.set_label(MAPS[fld]["label"])
    fig.supxlabel(f"slab of {a.slab} cells at cell {k0} along {'xyz'[a.axis]}; "
                  f"box {box:g} Mpc/h (comoving)", fontsize=8, color=MUTED)
    fig.savefig(os.path.join(a.out, "stills_grid.png"), dpi=110)
    plt.close(fig)

    # one image per redshift
    for p in picks:
        fig, axs = plt.subplots(1, len(fields), figsize=(4 * len(fields), 4.4), constrained_layout=True)
        for c, fld in enumerate(fields):
            im = panel(axs[c], frame_slab(a.run, p["frame"], fld, a.axis, a.slab, k0), fld, box)
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
