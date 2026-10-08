#!/usr/bin/env python3
"""Reionization movie from apps/radhydro movie frames (movie_dz output).

    make_movie.py RUN_DIR [--out RUN_DIR/movie/reion.mp4] [--fps 15]
                  [--slice-axis z] [--slice-index N/2] [--thick 1]
                  [--jobs 4] [--frames-only] [--dpi 100]

Needs movie_fields = fHI,temp,rhog in the run (rhog optional; without it
the density panel is left out). Each frame shows slices of x_HII, T and the
gas overdensity at the frame's redshift, and a panel with the volume- and
mass-weighted x_HII(z) and tau_es(z) up to that redshift (Planck 2018 band).
The PNGs go to RUN_DIR/movie/png/, then ffmpeg (libx264, yuv420p) makes the
MP4. Colour ranges are fixed over the movie so frames compare.
"""
import argparse
import os
import shutil
import subprocess
import sys
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "analysis"))
import reion  # noqa: E402

A = None  # options, set in main (inherited by the worker processes)
HIST = None


def load(run, frame, field):
    p = os.path.join(run, "movie", f"frame_{frame:05d}_{field}.bin")
    if not os.path.exists(p):
        return None
    v = np.fromfile(p, "<f4")
    n = round(len(v) ** (1 / 3))
    if n ** 3 != len(v):
        raise ValueError(f"{p}: not a cube")
    return v.reshape(n, n, n)  # [k][j][i], i fastest (Fortran order)


def take_slice(cube, axis, index, thick):
    n = cube.shape[0]
    idx = [(index + d) % n for d in range(thick)]
    ax = {"z": 0, "y": 1, "x": 2}[axis]
    return np.take(cube, idx, axis=ax).mean(axis=ax)


def render(k):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fr, z = int(HIST["frame"][k]), HIST["z"][k]
    fHI = load(A.run, fr, "fHI")
    temp = load(A.run, fr, "temp")
    rho = load(A.run, fr, "rhog")
    n = fHI.shape[0]
    si = A.slice_index if A.slice_index >= 0 else n // 2
    panels = [("x$_{\\rm HII}$", take_slice(1 - fHI, A.slice_axis, si, A.thick), "Blues_r", (0, 1), "lin")]
    if temp is not None:
        panels.append(("log$_{10}$ T [K]", np.log10(np.maximum(take_slice(temp, A.slice_axis, si, A.thick), 1)),
                       "inferno", (A.tmin, A.tmax), "lin"))
    if rho is not None:
        d = take_slice(rho, A.slice_axis, si, A.thick) / rho.mean()
        panels.append(("log$_{10}$ (1+δ) gas", np.log10(np.maximum(d, 1e-3)), "Greys_r", (-1, 1.5), "lin"))

    fig = plt.figure(figsize=(16, 9), dpi=A.dpi, facecolor="black")
    npan = len(panels)
    w = 0.92 / npan
    for i, (title, img, cmap, (vmin, vmax), _) in enumerate(panels):
        ax = fig.add_axes([0.04 + i * w, 0.38, w - 0.02, 0.56])
        im = ax.imshow(img, origin="lower", cmap=cmap, vmin=vmin, vmax=vmax,
                       extent=[0, A.box, 0, A.box], interpolation="nearest")
        ax.set_title(title, color="white", fontsize=14)
        ax.set_xticks([])
        ax.set_yticks([])
        cb = fig.colorbar(im, ax=ax, orientation="horizontal", fraction=0.05, pad=0.03)
        cb.ax.tick_params(colors="white", labelsize=10)
        cb.outline.set_edgecolor("white")
    fig.text(0.04, 0.955, f"z = {z:5.2f}", color="white", fontsize=22, weight="bold")
    fig.text(0.96, 0.955, f"{A.box:g} Mpc/h, {n}³ RT cells, slice {A.slice_axis} = {si}"
             f"{'' if A.thick == 1 else f' (+{A.thick-1})'}, ĉ = {A.c_label}",
             color="0.8", fontsize=12, ha="right")

    # History panel: x_HII(z) (left axis) and tau_es(z) (separate axis box)
    zall, xV, xM, tau = HIST["z"], HIST["xV"], HIST["xM"], HIST["tau"]
    past = zall >= z - 1e-9
    zr = (A.zlo, A.zhi)
    ax1 = fig.add_axes([0.06, 0.07, 0.40, 0.24], facecolor="black")
    ax1.plot(zall, xV, color="0.35", lw=1)
    ax1.plot(zall[past], xV[past], color="#6fa8ff", lw=2.5, label="volume")
    ax1.plot(zall[past], xM[past], color="#ffffff", lw=1.5, ls="--", label="mass")
    ax1.axvline(z, color="0.6", lw=1)
    ax1.set_xlim(zr[1], zr[0])
    ax1.set_ylim(0, 1.02)
    ax1.set_xlabel("z", color="white")
    ax1.set_ylabel("⟨x$_{\\rm HII}$⟩", color="white")
    ax1.legend(frameon=False, labelcolor="white", fontsize=10, loc="upper left")
    ax2 = fig.add_axes([0.55, 0.07, 0.40, 0.24], facecolor="black")
    ax2.axhspan(reion.PLANCK_TAU - reion.PLANCK_SIG, reion.PLANCK_TAU + reion.PLANCK_SIG,
                color="#ffb347", alpha=0.25, lw=0)
    ax2.axhline(reion.PLANCK_TAU, color="#ffb347", lw=1)
    ax2.text(zr[1] - 0.2, reion.PLANCK_TAU + reion.PLANCK_SIG + 0.002, "Planck 2018",
             color="#ffb347", fontsize=10)
    ax2.plot(zall, tau, color="0.35", lw=1)
    ax2.plot(zall[past], tau[past], color="#6fa8ff", lw=2.5)
    ax2.axvline(z, color="0.6", lw=1)
    ax2.set_xlim(zr[1], zr[0])
    ax2.set_ylim(0, max(0.08, 1.1 * tau.max()))
    ax2.set_xlabel("z", color="white")
    ax2.set_ylabel("τ$_{\\rm es}$(<z)", color="white")
    for ax in (ax1, ax2):
        ax.tick_params(colors="white")
        for sp in ax.spines.values():
            sp.set_color("0.6")
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    out = os.path.join(A.png_dir, f"f{k:05d}.png")
    fig.savefig(out, facecolor="black")
    plt.close(fig)
    return out


def main():
    global A, HIST
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--out")
    ap.add_argument("--fps", type=int, default=15)
    ap.add_argument("--slice-axis", default="z", choices=["x", "y", "z"])
    ap.add_argument("--slice-index", type=int, default=-1)
    ap.add_argument("--thick", type=int, default=1, help="cells averaged across the slice")
    ap.add_argument("--box", type=float, default=10.0)
    ap.add_argument("--c-label", default="0.01c")
    ap.add_argument("--tmin", type=float, default=2.0)
    ap.add_argument("--tmax", type=float, default=4.7)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--dpi", type=int, default=100)
    ap.add_argument("--frames-only", action="store_true")
    A = ap.parse_args()
    fr, z, xV, xM = reion.read_frames(A.run)
    zz, tau = reion.tau_history(z, xM)
    HIST = {"frame": fr, "z": z, "xV": xV, "xM": xM, "tau": np.interp(z, zz, tau)}
    A.zlo, A.zhi = z.min() - 0.1, z.max() + 0.1
    A.png_dir = os.path.join(A.run, "movie", "png")
    os.makedirs(A.png_dir, exist_ok=True)
    order = np.argsort(-z)  # high z first
    HIST = {k: v[order] for k, v in HIST.items()}
    with Pool(A.jobs) as p:
        outs = p.map(render, range(len(order)))
    print(f"{len(outs)} frames in {A.png_dir}; tau_es total = {tau[-1]:.4f}")
    if A.frames_only:
        return
    if not shutil.which("ffmpeg"):
        sys.exit("ffmpeg not found; frames are in " + A.png_dir)
    out = A.out or os.path.join(A.run, "movie", "reion.mp4")
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(A.fps),
                    "-i", os.path.join(A.png_dir, "f%05d.png"), "-c:v", "libx264", "-pix_fmt",
                    "yuv420p", "-crf", "18", "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", out],
                   check=True)
    print(out)


if __name__ == "__main__":
    main()
