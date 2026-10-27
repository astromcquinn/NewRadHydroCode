#!/usr/bin/env python3
"""Side-by-side movie of several runs that differ only in the RT (rays, M1,
OTVET; docs/rt_moments.md section 6, docs/movies.md).

    compare_movie.py OUT.mp4 RUN1[:label] RUN2[:label] ... [--axis z] [--pos P]
                     [--thick 1] [--zmax Z] [--zmin Z] [--every n]
                     [--fps 12] [--jobs 8] [--stills auto] [--diff] [--diff-range F]

RUN is a run directory (<problem>/run). Columns are runs; the rows are
slices of x_HII and log10 T, and of log10 Gamma_HI (fixed range 1e-14 -
1e-11 s^-1) when every run wrote gammaHI frames (--rows to choose). A time
panel shows x_HII,V(z) of every run up to the current frame; with the
gamma row a second one shows the median Gamma_HI of the ionised cells
(movie/frames.txt). The first run sets the frames; each other run
shows its movie frame nearest in redshift (the frames are written at the
same movie_dz redshifts, so they agree to within an RT substep; the
largest mismatch is printed). All runs must have the same N and box.
Full-cube and slab frames (movie_slab) both work; with slab frames the
slice must lie in the stored slabs of every run.
--diff shows, in every column after the first, the field minus the first
run's (x_HII and rho as differences, T and Gamma_HI in dex), on a diverging
scale of +-F (--diff-range, default 0.5 for x_HII and rho, 1 dex for Gamma
and T); the first column keeps the field itself. Uses the frame index, slab
reader and colour scales of make_movie.py.
"""
import argparse
import multiprocessing as mp
import os
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import make_movie as mm  # noqa: E402  (sets the Agg backend)
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import Normalize  # noqa: E402

COLORS = ["#4cc9f0", "#f72585", "#80ed99", "#ffd166"]
ROWS = ("xHII", "temp")  # default rows; "gamma" is added when every run has gammaHI frames
DIFF_RANGE = {"xHII": 0.5, "rhog": 0.5, "temp": 1.0, "gamma": 1.0}  # --diff half-ranges


class Renderer:
    def __init__(self, runs, labels, match, args, series):
        self.runs, self.labels, self.match, self.args, self.s = runs, labels, match, args, series
        r0 = runs[0]
        N, box = r0.N, r0.box
        self.N, self.cell = N, box / N
        self.ax = mm.AXES.index(args.axis)
        self.plane = [k for k in range(3) if k != self.ax]
        self.i0 = mm.plane_start(r0, args)  # also checks slab frames (movie_slab)
        for r in runs[1:]:
            mm.plane_start(r, args, self.i0)
        self.thick = max(1, args.thick)
        self.rows = rows = args.rows
        W, H = args.size
        plt.rcParams.update({"font.size": 12, "axes.titlesize": 13})
        self.fig = fig = plt.figure(figsize=(W / args.dpi, H / args.dpi), dpi=args.dpi,
                                    facecolor="black")
        nr = len(runs)
        # square slices on the left (rows: fields, columns: runs), the time
        # panel on the right; positions in pixels, converted to figure units
        gap, x0, y0 = 0.012 * W, 0.045 * W, 0.08 * H
        side = min((0.80 * H - (len(rows) - 1) * 2.5 * gap) / len(rows), 0.62 * W / nr)
        # centre the block of rows vertically when the width limits the panels (one row, many runs)
        y0 += 0.5 * (0.80 * H - (len(rows) * side + (len(rows) - 1) * 2.5 * gap))
        px = lambda x, y, w, h: [x / W, y / H, w / W, h / H]  # noqa: E731
        self.ims = {}
        hlab, vlab = mm.AXES[self.plane[0]], mm.AXES[self.plane[1]]
        self.diff = getattr(args, "diff", False)
        for i, fld in enumerate(rows):
            st = mm.STYLE[fld]
            lo, hi = mm.fixed_range(fld, args)
            dr = args.diff_range if getattr(args, "diff_range", None) else DIFF_RANGE.get(fld, 0.5)
            yb = y0 + (len(rows) - 1 - i) * (side + 2.5 * gap)
            for j in range(nr):
                a = fig.add_axes(px(x0 + j * (side + gap), yb, side, side))
                if self.diff and j > 0:
                    im = a.imshow(np.zeros((2, 2)), origin="lower", cmap="RdBu_r",
                                  norm=Normalize(-dr, dr), interpolation="nearest",
                                  extent=(0, box, 0, box))
                else:
                    im = a.imshow(np.zeros((2, 2)), origin="lower", cmap=st["cmap"],
                                  norm=Normalize(lo, hi), interpolation="nearest",
                                  extent=(0, box, 0, box))
                self._style(a)
                step = 1.0 if box <= 5 else 2.0 if box <= 12 else 10.0 ** np.floor(np.log10(box / 2))
                ticks = [t for t in np.arange(0, box, step)]
                a.set_xticks(ticks)
                a.set_yticks(ticks)
                a.set_xticks([]) if i < len(rows) - 1 else a.set_xlabel("%s [Mpc/h]" % hlab, color="white")
                a.set_yticks([]) if j > 0 else a.set_ylabel("%s [Mpc/h]" % vlab, color="white")
                if i == 0:
                    a.set_title(labels[j], color=COLORS[j % len(COLORS)], fontsize=16)
                self.ims[fld, j] = im
            cax = fig.add_axes(px(x0 + nr * (side + gap), yb, 0.012 * W, side))
            cb = fig.colorbar(im, cax=cax)
            cb.set_label((r"$\Delta$ " + st["label"] + " (run $-$ " + labels[0] + ")") if self.diff and nr > 1
                         else st["label"], color="white")
            cb.ax.tick_params(colors="white")
            cb.outline.set_edgecolor("white")
        xt = x0 + nr * (side + gap) + (0.11 if "gamma" in rows else 0.09) * W
        # time panel(s): x_HII,V, and with the gamma row the median Gamma_HI
        # of the ionised cells below it
        tw = W - xt - 0.03 * W
        tall = len(rows) * side + (len(rows) - 1) * 2.5 * gap
        z = series["z"]
        self.glines = []
        if "gamma" in rows:
            th = 0.42 * tall
            a = self.gax = fig.add_axes(px(xt, y0 + 0.05 * tall, tw, th))
            self._style(a)
            a.set_xlim(z.max(), z.min())
            a.set_ylim(*mm.fixed_range("gamma", args))
            a.set_xlabel("redshift", color="white")
            a.set_ylabel(r"$\log_{10}\Gamma_{\rm HI}$ (median, ionised)", color="white")
            for j, lab in enumerate(labels):
                c = COLORS[j % len(COLORS)]
                a.plot(series["zfull"][j], series["gfull"][j], color=c, lw=1, alpha=0.25)
                self.glines.append(a.plot([], [], color=c, lw=2.2)[0])
            self.gvline = a.axvline(z.max(), color="white", lw=1, alpha=0.7)
            a = self.tax = fig.add_axes(px(xt, y0 + 0.58 * tall, tw, 0.42 * tall))
        else:
            a = self.tax = fig.add_axes(px(xt, y0 + 0.5 * side, tw, tall - 0.5 * side))
        self._style(a)
        a.set_xlim(z.max(), z.min())
        a.set_ylim(0, 1.02)
        a.set_xlabel("redshift", color="white")
        a.set_ylabel(r"$x_{\rm HII,V}$", color="white")
        self.lines = []
        for j, lab in enumerate(labels):
            c = COLORS[j % len(COLORS)]
            a.plot(series["zfull"][j], series["xfull"][j], color=c, lw=1, alpha=0.25)
            self.lines.append(a.plot([], [], color=c, lw=2.2, label=lab)[0])
        self.vline = a.axvline(z.max(), color="white", lw=1, alpha=0.7)
        if nr <= 4:  # with more runs the coloured column titles are the legend
            a.legend(loc="lower right", facecolor="black", edgecolor="gray", labelcolor="white")
        self.title = fig.text(0.04, 0.955, "", color="white", fontsize=17, va="center")
        # the per-run x_HII,V readout: inside the panel for a few runs, above it for many
        if nr > 4:
            self.xtext = a.text(0.0, 1.02, "", transform=a.transAxes, color="white", va="bottom",
                                fontsize=10)
        else:
            self.xtext = a.text(0.02, 0.97, "", transform=a.transAxes, color="white", va="top",
                                fontsize=12)
        self.sub = fig.text(0.95, 0.955, "", color="white", fontsize=12, va="center", ha="right")

    @staticmethod
    def _style(a):
        a.set_facecolor("black")
        a.tick_params(colors="white")
        for s in a.spines.values():
            s.set_color("white")

    def render(self, k):
        s, args = self.s, self.args
        ref = {}
        for j, run in enumerate(self.runs):
            fr = run.frames[self.match[j][k]]
            for fld in self.rows:
                path, _ = run.file_for(fr, mm.FILE_OF[fld])
                im = self.ims[fld, j]
                if path is None:
                    im.set_data(np.full((2, 2), np.nan))
                    continue
                data, _ = mm.to_plot(fld, mm.read_slab(path, run.N, self.ax, self.i0, self.thick), args)
                if self.diff:
                    if j == 0:
                        ref[fld] = data
                    elif fld in ref:
                        data = data - ref[fld]
                im.set_data(data.T)
        z = s["z"][k]
        self.title.set_text("z = %.2f    t = %.0f Myr" % (z, s["t"][k]))
        if len(self.labels) > 4:
            self.xtext.set_text(r"$x_{\rm HII,V}$ = %.3f (%s)" % (s["xV"][0][k], self.labels[0]))
        else:
            self.xtext.set_text("\n".join(r"%s: $x_{\rm HII,V}$ = %.3f" % (lab, s["xV"][j][k])
                                          for j, lab in enumerate(self.labels)))
        proj = "slice" if self.thick == 1 else "projection over %d cells" % self.thick
        self.sub.set_text("%g Mpc/h, %d$^3$ — %s %s = %.2f Mpc/h" % (
            self.runs[0].box, self.N, proj, mm.AXES[self.ax], (self.i0 + 0.5 * self.thick) * self.cell))
        for j, ln in enumerate(self.lines):
            ln.set_data(s["z"][:k + 1], s["xV"][j][:k + 1])
        self.vline.set_xdata([z, z])
        for j, ln in enumerate(self.glines):
            ln.set_data(s["z"][:k + 1], s["g"][j][:k + 1])
        if self.glines:
            self.gvline.set_xdata([z, z])
        self.fig.canvas.draw()
        return np.asarray(self.fig.canvas.buffer_rgba())[..., :3].tobytes()


_R = None


def _init(*a):
    global _R
    _R = Renderer(*a)


def _render(k):
    return k, _R.render(k)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("out")
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--axis", choices=list(mm.AXES), default="z")
    ap.add_argument("--pos", type=float)
    ap.add_argument("--thick", type=int, default=1)
    ap.add_argument("--zmax", type=float)
    ap.add_argument("--zmin", type=float)
    ap.add_argument("--every", type=int, default=1)
    ap.add_argument("--range-xHII", type=float, nargs=2)
    ap.add_argument("--range-temp", type=float, nargs=2)
    ap.add_argument("--range-gamma", type=float, nargs=2, help="log10 Gamma_HI range")
    ap.add_argument("--rows", default="auto", help="comma list of xHII, temp, gamma, rhog; auto: "
                    "xHII,temp and gamma when every run has gammaHI frames")
    ap.add_argument("--size", type=int, nargs=2, default=[1920, 1080])
    ap.add_argument("--dpi", type=int, default=110)
    ap.add_argument("--fps", type=int, default=12)
    ap.add_argument("--crf", type=int, default=20)
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument("--diff", action="store_true",
                    help="columns after the first show the field minus the first run's (x_HII, rho: "
                         "differences; T, Gamma_HI: dex) on a diverging scale")
    ap.add_argument("--diff-range", type=float, default=None,
                    help="half-range of the --diff scale (default 0.5 for x_HII and rho, 1 for Gamma and T)")
    ap.add_argument("--stills", default="auto", help="auto (x_HII,V of the first run = 0.1, 0.5, "
                    "0.9, last), none, or frame indices")
    a = ap.parse_args()
    a.xhi_log = False
    runs, labels = [], []
    for r in a.runs:
        path, _, lab = r.partition(":")
        runs.append(mm.Run(path))
        labels.append(lab or os.path.basename(os.path.dirname(os.path.abspath(path))))
    if len({(r.N, r.box) for r in runs}) != 1:
        sys.exit("runs differ in N or box")
    r0 = runs[0]
    if a.rows == "auto":
        a.rows = list(ROWS) + (["gamma"] if all("gammaHI" in r.frames[-1]["files"] for r in runs)
                               else [])
    else:
        a.rows = a.rows.split(",")
        for fld in a.rows:
            if fld not in mm.STYLE or fld not in mm.FILE_OF:
                sys.exit("--rows: unknown field " + fld)
    fr = [f for f in r0.frames if (a.zmax is None or f["redshift"] <= a.zmax)
          and (a.zmin is None or f["redshift"] >= a.zmin)]
    # stop at the last redshift that every run reached
    zlast = max(min(f["redshift"] for f in r.frames) for r in runs)
    fr = [f for f in fr if f["redshift"] >= zlast - 1e-3][::max(1, a.every)]
    if not fr:
        sys.exit("no common frames")
    r0.frames = fr
    z = np.array([f["redshift"] for f in fr])
    match, xV, dzmax = [], [], 0.0
    for r in runs:
        zr = np.array([f["redshift"] for f in r.frames])
        idx = np.abs(zr[None, :] - z[:, None]).argmin(axis=1)
        dzmax = max(dzmax, np.abs(zr[idx] - z).max())
        match.append(idx)
        xV.append(np.array([r.frames[i]["xHII_V"] for i in idx]))
    lg = lambda v: np.log10(v) if v > 0 else np.nan  # noqa: E731
    series = dict(z=z, t=mm.cosmic_time_myr(z), xV=xV,
                  g=[np.array([lg(r.frames[i]["GammaHI_ion_med"]) for i in idx]) for r, idx in
                     zip(runs, match)],
                  zfull=[np.array([f["redshift"] for f in r.frames]) for r in runs],
                  xfull=[np.array([f["xHII_V"] for f in r.frames]) for r in runs],
                  gfull=[np.array([lg(f["GammaHI_ion_med"]) for f in r.frames]) for r in runs])
    # slab frames: check the slice here, not in the render workers
    i0 = mm.plane_start(r0, a)
    for r in runs[1:]:
        mm.plane_start(r, a, i0)
    print("%d frames, z = %.2f -> %.2f; largest redshift mismatch between runs %.4f"
          % (len(fr), z.max(), z.min(), dzmax))
    W, H = a.size
    a.size = [W - W % 2, H - H % 2]
    W, H = a.size
    if a.stills == "none":
        stills = set()
    elif a.stills == "auto":
        stills = {len(fr) - 1} | {int(np.argmax(xV[0] >= t)) for t in (0.1, 0.5, 0.9) if xV[0].max() >= t}
    else:
        stills = {int(x) for x in a.stills.split(",")}
    still_dir = os.path.splitext(a.out)[0] + "_stills"
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    ff = subprocess.Popen([mm.find_ffmpeg(), "-y", "-loglevel", "error", "-f", "rawvideo",
                           "-pix_fmt", "rgb24", "-s", "%dx%d" % (W, H), "-r", str(a.fps), "-i", "-",
                           "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", str(a.crf),
                           "-preset", "medium", "-movflags", "+faststart", a.out],
                          stdin=subprocess.PIPE)
    init = (runs, labels, match, a, series)
    if a.jobs > 1:
        pool = mp.get_context("fork").Pool(a.jobs, _init, init)
        it = pool.imap(_render, range(len(fr)), chunksize=1)
    else:
        _init(*init)
        it = map(_render, range(len(fr)))
    for k, rgb in it:
        ff.stdin.write(rgb)
        if k in stills:
            os.makedirs(still_dir, exist_ok=True)
            p = os.path.join(still_dir, "z%.2f.png" % z[k])
            plt.imsave(p, np.frombuffer(rgb, np.uint8).reshape(H, W, 3))
            print("wrote", p)
        if (k + 1) % 20 == 0:
            print("  frame %d/%d" % (k + 1, len(fr)), flush=True)
    if a.jobs > 1:
        pool.close()
        pool.join()
    ff.stdin.close()
    if ff.wait() != 0:
        sys.exit("ffmpeg failed")
    print("wrote %s (%.1f MB, %d frames)" % (a.out, os.path.getsize(a.out) / 1e6, len(fr)))


if __name__ == "__main__":
    main()
