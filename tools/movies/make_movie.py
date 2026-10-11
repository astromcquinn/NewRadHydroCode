#!/usr/bin/env python3
"""Reionisation movie from an apps/radhydro run directory (docs/movies.md).

    make_movie.py RUNDIR [-o OUT.mp4] [--style cosmo|planar]
                  [--axis z] [--pos 5.0] [--thick 1]
                  [--zoom X0,X1,Y0,Y1] [--zmax 10] [--zmin 5.5]
                  [--fps 12] [--jobs 4] [--stills auto]

RUNDIR is the output directory of `radhydro PROBLEM RUNDIR`. Frames come
from the `movie_dz` frames (RUNDIR/movie/frames.txt and
movie/frame_NNNNN_<field>.bin); a field that the movie frames do not carry
(`movie_fields`) is taken from the regular field output (step_NNNN/<field>.bin)
nearest in redshift, and labelled so. Without movie frames, every
step_NNNN with the fields is one frame.

Panels (style `cosmo`): slices or thin projections (mean over --thick
cells) of x_HII, T and rho_gas/<rho_gas> through the box, and a time panel
with the volume- and mass-weighted x_HII(z) and the Thomson optical depth
tau_es accumulated so far, with a marker at the current z. Style `planar`
(plane-parallel runs): a slice containing the ray direction --ray-axis,
and the I-front position (where the plane-averaged x_HII crosses 0.5)
against cosmic time.

Memory is bounded: each frame's fields are memory-mapped and only the slab
of the slice is read; rendered frames are piped to ffmpeg as raw RGB.

Field files are float32 N^3 in Fortran order (x fastest), or, with
movie_slab > 0, slab frames (frame_NNNNN_<field>.slab, geometry in
movie/slab.txt; tools/movies/movie_frames.py): then --axis must be one of
movie_slab_axes and the cells --pos .. --pos + --thick must lie in the
stored slab (by default the plane through the box centre, as for cubes).
rhog is in units of the mean baryon density. The cosmology is the original's
(cosmo_parameters.f90): Omega_m = 0.27, Omega_b = 0.045, h = 0.7,
X = 0.76, Y = 0.24 (single-precision literals).
"""
import argparse
import glob
import multiprocessing as mp
import os
import re
import shutil
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import movie_frames as mf  # noqa: E402

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import Normalize  # noqa: E402

# cosmo_parameters.f90 (real(4) literals, as include/radray/constants.hpp)
F = lambda x: float(np.float32(x))  # noqa: E731
OMEGAM, OMEGAB, H0, XHY, YHE, TCMB = F(0.27), F(0.045), F(0.70), F(0.76), F(0.24), F(2.725)
OMEGAR = float.fromhex("0x1.6559c34a8bc79p-14")
OMEGAL = 1 - OMEGAM
# physical constants (cgs)
SIGMA_T, C_CGS, MP, MPC, YR = 6.6524587e-25, 2.99792458e10, 1.67262192e-24, 3.0856776e24, 3.15576e7
H100 = 1e7 / MPC  # 100 km/s/Mpc in 1/s
RHOC100 = 3 * H100 ** 2 / (8 * np.pi * 6.6743e-8)  # g/cm^3 for h = 1

AXES = "xyz"
# Planck 2018 (TT,TE,EE+lowE+lensing) tau = 0.054 +- 0.007, quoted from
# memory (docs/movies.md): only a visual guide
PLANCK_TAU = (0.054, 0.007)
FIELDS = ("xHII", "temp", "rhog")
FILE_OF = {"xHII": "fHI", "temp": "temp", "rhog": "rhog", "gamma": "gammaHI"}
# perceptually uniform maps, fixed scales (override with --range-*)
STYLE = {
    "xHII": dict(cmap="cividis", label=r"$x_{\rm HII}$", log=False, range=(0.0, 1.0)),
    "xHI": dict(cmap="cividis_r", label=r"$\log_{10} x_{\rm HI}$", log=True, range=(-5.0, 0.0)),
    "temp": dict(cmap="inferno", label=r"$\log_{10} T$ [K]", log=True, range=(2.5, 4.7)),
    "rhog": dict(cmap="viridis", label=r"$\log_{10} \rho_{\rm gas}/\bar\rho_{\rm gas}$",
                 log=True, range=(-1.0, 2.0)),
    # Gamma_HI of the last RT substep (movie_fields gammaHI); 0 where no photons
    "gamma": dict(cmap="plasma", label=r"$\log_{10} \Gamma_{\rm HI}$ [s$^{-1}$]", log=True,
                  range=(-14.0, -11.0)),
}


# ----------------------------------------------------------------------------
# cosmology
# ----------------------------------------------------------------------------
def trapz(y, x):
    return float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(x)))


def hubble(z):
    """H(z) in 1/s, flat LCDM with radiation as in the original."""
    a = 1.0 / (1.0 + z)
    return H0 * H100 * np.sqrt(OMEGAM / a ** 3 + OMEGAR / a ** 4 + OMEGAL)


def cosmic_time_myr(z):
    """Age of the universe at z [Myr]: t = int_z^inf dz'/((1+z')H(z'))."""
    z = np.atleast_1d(np.asarray(z, dtype=float))
    out = np.empty_like(z)
    for i, zz in enumerate(z):
        # integrate in ln a from a = 0 to a(z): t = int da / (a H)
        lna = np.linspace(np.log(1e-8), np.log(1 / (1 + zz)), 4000)
        a = np.exp(lna)
        out[i] = trapz(1.0 / hubble(1 / a - 1), lna) / YR / 1e6
    return out


def nH0():
    """Comoving hydrogen number density [1/cm^3]."""
    return OMEGAB * RHOC100 * H0 ** 2 * XHY / MP


def tau_es_cumulative(z, xHII_M):
    """tau_es accumulated from the first frame (largest z) down to each z.

    n_e = n_H x_HII (1 + Y/(4X)): helium singly ionised like hydrogen (the
    movie frames carry no helium fields), mass-weighted x_HII.
    """
    order = np.argsort(-z)
    zs, xs = z[order], xHII_M[order]
    fe = 1 + YHE / (4 * XHY)
    integrand = SIGMA_T * C_CGS * nH0() * xs * fe * (1 + zs) ** 2 / hubble(zs)
    tau = np.zeros_like(zs)
    tau[1:] = np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * (zs[:-1] - zs[1:]))
    out = np.empty_like(tau)
    out[order] = tau
    return out


def tau_es_full_below(z):
    """tau_es from 0 to z for a fully ionised IGM (HeIII below z = 3)."""
    zz = np.linspace(0, z, 2000)
    fe = 1 + np.where(zz < 3, 2, 1) * YHE / (4 * XHY)
    return trapz(SIGMA_T * C_CGS * nH0() * fe * (1 + zz) ** 2 / hubble(zz), zz)


# ----------------------------------------------------------------------------
# run directory
# ----------------------------------------------------------------------------
def read_params(rundir):
    for p in (os.path.join(rundir, "params.txt"), os.path.join(rundir, "..", "params.txt")):
        if os.path.exists(p):
            d = {}
            for line in open(p):
                line = line.split("#")[0]
                if "=" in line:
                    k, v = line.split("=", 1)
                    d[k.strip()] = v.strip()
            return d
    return {}


def read_stats(path):
    d = {}
    for line in open(path):
        if "=" in line:
            k, v = line.split("=", 1)
            d[k.strip()] = float(v)
    return d


class Run:
    """Index of the frames and field outputs of a radhydro run directory."""

    def __init__(self, rundir, box=None):
        self.dir = rundir
        params = read_params(rundir)
        self.box = box if box else float(params.get("box", "nan"))
        if not np.isfinite(self.box):
            sys.exit("box size unknown: no params.txt in RUNDIR or RUNDIR/..; use --box")
        # regular outputs: step_NNNN/stats.txt (+ fields)
        self.steps = []
        for d in sorted(glob.glob(os.path.join(rundir, "step_[0-9][0-9][0-9][0-9]"))):
            st = os.path.join(d, "stats.txt")
            if os.path.exists(st):
                s = read_stats(st)
                s["dir"] = d
                self.steps.append(s)
        # movie frames
        self.frames = []
        self.skipped = 0
        ftxt = os.path.join(rundir, "movie", "frames.txt")
        if os.path.exists(ftxt):
            for line in open(ftxt):
                if line.startswith("#") or not line.strip():
                    continue
                f = line.split()
                fr = dict(frame=int(f[0]), istep=int(f[1]), istep_rt=int(f[2]),
                          redshift=float(f[3]), xHII_V=float(f[4]), xHII_M=float(f[5]))
                # whole-box Gamma_HI columns (runs since 2026-10-11; nan before)
                for i, key in enumerate(("GammaHI_V", "GammaHI_ion_med", "GammaHI_ion_mean",
                                         "fV_GammaHI")):
                    fr[key] = float(f[6 + i]) if len(f) > 6 + i else float("nan")
                fr["files"] = {}
                for fld in ("fHI", "temp", "rhog", "fHeI", "fHeII", "gamma", "gammaHI", "ngam"):
                    p = mf.frame_file(rundir, fr["frame"], fld)
                    if p:
                        fr["files"][fld] = p
                if not np.isfinite(fr["xHII_M"]):
                    # before the first RT step the rt_grid is not set up
                    # (fHI = rhog = 0): nothing to show
                    self.skipped += 1
                elif fr["files"]:  # a frame being written may be incomplete
                    self.frames.append(fr)
            self.source = "movie"
        if not self.frames:  # fall back to the field outputs
            for s in self.steps:
                files = {fld: os.path.join(s["dir"], fld + ".bin") for fld in ("fHI", "temp", "rhog")}
                files = {k: v for k, v in files.items() if os.path.exists(v)}
                if not np.isfinite(s["xHII_M"]):
                    self.skipped += 1
                elif files:
                    self.frames.append(dict(frame=int(s["istep"]), istep=int(s["istep"]),
                                            istep_rt=0, redshift=s["redshift"],
                                            xHII_V=s["xHII_V"], xHII_M=s["xHII_M"], files=files))
            self.source = "steps"
        if not self.frames:
            sys.exit("no movie frames and no field outputs in " + rundir)
        # slab frames (movie_slab > 0): geometry from movie/slab.txt
        self.slab = mf.run_slab_info(rundir) if self.source == "movie" else None
        p0 = next(iter(self.frames[0]["files"].values()))
        if mf.is_slab(p0):
            n = self.slab["N"] if self.slab else sys.exit("%s: no movie/slab.txt" % p0)
        else:
            n = mf.cube_N(p0)
            if n is None:
                sys.exit("%s is not a float32 cube" % p0)
        self.N = n
        # field outputs per field, for fields missing from the movie frames
        self.field_steps = {}
        for fld in ("fHI", "temp", "rhog", "gammaHI"):
            self.field_steps[fld] = [(s["redshift"], os.path.join(s["dir"], fld + ".bin"), int(s["istep"]))
                                     for s in self.steps
                                     if os.path.exists(os.path.join(s["dir"], fld + ".bin"))]

    def file_for(self, fr, fld):
        """(path, note): the movie frame's file, else the nearest field output."""
        if fld in fr["files"]:
            return fr["files"][fld], ""
        cands = self.field_steps.get(fld, [])
        if not cands:
            return None, "not written"
        z, path, istep = min(cands, key=lambda c: abs(c[0] - fr["redshift"]))
        return path, "step %d, z = %.2f" % (istep, z)


def open_cube(path, N):
    return np.memmap(path, dtype="<f4", mode="r", shape=(N, N, N), order="F")


def read_slab(path, N, axis, i0, thick):
    """Mean over cells i0..i0+thick-1 (periodic) along axis; returns a 2-D
    array indexed [first remaining axis, second remaining axis]. Reads full
    cubes and slab frames (movie_frames.read_plane)."""
    return mf.read_plane(path, axis, i0, thick, N=None if mf.is_slab(path) else N)


def plane_start(run, args, i0=None):
    """First cell of the slice/projection along args.axis: i0 if given,
    else from --pos, else the box centre (N//2 for full cubes; for slab
    frames the plane centred in the stored slab). Exits if slab frames do
    not hold it."""
    N, ax = run.N, AXES.index(args.axis)
    thick = max(1, args.thick)
    if i0 is not None:
        pass
    elif args.pos is not None:
        i0 = int(np.floor(args.pos / (run.box / N))) % N
    else:
        i0 = mf.default_i0(N, thick, run.slab)
    for fr in run.frames:
        for p in fr["files"].values():
            err = mf.plane_error(p, ax, i0, thick)
            if err:
                sys.exit(err)
    return i0


# ----------------------------------------------------------------------------
# rendering
# ----------------------------------------------------------------------------
def to_plot(fld, arr, args):
    if fld == "xHII":
        if args.xhi_log:
            return np.log10(np.clip(1 - arr, 1e-12, None)), "xHI"
        return 1 - arr, "xHII"
    return np.log10(np.clip(arr, 1e-30, None)), fld


def fixed_range(key, args):
    r = getattr(args, "range_" + key, None)
    return tuple(r) if r else STYLE[key]["range"]


class Renderer:
    """Builds the figure once and updates images and lines per frame."""

    def __init__(self, run, args, series):
        self.run, self.args, self.series = run, args, series
        N, box = run.N, run.box
        self.cell = box / N
        ax = AXES.index(args.axis)
        self.plane = [k for k in range(3) if k != ax]  # remaining axes (h, v)
        self.ax = ax
        self.i0 = plane_start(run, args)
        self.thick = max(1, args.thick)
        self.zoom = args.zoom or [0, box, 0, box]
        plt.rcParams.update({"font.size": 13, "axes.titlesize": 13})
        W, H = args.size
        self.fig = fig = plt.figure(figsize=(W / args.dpi, H / args.dpi), dpi=args.dpi,
                                    facecolor="black")
        self.txt_color = "white"
        if args.style == "cosmo":
            gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 0.42], left=0.05, right=0.94,
                                  bottom=0.08, top=0.87, wspace=0.22, hspace=0.32)
            self.slice_axes = [fig.add_subplot(gs[0, k]) for k in range(3)]
            self.panel_fields = list(FIELDS)
            self.tax = fig.add_subplot(gs[1, :])
        else:
            gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.5], left=0.05, right=0.94,
                                  bottom=0.08, top=0.87, wspace=0.15, hspace=0.32)
            self.slice_axes = [fig.add_subplot(gs[0, k]) for k in range(2)]
            self.panel_fields = ["xHII", "temp"]
            self.tax = fig.add_subplot(gs[1, :])
        self.ims, self.notes = [], []
        hlab, vlab = AXES[self.plane[0]], AXES[self.plane[1]]
        for a, fld in zip(self.slice_axes, self.panel_fields):
            key = "xHI" if (fld == "xHII" and args.xhi_log) else fld
            st = STYLE[key]
            lo, hi = fixed_range(key, args)
            im = a.imshow(np.zeros((2, 2)), origin="lower", cmap=st["cmap"],
                          norm=Normalize(lo, hi), interpolation="nearest",
                          extent=(0, box, 0, box))
            a.set_xlim(self.zoom[0], self.zoom[1])
            a.set_ylim(self.zoom[2], self.zoom[3])
            a.set_xlabel("%s [Mpc/h]" % hlab, color=self.txt_color)
            if not self.ims:
                a.set_ylabel("%s [Mpc/h]" % vlab, color=self.txt_color)
            a.set_title(st["label"], color=self.txt_color, fontsize=15)
            self._style_axis(a)
            cb = fig.colorbar(im, ax=a, fraction=0.046, pad=0.02)
            cb.ax.tick_params(colors=self.txt_color)
            cb.outline.set_edgecolor(self.txt_color)
            self.ims.append(im)
            self.notes.append(a.text(0.02, 0.02, "", transform=a.transAxes, color="white",
                                     fontsize=10, bbox=dict(fc="black", alpha=0.5, lw=0)))
        if args.style == "planar":
            d = "+" if args.source_face == "low" else "-"
            self.notes_ray = self.slice_axes[0].text(
                0.98, 0.98, "rays: %s%s →" % (d, args.ray_axis), transform=self.slice_axes[0].transAxes,
                ha="right", va="top", color="white", fontsize=11,
                bbox=dict(fc="black", alpha=0.5, lw=0))
        self.title = fig.text(0.04, 0.95, "", color=self.txt_color, fontsize=17, va="center")
        self.subtitle = fig.text(0.97, 0.95, "", color=self.txt_color, fontsize=13, va="center",
                                 ha="right")
        self._setup_time_panel()

    def _style_axis(self, a):
        a.set_facecolor("black")
        a.tick_params(colors=self.txt_color)
        for s in a.spines.values():
            s.set_color(self.txt_color)

    def _setup_time_panel(self):
        a, s = self.tax, self.series
        self._style_axis(a)
        if self.args.style == "cosmo":
            z = s["z"]
            a.set_xlim(z.max(), z.min())
            a.set_ylim(0, 1.02)
            a.set_xlabel("redshift", color=self.txt_color)
            a.set_ylabel(r"$x_{\rm HII}$", color=self.txt_color)
            self.lV, = a.plot([], [], color="#4cc9f0", lw=2.2, label=r"$x_{\rm HII}$ volume-weighted")
            self.lM, = a.plot([], [], color="#f72585", lw=2.2, ls="--",
                              label=r"$x_{\rm HII}$ mass-weighted")
            self.vline = a.axvline(z.max(), color="white", lw=1, alpha=0.7)
            self.dot, = a.plot([], [], "o", color="#4cc9f0", ms=7)
            a2 = a.twinx()
            self.tax2 = a2
            self._style_axis(a2)
            tmax = max(s["tau_low"].max(), PLANCK_TAU[0] + PLANCK_TAU[1], 1e-4)
            a2.set_ylim(0, 1.08 * tmax)
            a2.set_ylabel(r"$\tau_{\rm es}$", color="#ffd166")
            a2.tick_params(axis="y", colors="#ffd166")
            a2.axhspan(PLANCK_TAU[0] - PLANCK_TAU[1], PLANCK_TAU[0] + PLANCK_TAU[1],
                       color="#ffd166", alpha=0.12, lw=0)
            a2.text(0.995, PLANCK_TAU[0] + PLANCK_TAU[1], "Planck 2018",
                    transform=a2.get_yaxis_transform(), color="#ffd166", alpha=0.8,
                    ha="right", va="bottom", fontsize=10)
            self.lT, = a2.plot([], [], color="#ffd166", lw=2,
                               label=r"$\tau_{\rm es}$ accumulated since $z_{\rm start}$")
            self.lT2, = a2.plot([], [], color="#ffd166", lw=1.6, ls=":",
                                label=r"total $\tau_{\rm es}$ if fully ionised below $z$")
            hs = [self.lV, self.lM, self.lT, self.lT2]
            # legend on the twin axis, which is drawn on top
            a2.legend(hs, [h.get_label() for h in hs], loc="upper left", facecolor="black",
                      edgecolor="gray", labelcolor="white", fontsize=11, framealpha=0.9)
            self.tlabel = self.fig.text(0.05, 0.905, "", color=self.txt_color, fontsize=13,
                                        va="center")
        else:
            t, r = s["t"], s["front"]
            a.set_xlim(t.min(), t.max())
            ok = np.isfinite(r)
            a.set_ylim(0, self.run.box)
            a.set_xlabel("cosmic time [Myr]", color=self.txt_color)
            a.set_ylabel("I-front position [Mpc/h]", color=self.txt_color)
            self.lF, = a.plot([], [], color="#4cc9f0", lw=2.2,
                              label=r"I-front ($\langle x_{\rm HII}\rangle_\perp = 0.5$)")
            self.dot, = a.plot([], [], "o", color="#4cc9f0", ms=7)
            a.legend(loc="upper left", facecolor="black", edgecolor="gray", labelcolor="white")
            self.speed = self.fig.text(0.05, 0.905, "", color=self.txt_color, fontsize=13,
                                       va="center")
            if not ok.any():
                self.speed.set_text("no front found")

    def render(self, k):
        """Draw frame k; returns the RGB bytes."""
        run, args, s = self.run, self.args, self.series
        fr = run.frames[k]
        for im, note, fld in zip(self.ims, self.notes, self.panel_fields):
            path, src = run.file_for(fr, FILE_OF[fld])
            if path is None:
                im.set_data(np.full((2, 2), np.nan))
                note.set_text("%s not written" % FILE_OF[fld])
                continue
            sl = read_slab(path, run.N, self.ax, self.i0, self.thick)
            data, _ = to_plot(fld, sl, args)
            # imshow wants [v, h] with h = first remaining axis
            im.set_data(data.T)
            note.set_text(src)
        z = fr["redshift"]
        t = s["t"][k]
        proj = ("slice" if self.thick == 1 else "projection over %d cells (%.2f Mpc/h)"
                % (self.thick, self.thick * self.cell))
        self.title.set_text("z = %.2f    t = %.0f Myr" % (z, t))
        self.subtitle.set_text("%g Mpc/h box, %d$^3$ cells — %s %s = %.2f Mpc/h"
                               % (run.box, run.N, proj, AXES[self.ax],
                                  (self.i0 + 0.5 * self.thick) * self.cell))
        if args.style == "cosmo":
            sel = slice(0, k + 1)
            self.lV.set_data(s["z"][sel], s["xV"][sel])
            self.lM.set_data(s["z"][sel], s["xM"][sel])
            self.lT.set_data(s["z"][sel], s["tau"][sel])
            self.lT2.set_data(s["z"][sel], s["tau_low"][sel])
            self.vline.set_xdata([z, z])
            self.dot.set_data([z], [s["xV"][k]])
            self.tlabel.set_text(r"$x_{\rm HII,V}$ = %.3f    $x_{\rm HII,M}$ = %.3f    "
                                 r"$\tau_{\rm es}$ since $z_{\rm start}$ = %.4f    "
                                 r"total $\tau_{\rm es}$ if ionised below z = %.4f"
                                 % (s["xV"][k], s["xM"][k], s["tau"][k], s["tau_low"][k]))
        else:
            sel = slice(0, k + 1)
            self.lF.set_data(s["t"][sel], s["front"][sel])
            self.dot.set_data([t], [s["front"][k]])
            if (k > 0 and s["t"][k] > s["t"][k - 1] and np.isfinite(s["front"][k])
                    and np.isfinite(s["front"][k - 1])):
                dt = (s["t"][k] - s["t"][k - 1]) * 1e6 * YR
                v = (s["front"][k] - s["front"][k - 1]) / H0 * MPC / (1 + z) / dt
                self.speed.set_text("front speed (proper) ≈ %.0f km/s = %.2g c" % (v / 1e5, v / C_CGS))
        self.fig.canvas.draw()
        return np.asarray(self.fig.canvas.buffer_rgba())[..., :3].tobytes()


_R = None


def _init_worker(run, args, series):
    global _R
    _R = Renderer(run, args, series)


def _render(k):
    return k, _R.render(k)


# ----------------------------------------------------------------------------
# series
# ----------------------------------------------------------------------------
def front_positions(run, args):
    """I-front position along ray_axis per frame: first cell from the
    source face where the plane-averaged x_HII drops below 0.5 (streams
    one frame at a time)."""
    ax = AXES.index(args.ray_axis)
    cell = run.box / run.N
    out = np.full(len(run.frames), np.nan)
    noted = False
    for k, fr in enumerate(run.frames):
        path, _ = run.file_for(fr, "fHI")
        if path is None:
            continue
        mean, note = mf.axis_profile(path, ax, None if mf.is_slab(path) else run.N)
        if note and not noted:
            print("front position from the slab frames: x_HII averaged over the %s, "
                  "not the whole box" % note)
            noted = True
        prof = 1 - mean
        if args.source_face == "high":
            prof = prof[::-1]
        below = np.nonzero(prof < 0.5)[0]
        if prof[0] < 0.5:
            out[k] = 0.0
        elif below.size == 0:
            out[k] = run.box
        else:
            j = below[0]  # interpolate between cells j-1 (>= 0.5) and j
            f = (prof[j - 1] - 0.5) / (prof[j - 1] - prof[j])
            out[k] = (j - 0.5 + f) * cell
    return out


def pick_stills(run, series, spec):
    n = len(run.frames)
    if spec == "none":
        return []
    if spec != "auto":
        return sorted({min(n - 1, max(0, int(x))) for x in spec.split(",")})
    xv = series["xV"]
    ks = {n - 1}
    for target in (0.1, 0.5, 0.9):
        if xv.max() >= target:
            ks.add(int(np.argmax(xv >= target)))
    return sorted(ks)


def find_ffmpeg():
    exe = shutil.which("ffmpeg")
    if exe:
        return exe
    try:
        import imageio_ffmpeg
        return imageio_ffmpeg.get_ffmpeg_exe()
    except ImportError:
        sys.exit("ffmpeg not found: `pip install --user imageio-ffmpeg` (ships a static "
                 "ffmpeg), or `module load ffmpeg`, or --no-mp4 for PNG frames only")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("rundir")
    ap.add_argument("-o", "--out", help="MP4 path (default RUNDIR/movie/<style>.mp4)")
    ap.add_argument("--style", choices=["cosmo", "planar"], default="cosmo")
    ap.add_argument("--axis", choices=list(AXES), default="z", help="slice normal")
    ap.add_argument("--pos", type=float, help="slice position along --axis [Mpc/h] (default: centre)")
    ap.add_argument("--thick", type=int, default=1, help="cells averaged (thin projection)")
    ap.add_argument("--zoom", type=lambda s: [float(x) for x in s.split(",")],
                    help="H0,H1,V0,V1 [Mpc/h]: region of the slice plane shown")
    ap.add_argument("--zmax", type=float, help="first frame: z <= zmax")
    ap.add_argument("--zmin", type=float, help="last frame: z >= zmin")
    ap.add_argument("--every", type=int, default=1, help="use every n-th frame")
    ap.add_argument("--box", type=float, help="box [Mpc/h] (default from params.txt)")
    ap.add_argument("--xhi-log", action="store_true", help="show log10 x_HI instead of x_HII")
    ap.add_argument("--range-xHII", type=float, nargs=2)
    ap.add_argument("--range-xHI", type=float, nargs=2)
    ap.add_argument("--range-temp", type=float, nargs=2, help="log10 T range")
    ap.add_argument("--range-rhog", type=float, nargs=2, help="log10 rho/<rho> range")
    ap.add_argument("--ray-axis", choices=list(AXES), default="x", help="planar: ray direction")
    ap.add_argument("--source-face", choices=["low", "high"], default="low",
                    help="planar: rays enter at the low or high face of --ray-axis")
    ap.add_argument("--size", type=int, nargs=2, default=[1920, 1080], help="pixels W H")
    ap.add_argument("--dpi", type=int, default=120)
    ap.add_argument("--fps", type=int, default=12)
    ap.add_argument("--crf", type=int, default=20, help="x264 quality (lower: better, larger)")
    ap.add_argument("--jobs", type=int, default=1, help="parallel render processes")
    ap.add_argument("--stills", default="auto",
                    help="PNG stills: auto (x_HII,V = 0.1, 0.5, 0.9 and the last frame), "
                         "comma-separated frame indices, or none")
    ap.add_argument("--still-dir", help="default: next to the MP4")
    ap.add_argument("--no-mp4", action="store_true", help="only the stills")
    args = ap.parse_args()
    if args.style == "planar" and args.axis == args.ray_axis:
        # the slice must contain the ray direction
        args.axis = [a for a in AXES if a != args.ray_axis][-1]

    run = Run(args.rundir, args.box)
    fr = run.frames
    if args.zmax is not None:
        fr = [f for f in fr if f["redshift"] <= args.zmax]
    if args.zmin is not None:
        fr = [f for f in fr if f["redshift"] >= args.zmin]
    fr = fr[::max(1, args.every)]
    if not fr:
        sys.exit("no frames in the requested redshift range")
    # tau_es over all frames before selection would also count frames before
    # zmax; accumulate from the first frame of the whole run instead
    zall = np.array([f["redshift"] for f in run.frames])
    tau_all = tau_es_cumulative(zall, np.array([f["xHII_M"] for f in run.frames]))
    tau_of = dict(zip([f["frame"] for f in run.frames], tau_all))
    run.frames = fr
    z = np.array([f["redshift"] for f in fr])
    series = dict(z=z, t=cosmic_time_myr(z),
                  xV=np.array([f["xHII_V"] for f in fr]),
                  xM=np.array([f["xHII_M"] for f in fr]),
                  tau=np.array([tau_of[f["frame"]] for f in fr]))
    series["tau_low"] = series["tau"] + np.array([tau_es_full_below(zz) for zz in z])
    plane_start(run, args)  # slab frames: check the slice here, not in the render workers
    if args.style == "planar":
        series["front"] = front_positions(run, args)
    print("%s: N = %d, box = %g Mpc/h, %d frames (%s), z = %.2f -> %.2f; %d frames before "
          "the first RT step skipped" % (args.rundir, run.N, run.box, len(fr), run.source,
                                         z.max(), z.min(), run.skipped))

    out = args.out or os.path.join(args.rundir, "movie", args.style + ".mp4")
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    still_dir = args.still_dir or os.path.join(os.path.dirname(os.path.abspath(out)),
                                               os.path.splitext(os.path.basename(out))[0] + "_stills")
    stills = pick_stills(run, series, args.stills)
    W, H = args.size
    W, H = W - W % 2, H - H % 2
    args.size = [W, H]

    ff = None
    if not args.no_mp4:
        ff = subprocess.Popen([find_ffmpeg(), "-y", "-loglevel", "error",
                               "-f", "rawvideo", "-pix_fmt", "rgb24", "-s", "%dx%d" % (W, H),
                               "-r", str(args.fps), "-i", "-",
                               "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", str(args.crf),
                               "-preset", "medium", "-movflags", "+faststart", out],
                              stdin=subprocess.PIPE)
    ks = range(len(fr)) if not args.no_mp4 else stills
    if args.jobs > 1:
        pool = mp.get_context("fork").Pool(args.jobs, _init_worker, (run, args, series))
        it = pool.imap(_render, ks, chunksize=1)  # ordered, bounded by the pool
    else:
        _init_worker(run, args, series)
        it = map(_render, ks)
    for k, rgb in it:
        if len(rgb) != W * H * 3:
            sys.exit("rendered frame has the wrong size (%d bytes for %dx%d)" % (len(rgb), W, H))
        if ff:
            ff.stdin.write(rgb)
        if k in stills:
            os.makedirs(still_dir, exist_ok=True)
            img = np.frombuffer(rgb, np.uint8).reshape(H, W, 3)
            p = os.path.join(still_dir, "frame_%05d_z%.2f.png" % (fr[k]["frame"], fr[k]["redshift"]))
            plt.imsave(p, img)
            print("wrote", p)
        if (k + 1) % 20 == 0:
            print("  frame %d/%d" % (k + 1, len(fr)), flush=True)
    if args.jobs > 1:
        pool.close()
        pool.join()
    if ff:
        ff.stdin.close()
        if ff.wait() != 0:
            sys.exit("ffmpeg failed")
        print("wrote %s (%.1f MB, %d frames at %d fps)"
              % (out, os.path.getsize(out) / 1e6, len(fr), args.fps))


if __name__ == "__main__":
    main()
