#!/usr/bin/env python3
"""3D power spectrum of the HI photoionization rate in runs that differ in
the RT method (docs/rt_methods_512.md).

usage: gamma_power.py RUN[+RUN2...]:label [...] --outdir DIR [--zmin 4.8] [--zmax 5.6]
                      [--box 40] [--tol 0.03] [--nbins 24] [--workers 8] [--clip C] [--s0]

RUN is a run directory (<problem>/run); several directories of one run can
be joined with '+' (e.g. the run and its extra field cubes,
<problem>/run+<problem>/run_cubes, tools/hyak/field_cubes.slurm). Every
step_NNNN with gammaHI.bin (or gamma.bin, first component; --field for
another one-float-per-cell field such as gammaHI_avg, the substep-averaged
rate) and fHI.bin and zmin <= z <= zmax is used. The first run is the reference: for each of its
snapshots, every other run's snapshot nearest in redshift (within --tol) is
compared.

For each snapshot: delta = Gamma_HI/<Gamma_HI>_V - 1 over the whole box,
P(k) [(Mpc/h)^3] in log bins of |k| [h/Mpc] (fundamental to Nyquist), the
dimensionless Delta^2 = k^3 P / 2 pi^2, the ratio to the reference and the
cross-correlation coefficient with the reference, r(k) = P_ab / sqrt(P_aa P_bb).
Also the ionised fraction and <Gamma_HI>. Writes DIR/gamma_pk.txt and
DIR/gamma_pk.png.

--clip C caps Gamma_HI at C <Gamma_HI>_V before delta is formed (files
gamma_pk_clipC.*): the variance of delta is dominated by the source cells and
their proximity zones (in the 128^3 tests the top 100 of 2M cells hold
35-60 % of it), which add nearly white power; clipping isolates the
fluctuations of the background itself. The txt also gives the fraction of
the variance in the top 0.01 % of cells.

P(k_max), the power in the highest k bin, bounds a white-noise floor (for
the rays, the sampling noise of the discrete rays): its share of the power
at k is P(k_max)/P(k).

Note: Gamma_HI is the value of the last RT substep. With the rays it carries
the sampling noise of the discrete rays (roughly white at high k); the
moment methods have none, so the high-k ratio rays/moments contains it.

--s0 uses instead S0 of the tensor dump (radhydro tensor_dump = 1,
docs/otvet_closure.md), the Gamma_HI-weighted radiation averaged over the
hydro step's RT substeps (files gamma_pk_s0.*): for the rays it averages the
sampling noise down (to ~2 % per cell at 256^3, from the odd and even
substeps) and so compares the fields the moment methods have. The txt also
gives the quantiles of Gamma/<Gamma>_V in ionised cells (x_HII > 0.5).
"""
import argparse
import glob
import os
import sys

import numpy as np
import scipy.fft as sfft


def read_z(d):
    for line in open(os.path.join(d, "stats.txt")):
        if line.startswith("redshift"):
            return float(line.split("=")[1])
    return None


def snapshots(spec, zmin, zmax):
    """{z: step dir} over the '+'-joined run directories."""
    out = {}
    for run in spec.split("+"):
        for d in sorted(glob.glob(os.path.join(run, "step_[0-9][0-9][0-9][0-9]"))):
            if not os.path.exists(os.path.join(d, "fHI.bin")):
                continue
            if S0_MODE:
                if not os.path.exists(os.path.join(d, "tensor.bin")):
                    continue
            elif not (os.path.exists(os.path.join(d, FIELD + ".bin"))
                      or (FIELD == "gammaHI" and os.path.exists(os.path.join(d, "gamma.bin")))):
                continue
            z = read_z(d)
            if z is not None and zmin <= z <= zmax:
                out.setdefault(round(z, 5), d)
    return out


def box_of(spec):
    run = spec.split("+")[0].rstrip("/")
    p = os.path.join(os.path.dirname(run), "params.txt")
    if os.path.exists(p):
        for line in open(p):
            k, _, v = line.partition("=")
            if k.strip() == "box":
                return float(v.split("#")[0])
    return None


S0_MODE = False
FIELD = "gammaHI"  # --field: the one-float-per-cell Gamma_HI file of each step
QS = (0.001, 0.01, 0.16, 0.5, 0.84, 0.99, 0.999)


def load_gamma(d):
    if S0_MODE:
        m = np.memmap(os.path.join(d, "tensor.bin"), dtype=np.float32, mode="r").reshape(-1, 10)
        g = np.ascontiguousarray(m[:, 0])
        N = round(g.size ** (1 / 3))
        return g.reshape(N, N, N)
    p = os.path.join(d, FIELD + ".bin")
    if os.path.exists(p):
        g = np.fromfile(p, dtype=np.float32)
        N = round(g.size ** (1 / 3))
        return g.reshape(N, N, N)
    g = np.fromfile(os.path.join(d, "gamma.bin"), dtype=np.float32)
    N = round((g.size / 3) ** (1 / 3))
    return np.ascontiguousarray(g.reshape(N, N, N, 3)[..., 0])


class Binner:
    """|k| bins of the rfftn modes of an N^3 grid, with the half-plane weights."""

    def __init__(self, N, box, nbins):
        self.N, self.box = N, box
        kf = 2 * np.pi / box
        n = np.fft.fftfreq(N, 1.0 / N)
        nz = np.arange(N // 2 + 1)
        k2 = (n[:, None, None] ** 2 + n[None, :, None] ** 2) + nz[None, None, :] ** 2
        kmag = np.sqrt(k2.astype(np.float64)).ravel() * kf
        edges = np.logspace(np.log10(kf), np.log10(kf * N / 2 * np.sqrt(3)), nbins + 1)
        edges[0] *= 0.999
        self.idx = np.digitize(kmag, edges) - 1  # -1: the k = 0 mode
        w = np.full((N, N, N // 2 + 1), 2.0, dtype=np.float64)
        w[:, :, 0] = 1.0
        if N % 2 == 0:
            w[:, :, N // 2] = 1.0
        self.w = w.ravel()
        self.nb = nbins
        ok = (self.idx >= 0) & (self.idx < nbins)
        self.ok = ok
        self.cnt = np.bincount(self.idx[ok], weights=self.w[ok], minlength=nbins)
        self.k = np.bincount(self.idx[ok], weights=(self.w * kmag)[ok], minlength=nbins) / \
            np.maximum(self.cnt, 1)
        self.norm = box ** 3 / float(N) ** 6

    def power(self, fa, fb=None):
        p = (fa * np.conj(fb)).real if fb is not None else (fa.real ** 2 + fa.imag ** 2)
        p = p.ravel().astype(np.float64)
        s = np.bincount(self.idx[self.ok], weights=(self.w * p)[self.ok], minlength=self.nb)
        return s / np.maximum(self.cnt, 1) * self.norm


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--zmin", type=float, default=4.8)
    ap.add_argument("--zmax", type=float, default=5.6)
    ap.add_argument("--box", type=float, default=None)
    ap.add_argument("--tol", type=float, default=0.03)
    ap.add_argument("--nbins", type=int, default=24)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--clip", type=float, default=0.0)
    ap.add_argument("--zlist", default="", help="comma list: only the reference snapshots nearest "
                    "these redshifts")
    ap.add_argument("--s0", action="store_true",
                    help="use the step-averaged S0 of the tensor dump instead of gammaHI.bin")
    ap.add_argument("--field", default="gammaHI",
                    help="the Gamma_HI field file of each step (default gammaHI; gammaHI_avg for "
                         "the substep-averaged rate, fields = ...,gammaHI_avg)")
    a = ap.parse_args()
    global S0_MODE, FIELD
    S0_MODE = a.s0
    FIELD = a.field
    os.makedirs(a.outdir, exist_ok=True)
    specs, labels = [], []
    for r in a.runs:
        p, _, lab = r.partition(":")
        specs.append(p)
        labels.append(lab or os.path.basename(os.path.dirname(p.split("+")[0].rstrip("/"))))
    box = a.box or box_of(specs[0])
    if not box:
        sys.exit("box unknown: no params.txt next to the first run; use --box")
    snaps = [snapshots(s, a.zmin, a.zmax) for s in specs]
    if not snaps[0]:
        sys.exit("no snapshots with Gamma_HI in the reference for %g <= z <= %g" % (a.zmin, a.zmax))
    zref = sorted(snaps[0], reverse=True)
    if a.zlist:
        zref = sorted({min(zref, key=lambda z: abs(z - float(t))) for t in a.zlist.split(",")},
                      reverse=True)
    sets = []
    for z0 in zref:
        pick = []
        for sn in snaps:
            zn = min(sn, key=lambda z: abs(z - z0)) if sn else None
            pick.append((zn, sn[zn]) if zn is not None and abs(zn - z0) <= a.tol else None)
        sets.append((z0, pick))
    tag = ("gamma_pk" + ("_s0" if a.s0 else "") + ("_clip%g" % a.clip if a.clip > 0 else "")
           + ("_sel" if a.zlist else ""))
    rep = [("# Gamma_HI capped at %g <Gamma_HI>_V\n" % a.clip if a.clip > 0 else "") +
           "# Gamma_HI power spectra: delta = Gamma/<Gamma>_V - 1; k [h/Mpc], P [(Mpc/h)^3], "
           "Delta2 = k^3 P/2pi^2; ratio and r(k) relative to the first run (%s); box %g Mpc/h"
           % (labels[0], box)]
    results = []
    binner = None
    for z0, pick in sets:
        ref_f = None
        row = []
        for j, p in enumerate(pick):
            if p is None:
                row.append(None)
                continue
            zr, d = p
            g = load_gamma(d)
            if a.clip > 0:
                np.minimum(g, np.float32(a.clip * g.mean(dtype=np.float64)), out=g)
            if binner is None or binner.N != g.shape[0]:
                binner = Binner(g.shape[0], box, a.nbins)
            x = 1 - np.fromfile(os.path.join(d, "fHI.bin"), dtype=np.float32)
            gm = g.mean(dtype=np.float64)
            xv = x.mean(dtype=np.float64)
            gi = g.ravel()[x > 0.5]
            q = np.quantile(gi / gm, QS) if gi.size else np.full(len(QS), np.nan)
            del x, gi
            delta = (g / np.float32(gm) - np.float32(1)).astype(np.float32)
            del g
            d2 = np.square(delta, dtype=np.float64).ravel()
            ntop = max(1, d2.size // 10000)
            top = float(np.partition(d2, d2.size - ntop)[-ntop:].sum() / d2.sum())
            sd = float(np.sqrt(d2.mean()))
            del d2
            f = sfft.rfftn(delta, workers=a.workers)
            del delta
            P = binner.power(f)
            if j == 0:
                ref_f, Pref = f, P
                r = np.ones_like(P)
            else:
                r = binner.power(f, ref_f) / np.sqrt(P * Pref) if ref_f is not None else np.full_like(P, np.nan)
            row.append(dict(z=zr, dir=d, gmean=gm, xV=xv, P=P, r=r, sd=sd, top=top, q=q))
            if j > 0:
                del f
        if row[0] is None:
            continue
        results.append((z0, row))
        k = binner.k
        rep.append("\n# z_ref = %.4f" % z0)
        for lab, e in zip(labels, row):
            if e:
                rep.append("#   %-10s z = %.4f  x_HII,V = %.6f  <Gamma_HI>_V = %.4e  sigma_delta = %.3f  "
                           "var in top 0.01%% cells = %.2f  P(k_max) = %.3e  %s"
                           % (lab, e["z"], e["xV"], e["gmean"], e["sd"], e["top"], e["P"][-1], e["dir"]))
        rep.append("#   quantiles %s of Gamma/<Gamma>_V in ionised cells:" % " ".join("%g" % t for t in QS))
        for lab, e in zip(labels, row):
            if e:
                rep.append("#   %-10s q " % lab + " ".join("%.4f" % v for v in e["q"]))
        hdr = "# %10s" % "k" + "".join(" %11s %11s %7s %7s" % ("P_" + l, "D2_" + l, "rat_" + l[:3], "r_" + l[:3])
                                        for l, e in zip(labels, row) if e)
        rep.append(hdr)
        for i in range(len(k)):
            s = "  %10.4f" % k[i]
            for e in row:
                if e:
                    s += " %11.4e %11.4e %7.3f %7.3f" % (e["P"][i], k[i] ** 3 * e["P"][i] / (2 * np.pi ** 2),
                                                       e["P"][i] / row[0]["P"][i], e["r"][i])
            rep.append(s)
    with open(os.path.join(a.outdir, tag + ".txt"), "w") as fp:
        fp.write("\n".join(rep) + "\n")
    print("\n".join(rep[:3]))
    print("wrote", os.path.join(a.outdir, tag + ".txt"))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    colors = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#8b5cf6", "#d6336c", "#6b7280", "#0e9aa7"]
    nz = len(results)
    fig, ax = plt.subplots(3, nz, figsize=(3.6 * nz + 0.6, 9.0), squeeze=False, sharex=True)
    k = binner.k
    for c, (z0, row) in enumerate(results):
        for j, (lab, e) in enumerate(zip(labels, row)):
            if not e:
                continue
            col = colors[j % len(colors)]
            ax[0, c].loglog(k, k ** 3 * e["P"] / (2 * np.pi ** 2), "-", color=col, lw=1.6, label=lab)
            if j > 0:
                ax[1, c].semilogx(k, e["P"] / row[0]["P"], "-", color=col, lw=1.6, label=lab)
                ax[2, c].semilogx(k, e["r"], "-", color=col, lw=1.6, label=lab)
        ax[0, c].set_title("z = %.2f  (x$_{\\rm HII,V}$ = %.4f)" % (z0, row[0]["xV"]), fontsize=10)
        ax[1, c].axhline(1, color="0.6", lw=1)
        ax[1, c].set_yscale("log")
        ax[2, c].axhline(1, color="0.6", lw=1)
        ax[2, c].set_ylim(-0.1, 1.05)
        ax[2, c].set_xlabel(r"$k$ [$h$ Mpc$^{-1}$]")
        for r_ in range(3):
            ax[r_, c].spines[["top", "right"]].set_visible(False)
    ax[0, 0].set_ylabel(r"$\Delta^2_\Gamma(k)$, $\delta = \Gamma_{\rm HI}/\langle\Gamma_{\rm HI}\rangle - 1$")
    ax[1, 0].set_ylabel(r"$P_\Gamma / P_\Gamma$(%s)" % labels[0])
    ax[2, 0].set_ylabel(r"$r(k)$ with %s" % labels[0])
    ax[0, 0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    if a.clip > 0:
        fig.suptitle(r"$\Gamma_{\rm HI}$ capped at %g$\langle\Gamma_{\rm HI}\rangle$" % a.clip, fontsize=11)
    fig.savefig(os.path.join(a.outdir, tag + ".png"), dpi=110)
    print("wrote", os.path.join(a.outdir, tag + ".png"))


if __name__ == "__main__":
    main()
