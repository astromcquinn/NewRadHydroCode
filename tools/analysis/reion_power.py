#!/usr/bin/env python3
"""3D power spectra of the neutral hydrogen density and the ionised fraction
at the snapshots of runs that differ in the RT (docs/rt_methods_512.md).

usage: reion_power.py RUN:label [...] --outdir DIR [--box 40] [--nbins 24]

RUN is a run directory (<problem>/run) with step_NNNN/{fHI,rhog}.bin
snapshots; every snapshot after the sources switch on is used. The first
run is the reference: each other run's snapshot nearest in redshift (within
0.2; the snapshots land where a hydro step crosses the zdata1 redshift, so
they differ between runs, and the legend gives each run's z and x_HII,V) is
compared with it.

Fields (both dimensionless; Delta^2 of the field itself, not of its
contrast, the usual 21-cm convention):
  rho_HI = x_HI (1 + delta) = x_HI rho_gas/<rho_gas>: the neutral hydrogen
           density in units of the mean hydrogen density (the 21-cm
           brightness field up to T_0(z), without velocities and spin
           temperature);
  x_HII  = 1 - x_HI.
P(k) [(Mpc/h)^3] in log bins of |k| [h/Mpc], Delta^2 = k^3 P / 2 pi^2.
Writes DIR/reion_power_3d.txt and DIR/reion_power_3d.png (rows: the two
fields and their ratio to the reference; columns: snapshots).
"""
import argparse
import glob
import os
import sys

import numpy as np
import scipy.fft as sfft

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gamma_power as gp  # noqa: E402

COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#8b5cf6", "#d6336c", "#6b7280", "#0e9aa7"]
FIELDS = ("rhoHI", "xHII")
LABEL = {"rhoHI": r"$\rho_{\rm HI}=x_{\rm HI}(1+\delta)$", "xHII": r"$x_{\rm HII}$"}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--box", type=float, default=None)
    ap.add_argument("--nbins", type=int, default=24)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    runs, labels = [], []
    for r in a.runs:
        p, _, lab = r.partition(":")
        runs.append(p)
        labels.append(lab or os.path.basename(os.path.dirname(p.rstrip("/"))))
    box = a.box or gp.box_of(runs[0])
    if not box:
        sys.exit("box unknown; use --box")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # ---- 3D spectra at the snapshots ----------------------------------------
    outs = []
    for r in runs:
        o = {}
        for d in sorted(glob.glob(os.path.join(r, "step_[0-9][0-9][0-9][0-9]"))):
            if os.path.exists(os.path.join(d, "fHI.bin")) and os.path.exists(os.path.join(d, "rhog.bin")):
                z = gp.read_z(d)
                st = reion_stats(d)
                if z is not None and st > 0.001:  # after the sources switch on
                    o[round(z, 4)] = d
        outs.append(o)
    if not outs[0]:
        return
    b3 = None
    lines = ["# 3D power spectra at the snapshots: Delta2 = k^3 P/2pi^2; box %g Mpc/h" % box]
    sets = []
    for z0 in sorted(outs[0], reverse=True):
        row = []
        for o in outs:
            zn = min(o, key=lambda z: abs(z - z0)) if o else None
            row.append((zn, o[zn]) if zn is not None and abs(zn - z0) < 0.2 else None)
        sets.append((z0, row))
    res3 = []
    for z0, row in sets:
        rr = []
        for e in row:
            if e is None:
                rr.append(None)
                continue
            x = np.fromfile(os.path.join(e[1], "fHI.bin"), dtype=np.float32)
            N = round(x.size ** (1 / 3))
            x = x.reshape(N, N, N)
            if not x.max() > 0:  # before the first RT step the grid is not filled
                rr.append(None)
                continue
            rho = np.fromfile(os.path.join(e[1], "rhog.bin"), dtype=np.float32).reshape(N, N, N)
            if b3 is None:
                b3 = gp.Binner(N, box, a.nbins)
            rh = x * (rho / np.float32(rho.mean(dtype=np.float64)))
            P = {}
            for fld, f in (("rhoHI", rh), ("xHII", np.float32(1) - x)):
                P[fld] = b3.power(sfft.rfftn(f - np.float32(f.mean(dtype=np.float64)), workers=8))
            rr.append((e[0], float(1 - x.mean(dtype=np.float64)), P))
            del x, rho, rh
        res3.append((z0, rr))
    res3 = [r_ for r_ in res3 if any(r_[1])]  # drop snapshots before the RT
    kk = b3.k
    for z0, rr in res3:
        lines.append("\n# z_ref = %.4f; k [h/Mpc]: " % z0 + " ".join("%.4f" % v for v in kk))
        for lab, e in zip(labels, rr):
            if e:
                for fld in FIELDS:
                    lines.append("%-8s z = %.4f xV = %.5f %-5s " % (lab, e[0], e[1], fld) +
                                 " ".join("%.4e" % v for v in kk ** 3 * e[2][fld] / (2 * np.pi ** 2)))
    with open(os.path.join(a.outdir, "reion_power_3d.txt"), "w") as fp:
        fp.write("\n".join(lines) + "\n")
    fig, ax = plt.subplots(4, len(res3), figsize=(3.3 * len(res3) + 0.8, 10.5), squeeze=False,
                           sharex=True)
    for c, (z0, rr) in enumerate(res3):
        for i, fld in enumerate(FIELDS):
            top, rat = ax[2 * i, c], ax[2 * i + 1, c]
            for j, (lab, e) in enumerate(zip(labels, rr)):
                if e:
                    col = COLORS[j % len(COLORS)]
                    top.loglog(kk, kk ** 3 * e[2][fld] / (2 * np.pi ** 2), "-", color=col, lw=1.5,
                               label="%s (z = %.2f, $x_V$ = %.3f)" % (lab, e[0], e[1]))
                    if j > 0 and rr[0]:
                        rat.semilogx(kk, e[2][fld] / rr[0][2][fld], "-", color=col, lw=1.5)
            rat.axhline(1, color="0.6", lw=1)
            if c == 0:
                top.set_ylabel(r"$\Delta^2$ of " + LABEL[fld])
                rat.set_ylabel("ratio to " + labels[0])
            for x in (top, rat):
                x.spines[["top", "right"]].set_visible(False)
        xv0 = rr[0][1] if rr[0] else float("nan")
        ax[0, c].set_title(r"z = %.2f, $x_{\rm HII,V}$ = %.3f" % (z0, xv0), fontsize=9)
        ax[-1, c].set_xlabel(r"$k$ [$h$ Mpc$^{-1}$]")
        ax[0, c].legend(frameon=False, fontsize=6)
    fig.tight_layout()
    fig.savefig(os.path.join(a.outdir, "reion_power_3d.png"), dpi=110)
    print("wrote", os.path.join(a.outdir, "reion_power_3d.png"))


def reion_stats(d):
    for line in open(os.path.join(d, "stats.txt")):
        if line.startswith("xHII_V"):
            return float(line.split("=")[1])
    return float("nan")


if __name__ == "__main__":
    main()
