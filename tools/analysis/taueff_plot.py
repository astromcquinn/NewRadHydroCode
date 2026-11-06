#!/usr/bin/env python3
"""Figure 7 of the RT-closure paper: the distribution of the Lya effective
optical depth over 50 Mpc/h segments in every run of the 512^3 suite
(docs/taueff_21cmfast.md).

usage: taueff_plot.py JSON:label [JSON:label ...] --outdir DIR
                      [--zcdf 5.98,5.49,4.99] [--tol 0.06] [--seg seg|box]

JSON is the output of igm_stats.py --taueff L (lya.segments of every
snapshot). The first file is the reference (the case-B rays); colours and
line styles follow rt_suite_figures.py. Every run's snapshot nearest in
redshift to the reference's is paired with it (within --tol).

tau_eff = -ln <F>_segment; "raw" is the run's own Lya optical depth, "A"
the depth rescaled by the factor A that matches the Bosman et al. 2022 mean
flux at the snapshot's z (igm_stats.py), so the rescaled distributions
differ only in their shape. A segment with <F> = 0 in double precision
(tau_eff > 700) is kept at tau_eff = inf (it counts in the percentiles).

Writes to DIR:
  taueff.png    top: the cumulative distribution of the rescaled tau_eff at
                the --zcdf redshifts; bottom: the 16/50/84 % points of the
                rescaled tau_eff against z, their width (84 - 16) against z,
                and the raw median with -ln of the Bosman+22 mean flux
  taueff.txt    every number: per run and snapshot, z, A, <F>, the 5, 16,
                50, 84 and 95 % points of tau_eff (raw and rescaled) for the
                box-length and the --seg segments, and the CDF on a grid
"""
import argparse
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rt_suite_figures import COLORS, STYLES  # noqa: E402
from igm_stats import read_obs, DATA_DEFAULT  # noqa: E402

QS = (5, 16, 50, 84, 95)


def taueff(F):
    F = np.asarray(F, np.float64)
    with np.errstate(divide="ignore"):
        return np.where(F > 0, -np.log(np.maximum(F, 1e-300)), np.inf)


def load(spec):
    path, _, label = spec.partition(":")
    d = json.load(open(path))
    snaps = []
    for s in d["snapshots"]:
        S = s["lya"].get("segments")
        if S is None:
            continue
        e = {"z": s["z"], "xV": s.get("xHII_V"), "A": S.get("A", float("nan")),
             "F_mean": s["lya"]["F_mean"], "F_obs": s["lya"].get("F_obs", float("nan")),
             "seg": S["seg_mpch"], "box": S["box_mpch"]}
        for k in ("F_box_raw", "F_box_A", "F_seg_raw", "F_seg_A"):
            if k in S:
                e[k] = np.asarray(S[k])
        snaps.append(e)
    snaps.sort(key=lambda e: -e["z"])
    return label or d["name"], snaps


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--zcdf", default="5.98,5.49,4.99")
    ap.add_argument("--tol", type=float, default=0.06)
    ap.add_argument("--seg", default="seg", choices=("seg", "box"),
                    help="segments of the figure: the joined (--taueff L) or the box-length ones")
    ap.add_argument("--data", default=DATA_DEFAULT)
    ap.add_argument("--dpi", type=int, default=200)
    a = ap.parse_args()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(a.outdir, exist_ok=True)
    runs = [load(s) for s in a.runs]
    ref = runs[0][1]
    obs = read_obs(a.data)
    # pair every run's snapshots with the reference's
    pairs = []  # per reference snapshot: list per run of entry or None
    for e0 in ref:
        row = []
        for _, sn in runs:
            m = min(sn, key=lambda e: abs(e["z"] - e0["z"])) if sn else None
            row.append(m if m is not None and abs(m["z"] - e0["z"]) <= a.tol else None)
        pairs.append(row)
    seglen = ref[0]["seg"] if a.seg == "seg" else ref[0]["box"]

    # ---- table
    tgrid = np.linspace(0, 10, 101)
    with open(os.path.join(a.outdir, "taueff.txt"), "w") as f:
        f.write("# tau_eff = -ln <F> over segments (igm_stats.py --taueff; docs/taueff_21cmfast.md)\n")
        f.write(f"# box-length segments: {ref[0]['box']} Mpc/h; joined segments: {ref[0]['seg']} Mpc/h "
                "(a skewer + part of another; not one sightline; the box has no modes above its "
                "length, so the scatter is a lower bound)\n")
        f.write("# A: tau rescaled so that <F> of all pixels = Bosman+22 at the snapshot z\n")
        f.write("# columns: label z x_V A <F>_raw F_obs nseg then, for seg in (box, joined) and "
                "tau in (raw, A): tau_eff at " + ", ".join(f"{q}%" for q in QS) + "\n")
        for (label, sn) in runs:
            for e in sn:
                cols = [f"{label:14s}", f"{e['z']:.4f}", f"{e['xV']:.4f}", f"{e['A']:.4g}",
                        f"{e['F_mean']:.5g}", f"{e['F_obs']:.5g}", str(len(e["F_box_raw"]))]
                for sg in ("box", "seg"):
                    for t in ("raw", "A"):
                        k = f"F_{sg}_{t}"
                        if k in e:
                            cols += [f"{v:.4g}" for v in np.percentile(taueff(e[k]), QS)]
                        else:
                            cols += ["nan"] * len(QS)
                f.write(" ".join(cols) + "\n")
        f.write(f"\n# CDF of the rescaled tau_eff over the {seglen} Mpc/h segments on a grid of tau_eff "
                "(first line), one line per run and snapshot: label z then the CDF\n")
        f.write("tau_eff " + " ".join(f"{t:.2f}" for t in tgrid) + "\n")
        key = f"F_{a.seg}_A"
        for (label, sn) in runs:
            for e in sn:
                if key in e:
                    t = np.sort(taueff(e[key]))
                    c = np.searchsorted(t, tgrid, side="right") / len(t)
                    f.write(f"{label} {e['z']:.4f} " + " ".join(f"{v:.4f}" for v in c) + "\n")

    # ---- figure
    plt.rcParams.update({"font.size": 9, "axes.titlesize": 9, "axes.labelsize": 9,
                         "xtick.labelsize": 8, "ytick.labelsize": 8, "figure.dpi": a.dpi})
    zc = [float(x) for x in a.zcdf.split(",")]
    fig, axes = plt.subplots(2, 3, figsize=(10.0, 6.4))
    keyA, keyR = f"F_{a.seg}_A", f"F_{a.seg}_raw"
    for icol, zt in enumerate(zc):
        ax = axes[0, icol]
        ip = min(range(len(ref)), key=lambda i: abs(ref[i]["z"] - zt))
        for ir, (label, _) in enumerate(runs):
            e = pairs[ip][ir]
            if e is None or keyA not in e:
                continue
            t = np.sort(taueff(e[keyA]))
            c = np.arange(1, len(t) + 1) / len(t)
            ax.plot(t, c, color=COLORS[ir % len(COLORS)], ls=STYLES[ir % len(STYLES)],
                    lw=1.6 if ir == 0 else 1.1, label=label)
        e0 = ref[ip]
        ax.set_title(f"z = {e0['z']:.2f} ($x_V$ = {e0['xV']:.3f}), A = {e0['A']:.2f}")
        ax.set_xlabel(rf"$\tau_{{\rm eff}}$ ({seglen:g} Mpc/$h$, rescaled)")
        ax.set_ylim(0, 1)
        lo = np.nanpercentile(taueff(e0[keyA]), 0.5)
        hi = np.nanpercentile(taueff(e0[keyA]), 99.5)
        ax.set_xlim(max(0, lo - 0.4 * (hi - lo)), hi + 0.6 * (hi - lo))
    axes[0, 0].set_ylabel(r"$P(<\tau_{\rm eff})$")
    axes[0, 0].legend(frameon=False, fontsize=7, loc="lower right")
    nrun = len(runs)
    for ir, (label, sn) in enumerate(runs):
        c, ls = COLORS[ir % len(COLORS)], STYLES[ir % len(STYLES)]
        # snapshots outside the Bosman+22 range (A = nan) have no rescaled values
        sa = [e for e in sn if keyA in e]
        if sa:
            z = np.array([e["z"] for e in sa])
            dz = (ir - (nrun - 1) / 2) * 0.012  # offset so the error bars do not overlap
            pA = np.array([np.percentile(taueff(e[keyA]), (16, 50, 84)) for e in sa])
            axes[1, 0].errorbar(z + dz, pA[:, 1], yerr=[pA[:, 1] - pA[:, 0], pA[:, 2] - pA[:, 1]],
                                color=c, ls=ls, lw=1.0, marker="o", ms=2.5, capsize=1.5, label=label)
            axes[1, 1].plot(z, pA[:, 2] - pA[:, 0], color=c, ls=ls, marker="o", ms=2.5, lw=1.1)
        zr = np.array([e["z"] for e in sn])
        pR = np.array([np.percentile(taueff(e[keyR]), 50) for e in sn])
        axes[1, 2].plot(zr, pR, color=c, ls=ls, marker="o", ms=2.5, lw=1.1)
    b = obs.get("mean_flux")
    if b is not None:
        w = b["z"] <= 6.3
        tz, tv = b["z"][w], -np.log(b["v"][w])
        terr = [np.log(b["v"][w]) - np.log(b["v"][w] - b["em"][w]),
                np.log(b["v"][w] + b["ep"][w]) - np.log(b["v"][w])]
        for ax in (axes[1, 0], axes[1, 2]):
            ax.errorbar(tz, tv, yerr=[terr[1], terr[0]], color="0.45", fmt="s", ms=3, lw=0.8,
                        capsize=1.5, label=r"$-\ln\langle F\rangle$, Bosman+22", zorder=0)
    axes[1, 0].set_ylabel(r"$\tau_{\rm eff}$ 16/50/84 % (rescaled)")
    axes[1, 1].set_ylabel(r"$\tau_{\rm eff}$(84 %) $-$ $\tau_{\rm eff}$(16 %) (rescaled)")
    axes[1, 2].set_ylabel(r"median $\tau_{\rm eff}$ (raw)")
    axes[1, 0].legend(frameon=False, fontsize=6, ncol=2)
    for ax in axes[1]:
        ax.set_xlabel("z")
    for ax in axes.ravel():
        ax.grid(True, color="0.9", lw=0.6)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(a.outdir, "taueff.png"))
    plt.close(fig)
    print("wrote", a.outdir)


if __name__ == "__main__":
    main()
