#!/usr/bin/env python3
"""Batch driver for the ray-compression study (docs/ray_compression.md).

usage:
  study.py entropy RUN_DIR... --out entropy.json       per-snapshot entropies
  study.py codec   RUN_DIR... --out codec.json [--codec build/ray_codec]
  study.py figures --entropy entropy.json --codec codec.json
                   --precision precision.json --outdir docs/figures

RUN_DIR is a directory of ray dumps (rays_s*.npy + .json). The codec is
trained on the previous snapshot of the same run (the first snapshot is
coded with a model fitted on itself, flagged "oracle").
"""
import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

# Reference categorical palette (dataviz skill), fixed order
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"


def stems(run_dir):
    return [str(p)[:-4] for p in sorted(Path(run_dir).glob("rays_s*.npy"))]


def cmd_entropy(args):
    import entropy
    res = []
    for rd in args.runs:
        for st in stems(rd):
            r = entropy.analyse(st)
            r["run"] = Path(rd).parent.name
            res.append(r)
            print("%-28s %8d rays xHII=%.3f never=%.3f" % (Path(st).name, r["nrays"], r["xHII_vol"],
                                                           r["frac_never_merged"]), flush=True)
    json.dump(res, open(args.out, "w"), indent=1)


def cmd_codec(args):
    res = []
    variants = [("lossless", []), ("m10", ["--mbits", "10"]), ("m6", ["--mbits", "6"]), ("m6_o10", ["--mbits", "6", "--offbits", "10"]),
                ("no_origin", ["--no-origin"])]
    for rd in args.runs:
        ss = stems(rd)
        for i, st in enumerate(ss):
            train = ss[i - 1] if i > 0 else st
            for name, extra in variants:
                out = subprocess.run([args.codec, train, st, "--reps", str(args.reps)] + extra,
                                     check=True, capture_output=True, text=True).stdout
                r = json.loads(out)
                r["run"] = Path(rd).parent.name
                r["variant"] = name
                r["oracle"] = i == 0
                meta = json.load(open(st + ".json"))
                r["xHII_vol"] = meta["xHII_vol"]
                r["substep"] = meta["substep"]
                res.append(r)
                print("%-12s %-20s %-10s %7d rays  struct %.2f B  entropy %.2f B  enc %.0f ns dec %.0f ns" % (
                    r["run"], Path(st).name, name, r["nrays"], r["bytes_structural"], r["bytes_entropy"],
                    r["encode_ns_per_ray"], r["decode_ns_per_ray"]), flush=True)
    json.dump(res, open(args.out, "w"), indent=1)


# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------

def style(plt):
    plt.rcParams.update({
        "font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
        "ytick.color": INK2, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
        "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
        "lines.linewidth": 2, "savefig.dpi": 150, "figure.facecolor": "#fcfcfb",
        "axes.facecolor": "#fcfcfb"})


RUN_LABEL = {"snap_uniform_1src": "1 source (uniform)", "snap_lognormal_16src": "16 sources (lognormal)",
             "snap_ln_long": "16 sources, to overlap", "full64": "64³ full run"}


def run_group(r):
    return "16 sources (lognormal)" if r["run"] in ("snap_lognormal_16src", "snap_ln_long") else RUN_LABEL.get(r["run"], r["run"])


def fig_bytes(codec, outdir, plt):
    groups = {}
    for r in codec:
        groups.setdefault(run_group(r), []).append(r)
    fig, axes = plt.subplots(1, len(groups), figsize=(3.3 * len(groups), 3.0), sharey=True)
    for ax, (g, rs) in zip(np.atleast_1d(axes), groups.items()):
        for k, (lab, key, var) in enumerate([
                ("R-struct (fixed width)", "bytes_structural", "lossless"),
                ("rANS, float32 total", "bytes_entropy", "lossless"),
                ("rANS, 10-bit total", "bytes_entropy", "m10"),
                ("rANS, 6-bit total", "bytes_entropy", "m6"),
                ("rANS, 6-bit total + 10-bit offsets", "bytes_entropy", "m6_o10")]):
            pts = sorted([(r["xHII_vol"], r[key]) for r in rs if r["variant"] == var])
            if pts:
                x, y = zip(*pts)
                ax.plot(x, y, "-o", color=C[k], ms=4, label=lab)
        ax.axhline(80, color=INK2, lw=1, ls="--")
        ax.text(0.02, 80, " original: 80 B (+16 B links)", va="bottom", fontsize=7, color=INK2,
                transform=ax.get_yaxis_transform())
        ax.set_xscale("log")
        ax.set_title(g, fontsize=9, color=INK)
        ax.set_xlabel("volume-weighted x_HII of snapshot")
        ax.set_ylim(0, 90)
    np.atleast_1d(axes)[0].set_ylabel("bytes per ray")
    np.atleast_1d(axes)[0].legend(loc="center left", fontsize=7)
    fig.tight_layout()
    fig.savefig(Path(outdir) / "raycomp_bytes.png")
    plt.close(fig)


def fig_bits_by_group(codec, outdir, plt):
    # representative snapshots at the 6-bit total / 10-bit offset precision
    pick = [("full64", 1000, "64³ run, mid reionisation"), ("full64", 3000, "64³ run, after overlap"),
            ("snap_lognormal_16src", 720, "16 sources, x_HII 0.65"), ("snap_ln_long", 1440, "16 sources, after overlap")]
    reps = []
    for run, sub, lab in pick:
        for r in codec:
            if r["run"] == run and r["substep"] == sub and r["variant"] == "m6_o10":
                reps.append((lab, r))
    names = list(reps[0][1]["bits_by_group"].keys())
    fig, ax = plt.subplots(figsize=(6.5, 2.9))
    y = np.arange(len(reps))
    left = np.zeros(len(reps))
    for k, n in enumerate(names):
        v = np.array([r["bits_by_group"][n] for _, r in reps])
        ax.barh(y, v, left=left, color=C[k % 8], height=0.6, edgecolor="#fcfcfb", linewidth=1,
                hatch=None if k < 8 else "//", label=n)
        left += v
    for i, (_, r) in enumerate(reps):
        ax.text(left[i] + 1, i, "%.1f B" % r["bytes_entropy"], va="center", fontsize=7, color=INK2)
    ax.set_yticks(y)
    ax.set_yticklabels(["%s\n(x_HII=%.2f, %d rays)" % (lab, r["xHII_vol"], r["nrays"]) for lab, r in reps],
                       fontsize=7)
    ax.set_xlabel("bits per ray (rANS, 6-bit total, 10-bit offsets)")
    ax.legend(ncol=6, fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.28))
    ax.grid(axis="y", visible=False)
    ax.set_xlim(0, left.max() * 1.12)
    fig.tight_layout()
    fig.savefig(Path(outdir) / "raycomp_bits_by_field.png")
    plt.close(fig)


def fig_precision(prec, outdir, plt):
    m = prec["metrics"]
    probs = list(m.keys())
    metrics = [("xHII_mass", "x_HII (mass) max rel. err.", 0.01),
               ("T_mass", "T (mass) max rel. err.", 0.03),
               ("front_shift", "I-front shift (cells)", 0.5),
               ("photons_ionisations", "photoionisations rel. err.", 0.01)]
    fig, axes = plt.subplots(1, len(metrics), figsize=(11, 3.6))
    bits = [1, 2, 3, 4, 6, 8, 10]
    for ax, (key, lab, tgt) in zip(axes, metrics):
        for k, p in enumerate(probs):
            v = [m[p].get("d%d" % b, {}).get(key, np.nan) for b in bits]
            ax.plot(bits, np.maximum(v, 1e-6), "-o", ms=4, color=C[k], label=p + " (dithered)")
            vr = [m[p].get("r%d" % b, {}).get(key, np.nan) for b in (3, 6)]
            ax.plot([3, 6], np.maximum(vr, 1e-6), "s", ms=6, mfc="none", color=C[k],
                    label=p + " (round to nearest)")
            ref = m[p].get("nrmf32", {}).get(key, np.nan)
            if ref > 0:
                ax.axhline(ref, color=C[k], lw=1, ls=":")
        ax.axhline(tgt, color=INK, lw=1, ls="--")
        ax.set_yscale("log")
        ax.set_xlabel("photon mantissa bits")
        ax.set_title(lab, fontsize=9)
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, ncol=3, fontsize=7, loc="lower center", bbox_to_anchor=(0.5, 0.0))
    fig.text(0.5, 0.115, "dashed: target (docs/ROADMAP.md); dotted: baseline vs nray_max_factor 32 "
             "(the algorithm's own resolution sensitivity; zero for the uniform problems)",
             ha="center", fontsize=7, color=INK2)
    fig.tight_layout(rect=(0, 0.15, 1, 1))
    fig.savefig(Path(outdir) / "raycomp_precision.png")
    plt.close(fig)


def fig_history(prec, outdir, plt):
    h = prec["hist"]
    probs = list(h.keys())
    fig, axes = plt.subplots(1, len(probs), figsize=(3.5 * len(probs), 3.1))
    for ax, p in zip(np.atleast_1d(axes), probs):
        for k, v in enumerate(["d1", "d4", "struct_o10", "r6", "r3", "nrmf32"]):
            lab = {"d1": "1 bit, dithered", "d4": "4 bits, dithered", "struct_o10": "structural, 10-bit offsets",
                   "r6": "6 bits, round to nearest", "r3": "3 bits, round to nearest",
                   "nrmf32": "nray_max_factor 32 (reference)"}[v]
            if v not in h[p]:
                continue
            z = np.array(h[p][v]["z"])
            xb = np.array(h[p][v]["xm_b"])
            x = np.array(h[p][v]["xm"])
            ok = xb > 1e-3
            ax.plot(z[ok], x[ok] / xb[ok] - 1, color=C[k], lw=1.5, label=lab)
        ax.axhspan(-0.01, 0.01, color=GRID, alpha=0.6, lw=0)
        ax.set_ylim(-0.1, 0.1)
        ax.set_yscale("symlog", linthresh=1e-3)
        ax.invert_xaxis()
        ax.set_xlabel("z")
        ax.set_title(p, fontsize=9)
    np.atleast_1d(axes)[0].set_ylabel("x_HII,mass / baseline - 1")
    h, l = np.atleast_1d(axes)[0].get_legend_handles_labels()
    fig.legend(h, l, ncol=6, fontsize=7, loc="lower center")
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    fig.savefig(Path(outdir) / "raycomp_xhii_history.png")
    plt.close(fig)


def fig_determinism(ent, outdir, plt, snapshot_stems):
    import snap
    fig, axes = plt.subplots(1, 3, figsize=(10, 2.7))
    for k, st in enumerate(snapshot_stems):
        a, meta = snap.load(st)
        d = snap.derive(a, meta)
        ex = d["never_merged"] & ~d["clamped"] & (a["src"] >= 0)
        lab = "%s (%d rays)" % (RUN_LABEL.get(Path(st).parent.parent.name, ""), ex.sum())
        res = np.abs(d["x_resid_cen"][ex]).max(1)
        axes[0].hist(np.log10(np.maximum(res, 1e-9)), bins=60, histtype="step", color=C[k], lw=1.5, label=lab)
        dev = np.hypot(d["dev_th"], d["dev_ph"])
        axes[1].hist(np.log10(np.maximum(dev[ex], 1e-9)), bins=60, histtype="step", color=C[k], lw=1.5)
        axes[1].hist(np.log10(np.maximum(dev[~d["never_merged"]], 1e-9)), bins=60, histtype="step",
                     color=C[k], lw=1, ls=":")
        age = meta["substep"] - a["birth"]
        axes[2].hist((d["r"] - age)[ex], bins=np.linspace(0.999, 1.001, 41), histtype="step", color=C[k], lw=1.5)
    axes[0].set_xlabel("log10 |x - (x_src + r·centre(ipix))|  (cells)")
    axes[1].set_xlabel("log10 |dx - centre(ipix)|  (pixel sizes)\nsolid: never merged, dotted: merged")
    axes[2].set_xlabel("r - (substep - birth substep)  (cells)")
    for ax in axes:
        ax.set_ylabel("rays")
        ax.set_yscale("log")
    h, l = axes[0].get_legend_handles_labels()
    axes[2].legend(h, l, fontsize=7, loc="upper left")
    fig.tight_layout()
    fig.savefig(Path(outdir) / "raycomp_determinism.png")
    plt.close(fig)


def cmd_figures(args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    style(plt)
    Path(args.outdir).mkdir(parents=True, exist_ok=True)
    if args.codec:
        codec = json.load(open(args.codec))
        fig_bytes(codec, args.outdir, plt)
        fig_bits_by_group(codec, args.outdir, plt)
    if args.precision:
        prec = json.load(open(args.precision))
        fig_precision(prec, args.outdir, plt)
        fig_history(prec, args.outdir, plt)
    if args.det:
        fig_determinism(None, args.outdir, plt, args.det)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    e = sub.add_parser("entropy")
    e.add_argument("runs", nargs="+")
    e.add_argument("--out", required=True)
    c = sub.add_parser("codec")
    c.add_argument("runs", nargs="+")
    c.add_argument("--out", required=True)
    c.add_argument("--codec", default=str(Path(__file__).resolve().parents[2] / "build" / "ray_codec"))
    c.add_argument("--reps", type=int, default=3)
    f = sub.add_parser("figures")
    f.add_argument("--codec")
    f.add_argument("--precision")
    f.add_argument("--det", nargs="*")
    f.add_argument("--outdir", default="docs/figures")
    args = ap.parse_args()
    {"entropy": cmd_entropy, "codec": cmd_codec, "figures": cmd_figures}[args.cmd](args)


if __name__ == "__main__":
    main()
