#!/usr/bin/env python3
"""Compare runs measured by igm_stats.py with the z = 5-6 IGM data.

    igm_plot.py igm_RUN1.json[:label] igm_RUN2.json[:label] ...
                [--data data/igm/igm_z5_6.txt] [--out igm.png] [--table igm.txt]
                [--zlim 4.7,7.1]

Four panels against redshift, with the data of data/igm/igm_z5_6.txt and
their asymmetric 68% errors:
  Gamma_HI   solid: volume average over cells with x_HI < 1e-2 (the
             estimator compared with Gaikwad et al. 2023); dotted:
             transmission-weighted (Lya skewers); crosses: the Gamma implied
             by the Bosman mean flux, Gamma(x_HI < 1e-2) / A, where
             <exp(-A tau)> = <F>_obs (photoionization equilibrium, tau ~ 1/Gamma).
  lambda_mfp solid: lambda_e (<exp(-tau_912)> = 1/e) from starts in ionized
             cells; dashed: from source cells; the grey dotted line is the
             proper box size of the first run.
  <F_Lya>    mean transmitted flux over skewers along the three axes.
  ndot       escaped ionizing photons s^-1 cMpc^-3 at every hydro step
             (step_NNNN/stats.txt), with the Becker et al. 2021 values and
             the ndot derived from the Gaikwad Gamma and the Zhu lambda.
The text table (default: --out with .txt) lists, per run and snapshot, the
measured values, the data interpolated to the same z and their ratios.
"""
import argparse
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from igm_stats import DATA_DEFAULT, obs_interp, read_obs  # noqa: E402

# reference categorical palette (dataviz skill), fixed order; data in ink
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
INK, MUTED = "#0b0b0b", "#8a8984"
REF_STYLE = {"Becker21": ("s", "Becker+21"), "Zhu23": ("o", "Zhu+23"),
             "Gaikwad23": ("o", "Gaikwad+23"), "Bosman22": ("o", "Bosman+22"),
             "Gaikwad23_Gamma+Zhu23_lambda": ("D", r"$\Gamma$ Gaikwad+23, $\lambda$ Zhu+23")}


def load(spec):
    path, label = (spec.rsplit(":", 1) + [None])[:2] if ":" in spec else (spec, None)
    d = json.load(open(path))
    d["label"] = label or d["name"]
    return d


def thr_index(snap, t=1e-2):
    return snap["gamma"]["thresholds"].index(t)


def data_points(ax, b, scale=1.0, log10=False):
    for ref in dict.fromkeys(b["ref"]):
        m = np.array([r == ref for r in b["ref"]])
        z, v, ep, em = b["z"][m], b["v"][m], b["ep"][m], b["em"][m]
        if log10:
            y = 10 ** v
            lo, hi = y - 10 ** (v - em), 10 ** (v + ep) - y
        else:
            y, lo, hi = v * scale, em * scale, ep * scale
        mk, lab = REF_STYLE.get(ref, ("o", ref))
        ax.errorbar(z, y, yerr=[lo, hi], fmt=mk, ms=4.5, color=INK, mfc="white" if mk == "D"
                    else INK, ecolor=INK, elinewidth=1, capsize=2, label=lab, zorder=5)


def ndot_at(run, z):
    nz, nd = np.array(run["ndot"]["z"]), np.array(run["ndot"]["ndot"])
    if len(nz) == 0:
        return float("nan")
    return float(nd[np.argmin(np.abs(nz - z))])


def table(runs, obs, path):
    cols = ["z", "xHII_V", "G_V", "G_ion", "G_Fw", "G_reqF", "G_obs", "G/obs",
            "lam_e_rnd", "lam_e_ion", "lam_fit_ion", "lam_e_src", "lam_obs", "lam/obs", "box",
            "F", "F_obs", "F/obs", "A", "lg_ndot", "lg_ndot_B21", "lg_ndot_GZ"]
    lines = ["# igm_plot.py: Gamma in 1e-12 s^-1 (G_ion: x_HI < 1e-2; G_Fw: transmission-"
             "weighted; G_reqF = G_ion/A), lambda and box in pMpc (lam_e_*: <F> = 1/e from "
             "random / ionized / source starts; lam_fit_ion: exp fit), A = tau rescaling for "
             "Bosman <F>, ndot in log10 s^-1 cMpc^-3 at the nearest hydro step; *_obs: data "
             "interpolated to z (lam_obs: Becker+21 and Zhu+23 pooled; G_obs: Gaikwad+23; lg_ndot_GZ: ndot_gz block)"]
    for r in runs:
        lines.append(f"## {r['label']}  ({r['run']}; N = {r['N']}, box = {r['box']} Mpc/h, "
                     f"f_esc = {r['f_esc']}, c_reduced = {r['c_reduced']})")
        lines.append(" ".join(f"{c:>11s}" for c in cols))
        for s in r["snapshots"]:
            z = s["z"]
            g, m, ly = s["gamma"], s["mfp"], s["lya"]
            i2 = thr_index(s)
            gi = g["ion_V"][i2] / 1e-12
            go = obs_interp(obs, "gamma", z)
            lo = obs_interp(obs, "lambda_mfp", z)
            li = m.get("ionized", m["random"])["lambda_e"]
            nd = ndot_at(r, z)
            vals = [z, s["xHII_V"], g["V"] / 1e-12, gi, ly["gamma_Fweighted"] / 1e-12,
                    g.get("req_from_F", float("nan")) / 1e-12, go, gi / go,
                    m["random"]["lambda_e"], li, m.get("ionized", m["random"])["lambda_fit"],
                    m.get("sources", {}).get("lambda_e", float("nan")), lo, li / lo,
                    s["box_pmpc"], ly["F_mean"], ly.get("F_obs", float("nan")),
                    ly["F_mean"] / ly.get("F_obs", float("nan")), ly.get("A", float("nan")),
                    math.log10(nd) if nd > 0 else float("nan"),
                    obs_interp(obs, "ndot", z), obs_interp(obs, "ndot_gz", z)]
            lines.append(" ".join(f"{v:11.4g}" for v in vals))
    open(path, "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print("#", path)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--data", default=DATA_DEFAULT)
    ap.add_argument("--out", default="igm.png")
    ap.add_argument("--table", default=None)
    ap.add_argument("--zlim", default="4.7,7.1")
    a = ap.parse_args()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    obs = read_obs(a.data)
    runs = [load(s) for s in a.runs]
    zlo, zhi = [float(x) for x in a.zlim.split(",")]
    fig, axs = plt.subplots(2, 2, figsize=(10, 7.6), sharex=True, constrained_layout=True)
    (ag, al), (af, an) = axs
    data_points(ag, obs["gamma"])
    data_points(al, obs["lambda_mfp"])
    data_points(af, obs["mean_flux"])
    data_points(an, obs["ndot"], log10=True)
    data_points(an, obs["ndot_gz"], log10=True)
    handles = []
    for i, r in enumerate(runs):
        c = COLORS[i % len(COLORS)]
        sn = sorted(r["snapshots"], key=lambda s: s["z"])
        z = np.array([s["z"] for s in sn])
        if len(sn):
            gi = np.array([s["gamma"]["ion_V"][thr_index(s)] for s in sn]) / 1e-12
            gf = np.array([s["lya"]["gamma_Fweighted"] for s in sn]) / 1e-12
            gr = np.array([s["gamma"].get("req_from_F", np.nan) for s in sn]) / 1e-12
            ag.plot(z, gi, "-o", color=c, lw=2, ms=4)
            ag.plot(z, gf, ":", color=c, lw=1.5)
            ag.plot(z, gr, "x", color=c, ms=7, mew=1.8)
            li = np.array([s["mfp"].get("ionized", s["mfp"]["random"])["lambda_e"] for s in sn])
            le = np.array([s["mfp"].get("ionized", s["mfp"]["random"])["lambda_e_err"]
                           for s in sn])
            ls = np.array([s["mfp"].get("sources", {}).get("lambda_e", np.nan) for s in sn])
            al.errorbar(z, li, yerr=le, fmt="-o", color=c, lw=2, ms=4, capsize=0)
            al.plot(z, ls, "--", color=c, lw=1.5)
            af.plot(z, [s["lya"]["F_mean"] for s in sn], "-o", color=c, lw=2, ms=4)
        nz, nd = np.array(r["ndot"]["z"]), np.array(r["ndot"]["ndot"])
        an.plot(nz, nd, "-", color=c, lw=2)
        handles.append(Line2D([], [], color=c, lw=2, marker="o", ms=4,
                              label=f"{r['label']} (f_esc = {r['f_esc']:g})"))
    if runs and runs[0]["snapshots"]:
        sn = sorted(runs[0]["snapshots"], key=lambda s: s["z"])
        al.plot([s["z"] for s in sn], [s["box_pmpc"] for s in sn], ":", color=MUTED, lw=1.5)
        al.annotate("box (proper)", (sn[0]["z"], sn[0]["box_pmpc"]), xytext=(4, 4),
                    textcoords="offset points", color=MUTED, fontsize=8)
    ag.set_ylabel(r"$\Gamma_{\rm HI}$ [$10^{-12}$ s$^{-1}$]")
    ag.set_yscale("log")
    al.set_ylabel(r"$\lambda_{\rm mfp}^{912}$ [pMpc]")
    al.set_yscale("log")
    af.set_ylabel(r"$\langle F_{\rm Ly\alpha}\rangle$")
    af.set_yscale("log")
    an.set_ylabel(r"$\dot n_{\rm ion}$ [s$^{-1}$ cMpc$^{-3}$]")
    an.set_yscale("log")
    for ax in axs.ravel():
        ax.set_xlim(zlo, zhi)
        ax.grid(alpha=0.25, lw=0.5, which="both")
        ax.legend(frameon=False, fontsize=7.5, loc="best")
    for ax in axs[1]:
        ax.set_xlabel("redshift")
    style = [Line2D([], [], color=MUTED, lw=2, marker="o", ms=4, label=r"$x_{\rm HI}<10^{-2}$ cells"),
             Line2D([], [], color=MUTED, lw=1.5, ls=":", label="transmission-weighted"),
             Line2D([], [], color=MUTED, ls="", marker="x", mew=1.8, label=r"implied by $\langle F\rangle_{\rm obs}$")]
    ag.legend(handles=ag.get_legend_handles_labels()[0] + style, frameon=False, fontsize=7.5)
    al.legend(handles=al.get_legend_handles_labels()[0] + [
        Line2D([], [], color=MUTED, lw=2, marker="o", ms=4, label="from ionized cells"),
        Line2D([], [], color=MUTED, lw=1.5, ls="--", label="from sources")],
        frameon=False, fontsize=7.5)
    fig.legend(handles=handles, loc="outside upper center", ncol=min(4, len(handles)),
               frameon=False, fontsize=9)
    fig.savefig(a.out, dpi=140)
    print("#", a.out)
    table(runs, obs, a.table or os.path.splitext(a.out)[0] + ".txt")


if __name__ == "__main__":
    main()
