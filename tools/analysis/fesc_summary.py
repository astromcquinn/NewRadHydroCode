#!/usr/bin/env python3
"""Runs with an escape fraction that falls with time (docs/fesc_evolution.md):
the source model, the reionization history and the IGM at z = 5-6 against
the data, in one table and one figure.

    fesc_summary.py RUN1[:label] RUN2[:label] ... --igm-dir DIR
                    [--zref 5.9,5.4,5.0] [--out fig.png] [--table t.md]
                    [--data data/igm/igm_z5_6.txt]

RUN is a problem directory (params.txt and run/). DIR holds the
igm_stats.py output of each run as igm_<basename of RUN>.json (the naming
of tools/hyak/igm_diag.slurm); a missing file leaves the IGM columns empty.
Reads:
  params.txt          f_esc, f_esc_beta, f_esc_z0, f_esc_max:
                      f_esc(z) = min(f_esc_max, f_esc ((1+z)/(1+z0))^beta)
                      (src/halo_sources.cpp f_esc_at; beta = 0: constant);
  run/movie/frames.txt  x_HII,V and x_HII,M: z at 10 ... 99.9 per cent and
                      tau_es (reion_history.py);
  the igm json        Gamma_HI(x_HI < 1e-2), lambda_mfp (exponential fit from
                      ionized starts, lam_fit_ion of igm_plot.py), <F_Lya>, A
                      and ndot_ion at every hydro step.
The IGM values are interpolated to the --zref redshifts between the two
bracketing snapshots (linear in z of the log; A linear), so runs whose
snapshots fall at slightly different z compare at the same z; nothing is
extrapolated. Ratios to the data use the same interpolation of
data/igm/igm_z5_6.txt as igm_plot.py.

Figure: f_esc(z); ndot_ion with Becker+21 and the Gaikwad-Zhu values;
1 - x_HII,V; Gamma_HI with Gaikwad+23; lambda_mfp with Becker+21 and Zhu+23;
<F_Lya> with Bosman+22.
"""
import argparse
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from igm_plot import COLORS, INK, MUTED, data_points, thr_index  # noqa: E402
from igm_stats import DATA_DEFAULT, obs_interp, read_obs  # noqa: E402
from reion_history import LEVELS, read as read_frames, tau_es, z_at  # noqa: E402

PLANCK_TAU, PLANCK_ERR = 0.054, 0.007


def read_params(path):
    p = {}
    for line in open(path):
        s = line.split("#", 1)[0].strip()
        if "=" in s:
            k, v = s.split("=", 1)
            p[k.strip()] = v.strip()
    return p


def fesc_of(p):
    f6 = float(p.get("f_esc", 0.2))
    beta = float(p.get("f_esc_beta", 0))
    z0 = float(p.get("f_esc_z0", 6))
    fmax = float(p.get("f_esc_max", 1))
    if beta == 0:
        return f6, beta, lambda z: np.full_like(np.asarray(z, float), f6)
    return f6, beta, lambda z: np.minimum(fmax, f6 * ((1 + np.asarray(z, float)) / (1 + z0)) ** beta)


def igm_series(d):
    """Snapshots sorted by z: z, x_V, Gamma_ion, lambda_fit_ion, <F>, A."""
    sn = sorted(d["snapshots"], key=lambda s: s["z"])
    out = {"z": [], "xV": [], "G": [], "lam": [], "F": [], "A": []}
    for s in sn:
        m = s["mfp"].get("ionized", s["mfp"]["random"])
        out["z"].append(s["z"])
        out["xV"].append(s["xHII_V"])
        out["G"].append(s["gamma"]["ion_V"][thr_index(s)])
        out["lam"].append(m.get("lambda_fit", float("nan")))
        out["F"].append(s["lya"]["F_mean"])
        out["A"].append(s["lya"].get("A", float("nan")))
    return {k: np.array(v, float) for k, v in out.items()}


def interp_z(zs, v, z, log=True):
    """Value at z between the bracketing snapshots; nan outside or if either is nan."""
    if len(zs) < 2 or z < zs[0] or z > zs[-1]:
        return float("nan")
    i = min(max(np.searchsorted(zs, z), 1), len(zs) - 1)
    a, b = v[i - 1], v[i]
    if not (np.isfinite(a) and np.isfinite(b)) or (log and (a <= 0 or b <= 0)):
        return float("nan")
    w = (z - zs[i - 1]) / (zs[i] - zs[i - 1])
    return float(np.exp((1 - w) * np.log(a) + w * np.log(b))) if log else float((1 - w) * a + w * b)


def ndot_at(d, z):
    nz, nd = np.array(d["ndot"]["z"]), np.array(d["ndot"]["ndot"])
    if len(nz) == 0:
        return float("nan")
    o = np.argsort(nz)
    return float(np.interp(z, nz[o], nd[o]))


def fmt(v, f="{:.2f}"):
    return "—" if v is None or not np.isfinite(v) else f.format(v)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--igm-dir", required=True)
    ap.add_argument("--zref", default="5.9,5.4,5.0")
    ap.add_argument("--data", default=DATA_DEFAULT)
    ap.add_argument("--out")
    ap.add_argument("--table")
    a = ap.parse_args()
    zref = [float(x) for x in a.zref.split(",")]
    obs = read_obs(a.data)

    runs = []
    for spec in a.runs:
        path, _, label = spec.partition(":")
        path = os.path.normpath(path)
        name = os.path.basename(path)
        p = read_params(os.path.join(path, "params.txt"))
        f6, beta, ff = fesc_of(p)
        z, xv, xm = read_frames(os.path.join(path, "run"))
        js = os.path.join(a.igm_dir, f"igm_{name}.json")
        d = json.load(open(js)) if os.path.exists(js) else None
        runs.append(dict(name=name, label=label or name, f6=f6, beta=beta, f=ff,
                         z=z, xv=xv, xm=xm, igm=d, ser=igm_series(d) if d else None))

    out = []
    out.append("**Sources.** f_esc(z) and the escaped ṅ_ion [log s⁻¹ cMpc⁻³] at the "
               "hydro step nearest z (Becker+21: 51.48 at z = 5.9, 51.20 at 5.4).")
    out.append("")
    zf = [10, 8, 7, 6, 5.4, 5]
    out.append("| run | f6 | β | " + " | ".join(f"f(z={x:g})" for x in zf) + " | "
               + " | ".join(f"ṅ(z={x:g})" for x in zref) + " |")
    out.append("|---" * (3 + len(zf) + len(zref)) + "|")
    for r in runs:
        nd = [math.log10(ndot_at(r["igm"], x)) if r["igm"] else float("nan") for x in zref]
        out.append(f"| {r['label']} | {r['f6']:g} | {r['beta']:g} | "
                   + " | ".join(f"{float(r['f'](x)):.3f}" for x in zf) + " | "
                   + " | ".join(fmt(v) for v in nd) + " |")
    out.append("")
    out.append(f"**Reionization history** (volume-weighted x_HII; τ_es with x = 1 below the "
               f"last frame; Planck 2018: {PLANCK_TAU} ± {PLANCK_ERR}).")
    out.append("")
    out.append("| run | " + " | ".join(f"z{l * 100:g}" for l in LEVELS) + " | x_V(end) | τ_es | (τ−τ_P)/σ |")
    out.append("|---" * (len(LEVELS) + 4) + "|")
    for r in runs:
        t = tau_es(r["z"], r["xm"])
        out.append(f"| {r['label']} | " + " | ".join(fmt(z_at(r["z"], r["xv"], l)) for l in LEVELS)
                   + f" | {r['xv'][-1]:.4f} | {t:.4f} | {(t - PLANCK_TAU) / PLANCK_ERR:+.1f} |")
    out.append("")
    out.append("**IGM** at fixed z (interpolated between snapshots): Γ_HI(x_HI < 10⁻²) "
               "[10⁻¹² s⁻¹], λ_mfp [pMpc], ⟨F⟩ and A, with the ratio to the data in "
               "brackets (Gaikwad+23, Becker+21/Zhu+23 pooled, Bosman+22).")
    for z in zref:
        go = obs_interp(obs, "gamma", z)
        lo = obs_interp(obs, "lambda_mfp", z)
        fo = obs_interp(obs, "mean_flux", z)
        out.append("")
        out.append(f"z = {z:g} (data: Γ {fmt(go * 1, '{:.3f}')}, λ {fmt(lo)}, ⟨F⟩ {fmt(fo, '{:.4f}')})")
        out.append("")
        out.append("| run | x_V | Γ_HI | λ_mfp | ⟨F⟩ | A |")
        out.append("|---|---|---|---|---|---|")
        for r in runs:
            s = r["ser"]
            if s is None:
                out.append(f"| {r['label']} | — | — | — | — | — |")
                continue
            xv = float(np.interp(z, r["z"][::-1], r["xv"][::-1]))
            G = interp_z(s["z"], s["G"], z) / 1e-12
            L = interp_z(s["z"], s["lam"], z)
            F = interp_z(s["z"], s["F"], z)
            A = interp_z(s["z"], s["A"], z, log=False)
            out.append(f"| {r['label']} | {xv:.4f} | {fmt(G, '{:.3f}')} ({fmt(G / go, '{:.2f}×')}) | "
                       f"{fmt(L, '{:.1f}')} ({fmt(L / lo, '{:.1f}×')}) | "
                       f"{fmt(F, '{:.3f}')} ({fmt(F / fo, '{:.1f}×')}) | {fmt(A)} |")
    text = "\n".join(out)
    print(text)
    if a.table:
        open(a.table, "w").write(text + "\n")
    if a.out:
        plot(runs, obs, a.out)


def plot(runs, obs, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2, 3, figsize=(13, 7.4), sharex=True, constrained_layout=True)
    (af, an, ax), (ag, al, aF) = axs
    data_points(an, obs["ndot"], log10=True)
    data_points(an, obs["ndot_gz"], log10=True)
    data_points(ag, obs["gamma"])
    data_points(al, obs["lambda_mfp"])
    data_points(aF, obs["mean_flux"])
    zz = np.linspace(4.9, 12, 300)
    for i, r in enumerate(runs):
        c = COLORS[i % len(COLORS)]
        lab = f"{r['label']}"
        af.plot(zz, r["f"](zz), color=c, lw=2, label=lab)
        ax.semilogy(r["z"], np.clip(1 - r["xv"], 1e-5, 1), color=c, lw=2)
        if r["igm"]:
            nz, nd = np.array(r["igm"]["ndot"]["z"]), np.array(r["igm"]["ndot"]["ndot"])
            an.plot(nz, nd, color=c, lw=2)
            s = r["ser"]
            ag.plot(s["z"], s["G"] / 1e-12, "-o", color=c, lw=2, ms=4)
            al.plot(s["z"], s["lam"], "-o", color=c, lw=2, ms=4)
            aF.plot(s["z"], s["F"], "-o", color=c, lw=2, ms=4)
    if runs and runs[0]["ser"] is not None:
        sn = sorted(runs[0]["igm"]["snapshots"], key=lambda s: s["z"])
        al.plot([s["z"] for s in sn], [s["box_pmpc"] for s in sn], ":", color=MUTED, lw=1.5)
        al.annotate("box (proper)", (sn[0]["z"], sn[0]["box_pmpc"]), xytext=(4, 4),
                    textcoords="offset points", color=MUTED, fontsize=8)
    af.set_ylabel(r"$f_{\rm esc}(z)$")
    af.set_yscale("log")
    an.set_ylabel(r"$\dot n_{\rm ion}$ [s$^{-1}$ cMpc$^{-3}$]")
    an.set_yscale("log")
    ax.set_ylabel(r"$1 - x_{\rm HII,V}$")
    ag.set_ylabel(r"$\Gamma_{\rm HI}$ [$10^{-12}$ s$^{-1}$]")
    ag.set_yscale("log")
    al.set_ylabel(r"$\lambda_{\rm mfp}^{912}$ [pMpc]")
    al.set_yscale("log")
    aF.set_ylabel(r"$\langle F_{\rm Ly\alpha}\rangle$")
    aF.set_yscale("log")
    aF.set_ylim(bottom=1e-3)
    for a_ in axs.ravel():
        a_.set_xlim(4.8, 9)
        a_.invert_xaxis()
        a_.grid(alpha=0.25, lw=0.5, which="both")
        a_.tick_params(colors=MUTED)
        a_.yaxis.label.set_color(INK)
        for s_ in ("top", "right"):
            a_.spines[s_].set_visible(False)
    for a_ in (an, ag, al, aF):
        a_.legend(frameon=False, fontsize=7.5, loc="best")
    for a_ in axs[1]:
        a_.set_xlabel("redshift")
    af.legend(frameon=False, fontsize=8)
    fig.savefig(path, dpi=140)
    print("#", path)


if __name__ == "__main__":
    main()
