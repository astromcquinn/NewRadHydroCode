#!/usr/bin/env python3
"""Figures and tables of the convergence series (docs/rt_moments.md section 8)
from the output of rt_convergence_cmp.py and igm_stats.py.

usage: rt_convergence_fig.py CMP.json [CMP2.json ...] --igm DIR[,DIR2] --outdir docs/figures

Writes rt_convergence_res.{png,txt} (M1 and M1 with moment_inject = 0 against
the rays at N = 128, 256, 512, at matched x_HII,V) and
rt_convergence_rays.{png,txt} (the ray budgets against the largest).
"""
import argparse
import glob
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# the snapshots are few (zdata1) and the histories differ with N: the three
# epochs before overlap are those every N has a pair of snapshots near
TARGETS = [(0.57, r"$x_{\rm HII,V}\approx0.57$"), (0.77, r"$\approx0.77$"), (0.89, r"$\approx0.86$-$0.92$"),
           (2.0, "after overlap")]
KTAB = (0.3, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0)
COL = {128: "C0", 256: "C1", 512: "C3"}


def load(paths):
    groups, box = {}, None
    for p in paths:
        j = json.load(open(p))
        box = j["box"]
        groups.update(j["groups"])
    return groups, box


def load_igm(dirs):
    out = {}
    for d in dirs:
        for f in glob.glob(os.path.join(d, "igm_*.json")):
            j = json.load(open(f))
            out[os.path.basename(f)[4:-5]] = [(s["z"], s["mfp"]["random"]["lambda_e"]) for s in j["snapshots"]
                                              if "mfp" in s]
    return out


def mfp(igm, sn):
    name = sn["dir"].split(os.sep)[0]
    if name not in igm or not igm[name]:
        return float("nan")
    z, lam = min(igm[name], key=lambda t: abs(t[0] - sn["z"]))
    return lam if abs(z - sn["z"]) < 0.02 else float("nan")


def pick_ref(g, target, lab=None):
    """Index of the reference snapshot for a target x_V (2.0: the last one, if
    past overlap); with lab, the one whose paired snapshot of run lab is also
    near the target."""
    ref = g["runs"][g["reference"]]
    ok = [i for i, s in enumerate(ref) if "xV" in s]
    if target > 1:
        i = min(ok, key=lambda q: ref[q]["z"])
        return i if ref[i]["xV"] > 0.995 else None

    def cost(q):
        c = abs(ref[q]["xV"] - target)
        if lab in g["runs"]:
            pr = pair_of(g, lab, q)
            c += abs(pr[0]["xV"] - ref[q]["xV"]) if pr else 1.0
        return c
    i = min(ok, key=cost)
    return i if abs(ref[i]["xV"] - target) < 0.08 else None


def pair_of(g, lab, iref):
    """The pair (run snapshot, cross spectra) of run lab with reference snapshot iref."""
    ref = g["runs"][g["reference"]][iref]
    c = [p for p in g["pairs"] if p["run"] == lab and p["ref_i"] == iref]
    if not c:
        return None
    if ref["xV"] > 0.995:
        p = min(c, key=lambda p: abs(g["runs"][lab][p["i"]]["z"] - ref["z"]))
    else:
        p = min(c, key=lambda p: abs(g["runs"][lab][p["i"]]["xV"] - ref["xV"]))
    return g["runs"][lab][p["i"]], p


def spectra(sn, ref, p, f):
    """k, P/P_ref, r(k) of field f ('gamma' or 'x') where both have modes."""
    k = np.array(ref["k"])
    ok = np.array(ref["nmodes"]) > 0
    a, b, c = (np.array(sn["P_" + f]), np.array(ref["P_" + f]), np.array(p["C_" + f]))
    with np.errstate(all="ignore"):
        return k[ok], (a / b)[ok], (c / np.sqrt(a * b))[ok]


def at(k, y, kt):
    m = np.isfinite(y) & (y > -9)
    if kt < k[m][0] * 0.9 or kt > k[m][-1] * 1.1:
        return float("nan")
    return float(np.interp(np.log(kt), np.log(k[m]), y[m]))


def fmt_row(vals):
    return " ".join("%6.2f" % v if np.isfinite(v) else "     -" for v in vals)


def res_series(groups, box, igm, outdir):
    Ns = sorted(int(g[1:]) for g in groups if g.startswith("N"))
    fig, ax = plt.subplots(4, len(TARGETS), figsize=(4.1 * len(TARGETS), 12.4), sharex=True, squeeze=False)
    rep = ["# Resolution series (tools/hyak/rt_convergence.sh; 20 Mpc/h, 512^3 particles, gas and RT on N^3).",
           "# P_Gamma of delta = Gamma_HI/<Gamma_HI>_V - 1 (gammaHI_avg), ratio to the ray run of the same N and",
           "# cross-correlation r(k) with it, at snapshots matched in x_HII,V (after overlap: in z).",
           "# lambda: the mean free path at 912 A, lambda_e of igm_stats.py (random starts) [pMpc].",
           "# k [h/Mpc] columns: " + " ".join("%6.1f" % k for k in KTAB)]
    raysP = {}
    for it, (tx, tlab) in enumerate(TARGETS):
        rep.append("\n## " + ("after overlap (last snapshot)" if tx > 1 else "x_HII,V ~ %.2f" % tx))
        for N in Ns:
            g = groups["N%d" % N]
            ref = g["runs"][g["reference"]]
            ir = pick_ref(g, tx, "M1")
            if ir is None:
                rep.append("N = %d: no %s snapshot at this epoch" % (N, g["reference"]))
                continue
            r0 = ref[ir]
            raysP[(it, N)] = r0
            rep.append("N = %d  reference %s: z = %.2f, x_V = %.3f, <Gamma>_V = %.3g, median Gamma(ionised) = %.3g, "
                       "lambda = %.2f" % (N, g["reference"], r0["z"], r0["xV"], r0["gamma_V"],
                                          r0["gamma_med_ion"], mfp(igm, r0)))
            dx = box / N
            for lab, ls in (("M1", "-"), ("M1_inject0", "--"), ("OTVET", "-."), ("M1_order2", (0, (4, 1, 1, 1, 1, 1))),
                            ("rays", ":")):
                if lab not in g["runs"] or lab == g["reference"]:
                    continue
                pr = pair_of(g, lab, ir)
                if pr is None:
                    rep.append("   %-10s no snapshot paired with this one" % lab)
                    continue
                sn, p = pr
                k, ratio, r = spectra(sn, r0, p, "gamma")
                kx, ratx, rx = spectra(sn, r0, p, "x")
                rep.append("   %-10s z = %.2f, x_V = %.3f, <Gamma>_V = %.3g, median Gamma(ionised) = %.3g, lambda = %.2f"
                           % (lab, sn["z"], sn["xV"], sn["gamma_V"], sn["gamma_med_ion"], mfp(igm, sn)))
                rep.append("      P_Gamma/ref " + fmt_row([at(k, ratio, kt) for kt in KTAB]))
                rep.append("      r_Gamma     " + fmt_row([at(k, r, kt) for kt in KTAB]))
                rep.append("      P_xHII/ref  " + fmt_row([at(kx, ratx, kt) for kt in KTAB]))
                rep.append("      r_xHII      " + fmt_row([at(kx, rx, kt) for kt in KTAB]))
                x1 = k[np.argmax((ratio > 1) & (k > 2))] if ((ratio > 1) & (k > 2)).any() else float("nan")
                if lab == "M1":
                    # the first k above 2 h/Mpc where M1 has more power than the rays, log-interpolated
                    j = np.where((ratio[1:] > 1) & (ratio[:-1] <= 1) & (k[1:] > 2))[0]
                    if len(j):
                        j = j[0]
                        x1 = np.exp(np.interp(0.0, np.log(ratio[j:j + 2]), np.log(k[j:j + 2])))
                        rep.append("      M1 overtakes the reference at k = %.1f h/Mpc, k dx = %.2f" % (x1, x1 * dx))
                    else:
                        rep.append("      M1 does not cross the reference from below above k = 2 h/Mpc")
                if lab == "rays":
                    lw, al = 1.0, 0.8
                else:
                    lw, al = 1.6, 1.0
                kw = dict(color=COL.get(N, "k"), ls=ls, lw=lw, alpha=al)
                ax[0, it].loglog(k, np.maximum(ratio, 1e-4), label="N = %d %s (x = %.2f | %.2f)"
                                 % (N, lab.replace("_inject0", " inj 0").replace("_order2", " order 2"), sn["xV"], r0["xV"]), **kw)
                ax[1, it].semilogx(k, r, **kw)
            ax[0, it].axvline(1.3 / dx, color=COL.get(N, "k"), lw=0.6, alpha=0.5)
        # the references themselves against N at fixed k
        base = raysP.get((it, 256))
        for N in Ns:
            s = raysP.get((it, N))
            if s is None:
                continue
            k = np.array(s["k"])
            ok = np.array(s["nmodes"]) > 0
            d2 = k ** 3 * np.array(s["P_gamma"]) / (2 * np.pi ** 2)
            g = groups["N%d" % N]
            ax[2, it].loglog(k[ok], d2[ok], color=COL.get(N, "k"),
                             label="N = %d %s (z = %.2f, x = %.2f)" % (N, g["reference"], s["z"], s["xV"]))
            if base is not None and N != 256:
                kb = np.array(base["k"])
                okb = (np.array(base["nmodes"]) > 0) & ok
                with np.errstate(all="ignore"):
                    rr = (np.array(s["P_gamma"]) / np.array(base["P_gamma"]))[okb]
                rep.append("N = %d reference / N = 256 reference, P_Gamma: " % N + fmt_row([at(kb[okb], rr, kt) for kt in KTAB]))
        # M1 itself against N (its own snapshot nearest the epoch: the 512^3
        # ray run may not have reached it)
        for N in Ns:
            sns = [s for s in groups["N%d" % N]["runs"].get("M1", []) if "xV" in s]
            if not sns:
                continue
            if tx > 1:
                s = min(sns, key=lambda q: q["z"])
                s = s if s["xV"] > 0.995 else None
            else:
                s = min(sns, key=lambda q: abs(q["xV"] - tx))
                s = s if abs(s["xV"] - tx) < 0.08 else None
            if s is None:
                continue
            k = np.array(s["k"])
            ok = np.array(s["nmodes"]) > 0
            d2 = k ** 3 * np.array(s["P_gamma"]) / (2 * np.pi ** 2)
            ax[3, it].loglog(k[ok], d2[ok], color=COL.get(N, "k"),
                             label="N = %d M1 (z = %.2f, x = %.2f)" % (N, s["z"], s["xV"]))
            rep.append("N = %d M1 alone: z = %.2f, x_V = %.3f, lambda = %.2f, Delta2_Gamma: " % (N, s["z"], s["xV"], mfp(igm, s))
                       + " ".join("%8.3g" % at(k[ok], d2[ok], kt) for kt in KTAB))
        ax[3, it].legend(fontsize=6, loc="upper left")
        ax[3, it].set_xlabel(r"$k$ [$h$/Mpc]")
        ax[0, it].set_title(tlab)
        ax[0, it].axhline(1, color="k", lw=0.5)
        ax[0, it].set_ylim(3e-3, 30)
        ax[1, it].set_ylim(-0.3, 1.05)
        ax[1, it].axhline(1, color="k", lw=0.5)
        ax[0, it].legend(fontsize=6, loc="lower left")
        ax[2, it].legend(fontsize=6, loc="upper left")
    ax[0, 0].set_ylabel(r"$P_\Gamma$ / $P_\Gamma$(rays, same $N$)")
    ax[1, 0].set_ylabel(r"$r(k)$ with the rays")
    ax[2, 0].set_ylabel(r"$\Delta^2_\Gamma$ of the ray runs")
    ax[3, 0].set_ylabel(r"$\Delta^2_\Gamma$ of the M1 runs")
    fig.suptitle("M1 (solid), M1 inject 0 (dashed), OTVET (dash-dot), M1 moment_order 2 (dash-dot-dot), "
                 "f = 15 rays where not the reference (dotted) against the "
                 r"ray reference; thin vertical lines: $k\Delta x = 1.3$", fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "rt_convergence_res.png"), dpi=130)
    # global histories
    rep.append("\n## histories: z, x_HII,V, <Gamma>_V, median Gamma(ionised) [1/s], lambda [pMpc]; "
               "cells above 100 x mean density: number, x_HI (mass-weighted), mean Gamma")
    for N in Ns:
        g = groups["N%d" % N]
        for lab, sns in g["runs"].items():
            rep.append("N = %d %s" % (N, lab))
            for s in sns:
                if "xV" not in s:
                    continue
                d = s["dense"]
                rep.append("   %.2f %.4f %.3e %.3e %6.2f | %7d %s"
                           % (s["z"], s["xV"], s["gamma_V"], s["gamma_med_ion"] or float("nan"), mfp(igm, s), d["n"],
                              ("%.3e %.3e" % (d["xHI_M"], d["gamma_mean"])) if d["n"] else "-"))
    open(os.path.join(outdir, "rt_convergence_res.txt"), "w").write("\n".join(rep) + "\n")


def ray_series(groups, igm, outdir):
    if "rays" not in groups:
        return
    g = groups["rays"]
    ref = g["runs"][g["reference"]]
    labs = [l for l in g["runs"] if l != g["reference"]]
    cols = {l: c for l, c in zip(sorted(labs, key=lambda s: int(s[1:])), ("C3", "C1", "C2", "C0", "C4"))}
    fig, ax = plt.subplots(4, len(TARGETS), figsize=(4.1 * len(TARGETS), 11), sharex=True, squeeze=False)
    rep = ["# Ray-budget series at 256^3 (nray_max_factor f; reference %s): ratio of P(k) to the reference and" % g["reference"],
           "# r(k) with it, for delta_Gamma (gammaHI_avg) and x_HII, at the reference's snapshots.",
           "# k [h/Mpc] columns: " + " ".join("%6.1f" % k for k in KTAB)]
    for it, (tx, tlab) in enumerate(TARGETS):
        ir = pick_ref(g, tx)
        if ir is None:
            continue
        r0 = ref[ir]
        rep.append("\n## %s: reference z = %.2f, x_V = %.4f, median Gamma(ionised) = %.3e, lambda = %.2f pMpc"
                   % ("after overlap" if tx > 1 else "x_HII,V ~ %.2f" % tx, r0["z"], r0["xV"], r0["gamma_med_ion"],
                      mfp(igm, r0)))
        for lab in sorted(labs, key=lambda s: int(s[1:])):
            pr = pair_of(g, lab, ir)
            if pr is None:
                continue
            sn, p = pr
            k, ratio, r = spectra(sn, r0, p, "gamma")
            kx, ratx, rx = spectra(sn, r0, p, "x")
            rep.append("  %-4s z = %.2f, x_V = %.4f, median Gamma(ionised) = %.3e, lambda = %.2f"
                       % (lab, sn["z"], sn["xV"], sn["gamma_med_ion"], mfp(igm, sn)))
            rep.append("      P_Gamma/ref " + fmt_row([at(k, ratio, kt) for kt in KTAB]))
            rep.append("      r_Gamma     " + fmt_row([at(k, r, kt) for kt in KTAB]))
            rep.append("      P_xHII/ref  " + fmt_row([at(kx, ratx, kt) for kt in KTAB]))
            rep.append("      r_xHII      " + fmt_row([at(kx, rx, kt) for kt in KTAB]))
            dev = np.abs(ratio - 1)
            bad = k[dev > 0.1]
            rep.append("      |P_Gamma/ref - 1| max at k <= 2 h/Mpc: %.3f; first k with more than 10 %%: %s"
                       % (dev[k <= 2].max(), "%.1f" % bad[0] if len(bad) else "none up to %.0f" % k[-1]))
            kw = dict(color=cols[lab], label=lab)
            ax[0, it].semilogx(k, ratio, **kw)
            ax[1, it].semilogx(k, 1 - r, **kw)
            ax[2, it].semilogx(kx, ratx, **kw)
            ax[3, it].semilogx(kx, 1 - rx, **kw)
        ax[0, it].set_title("%s (z = %.2f, x = %.3f)" % (tlab, r0["z"], r0["xV"]), fontsize=9)
        for a in (ax[0, it], ax[2, it]):
            a.axhline(1, color="k", lw=0.5)
            a.axhspan(0.9, 1.1, color="k", alpha=0.07)
            a.set_ylim(0.6, 1.6)
        for a in (ax[1, it], ax[3, it]):
            a.set_yscale("symlog", linthresh=1e-3)
            a.set_ylim(-1e-3, 1)
        ax[3, it].set_xlabel(r"$k$ [$h$/Mpc]")
    ax[0, 0].legend(fontsize=7)
    ax[0, 0].set_ylabel(r"$P_\Gamma(f)/P_\Gamma$(%s)" % g["reference"])
    ax[1, 0].set_ylabel(r"$1 - r_\Gamma(k)$")
    ax[2, 0].set_ylabel(r"$P_{x_{\rm HII}}(f)/P_{x_{\rm HII}}$(%s)" % g["reference"])
    ax[3, 0].set_ylabel(r"$1 - r_{x_{\rm HII}}(k)$")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "rt_convergence_rays.png"), dpi=130)
    rep.append("\n## histories: z, x_HII,V, <Gamma>_V, median Gamma(ionised) [1/s], lambda [pMpc]")
    for lab, sns in g["runs"].items():
        rep.append(lab)
        for s in sns:
            if "xV" in s:
                rep.append("   %.2f %.4f %.3e %.3e %6.2f" % (s["z"], s["xV"], s["gamma_V"], s["gamma_med_ion"], mfp(igm, s)))
    open(os.path.join(outdir, "rt_convergence_rays.txt"), "w").write("\n".join(rep) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("cmp", nargs="+")
    ap.add_argument("--igm", default="")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    groups, box = load(a.cmp)
    igm = load_igm([d for d in a.igm.split(",") if d])
    os.makedirs(a.outdir, exist_ok=True)
    res_series(groups, box, igm, a.outdir)
    ray_series(groups, igm, a.outdir)


if __name__ == "__main__":
    main()
