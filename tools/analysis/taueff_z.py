#!/usr/bin/env python3
"""tau_eff against z (and against lambda_mfp) for the RT methods, raw (no rescaling).

Reads what taueff_slabs.py and igm_stats.py already wrote; computes nothing new:

  SLABS/<label>_slabs.txt   50 Mpc/h segment tau_eff from the movie slabs,
                            every dz = 0.025 (real space; the raw columns)
  IGM/<json>                igm_stats.py of the snapshots: -ln<F> of all
                            pixels (redshift space) and lambda_mfp

Usage:
  taueff_z.py --slabs DIR --outdir DIR label:igm.json [label:igm.json ...]

Panels:
  (a) median of the 50 Mpc/h tau_eff against z, the 16-84 % band of the
      first run shaded, -ln<F> of the snapshots as markers;
  (b) the median over the first run's median;
  (c) the 84-16 % width;
  (d) the median against lambda_mfp (random starts, <F> = 1/e;
      docs/igm_diagnostics.md), ln lambda interpolated linearly in z between
      the snapshots. Where lambda_e is undefined (it exceeds the 40 pMpc
      rays at z ~ 5) the exponential fit from ionized starts is used.

The slab curves are smoothed with a [1, 2, 1]/4 kernel over frames (end
frames dropped): rho and T in the slabs change only at hydro steps, about
every second frame, so the raw curves step in pairs, and two runs whose hydro
steps fall on different frames zigzag against each other in (b).

Writes taueff_z.png and taueff_z.txt (docs/taueff_21cmfast.md A.6), and
taueff_mfp_width.png: lambda_mfp of every run at the snapshots (left; open
symbols where the exponential fit stands in) beside the 84-16 % width of
the raw tau_eff against z (right).

--full FILE (taueff_plot.py's taueff.txt): the full-box values at the
snapshots (3 x 8192 random skewers, redshift space, 50 Mpc/h joined
segments) as filled markers on every panel, for comparison with the slabs.

Without --slabs (only --full): the snapshots alone, connected by lines; no
slab curves, no smoothing, lambda_mfp exact at each snapshot.

--norm: the same figures (suffix _A) for tau rescaled per frame by A so that
<F> of the slab pixels equals Bosman+22's mean transmission at the frame's z
(data/igm/igm_z5_6.txt, Table 4; taueff_slabs.py; z <= 6.2 only). The
markers of (a) are then the observed -ln<F>.
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def read_slabs(path, norm=False):
    """z, tau_eff at 5, 16, 50, 84, 95 % (raw, or rescaled by A if norm)."""
    sl = np.loadtxt(path, ndmin=2)
    z, p = sl[:, 1], (sl[:, 5:10] if norm else sl[:, 11:16])
    o = np.argsort(z)
    return z[o], p[o]


def read_igm(path):
    """Per snapshot: z, -ln<F>, lambda_mfp in pMpc (and whether it is the fit)."""
    d = json.load(open(path))
    z, tau, lam, fit = [], [], [], []
    for s in d["snapshots"]:
        m = s.get("mfp", {})
        le = m.get("random", {}).get("lambda_e", np.nan)
        isfit = not np.isfinite(le)
        if isfit:
            le = m.get("ionized", {}).get("lambda_fit", np.nan)
        z.append(s["z"])
        tau.append(s["lya"]["tau_eff"])
        lam.append(le)
        fit.append(isfit)
    o = np.argsort(z)
    return np.array(z)[o], np.array(tau)[o], np.array(lam)[o], np.array(fit)[o]


def read_full(path, norm):
    """taueff_plot.py's taueff.txt: label -> (z, tau at 5..95 % of the 50 Mpc/h
    segments, raw or rescaled)."""
    out = {}
    if not path:
        return out
    for line in open(path):
        if line.startswith("#"):
            continue
        t = line.split()
        if len(t) != 27:          # the CDF block that follows has other columns
            continue
        p = [float(v) for v in (t[22:27] if norm else t[17:22])]
        out.setdefault(t[0], []).append([float(t[1])] + p)
    return {k: np.array(sorted(v)) for k, v in out.items()}


def smooth121(z, p):
    """[1, 2, 1]/4 over frames along axis 0, the end frames dropped."""
    return z[1:-1], 0.25 * p[:-2] + 0.5 * p[1:-1] + 0.25 * p[2:]


def interp_ext(x, xp, fp):
    """np.interp, extended linearly from the two end points on either side."""
    o = np.argsort(xp)
    xp, fp = np.asarray(xp)[o], np.asarray(fp)[o]
    y = np.interp(x, xp, fp)
    lo, hi = x < xp[0], x > xp[-1]
    y[lo] = fp[0] + (x[lo] - xp[0]) * (fp[1] - fp[0]) / (xp[1] - xp[0])
    y[hi] = fp[-1] + (x[hi] - xp[-1]) * (fp[-1] - fp[-2]) / (xp[-1] - xp[-2])
    return y


def lam_at(z, zs, lam):
    ok = np.isfinite(lam) & (lam > 0)
    out = np.exp(np.interp(z, zs[ok], np.log(lam[ok])))
    out[(z < zs[ok].min() - 5e-3) | (z > zs[ok].max() + 5e-3)] = np.nan   # 4-digit z of the tables
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("runs", nargs="+", help="label:igm_stats json (the label names SLABS/<label>_slabs.txt)")
    ap.add_argument("--slabs", default=None, help="the slab curves; omit to plot the snapshots alone (--full)")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--dpi", type=int, default=200)
    ap.add_argument("--norm", action="store_true", help="tau rescaled to Bosman+22's mean flux")
    ap.add_argument("--full", default=None, help="taueff_plot.py's taueff.txt: full-box snapshot markers")
    a = ap.parse_args()
    sfx = "_A" if a.norm else ""
    kind = "rescaled to $\\langle F\\rangle_{\\rm obs}$" if a.norm else "raw"
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from rt_suite_figures import COLORS, STYLES

    os.makedirs(a.outdir, exist_ok=True)
    full = read_full(a.full, a.norm)
    snap_only = a.slabs is None
    if snap_only and not full:
        ap.error("without --slabs, --full is required")
    runs = []
    for spec in a.runs:
        lab, _, js = spec.partition(":")
        if snap_only:
            fb = full[lab]
            fb = fb[np.all(np.isfinite(fb[:, 1:6]), axis=1)]
            zs, ps = fb[:, 0], fb[:, 1:6]
        else:
            zs, ps = read_slabs(os.path.join(a.slabs, f"{lab}_slabs.txt"), a.norm)
        runs.append((lab, zs, ps, read_igm(js)))
    sm = (lambda z, p: (z, p)) if snap_only else smooth121
    if snap_only:
        sfx = "_snap" + sfx
        full = {}
    fr0 = full.get(runs[0][0])

    fig, ax = plt.subplots(2, 2, figsize=(9.0, 7.0))
    (a0, a1), (a2, a3) = ax
    z0, p0 = sm(*[x[np.all(np.isfinite(runs[0][2][:, 1:4]), axis=1)] for x in runs[0][1:3]])
    with open(os.path.join(a.outdir, f"taueff_z{sfx}.txt"), "w") as f:
        f.write(f"# taueff_z.py: {'tau rescaled to the Bosman+22 mean flux' if a.norm else 'raw'}: tau_eff over 50 Mpc/h segments "
                + ("at the snapshots (full box, 3 x 8192 skewers, redshift space)\n# " if snap_only else
                   "from the movie slabs (real space)\n# smoothed [1,2,1]/4 over frames; ")
                + "columns: label z tau16 tau50 tau84 width(84-16) lambda_mfp[pMpc] (interpolated)\n")
        for ir, (lab, zs, ps, (zi, ti, li, fi)) in enumerate(runs):
            c, ls = COLORS[ir % len(COLORS)], STYLES[ir % len(STYLES)]
            lw = 1.8 if ir == 0 else 1.1
            m = np.all(np.isfinite(ps[:, 1:4]), axis=1)
            zs, ps = sm(zs[m], ps[m])
            m = np.ones(len(zs), bool)
            mko = dict(marker="os^vDph"[ir % 7], ms=4.5, mec="white", mew=0.5) if snap_only else {}
            z, med, w = zs, ps[:, 2], ps[:, 3] - ps[:, 1]
            lam = lam_at(z, zi, li)
            for row in zip(z, ps[m, 1], med, ps[m, 3], w, lam):
                f.write(f"{lab:10s} " + " ".join(f"{v:.4g}" for v in row) + "\n")
            if ir == 0:
                a0.fill_between(z, ps[m, 1], ps[m, 3], color=c, alpha=0.12, lw=0,
                                label=f"{lab} 16-84 %")
            a0.plot(z, med, color=c, ls=ls, lw=lw, label=lab, **mko)
            k = (zi >= z.min() - 0.02) & (zi <= z.max() + 0.02)
            if not a.norm and not full and not snap_only:
                a0.plot(zi[k], ti[k], color=c, marker="o", ms=3.5, mfc="none", ls="none")
            fb = full.get(lab)
            if fb is not None:
                fb = fb[np.all(np.isfinite(fb[:, 1:]), axis=1)]
                mk = dict(color=c, marker="os^vDph"[ir % 7], ms=4.5, ls="none", mec="white", mew=0.5, zorder=5)
                a0.plot(fb[:, 0], fb[:, 3], **mk)
                if fr0 is not None:
                    g = fr0[np.isfinite(fr0[:, 3])]
                    ref = np.exp(np.interp(fb[:, 0], g[:, 0], np.log(g[:, 3]), left=np.nan, right=np.nan))
                    a1.plot(fb[:, 0], fb[:, 3] / ref, **mk)
                a2.plot(fb[:, 0], fb[:, 4] - fb[:, 2], **mk)
                a3.plot(lam_at(fb[:, 0], zi, li), fb[:, 3], **mk)
            ref = np.exp(interp_ext(z, z0, np.log(p0[:, 2])))   # log-linear past the ends, dz <= 0.03
            ref[(z < z0.min() - 0.03) | (z > z0.max() + 0.03)] = np.nan
            a1.plot(z, med / ref, color=c, ls=ls, lw=lw, **mko)
            a2.plot(z, w, color=c, ls=ls, lw=lw, **mko)
            a3.plot(lam, med, color=c, ls=ls, lw=lw, **mko)
        f.write("# -ln<F> of all pixels at the snapshots (redshift space) and lambda_mfp:\n"
                "# label z tau_eff lambda_mfp[pMpc] lambda_is_fit\n")
        for lab, _, _, (zi, ti, li, fi) in runs:
            for row in zip(zi, ti, li, fi):
                f.write(f"#S {lab:10s} {row[0]:.4f} {row[1]:.4g} {row[2]:.4g} {int(row[3])}\n")

    if a.norm:
        import igm_stats as ig
        ob = ig.read_obs()["mean_flux"]
        k = (ob["z"] >= 4.95) & (ob["z"] <= 6.25)
        a0.errorbar(ob["z"][k], -np.log(ob["v"][k]),
                    yerr=[np.log(ob["v"][k] + ob["ep"][k]) - np.log(ob["v"][k]),
                          np.log(ob["v"][k]) - np.log(np.maximum(ob["v"][k] - ob["em"][k], 1e-4))],
                    color="0.4", marker="s", ms=3, mfc="none", ls="none", lw=0.8,
                    label=r"$-\ln\langle F\rangle$ Bosman+22")
    elif not full and not snap_only:
        a0.plot([], [], color="0.4", marker="o", mfc="none", ls="none", label=r"$-\ln\langle F\rangle$ (snapshots)")
    if full:
        a0.plot([], [], color="0.4", marker="o", ms=4.5, ls="none",
                label="full box (snapshots, z-space)")
        a0.plot([], [], color="0.4", ls="-", label="slabs (real space)")
    a0.set_yscale("log")
    a0.set_ylabel(r"median $\tau_{\rm eff}$ (50 Mpc/$h$, " + kind + ")")
    a0.legend(frameon=False, fontsize=7, ncol=2)
    a1.axhline(1, color="0.6", lw=0.8)
    a1.set_ylabel(f"median / median({runs[0][0]})")
    a2.set_yscale("log")
    a2.set_ylabel(r"$\tau_{\rm eff}$(84 %) $-$ $\tau_{\rm eff}$(16 %)")
    for x in (a0, a1, a2):
        x.set_xlabel("z")
    a3.set_xscale("log")
    a3.set_yscale("log")
    a3.set_xlabel(r"$\lambda_{\rm mfp}$ [pMpc] (912 Å" + (")" if snap_only else ", interpolated between snapshots)"))
    a3.set_ylabel(r"median $\tau_{\rm eff}$ (50 Mpc/$h$, " + kind + ")")
    from matplotlib.ticker import ScalarFormatter, NullFormatter
    for x in (a0, a2, a3):
        x.yaxis.set_major_formatter(ScalarFormatter())
        x.yaxis.set_minor_formatter(NullFormatter())
    a0.set_yticks([1.5, 2, 3, 4, 6])
    a3.set_yticks([1.5, 2, 3, 4, 6])
    a2.set_yticks([0.5, 0.7, 1, 1.5, 2, 3])
    a3.xaxis.set_major_formatter(ScalarFormatter())
    a3.xaxis.set_minor_formatter(NullFormatter())
    a3.set_xticks([1, 2, 5, 10, 20, 50])
    for x in ax.flat:
        x.grid(True, color="0.9", lw=0.6)
        for sp in ("top", "right"):
            x.spines[sp].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(a.outdir, f"taueff_z{sfx}.png"), dpi=a.dpi)
    plt.close(fig)
    mfp_width_figure(runs, a.outdir, a.dpi, sfx, kind, full, snap_only)


def mfp_width_figure(runs, outdir, dpi, sfx="", kind="raw", full=None, snap_only=False):
    """Left: lambda_mfp against z at the snapshots; right: the raw 84-16 % width."""
    import matplotlib.pyplot as plt
    from matplotlib.ticker import ScalarFormatter, NullFormatter
    from rt_suite_figures import COLORS, STYLES
    fig, (b0, b1) = plt.subplots(1, 2, figsize=(9.0, 3.8), sharex=True)
    zlo = min(r[1].min() for r in runs) - 0.03
    zhi = max(r[1].max() for r in runs) + 0.03
    for ir, (lab, zs, ps, (zi, ti, li, fi)) in enumerate(runs):
        c, ls = COLORS[ir % len(COLORS)], STYLES[ir % len(STYLES)]
        lw = 1.8 if ir == 0 else 1.1
        k = (zi >= zlo) & (zi <= zhi) & np.isfinite(li)
        b0.plot(zi[k], li[k], color=c, ls=ls, lw=lw, label=lab)
        b0.plot(zi[k & ~fi], li[k & ~fi], color=c, marker="o", ms=3.5, ls="none")
        b0.plot(zi[k & fi], li[k & fi], color=c, marker="o", ms=3.5, mfc="none", ls="none")
        m = np.all(np.isfinite(ps[:, 1:4]), axis=1)
        z, p = (zs[m], ps[m]) if snap_only else smooth121(zs[m], ps[m])
        mko = dict(marker="os^vDph"[ir % 7], ms=4.5, mec="white", mew=0.5) if snap_only else {}
        b1.plot(z, p[:, 3] - p[:, 1], color=c, ls=ls, lw=lw, **mko)
        fb = (full or {}).get(lab)
        if fb is not None:
            fb = fb[np.all(np.isfinite(fb[:, 1:]), axis=1)]
            b1.plot(fb[:, 0], fb[:, 4] - fb[:, 2], color=c, marker="os^vDph"[ir % 7], ms=4.5, ls="none",
                    mec="white", mew=0.5, zorder=5)
    if full:
        b1.plot([], [], color="0.4", marker="o", ms=4.5, ls="none", label="full box (snapshots, z-space)")
        b1.plot([], [], color="0.4", ls="-", label="slabs (real space)")
        b1.legend(frameon=False, fontsize=7)
    b0.set_yscale("log")
    b0.set_ylabel(r"$\lambda_{\rm mfp}$ (912 Å) [pMpc]")
    b0.set_yticks([1, 2, 5, 10, 20, 50])
    b0.plot([], [], color="0.4", marker="o", mfc="none", ls="none", label="exp. fit (λ > rays)")
    b0.legend(frameon=False, fontsize=7, ncol=2)
    b1.set_yscale("log")
    b1.set_ylabel(r"$\tau_{\rm eff}$(84 %) $-$ $\tau_{\rm eff}$(16 %)" + "\n50 Mpc/$h$, " + kind)
    b1.set_yticks([0.5, 0.7, 1, 1.5, 2, 3])
    for x in (b0, b1):
        x.yaxis.set_major_formatter(ScalarFormatter())
        x.yaxis.set_minor_formatter(NullFormatter())
        x.set_xlabel("z")
        x.set_xlim(zlo, zhi)
        x.grid(True, color="0.9", lw=0.6)
        for sp in ("top", "right"):
            x.spines[sp].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, f"taueff_mfp_width{sfx}.png"), dpi=dpi)
    plt.close(fig)


if __name__ == "__main__":
    main()
