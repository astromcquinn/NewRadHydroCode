#!/usr/bin/env python3
"""FoF mass functions of several runs at matched redshifts (docs/pm_refine.md).

usage: hmf_runs_compare.py --run LABEL:CATDIR:BOX:NDM [--run ...] --ref LABEL
                           [--z 5.5 6 7 8 10] [--dlogm 0.25] [--nmin 20]
                           [--hmf-table build/hmf_table] [--out fig.png]
                           [--table out.txt]

Each run is a directory of on-the-fly catalogs halo.z=XX.XXXX (radhydro
write_catalogs, read as in hmf_compare.py), its box in Mpc/h and Ndm. For
each requested z every run uses its nearest catalog (the redshifts differ by
the time steps; the difference is listed); with --interp each run's dn/dlogM
is interpolated in log n, linearly in z, between the two catalogs that
bracket z (the counts shown are the nearest catalog's). Bins of dlogm dex
start at the
--nmin-particle mass of the finest run. Prints dn/dlog10M, the ratio of each
run to the --ref run with the Poisson error of both, the number of halos
per bin, and the ratio to Sheth-Tormen (apps/hmf_table). Writes a figure:
dn/dlogM (top) and the ratio to the reference (bottom) per redshift.
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hmf_compare import analytic, catalogs, read_catalog  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "ics"))
import pmcosmo as pc  # noqa: E402

MSUN_PER_RHOC = 2.775e11  # rho_crit h^-2 Msun Mpc^-3 (only for the particle mass label)


def particle_mass(box, ndm):
    # cosmo_tools.f90 mp_unit (Msun): rhoc om (box/h)^3 / Npart, as in halos.cpp
    h = pc.HUBBLE0
    return MSUN_PER_RHOC * h * h * pc.OMEGAM * (box / h) ** 3 / ndm ** 3


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", action="append", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--z", type=float, nargs="+", default=[5.5, 6.0, 7.0, 8.0, 10.0])
    ap.add_argument("--dlogm", type=float, default=0.25)
    ap.add_argument("--nmin", type=int, default=20)
    ap.add_argument("--hmf-table", default="build/hmf_table")
    ap.add_argument("--out", default=None)
    ap.add_argument("--table", default=None)
    ap.add_argument("--interp", action="store_true",
                    help="interpolate each run to the requested z (see above)")
    a = ap.parse_args()

    runs = []
    for spec in a.run:
        lab, d, box, ndm = spec.split(":")
        runs.append(dict(label=lab, dir=d, box=float(box), ndm=int(ndm), cats=catalogs(d),
                         mp=particle_mass(float(box), int(ndm))))
    ref = [r for r in runs if r["label"] == a.ref][0]
    mp_min = min(r["mp"] for r in runs)
    lo = np.log10(a.nmin * mp_min)
    lines = []

    def out(s=""):
        print(s)
        lines.append(s)

    for r in runs:
        out(f"# {r['label']}: {r['dir']}, box {r['box']:g} Mpc/h, Ndm {r['ndm']}, "
            f"m_p = {r['mp']:.3e} Msun, {a.nmin} particles = {a.nmin * r['mp']:.3e} Msun, "
            f"{len(r['cats'])} catalogs")
    results = {}
    for z in a.z:
        sel = {}
        brack = {}
        for r in runs:
            if not r["cats"]:
                continue
            zc, p = min(r["cats"], key=lambda c: abs(c[0] - z))
            h = read_catalog(p)
            m = h["M"].astype(float)
            m = m[h["N"] >= a.nmin]
            sel[r["label"]] = (zc, m, r)
            lo_c = [c for c in r["cats"] if c[0] <= z]
            hi_c = [c for c in r["cats"] if c[0] >= z]
            if a.interp and lo_c and hi_c:
                c1, c2 = max(lo_c), min(hi_c)
                ms = []
                for zz, pp in (c1, c2):
                    hh = read_catalog(pp)
                    ms.append((zz, hh["M"].astype(float)[hh["N"] >= a.nmin]))
                brack[r["label"]] = ms
        if ref["label"] not in sel:
            continue
        mmax = max(np.log10(s[1].max()) for s in sel.values() if s[1].size)
        edges = np.arange(lo, mmax + a.dlogm, a.dlogm)
        cen = 0.5 * (edges[1:] + edges[:-1])
        zref = sel[ref["label"]][0]
        st = analytic(a.hmf_table, zref, edges[0], edges[-1], 0.01)
        st_d = np.interp(cen, st[:, 0], st[:, 2])  # dn/dlog10M, Sheth-Tormen
        out()
        out(f"z = {z}: catalogs at " + ", ".join(f"{k} z={v[0]:.4f}" for k, v in sel.items())
            + ("  (interpolated to z: " + ", ".join(
                f"{k} {v[0][0]:.3f}-{v[1][0]:.3f}" for k, v in brack.items()) + ")" if brack else ""))
        hdr = " log M   " + "".join(f"{k:>24s}" for k in sel) + "".join(
            f"   {k}/{ref['label']:<10s}" for k in sel if k != ref["label"])
        out(hdr + "      ST")
        res = {}
        for k, (zc, m, r) in sel.items():
            n, _ = np.histogram(np.log10(m), edges)
            vol = (r["box"] / pc.HUBBLE0) ** 3  # Mpc^3 (h from cosmo_parameters)
            d = n / vol / a.dlogm
            if k in brack:
                (z1, m1), (z2, m2) = brack[k]
                n1 = np.histogram(np.log10(m1), edges)[0].astype(float)
                n2 = np.histogram(np.log10(m2), edges)[0].astype(float)
                w = 0.0 if z2 == z1 else (z - z1) / (z2 - z1)
                with np.errstate(divide="ignore", invalid="ignore"):
                    ni = np.where((n1 > 0) & (n2 > 0), n1 ** (1 - w) * n2 ** w, 0.0)
                d = ni / vol / a.dlogm
            res[k] = (n, d)
        nr, dr = res[ref["label"]]
        for b in range(cen.size):
            s = f" {cen[b]:6.3f} "
            for k in sel:
                n, d = res[k]
                s += f"  {d[b]:10.3e} ({n[b]:6d})   "
            for k in sel:
                if k == ref["label"]:
                    continue
                n, d = res[k]
                if nr[b] > 0 and n[b] > 0:
                    rat = d[b] / dr[b]
                    err = rat * np.sqrt(1 / n[b] + 1 / nr[b])
                    s += f"   {rat:6.3f} +- {err:5.3f}    "
                else:
                    s += "        --             "
            s += f"  {st_d[b]:10.3e}"
            out(s)
        results[z] = dict(cen=cen, res=res, st=st_d, sel={k: v[0] for k, v in sel.items()})

    if a.table:
        with open(a.table, "w") as f:
            f.write("\n".join(lines) + "\n")
    if a.out:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        zs = list(results)
        fig, ax = plt.subplots(2, len(zs), figsize=(3.6 * len(zs), 6), sharex=True,
                               gridspec_kw=dict(height_ratios=[2, 1]), squeeze=False)
        for j, z in enumerate(zs):
            R = results[z]
            cen = R["cen"]
            nr, dr = R["res"][ref["label"]]
            ax[0, j].plot(cen, R["st"], "k:", lw=1, label="Sheth-Tormen")
            for k, (n, d) in R["res"].items():
                ok = n > 0
                ax[0, j].errorbar(cen[ok], d[ok], yerr=d[ok] / np.sqrt(n[ok]), fmt="o-", ms=3,
                                  label=f"{k} (z={R['sel'][k]:.2f})")
                if k != ref["label"]:
                    ok2 = ok & (nr > 0)
                    rat = d[ok2] / dr[ok2]
                    err = rat * np.sqrt(1 / n[ok2] + 1 / nr[ok2])
                    ax[1, j].errorbar(cen[ok2], rat, yerr=err, fmt="o-", ms=3, label=k)
            ax[0, j].set_yscale("log")
            ax[0, j].set_title(f"z = {z}")
            ax[1, j].axhline(1, color="k", lw=0.8)
            ax[1, j].set_ylim(0, 2)
            ax[1, j].set_xlabel(r"log$_{10}$ M [M$_\odot$]")
            for r in runs:
                for axx in ax[:, j]:
                    axx.axvline(np.log10(a.nmin * r["mp"]), color="0.7", lw=0.6, ls="--")
        ax[0, 0].set_ylabel(r"dn/dlog$_{10}$M [Mpc$^{-3}$]")
        ax[1, 0].set_ylabel(f"ratio to {ref['label']}")
        ax[0, 0].legend(fontsize=7)
        ax[1, 0].legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(a.out, dpi=110)
        print("wrote", a.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
