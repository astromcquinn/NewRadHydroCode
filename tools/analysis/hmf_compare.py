#!/usr/bin/env python3
"""Halo mass function of on-the-fly FoF catalogs vs Sheth-Tormen and Tinker.

usage: hmf_compare.py CATALOG_DIR --box 10 --N 128 [--z 6 8 10 12]
                      [--hmf-table build/hmf_table] [--out fig.png]

CATALOG_DIR holds the legacy-format catalogs halo.z=XX.XXXX written by
apps/radhydro with write_catalogs (docs/halos_and_sources.md): an int64
Nhalo followed by 56-byte source_type records (M in Msun). For each
requested z the nearest catalog is used. The analytic curves come from
apps/hmf_table (the same code as the abundance matching). Prints a table of
the binned counts and the ratios to ST and Tinker, and writes a figure.
"""
import argparse
import glob
import os
import re
import subprocess
import sys

import numpy as np

trapezoid = np.trapezoid if hasattr(np, "trapezoid") else np.trapz  # numpy < 2.0

REC = np.dtype([("x", "<f4", 3), ("v", "<f4", 3), ("N", "<i4"), ("M", "<f4"), ("R", "<f4"),
                ("T", "<f4"), ("L", "<f4"), ("sfr", "<f4"), ("mstar", "<f4"), ("fgas", "<f4")])
assert REC.itemsize == 56


def read_catalog(path):
    with open(path, "rb") as f:
        n = int(np.frombuffer(f.read(8), "<i8")[0])
        h = np.frombuffer(f.read(56 * n), REC)
    return h


def catalogs(d):
    out = []
    for p in glob.glob(os.path.join(d, "halo.z=*")):
        m = re.search(r"halo\.z=([0-9.]+)$", p)
        if m:
            out.append((float(m.group(1)), p))
    return sorted(out)


def analytic(exe, z, lo, hi, d=0.05):
    r = subprocess.run([exe, str(z), str(lo), str(hi), str(d)], check=True, capture_output=True,
                       text=True)
    a = np.loadtxt(r.stdout.splitlines())
    return a  # logM, dn/dlog10M PS ST T, n(>M) PS ST T


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("catdir")
    ap.add_argument("--box", type=float, required=True, help="Mpc/h")
    ap.add_argument("--N", type=int, required=True, help="Ndm per side")
    ap.add_argument("--z", type=float, nargs="+", default=[6, 8, 10, 12])
    ap.add_argument("--hmf-table", default="build/hmf_table")
    ap.add_argument("--nmin", type=int, default=20)
    ap.add_argument("--dlog", type=float, default=0.25)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    h = 0.7
    vol = (a.box / h) ** 3  # Mpc^3
    # particle mass, Msun: rho_crit Omega_m (box)^3 / N^3 (all matter, docs)
    mp = 2.775e11 * 0.27 * a.box ** 3 / a.N ** 3 / h
    cats = catalogs(a.catdir)
    if not cats:
        sys.exit("no catalogs in " + a.catdir)
    rows = []
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    nz = len(a.z)
    fig, axes = plt.subplots(2, nz, figsize=(3.6 * nz, 6), sharex=True,
                             gridspec_kw={"height_ratios": [2, 1]}, constrained_layout=True)
    axes = np.atleast_2d(axes).reshape(2, nz)
    print(f"# particle mass {mp:.3e} Msun, volume {vol:.1f} Mpc^3, N >= {a.nmin}")
    print("# z_cat logM_lo logM_hi count dn/dlog10M  ratio_ST  ratio_Tinker  ratio_ST_corr")
    for iz, zt in enumerate(a.z):
        zc, path = min(cats, key=lambda c: abs(c[0] - zt))
        cat = read_catalog(path)
        cat = cat[cat["N"] >= a.nmin]
        logM = np.log10(cat["M"].astype(np.float64))
        # Warren et al. (2006) FoF low-N correction N -> N (1 - N^-0.6), shown as an alternative
        logMc = np.log10(cat["N"] * (1 - cat["N"].astype(float) ** -0.6) * mp)
        lo = np.log10(a.nmin * mp)
        edges = np.arange(lo, logM.max() + a.dlog, a.dlog) if len(logM) else np.array([lo, lo + a.dlog])
        cnt, _ = np.histogram(logM, edges)
        cntc, _ = np.histogram(logMc, edges)
        mid = 0.5 * (edges[1:] + edges[:-1])
        an = analytic(a.hmf_table, zc, lo - 0.3, max(edges[-1], lo + 1) + 0.3)
        # Bin-averaged analytic counts (integrate dn/dlog10M over the bin)
        def binavg(col):
            out = []
            for e0, e1 in zip(edges[:-1], edges[1:]):
                s = (an[:, 0] >= e0 - 1e-9) & (an[:, 0] <= e1 + 1e-9)
                out.append(trapezoid(an[s, col], an[s, 0]) / (e1 - e0))
            return np.array(out)
        st, ti = binavg(2), binavg(3)
        dn = cnt / vol / a.dlog
        dnc = cntc / vol / a.dlog
        err = np.sqrt(np.maximum(cnt, 1)) / vol / a.dlog
        for k in range(len(cnt)):
            if cnt[k] == 0:
                continue
            print(f"{zc:.3f} {edges[k]:.2f} {edges[k+1]:.2f} {cnt[k]:6d} {dn[k]:.3e} "
                  f"{dn[k]/st[k]:.3f} {dn[k]/ti[k]:.3f} {dnc[k]/st[k]:.3f}")
            rows.append((zc, edges[k], cnt[k], dn[k] / st[k], dn[k] / ti[k]))
        ax, rx = axes[0, iz], axes[1, iz]
        ax.plot(an[:, 0], np.log10(an[:, 2]), "-", color="0.25", lw=1.5, label="Sheth–Tormen")
        ax.plot(an[:, 0], np.log10(an[:, 3]), "--", color="0.25", lw=1.5, label="Tinker 08")
        ok = cnt > 0
        ax.errorbar(mid[ok], np.log10(dn[ok]),
                    yerr=[np.log10(dn[ok]) - np.log10(np.maximum(dn[ok] - err[ok], 1e-12)),
                          np.log10(dn[ok] + err[ok]) - np.log10(dn[ok])],
                    fmt="o", ms=5, color="#2a6fbd", label=f"FoF, {a.N}³")
        ax.set_title(f"z = {zc:.2f}  ({len(cat)} halos)", fontsize=10)
        ax.set_ylabel(r"log$_{10}$ dn/dlog$_{10}$M [Mpc$^{-3}$]" if iz == 0 else "")
        rx.axhline(1, color="0.6", lw=1)
        rx.errorbar(mid[ok], dn[ok] / st[ok], yerr=err[ok] / st[ok], fmt="o", ms=4, color="#2a6fbd",
                    label="FoF / ST")
        rx.errorbar(mid[ok] + 0.03, dn[ok] / ti[ok], yerr=err[ok] / ti[ok], fmt="s", ms=4, mfc="none",
                    color="#2a6fbd", label="FoF / Tinker")
        rx.set_ylim(0, 2.5)
        rx.set_xlabel(r"log$_{10}$ M [M$_\odot$]")
        rx.set_ylabel("ratio" if iz == 0 else "")
        if iz == 0:
            ax.legend(frameon=False, fontsize=8)
            rx.legend(frameon=False, fontsize=8)
        for x in (ax, rx):
            x.grid(alpha=0.2, lw=0.5)
            for sp in ("top", "right"):
                x.spines[sp].set_visible(False)
    out = a.out or os.path.join(a.catdir, "hmf_compare.png")
    fig.savefig(out, dpi=130)
    print("#", out)


if __name__ == "__main__":
    main()
