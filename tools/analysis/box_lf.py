#!/usr/bin/env python3
"""UV luminosity function of the simulated box vs the observed LFs.

    box_lf.py CATALOG_DIR AM_PARAMS LF_TABLE --box 10 --N 128 [--out fig.png]

Applies the fitted L_UV(M_halo, z) of AM_PARAMS (apps/am_fit output; the
formula of include/radray/abundance_matching.hpp LuvModel) to the FoF halos
of each catalog (N >= 20 particles), and compares the binned LF of the box
with the LF_TABLE rows (data/uv_lf/lf_params.txt) nearest in redshift.
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hmf_compare import catalogs, read_catalog  # noqa: E402


def read_params(path):
    p = {}
    for l in open(path):
        l = l.split("#")[0].strip()
        if "=" in l:
            k, v = [x.strip() for x in l.split("=", 1)]
            p[k] = v
    return p


def log_luv(p, logM, z):
    g = lambda k, d=0.0: float(p.get(k, d))
    x = np.clip(z, g("am_zlo", 0), g("am_zhi", 1e30)) - g("am_z0", 8)
    lLp = g("am_l0") + g("am_l1") * x + g("am_l2") * x * x
    a = g("am_a0") + g("am_a1") * x
    lMp = g("am_m0") + g("am_m1") * x
    b = g("am_b")
    u = logM - lMp
    ea, eb = -a * u, -b * u
    mx = np.maximum(ea, eb)
    out = lLp + np.log10(2) - mx - np.log10(10 ** (ea - mx) + 10 ** (eb - mx))
    zh = g("am_zhi", 1e30)
    if "am_norm_z" in p and z > zh:
        kz = [float(v) for v in p["am_norm_z"].split(",")]
        kd = [float(v) for v in p["am_norm_dl"].split(",")]
        out = out + np.interp(z, [zh] + kz, [0.0] + kd)
    return out


def lf_row(r, M):
    z, form, Ms, ps, al, be = r[:6]
    dM = M - Ms
    if form == "schechter":
        xx = 10 ** (-0.4 * dM)
        return 0.4 * np.log(10) * ps * xx ** (al + 1) * np.exp(-xx)
    return ps / (10 ** (0.4 * (al + 1) * dM) + 10 ** (0.4 * (be + 1) * dM))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("catdir")
    ap.add_argument("am_params")
    ap.add_argument("lf_table")
    ap.add_argument("--box", type=float, required=True)
    ap.add_argument("--N", type=int, required=True)
    ap.add_argument("--nmin", type=int, default=20)
    ap.add_argument("--out")
    a = ap.parse_args()
    p = read_params(a.am_params)
    rows = []
    for l in open(a.lf_table):
        t = l.split("#")[0].split()
        if len(t) >= 11 and float(t[8]) > 0:
            rows.append((float(t[0]), t[1], float(t[2]), float(t[3]), float(t[4]), float(t[5]),
                         float(t[6]), float(t[7]), t[9]))
    vol = (a.box / 0.7) ** 3
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    cats = [c for c in catalogs(a.catdir) if c[0] < 11.5]
    sel = [min(cats, key=lambda c: abs(c[0] - r[0])) for r in rows]
    pairs = [(r, c) for r, c in zip(rows, sel) if abs(r[0] - c[0]) < 0.6]
    fig, axes = plt.subplots(1, len(pairs), figsize=(3.6 * len(pairs), 3.8), constrained_layout=True)
    axes = np.atleast_1d(axes)
    print("# z_LF ref z_cat  M_UV  phi_box  phi_obs  ratio  count")
    for ax, (r, (zc, path)) in zip(axes, pairs):
        h = read_catalog(path)
        h = h[h["N"] >= a.nmin]
        lL = log_luv(p, np.log10(h["M"].astype(float)), zc)
        Muv = 51.595 - 2.5 * lL
        edges = np.arange(-23, -9.99, 0.5)
        cnt, _ = np.histogram(Muv, edges)
        mid = 0.5 * (edges[1:] + edges[:-1])
        phi = cnt / vol / 0.5
        ok = cnt > 0
        Mg = np.linspace(-23, -10, 200)
        ax.plot(Mg, np.log10(lf_row(r, Mg)), "--", color="0.3", lw=1.5, label=f"{r[8]} z={r[0]:g}")
        obs = (Mg >= r[6]) & (Mg <= r[7])
        ax.plot(Mg[obs], np.log10(lf_row(r, Mg[obs])), "-", color="0.3", lw=3, alpha=0.4)
        ax.errorbar(mid[ok], np.log10(phi[ok]), yerr=0.434 / np.sqrt(cnt[ok]), fmt="o", ms=4,
                    color="#2a6fbd", label=f"box, z={zc:.2f}")
        for m, ph, c in zip(mid[ok], phi[ok], cnt[ok]):
            print(f"{r[0]:g} {r[8]} {zc:.2f} {m:.2f} {ph:.3e} {lf_row(r, m):.3e} {ph / lf_row(r, m):.3f} {c}")
        ax.set_xlim(-23, -10)
        ax.set_ylim(-6, 0.5)
        ax.invert_xaxis()
        ax.set_xlabel("M$_{UV}$")
        ax.set_ylabel(r"log$_{10}$ φ [Mpc$^{-3}$ mag$^{-1}$]")
        ax.legend(frameon=False, fontsize=8, loc="lower left")
        ax.grid(alpha=0.2, lw=0.5)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    out = a.out or os.path.join(a.catdir, "box_lf.png")
    fig.savefig(out, dpi=130)
    print("#", out)


if __name__ == "__main__":
    main()
