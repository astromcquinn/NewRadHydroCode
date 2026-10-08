#!/usr/bin/env python3
"""Science sanity checks of an apps/radhydro run (docs/full_run.md).

usage: full_run_report.py RUN_DIR PROBLEM_DIR [--figdir DIR] [--prefix NAME]

Reads RUN_DIR/step_*/stats.txt, timing.txt, Pk/pk.*, the on-the-fly halo
catalogs (write_catalogs directory, legacy format) and the fields of the
last step, and makes (if matplotlib is available) small PNG figures:
  <prefix>_xhii.png    global ionized fraction (volume/mass weighted) vs z
  <prefix>_hmf.png     FoF halo mass function at a few z vs Sheth-Tormen
  <prefix>_pk.png      DM power spectrum vs linear theory at a few z
  <prefix>_trho.png    T - rho/<rho> of the gas at the last step
  <prefix>_timing.png  time per step by component
and prints the numbers. Linear theory: Eisenstein-Hu no-wiggle P(k) and the
growth factor of tools/ics/pmcosmo.py (the ICs' own cosmology).
"""
import argparse
import glob
import os
import struct
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "ics"))
import pmcosmo as pc  # noqa: E402

MSUN_G = 1.98900e33
MPC_CM = 3.0856775807e24
RHOC = 1.87890e-29  # h^2 g/cm^3 (physical_constants.f90 rhoc_cgs)


def kv(fn):
    d = {}
    for line in open(fn):
        if "=" in line:
            k, v = line.split("=", 1)
            d[k.strip()] = v.strip()
    return d


def params(problem):
    p = kv(os.path.join(problem, "params.txt"))
    return {k: v.split("#")[0].strip() for k, v in p.items()}


def read_stats(run):
    rows = []
    for d in sorted(glob.glob(os.path.join(run, "step_*"))):
        s = kv(os.path.join(d, "stats.txt"))
        rows.append({k: float(v.split()[0].replace("nan", "nan")) if v else np.nan
                     for k, v in s.items()})
    return rows


def read_catalog(fn):
    data = open(fn, "rb").read()
    n = struct.unpack("<q", data[:8])[0]
    recs = []
    for i in range(n):
        r = struct.unpack("<6fi7f", data[8 + 56 * i: 8 + 56 * (i + 1)])
        recs.append(r)
    return recs


def sigma_of_M(M, rho_m):
    """sigma(M) at z = 0, M in Msun/h, rho_m in (Msun/h)/(Mpc/h)^3."""
    R = (3 * M / (4 * np.pi * rho_m)) ** (1 / 3.0)
    return np.array([pc.sigma_R(r, pc.power_z0, n=20000) for r in np.atleast_1d(R)])


def sheth_tormen_dndlnM(M, z, g):
    """dn/dlnM in (h/Mpc)^3, M in Msun/h."""
    rho_m = 2.77536627e11 * pc.OMEGAM
    D = g(1 / (1 + z)) / g(1.0)
    M = np.atleast_1d(M)
    s0 = sigma_of_M(M, rho_m)
    s0b = sigma_of_M(M * 1.05, rho_m)
    dlns_dlnM = (np.log(s0b) - np.log(s0)) / np.log(1.05)
    nu = (1.686 / (s0 * D)) ** 2
    A, a, p = 0.3222, 0.707, 0.3
    fnu = A * np.sqrt(2 * a * nu / np.pi) * (1 + (a * nu) ** -p) * np.exp(-a * nu / 2)
    return rho_m / M * fnu * np.abs(dlns_dlnM)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("problem")
    ap.add_argument("--figdir", default=None)
    ap.add_argument("--prefix", default="full")
    ap.add_argument("--catalogs", default=None, help="on-the-fly catalog directory")
    ap.add_argument("--hmf-z", default="", help="comma-separated catalog redshifts to plot")
    a = ap.parse_args()
    P = params(a.problem)
    box = float(P.get("box", 5))
    N = int(P.get("N", 32))
    ndm = int(P.get("Ndm", N // 2))
    st = read_stats(a.run)
    g = pc.Growth(1.0)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        plt = None
    figdir = a.figdir
    if figdir:
        os.makedirs(figdir, exist_ok=True)

    def save(fig, name):
        if figdir:
            fig.savefig(os.path.join(figdir, f"{a.prefix}_{name}.png"), dpi=80)
            plt.close(fig)

    # 1. ionized fraction
    z = np.array([s["redshift"] for s in st])
    xv = np.array([s.get("xHII_V", np.nan) for s in st])
    xm = np.array([s.get("xHII_M", np.nan) for s in st])
    ns = np.array([s.get("Nsource", 0) for s in st])
    rt_on = np.array([s.get("idata1", 1) for s in st]) > 1
    print("# ionized fraction (steps with RT on)")
    print("#   z        xHII_V      xHII_M     Nsource   T_V")
    for s, on in zip(st, rt_on):
        if on:
            print(f"  {s['redshift']:8.4f}  {s['xHII_V']:.4e}  {s['xHII_M']:.4e}  "
                  f"{int(s['Nsource']):6d}  {s['T_V']:.4e}")
    if plt and figdir:
        fig, ax = plt.subplots(figsize=(4.5, 3.2))
        m = rt_on & np.isfinite(xv)
        ax.semilogy(z[m], xv[m], "o-", ms=3, label="volume weighted")
        ax.semilogy(z[m], xm[m], "s--", ms=3, label="mass weighted")
        ax.set_xlabel("z")
        ax.set_ylabel("x_HII")
        ax.invert_xaxis()
        ax.legend(fontsize=8)
        ax.set_title(f"N={N}, box={box:g} Mpc/h", fontsize=9)
        fig.tight_layout()
        save(fig, "xhii")

    # 2. halo mass function
    cats = sorted(glob.glob(os.path.join(a.catalogs or "", "halo.z=*"))) if a.catalogs else []
    zsel = [float(x) for x in a.hmf_z.split(",") if x]
    mpart = RHOC * pc.OMEGAM * (box * MPC_CM) ** 3 / pc.HUBBLE0 / ndm ** 3 / MSUN_G  # Msun
    vol = box ** 3  # (Mpc/h)^3
    if cats:
        byz = {float(os.path.basename(c).split("=")[1]): c for c in cats}
        if not zsel:
            zs = sorted(byz)
            zsel = [zs[0], zs[len(zs) // 2], zs[-1]] if len(zs) > 2 else zs
        if plt and figdir:
            fig, ax = plt.subplots(figsize=(4.5, 3.2))
        print("# halo mass function: z, N_halos, M range [Msun], counts per dex vs ST")
        for zc in zsel:
            zc = min(byz, key=lambda q: abs(q - zc))
            recs = read_catalog(byz[zc])
            M = np.array([r[7] for r in recs]) * pc.HUBBLE0  # Msun/h
            if len(M) == 0:
                print(f"  z={zc}: no halos")
                continue
            bins = np.arange(np.floor(np.log10(M.min())), np.log10(M.max()) + 0.25, 0.25)
            cnt, _ = np.histogram(np.log10(M), bins)
            mc = 10 ** (0.5 * (bins[1:] + bins[:-1]))
            st_n = sheth_tormen_dndlnM(mc, zc, g) * np.log(10) * 0.25 * vol
            print(f"  z={zc:.4f}: {len(M)} halos, {M.min() / pc.HUBBLE0:.3g} - "
                  f"{M.max() / pc.HUBBLE0:.3g} Msun")
            for c, s_, m_ in zip(cnt, st_n, mc):
                print(f"     log10 M[Msun/h] = {np.log10(m_):5.2f}: FoF {c:4d}  ST {s_:8.2f}")
            if plt and figdir:
                ok = cnt > 0
                line = ax.errorbar(np.log10(mc[ok]), cnt[ok] / (0.25 * vol),
                                   yerr=np.sqrt(cnt[ok]) / (0.25 * vol), fmt="o", ms=3,
                                   label=f"FoF z={zc:.2f}")
                mm = np.logspace(bins[0], bins[-1], 30)
                ax.plot(np.log10(mm), sheth_tormen_dndlnM(mm, zc, g) * np.log(10), "-",
                        color=line[0].get_color(), lw=1)
        if plt and figdir:
            ax.set_yscale("log")
            ax.set_xlabel("log10 M [Msun/h]")
            ax.set_ylabel("dn/dlog10M [(h/Mpc)^3]")
            ax.axvline(np.log10(20 * mpart * pc.HUBBLE0), color="gray", ls=":", lw=0.8)
            ax.legend(fontsize=7)
            ax.set_title("FoF (points) vs Sheth-Tormen (lines)", fontsize=9)
            fig.tight_layout()
            save(fig, "hmf")

    # 3. power spectrum
    pks = sorted(glob.glob(os.path.join(a.run, "Pk", "pk.z=*")))
    if pks:
        print("# DM power spectrum Delta^2(k) (CIC-deconvolved) / linear, lowest 4 bins")
        sel = [pks[0], pks[len(pks) // 2], pks[-1]] if len(pks) > 2 else pks
        if plt and figdir:
            fig, ax = plt.subplots(figsize=(4.5, 3.2))
        for f in sel:
            zz = float(os.path.basename(f).split("=")[1])
            d = np.loadtxt(f)
            k = d[:, 0]
            lin = k ** 3 * pc.power_z0(k) / (2 * np.pi ** 2) * (g(1 / (1 + zz)) / g(1.0)) ** 2
            print(f"  z={zz:.4f}: " + " ".join(f"k={kk:.3g}:{r:.3f}"
                                               for kk, r in zip(k[:4], d[:4, 3] / lin[:4])))
            if plt and figdir:
                line = ax.loglog(k, d[:, 3], "o", ms=3, label=f"DM z={zz:.2f}")
                ax.loglog(k, d[:, 5], "x", ms=3, color=line[0].get_color())
                ax.loglog(k, lin, "-", lw=1, color=line[0].get_color())
        if plt and figdir:
            ax.set_xlabel("k [h/Mpc]")
            ax.set_ylabel("Delta^2(k)")
            ax.legend(fontsize=7)
            ax.set_title("DM (o), gas (x), linear (lines)", fontsize=9)
            fig.tight_layout()
            save(fig, "pk")

    # 4. T - rho at the last step with fields
    last = None
    for d in sorted(glob.glob(os.path.join(a.run, "step_*")))[::-1]:
        if os.path.exists(os.path.join(d, "temp.bin")):
            last = d
            break
    if last:
        T = np.fromfile(os.path.join(last, "temp.bin"), dtype="<f4")
        rho = np.fromfile(os.path.join(last, "rhog.bin"), dtype="<f4")
        fHI = np.fromfile(os.path.join(last, "fHI.bin"), dtype="<f4")
        s = kv(os.path.join(last, "stats.txt"))
        delta = rho / rho.mean()
        ion = fHI < 0.5
        print(f"# T-rho at z = {float(s['redshift']):.4f} ({os.path.basename(last)}): "
              f"ionized cells {ion.mean():.3f}")
        for lab, m in (("ionized (fHI<0.5)", ion), ("neutral", ~ion)):
            if m.sum() > 10:
                x = np.log10(delta[m])
                y = np.log10(T[m])
                sel = (x > -0.5) & (x < 0.5)
                if sel.sum() > 10:
                    gam = np.polyfit(x[sel], y[sel], 1)
                    print(f"  {lab}: T0 = {10 ** gam[1]:.4g} K, gamma-1 = {gam[0]:.3f} "
                          f"({sel.sum()} cells with |log delta| < 0.5)")
        if plt and figdir:
            fig, ax = plt.subplots(figsize=(4.5, 3.2))
            ax.hist2d(np.log10(delta), np.log10(np.maximum(T, 1)), bins=60, cmap="viridis",
                      cmin=1)
            ax.set_xlabel("log10 rho_gas/<rho_gas>")
            ax.set_ylabel("log10 T [K]")
            ax.set_title(f"z = {float(s['redshift']):.2f}", fontsize=9)
            fig.tight_layout()
            save(fig, "trho")

    # 5. timing
    tf = os.path.join(a.run, "timing.txt")
    if os.path.exists(tf):
        t = np.loadtxt(tf, comments="#")
        t = np.atleast_2d(t)
        names = ["hydro", "pm", "halos", "rt", "summaries", "io"]
        tot = t[:, 5:11].sum(axis=0)
        print("# timing totals [s]: " + ", ".join(f"{n} {v:.1f}" for n, v in zip(names, tot))
              + f", all {t[:, 11].sum():.1f}")
        rton = t[:, 8] > 0
        if rton.any():
            sub = t[rton]
            print(f"# RT: {int(sub[:, 2].sum())} substeps, {sub[:, 8].sum():.1f} s, "
                  f"{1e3 * sub[:, 8].sum() / sub[:, 2].sum():.1f} ms per substep; "
                  f"rays at the end of the RT steps: "
                  + ", ".join(str(int(x)) for x in sub[:, 4]))
        if plt and figdir:
            fig, ax = plt.subplots(figsize=(4.5, 3.2))
            bottom = np.zeros(len(t))
            for i, n in enumerate(names[:4]):
                ax.bar(t[:, 0], t[:, 5 + i], bottom=bottom, label=n)
                bottom += t[:, 5 + i]
            ax.set_yscale("symlog", linthresh=0.1)
            ax.set_xlabel("step")
            ax.set_ylabel("time per step [s]")
            ax.legend(fontsize=7)
            fig.tight_layout()
            save(fig, "timing")


if __name__ == "__main__":
    main()
