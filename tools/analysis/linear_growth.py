#!/usr/bin/env python3
"""Linear growth of dark matter and gas from MUSIC ICs against CAMB.

    linear_growth.py ROOT [--z 100 50 30] [--out ROOT/linear_growth.png] [--json ROOT/linear_growth.json]

The test of tools/hyak/linear_growth.slurm (docs/thermal_history.md, "Linear
growth from MUSIC ICs"): a uniform box from MUSIC with CAMB's separate CDM and
baryon transfer functions at z_init, run with hydro and gravity only
(thermal_init = recfast) to the zdata1 outputs. Reads

  ROOT/prob/ics/xv_dm.init, u_gas.init     the ICs (z = z_init)
  ROOT/run/DM/xv_dm.<ztag>, Gas/u_gas.<ztag>   the snapshots
  ROOT/run/Pk/pk.<ztag>                        the code's own P(k) (cross-check)
  ROOT/prob/params.txt                         N, Ndm, box, z_init, cosmology

and computes P_c(k) and P_b(k) with one estimator for the ICs and the
snapshots: the DM by CIC on an N^3 mesh with the window deconvolved, the gas
from the cell densities; bins of unit width in |n| (k = 2 pi n / L), as the
code's calc_pk_3d. CAMB (the run's cosmology, Tcmb 2.725 K, YHe 0.24,
N_eff 3.046, no massive neutrinos, sigma_8 normalised at z = 0) gives the
linear P_cdm and P_b at z_init and the outputs.

Printed per output z and k bin:

  G_c, G_b    the growth since z_init of each species over CAMB's,
              [P_sim(z)/P_sim(z_init)] / [P_CAMB(z)/P_CAMB(z_init)]: the
              same modes, so the box's cosmic variance cancels; 1 for exact
              linear growth;
  R_bc        (P_b/P_c)_sim / (P_b/P_c)_CAMB: the baryons' catch-up and
              their pressure (Jeans) suppression against CAMB's, also free of
              the ICs' normalisation;
  A_c, A_b    P_sim/P_CAMB (the absolute amplitude; carries the box's
              cosmic variance and the ICs' normalisation convention);
  D2_c        CAMB's Delta^2_cdm (k^3 P / 2 pi^2): linear theory holds where
              it is << 1.
"""
import argparse
import glob
import json
import os
import re

import numpy as np


def params(path):
    p = {}
    with open(path) as f:
        for line in f:
            line = line.split("#")[0]
            if "=" in line:
                k, v = line.split("=", 1)
                p[k.strip()] = v.strip()
    return p


def ztag_files(run, sub, stem):
    out = {}
    for f in glob.glob(os.path.join(run, sub, stem + ".z=*")):
        out[float(re.search(r"z=\s*([0-9.]+)", f).group(1))] = f
    return out


def find_z(files, z):
    for zz, f in files.items():
        if abs(zz - z) < 1e-3 * (1 + z):
            return f
    raise SystemExit(f"no file at z = {z}: {sorted(files)}")


def cic(pos, n):
    """Mass assignment of positions in mesh cells (periodic) to an n^3 mesh."""
    rho = np.zeros(n * n * n)
    x = np.mod(pos, n)
    i0 = np.floor(x).astype(np.int64)
    d = x - i0
    for dx in (0, 1):
        wx = d[:, 0] if dx else 1 - d[:, 0]
        ix = (i0[:, 0] + dx) % n
        for dy in (0, 1):
            wy = d[:, 1] if dy else 1 - d[:, 1]
            iy = (i0[:, 1] + dy) % n
            for dz in (0, 1):
                wz = d[:, 2] if dz else 1 - d[:, 2]
                iz = (i0[:, 2] + dz) % n
                rho += np.bincount(ix + n * (iy + n * iz), weights=wx * wy * wz, minlength=n ** 3)
    return rho.reshape(n, n, n)  # [k][j][i]


def pk(delta, box, deconvolve):
    """P(k) [(Mpc/h)^3] in bins of unit |n| (1 .. n/2), and the bins' mean k."""
    n = delta.shape[0]
    dk = np.fft.rfftn(delta) / n ** 3
    f = np.fft.fftfreq(n) * n
    fz = np.fft.rfftfreq(n) * n
    nx, ny, nz = np.meshgrid(f, f, fz, indexing="ij")
    if deconvolve:  # the CIC window per axis: sinc^2(pi n_i / n)
        w = (np.sinc(nx / n) * np.sinc(ny / n) * np.sinc(nz / n)) ** 2
        dk = dk / w
    p = np.abs(dk) ** 2 * box ** 3
    nn = np.sqrt(nx ** 2 + ny ** 2 + nz ** 2)
    # the rfft half: the planes 0 < nz < n/2 stand for two modes each
    mult = np.where((nz == 0) | (nz == n // 2), 1.0, 2.0)
    b = np.rint(nn).astype(int)
    nb = n // 2
    sel = (b >= 1) & (b <= nb)
    cnt = np.bincount(b[sel], weights=mult[sel], minlength=nb + 1)[1:]
    ps = np.bincount(b[sel], weights=(mult * p)[sel], minlength=nb + 1)[1:] / cnt
    ks = np.bincount(b[sel], weights=(mult * nn)[sel], minlength=nb + 1)[1:] / cnt * 2 * np.pi / box
    return ks, ps


def read_dm(path, ndm, nm):
    xv = np.fromfile(path, dtype=np.float32).reshape(-1, 6).astype(np.float64)
    assert xv.shape[0] == ndm ** 3, (path, xv.shape)
    return cic(xv[:, :3] * (nm / ndm), nm)


def read_gas(path, n):
    u = np.fromfile(path, dtype=np.float32).reshape(n, n, n, 5)
    return u[..., 0].astype(np.float64)


def delta_of(rho):
    return rho / rho.mean() - 1


def camb_pk(cos, zs, kmax):
    import camb
    h = cos["hubble0"]
    p = camb.CAMBparams()
    p.set_cosmology(H0=100 * h, ombh2=cos["omegab"] * h * h, omch2=(cos["omegam"] - cos["omegab"]) * h * h,
                    TCMB=2.725, YHe=0.24, mnu=0, nnu=3.046, omk=0, tau=0.054)
    p.InitPower.set_params(ns=cos["ns"], As=2.1e-9)
    p.set_matter_power(redshifts=sorted(set(list(zs) + [0.0]), reverse=True), kmax=kmax * h * 1.2,
                       k_per_logint=30)
    p.NonLinear = camb.model.NonLinear_none
    r = camb.get_results(p)
    norm = (cos["sigma8"] / r.get_sigma8_0()) ** 2
    out = {}
    for var, tag in (("delta_cdm", "c"), ("delta_baryon", "b")):
        f = r.get_matter_power_interpolator(nonlinear=False, var1=var, var2=var, hubble_units=True,
                                            k_hunit=True, kmax=kmax * 1.1)
        out[tag] = lambda z, k, f=f: norm * f.P(z, k)
    return out, camb.__version__


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root")
    ap.add_argument("--z", type=float, nargs="+", default=[100, 50, 30])
    ap.add_argument("--out")
    ap.add_argument("--json")
    a = ap.parse_args()
    prob, run = os.path.join(a.root, "prob"), os.path.join(a.root, "run")
    P = params(os.path.join(prob, "params.txt"))
    n, ndm, box, z0 = int(P["N"]), int(P.get("Ndm", P["N"])), float(P["box"]), float(P["z_init"])
    cos = {k: float(P[k]) for k in ("omegam", "omegab", "hubble0", "sigma8", "ns")}
    ics = os.path.join(prob, P.get("ics", "ics"))
    kmax = np.pi * n / box * np.sqrt(3)
    camb, cver = camb_pk(cos, [z0] + a.z, kmax)

    def sim(dm_path, gas_path):
        kc, pc_ = pk(delta_of(read_dm(dm_path, ndm, n)), box, True)
        kb, pb_ = pk(delta_of(read_gas(gas_path, n)), box, False)
        return kc, pc_, pb_

    k, pc0, pb0 = sim(os.path.join(ics, "xv_dm.init"), os.path.join(ics, "u_gas.init"))
    lin0 = {s: camb[s](z0, k) for s in "cb"}
    dmf, gasf, pkf = ztag_files(run, "DM", "xv_dm"), ztag_files(run, "Gas", "u_gas"), ztag_files(run, "Pk", "pk")
    res = {"camb": cver, "N": n, "box": box, "z_init": z0, "cosmology": cos, "k": k.tolist(),
           "ics": {"A_c": (pc0 / lin0["c"]).tolist(), "A_b": (pb0 / lin0["b"]).tolist()}, "z": {}}
    print(f"N = {n}, box = {box} Mpc/h, z_init = {z0}, CAMB {cver}; cosmology {cos}")
    print(f"ICs (z = {z0}): A_c = P_c/P_CAMB, A_b, R_bc at k = {k[0]:.0f}, {k[n // 8]:.0f}, {k[n // 4]:.0f}, {k[-1]:.0f} h/Mpc: "
          + "  ".join(f"{pc0[i] / lin0['c'][i]:.3f}/{pb0[i] / lin0['b'][i]:.3f}/{(pb0[i] / pc0[i]) / (lin0['b'][i] / lin0['c'][i]):.3f}"
                      for i in (0, n // 8, n // 4, n // 2 - 1)))
    show = sorted(set([0, 1, 3, 7, 15, 31, n // 4, 3 * n // 8, n // 2 - 1]) & set(range(n // 2)))
    for z in a.z:
        kz, pc, pb = sim(find_z(dmf, z), find_z(gasf, z))
        lin = {s: camb[s](z, k) for s in "cb"}
        Gc = (pc / pc0) / (lin["c"] / lin0["c"])
        Gb = (pb / pb0) / (lin["b"] / lin0["b"])
        Rbc = (pb / pc) / (lin["b"] / lin["c"])
        D2 = k ** 3 * lin["c"] / (2 * np.pi ** 2)
        d = {"G_c": Gc.tolist(), "G_b": Gb.tolist(), "R_bc": Rbc.tolist(), "A_c": (pc / lin["c"]).tolist(),
             "A_b": (pb / lin["b"]).tolist(), "D2_c": D2.tolist(),
             "growth_camb_c": (lin["c"] / lin0["c"]).tolist(), "growth_camb_b": (lin["b"] / lin0["b"]).tolist()}
        # the code's own P(k) (DM deconvolved, gas) against this estimator
        try:
            c = np.loadtxt(find_z(pkf, z))
            kk = c[:, 0]
            m = min(len(kk), len(k))
            d["code_pk_over_this_c"] = (c[:m, 3] / (kk[:m] ** 3 * pc[:m] / (2 * np.pi ** 2))).tolist()
            d["code_pk_over_this_b"] = (c[:m, 5] / (kk[:m] ** 3 * pb[:m] / (2 * np.pi ** 2))).tolist()
        except SystemExit:
            pass
        res["z"][f"{z:g}"] = d
        print(f"\nz = {z:g}: growth since z_init, CAMB: c x{lin['c'][0] / lin0['c'][0]:.1f}, b x{lin['b'][0] / lin0['b'][0]:.1f} (lowest k)")
        print("   k [h/Mpc]    G_c     G_b    R_bc    A_c     A_b    D2_c(CAMB)" +
              ("   code/this c, b" if "code_pk_over_this_c" in d else ""))
        for i in show:
            extra = ""
            if "code_pk_over_this_c" in d and i < len(d["code_pk_over_this_c"]):
                extra = f"   {d['code_pk_over_this_c'][i]:.4f} {d['code_pk_over_this_b'][i]:.4f}"
            print(f"  {k[i]:9.1f}  {Gc[i]:6.3f}  {Gb[i]:6.3f}  {Rbc[i]:6.3f}  {pc[i] / lin['c'][i]:6.3f}  "
                  f"{pb[i] / lin['b'][i]:6.3f}  {D2[i]:9.2e}{extra}")
        lo = D2 < 0.1
        if lo.any():
            print(f"  where D2_c < 0.1 (k <= {k[lo].max():.0f} h/Mpc): G_c {Gc[lo].min():.3f} .. {Gc[lo].max():.3f}, "
                  f"G_b {Gb[lo].min():.3f} .. {Gb[lo].max():.3f}, R_bc {Rbc[lo].min():.3f} .. {Rbc[lo].max():.3f}")
    if a.json:
        with open(a.json, "w") as f:
            json.dump(res, f, indent=1)
    if a.out:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 3, figsize=(13, 4), sharex=True)
        for zs, d in res["z"].items():
            for j, key in enumerate(("G_c", "G_b", "R_bc")):
                ax[j].semilogx(k, d[key], label=f"z = {zs}")
        for j, t in enumerate(("CDM growth / CAMB", "baryon growth / CAMB", "(P_b/P_c) / CAMB")):
            ax[j].set_title(t)
            ax[j].axhline(1, color="k", lw=0.5)
            ax[j].set_xlabel("k [h/Mpc]")
        ax[0].legend()
        fig.suptitle(f"{n}^3, {box} Mpc/h, MUSIC + CAMB ICs at z = {z0:g}, hydro + gravity")
        fig.tight_layout()
        fig.savefig(a.out, dpi=120)
        print("wrote", a.out)


if __name__ == "__main__":
    main()
