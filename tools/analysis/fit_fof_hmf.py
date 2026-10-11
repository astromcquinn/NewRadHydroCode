#!/usr/bin/env python3
"""Fit the simulated FoF mass function as Sheth-Tormen times a smooth
resolution/normalisation factor, for abundance matching against the
simulated halos (am_fit --sim-hmf, docs/abundance_matching.md).

    fit_fof_hmf.py CATALOG_DIR --box 10 --N 128 [--out data/hmf/fof_128_10mpc.txt]

Model, with N = M / m_particle:
    dn/dlog10M (sim) = dn/dlog10M (ST) * R(N, z),
    R = r_inf / (1 + (N / N_c)^-p),   N_c = N_c0 ((1+z)/7)^q.
Fitted by maximum Poisson likelihood to the counts of all catalogs (halo.z=*,
N >= nmin) in bins of 0.1 dex, including catalogs with few or no halos.
Writes key = value lines read by am_fit --sim-hmf.
"""
import argparse
import os
import subprocess
import sys

import numpy as np

trapezoid = np.trapezoid if hasattr(np, "trapezoid") else np.trapz  # numpy < 2.0

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hmf_compare import catalogs, read_catalog  # noqa: E402

KEYS = ["fof_r_inf", "fof_p", "fof_Nc0", "fof_q"]


def ratio(par, N, z):
    lr, p, lNc0, q = par
    Nc = np.exp(lNc0) * ((1 + z) / 7.0) ** q
    return np.exp(lr) / (1 + (N / Nc) ** (-p))


def nelder_mead(f, x0, step, iters=4000):
    n = len(x0)
    s = [np.array(x0, float)] + [np.array(x0, float) + np.eye(n)[i] * step[i] for i in range(n)]
    fs = [f(x) for x in s]
    for _ in range(iters):
        o = np.argsort(fs)
        s, fs = [s[i] for i in o], [fs[i] for i in o]
        if abs(fs[-1] - fs[0]) < 1e-10:
            break
        c = np.mean(s[:-1], axis=0)
        xr = c + (c - s[-1]); fr = f(xr)
        if fr < fs[0]:
            xe = c + 2 * (c - s[-1]); fe = f(xe)
            s[-1], fs[-1] = (xe, fe) if fe < fr else (xr, fr)
        elif fr < fs[-2]:
            s[-1], fs[-1] = xr, fr
        else:
            xc = c + 0.5 * (s[-1] - c) if fr >= fs[-1] else c + 0.5 * (xr - c)
            fc = f(xc)
            if fc < min(fr, fs[-1]):
                s[-1], fs[-1] = xc, fc
            else:
                s = [s[0]] + [s[0] + 0.5 * (x - s[0]) for x in s[1:]]
                fs = [fs[0]] + [f(x) for x in s[1:]]
    i = int(np.argmin(fs))
    return s[i], fs[i]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("catdir")
    ap.add_argument("--box", type=float, required=True)
    ap.add_argument("--N", type=int, required=True)
    ap.add_argument("--nmin", type=int, default=20)
    ap.add_argument("--hmf-table", default="build/hmf_table")
    ap.add_argument("--out")
    a = ap.parse_args()
    h = 0.7
    vol = (a.box / h) ** 3
    mp = 2.775e11 * 0.27 * a.box ** 3 / a.N ** 3 / h
    lo = np.log10(a.nmin * mp)
    d = 0.1
    data = []  # (z, N_mid array, expected-ST counts per bin, observed counts)
    for zc, path in catalogs(a.catdir):
        cat = read_catalog(path)
        cat = cat[cat["N"] >= a.nmin]
        edges = np.arange(lo, 12.5 + d, d)
        cnt, _ = np.histogram(np.log10(cat["M"].astype(float)), edges)
        r = subprocess.run([a.hmf_table, str(zc), str(lo), str(edges[-1]), "0.02"], check=True,
                           capture_output=True, text=True)
        t = np.loadtxt(r.stdout.splitlines())
        # fine grid inside each bin for the ST integral times R
        data.append((zc, edges, t, cnt))

    def model_counts(par, zc, edges, t):
        lm, st = t[:, 0], t[:, 2]
        R = ratio(par, 10 ** lm / mp, zc)
        y = st * R
        out = []
        for e0, e1 in zip(edges[:-1], edges[1:]):
            s = (lm >= e0 - 1e-9) & (lm <= e1 + 1e-9)
            out.append(trapezoid(y[s], lm[s]) * vol)
        return np.array(out)

    def nll(par):
        if par[1] <= 0:
            return 1e30
        tot = 0.0
        for zc, edges, t, cnt in data:
            mu = np.maximum(model_counts(par, zc, edges, t), 1e-300)
            tot += np.sum(mu - cnt * np.log(mu))
        return tot

    best, f = nelder_mead(nll, [np.log(1.2), 1.5, np.log(60.0), 1.0], [0.2, 0.5, 0.5, 0.5])
    best, f = nelder_mead(nll, best, [0.05, 0.2, 0.2, 0.2])
    lr, p, lNc0, q = best
    vals = [np.exp(lr), p, np.exp(lNc0), q]
    lines = [f"# fit_fof_hmf.py {a.catdir} --box {a.box:g} --N {a.N} (nmin {a.nmin})",
             "# dn/dlog10M(sim) = dn/dlog10M(ST) * fof_r_inf / (1 + (N/N_c)^-fof_p),",
             "# N = M/mpart_msun, N_c = fof_Nc0 ((1+z)/7)^fof_q",
             f"mpart_msun = {mp:.6e}"] + [f"{k} = {v:.6g}" for k, v in zip(KEYS, vals)]
    # Goodness: observed vs model totals per catalog
    lines.append("# z  N_obs  N_model  (N >= nmin)")
    for zc, edges, t, cnt in data:
        lines.append(f"# {zc:.3f} {cnt.sum():5d} {model_counts(best, zc, edges, t).sum():9.1f}")
    txt = "\n".join(lines) + "\n"
    if a.out:
        os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
        open(a.out, "w").write(txt)
    print(txt)


if __name__ == "__main__":
    main()
