#!/usr/bin/env python3
"""Fit the Sheth-Tormen coefficients (A, a, p) to on-the-fly FoF catalogs
(docs/abundance_matching.md: the production abundance matching uses this
fitted ST, so the L_UV(M) relation does not depend on the box or the
resolution).

usage: fit_st_fof.py --run CATDIR:BOX:NDM [--run ...] [--nmin 300]
                     [--zmin 5.4] [--zmax 12] [--dlogm 0.1]
                     [--hmf-table build/hmf_table] [--out data/hmf/st_fof_X.txt]
                     [--check CATDIR:BOX:NDM:NMIN ...] [--fix-p P]

f(nu) = A sqrt(2a/pi) (1 + (a nu^2)^-p) nu exp(-a nu^2/2), with
nu = 1.686/sigma(M,z) in the simulation cosmology (apps/hmf_table). The fit
maximises the Poisson likelihood of the halo counts in dlogm bins above
--nmin particles, over every catalog with zmin <= z <= zmax. Consecutive
catalogs share most of their halos, so the likelihood is not that of
independent data: the best fit is meaningful, its formal errors are not
(not printed). --check runs (other boxes / resolutions) are only compared
with the fit, each above its own NMIN particles.
"""
import argparse
import os
import subprocess
import sys

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hmf_compare import catalogs, read_catalog  # noqa: E402

H = 0.7


def mpart(box, ndm):  # Msun, cosmo_tools.f90 mp_unit (as hmf_compare.py)
    return 2.775e11 * 0.27 * box ** 3 / ndm ** 3 / H


def table(exe, z, lo, hi, d, st):
    r = subprocess.run([exe, str(z), str(lo), str(hi), str(d), "--st-params",
                        ",".join(str(x) for x in st)], check=True, capture_output=True, text=True)
    a = np.loadtxt([l for l in r.stdout.splitlines() if not l.startswith("#")])
    return a[:, 0], a[:, 2]  # logM, dn/dlog10M (ST)


def nu_base(exe, z, lo, hi, d):
    """nu(logM) and the prefactor b with dn/dlog10M = b f(nu), from two ST
    tables with p = 0: f(A=1,a=1) = 2 sqrt(2/pi) nu e^{-nu^2/2},
    f(A=1,a=1/2) = 2 sqrt(1/pi) nu e^{-nu^2/4}."""
    lm, d1 = table(exe, z, lo, hi, d, (1, 1, 0))
    _, d2 = table(exe, z, lo, hi, d, (1, 0.5, 0))
    nu = np.sqrt(4 * np.log(np.sqrt(2) * d2 / d1))
    f1 = 2 * np.sqrt(2 / np.pi) * nu * np.exp(-0.5 * nu * nu)
    return lm, nu, d1 / f1


def f_st(nu, A, a, p):
    an2 = a * nu * nu
    return A * np.sqrt(2 * a / np.pi) * (1 + an2 ** (-p)) * nu * np.exp(-0.5 * an2)


class Cat:
    def __init__(self, exe, z, path, box, ndm, nmin, dlogm, fine=0.01):
        h = read_catalog(path)
        m = h["M"].astype(float)[h["N"] >= nmin]
        self.z, self.vol = z, (box / H) ** 3
        lo = np.log10(nmin * mpart(box, ndm))
        hi = max(np.log10(m.max()) if m.size else lo, lo) + dlogm
        self.edges = np.arange(lo, hi + 1e-9, dlogm)
        self.n = np.histogram(np.log10(m), self.edges)[0]
        self.lm, self.nu, self.base = nu_base(exe, z, lo, self.edges[-1] + fine, fine)
        # trapezoid weights of each fine point in each bin
        nb = self.edges.size - 1
        self.W = np.zeros((nb, self.lm.size))
        for b in range(nb):
            s = (self.lm >= self.edges[b] - 1e-9) & (self.lm <= self.edges[b + 1] + 1e-9)
            idx = np.nonzero(s)[0]
            w = np.full(idx.size, fine)
            w[0] = w[-1] = fine / 2
            self.W[b, idx] = w

    def mu(self, A, a, p):
        return self.vol * (self.W @ (self.base * f_st(self.nu, A, a, p)))


def load(spec, exe, nmin, dlogm, zmin, zmax):
    parts = spec.split(":")
    d, box, ndm = parts[0], float(parts[1]), int(parts[2])
    if len(parts) > 3:
        nmin = int(parts[3])
    out = []
    for z, p in catalogs(d):
        if zmin <= z <= zmax:
            out.append(Cat(exe, z, p, box, ndm, nmin, dlogm))
    return out


def report(cats, A, a, p, label, out):
    out(f"# {label}: data / fitted ST per bin (counts), selected redshifts")
    zs = sorted({round(c.z, 1) for c in cats})
    pick = [min(cats, key=lambda c: abs(c.z - zt)) for zt in (5.5, 6, 7, 8, 10, 12)
            if min(abs(c.z - zt) for c in cats) < 0.3]
    for c in pick:
        mu = c.mu(A, a, p)
        cells = []
        for b in range(c.n.size):
            if c.n[b] >= 10:
                cells.append(f"{0.5 * (c.edges[b] + c.edges[b + 1]):.2f}:{c.n[b] / mu[b]:.2f}({c.n[b]})")
        out(f"  z={c.z:.2f} total {c.n.sum()} / {mu.sum():.0f} = {c.n.sum() / mu.sum():.3f} | "
            + " ".join(cells))
    _ = zs


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", action="append", required=True)
    ap.add_argument("--check", action="append", default=[])
    ap.add_argument("--nmin", type=int, default=300)
    ap.add_argument("--zmin", type=float, default=5.4)
    ap.add_argument("--zmax", type=float, default=12.0)
    ap.add_argument("--dlogm", type=float, default=0.1)
    ap.add_argument("--hmf-table", default="build/hmf_table")
    ap.add_argument("--fix-p", type=float, default=None)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    cats = []
    for spec in a.run:
        cats += load(spec, a.hmf_table, a.nmin, a.dlogm, a.zmin, a.zmax)
    lines = []

    def out(s):
        print(s)
        lines.append(s)

    def nll(x):
        A, aa, p = np.exp(x[0]), np.exp(x[1]), (a.fix_p if a.fix_p is not None else x[2])
        s = 0.0
        for c in cats:
            mu = np.maximum(c.mu(A, aa, p), 1e-300)
            s += np.sum(mu - c.n * np.log(mu))
        return s

    x0 = [np.log(0.3222), np.log(0.707), 0.3]
    r = minimize(nll, x0, method="Nelder-Mead", options=dict(xatol=1e-6, fatol=1e-6, maxiter=4000))
    A, aa = np.exp(r.x[0]), np.exp(r.x[1])
    p = a.fix_p if a.fix_p is not None else r.x[2]
    ntot = sum(c.n.sum() for c in cats)
    out("# Sheth-Tormen fitted to FoF halo counts (tools/analysis/fit_st_fof.py)")
    out("# runs: " + " ".join(a.run))
    out(f"# {len(cats)} catalogs, {a.zmin} <= z <= {a.zmax}, N >= {a.nmin} particles, "
        f"{ntot} halos (summed over catalogs), dlogM {a.dlogm}")
    out(f"# -lnL: fit {r.fun:.1f}, standard ST {nll(x0):.1f}")
    out(f"st_A = {A:.6g}")
    out(f"st_a = {aa:.6g}")
    out(f"st_p = {p:.6g}")
    out(f"# am_fit --hmf st --st-params {A:.6g},{aa:.6g},{p:.6g}")
    report(cats, A, aa, p, "fit runs", out)
    for spec in a.check:
        cc = load(spec, a.hmf_table, a.nmin, a.dlogm, a.zmin, a.zmax)
        report(cc, A, aa, p, "check " + spec, out)
    if a.out:
        with open(a.out, "w") as f:
            f.write("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
