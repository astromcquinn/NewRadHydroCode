#!/usr/bin/env python3
"""Check the spike hydro kernel against the exact Sod solution.

    hydro_serial --n 64 --steps S --ic sod --out sod64.bin   (prints t)
    python3 spike/bench/check_sod.py sod64.bin 64 <t>

--ic sod puts Sod tubes (rho,p = 1,1 | 0.125,0.1, gamma = 5/3) at x = 0.25
(left state at x < 0.25) and its mirror image at x = 0.75 in a periodic box.
Prints the L1 density error along x (averaged over y, z) against the exact
Riemann solution (Toro 2009, chapter 4).
"""
import sys

import numpy as np

G = 5.0 / 3.0


def exact_sod(x, t, x0, rl=1.0, pl=1.0, rr=0.125, pr=0.1, ul=0.0, ur=0.0):
    cl, cr = np.sqrt(G * pl / rl), np.sqrt(G * pr / rr)

    def f(p, rk, pk, ck):
        if p > pk:
            A, B = 2 / ((G + 1) * rk), (G - 1) / (G + 1) * pk
            return (p - pk) * np.sqrt(A / (p + B))
        return 2 * ck / (G - 1) * ((p / pk) ** ((G - 1) / (2 * G)) - 1)

    lo, hi = 1e-8, 10.0
    for _ in range(200):
        pm = 0.5 * (lo + hi)
        if f(pm, rl, pl, cl) + f(pm, rr, pr, cr) + (ur - ul) > 0:
            hi = pm
        else:
            lo = pm
    um = 0.5 * (ul + ur) + 0.5 * (f(pm, rr, pr, cr) - f(pm, rl, pl, cl))
    rho = np.empty_like(x)
    for n, xx in enumerate(x):
        s = (xx - x0) / t
        if s < um:  # left of contact
            if pm <= pl:  # rarefaction
                shl, cml = ul - cl, cl * (pm / pl) ** ((G - 1) / (2 * G))
                stl = um - cml
                if s < shl:
                    rho[n] = rl
                elif s > stl:
                    rho[n] = rl * (pm / pl) ** (1 / G)
                else:
                    rho[n] = rl * (2 / (G + 1) + (G - 1) / ((G + 1) * cl) * (ul - s)) ** (2 / (G - 1))
            else:
                rho[n] = rl  # (shock to the left: not needed for Sod)
        else:  # right of contact: shock
            sr = ur + cr * np.sqrt((G + 1) / (2 * G) * pm / pr + (G - 1) / (2 * G))
            if s > sr:
                rho[n] = rr
            else:
                rho[n] = rr * ((pm / pr + (G - 1) / (G + 1)) / ((G - 1) / (G + 1) * pm / pr + 1))
    return rho


def main():
    path, n, t = sys.argv[1], int(sys.argv[2]), float(sys.argv[3])
    U = np.fromfile(path, dtype="<f8").reshape((n, n, n, 5))  # (k, j, i, comp)
    rho = U[..., 0].mean(axis=(0, 1))
    x = (np.arange(n) + 0.5) / n
    ex = np.where(x < 0.5, exact_sod(x, t, 0.25), exact_sod(1 - x, t, 0.25))
    l1 = np.mean(np.abs(rho - ex))
    spread = np.abs(U[..., 0] - rho[None, None, :]).max()
    print(f"n = {n}  t = {t:.6g}  L1(rho) = {l1:.4e}  max |rho - <rho>_yz| = {spread:.2e}")


if __name__ == "__main__":
    main()
