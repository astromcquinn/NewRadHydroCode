#!/usr/bin/env python3
"""Constant-opacity transfer function of the moment RT against the analytic forms.

Reads the .L / .N cubes of tools/analysis/moment_transfer.cpp and prints (and
plots) T(k) = (P_N / P_L)^(1/2) of delta = X / <X> - 1 with

  exact      arctan(k lam) / (k lam)                    (Gnedin 2025, eq. 2)
  M1 lin     1 / (1 + k^2 lam^2 / 3)                    (P = N I/3: Wu et al. 2021, eq. 5.14
                                                         without the opacity response)
  OTVET lin  (1 + 2 pi^2 k lam^2 / (3 L_O)) / (1 + k^2 lam^2 / 3), L_O = 7.66 L (their eq. 5.15)
  + GLF      the same two with the scheme's numerical diffusion: per unit c/dx the
             first-order GLF flux adds (1/2) Lap to both moment equations, so that with
             s_d = sin k_d, l = sum_d (2 cos k_d - 2), kap = 1/lam, q = kap - l/2:
               M1     T = kap / (q + |s|^2 / (3 q))
               OTVET  T = (kap + |s|^2 2 pi^2 / (3 q k L_O)) / (q + |s|^2 / (3 q))

  + MUSCL    the same for moment_order = 2 where the limiter is inactive (central
             slopes): the dissipation is a fourth difference, q = kap + sum_d
             (2 - 2 cos k_d)^2 / 8, and the flux difference has s_d = sin k_d
             (3 - cos k_d) / 2 in place of sin k_d

usage: moment_transfer.py N lam prefix [prefix ...] [--png out.png] [--title t]
"""
import sys

import numpy as np

LO_PER_L = 7.66  # integral of d^3x / r^2 over the minimum-image cube, in box lengths


def load(prefix, N):
    L = np.fromfile(prefix + ".L", dtype=np.float32).astype(np.float64).reshape(N, N, N)
    n = np.fromfile(prefix + ".N", dtype=np.float32).astype(np.float64).reshape(N, N, N)
    return L, n


def kgrid(N):
    k1 = 2 * np.pi * np.fft.fftfreq(N)  # rad per cell
    kz, ky, kx = np.meshgrid(k1, k1, k1, indexing="ij")
    return kx, ky, kz


def shells(N, nbin=14):
    kx, ky, kz = kgrid(N)
    k = np.sqrt(kx**2 + ky**2 + kz**2)
    edges = np.geomspace(2 * np.pi / N * 0.99, np.pi * 1.001, nbin + 1)
    idx = np.digitize(k.ravel(), edges) - 1
    ok = (idx >= 0) & (idx < nbin)
    return k, idx, ok, nbin


def binned(q, idx, ok, nbin):
    s = np.bincount(idx[ok], weights=q.ravel()[ok], minlength=nbin)
    c = np.bincount(idx[ok], minlength=nbin)
    return s / np.maximum(c, 1)


def theory(N, lam):
    kx, ky, kz = kgrid(N)
    k = np.sqrt(kx**2 + ky**2 + kz**2)
    k[0, 0, 0] = 1e-30
    kap = 1.0 / lam
    s2 = np.sin(kx) ** 2 + np.sin(ky) ** 2 + np.sin(kz) ** 2
    lap = 2 * (np.cos(kx) + np.cos(ky) + np.cos(kz)) - 6
    q = kap - lap / 2
    LO = LO_PER_L * N
    t = {}
    t["exact"] = np.arctan(k * lam) / (k * lam)
    t["M1 lin"] = 1 / (1 + (k * lam) ** 2 / 3)
    t["M1 lin+GLF"] = kap / (q + s2 / (3 * q))
    t["OTVET lin"] = (1 + 2 * np.pi**2 * k * lam**2 / (3 * LO)) / (1 + (k * lam) ** 2 / 3)
    t["OTVET lin+GLF"] = (kap + s2 * 2 * np.pi**2 / (3 * q * k * LO)) / (q + s2 / (3 * q))
    q2 = kap + ((2 - 2 * np.cos(kx)) ** 2 + (2 - 2 * np.cos(ky)) ** 2 + (2 - 2 * np.cos(kz)) ** 2) / 8
    t2 = sum((np.sin(x) * (3 - np.cos(x)) / 2) ** 2 for x in (kx, ky, kz))
    t["M1 lin+MUSCL"] = kap / (q2 + t2 / (3 * q2))
    t["OTVET lin+MUSCL"] = (kap + t2 * 2 * np.pi**2 / (3 * q2 * k * LO)) / (q2 + t2 / (3 * q2))
    return t


def main():
    a = sys.argv[1:]
    png = title = None
    if "--png" in a:
        i = a.index("--png")
        png = a[i + 1]
        del a[i : i + 2]
    if "--title" in a:
        i = a.index("--title")
        title = a[i + 1]
        del a[i : i + 2]
    N, lam = int(a[0]), float(a[1])
    prefixes = a[2:]
    k, idx, ok, nbin = shells(N)
    kb = binned(k, idx, ok, nbin)
    th = {name: np.sqrt(binned(v**2, idx, ok, nbin)) for name, v in theory(N, lam).items()}
    cols = {}
    for p in prefixes:
        L, n = load(p, N)
        dL = np.fft.fftn(L / L.mean() - 1)
        dN = np.fft.fftn(n / n.mean() - 1)
        PL = binned(np.abs(dL) ** 2, idx, ok, nbin)
        PN = binned(np.abs(dN) ** 2, idx, ok, nbin)
        X = binned((dN * np.conj(dL)).real, idx, ok, nbin)
        name = p.split("/")[-1]
        cols[name] = (np.sqrt(PN / PL), X / np.sqrt(PN * PL))
        print(f"# {name}: <N>/(<L> lam) = {n.mean() / (L.mean() * lam):.4f}")
    names = list(cols)
    print(f"# N = {N}, lambda = {lam} cells (lambda / L = {lam / N:.3f}); T(k), then r(k) with L in brackets")
    print("# k*lam  " + "  ".join(f"{s:>13}" for s in list(th) + names))
    for b in range(nbin):
        row = f"{kb[b] * lam:7.2f}  " + "  ".join(f"{th[s][b]:13.4e}" for s in th)
        row += "  " + "  ".join(f"{cols[s][0][b]:.3e} ({cols[s][1][b]:+.2f})" for s in names)
        print(row)
    if png:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(6.4, 4.8))
        style = {"exact": "k:", "M1 lin": "C1--", "M1 lin+GLF": "C1-.", "OTVET lin": "C0--", "OTVET lin+GLF": "C0-.",
                 "M1 lin+MUSCL": "C1:", "OTVET lin+MUSCL": "C0:"}
        for s, v in th.items():
            ax.loglog(kb * lam, v, style[s], lw=1, label=s)
        for s in names:
            ax.loglog(kb * lam, cols[s][0], "o-", ms=3, lw=1.5, label=s)
        ax.set_xlabel(r"$k\lambda$")
        ax.set_ylabel(r"$T(k) = (P_N/P_L)^{1/2}$")
        ax.set_title(title or f"N = {N}, lambda = {lam:g} cells")
        ax.legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(png, dpi=130)


if __name__ == "__main__":
    main()
