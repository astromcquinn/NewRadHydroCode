#!/usr/bin/env python3
"""The thermal history of the gas before any heating (docs/thermal_history.md):
the simulated temperature of mean-density gas and of voids against the
recombination-era reference (CAMB + RECFAST for the simulation cosmology,
data/thermal/recfast_sim_cosmology.txt), the code's temperature floor
temp_min(z) and T_CMB.

usage: thermal_history.py OUT.png RUN:N:BOX:LABEL [RUN:N:BOX:LABEL ...]

RUN is an apps/radhydro output directory with u_gas (and, for hydro_scheme =
unsplit, e_gas) in its step_* fields; N the gas grid, BOX in Mpc/h. The
temperature is that of neutral gas, (gamma-1) e/rho / k * mu mH with
mu = 1/(XHy + YHe/4) (calc_sim_var's temp1_unit). Prints a table per run.
"""
import glob
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REF = os.path.join(HERE, "..", "..", "data", "thermal", "recfast_sim_cosmology.txt")

f32 = lambda x: float(np.float32(x))
OM, OB, H0H, XHY, YHE, GAM = f32(0.27), f32(0.045), f32(0.70), f32(0.76), f32(0.24), f32(5 / 3)
H0, KB, MH, MPC, TCMB0 = 3.24086e-18, 1.38066e-16, 1.67223e-24, 3.08560e24, 2.725
ZDEC1 = 1e3 * (OB * H0H**2) ** 0.4  # calc_sim_var: 1 + z of the floor's decoupling


def temp1_unit(a, box, N):
    x = a * box / H0H / N * MPC
    t = 2 / (3 * H0 * H0H * np.sqrt(OM)) * a * a
    return (x / t) ** 2 * (GAM - 1) / KB * MH / (XHY + YHE / 4)


def floor(z):
    z = np.asarray(z, dtype=float)
    return np.where(1 + z > ZDEC1, TCMB0 * (1 + z), TCMB0 * (1 + z) ** 2 / ZDEC1)


def reference():
    d = np.loadtxt(REF)
    z, T = d[::-1, 0], d[::-1, 2]
    return lambda zz: np.exp(np.interp(np.log1p(zz), np.log1p(z), np.log(T)))


def run_history(run, N, box):
    rows = []
    for d in sorted(glob.glob(os.path.join(run, "step_*"))):
        if not os.path.exists(os.path.join(d, "u_gas.bin")):
            continue
        st = dict(l.split(" = ", 1) for l in open(os.path.join(d, "stats.txt")) if " = " in l)
        z = float(st["redshift"])
        u = np.fromfile(os.path.join(d, "u_gas.bin"), dtype=np.float32).reshape(-1, 5).astype(np.float64)
        rho = u[:, 0]
        pe = os.path.join(d, "e_gas.bin")
        if os.path.exists(pe):
            eth = np.fromfile(pe, dtype=np.float32).astype(np.float64)
        else:
            eth = u[:, 4] - (u[:, 1] ** 2 + u[:, 2] ** 2 + u[:, 3] ** 2) / rho / 2
        T = eth / rho * temp1_unit(1 / (1 + z), box, N)
        D = rho * OM / OB
        mean = (D > 0.9) & (D < 1.1)
        void = (D > 0.3) & (D < 0.4)
        rows.append((z, np.median(T[mean]), np.median(T[void]) if void.sum() > 10 else np.nan,
                     np.mean(T <= 1.01)))
    return np.array(rows)


def main():
    out = sys.argv[1]
    Tref = reference()
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
    zz = np.linspace(5.5, 100, 400)
    for a_ in ax:
        a_.grid(True, color="#e4e3de", lw=0.6)
    ax[0].semilogy(zz, Tref(zz), color="#222222", lw=2, label="RECFAST (no heating)")
    ax[0].semilogy(zz, floor(zz), color="#6b6b66", lw=1.2, ls=":", label="code floor temp_min(z)")
    ax[0].semilogy(zz, TCMB0 * (1 + zz), color="#6b6b66", lw=1.2, ls="--", label="T_CMB")
    ax[1].axhline(1, color="#222222", lw=1.5)
    ax[1].plot(zz, floor(zz) / Tref(zz), color="#6b6b66", lw=1.2, ls=":")
    colors = ["#6da7ec", "#104281", "#eb6834", "#a33d12", "#1baf7a", "#0d6e4a"]
    for i, spec in enumerate(sys.argv[2:]):
        run, N, box, label = spec.split(":", 3)
        h = run_history(run, int(N), float(box))
        if h.size == 0:
            print(f"# {label}: no step_* with u_gas in {run}, skipped")
            continue
        c = colors[i % len(colors)]
        print(f"# {label} ({run})\n#      z  T_mean[K]  T/RECFAST  T_void[K]  frac(T<=1.01K)")
        for z, Tm, Tv, f1 in h:
            print(f"{z:8.2f} {Tm:10.3f} {Tm / Tref(z):9.3f} {Tv:10.3f} {f1:8.3f}")
        ax[0].semilogy(h[:, 0], h[:, 1], "o-", color=c, ms=4, lw=1.5, label=label)
        ax[1].plot(h[:, 0], h[:, 1] / Tref(h[:, 0]), "o-", color=c, ms=4, lw=1.5, label=label)
    ax[0].set_xlabel("z")
    ax[0].set_ylabel("T of mean-density gas [K]")
    ax[0].set_xlim(100, 5.5)
    ax[0].set_ylim(0.3, 400)
    ax[0].legend(fontsize=8, frameon=False)
    ax[1].set_xlabel("z")
    ax[1].set_ylabel("T / T_RECFAST")
    ax[1].set_xlim(100, 5.5)
    ax[1].set_yscale("log")
    ax[1].set_ylim(0.4, 10)
    ax[1].set_title("dotted: the code's floor", fontsize=9, color="#6b6b66")
    fig.tight_layout()
    fig.savefig(out, dpi=130)
    print("wrote", out)


if __name__ == "__main__":
    main()
