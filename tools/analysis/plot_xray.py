#!/usr/bin/env python3
"""The X-ray preheating model of xray_heating = 1 without a simulation
(docs/xray_heating.md): runs build/xray_bkg for several L_X/SFR and plots
the SFR density of the source model, the X-ray heating rate per H atom of
the mean-density neutral IGM, its one-zone temperature and ionised fraction.

usage: plot_xray.py XRAY_BKG AM_PARAMS OUT.png [--lx 39.5,40.5,41.5] [--set key=value ...]
"""
import subprocess
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ordinal single-hue ramp (light -> dark = low -> high L_X/SFR)
RAMP = ["#6da7ec", "#2a78d6", "#104281"]
INK, MUTED, GRID = "#222222", "#6b6b66", "#e4e3de"


def run(exe, am, lx, sets):
    cmd = [exe, am, "--set", f"xray_log_lx_sfr={lx}", "--zmin", "5.5", "--zout-step", "0.1"] + sets
    out = subprocess.run(cmd, check=True, capture_output=True, text=True).stdout
    rows = [l.split() for l in out.splitlines() if l and not l.startswith("#")]
    return np.array(rows, dtype=float)


def main():
    exe, am, png = sys.argv[1:4]
    lxs, sets = [39.5, 40.5, 41.5], []
    a = 4
    while a < len(sys.argv):
        if sys.argv[a] == "--lx":
            lxs = [float(x) for x in sys.argv[a + 1].split(",")]
        elif sys.argv[a] == "--set":
            sets += ["--set", sys.argv[a + 1]]
        a += 2
    data = {lx: run(exe, am, lx, sets) for lx in lxs}
    d0 = data[lxs[0]]
    # columns: z M_min sfrd eps gamma_HI phi_HI Hx_per_H x_e T T_cmb T_ad heat
    fig, ax = plt.subplots(2, 2, figsize=(9, 7), sharex=True)
    for a_ in ax.flat:
        a_.grid(True, color=GRID, lw=0.6)
        a_.tick_params(colors=MUTED, labelcolor=INK)
        for sp in a_.spines.values():
            sp.set_color(GRID)
    ax[0, 0].semilogy(d0[:, 0], d0[:, 2], color=RAMP[1], lw=2)
    ax[0, 0].set_ylabel(r"$\dot\rho_*$ [M$_\odot$ yr$^{-1}$ Mpc$^{-3}$]")
    ax[0, 0].set_title("SFR density, halos with $T_{vir} \\geq 10^4$ K", fontsize=10, color=INK)
    for i, lx in enumerate(lxs):
        d = data[lx]
        lab = rf"$\log L_X/{{\rm SFR}} = {lx:g}$"
        c = RAMP[i % len(RAMP)]
        ax[0, 1].semilogy(d[:, 0], d[:, 6], color=c, lw=2, label=lab)
        ax[1, 0].semilogy(d[:, 0], d[:, 8], color=c, lw=2, label=lab)
        ax[1, 1].semilogy(d[:, 0], d[:, 7], color=c, lw=2, label=lab)
    ax[1, 0].semilogy(d0[:, 0], d0[:, 9], color=MUTED, lw=1.2, ls="--", label=r"$T_{\rm CMB}$")
    ax[1, 0].semilogy(d0[:, 0], d0[:, 10], color=MUTED, lw=1.2, ls=":", label="adiabatic (code's floor)")
    ax[0, 1].set_ylabel("X-ray heat per H atom [erg/s]")
    ax[0, 1].set_title("heating of the mean-density neutral IGM", fontsize=10, color=INK)
    ax[1, 0].set_ylabel("T [K]")
    ax[1, 0].set_title("one-zone temperature", fontsize=10, color=INK)
    ax[1, 1].set_ylabel(r"$x_e$")
    ax[1, 1].set_title("one-zone ionised fraction", fontsize=10, color=INK)
    ax[0, 1].set_ylim(1e-35, None)
    ax[0, 0].set_ylim(1e-7, 1e-1)
    for a_ in ax[1]:
        a_.set_xlabel("z")
    ax[0, 0].set_xlim(d0[:, 0].max() if d0[:, 0].max() < 25 else 25, d0[:, 0].min())
    ax[1, 0].legend(fontsize=8, frameon=False, labelcolor=INK)
    ax[0, 1].legend(fontsize=8, frameon=False, labelcolor=INK)
    fig.tight_layout()
    fig.savefig(png, dpi=130)
    print("wrote", png)


if __name__ == "__main__":
    main()
