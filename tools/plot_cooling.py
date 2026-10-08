#!/usr/bin/env python3
"""Cooling-function components for legacy_physics = 1 and 0.

usage: plot_cooling.py <rates.csv> <out.png> [--nH 1.0]

<rates.csv> comes from `build/tests/test_rates --csv rates.csv` (the code's
own rate tables, both modes). The figure has two panels:

 (a) the collisional-ionisation-equilibrium (CIE) cooling function
     Lambda/nH^2 of H+He gas, split into recombination (slots 4-6),
     collisional ionisation (7-9), collisional excitation (10-12),
     HeII dielectronic recombination (15) and bremsstrahlung (13), plus the
     total. Inverse Compton (14) depends on redshift and is left out.
     The ionisation state uses the code's own rates: case B recombination
     (alpha_rec(4:6), which include the dielectronic term in the new mode)
     against collisional ionisation (gamma_col), no radiation. The HeI
     (2^3S) excitation term (slot 11) scales as ne^2 nHeII, so it depends on
     density: it is evaluated at --nH (default 1 cm^-3).
 (b) the ratio new/legacy of each coefficient changed by P1, P4-P8.
 (c) the table-lookup error (P9) for a few steep rates: the worst case of
     the legacy nearest-node lookup (half an interval off, |ln(R_i+1/R_i)|/2)
     and of the log-log interpolation at mid-interval (|second difference
     of ln R|/8), both from the tables themselves; test_rates checks the
     latter against direct formula evaluation.

It also prints the components at a few temperatures.
"""
import argparse

import numpy as np

XHY, YHE = 0.76, 0.24  # (double)0.76f etc. are irrelevant at plot precision
Y = YHE / (4 * XHY)    # nHe/nH, mHe = 4 mH as in the code

# Categorical slots 1-5 (reference palette, fixed order) and text inks
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"


def load(fn):
    a = np.genfromtxt(fn, delimiter=",", names=True)
    return a


def cie(a, m, nH):
    """CIE abundances and cooling components (per nH^2) for mode m ('L'/'N')."""
    g = lambda k: a[f"{m}_{k}"]  # noqa: E731
    with np.errstate(divide="ignore", invalid="ignore"):
        rH = np.where(g("alpha4") > 0, g("gcol1") / g("alpha4"), 0)       # nHII/nHI
        r1 = np.where(g("alpha5") > 0, g("gcol2") / g("alpha5"), 0)       # nHeII/nHeI
        r2 = np.where(g("alpha6") > 0, g("gcol3") / g("alpha6"), 0)       # nHeIII/nHeII
    xHI = 1 / (1 + rH)
    xHII = 1 - xHI
    fHeI = 1 / (1 + r1 + r1 * r2)
    fHeII = r1 * fHeI
    fHeIII = r2 * fHeII
    nHI, nHII = xHI * nH, xHII * nH
    nHeI, nHeII, nHeIII = (f * Y * nH for f in (fHeI, fHeII, fHeIII))
    ne = nHII + nHeII + 2 * nHeIII
    nff = nHII + nHeII + 4 * nHeIII
    c = lambda n: g(f"cool{n}")  # noqa: E731
    comp = {
        "recombination": c(4) * ne * nHII + c(5) * ne * nHeII + c(6) * ne * nHeIII,
        "collisional ionisation": c(7) * ne * nHI + c(8) * ne * nHeI + c(9) * ne * nHeII,
        "collisional excitation": c(10) * ne * nHI
        + (c(11) * ne * ne * nHeII if m == "N" else c(11) * ne * nHeI)
        + c(12) * ne * nHeII,
        "dielectronic (HeII)": c(15) * ne * nHeII,
        "bremsstrahlung": c(13) * ne * nff,
    }
    comp = {k: v / nH**2 for k, v in comp.items()}
    comp["total"] = sum(comp.values())
    return comp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("out")
    ap.add_argument("--nH", type=float, default=1.0)
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    a = load(args.csv)
    lgT = a["logT"]
    T = 10**lgT
    sel = (lgT >= 4.0) & (lgT <= 8.0)
    L, N = cie(a, "L", args.nH), cie(a, "N", args.nH)

    plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                         "xtick.color": INK2, "ytick.color": INK2, "axes.linewidth": 0.8})
    fig, (ax, bx, cx) = plt.subplots(1, 3, figsize=(16, 4.6), constrained_layout=True)

    names = ["recombination", "collisional ionisation", "collisional excitation",
             "dielectronic (HeII)", "bremsstrahlung"]
    for k, col in zip(names, C):
        ax.plot(T[sel], N[k][sel], color=col, lw=2, label=k)
        if np.any(L[k][sel] > 0):
            ax.plot(T[sel], L[k][sel], color=col, lw=1.2, ls=(0, (4, 3)))
    ax.plot(T[sel], N["total"][sel], color=INK, lw=2, label="total")
    ax.plot(T[sel], L["total"][sel], color=INK, lw=1.2, ls=(0, (4, 3)))
    ax.plot([], [], color=INK2, lw=2, label="legacy_physics = 0 (solid)")
    ax.plot([], [], color=INK2, lw=1.2, ls=(0, (4, 3)), label="legacy_physics = 1 (dashed)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(1e-25, 3e-21)
    ax.set_xlim(1e4, 1e8)
    ax.set_xlabel("T [K]")
    ax.set_ylabel(r"$\Lambda / n_H^2$  [erg cm$^3$ s$^{-1}$]")
    ax.set_title(f"(a) CIE cooling function components, H+He, $n_H$ = {args.nH:g} cm$^{{-3}}$",
                 loc="left", color=INK, fontsize=9.5)
    ax.grid(True, which="major", color=GRID, lw=0.6)
    ax.legend(frameon=False, fontsize=8, loc="upper right")

    ratios = [
        ("HII rec. cooling B (P1)", a["N_cool4"] / a["L_cool4"]),
        ("HeII rec. cooling B (P1)", a["N_cool5"] / a["L_cool5"]),
        ("HeIII rec. cooling B (P1)", a["N_cool6"] / a["L_cool6"]),
        ("HI coll. ion. cooling (P5)", a["N_cool7"] / np.where(a["L_cool7"] > 0, a["L_cool7"], np.nan)),
        ("HeIII rec. rate B (P4)", a["N_alpha6"] / a["L_alpha6"]),
    ]
    with np.errstate(divide="ignore", invalid="ignore"):
        ratios.append(("HeII rec. rate B, +dielectronic (P6)", a["N_alpha5"] / a["L_alpha5"]))
    ratios.append(("HeII rec. cooling A (P7)", a["N_cool2"] / a["L_cool2"]))
    ratios.append(("bremsstrahlung, g_ff(T) (P8)", a["N_cool13"] / a["L_cool13"]))
    cols = C + [INK2, C[0], C[1]]
    styles = ["-"] * 6 + [(0, (1, 1.5)), (0, (1, 1.5))]
    for (lab, r), col, ls in zip(ratios, cols, styles):
        bx.plot(T[sel], r[sel], color=col, lw=2, ls=ls, label=lab)
    bx.axhline(1, color=INK2, lw=0.8)
    bx.set_xscale("log")
    bx.set_yscale("log")
    bx.set_xlim(1e4, 1e8)
    bx.set_xlabel("T [K]")
    bx.set_ylabel("new / legacy")
    bx.set_title("(b) coefficients changed by legacy_physics = 0", loc="left", color=INK,
                 fontsize=9.5)
    bx.grid(True, which="major", color=GRID, lw=0.6)
    bx.legend(frameon=False, fontsize=7.5, loc="lower left")

    # (c) lookup error, from the legacy_physics = 0 tables
    lk = [("HI coll. excitation (10)", "N_cool10"), ("HeII coll. ionisation (9)", "N_cool9"),
          ("HII rec. rate B (alpha 4)", "N_alpha4")]
    Tm = 10 ** (0.5 * (lgT[1:] + lgT[:-1]))
    selm = (Tm >= 1e4) & (Tm <= 1e8)
    for (lab, key), col in zip(lk, C):
        with np.errstate(divide="ignore", invalid="ignore"):
            lr = np.log(a[key])  # -inf below the code's cut-offs
            near = np.abs(np.diff(lr)) / 2
            d2 = np.full(lr.size, np.nan)
            d2[1:-1] = lr[:-2] - 2 * lr[1:-1] + lr[2:]
            interp = np.abs(0.5 * (d2[:-1] + d2[1:])) / 8
        cx.plot(Tm[selm], near[selm], color=col, lw=1.2, ls=(0, (4, 3)))
        cx.plot(Tm[selm], interp[selm], color=col, lw=2, label=lab)
    cx.plot([], [], color=INK2, lw=2, label="log-log interpolation (solid)")
    cx.plot([], [], color=INK2, lw=1.2, ls=(0, (4, 3)), label="nearest node, legacy (dashed)")
    cx.set_xscale("log")
    cx.set_yscale("log")
    cx.set_xlim(1e4, 1e8)
    cx.set_ylim(1e-8, 1)
    cx.set_xlabel("T [K]")
    cx.set_ylabel("max relative error of the rate lookup")
    cx.set_title("(c) table lookup (P9), $\\Delta\\log T$ = 0.01", loc="left", color=INK,
                 fontsize=9.5)
    cx.grid(True, which="major", color=GRID, lw=0.6)
    cx.legend(frameon=False, fontsize=7.5, loc="lower left")
    fig.savefig(args.out, dpi=150)

    print(f"CIE cooling components, Lambda/nH^2 [erg cm^3/s], nH = {args.nH:g}")
    print(f"{'T':>8} {'mode':>4} " + " ".join(f"{k[:14]:>14}" for k in names + ["total"]))
    for lt in (4.2, 4.5, 5.0, 5.5, 6.0, 7.0):
        i = int(np.argmin(np.abs(lgT - lt)))
        for m, d in (("L", L), ("N", N)):
            print(f"{T[i]:8.3g} {m:>4} " + " ".join(f"{d[k][i]:14.3e}" for k in names + ["total"]))


if __name__ == "__main__":
    main()
