#!/usr/bin/env python3
"""Over-cooled and over-heated cells of the unheated IGM in the unsplit hydro
(docs/hydro_unsplit.md section 7, U7; docs/thermal_history.md section 4).

usage: cold_cells.py OUT.png RUN:N:BOX:REF:LABEL [RUN:N:BOX:REF:LABEL ...]

RUN is an apps/radhydro output directory of a run without heating with
u_gas and e_gas in its step_* fields, N the gas grid, BOX in Mpc/h, REF the
reference temperature of mean-density gas: floor (the code's temp_min(z),
thermal_init = floor) or recfast (thermal_init = recfast). A cell's
adiabat is T_ref(z) Delta^(2/3). Prints per run and step the fraction of
cells colder than 0.1 x their adiabat (over-cooled), the percentiles of
T/adiabat over the cells with Delta < 3 and the median T/T_ref of the
mean-density cells (0.9 < Delta < 1.1); plots the first and the
percentiles against z.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import thermal_history as th  # noqa: E402  (also sets the Agg backend)
import matplotlib.pyplot as plt  # noqa: E402

PCT = (1, 10, 50, 90, 99)


def snapshots(run):
    for d in sorted(os.listdir(run)):
        p = os.path.join(run, d)
        if not (d.startswith("step_") and os.path.exists(os.path.join(p, "e_gas.bin"))):
            continue
        st = dict(l.split(" = ", 1) for l in open(os.path.join(p, "stats.txt")) if " = " in l)
        u = np.fromfile(os.path.join(p, "u_gas.bin"), dtype=np.float32).reshape(-1, 5).astype(np.float64)
        e = np.fromfile(os.path.join(p, "e_gas.bin"), dtype=np.float32).astype(np.float64)
        yield float(st["redshift"]), u[:, 0], e


def history(run, N, box, ref):
    Tref = th.floor if ref == "floor" else th.reference()
    rows = []
    for z, rho, e in snapshots(run):
        T = e / rho * th.temp1_unit(1 / (1 + z), box, N)
        D = rho * th.OM / th.OB
        r = T / (Tref(z) * D ** (2 / 3))
        mean = (D > 0.9) & (D < 1.1)
        rows.append([z, np.mean(r < 0.1), *np.percentile(r[D < 3], PCT), np.median(T[mean]) / Tref(z)])
    return np.array(rows)


def main():
    out = sys.argv[1]
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
    colors = ["#6da7ec", "#104281", "#eb6834", "#a33d12", "#1baf7a", "#0d6e4a"]
    for i, spec in enumerate(sys.argv[2:]):
        run, N, box, ref, label = spec.split(":", 4)
        h = history(run, int(N), float(box), ref)
        c = colors[i % len(colors)]
        print(f"# {label} ({run})\n#      z  f_cold  " + "  ".join(f"p{p:<4d}" for p in PCT) + "  Tmean/Tref")
        for row in h:
            print(f"{row[0]:8.2f} {row[1]:7.4f} " + " ".join(f"{x:6.3f}" for x in row[2:7]) + f" {row[7]:8.3f}")
        ls = "--" if "self" not in label else "-"
        ax[0].plot(h[:, 0], 100 * h[:, 1], ls, color=c, lw=1.5, label=label)
        ax[1].fill_between(h[:, 0], h[:, 3], h[:, 5], color=c, alpha=0.15, lw=0)
        ax[1].plot(h[:, 0], h[:, 4], ls, color=c, lw=1.5, label=label)
        ax[1].plot(h[:, 0], h[:, 2], ls, color=c, lw=0.8)
        ax[1].plot(h[:, 0], h[:, 6], ls, color=c, lw=0.8)
    for a in ax:
        a.grid(True, color="#e4e3de", lw=0.6)
        a.set_xlabel("z")
        a.set_xlim(60, 5.5)
    ax[0].set_ylabel("cells below 0.1 x their adiabat [%]")
    ax[0].legend(fontsize=8, frameon=False)
    ax[1].set_yscale("log")
    ax[1].set_ylim(1e-3, 1e3)
    ax[1].set_ylabel("T / adiabat, cells with Delta < 3")
    ax[1].set_title("median (thick), 10-90 % (band), 1 % and 99 % (thin)", fontsize=9, color="#6b6b66")
    fig.tight_layout()
    fig.savefig(out, dpi=130)
    print("wrote", out)


if __name__ == "__main__":
    main()
