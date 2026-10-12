#!/usr/bin/env python3
"""Figures of docs/rt_multicell.md from the JSON outputs of
rt_multicell_ifront.py and rt_multicell_fields.py.

usage: rt_multicell_figures.py ifront --json A.json B.json C.json
                               --labels "share" "share+pc" "exposure+pc" --out PNG
       rt_multicell_figures.py profiles --json A.json [B.json] --labels ... --c 1
                               --out PNG
       rt_multicell_figures.py fields --json f1.json f2.json ... --labels ... --out PNG

ifront: R - R(k=1) against t for every k, one column per JSON (scheme), one
row per c_reduced. profiles: the shell-mean x_HI and T at the profile steps
of each JSON for k = 1, 2, 4, 8 (one row pair per JSON). fields: the
acceptance metrics of rt_multicell_fields.py against k.
"""
import argparse
import json

import numpy as np

COLS = {1: "#1f77b4", 2: "#ff7f0e", 4: "#2ca02c", 8: "#d62728", 16: "#9467bd"}


def load(fn):
    return json.load(open(fn))


def fig_ifront(a):
    import matplotlib.pyplot as plt
    ds = [load(f) for f in a.json]
    cs = sorted({r["c"] for r in ds[0]["runs"]})
    fig, axs = plt.subplots(len(cs), len(ds), figsize=(4.2 * len(ds), 3.2 * len(cs)),
                            sharex=True, squeeze=False)
    for j, (d, lab) in enumerate(zip(ds, a.labels)):
        for i, c in enumerate(cs):
            ax = axs[i, j]
            r1 = [r for r in d["runs"] if r["c"] == c and r["k"] == 1][0]
            for r in d["runs"]:
                if r["c"] != c:
                    continue
                dR = np.array(r["R"]) - np.array(r1["R"])
                ax.plot(r["t_myr"], dR, "o-", ms=3, color=COLS.get(r["k"]), label=f"k = {r['k']}")
            ax.axhspan(-0.5, 0.5, color="0.92")
            ax.axhline(0, color="k", lw=0.5)
            ax.set_ylim(-3, 1.5)
            ax.set_title(f"{lab}, c_reduced = {c:g}", fontsize=9)
            if j == 0:
                ax.set_ylabel("R - R(k=1) [cells]")
            if i == len(cs) - 1:
                ax.set_xlabel("t [Myr]")
            if i == 0 and j == 0:
                ax.legend(fontsize=8)
    fig.suptitle(f"uniform I-front, N = {ds[0]['N']}: front radius against k = 1 "
                 "(grey: half a cell)", fontsize=10)
    fig.tight_layout()
    fig.savefig(a.out, dpi=110)


def fig_profiles(a):
    import matplotlib.pyplot as plt
    ds = [load(f) for f in a.json]
    steps = sorted({int(s) for r in ds[0]["runs"] for s in r["profiles"]})
    fig, axs = plt.subplots(2 * len(ds), len(steps), figsize=(3.6 * len(steps), 2.6 * 2 * len(ds)),
                            sharex=True, squeeze=False)
    for i, (d, lab) in enumerate(zip(ds, a.labels)):
        for j, st in enumerate(steps):
            ax, axT = axs[2 * i, j], axs[2 * i + 1, j]
            for r in d["runs"]:
                if r["c"] != a.c or str(st) not in r["profiles"]:
                    continue
                p = r["profiles"][str(st)]
                rr = np.arange(len(p["xHI"])) + 0.5
                ax.semilogy(rr, np.maximum(p["xHI"], 1e-7), "-", color=COLS.get(r["k"]),
                            label=f"k = {r['k']}")
                axT.plot(rr, np.array(p["T"]) / 1e4, "-", color=COLS.get(r["k"]))
                t = p["t"] / 3.156e13
            ax.set_ylim(1e-6, 2)
            ax.set_title(f"{lab}: step {st}, t = {t:.1f} Myr", fontsize=8)
            if j == 0:
                ax.set_ylabel("x_HI (shell mean)")
                axT.set_ylabel("T [10^4 K]")
            if i == len(ds) - 1:
                axT.set_xlabel("r [RT cells]")
            if i == 0 and j == 0:
                ax.legend(fontsize=7)
    fig.suptitle(f"uniform I-front profiles, c_reduced = {a.c:g}, N = {ds[0]['N']}", fontsize=10)
    fig.tight_layout()
    fig.savefig(a.out, dpi=110)


def fig_fields(a):
    import matplotlib.pyplot as plt
    ds = [load(f) for f in a.json]
    keys = [("front_shift_max", "I-front shift, max over steps [cells]"),
            ("xHII_M_last_rel", "|x_HII,M / x_HII,M(k=1) - 1|, last step"),
            ("T_ion_rms", "rms T/T(k=1) - 1, ionised cells, last step")]
    fig, axs = plt.subplots(1, len(keys), figsize=(4.3 * len(keys), 3.4))
    for d, lab in zip(ds, a.labels):
        ks, vals = [], {k: [] for k, _ in keys}
        base = d["runs"]["1"]
        for k, e in sorted(d["runs"].items(), key=lambda x: int(x[0])):
            if k == "1":
                continue
            ks.append(int(k))
            vals["front_shift_max"].append(e["front_shift_max"])
            vals["xHII_M_last_rel"].append(abs(e["xHII_M_last"] / base["xHII_M_last"] - 1))
            vals["T_ion_rms"].append(e["rms_max_last"]["T_ion"]["rms"])
        for ax, (key, _) in zip(axs, keys):
            ax.plot(ks, vals[key], "o-", label=lab)
    for ax, (key, t) in zip(axs, keys):
        ax.set_xscale("log", base=2)
        ax.set_yscale("log")
        ax.set_xlabel("k (RT cells per substep)")
        ax.set_title(t, fontsize=8)
    axs[0].axhline(0.5, color="0.5", ls=":")
    axs[1].axhline(0.01, color="0.5", ls=":")
    axs[2].axhline(0.03, color="0.5", ls=":")
    axs[0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(a.out, dpi=110)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("kind", choices=("ifront", "profiles", "fields"))
    ap.add_argument("--json", nargs="+", required=True)
    ap.add_argument("--labels", nargs="+", required=True)
    ap.add_argument("--c", type=float, default=1.0)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    import matplotlib
    matplotlib.use("Agg")
    {"ifront": fig_ifront, "profiles": fig_profiles, "fields": fig_fields}[a.kind](a)


if __name__ == "__main__":
    main()
