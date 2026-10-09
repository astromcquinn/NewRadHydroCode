#!/usr/bin/env python3
"""Compare two apps/radhydro runs of the same problem, e.g. the legacy and
the unsplit hydro scheme (docs/hydro_unsplit.md section 8).

    compare_schemes.py RUN_A:label RUN_B:label [--out FIG.png] [--json OUT.json]

Reads every step_*/stats.txt (x_HII(z), T_V, T_M, nstep_rt), timing.txt
(hydro step counts, RT substeps, wall per hydro step, the unsplit
counters), the field outputs at the data1 steps (temp.bin, rhog.bin on the
RT grid: temperature distributions and T-delta medians, before and after
reionisation) and Pk/pk.z=* (gas/DM power ratio). Prints the comparison
tables and draws the figure.
"""
import argparse
import glob
import json
import os
import re
import sys

import numpy as np


def read_stats(run):
    rows = []
    for d in sorted(glob.glob(os.path.join(run, "step_*"))):
        p = os.path.join(d, "stats.txt")
        if not os.path.exists(p):
            continue
        s = {"dir": d}
        for line in open(p):
            if "=" in line:
                k, v = line.split("=", 1)
                try:
                    s[k.strip()] = float(v)
                except ValueError:
                    pass
        rows.append(s)
    return rows


def read_timing(run):
    p = os.path.join(run, "timing.txt")
    steps, extra = [], {}
    for line in open(p):
        if line.startswith("#"):
            if line.startswith("# unsplit"):
                toks = line[1:].split()
                extra = {toks[i]: float(toks[i + 1]) for i in range(1, len(toks) - 1, 2)}
            elif line.startswith("# total"):
                toks = line[1:].split()
                extra.update({"total_" + toks[i]: float(toks[i + 1])
                              for i in range(0, len(toks) - 1, 2)})
            continue
        v = line.split()
        if len(v) >= 12:
            steps.append([float(x) for x in v])
    a = np.array(steps)
    names = ["istep", "z", "nstep_rt", "Nsource", "nray", "t_hydro", "t_pm", "t_halos", "t_rt",
             "t_sum", "t_io", "t_total"]
    return {n: a[:, i] for i, n in enumerate(names)}, extra


def read_field(d, name, N):
    p = os.path.join(d, name + ".bin")
    if not os.path.exists(p):
        return None
    return np.fromfile(p, dtype="<f4").reshape((N, N, N), order="F").astype(np.float64)


def read_pk(run):
    out = {}
    for f in glob.glob(os.path.join(run, "Pk", "pk.z=*")):
        z = float(re.search(r"z=([0-9.]+)", f).group(1))
        out[round(z, 3)] = np.loadtxt(f)
    return out


def N_of(run):
    with open(os.path.join(os.path.dirname(run.rstrip("/")), "params.txt")) as f:
        for line in f:
            if line.split("=")[0].strip() == "N":
                return int(line.split("=")[1].split("#")[0])
    raise SystemExit("N not found in params.txt")


def box_of(run):
    with open(os.path.join(os.path.dirname(run.rstrip("/")), "params.txt")) as f:
        for line in f:
            if line.split("=")[0].strip() == "box":
                return float(line.split("=")[1].split("#")[0])
    return 10.0


def neutral_temperature(d, N, a, box):
    """T [K] of neutral gas from u_gas (and e_gas, the dual-energy field, if
    present): (e/rho) * temp1_unit, as the hydro suite does; for the field
    outputs before the RT starts, where rt_grid temp is not yet set."""
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..",
                                    "tests", "hydro_suite"))
    import cosmo_gen as cg
    u = np.fromfile(os.path.join(d, "u_gas.bin"), dtype="<f4").reshape((5, N, N, N), order="F")
    u = u.astype(np.float64)
    rho = u[0]
    ef = os.path.join(d, "e_gas.bin")
    if os.path.exists(ef):
        e = np.fromfile(ef, dtype="<f4").reshape((N, N, N), order="F").astype(np.float64)
    else:
        e = u[4] - 0.5 * (u[1:4] ** 2).sum(axis=0) / rho
    T = e / rho * cg.units(a, box, N)["temp1_unit"]
    return np.maximum(T, 1.0), rho / rho.mean()


def tdist(d, N, a=None, box=None):
    """temperature distribution on the RT grid at one field output (rt_grid
    temp; before the RT starts it is zero and the neutral-gas temperature
    from u_gas/e_gas is used instead, flagged by 'source')"""
    T = read_field(d, "temp", N)
    rho = read_field(d, "rhog", N)
    if T is None or rho is None:
        return None
    source = "rt_grid"
    if not np.any(T > 0) and os.path.exists(os.path.join(d, "u_gas.bin")) and a is not None:
        T, rho = neutral_temperature(d, N, a, box)
        source = "u_gas/e_gas (neutral)"
    T, rho = T.ravel(), rho.ravel()
    lt = np.log10(T)
    ld = np.log10(rho)
    bins = np.linspace(-1.5, 3, 19)
    idx = np.digitize(ld, bins) - 1
    med = []
    for b in range(len(bins) - 1):
        s = lt[idx == b]
        med.append(float(np.median(s)) if s.size >= 10 else float("nan"))
    hist, edges = np.histogram(lt, bins=np.linspace(-0.5, 5, 56), weights=rho)
    hist = hist / hist.sum()
    under = rho < 1
    return {"source": source, "T_med": float(np.median(T)), "T_med_under": float(np.median(T[under])),
            "T_mw": float(np.sum(rho * T) / np.sum(rho)), "T_vw": float(np.mean(T)),
            "f_T_lt_10": float(np.mean(T < 10)), "f_T_gt_1e4": float(np.mean(T > 1e4)),
            "trho_logd": (0.5 * (bins[1:] + bins[:-1])).tolist(), "trho_logT_med": med,
            "hist": hist.tolist(), "hist_edges": edges.tolist()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs=2)
    ap.add_argument("--out")
    ap.add_argument("--json")
    a = ap.parse_args()
    R = {}
    for spec in a.runs:
        run, label = (spec.split(":") + [None])[:2]
        label = label or os.path.basename(run.rstrip("/"))
        N = N_of(run)
        st = read_stats(run)
        tm, extra = read_timing(run)
        dists = {}
        box = box_of(run)
        for s in st:
            if os.path.exists(os.path.join(s["dir"], "temp.bin")):
                dists[round(s["redshift"], 3)] = tdist(s["dir"], N, s["scalefactor"], box)
        R[label] = {"run": run, "N": N, "stats": st, "timing": tm, "extra": extra,
                    "dists": dists, "pk": read_pk(run)}

    labels = list(R)
    print("== steps and cost")
    print(f"{'':12s} {'hydro steps':>12s} {'RT substeps':>12s} {'wall/hydro step [s]':>20s} "
          f"{'t_hydro/step [s]':>17s} {'t_rt/step [s]':>14s} {'total [s]':>10s}")
    for l in labels:
        tm = R[l]["timing"]
        rt_on = tm["t_rt"] > 0
        print(f"{l:12s} {len(tm['istep']):12d} {int(tm['nstep_rt'][rt_on].sum()):12d} "
              f"{tm['t_total'].mean():20.3f} {tm['t_hydro'].mean():17.3f} "
              f"{tm['t_rt'].mean():14.3f} {tm['t_total'].sum():10.1f}")
        if R[l]["extra"]:
            print("   ", R[l]["extra"])

    print("== ionisation history x_HII (mass-weighted) and temperatures")
    zs = sorted({round(s["redshift"], 2) for l in labels for s in R[l]["stats"]
                 if s["redshift"] <= 15.5}, reverse=True)
    ref = R[labels[0]]["stats"]
    print(f"{'z':>7s} " + " ".join(f"{l + ' xHII_M':>14s} {l + ' T_M':>12s} {l + ' T_V':>12s}"
                                   for l in labels))
    for z in [15, 14, 13, 12, 11, 10, 9, 8, 7.5, 7, 6.5, 6, 5.5]:
        row = f"{z:7.2f} "
        for l in labels:
            st = R[l]["stats"]
            s = min(st, key=lambda s: abs(s["redshift"] - z))
            if abs(s["redshift"] - z) > 0.3:
                row += f"{'-':>14s} {'-':>12s} {'-':>12s} "
                continue
            row += f"{s['xHII_M']:14.4f} {s['T_M']:12.4g} {s['T_V']:12.4g} "
        print(row)

    print("== temperature distribution at the field outputs (RT grid temp, mass-weighted)")
    print(f"{'label':12s} {'z':>6s} {'median T':>9s} {'median T(d<1)':>14s} {'T_mw':>9s} "
          f"{'f(T<10K)':>9s} {'f(T>1e4)':>9s}  source")
    for l in labels:
        for z, dd in sorted(R[l]["dists"].items(), reverse=True):
            if dd:
                print(f"{l:12s} {z:6.2f} {dd['T_med']:9.3g} {dd['T_med_under']:14.3g} "
                      f"{dd['T_mw']:9.3g} {dd['f_T_lt_10']:9.3f} {dd['f_T_gt_1e4']:9.3f}  "
                      f"{dd['source']}")

    print("== gas / DM power (Delta^2 ratio, CIC-deconvolved DM) at k = 2, 5, 10, 20 h/Mpc")
    for l in labels:
        for z, pk in sorted(R[l]["pk"].items(), reverse=True):
            k, dm, gas = pk[:, 0], pk[:, 3], pk[:, 5]
            r = [float(np.interp(kk, k, gas / dm)) for kk in (2, 5, 10, 20)]
            print(f"{l:12s} z={z:6.2f}  " + "  ".join(f"{x:.3f}" for x in r))

    if a.json:
        out = {l: {"N": R[l]["N"], "extra": R[l]["extra"],
                   "steps": len(R[l]["timing"]["istep"]),
                   "rt_substeps": int(R[l]["timing"]["nstep_rt"][R[l]["timing"]["t_rt"] > 0].sum()),
                   "wall_per_step": float(R[l]["timing"]["t_total"].mean()),
                   "t_hydro_per_step": float(R[l]["timing"]["t_hydro"].mean()),
                   "dists": R[l]["dists"],
                   "xHII": [(s["redshift"], s["xHII_M"], s["T_M"], s["T_V"]) for s in R[l]["stats"]]}
               for l in labels}
        with open(a.json, "w") as f:
            json.dump(out, f, indent=1)

    if a.out:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(2, 3, figsize=(12, 7))
        C = ["#1f5fa8", "#d1495b"]
        for i, l in enumerate(labels):
            st = R[l]["stats"]
            z = np.array([s["redshift"] for s in st])
            ax[0, 0].plot(z, [s["xHII_M"] for s in st], color=C[i], label=l)
            ax[0, 1].plot(z, [s["T_M"] for s in st], color=C[i], label=l + " T_M")
            ax[0, 1].plot(z, [s["T_V"] for s in st], color=C[i], ls="--", label=l + " T_V")
            tm = R[l]["timing"]
            ax[0, 2].plot(tm["z"], tm["nstep_rt"], color=C[i], label=l)
            dists = R[l]["dists"]
            if dists:
                zz = sorted(dists, reverse=True)
                for zsel, ls in ((zz[0], "--"), (zz[-1], "-")):
                    dd = dists[zsel]
                    if dd:
                        e = np.array(dd["hist_edges"])
                        ax[1, 0].step(0.5 * (e[1:] + e[:-1]), dd["hist"], where="mid", color=C[i],
                                      ls=ls, label=f"{l} z={zsel:.1f}")
                        ax[1, 1].plot(dd["trho_logd"], dd["trho_logT_med"], color=C[i], ls=ls,
                                      label=f"{l} z={zsel:.1f}")
            pks = R[l]["pk"]
            if pks:
                zl = min(pks)
                pk = pks[zl]
                ax[1, 2].plot(pk[:, 0], pk[:, 5] / pk[:, 3], color=C[i], label=f"{l} z={zl:.1f}")
        ax[0, 0].set(xlabel="z", ylabel="x_HII (mass)", xlim=(16, 5))
        ax[0, 1].set(xlabel="z", ylabel="T [K]", yscale="log", xlim=(16, 5))
        ax[0, 2].set(xlabel="z", ylabel="RT substeps per hydro step", xlim=(16, 5))
        ax[1, 0].set(xlabel="log10 T [K]", ylabel="mass fraction per bin")
        ax[1, 1].set(xlabel="log10 rho/rho_mean", ylabel="median log10 T [K]")
        ax[1, 2].set(xlabel="k [h/Mpc]", ylabel="P_gas / P_dm", xscale="log", ylim=(0, 1.3))
        for x in ax.ravel():
            x.legend(fontsize=7)
            x.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(a.out, dpi=110)
        print("wrote", a.out)


if __name__ == "__main__":
    main()
