#!/usr/bin/env python3
"""Compare full radhydro runs that differ in rt_cells_per_substep
(docs/rt_multicell.md, step 4: the production configuration).

usage: rt_multicell_compare.py ROOT [--N 64] [--base base] [--json OUT]
                               [--fig PREFIX] [--slice-z 8]

ROOT holds one radhydro output directory per run: `base` (k = 1, the
reference) and the variants (k2, k4, ...; any other directory with a
movie/frames.txt). Reported per variant against the base:
  * the acceptance metrics of docs/ROADMAP.md with
    tools/ray_compression/precision_metrics.py at the zdata1 field outputs
    (x_HII,V and x_HII,M relative error, T_M, rms T and x_HII per ionised
    cell, mean I-front displacement in cells);
  * the x_HII,V(z) and x_HII,M(z) histories of movie/frames.txt (every
    movie_dz), the variant interpolated to the base frames: maximum
    relative error where x > 1e-3 (and > 0.1), maximum absolute error,
    the redshifts of x_HII,V = 0.5 and 0.99;
  * tau_es (tools/analysis/reion.py) and its relative error;
  * T_M(z) of every hydro step (stats.txt), interpolated, maximum relative
    error at z <= the first zdata1 output;
  * the cost (timing.txt): wall time, RT time, RT time per hydro step with
    RT, RT substeps, ray-substeps (sum over substeps of the live rays, from
    Log/ray.log) and the RT time per ray-substep.
"""
import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "ray_compression"))
import reion  # noqa: E402
import precision_metrics as pm  # noqa: E402


def read_kv(path):
    d = {}
    for line in open(path):
        if "=" in line:
            k, v = line.split("=", 1)
            try:
                d[k.strip()] = float(v.split()[0])
            except (ValueError, IndexError):
                pass
    return d


def step_stats(run):
    rows = []
    for s in sorted(Path(run).glob("step_*")):
        if (s / "stats.txt").exists():
            rows.append(read_kv(s / "stats.txt"))
    return rows


def timing(run):
    f = Path(run) / "timing.txt"
    a = np.loadtxt(f, ndmin=2)
    # istep redshift nstep_rt Nsource nray_active t_hydro t_pm t_halos t_rt t_summaries t_io t_total
    rt = a[:, 8]
    with_rt = a[:, 2] > 0
    return {"wall": float(a[:, 11].sum()), "t_rt": float(rt.sum()),
            "t_rt_per_step": float(rt[with_rt & (a[:, 3] > 0)].mean()) if np.any(with_rt & (a[:, 3] > 0)) else 0.0,
            "substeps": float(a[a[:, 3] > 0, 2].sum()), "hydro_steps": int(len(a)),
            "rt_steps": int(np.count_nonzero(a[:, 3] > 0))}


def ray_substeps(run):
    """sum over RT substeps of the live rays (ray.log 'Rays' lines)"""
    f = Path(run) / "Log" / "ray.log"
    tot, n, mx = 0.0, 0, 0.0
    if not f.exists():
        return 0.0, 0, 0.0
    for line in open(f):
        if line.lstrip().startswith("Rays") and ":" in line:
            vals = line.split(":")[1].split()
            if len(vals) != 2:  # calc_ngam_ray threshold lines
                continue
            v = float(vals[0])
            tot += v
            n += 1
            mx = max(mx, v)
    return tot, n, mx


def ray_counts(run):
    """per-substep casts, splits, merges and live rays (ray.log)"""
    f = Path(run) / "Log" / "ray.log"
    out = {"Rays": [], "Cast": [], "Split": [], "Merge": []}
    if not f.exists():
        return out
    for line in open(f):
        p = line.split(":")
        if len(p) == 2 and p[0].strip() in out:
            v = p[1].split()
            if p[0].strip() == "Rays" and len(v) != 2:
                continue
            out[p[0].strip()].append(float(v[0]))
    return {k: np.array(v) for k, v in out.items()}


def zcross(z, x, lev):
    o = np.argsort(z)[::-1]  # descending z
    zs, xs = z[o], x[o]
    above = np.nonzero(xs >= lev)[0]
    if not len(above):
        return float("nan")
    k = above[0]
    if k == 0:
        return float(zs[0])
    return float(np.interp(lev, [xs[k - 1], xs[k]], [zs[k - 1], zs[k]]))


def history_metrics(base, var):
    _, zb, xVb, xMb = reion.read_frames(base, valid_only=True)
    _, zv, xVv, xMv = reion.read_frames(var, valid_only=True)
    ob, ov = np.argsort(zb), np.argsort(zv)
    m = {}
    for lab, xb, xv in (("V", xVb, xVv), ("M", xMb, xMv)):
        xi = np.interp(zb[ob], zv[ov], xv[ov])
        xbb = xb[ob]
        inside = (zb[ob] >= zv.min()) & (zb[ob] <= zv.max())
        for thr, key in ((1e-3, ""), (0.1, "_x01")):
            sel = inside & (xbb > thr)
            m[f"xHII_{lab}_rel_max{key}"] = float(np.max(np.abs(xi[sel] / xbb[sel] - 1))) if sel.any() else 0.0
        m[f"xHII_{lab}_abs_max"] = float(np.max(np.abs(xi[inside] - xbb[inside])))
    for lev in (0.1, 0.5, 0.9, 0.99):
        m[f"z_x{lev}_base"] = zcross(zb, xVb, lev)
        m[f"z_x{lev}"] = zcross(zv, xVv, lev)
    zz, tau = reion.tau_history(zb, xMb)
    zz2, tau2 = reion.tau_history(zv, xMv)
    m["tau_es_base"] = float(tau[-1])
    m["tau_es"] = float(tau2[-1])
    m["tau_es_rel"] = float(tau2[-1] / tau[-1] - 1)
    return m


def tm_metrics(base, var, zmax):
    sb, sv = step_stats(base), step_stats(var)
    zb = np.array([r["redshift"] for r in sb])
    zv = np.array([r["redshift"] for r in sv])
    out = {}
    for key in ("T_M", "T_V"):
        tb = np.array([r.get(key, np.nan) for r in sb])
        tv = np.array([r.get(key, np.nan) for r in sv])
        ob, ov = np.argsort(zb), np.argsort(zv)
        ti = np.interp(zb[ob], zv[ov], tv[ov])
        sel = (zb[ob] <= zmax) & (zb[ob] >= zv.min()) & np.isfinite(ti) & np.isfinite(tb[ob])
        out[f"{key}_rel_max"] = float(np.max(np.abs(ti[sel] / tb[ob][sel] - 1))) if sel.any() else 0.0
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root")
    ap.add_argument("--N", type=int, default=64)
    ap.add_argument("--base", default="base")
    ap.add_argument("--json")
    ap.add_argument("--fig")
    ap.add_argument("--zmax", type=float, default=15.0, help="T_M compared at z <= zmax (RT on)")
    a = ap.parse_args()
    root = Path(a.root)
    base = root / a.base
    variants = sorted(p for p in root.iterdir()
                      if p.is_dir() and p.name != a.base and (p / "movie" / "frames.txt").exists())
    bfields = pm.load_run(base, a.N)
    res = {"base": {"timing": timing(base)}}
    rs, nsub, rmax = ray_substeps(base)
    res["base"]["timing"].update(ray_substeps=rs, ray_log_substeps=nsub, rays_max=rmax)
    res["base"]["history"] = history_metrics(base, base)
    for v in variants:
        r = {}
        vf = pm.load_run(v, a.N)
        r["fields"], _ = pm.metrics(bfields, vf)
        r["history"] = history_metrics(base, v)
        r["T"] = tm_metrics(base, v, a.zmax)
        r["timing"] = timing(v)
        rs, nsub, rmax = ray_substeps(v)
        r["timing"].update(ray_substeps=rs, ray_log_substeps=nsub, rays_max=rmax)
        res[v.name] = r
    # table
    b = res["base"]
    print(f"base: tau_es {b['history']['tau_es_base']:.5f}, z(x_V = 0.5) "
          f"{b['history']['z_x0.5_base']:.3f}, z(0.99) {b['history']['z_x0.99_base']:.3f}; "
          f"wall {b['timing']['wall']:.0f} s, RT {b['timing']['t_rt']:.0f} s, "
          f"{b['timing']['substeps']:.0f} substeps")
    cols = ("xHII_V_rel_max", "xHII_M_rel_max", "xHII_V_rel_max_x01", "xHII_V_abs_max", "tau_es_rel")
    print("variant " + " ".join(f"{c:>18s}" for c in cols) +
          "   T_M     Tion_rms  front  dz(0.5) dz(0.99)  wall  t_rt  t_rt/step")
    for v in variants:
        r = res[v.name]
        h, f, t, T = r["history"], r["fields"], r["timing"], r["T"]
        print(f"{v.name:7s} " + " ".join(f"{h[c]:18.3e}" for c in cols) +
              f"  {T['T_M_rel_max']:.3e} {f.get('T_ion_cell_rms', float('nan')):.3e}"
              f" {f['front_shift']:.3f} {h['z_x0.5'] - h['z_x0.5_base']:+.3f}"
              f" {h['z_x0.99'] - h['z_x0.99_base']:+.3f}"
              f" {t['wall']:6.0f} {t['t_rt']:5.0f} {t['t_rt_per_step']:6.2f}")
    if a.json:
        with open(a.json, "w") as fjs:
            json.dump(res, fjs, indent=1)
    if a.fig:
        figures(a, base, variants, bfields)


def figures(a, base, variants, bfields):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    runs = [base] + variants
    cols = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"]
    fig, axs = plt.subplots(3, 1, figsize=(6, 8.5), sharex=True,
                            gridspec_kw={"height_ratios": [2, 1, 1]})
    _, zb, xVb, xMb = reion.read_frames(base, valid_only=True)
    ob = np.argsort(zb)
    sb = step_stats(base)
    zsb = np.array([r["redshift"] for r in sb])
    tmb = np.array([r["T_M"] for r in sb])
    for i, run in enumerate(runs):
        _, z, xV, xM = reion.read_frames(run, valid_only=True)
        o = np.argsort(z)
        lab = "k = 1" if run == base else run.name.replace("k", "k = ")
        axs[0].plot(z[o], xV[o], "-", color=cols[i % 6], label=lab + " (vol)")
        axs[0].plot(z[o], xM[o], "--", color=cols[i % 6], lw=0.8)
        xi = np.interp(zb[ob], z[o], xV[o])
        sel = xVb[ob] > 1e-3
        axs[1].plot(zb[ob][sel], xi[sel] / xVb[ob][sel] - 1, "-", color=cols[i % 6])
        st = step_stats(run)
        zs = np.array([r["redshift"] for r in st])
        tm = np.array([r["T_M"] for r in st])
        os_ = np.argsort(zs)
        osb = np.argsort(zsb)
        ti = np.interp(zsb[osb], zs[os_], tm[os_])
        sel = zsb[osb] <= a.zmax
        axs[2].plot(zsb[osb][sel], ti[sel] / tmb[osb][sel] - 1, "-", color=cols[i % 6])
    axs[0].set_ylabel("x_HII (solid: volume, dashed: mass)")
    axs[0].legend(fontsize=8)
    axs[1].axhspan(-0.01, 0.01, color="0.9")
    axs[1].set_ylabel("x_HII,V / x_HII,V(k=1) - 1")
    axs[1].set_ylim(-0.1, 0.1)
    axs[2].axhspan(-0.03, 0.03, color="0.9")
    axs[2].set_ylabel("T_M / T_M(k=1) - 1")
    axs[2].set_xlabel("z")
    axs[2].invert_xaxis()
    fig.suptitle(f"{a.N}^3 production configuration: rays crossing k cells per substep", fontsize=10)
    fig.tight_layout()
    fig.savefig(a.fig + "_history.png", dpi=120)
    plt.close(fig)
    # slices at the zdata1 output nearest x_HII,V = 0.5
    xs = [1 - b["fHI"].mean() for b in bfields]
    i5 = int(np.argmin(np.abs(np.array(xs) - 0.5)))
    vfs = [pm.load_run(v, a.N) for v in variants]
    n = 1 + len(variants)
    fig, axs = plt.subplots(2, n, figsize=(3.2 * n, 6.2), squeeze=False)
    k0 = a.N // 2
    for j, (lab, fl) in enumerate([("k = 1", bfields)] + [(v.name, f) for v, f in zip(variants, vfs)]):
        if i5 >= len(fl):
            continue
        s = fl[i5]
        im = axs[0, j].imshow(np.log10(np.maximum(s["fHI"][k0], 1e-6)).T, origin="lower",
                              vmin=-6, vmax=0, cmap="viridis")
        axs[0, j].set_title(f"{lab}: z = {s['z']:.2f}, x_V = {1 - s['fHI'].mean():.3f}", fontsize=8)
        im2 = axs[1, j].imshow(np.log10(np.maximum(s["T"][k0], 1)).T, origin="lower",
                               vmin=1, vmax=4.7, cmap="inferno")
        for ax in axs[:, j]:
            ax.set_xticks([])
            ax.set_yticks([])
    fig.colorbar(im, ax=axs[0, :].tolist(), shrink=0.8, label="log x_HI")
    fig.colorbar(im2, ax=axs[1, :].tolist(), shrink=0.8, label="log T [K]")
    fig.savefig(a.fig + "_slices.png", dpi=110)
    plt.close(fig)


if __name__ == "__main__":
    main()
