#!/usr/bin/env python3
"""Per-cell fields and ray statistics of an RT problem for several
rt_cells_per_substep (docs/rt_multicell.md, step 3).

usage: rt_multicell_fields.py RT_DRIVER PROBLEM_DIR WORKDIR [--nsteps 100]
           [--k 1 2 4] [--c 1] [--json OUT] [--fig PREFIX] [--set key=value ...]

Runs rt_driver on PROBLEM_DIR (legacy_rt = legacy_physics = 0,
nray_max_factor = 16 as tools/ray_compression/run_precision.sh, and
rt_substep_carry = 1) for every k, then compares each k > 1 with k = 1:
  * tools/compare/compare.py on the last step (--json): the maximum
    relative difference |a - b| / max(|a|, |b|, 1e-6 max|a|) of every field;
  * the rms of the same relative difference over all cells and over the
    cells ionised in the k = 1 run (x_HII > 0.5), for x_HI (fHI), x_HII,
    T (temp), the photon counts (ngam, summed over the 5 groups) and the
    photo-ionisation rate (gamma, HI);
  * the volume- and mass-weighted x_HII, T_M and the mean I-front
    displacement (cells whose x_HII > 0.5 state differs / ionised-neutral
    faces of k = 1), at every step;
  * the ray statistics per substep from ray.log (live rays, casts, splits,
    merges; their run means, and the totals over the run, i.e. the work).
"""
import argparse
import json
import os
import shutil
import subprocess
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
COMPARE = HERE.parent / "compare" / "compare.py"


def parse_kv(fn):
    d = {}
    for line in open(fn):
        if "=" in line:
            k, v = line.split("=", 1)
            try:
                d[k.strip()] = [float(x) for x in v.split()]
            except ValueError:
                pass
    return d


def raylog(fn):
    out = {"Rays": [], "Cast": [], "Split": [], "Merge": [], "Delete": []}
    for line in open(fn):
        p = line.split(":")
        if len(p) == 2 and p[0].strip() in out:
            out[p[0].strip()].append(float(p[1].split()[0]))
    return {k: np.array(v) for k, v in out.items()}


def field(step_dir, name, N):
    a = np.fromfile(step_dir / f"{name}.bin", dtype=np.float32).astype(np.float64)
    nc = a.size // N ** 3
    return a.reshape(nc, N ** 3) if nc > 1 else a


def relrms(a, b, sel=None):
    floor = 1e-6 * np.max(np.abs(a))
    rel = np.abs(a - b) / np.maximum(np.maximum(np.abs(a), np.abs(b)), max(floor, 1e-300))
    if sel is not None:
        rel = rel[sel]
    return float(np.sqrt(np.mean(rel ** 2))) if rel.size else 0.0, float(rel.max()) if rel.size else 0.0


def faces(ion, N):
    v = ion.reshape(N, N, N)
    return sum(int(np.count_nonzero(v != np.roll(v, 1, axis=a))) for a in range(3))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("driver")
    ap.add_argument("problem")
    ap.add_argument("work")
    ap.add_argument("--nsteps", type=int, default=100)
    ap.add_argument("--k", type=int, nargs="+", default=[1, 2, 4])
    ap.add_argument("--c", type=float, default=1.0)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--json")
    ap.add_argument("--fig")
    ap.add_argument("--set", action="append", default=[])
    a = ap.parse_args()
    if 1 not in a.k:
        a.k = [1] + a.k
    work = Path(a.work)
    work.mkdir(parents=True, exist_ok=True)
    runs = {}
    for k in a.k:
        out = work / f"k{k}"
        if out.is_dir():
            shutil.rmtree(out)
        cmd = [a.driver, a.problem, str(out), "--nsteps", str(a.nsteps), "--set", "legacy_rt=0",
               "--set", "legacy_physics=0", "--set", "nray_max_factor=16", "--set",
               "rt_substep_carry=1", "--set", f"c_reduced={a.c}", "--set",
               f"rt_cells_per_substep={k}", "--set", f"threads={a.threads}"]
        cmd += sum((["--set", kv] for kv in a.set), [])
        t0 = time.time()
        with open(str(out) + ".stdout", "w") as log:
            subprocess.run(cmd, check=True, stdout=log, stderr=subprocess.STDOUT)
        wall = time.time() - t0
        txt = open(str(out) + ".stdout").read()
        rt_line = [l for l in txt.splitlines() if "ray_tracing" in l or "RT phases" in l]
        runs[k] = {"dir": out, "wall": wall, "ray": raylog(out / "ray.log"), "timers": rt_line[-1:]}
        print(f" k = {k}: {wall:.1f} s", flush=True)
    N = int(parse_kv(Path(a.problem) / "params.txt")["N"][0])
    last = f"step_{a.nsteps:04d}"
    b = runs[1]["dir"]
    res = {"problem": os.path.basename(a.problem.rstrip("/")), "nsteps": a.nsteps, "c": a.c,
           "N": N, "runs": {}}
    hist = {}
    for k in a.k:
        d = runs[k]["dir"]
        xs = []
        for n in range(1, a.nsteps + 1):
            sd = d / f"step_{n:04d}"
            fHI = field(sd, "fHI", N)
            T = field(sd, "temp", N)
            rho = field(sd, "rhog", N)
            st = parse_kv(sd / "stats.txt")
            xs.append(((1 - fHI).mean(), ((1 - fHI) * rho).sum() / rho.sum(),
                       (T * rho).sum() / rho.sum(), st["nstep_rt"][0]))
        hist[k] = np.array(xs)
    for k in a.k:
        r = runs[k]
        ray = r["ray"]
        ent = {"wall": r["wall"], "substeps": float(hist[k][:, 3].sum()),
               "rays_mean": float(ray["Rays"].mean()), "rays_max": float(ray["Rays"].max()),
               "cast_mean": float(ray["Cast"].mean()), "split_mean": float(ray["Split"].mean() / 4),
               "merge_mean": float(ray["Merge"].mean()), "cast_total": float(ray["Cast"].sum()),
               "split_total": float(ray["Split"].sum() / 4), "merge_total": float(ray["Merge"].sum()),
               "xHII_V_last": float(hist[k][-1, 0]), "xHII_M_last": float(hist[k][-1, 1]),
               "T_M_last": float(hist[k][-1, 2])}
        if k != 1:
            d = r["dir"]
            js = work / f"compare_k{k}.json"
            subprocess.run(["python3", str(COMPARE), str(b), str(d), "--steps", str(a.nsteps),
                            "--json", str(js), "-q"], stdout=subprocess.DEVNULL)
            rep = json.load(open(js))
            ent["compare_max_rel"] = {f["field"]: f["max_rel"] for f in rep["fields"]}
            sb, sd = b / last, d / last
            fb, fd = field(sb, "fHI", N), field(sd, "fHI", N)
            ion = (1 - fb) > 0.5
            m = {}
            m["xHI"] = relrms(fb, fd)
            m["xHI_ion"] = relrms(fb, fd, ion)
            m["xHII"] = relrms(1 - fb, 1 - fd)
            Tb, Td = field(sb, "temp", N), field(sd, "temp", N)
            m["T"] = relrms(Tb, Td)
            m["T_ion"] = relrms(Tb, Td, ion)
            nb, nd = field(sb, "ngam", N).sum(axis=0), field(sd, "ngam", N).sum(axis=0)
            m["ngam"] = relrms(nb, nd)
            m["ngam_ion"] = relrms(nb, nd, ion)
            gb, gd = field(sb, "gamma", N)[0], field(sd, "gamma", N)[0]
            m["gamma_ion"] = relrms(gb, gd, ion)
            ent["rms_max_last"] = {kk: {"rms": v[0], "max": v[1]} for kk, v in m.items()}
            # histories
            xb, xk = hist[1], hist[k]
            sel = xb[:, 0] > 1e-3
            ent["xHII_V_rel_max"] = float(np.max(np.abs(xk[sel, 0] / xb[sel, 0] - 1)))
            ent["xHII_M_rel_max"] = float(np.max(np.abs(xk[sel, 1] / xb[sel, 1] - 1)))
            ent["T_M_rel_max"] = float(np.max(np.abs(xk[sel, 2] / xb[sel, 2] - 1)))
            shift = 0.0
            for n in range(1, a.nsteps + 1, max(1, a.nsteps // 20)):
                ib = (1 - field(b / f"step_{n:04d}", "fHI", N)) > 0.5
                ik = (1 - field(d / f"step_{n:04d}", "fHI", N)) > 0.5
                fc = faces(ib, N)
                if fc:
                    shift = max(shift, np.count_nonzero(ib != ik) / fc)
            ib = (1 - fb) > 0.5
            ik = (1 - fd) > 0.5
            ent["front_shift_last"] = float(np.count_nonzero(ib != ik) / max(faces(ib, N), 1))
            ent["front_shift_max"] = float(max(shift, ent["front_shift_last"]))
        res["runs"][k] = ent
    # report
    print(f"\n{res['problem']}: {a.nsteps} steps, c_reduced = {a.c:g}")
    print(" k  substeps  rays/substep  casts/substep  splits/substep  merges/substep  "
          "casts total   wall")
    for k in a.k:
        e = res["runs"][k]
        print(f" {k}  {e['substeps']:8.0f}  {e['rays_mean']:12.0f}  {e['cast_mean']:13.0f}  "
              f"{e['split_mean']:14.0f}  {e['merge_mean']:14.0f}  {e['cast_total']:11.3e}  "
              f"{e['wall']:5.1f}")
    print(" k  x_V rel  x_M rel  T_M rel  front(last/max)  rms|max rel at the last step: "
          "xHI xHI(ion) T T(ion) ngam gamma(ion)")
    for k in a.k:
        if k == 1:
            continue
        e = res["runs"][k]
        m = e["rms_max_last"]
        print(f" {k}  {e['xHII_V_rel_max']:.2e} {e['xHII_M_rel_max']:.2e} {e['T_M_rel_max']:.2e}"
              f"  {e['front_shift_last']:.3f}/{e['front_shift_max']:.3f}  " +
              "  ".join(f"{m[q]['rms']:.2e}|{m[q]['max']:.2e}"
                        for q in ("xHI", "xHI_ion", "T", "T_ion", "ngam", "gamma_ion")))
    if a.json:
        with open(a.json, "w") as f:
            json.dump(res, f, indent=1)
    if a.fig:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        cols = {1: "#1f77b4", 2: "#ff7f0e", 4: "#2ca02c", 8: "#d62728"}
        fig, axs = plt.subplots(1, 3, figsize=(13, 3.6))
        for k in a.k:
            ray = runs[k]["ray"]
            # substeps on a common time axis: substep index * k
            x = np.arange(len(ray["Rays"])) * k
            axs[0].plot(x, ray["Rays"], color=cols.get(k), label=f"k = {k}")
            axs[1].plot(x, ray["Split"] / 4, color=cols.get(k), lw=0.7)
            axs[2].plot(x, ray["Merge"], color=cols.get(k), lw=0.7)
        axs[0].set_ylabel("live rays")
        axs[1].set_ylabel("splits per substep")
        axs[2].set_ylabel("merges per substep")
        for ax in axs:
            ax.set_xlabel("time [k = 1 substeps]")
        axs[0].legend(fontsize=8)
        fig.suptitle(f"{res['problem']}: ray statistics per substep", fontsize=10)
        fig.tight_layout()
        fig.savefig(a.fig + "_rays.png", dpi=110)
        plt.close(fig)
        # slices at the last step
        n = len(a.k)
        fig, axs = plt.subplots(2, n, figsize=(3.3 * n, 6.4), squeeze=False)
        kk = N // 2
        for j, k in enumerate(a.k):
            sd = runs[k]["dir"] / last
            fHI = field(sd, "fHI", N).reshape(N, N, N)
            T = field(sd, "temp", N).reshape(N, N, N)
            im = axs[0, j].imshow(np.log10(np.maximum(fHI[kk], 1e-6)).T, origin="lower", vmin=-6,
                                  vmax=0, cmap="viridis")
            im2 = axs[1, j].imshow(np.log10(np.maximum(T[kk], 1)).T, origin="lower", vmin=2,
                                   vmax=4.7, cmap="inferno")
            axs[0, j].set_title(f"k = {k}", fontsize=9)
            for ax in axs[:, j]:
                ax.set_xticks([])
                ax.set_yticks([])
        fig.colorbar(im, ax=axs[0, :].tolist(), shrink=0.8, label="log x_HI")
        fig.colorbar(im2, ax=axs[1, :].tolist(), shrink=0.8, label="log T [K]")
        fig.suptitle(f"{res['problem']}, step {a.nsteps}, slice k = {kk}", fontsize=10)
        fig.savefig(a.fig + "_slices.png", dpi=110)
        plt.close(fig)


if __name__ == "__main__":
    main()
