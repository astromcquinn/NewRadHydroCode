#!/usr/bin/env python3
"""I-front in uniform gas with rays crossing k RT cells per substep
(rt_cells_per_substep, docs/rt_multicell.md).

usage: rt_multicell_ifront.py RT_DRIVER UNIFORM_PROBLEM_DIR WORKDIR
           [--N 32] [--da 2e-4] [--nsteps 15] [--luv 8e28]
           [--k 1 2 4] [--c 0.1 1] [--carry 1] [--threads 1]
           [--json OUT.json] [--plot PREFIX] [--profile-steps 2 7 15]
           [--set key=value ...] [--check]

Builds the problem of tests/ifront_rsla_check.py from UNIFORM_PROBLEM_DIR
(uniform neutral gas at the mean density, z = 10, one source of Luv in the
centre cell; N other than 32 replicates the uniform cell), runs rt_driver
(legacy_rt = legacy_physics = 0) for every (c_reduced, k) and measures
after every hydro step:
  * the I-front radius R = (3 V_ion / 4 pi)^(1/3) in RT cells, V_ion = sum
    of x_HII over the cells;
  * the finite-c solution R = R_inf(t - R/c_hat) (Shapiro et al. 2006),
    R_inf from the c_reduced = 1, k = 1 run, c_hat = k/dtsec_rt cells/s
    (what the code does: k cells per substep);
  * the spherically averaged x_HI and T profiles at --profile-steps;
  * the ray statistics per substep from ray.log (live rays, casts, splits,
    merges) and the emitted photons.
rt_substep_carry = 1 (--carry) makes the total light travel exact over the
run (with k cells per substep and few cells of light per hydro step the
rounding of nstep_rt is otherwise large); light must cross at least k cells
per hydro step, or the nstep_rt >= 1 guard makes it faster than c_reduced
(such runs are flagged).

--check (the multicell CTest): exits 1 unless, for every k > 1 and c, the
front stays within 0.5 cells of k = 1 for R >= 2 cells, matches the finite-c
solution as tests/ifront_rsla_check.py requires (12 %, rms 6 %, R >= 2),
is causal (R <= c_hat t + 1.5), and the emitted photons agree to 1e-6.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import time

import numpy as np

MYR = 3.156e13


def make_problem(src, dst, N, da, nsteps, luv):
    os.makedirs(dst, exist_ok=True)
    g = np.fromfile(os.path.join(src, "gas.bin"), dtype=np.float32).reshape(-1, 9)
    if not np.all(g == g[0]):
        raise SystemExit("the problem gas is not uniform: " + src)
    np.tile(g[0], N ** 3).astype(np.float32).tofile(os.path.join(dst, "gas.bin"))
    with open(os.path.join(dst, "params.txt"), "w") as f:
        f.write(f"N = {N}\nz_init = 10.0\nda_sim = {da}\nnsteps = {nsteps}\n"
                "legacy_rt = 0\nlegacy_physics = 0\n")
    c = N // 2
    with open(os.path.join(dst, "sources.txt"), "w") as f:
        f.write(f"{c} {c} {c} {luv:e} 1.0\n")
    return c  # 1-based source cell


def parse(fn):
    d = {}
    for line in open(fn):
        if "=" in line:
            k, v = line.split("=", 1)
            try:
                d[k.strip()] = [float(x) for x in v.split()]
            except ValueError:
                pass
    return d


def parse_raylog(fn):
    """Per-substep ray statistics of ray.log: live rays, casts, merges,
    splits (the first number of each line)."""
    out = {"rays": [], "cast": [], "merge": [], "split": [], "delete": []}
    key = {"Rays": "rays", "Cast": "cast", "Merge": "merge", "Split": "split",
           "Delete": "delete"}
    for line in open(fn):
        parts = line.split(":")
        if len(parts) != 2:
            continue
        name = parts[0].strip()
        if name in key:
            out[key[name]].append(float(parts[1].split()[0]))
    return {k: np.array(v) for k, v in out.items()}


def radial_profiles(sd, N, src, nbin):
    fHI = np.fromfile(os.path.join(sd, "fHI.bin"), dtype=np.float32).astype(np.float64)
    T = np.fromfile(os.path.join(sd, "temp.bin"), dtype=np.float32).astype(np.float64)
    # Fortran order (i fastest): index = i + N j + N^2 k
    i = np.arange(N) + 0.5
    c0 = src - 0.5
    d = np.abs(i - c0)
    d = np.minimum(d, N - d)
    r = np.sqrt(d[None, None, :] ** 2 + d[None, :, None] ** 2 + d[:, None, None] ** 2).ravel()
    b = np.minimum((r).astype(int), nbin - 1)
    cnt = np.bincount(b, minlength=nbin)
    xhi = np.bincount(b, weights=fHI, minlength=nbin) / np.maximum(cnt, 1)
    tt = np.bincount(b, weights=T, minlength=nbin) / np.maximum(cnt, 1)
    return xhi, tt, fHI, T


def run(driver, prob, out, c, k, nsteps, N, src, carry, threads, extra, profile_steps):
    if os.path.isdir(out):
        shutil.rmtree(out)
    os.makedirs(out)
    cmd = [driver, prob, out, "--set", f"c_reduced={c}", "--set", f"rt_cells_per_substep={k}",
           "--set", f"rt_substep_carry={carry}", "--set", f"threads={threads}"]
    cmd += sum((["--set", kv] for kv in extra), [])
    t0 = time.time()
    with open(out + ".stdout", "w") as log:
        subprocess.run(cmd, check=True, stdout=log, stderr=subprocess.STDOUT)
    wall = time.time() - t0
    guard = "warning: nstep_rt" in open(out + ".stdout").read()
    res = {"c": c, "k": k, "t": [], "R": [], "chat": [], "nstep": [], "ngam": 0.0,
           "wall": wall, "guard": guard, "profiles": {}}
    sec = 0.0
    nbin = N // 2
    for n in range(1, nsteps + 1):
        sd = os.path.join(out, f"step_{n:04d}")
        st = parse(os.path.join(sd, "rt_stats.txt"))
        ns = parse(os.path.join(sd, "stats.txt"))
        dtsec_rt = st["dtsec_rt"][0]
        nstep = ns["nstep_rt"][0]
        sec += dtsec_rt * nstep
        xhi, tt, fHI, T = radial_profiles(sd, N, src, nbin)
        V = np.sum(1 - fHI)
        res["t"].append(sec)
        res["R"].append((3 * V / (4 * np.pi)) ** (1 / 3))
        res["chat"].append(k / dtsec_rt)  # cells per second: k cells per substep
        res["nstep"].append(nstep)
        res["ngam"] += st["ngam"][0]
        if n in profile_steps:
            res["profiles"][n] = {"xHI": xhi.tolist(), "T": tt.tolist(), "t": sec}
    rl = parse_raylog(os.path.join(out, "ray.log"))
    res["raystats"] = {kk: v.tolist() for kk, v in rl.items()}
    for kk in ("t", "R", "chat", "nstep"):
        res[kk] = np.array(res[kk])
    return res


def retarded(t, chat, Rinf_t, Rinf_R):
    """Solve R = R_inf(t - R/chat) by bisection (R_inf^3 interpolated in t)."""
    out = []
    for ti, ci in zip(t, chat):
        lo, hi = 0.0, ci * ti
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            Rr = np.interp(ti - mid / ci, Rinf_t, Rinf_R ** 3, left=0.0) ** (1 / 3)
            if mid > Rr:
                hi = mid
            else:
                lo = mid
        out.append(0.5 * (lo + hi))
    return np.array(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("driver")
    ap.add_argument("uniform")
    ap.add_argument("work")
    ap.add_argument("--N", type=int, default=32)
    ap.add_argument("--da", type=float, default=2e-4)
    ap.add_argument("--nsteps", type=int, default=15)
    ap.add_argument("--luv", type=float, default=8e28)
    ap.add_argument("--k", type=int, nargs="+", default=[1, 2, 4])
    ap.add_argument("--c", type=float, nargs="+", default=[0.1, 1.0])
    ap.add_argument("--carry", type=int, default=1)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--profile-steps", type=int, nargs="+", default=None)
    ap.add_argument("--rmin", type=float, default=2.0)
    ap.add_argument("--json")
    ap.add_argument("--plot")
    ap.add_argument("--set", action="append", default=[])
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    if 1 not in a.k:
        a.k = [1] + a.k
    if 1.0 not in a.c:
        a.c = a.c + [1.0]  # R_inf
    psteps = a.profile_steps or sorted({2, (a.nsteps + 1) // 2, a.nsteps})
    prob = os.path.join(a.work, "problem")
    src = make_problem(a.uniform, prob, a.N, a.da, a.nsteps, a.luv)
    runs = {}
    for c in a.c:
        for k in a.k:
            tag = f"c{c:g}_k{k}"
            runs[(c, k)] = run(a.driver, prob, os.path.join(a.work, tag), c, k, a.nsteps, a.N,
                               src, a.carry, a.threads, a.set, psteps)
            r = runs[(c, k)]
            print(f" ran {tag}: {r['wall']:.1f} s, substeps/step {r['nstep'].min():g}-"
                  f"{r['nstep'].max():g}{'  (nstep_rt >= 1 guard active)' if r['guard'] else ''}",
                  flush=True)
    base = runs[(1.0, 1)]
    Rinf_t = np.concatenate([[0.0], base["t"] - base["R"] / base["chat"]])
    Rinf_R = np.concatenate([[0.0], base["R"]])
    ok = True
    summary = {"N": a.N, "da": a.da, "nsteps": a.nsteps, "luv": a.luv, "carry": a.carry,
               "runs": []}
    for c in a.c:
        r1 = runs[(c, 1)]
        print(f"\n c_reduced = {c:g}: R [cells] per step (pred: finite-c solution)")
        hdr = " step  t[Myr]  " + "  ".join(f"R(k={k})" for k in a.k) + "   pred(k=1)"
        print(hdr)
        preds = {k: retarded(runs[(c, k)]["t"], runs[(c, k)]["chat"], Rinf_t, Rinf_R) for k in a.k}
        for n in range(a.nsteps):
            print(f" {n + 1:4d} {r1['t'][n] / MYR:7.2f}  " +
                  "  ".join(f"{runs[(c, k)]['R'][n]:7.3f}" for k in a.k) +
                  f"   {preds[1][n]:8.3f}")
        for k in a.k:
            r = runs[(c, k)]
            pred = preds[k]
            sel = r["R"] >= a.rmin
            rel = r["R"] / pred - 1
            dR = r["R"] - r1["R"]
            rs = r["raystats"]
            ent = {
                "c": c, "k": k, "wall": r["wall"], "guard": r["guard"],
                "nstep_total": float(r["nstep"].sum()),
                "chat_over_c": float(r["chat"][0] / (base["chat"][0] / 1.0)) if c != 1.0 or k != 1 else 1.0,
                "R": r["R"].tolist(), "t_myr": (r["t"] / MYR).tolist(), "pred": pred.tolist(),
                "rel_err_max": float(np.max(np.abs(rel[sel]))) if sel.any() else 0.0,
                "rel_err_rms": float(np.sqrt(np.mean(rel[sel] ** 2))) if sel.any() else 0.0,
                "dR_vs_k1_max": float(np.max(np.abs(dR[r1["R"] >= a.rmin]))) if k != 1 else 0.0,
                "dR_vs_k1_last": float(dR[-1]),
                "ngam_emitted": r["ngam"],
                "rays_mean": float(np.mean(rs["rays"])) if len(rs["rays"]) else 0.0,
                "rays_last": float(rs["rays"][-1]) if len(rs["rays"]) else 0.0,
                "cast_per_substep": float(np.mean(rs["cast"])) if len(rs["cast"]) else 0.0,
                "split_per_substep": float(np.mean(rs["split"])) / 4 if len(rs["split"]) else 0.0,
                "merge_per_substep": float(np.mean(rs["merge"])) if len(rs["merge"]) else 0.0,
                "cast_total": float(np.sum(rs["cast"])),
                "profiles": r["profiles"],
            }
            summary["runs"].append(ent)
            print(f"   k={k}: finite-c max|rel| {ent['rel_err_max']:.3f} rms {ent['rel_err_rms']:.3f}"
                  f"; max|R-R(k=1)| {ent['dR_vs_k1_max']:.3f} (last {ent['dR_vs_k1_last']:+.3f})"
                  f"; casts {ent['cast_total']:.3g}; wall {r['wall']:.1f} s"
                  f"{'; GUARD' if r['guard'] else ''}")
            if a.check and k != 1:
                if ent["dR_vs_k1_max"] > 0.5:
                    print(f"  FAIL: c={c:g} k={k}: front more than 0.5 cells from k = 1")
                    ok = False
                if ent["rel_err_max"] > 0.12 or ent["rel_err_rms"] > 0.06:
                    print(f"  FAIL: c={c:g} k={k}: finite-c solution missed")
                    ok = False
                if np.any(r["R"] > r["chat"] * r["t"] + 1.5):
                    print(f"  FAIL: c={c:g} k={k}: front faster than c_hat")
                    ok = False
                if abs(r["ngam"] / r1["ngam"] - 1) > 1e-6:
                    print(f"  FAIL: c={c:g} k={k}: emitted photons differ")
                    ok = False
                if r["guard"]:
                    print(f"  FAIL: c={c:g} k={k}: nstep_rt >= 1 guard (light < k cells per step)")
                    ok = False
    if a.json:
        with open(a.json, "w") as f:
            json.dump(summary, f, indent=1)
    if a.plot:
        plot(a, runs, summary, psteps, Rinf_t, Rinf_R)
    if a.check:
        print("PASS" if ok else "FAIL")
        sys.exit(0 if ok else 1)


def plot(a, runs, summary, psteps, Rinf_t, Rinf_R):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    cols = {1: "#1f77b4", 2: "#ff7f0e", 4: "#2ca02c", 8: "#d62728", 16: "#9467bd"}
    # Front radius and residuals
    fig, axs = plt.subplots(2, len(a.c), figsize=(5 * len(a.c), 6.5), sharex="col",
                            gridspec_kw={"height_ratios": [2, 1]}, squeeze=False)
    for j, c in enumerate(a.c):
        ax, axr = axs[0, j], axs[1, j]
        r1 = runs[(c, 1)]
        for k in a.k:
            r = runs[(c, k)]
            pred = retarded(r["t"], r["chat"], Rinf_t, Rinf_R)
            ax.plot(r["t"] / MYR, r["R"], "o-", ms=3, color=cols.get(k), label=f"k = {k}")
            ax.plot(r["t"] / MYR, pred, "--", lw=0.8, color=cols.get(k))
            axr.plot(r["t"] / MYR, r["R"] - r1["R"], "o-", ms=3, color=cols.get(k))
        ax.plot(r1["t"] / MYR, r1["chat"] * r1["t"], ":", color="gray", label=r"$\hat c t$")
        ax.set_ylim(0, 1.15 * max(runs[(1.0, 1)]["R"].max(), r1["R"].max()))
        ax.set_title(f"c_reduced = {c:g}  (dashed: R = R_inf(t - R/c_hat))", fontsize=9)
        ax.set_ylabel("I-front radius [RT cells]")
        ax.legend(fontsize=8)
        axr.axhspan(-0.5, 0.5, color="0.9")
        axr.axhline(0, color="k", lw=0.5)
        axr.set_ylabel("R - R(k=1) [cells]")
        axr.set_xlabel("t [Myr]")
    fig.suptitle(f"uniform I-front, N = {a.N}, da = {a.da:g}: rays crossing k cells per substep",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(a.plot + "_radius.png", dpi=120)
    plt.close(fig)
    # Profiles of x_HI and T at three times
    for c in a.c:
        fig, axs = plt.subplots(2, len(psteps), figsize=(4 * len(psteps), 5.5), sharex=True,
                                squeeze=False)
        for j, n in enumerate(psteps):
            for k in a.k:
                p = runs[(c, k)]["profiles"].get(n)
                if not p:
                    continue
                rr = np.arange(len(p["xHI"])) + 0.5
                axs[0, j].semilogy(rr, np.maximum(p["xHI"], 1e-7), "-", color=cols.get(k),
                                   label=f"k = {k}")
                axs[1, j].plot(rr, np.array(p["T"]) / 1e4, "-", color=cols.get(k))
            t = runs[(c, 1)]["profiles"][n]["t"] / MYR
            axs[0, j].set_title(f"step {n}, t = {t:.1f} Myr", fontsize=9)
            axs[0, j].set_ylim(1e-6, 2)
            axs[1, j].set_xlabel("r [RT cells]")
        axs[0, 0].set_ylabel("x_HI (shell mean)")
        axs[1, 0].set_ylabel("T [10^4 K] (shell mean)")
        axs[0, 0].legend(fontsize=8)
        fig.suptitle(f"uniform I-front profiles, c_reduced = {c:g}, N = {a.N}", fontsize=10)
        fig.tight_layout()
        fig.savefig(a.plot + f"_profiles_c{c:g}.png", dpi=120)
        plt.close(fig)


if __name__ == "__main__":
    main()
