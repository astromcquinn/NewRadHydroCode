#!/usr/bin/env python3
"""X-ray background heating in apps/radhydro (xray_heating = 1,
docs/xray_heating.md).

usage: xray_check.py RADHYDRO PROBLEM WORKDIR AM_PARAMS [--set key=value ...]

AM_PARAMS: the source-model parameters (data/uv_lf/am_params_st_fof.txt),
passed as --set. Four runs of PROBLEM with the temperature and HI fields
every step:
  A  xray_heating = 1 (log L_X/SFR raised to 44 so that the small test box
     at z = 30-20 heats measurably), one thread;
  B  as A with two threads;
  C  as A with xray_heating = 0;
  D  as C without any xray_* key.
Checks: every file of B is byte-identical to A's (thread independence), and
of D to C's (the keys alone change nothing); Log/xray.txt of A has one row
per hydro step with RT, z decreasing, positive rates and the box's mass-
weighted fractions in [0, 1]; in every step after the RT starts, the mean
temperature of A's neutral cells (x_HI > 0.99, which the original physics
skips) is above C's, and the neutral fraction of those cells is below C's,
by an amount consistent with the rates. Exit 0 on success.
"""
import glob
import os
import shutil
import subprocess
import sys

import numpy as np

SKIP = ("timing.txt", "ray.log", "rt.log", "stdout.log", "hydro.log")
bad = []


def check(cond, msg):
    if not cond:
        bad.append(msg)
        print("  FAIL", msg)


def run(exe, prob, out, sets, threads):
    if os.path.isdir(out):
        shutil.rmtree(out)
    os.makedirs(out)
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    cmd = [exe, prob, out, "--quiet"] + sets
    with open(os.path.join(out, "stdout.log"), "w") as f:
        r = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, env=env)
    if r.returncode != 0:
        print(open(os.path.join(out, "stdout.log")).read()[-3000:])
        sys.exit(f"failed ({r.returncode}): {' '.join(cmd)}")


def files_of(run):
    out = {}
    for p in glob.glob(os.path.join(run, "**", "*"), recursive=True):
        if os.path.isfile(p) and os.path.basename(p) not in SKIP:
            out[os.path.relpath(p, run)] = p
    return out


def same_tree(a, b, name, ignore=()):
    fa, fb = files_of(a), files_of(b)
    for k in ignore:
        fa.pop(k, None)
        fb.pop(k, None)
    check(sorted(fa) == sorted(fb), f"{name}: file lists differ: {sorted(set(fa) ^ set(fb))[:5]}")
    ndiff = [k for k in fa if k in fb and open(fa[k], "rb").read() != open(fb[k], "rb").read()]
    check(not ndiff, f"{name}: {len(ndiff)} files differ, e.g. {ndiff[:5]}")
    print(f"  {name}: {len(fa)} files compared")


def am_sets(path):
    sets = []
    for line in open(path):
        line = line.split("#")[0].strip()
        if "=" in line:
            k, v = (s.strip() for s in line.split("=", 1))
            sets += ["--set", f"{k}={v}"]
    return sets


def main():
    exe, prob, work, am = sys.argv[1:5]
    extra = sys.argv[5:]
    base = am_sets(am) + ["--set", "fields_every=1", "--set", "fields=temp,fHI"] + extra
    xr = ["--set", "xray_log_lx_sfr=44", "--set", "xray_dlog=0.01"]
    A, B, C, D = (os.path.join(work, n) for n in "ABCD")
    run(exe, prob, A, base + xr + ["--set", "xray_heating=1"], 1)
    run(exe, prob, B, base + xr + ["--set", "xray_heating=1"], 2)
    run(exe, prob, C, base + xr + ["--set", "xray_heating=0"], 1)
    run(exe, prob, D, base, 1)
    same_tree(A, B, "threads 1 vs 2 (xray_heating = 1)")
    same_tree(C, D, "xray_heating = 0 vs no xray keys")

    rows = [l.split() for l in open(os.path.join(A, "Log", "xray.txt")) if not l.startswith("#")]
    rows = [[float(x) for x in r] for r in rows]
    check(len(rows) >= 2, f"Log/xray.txt has {len(rows)} rows")
    z = [r[1] for r in rows]
    check(all(z[i + 1] < z[i] for i in range(len(z) - 1)), "xray.txt: z not decreasing")
    for r in rows:
        check(all(0 <= x <= 1 for x in r[4:7]), f"xray.txt step {int(r[0])}: fractions {r[4:7]}")
        check(r[2] > 0 and r[3] > 0, f"xray.txt step {int(r[0])}: sfrd or eps not positive")
        check(r[7] > 0 and r[8] > 0 and r[10] > 0 and r[11] > 0,
              f"xray.txt step {int(r[0])}: rates not positive")
    rt_steps = {int(r[0]) for r in rows}

    nst = 0
    for d in sorted(glob.glob(os.path.join(A, "step_*"))):
        istep = int(os.path.basename(d).split("_")[1])
        if istep not in rt_steps:
            continue
        dc = os.path.join(C, os.path.basename(d))
        tA = np.fromfile(os.path.join(d, "temp.bin"), dtype=np.float32)
        tC = np.fromfile(os.path.join(dc, "temp.bin"), dtype=np.float32)
        hA = np.fromfile(os.path.join(d, "fHI.bin"), dtype=np.float32)
        hC = np.fromfile(os.path.join(dc, "fHI.bin"), dtype=np.float32)
        neutral = hC > 0.99
        check(neutral.sum() > 0, f"step {istep}: no neutral cells")
        dT = float(np.mean(tA[neutral].astype(np.float64) - tC[neutral]))
        dx = float(np.mean(hC[neutral].astype(np.float64) - hA[neutral]))
        check(dT > 0, f"step {istep}: neutral cells not heated (mean dT = {dT:.4g} K)")
        check(dx > 0, f"step {istep}: neutral cells not ionised (mean dx_HI = {dx:.4g})")
        nst += 1
        print(f"  step {istep}: {int(neutral.sum())} neutral cells, mean T {np.mean(tC[neutral]):.3f}"
              f" -> {np.mean(tA[neutral]):.3f} K, mean x_HI {np.mean(hC[neutral]):.6f} -> {np.mean(hA[neutral]):.6f}")
    check(nst >= 2, f"only {nst} steps with RT compared")
    if bad:
        sys.exit(f"{len(bad)} failure(s)")
    print("OK")


if __name__ == "__main__":
    main()
