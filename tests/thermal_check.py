#!/usr/bin/env python3
"""thermal_init = recfast in apps/radhydro (docs/thermal_history.md).

usage: thermal_check.py RADHYDRO PROBLEM WORKDIR [--set key=value ...]

Four runs of PROBLEM without RT (one zdata1, 12), hydro_scheme = unsplit,
u_gas and e_gas every step:
  A  thermal_init = recfast, one thread;
  B  as A with two threads;
  C  thermal_init = floor, thermal_floor = 1 (the defaults, set explicitly);
  D  no thermal_* key.
Checks: B byte-identical to A, D to C; at every step with z >= 30 the
median temperature of the mean-density cells (0.9 < Delta < 1.1) within 1 %
of RECFAST (data/thermal/recfast_sim_cosmology.txt) in A, and within 2 % of
the code's floor curve temp_min(z) (the original) in C. (Later the small
test box is non-linear: its mean-density cells are compressed and shocked.)
Exit 0 on success.
"""
import glob
import os
import shutil
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "tools", "analysis"))
import thermal_history as th  # noqa: E402

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


def same_tree(a, b, name):
    fa, fb = files_of(a), files_of(b)
    check(sorted(fa) == sorted(fb), f"{name}: file lists differ")
    nd = [k for k in fa if k in fb and open(fa[k], "rb").read() != open(fb[k], "rb").read()]
    check(not nd, f"{name}: {len(nd)} files differ, e.g. {nd[:5]}")
    print(f"  {name}: {len(fa)} files compared")


def params(prob):
    p = {}
    for line in open(os.path.join(prob, "params.txt")):
        line = line.split("#")[0]
        if "=" in line:
            k, v = (s.strip() for s in line.split("=", 1))
            p[k] = v
    return int(p["N"]), float(p["box"])


def main():
    exe, prob, work = sys.argv[1:4]
    extra = sys.argv[4:]
    N, box = params(prob)
    base = ["--set", "zdata1=12", "--set", "hydro_scheme=unsplit", "--set", "fields_every=1",
            "--set", "fields=u_gas,e_gas"] + extra
    A, B, C, D = (os.path.join(work, n) for n in "ABCD")
    run(exe, prob, A, base + ["--set", "thermal_init=recfast"], 1)
    run(exe, prob, B, base + ["--set", "thermal_init=recfast"], 2)
    run(exe, prob, C, base + ["--set", "thermal_init=floor", "--set", "thermal_floor=1"], 1)
    run(exe, prob, D, base, 1)
    same_tree(A, B, "recfast: threads 1 vs 2")
    same_tree(C, D, "thermal_init = floor, thermal_floor = 1 vs no key")
    Tref = th.reference()
    hA = th.run_history(A, N, box)
    hC = th.run_history(C, N, box)
    check(len(hA) >= 10, f"only {len(hA)} steps")
    worstA = max(abs(Tm / Tref(z) - 1) for z, Tm, _, _ in hA if z >= 30)
    worstC = max(abs(Tm / th.floor(z) - 1) for z, Tm, _, _ in hC if z >= 30)
    print(f"  recfast: {len(hA)} steps, z = {hA[0, 0]:.1f} -> {hA[-1, 0]:.1f}, max |T/T_RECFAST - 1| = {worstA:.4f} at z >= 30")
    print(f"  floor:   max |T/temp_min - 1| = {worstC:.4f}")
    check(worstA < 0.01, f"recfast: mean-density T off RECFAST by {worstA:.4f}")
    check(worstC < 0.02, f"floor: mean-density T off the floor curve by {worstC:.4f}")
    if bad:
        sys.exit(f"{len(bad)} failure(s)")
    print("OK")


if __name__ == "__main__":
    main()
