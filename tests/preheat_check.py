#!/usr/bin/env python3
"""The preheating floor T_min / z_floor in the hydro (docs/thermal_history.md section 5).

usage: preheat_check.py RADHYDRO PROBLEM WORKDIR [--set key=value ...]

Two runs of PROBLEM without RT (one zdata1, 20), hydro_scheme = unsplit,
u_gas and e_gas every step:
  A  T_min = 1000, z_floor = 30;
  B  no T_min / z_floor key.
Checks: every step with z >= 30 byte-identical between A and B (the keys
are inert above z_floor); at every step with z < 30 the median
temperature of the mean-density cells (0.9 < Delta < 1.1) in A is between
1000 K cooled adiabatically over that one step (the floor acts at the
step's start, the expansion then cools the gas by ((1+z)/(1+z_prev))^2)
and 1000 K + 1 % (lifted to the floor in the first step after z_floor
and held there), while in B it stays below a tenth of that (the code's
floor curve). Also two threads byte-identical to one in A. Exit 0 on
success.
"""
import glob
import os
import shutil
import subprocess
import sys

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


def same_files(a, b, keys, name):
    fa, fb = files_of(a), files_of(b)
    missing = [k for k in keys if k not in fa or k not in fb]
    check(not missing, f"{name}: missing {missing[:5]}")
    nd = [k for k in keys if k in fa and k in fb and open(fa[k], "rb").read() != open(fb[k], "rb").read()]
    check(not nd, f"{name}: {len(nd)} files differ, e.g. {nd[:5]}")
    print(f"  {name}: {len(keys)} files compared")


def params(prob):
    p = {}
    for line in open(os.path.join(prob, "params.txt")):
        line = line.split("#")[0]
        if "=" in line:
            k, v = (s.strip() for s in line.split("=", 1))
            p[k] = v
    return int(p["N"]), float(p["box"])


def step_dirs(run):  # the steps of th.run_history: those with the gas fields
    return [d for d in sorted(glob.glob(os.path.join(run, "step_*"))) if os.path.exists(os.path.join(d, "u_gas.bin"))]


def main():
    exe, prob, work = sys.argv[1:4]
    extra = sys.argv[4:]
    N, box = params(prob)
    T_pre, z_pre = 1000.0, 30.0
    base = ["--set", "zdata1=20", "--set", "hydro_scheme=unsplit", "--set", "fields_every=1",
            "--set", "fields=u_gas,e_gas"] + extra
    pre = ["--set", f"T_min={T_pre:g}", "--set", f"z_floor={z_pre:g}"]
    A, A2, B = (os.path.join(work, n) for n in ("A", "A2", "B"))
    run(exe, prob, A, base + pre, 1)
    run(exe, prob, A2, base + pre, 2)
    run(exe, prob, B, base, 1)
    hA = th.run_history(A, N, box)
    hB = th.run_history(B, N, box)
    check(len(hA) >= 10, f"only {len(hA)} steps")
    check(len(hA) == len(hB), f"{len(hA)} steps with the floor, {len(hB)} without")
    # the keys are inert above z_floor: the steps there byte-identical
    sa = step_dirs(A)
    above = [i for i, (z, *_) in enumerate(hA) if z >= z_pre]
    below = [i for i, (z, *_) in enumerate(hA) if z < z_pre - 0.5]  # the step crossing z_pre starts above it
    check(len(above) >= 3 and len(below) >= 3, f"{len(above)} steps above z = {z_pre:g}, {len(below)} below")
    keys = []
    for i in above:
        d = os.path.relpath(sa[i], A)
        keys += [os.path.join(d, os.path.basename(p)) for p in glob.glob(os.path.join(sa[i], "*"))
                 if os.path.basename(p) not in SKIP]
    same_files(A, B, keys, f"T_min vs no key, the {len(above)} steps with z >= {z_pre:g}")
    fa = files_of(A)
    same_files(A, A2, sorted(fa), "T_min: threads 1 vs 2")
    # the floor acts on the units of the step's start; the output is after
    # the step's adiabatic expansion: T_pre ((1+z)/(1+z_prev))^2 <= T <= 1.01 T_pre
    lo = min(hA[i, 1] / (T_pre * ((1 + hA[i, 0]) / (1 + hA[i - 1, 0])) ** 2) for i in below)
    hi = max(hA[i, 1] / T_pre for i in below)
    ratioB = max(hB[i, 1] / T_pre for i in below)
    print(f"  with the floor: {len(hA)} steps, z = {hA[0, 0]:.1f} -> {hA[-1, 0]:.1f}; at z < {z_pre:g} "
          f"T / (floor cooled over the step) >= {lo:.4f}, T / floor <= {hi:.4f}")
    print(f"  without:        T / floor <= {ratioB:.4f} at z < {z_pre:g}")
    check(lo > 0.99, f"T_min: mean-density T below the floor cooled over the step ({lo:.4f})")
    check(hi < 1.01, f"T_min: mean-density T above the floor ({hi:.4f})")
    check(ratioB < 0.1, f"no key: mean-density T at {ratioB:.4f} of the floor (the keys leaked)")
    if bad:
        sys.exit(f"{len(bad)} failure(s)")
    print("OK")


if __name__ == "__main__":
    main()
