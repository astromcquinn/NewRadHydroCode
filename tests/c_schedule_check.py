#!/usr/bin/env python3
"""c_reduced_schedule (apps/radhydro, docs/production_run.md section 1).

usage: c_schedule_check.py [--launch LAUNCHER FLAGS NPFLAG NP] RADHYDRO PROBLEM WORKDIR
                           [--set key=value ...]

With --launch (as backend_pair_check.py) the runs go through MPI: on the
AMReX backends every rank must take the speed rank 0 takes from x.

1. Identity: a one-speed schedule ("C x C") gives the run of the constant
   c_reduced = C byte for byte: every file the same, except that stats.txt
   gains its c_reduced line.
2. The rule: with "0.3 0.005 0.6 0.01 1" every step's c_reduced (stats.txt)
   is c0 = 0.3 below x = 0.005, 0.6 up to 0.01 and 1 above, where x is the
   x_HII,V of the previous step's stats.txt, or 0 while the RT has not yet
   written its grid (x_HII,V = 1 exactly there: fHI is still zero). At least
   two switches must happen (in the smoke problem at steps 14 and 15), and
   every change is announced on stdout.
Exit 0 on success.
"""
import glob
import os
import shutil
import subprocess
import sys

SCHED = [0.3, 0.005, 0.6, 0.01, 1.0]


def run(cmd, log):
    with open(log, "w") as f:
        r = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
    if r.returncode != 0:
        print(open(log).read()[-3000:])
        sys.exit(f"failed: {' '.join(cmd)}")
    return open(log).read()


def files(d):
    out = {}
    for p in glob.glob(os.path.join(d, "**", "*"), recursive=True):
        if os.path.isfile(p) and os.path.relpath(p, d) != "timing.txt":  # wall times
            out[os.path.relpath(p, d)] = p
    return out


def stats(path):
    d = {}
    for line in open(path):
        if " = " in line:
            k, v = line.split(" = ", 1)
            d[k.strip()] = v.strip()
    return d


def rule(x):
    c = SCHED[0]
    for i in range(1, len(SCHED), 2):
        if x >= SCHED[i]:
            c = SCHED[i + 1]
    return c


def main():
    args = sys.argv[1:]
    launch = []
    if args and args[0] == "--launch":
        launcher, flags, npflag, np = args[1:5]
        launch = [launcher] + [f for f in flags.split(";") if f] + [npflag, np]
        args = args[5:]
    exe, prob, work = args[:3]
    sets = args[3:]
    if os.path.isdir(work):
        shutil.rmtree(work)
    os.makedirs(work)
    common = launch + [exe, prob]
    opts = ["--quiet", "--set", "fields_every=1"] + sets
    a, b, c = (os.path.join(work, n) for n in ("const", "sched1", "sched"))
    run(common + [a] + opts + ["--set", "c_reduced=0.3"], a + ".log")
    run(common + [b] + opts + ["--set", "c_reduced_schedule=0.3 0.5 0.3"], b + ".log")
    fa, fb = files(a), files(b)
    bad = sorted(set(fa) ^ set(fb))
    if bad:
        sys.exit(f"identity: different file sets: {bad[:10]}")
    for k in sorted(fa):
        da, db = open(fa[k], "rb").read(), open(fb[k], "rb").read()
        if os.path.basename(k) == "stats.txt":
            lb = [l for l in db.decode().splitlines() if not l.startswith("c_reduced = ")]
            if lb != da.decode().splitlines():
                sys.exit(f"identity: {k} differs beyond the c_reduced line")
        elif da != db:
            sys.exit(f"identity: {k} differs")
    print(f"identity: {len(fa)} files byte-identical (stats.txt apart from c_reduced)")

    log = run(common + [c] + opts + ["--set", "c_reduced_schedule=" + " ".join(map(str, SCHED))],
              c + ".log")
    steps = sorted(glob.glob(os.path.join(c, "step_*", "stats.txt")))
    st = [stats(p) for p in steps]
    used, prev_x = [], None
    for s in st:
        istep = int(s["istep"])
        if istep == 0:
            prev_x = float(s["xHII_V"])
            continue
        x = 0.0 if prev_x == 1.0 else prev_x
        want = rule(x)
        got = float(s["c_reduced"])
        if got != want:
            sys.exit(f"step {istep}: c_reduced = {got}, rule gives {want} (x = {x})")
        used.append(got)
        prev_x = float(s["xHII_V"])
    switches = sum(1 for i in range(1, len(used)) if used[i] != used[i - 1])
    if switches < 2:
        sys.exit(f"only {switches} switches: {used}")
    ann = [l for l in log.splitlines() if "c_reduced_schedule:" in l]
    if len(ann) != switches + 1:
        sys.exit(f"{len(ann)} announcements for {switches} switches")
    print(f"rule: {len(used)} steps, speeds {sorted(set(used))}, {switches} switches")


if __name__ == "__main__":
    main()
