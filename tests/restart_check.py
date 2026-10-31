#!/usr/bin/env python3
"""Bit-identical restart of apps/radhydro.

usage: restart_check.py [--launch LAUNCHER FLAGS NPFLAG NP [--restart-np NPB]] RADHYDRO PROBLEM WORKDIR STEP
                        [--set key=value ...]

(--launch: both runs under the MPI launcher; FLAGS a ;-separated list, as
CMake passes it. --restart-np: the restart on NPB ranks instead of NP, a
checkpoint read on another rank count.)

Run A writes a checkpoint after every step (restart_every = 1,
restart_keep_latest = 0) and all fields every step. Run B restarts from A's
checkpoint of step STEP into a fresh directory. Every file B writes (step_*
directories, DM/, Gas/, RT/, Ray/, Pk/, movie frames) must be byte-identical
to A's file of the same name, and B's Log/hydro.log must equal the tail of
A's (also movie/frames.txt). Exit 0 on success.
"""
import glob
import os
import shutil
import subprocess
import sys


def run(cmd, log):
    with open(log, "w") as f:
        r = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
    if r.returncode != 0:
        print(open(log).read()[-3000:])
        sys.exit(f"failed: {' '.join(cmd)}")


def main():
    args = sys.argv[1:]
    launch = launch_b = []
    if args and args[0] == "--launch":
        launcher, flags, npflag, np = args[1:5]
        launch = [launcher] + [f for f in flags.split(";") if f] + [npflag, np]
        launch_b = launch
        args = args[5:]
        if args and args[0] == "--restart-np":
            launch_b = launch[:-1] + [args[1]]
            args = args[2:]
    exe, prob, work, step = args[0], args[1], args[2], int(args[3])
    sets = args[4:]
    a, b = os.path.join(work, "A"), os.path.join(work, "B")
    for d in (a, b):
        if os.path.isdir(d):
            shutil.rmtree(d)
        os.makedirs(d)
    common = sets + ["--set", "fields_every=1"]
    run(launch + [exe, prob, a, "--quiet", "--set", "restart_every=1", "--set", "restart_keep_latest=0"]
        + common, a + ".log")
    ck = os.path.join(a, "restart", f"checkpoint_{step:04d}.bin")
    run(launch_b + [exe, prob, b, "--quiet", "--restart", ck] + common, b + ".log")

    bad, n = [], 0
    for f in sorted(glob.glob(os.path.join(b, "**", "*"), recursive=True)):
        rel = os.path.relpath(f, b)
        if (os.path.isdir(f) or rel in ("timing.txt", "movie/frames.txt")
                or rel.startswith("Log/")):
            continue
        g = os.path.join(a, rel)
        n += 1
        if not os.path.exists(g):
            bad.append(f"{rel}: not in the uninterrupted run")
        elif open(f, "rb").read() != open(g, "rb").read():
            bad.append(f"{rel}: differs")
    for rel in ("Log/hydro.log", "movie/frames.txt"):
        fa, fb = os.path.join(a, rel), os.path.join(b, rel)
        if not os.path.exists(fa):
            continue
        ha = open(fa).read().splitlines()
        hb = [x for x in open(fb).read().splitlines() if not x.startswith("#")]
        n += 1
        if hb and ha[-len(hb):] != hb:
            bad.append(f"{rel}: differs from the tail of the uninterrupted run")
    nsteps = len(glob.glob(os.path.join(b, "step_*")))
    ranks = f" ({launch[-1]} ranks, the restart on {launch_b[-1]})" if launch_b != launch else ""
    print(f"restart from step {step}{ranks}: {nsteps} steps, {n} files compared")
    for x in bad[:20]:
        print("  DIFF", x)
    if nsteps == 0:
        bad.append("the restarted run wrote no steps")
    print("IDENTICAL" if not bad else f"DIFFERENT ({len(bad)} items)")
    sys.exit(0 if not bad else 1)


if __name__ == "__main__":
    main()
