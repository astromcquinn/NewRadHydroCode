#!/usr/bin/env python3
"""Two runs of a driver that differ only in a set of params.txt overrides
must give byte-identical outputs (e.g. hydro_backend = serial against
amrex: the unsplit step on the AMReX level is the same arithmetic).

usage: backend_pair_check.py [--launch LAUNCHER FLAGS NPFLAG NP] WORKDIR DRIVER PROBLEM
                             "k=v,k=v" "k=v,k=v" [common args ...]

Each run is `[LAUNCHER FLAGS NPFLAG NP] DRIVER PROBLEM <WORKDIR>/<n> [common args]
--set k=v ...` (FLAGS: a ;-separated list, as CMake passes it). Timing
files and logs with wall-clock times are skipped. Exit 0 if identical.
"""
import glob
import os
import shutil
import subprocess
import sys

SKIP = ("timing.txt", "ray.log", "rt.log", "stdout.log")


def main():
    args = sys.argv[1:]
    launch = []
    if args and args[0] == "--launch":
        launcher, flags, npflag, np = args[1:5]
        launch = [launcher] + [f for f in flags.split(";") if f] + [npflag, np]
        args = args[5:]
    work, drv, prob, seta, setb = args[:5]
    common = args[5:]
    outs = []
    for n, sets in enumerate((seta, setb)):
        o = os.path.join(work, f"run{n}")
        if os.path.isdir(o):
            shutil.rmtree(o)
        os.makedirs(o)
        cmd = launch + [drv, prob, o] + common
        for kv in sets.split(","):
            if kv:
                cmd += ["--set", kv]
        with open(o + ".log", "w") as f:
            r = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
        if r.returncode != 0:
            print(open(o + ".log").read()[-2000:])
            sys.exit(f"driver failed with {sets}")
        outs.append(o)
    files = sorted(os.path.relpath(f, outs[0])
                   for f in glob.glob(os.path.join(outs[0], "**", "*"), recursive=True)
                   if os.path.isfile(f) and os.path.basename(f) not in SKIP
                   and not os.path.relpath(f, outs[0]).startswith("Log/"))
    if not files:
        sys.exit("no output files")
    bad = []
    for rel in files:
        a = open(os.path.join(outs[0], rel), "rb").read()
        fb = os.path.join(outs[1], rel)
        if not os.path.exists(fb) or open(fb, "rb").read() != a:
            bad.append(rel)
    if bad:
        print("files that differ between", seta, "and", setb, ":")
        for b in bad[:40]:
            print(" ", b)
        sys.exit(f"{len(bad)} of {len(files)} files differ")
    print(f"{len(files)} files identical between '{seta}' and '{setb}'")


if __name__ == "__main__":
    main()
