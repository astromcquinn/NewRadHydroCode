#!/usr/bin/env python3
"""Rank-count independence of apps/radhydro with rt_backend = amrex: run the
driver under the MPI launcher with each rank count and require byte-identical
outputs (step_*, DM/Gas/RT files, movie frames, checkpoints).

usage: rank_check.py WORKDIR LAUNCH "LAUNCH_FLAGS" NPFLAG "NP1 NP2 ..." DRIVER PROBLEM [args ...]
"""
import glob
import os
import shutil
import subprocess
import sys

SKIP = ("timing.txt", "ray.log", "rt.log", "stdout.log")


def main():
    work, launch, flags, npflag, nps = sys.argv[1:6]
    drv, prob = sys.argv[6], sys.argv[7]
    args = sys.argv[8:]
    flags = [f for f in flags.replace(";", " ").split() if f]
    # srun inherits the allocation's --cpus-per-task; give it the test's threads
    if launch.endswith("srun") and not any("cpus-per-task" in f for f in flags):
        flags.append("--cpus-per-task=" + os.environ.get("OMP_NUM_THREADS", "1"))
    outs = {}
    for npr in [int(x) for x in nps.split()]:
        o = os.path.join(work, f"np{npr}")
        if os.path.isdir(o):
            shutil.rmtree(o)
        os.makedirs(o)
        cmd = ([launch] + flags + [npflag, str(npr)] if (launch and npr > 1) else []) + [drv, prob, o] + args
        with open(o + ".log", "w") as f:
            r = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
        if r.returncode != 0:
            print(open(o + ".log").read()[-2000:])
            sys.exit(f"driver failed with {npr} ranks")
        outs[npr] = o
    base = min(outs)
    files = sorted(os.path.relpath(f, outs[base])
                   for f in glob.glob(os.path.join(outs[base], "**", "*"), recursive=True)
                   if os.path.isfile(f) and os.path.basename(f) not in SKIP
                   and not os.path.relpath(f, outs[base]).startswith("Log/"))
    bad = []
    for npr, o in outs.items():
        if npr == base:
            continue
        for rel in files:
            g = os.path.join(o, rel)
            if not os.path.exists(g) or open(g, "rb").read() != open(os.path.join(outs[base], rel), "rb").read():
                bad.append(f"{npr} ranks: {rel}")
    print(f"{len(files)} files compared for ranks {sorted(outs)}")
    for x in bad[:20]:
        print("  DIFF", x)
    print("IDENTICAL" if not bad and files else "DIFFERENT")
    sys.exit(0 if not bad and files else 1)


if __name__ == "__main__":
    main()
