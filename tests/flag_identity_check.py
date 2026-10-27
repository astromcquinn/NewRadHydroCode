#!/usr/bin/env python3
"""A flag that must change nothing in a given configuration: run apps/radhydro
with FLAG=0 and FLAG=1 (under the MPI launcher) and require byte-identical
outputs (step_*, DM/Gas/RT files, movie frames, checkpoints).

usage: flag_identity_check.py WORKDIR LAUNCH "LAUNCH_FLAGS" NPFLAG NP DRIVER PROBLEM FLAG [args ...]
"""
import filecmp
import os
import shutil
import subprocess
import sys

SKIP = ("timing.txt", "ray.log", "rt.log", "stdout.log")


def main():
    work, launch, flags, npflag, np_, drv, prob, flag = sys.argv[1:9]
    args = sys.argv[9:]
    flags = [f for f in flags.replace(";", " ").split() if f]
    if launch.endswith("srun") and not any("cpus-per-task" in f for f in flags):
        flags.append("--cpus-per-task=" + os.environ.get("OMP_NUM_THREADS", "1"))
    outs = {}
    for v in (0, 1):
        o = os.path.join(work, f"{flag}{v}")
        if os.path.isdir(o):
            shutil.rmtree(o)
        os.makedirs(o)
        cmd = [launch] + flags + [npflag, np_, drv, prob, o] + args + ["--set", f"{flag}={v}"]
        log = os.path.join(work, f"{flag}{v}.log")
        with open(log, "w") as f:
            r = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
        if r.returncode != 0:
            print(open(log).read()[-3000:])
            sys.exit(f"failed: {' '.join(cmd)}")
        files = {}
        for root, _, names in os.walk(o):
            for n in names:
                if n in SKIP:
                    continue
                p = os.path.join(root, n)
                files[os.path.relpath(p, o)] = p
        outs[v] = files
    if set(outs[0]) != set(outs[1]):
        sys.exit(f"different file sets: {sorted(set(outs[0]) ^ set(outs[1]))}")
    bad = [k for k in sorted(outs[0]) if not filecmp.cmp(outs[0][k], outs[1][k], shallow=False)]
    print(f"{len(outs[0])} files compared for {flag} = 0 and 1")
    if bad:
        sys.exit("DIFFERENT: " + " ".join(bad[:20]))
    print("IDENTICAL")


if __name__ == "__main__":
    main()
