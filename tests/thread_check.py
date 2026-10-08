#!/usr/bin/env python3
"""Thread-count independence: run a driver with 1, 2 and 4 OpenMP threads
and require byte-identical outputs.

usage: thread_check.py WORKDIR DRIVER PROBLEM [driver args ...]

DRIVER is called as `DRIVER PROBLEM <WORKDIR>/tN [driver args]` with
OMP_NUM_THREADS = N. Every output file (step_*/stats.txt, rt_stats.txt,
*.bin, Log/hydro.log, Pk, DM, Gas, RT, Ray, movie frames) must be
byte-identical between the thread counts; timing files, logs with wall-clock
times and stdout are skipped. Exit 0 if identical.
"""
import glob
import os
import shutil
import subprocess
import sys

SKIP = ("timing.txt", "ray.log", "rt.log", "stdout.log")


def main():
    work, drv, prob = sys.argv[1], sys.argv[2], sys.argv[3]
    args = sys.argv[4:]
    outs = {}
    for t in (1, 2, 4):
        o = os.path.join(work, f"t{t}")
        if os.path.isdir(o):
            shutil.rmtree(o)
        os.makedirs(o)
        env = dict(os.environ, OMP_NUM_THREADS=str(t))
        with open(o + ".log", "w") as f:
            r = subprocess.run([drv, prob, o] + args, env=env, stdout=f,
                               stderr=subprocess.STDOUT)
        if r.returncode != 0:
            print(open(o + ".log").read()[-2000:])
            sys.exit(f"driver failed with {t} threads")
        outs[t] = o
    files = sorted(os.path.relpath(f, outs[1])
                   for f in glob.glob(os.path.join(outs[1], "**", "*"), recursive=True)
                   if os.path.isfile(f) and os.path.basename(f) not in SKIP)
    bad = []
    for t in (2, 4):
        for rel in files:
            g = os.path.join(outs[t], rel)
            if not os.path.exists(g) or open(g, "rb").read() != open(
                    os.path.join(outs[1], rel), "rb").read():
                bad.append(f"{t} threads: {rel}")
    print(f"{len(files)} files compared for 1, 2 and 4 threads")
    for x in bad[:20]:
        print("  DIFF", x)
    print("IDENTICAL" if not bad and files else "DIFFERENT")
    sys.exit(0 if not bad and files else 1)


if __name__ == "__main__":
    main()
