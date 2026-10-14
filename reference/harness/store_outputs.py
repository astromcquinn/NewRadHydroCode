#!/usr/bin/env python3
"""Store the compact subset of RT-harness runs that the CTests compare with.

usage: store_outputs.py RUN_ROOT OUT_ROOT [problem ...]

For every problem directory RUN_ROOT/<p> (default: all of them) OUT_ROOT/<p>
receives the subset README.md describes: stats.txt and rt_stats.txt of
every step, the binary fields of the last step only (none for
uniform_8src_smallbuf, whose fields equal uniform_8src's), ray.log, rt.log,
stdout.log and Log/. OUT_ROOT/<p> is replaced. The reference sets live in
outputs/, outputs_legacy_rt0/, outputs_legacy_rt0_phys0/ and
outputs_c0.1_fixed/; .github/workflows/regen_references.yml regenerates them
on the CI platform (docs/cluster.md section 3).
"""
import glob
import os
import shutil
import sys

NO_FIELDS = {"uniform_8src_smallbuf"}


def store(run, out, name):
    shutil.rmtree(out, ignore_errors=True)
    os.makedirs(out)
    steps = sorted(glob.glob(os.path.join(run, "step_*")))
    for i, sd in enumerate(steps):
        dd = os.path.join(out, os.path.basename(sd))
        os.makedirs(dd)
        for f in ("stats.txt", "rt_stats.txt"):
            shutil.copy(os.path.join(sd, f), dd)
        if i == len(steps) - 1 and name not in NO_FIELDS:
            for b in sorted(glob.glob(os.path.join(sd, "*.bin"))):
                shutil.copy(b, dd)
    for f in ("ray.log", "rt.log", "stdout.log"):
        if os.path.exists(os.path.join(run, f)):
            shutil.copy(os.path.join(run, f), out)
    if os.path.isdir(os.path.join(run, "Log")):
        shutil.copytree(os.path.join(run, "Log"), os.path.join(out, "Log"))
    print(f"{run}: {len(steps)} steps -> {out}")


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    run_root, out_root = sys.argv[1], sys.argv[2]
    names = sys.argv[3:] or sorted(
        d for d in os.listdir(run_root) if os.path.isdir(os.path.join(run_root, d, "step_0001")))
    for p in names:
        store(os.path.join(run_root, p), os.path.join(out_root, p), p)


if __name__ == "__main__":
    main()
