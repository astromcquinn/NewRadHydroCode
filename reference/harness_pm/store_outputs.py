#!/usr/bin/env python3
"""Store a compact copy of harness runs as the reference for CTest.

usage: store_outputs.py RUNS_ROOT OUT_ROOT [problem ...]

For each problem, OUT_ROOT/<p>/ receives
  step_NNNN/stats.txt     every step
  step_NNNN/sha256.txt    SHA-256 of every .bin file of that step
  step_NNNN/{xv_dm,phi}.bin  last step only (for compare.py diagnostics)
  Pk/                     all power spectra (text)
  stdout.log              harness log
tests/compare_pm_reference.py checks a C++ run against this.
"""
import glob
import hashlib
import os
import shutil
import sys


def main():
    runs, out = sys.argv[1], sys.argv[2]
    problems = sys.argv[3:] or sorted(os.listdir(runs))
    for p in problems:
        src, dst = os.path.join(runs, p), os.path.join(out, p)
        shutil.rmtree(dst, ignore_errors=True)
        os.makedirs(dst)
        steps = sorted(glob.glob(os.path.join(src, "step_*")))
        for i, sd in enumerate(steps):
            dd = os.path.join(dst, os.path.basename(sd))
            os.makedirs(dd)
            shutil.copy(os.path.join(sd, "stats.txt"), dd)
            with open(os.path.join(dd, "sha256.txt"), "w") as f:
                for b in sorted(glob.glob(os.path.join(sd, "*.bin"))):
                    h = hashlib.sha256(open(b, "rb").read()).hexdigest()
                    f.write(f"{os.path.basename(b)} {h}\n")
            if i == len(steps) - 1:
                for n in ("xv_dm.bin", "phi.bin"):
                    shutil.copy(os.path.join(sd, n), dd)
        if os.path.isdir(os.path.join(src, "Pk")) and os.listdir(os.path.join(src, "Pk")):
            shutil.copytree(os.path.join(src, "Pk"), os.path.join(dst, "Pk"))
        shutil.copy(os.path.join(src, "stdout.log"), dst)
        print(f"{p}: {len(steps)} steps -> {dst}")


if __name__ == "__main__":
    main()
