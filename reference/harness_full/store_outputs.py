#!/usr/bin/env python3
"""Store a compact copy of full-harness runs as the CTest reference.

usage: store_outputs.py RUN_DIR OUT_DIR [--fields]

OUT_DIR receives
  step_NNNN/stats.txt     every step
  step_NNNN/sha256.txt    SHA-256 of every .bin field file of that step
  step_NNNN/*.bin         last step only, with --fields (compare.py diagnostics)
  Log/hydro.log           summarize_gas records (all steps)
  Pk/                     power spectra (text)
  DM/sha256.txt           SHA-256 of the write_data_1/2 particle files
tests/compare_full_reference.py checks a C++ run against this.
"""
import glob
import hashlib
import os
import shutil
import sys


def sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def main():
    src, dst = sys.argv[1], sys.argv[2]
    fields = "--fields" in sys.argv[3:]
    shutil.rmtree(dst, ignore_errors=True)
    os.makedirs(dst)
    steps = sorted(glob.glob(os.path.join(src, "step_*")))
    for i, sd in enumerate(steps):
        dd = os.path.join(dst, os.path.basename(sd))
        os.makedirs(dd)
        shutil.copy(os.path.join(sd, "stats.txt"), dd)
        bins = sorted(glob.glob(os.path.join(sd, "*.bin")))
        if bins:
            with open(os.path.join(dd, "sha256.txt"), "w") as f:
                for b in bins:
                    f.write(f"{os.path.basename(b)} {sha(b)}\n")
        if fields and i == len(steps) - 1:
            for b in bins:
                shutil.copy(b, dd)
    os.makedirs(os.path.join(dst, "Log"))
    shutil.copy(os.path.join(src, "Log", "hydro.log"), os.path.join(dst, "Log"))
    if glob.glob(os.path.join(src, "Pk", "*")):
        shutil.copytree(os.path.join(src, "Pk"), os.path.join(dst, "Pk"))
    dm = sorted(glob.glob(os.path.join(src, "DM", "xv_dm.*")))
    if dm:
        os.makedirs(os.path.join(dst, "DM"))
        with open(os.path.join(dst, "DM", "sha256.txt"), "w") as f:
            for b in dm:
                f.write(f"{os.path.basename(b)} {sha(b)}\n")
    print(f"{src}: {len(steps)} steps -> {dst}")


if __name__ == "__main__":
    main()
