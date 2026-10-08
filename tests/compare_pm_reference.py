#!/usr/bin/env python3
"""Strict comparison of a PM run with the stored Fortran reference.

usage: compare_pm_reference.py REF_DIR NEW_DIR

REF_DIR is reference/harness_pm/outputs/<problem> (store_outputs.py):
  step_NNNN/stats.txt   every value must be the identical double
  step_NNNN/sha256.txt  every listed .bin of NEW must have this SHA-256
  step_NNNN/*.bin       byte-identical
  Pk/*                  byte-identical
Exit 0 if identical, 1 if anything differs, 2 on missing files.
"""
import glob
import hashlib
import os
import sys


def parse(fn):
    d = {}
    for line in open(fn):
        if "=" in line:
            k, v = line.split("=", 1)
            d[k.strip()] = float(v.replace("D", "e"))
    return d


def main():
    ref, new = sys.argv[1], sys.argv[2]
    steps = sorted(glob.glob(os.path.join(ref, "step_*")))
    if not steps:
        print("no steps in", ref)
        return 2
    bad, nhash, nval, nfile = [], 0, 0, 0

    def need(path):
        if not os.path.exists(path):
            print("missing", path)
            sys.exit(2)
        return path

    for sd in steps:
        step = os.path.basename(sd)
        nd = os.path.join(new, step)
        a, b = parse(os.path.join(sd, "stats.txt")), parse(need(os.path.join(nd, "stats.txt")))
        for k, v in a.items():
            nval += 1
            if k not in b or b[k] != v:
                bad.append(f"{step}/stats.txt {k}: ref {v!r} new {b.get(k)!r}")
        for line in open(os.path.join(sd, "sha256.txt")):
            name, h = line.split()
            nhash += 1
            hn = hashlib.sha256(open(need(os.path.join(nd, name)), "rb").read()).hexdigest()
            if hn != h:
                bad.append(f"{step}/{name}: SHA-256 differs")
        for f in glob.glob(os.path.join(sd, "*.bin")):
            nfile += 1
            if open(f, "rb").read() != open(need(os.path.join(nd, os.path.basename(f))), "rb").read():
                bad.append(f"{step}/{os.path.basename(f)}: not byte-identical")
    for f in sorted(glob.glob(os.path.join(ref, "Pk", "*"))):
        nfile += 1
        g = need(os.path.join(new, "Pk", os.path.basename(f)))
        if open(f, "rb").read() != open(g, "rb").read():
            bad.append(f"Pk/{os.path.basename(f)}: differs")
    print(f"{ref}: {len(steps)} steps, {nval} stat values, {nhash} field hashes, "
          f"{nfile} files compared")
    for x in bad[:20]:
        print("  DIFF", x)
    print("IDENTICAL" if not bad else f"DIFFERENT ({len(bad)} items)")
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
