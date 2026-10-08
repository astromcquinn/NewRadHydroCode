#!/usr/bin/env python3
"""Strict comparison of a C++ run with Fortran harness outputs.

usage: compare_reference.py REF_DIR NEW_DIR [rtol]

For every step_* directory of REF_DIR and every file in it:
  *.bin                     must be byte-identical
  stats.txt, rt_stats.txt   every numeric value must agree to rtol
                            (default 0: identical doubles; both sides
                            print 17 significant digits)
Files present only in NEW_DIR are ignored (the stored reference subset in
reference/harness/outputs has fields for the last step only).
Exit 0 if everything agrees, 1 otherwise, 2 on missing files.
"""
import glob
import os
import sys


def parse(fn):
    d = {}
    for line in open(fn):
        if "=" not in line:
            continue
        k, v = line.split("=", 1)
        d[k.strip()] = [float(x.replace("D", "e")) for x in v.split()]
    return d


def main():
    ref, new = sys.argv[1], sys.argv[2]
    rtol = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0
    steps = sorted(glob.glob(os.path.join(ref, "step_*")))
    if not steps:
        print("no steps in", ref)
        sys.exit(2)
    nbin = nval = 0
    bad = []
    worst = {}
    for sd in steps:
        step = os.path.basename(sd)
        for fn in sorted(os.listdir(sd)):
            a = os.path.join(sd, fn)
            b = os.path.join(new, step, fn)
            if not os.path.exists(b):
                print("missing", b)
                sys.exit(2)
            if fn.endswith(".bin"):
                nbin += 1
                if open(a, "rb").read() != open(b, "rb").read():
                    bad.append(f"{step}/{fn}: not byte-identical")
            elif fn.endswith(".txt"):
                pa, pb = parse(a), parse(b)
                for k, va in pa.items():
                    vb = pb.get(k)
                    if vb is None or len(vb) != len(va):
                        bad.append(f"{step}/{fn}: key {k} missing/size")
                        continue
                    for x, y in zip(va, vb):
                        nval += 1
                        r = 0.0 if x == y else abs(x - y) / max(abs(x), abs(y))
                        key = f"{fn}:{k}"
                        worst[key] = max(worst.get(key, 0.0), r)
                        if r > rtol:
                            bad.append(f"{step}/{fn}: {k} ref={x!r} new={y!r} rel={r:.3e}")
    print(f"{ref}: {len(steps)} steps, {nbin} binary files, {nval} stat values compared")
    nonzero = {k: v for k, v in worst.items() if v > 0}
    for k, v in sorted(nonzero.items()):
        print(f"  {k:28s} max_rel = {v:.3e}")
    for b in bad[:30]:
        print("  MISMATCH", b)
    if bad:
        print(f"FAILED ({len(bad)} mismatches)")
        sys.exit(1)
    print("IDENTICAL" if not nonzero else "OK (within rtol)")


if __name__ == "__main__":
    main()
