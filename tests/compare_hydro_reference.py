#!/usr/bin/env python3
"""Strict comparison of a hydro run with Fortran hydro-harness outputs.

usage: compare_hydro_reference.py REF_DIR NEW_DIR

For every file of every step_* directory of REF_DIR (the stored subset in
reference/harness_hydro/outputs, or a full harness run):
  *.bin              byte-identical
  stats.txt          every number identical (both sides print 17 digits)
  gas_summary.txt    byte-identical (the summarize_gas "(50es13.5)" record)
and REF_DIR/domain_gas.txt byte-identical.  Files only in NEW_DIR are ignored.
Exit 0 if identical, 1 on a mismatch, 2 on missing files.
"""
import math
import os
import sys


def parse(fn):
    d = {}
    for line in open(fn):
        if "=" in line:
            k, v = line.split("=", 1)
            d[k.strip()] = [float(x.replace("D", "E")) for x in v.split()]
    return d


def same(x, y):
    return x == y or (math.isnan(x) and math.isnan(y))


def main():
    ref, new = sys.argv[1], sys.argv[2]
    steps = sorted(d for d in os.listdir(ref) if d.startswith("step_"))
    if not steps:
        print("no steps in", ref)
        return 2
    bad, nbin, nval, ntxt = [], 0, 0, 0
    pairs = [("domain_gas.txt", os.path.join(ref, "domain_gas.txt"),
              os.path.join(new, "domain_gas.txt"))]
    for s in steps:
        for fn in sorted(os.listdir(os.path.join(ref, s))):
            pairs.append((f"{s}/{fn}", os.path.join(ref, s, fn), os.path.join(new, s, fn)))
    for name, a, b in pairs:
        if not os.path.exists(a):
            continue
        if not os.path.exists(b):
            print("missing", b)
            return 2
        if name.endswith("stats.txt"):
            pa, pb = parse(a), parse(b)
            for k, va in pa.items():
                vb = pb.get(k)
                if vb is None or len(va) != len(vb):
                    bad.append(f"{name}: key {k} missing or wrong length")
                    continue
                for x, y in zip(va, vb):
                    nval += 1
                    if not same(x, y):
                        bad.append(f"{name}: {k} ref={x!r} new={y!r}")
        else:
            if name.endswith(".bin"):
                nbin += 1
            else:
                ntxt += 1
            if open(a, "rb").read() != open(b, "rb").read():
                bad.append(f"{name}: not byte-identical")
    print(f"{ref}: {len(steps)} steps, {nbin} binary files, {ntxt} text files, "
          f"{nval} stat values compared")
    for b in bad[:30]:
        print("  MISMATCH", b)
    if bad:
        print(f"FAILED ({len(bad)} mismatches)")
        return 1
    print("IDENTICAL")
    return 0


if __name__ == "__main__":
    sys.exit(main())
