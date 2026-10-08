#!/usr/bin/env python3
"""Strict comparison of a full radhydro run with a stored Fortran reference.

usage: compare_full_reference.py REF_DIR NEW_DIR

REF_DIR is reference/harness_full/outputs/<case> (store_outputs.py) or a
complete run directory:
  step_NNNN/stats.txt    every value must be the identical number (NaN = NaN)
  step_NNNN/sha256.txt   every listed .bin of NEW must have this SHA-256
  step_NNNN/*.bin        byte-identical
  Log/hydro.log          byte-identical (summarize_gas unit-100 records)
  Pk/*, DM/*             byte-identical (Pk) / SHA-256 (DM, in DM/sha256.txt)
  Gas/*, RT/*            byte-identical (write_data_1/2 with the F1 fix; only
                         when REF has them, i.e. a live Fortran run)
  Ray/*                  the C++ file is a prefix of the Fortran one (which
                         holds all Nray_max slots) and the Fortran remainder
                         is the untouched tail: zeros (ray) or the identity
                         (ray_index); htoc_ray byte-identical
NEW_DIR must have the same steps as REF_DIR (no more, no fewer).
Exit 0 if identical, 1 if anything differs, 2 on missing files.
"""
import glob
import hashlib
import math
import os
import sys


def parse(fn):
    d = {}
    for line in open(fn):
        if "=" in line:
            k, v = line.split("=", 1)
            d[k.strip()] = [float(x.replace("D", "e")) for x in v.split()]
    return d


def same(a, b):
    return len(a) == len(b) and all(x == y or (math.isnan(x) and math.isnan(y))
                                    for x, y in zip(a, b))


def sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


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

    nsteps_new = len(glob.glob(os.path.join(new, "step_*")))
    if nsteps_new != len(steps):
        bad.append(f"number of steps: ref {len(steps)}, new {nsteps_new}")
    for sd in steps:
        step = os.path.basename(sd)
        nd = os.path.join(new, step)
        a, b = parse(os.path.join(sd, "stats.txt")), parse(need(os.path.join(nd, "stats.txt")))
        for k, v in a.items():
            nval += 1
            if k not in b or not same(v, b[k]):
                bad.append(f"{step}/stats.txt {k}: ref {v!r} new {b.get(k)!r}")
        hfile = os.path.join(sd, "sha256.txt")
        if os.path.exists(hfile):
            for line in open(hfile):
                name, h = line.split()
                nhash += 1
                if sha(need(os.path.join(nd, name))) != h:
                    bad.append(f"{step}/{name}: SHA-256 differs")
        for f in glob.glob(os.path.join(sd, "*.bin")):
            nfile += 1
            if open(f, "rb").read() != open(need(os.path.join(nd, os.path.basename(f))), "rb").read():
                bad.append(f"{step}/{os.path.basename(f)}: not byte-identical")
    for sub in ("Log/hydro.log",):
        f = os.path.join(ref, sub)
        if os.path.exists(f):
            nfile += 1
            if open(f, "rb").read() != open(need(os.path.join(new, sub)), "rb").read():
                bad.append(f"{sub}: differs")
    for f in sorted(glob.glob(os.path.join(ref, "Pk", "*"))):
        nfile += 1
        g = need(os.path.join(new, "Pk", os.path.basename(f)))
        if open(f, "rb").read() != open(g, "rb").read():
            bad.append(f"Pk/{os.path.basename(f)}: differs")
    dmh = os.path.join(ref, "DM", "sha256.txt")
    dms = ([l.split() for l in open(dmh)] if os.path.exists(dmh) else
           [(os.path.basename(f), sha(f)) for f in glob.glob(os.path.join(ref, "DM", "xv_dm.*"))])
    for name, h in dms:
        nhash += 1
        if sha(need(os.path.join(new, "DM", name))) != h:
            bad.append(f"DM/{name}: SHA-256 differs")
    for sub in ("Gas", "RT"):
        for f in sorted(glob.glob(os.path.join(ref, sub, "*"))):
            nfile += 1
            g = need(os.path.join(new, sub, os.path.basename(f)))
            if open(f, "rb").read() != open(g, "rb").read():
                bad.append(f"{sub}/{os.path.basename(f)}: differs")
    for f in sorted(glob.glob(os.path.join(ref, "Ray", "*"))):
        nfile += 1
        name = os.path.basename(f)
        a = open(f, "rb").read()
        b = open(need(os.path.join(new, "Ray", name)), "rb").read()
        if name.startswith("htoc_ray"):
            ok = a == b
        else:
            ok = len(b) <= len(a) and a[:len(b)] == b
            tail = a[len(b):]
            if ok and name.startswith("ray."):
                ok = tail.count(0) == len(tail)
            elif ok and name.startswith("ray_index."):
                import struct
                n0 = len(b) // 8
                vals = struct.unpack(f"<{len(tail) // 8}q", tail)
                ok = all(v == n0 + 1 + i for i, v in enumerate(vals))
        if not ok:
            bad.append(f"Ray/{name}: differs")
    print(f"{ref}: {len(steps)} steps, {nval} stat values, {nhash} hashes, "
          f"{nfile} files compared")
    for x in bad[:20]:
        print("  DIFF", x)
    print("IDENTICAL" if not bad else f"DIFFERENT ({len(bad)} items)")
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
