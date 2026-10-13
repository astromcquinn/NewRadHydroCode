#!/usr/bin/env python3
"""Fields of a radhydro run on the AMReX backend against the serial run of
the same step (tools/hyak/amrex_scaling.slurm): for every field of the last
common step_* directory the rms and the maximum relative difference, the
fraction of cells beyond --rtol, and the stats.txt lines that differ.

usage: amrex_scaling_compare.py SERIAL_DIR AMREX_DIR [--rtol 1e-3]

Round-off between the two codes (the deposit and merge order,
docs/stage_d.md section 2.3) shows as rms differences of 1e-6 to 1e-4 and
single cells (I-front, removal threshold) at 1e-3; anything larger is a
real difference. Exit status 1 when a field's rms relative difference
exceeds --rtol.
"""
import argparse
import os
import sys

import numpy as np


def steps(d):
    return sorted(x for x in os.listdir(d) if x.startswith("step_") and os.path.isdir(os.path.join(d, x)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("serial")
    ap.add_argument("amrex")
    ap.add_argument("--rtol", type=float, default=1e-3, help="rms tolerance and per-cell threshold")
    a = ap.parse_args()
    common = sorted(set(steps(a.serial)) & set(steps(a.amrex)))
    if not common:
        sys.exit("no common step directory")
    last = common[-1]
    print(f"step {last}: field  rms_rel  max_rel  frac_cells>{a.rtol:g}  (n cells)")
    bad = []
    sd, ad = os.path.join(a.serial, last), os.path.join(a.amrex, last)
    for f in sorted(x for x in os.listdir(sd) if x.endswith(".bin")):
        if not os.path.exists(os.path.join(ad, f)):
            print(f"  {f}: missing in {ad}")
            continue
        x = np.fromfile(os.path.join(sd, f), dtype="<f4").astype(float)
        y = np.fromfile(os.path.join(ad, f), dtype="<f4").astype(float)
        if x.size != y.size:
            print(f"  {f}: sizes differ {x.size} vs {y.size}")
            bad.append(f)
            continue
        scale = np.sqrt(np.mean(x ** 2))
        if scale == 0:
            rms = np.sqrt(np.mean((y - x) ** 2))
            rel = np.abs(y - x)
        else:
            rms = np.sqrt(np.mean((y - x) ** 2)) / scale
            rel = np.abs(y - x) / np.maximum(np.abs(x), 1e-6 * np.max(np.abs(x)))
        frac = np.mean(rel > a.rtol)
        print(f"  {f[:-4]:8s} {rms:9.3e} {rel.max():9.3e} {frac:9.3e}  ({x.size})")
        if rms > a.rtol:
            bad.append(f)
    # stats.txt: the lines that differ
    try:
        s1 = open(os.path.join(sd, "stats.txt")).read().splitlines()
        s2 = open(os.path.join(ad, "stats.txt")).read().splitlines()
        diff = [(u, v) for u, v in zip(s1, s2) if u != v]
        print(f"stats.txt: {len(diff)} of {len(s1)} lines differ")
        for u, v in diff[:20]:
            print(f"  serial: {u}\n  amrex:  {v}")
    except OSError as e:
        print(f"stats.txt: {e}")
    print("PASS" if not bad else "FAIL: " + " ".join(bad))
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
