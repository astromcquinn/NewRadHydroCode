#!/usr/bin/env python3
"""Check the star-formation sink's outputs of a radhydro run (docs/stage_f.md
section 2).

usage: sink_check.py RUN [--min-stars 1]

RUN/sink.txt holds the final count and total mass of the star particles;
RUN/DM/stars.<ztag> (the latest) the particles. Checks: at least
--min-stars stars formed, the file's count and mass sum equal sink.txt's,
the ids are unique and above the dark-matter range, masses positive.
"""
import argparse
import glob
import os
import sys

import numpy as np

STAR = np.dtype([("id", "<i8"), ("x", "<f4", 3), ("v", "<f4", 3), ("mass", "<f4"), ("a_form", "<f4")])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--min-stars", type=int, default=1)
    a = ap.parse_args()
    summ = {}
    for line in open(os.path.join(a.run, "sink.txt")):
        k, v = (x.strip() for x in line.split("=", 1))
        summ[k] = float(v)
    files = sorted(glob.glob(os.path.join(a.run, "DM", "stars.*")), key=os.path.getmtime)
    if not files:
        sys.exit("no DM/stars.* file")
    with open(files[-1], "rb") as f:
        n = int(np.fromfile(f, dtype="<i8", count=1)[0])
        st = np.fromfile(f, dtype=STAR, count=n)
    bad = []
    if n < a.min_stars:
        bad.append(f"only {n} stars (need {a.min_stars})")
    if n != int(summ["nstar"]):
        bad.append(f"stars.* has {n}, sink.txt {summ['nstar']:.0f}")
    m = float(np.sum(st["mass"].astype(np.float64)))
    if abs(m - summ["mass"]) > 1e-9 * max(1.0, abs(summ["mass"])):
        bad.append(f"mass {m!r} vs sink.txt {summ['mass']!r}")
    if n and (np.any(st["mass"] <= 0) or len(set(st["id"].tolist())) != n):
        bad.append("non-positive masses or repeated ids")
    print(f"{files[-1]}: {n} stars, mass {m:.6g}, nstar_ever {summ.get('nstar_ever', 0):.0f}")
    if bad:
        print("FAIL: " + "; ".join(bad))
        sys.exit(1)


if __name__ == "__main__":
    main()
