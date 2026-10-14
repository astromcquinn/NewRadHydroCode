#!/usr/bin/env python3
"""The I-front of a radhydro run from its fHI fields: the radius of the
ionised volume, R = (3 V/4 pi)^(1/3) with V = sum(1 - x_HI) cells, per step,
in level-0 cells (the refined runs' level-0 field is the composite).

usage: ifront_full.py RUN [RUN ...] [--N n] [--check REF BASE NEW tol]

Prints R per step for every run. --check REF BASE NEW tol: the NEW run
(refined, base grid N) must be at least as close to REF (uniform 2N, in
units of the base cell) as BASE (uniform N) is, within tol base cells, at
every common step: |R_new - R_ref| <= |R_base - R_ref| + tol. Exit 1 if not.
"""
import argparse
import glob
import os
import sys

import numpy as np


def radii(run, N):
    out = {}
    for d in sorted(glob.glob(os.path.join(run, "step_*"))):
        f = os.path.join(d, "fHI.bin")
        if not os.path.exists(f):
            continue
        x = np.fromfile(f, dtype=np.float32)
        if x.size != N ** 3:
            raise SystemExit(f"{f}: {x.size} values, expected {N ** 3}")
        V = float(np.sum(1.0 - x.astype(np.float64)))
        if V >= 0.999 * N ** 3:
            continue  # before the RT starts the field is unset (x_HI = 0 everywhere)
        out[os.path.basename(d)] = (3 * V / (4 * np.pi)) ** (1 / 3)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="*")
    ap.add_argument("--N", type=int, default=32)
    ap.add_argument("--check", nargs=4, metavar=("REF", "BASE", "NEW", "TOL"))
    a = ap.parse_args()
    for r in a.runs:
        R = radii(r, a.N)
        print(r, " ".join(f"{k[5:]}:{v:.3f}" for k, v in R.items()))
    if a.check:
        ref, base, new, tol = a.check
        tol = float(tol)
        Rr = radii(ref, 2 * a.N)
        Rb = radii(base, a.N)
        Rn = radii(new, a.N)
        bad = 0
        for k in sorted(set(Rb) & set(Rn) & set(Rr)):
            rr = Rr[k] / 2  # the 2N run's radius in base cells
            eb, en = abs(Rb[k] - rr), abs(Rn[k] - rr)
            flag = "" if en <= eb + tol else "  <-- worse than the base run"
            if flag:
                bad += 1
            print(f"{k}: R_ref {rr:.3f} base {Rb[k]:.3f} (err {eb:.3f}) refined {Rn[k]:.3f} (err {en:.3f}){flag}")
        print("PASS" if not bad else f"FAIL: {bad} steps")
        sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
