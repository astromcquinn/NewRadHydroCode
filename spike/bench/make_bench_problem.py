#!/usr/bin/env python3
"""Benchmark RT problems for the framework spike.

    python3 spike/bench/make_bench_problem.py --N 64 --nsrc 64 --nsteps 4 \
        --nray-max-factor 16 --out /path/bench_N64

Log-normal density field and sources in the densest cells, exactly as
reference/harness/make_test_problem.py builds lognormal_16src, but with a
free grid size and source count. legacy_rt = 0 is written into params.txt,
since the spike implements only the sorted-array (legacy_rt = 0) algorithm.
"""
import argparse
import os
import sys

import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, "..", "..", "reference", "harness"))
import make_test_problem as mtp  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--N", type=int, default=64)
    ap.add_argument("--nsrc", type=int, default=64)
    ap.add_argument("--nsteps", type=int, default=4)
    ap.add_argument("--nray-max-factor", type=int, default=16)
    ap.add_argument("--luv", type=float, default=1.0e28)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    n = a.N
    d = mtp.lognormal_delta(n)
    order = np.argsort(d.ravel(order="F"), kind="stable")[::-1][: a.nsrc]
    sources = []
    for r, idx in enumerate(order):
        i, j, k = np.unravel_index(idx, (n, n, n), order="F")
        luv = a.luv * (1.0 + 0.1 * (r % 16))
        sources.append((int(i) + 1, int(j) + 1, int(k) + 1, luv, luv / 2.0e28))
    gas, temp_min = mtp.make_gas(d, 10.0, n)
    note = f"Spike benchmark: log-normal field, N={n}, {a.nsrc} sources in the densest cells.\n"
    mtp.write_problem(a.out, n, gas, sources, 10.0, 1e-4, a.nsteps, note, a.nray_max_factor)
    with open(os.path.join(a.out, "params.txt"), "a") as f:
        f.write("legacy_rt = 0\n")
    print(f"written {a.out}")


if __name__ == "__main__":
    main()
