#!/usr/bin/env python3
"""Copy the stored subset of harness runs into reference/harness_hydro/outputs.

usage: store_reference.py [runs=reference/harness_hydro/runs] [outputs=.../outputs]

Kept per problem: stats.txt and gas_summary.txt of every step,
domain_gas.txt, stdout.log, and the binary fields (u_gas, w_gas, v_gas) of
the LAST step only.  sod_y and sod_z keep no binaries (their fields are
covered by the live Fortran-vs-C++ comparison when gfortran is available;
the stored stats include all conserved totals).
"""
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
NO_FIELDS = {"sod_y", "sod_z"}


def main():
    runs = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "runs")
    outs = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "outputs")
    for name in sorted(os.listdir(runs)):
        src = os.path.join(runs, name)
        if not os.path.isdir(src):
            continue
        dst = os.path.join(outs, name)
        shutil.rmtree(dst, ignore_errors=True)
        os.makedirs(dst)
        for fn in ("domain_gas.txt", "stdout.log"):
            if os.path.isfile(os.path.join(src, fn)):
                shutil.copy2(os.path.join(src, fn), dst)
        steps = sorted(d for d in os.listdir(src) if d.startswith("step_"))
        for s in steps:
            os.makedirs(os.path.join(dst, s))
            for fn in os.listdir(os.path.join(src, s)):
                keep = fn.endswith(".txt") or (s == steps[-1] and name not in NO_FIELDS)
                if keep:
                    shutil.copy2(os.path.join(src, s, fn), os.path.join(dst, s, fn))
        print(f"{name}: {len(steps)} steps stored")


if __name__ == "__main__":
    main()
