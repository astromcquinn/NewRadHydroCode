#!/usr/bin/env python3
"""Conservation of a refined hydro_driver run against its uniform twin.

usage: amr_conservation_check.py REF_DIR NEW_DIR [--rtol 1e-9] [--levels n]

REF_DIR and NEW_DIR are hydro_driver output directories (step_NNNN/stats.txt
with the conserved totals of the level-0 state). With refinement
(amrex_max_level > 0, docs/stage_f.md) the level-0 fields differ from the
uniform run, but the box is periodic and the fine levels are refluxed and
averaged down, so the total mass, momentum and energy of every step must
agree with the uniform run to round-off, and the species totals too.
--levels n additionally requires the driver's stdout (NEW_DIR.stdout, as
run_amrex.cmake writes it) to report that many levels in use.
"""
import argparse
import glob
import os
import sys


def totals(path):
    t = {}
    for line in open(path):
        if "=" not in line:
            continue
        k, v = (x.strip() for x in line.split("=", 1))
        if k in ("mass", "energy", "energy_grid"):
            t[k] = [float(v)]
        elif k in ("momentum", "w_total"):
            t[k] = [float(x) for x in v.split()]
    return t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ref")
    ap.add_argument("new")
    ap.add_argument("--rtol", type=float, default=1e-9)
    ap.add_argument("--levels", type=int, default=0)
    a = ap.parse_args()
    steps = sorted(os.path.basename(d) for d in glob.glob(os.path.join(a.new, "step_*")))
    if not steps:
        sys.exit(f"no steps in {a.new}")
    bad = 0
    for s in steps:
        r = totals(os.path.join(a.ref, s, "stats.txt"))
        n = totals(os.path.join(a.new, s, "stats.txt"))
        for k in ("mass", "momentum", "energy", "w_total"):
            if k not in r or k not in n:
                continue
            scale = max(abs(x) for x in r[k]) or 1.0
            if k == "momentum":
                scale = max(scale, abs(r["mass"][0]) * 1e-3)  # momenta sum to ~0: relative to a velocity scale
            for x, y in zip(r[k], n[k]):
                if abs(x - y) > a.rtol * scale:
                    print(f"{s}: {k} {x!r} vs {y!r} (diff {abs(x - y):.3e}, allowed {a.rtol * scale:.3e})")
                    bad += 1
    if a.levels > 0:
        out = a.new.rstrip("/") + ".stdout"
        txt = open(out).read() if os.path.exists(out) else ""
        if f"levels {a.levels}" not in txt:
            print(f"{out}: expected 'levels {a.levels}' in the AMReX hydro timers")
            bad += 1
    print(f"{len(steps)} steps compared, {bad} violations")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
