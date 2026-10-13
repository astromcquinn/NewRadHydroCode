#!/usr/bin/env python3
"""A random3d hydro problem (reference/harness_hydro/make_hydro_problems.py)
at any N for the hydro timing and scaling runs (hydro_driver,
hydro_scheme = unsplit; docs/stage_e_hydro.md).

usage: make_hydro_problem.py OUTDIR [--N 256] [--steps 2] [--da 0.002]
"""
import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "reference", "harness_hydro"))
import make_hydro_problems as m  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("outdir")
    ap.add_argument("--N", type=int, default=256)
    ap.add_argument("--steps", type=int, default=2)
    ap.add_argument("--da", type=float, default=0.002)
    a = ap.parse_args()
    rho, vel, eth, vg, w, readme = m.random3d(N=a.N)
    root, name = os.path.split(os.path.abspath(a.outdir))
    m.write_problem(root, name, {"z_init": 10, "nsteps": a.steps, "timestep": "fixed", "da_sim": a.da,
                                 "hydro_scheme": "unsplit"},
                    rho, vel, eth, vg=None, w=w, readme=readme, force=True)


if __name__ == "__main__":
    main()
