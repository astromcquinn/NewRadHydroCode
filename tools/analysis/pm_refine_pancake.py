#!/usr/bin/env python3
"""Zel'dovich pancake with gas and dark matter on the refined PM mesh
(docs/pm_refine.md; the problem of docs/hydro_tests.md section 6).

Runs the hydro-suite driver on the pancake with pm_refine = 1 (Ndm = N, the
original) and pm_refine = 2 (Ndm = 2N on the 2N mesh), for both hydro
schemes (legacy: calc_kick; unsplit: calc_kick_unsplit), and measures the
errors against the exact Zel'dovich solution before the caustic (z_c = 5):
  gas  L1 of rho and v (analyze.py analyze_pancake);
  DM   x and v of every particle along the wave, max error relative to
       the displacement / velocity amplitude, and the error of the fitted
       sin(kq) mode.

usage: pm_refine_pancake.py RUNDIR [--N 32 64] [--driver build/hydro_suite_driver]
                            [--jobs 4] [--check] [--json out.json]
--check: pm_refine = 2 must be within 20% (+ 2e-3 absolute) of pm_refine = 1
in every measure, at every N, scheme and redshift (exit status 1 if not).
"""
import argparse
import json
import os
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SUITE = os.path.join(HERE, "..", "..", "tests", "hydro_suite")
sys.path.insert(0, SUITE)
sys.path.insert(0, os.path.join(HERE, "..", "ics"))
import analyze as an  # noqa: E402
import cosmo_gen as cg  # noqa: E402
import hs_common as hc  # noqa: E402
import pmcosmo as pc  # noqa: E402

Z_OUT = (20, 10, 7, 5.5)


def dm_errors(d, meta, Ndm):
    xv = np.fromfile(os.path.join(d, "xv_dm.bin"), "<f4").reshape(6, -1, order="F").astype(float)
    info = hc.read_state(d, meta["N"])[3]
    a = info["scalefactor"]
    G = pc.Growth(1.0)
    a_c = 1 / (1 + meta["z_c"])
    A = G(a) / G(a_c)
    Ad = G.dtau(a) / G(a_c)
    k = 2 * np.pi / Ndm
    c = Ndm / 2.0
    q1 = np.arange(Ndm) + 0.5
    q = np.meshgrid(q1, q1, q1, indexing="ij")[0].ravel(order="F")
    s = np.sin(k * (q - c)) / k
    dx = (xv[0] - (q - A * s) + Ndm / 2) % Ndm - Ndm / 2
    dv = xv[3] + Ad * s
    disp = (xv[0] - q + Ndm / 2) % Ndm - Ndm / 2
    Am = -(disp * s).sum() / (s * s).sum()
    Vm = -(xv[3] * s).sum() / (s * s).sum()
    return {"max_x": float(np.abs(dx).max() / (A / k)), "max_v": float(np.abs(dv).max() / (Ad / k)),
            "mode_x": float(abs(Am / A - 1)), "mode_v": float(abs(Vm / Ad - 1))}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("rundir")
    ap.add_argument("--N", type=int, nargs="+", default=[32, 64])
    ap.add_argument("--schemes", nargs="+", default=["legacy", "unsplit"])
    ap.add_argument("--driver", default=None)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--json", default=None)
    a = ap.parse_args()
    rundir = os.path.abspath(a.rundir)
    cases = []
    for N in a.N:
        for r in (1, 2):
            pdir = os.path.join(rundir, "problems", f"pancake_N{N}_r{r}")
            meta = cg.gen_pancake(pdir, N=N, Ndm=r * N, T_init=1000.0, z_out=Z_OUT)
            for sch in a.schemes:
                cases.append((N, r, sch, pdir, os.path.join(rundir, f"pancake_N{N}_r{r}_{sch}"),
                              meta))

    def run(c):
        N, r, sch, pdir, odir, meta = c
        sets = [f"pm_refine={r}"] + ([f"hydro_scheme={sch}"] if sch != "legacy" else [])
        if r == 1:
            sets = sets[1:]  # the default: no pm_refine line at all
        w = hc.run_driver(pdir, odir, sets, driver=a.driver)
        return w

    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        walls = list(ex.map(run, cases))

    res = {}
    print(" N  scheme   r   z     gas L1 rho  gas L1 v   DM max x   DM max v   DM mode x  DM mode v  wall")
    for c, w in zip(cases, walls):
        N, r, sch, pdir, odir, meta = c
        g = an.analyze_pancake(odir, meta, None)
        outs = hc.outputs(odir)[1:]
        for o, d in zip(g["outputs"], outs):
            zd = o.get("zeldovich")
            if not zd:
                continue
            dm = dm_errors(d, meta, r * N)
            key = f"N{N}_{sch}_z{o['z']:.2f}"
            res.setdefault(key, {})[f"r{r}"] = {"l1_rho": zd["l1_rho"], "l1_v": zd["l1_v"], **dm}
            print(f"{N:3d}  {sch:8s} {r:d}  {o['z']:5.2f}  {zd['l1_rho']:9.3e}  {zd['l1_v']:9.3e}  "
                  f"{dm['max_x']:9.3e}  {dm['max_v']:9.3e}  {dm['mode_x']:9.3e}  "
                  f"{dm['mode_v']:9.3e}  {w:6.1f} s")
    fails = []
    for key, rr in res.items():
        if "r1" not in rr or "r2" not in rr:
            continue
        for m, v1 in rr["r1"].items():
            v2 = rr["r2"][m]
            if not (v2 <= 1.2 * v1 + 2e-3):
                fails.append(f"{key} {m}: r2 {v2:.3e} vs r1 {v1:.3e}")
    if a.json:
        with open(a.json, "w") as f:
            json.dump(res, f, indent=1)
    if a.check:
        for m in fails:
            print("FAIL:", m)
        print("RESULT:", "PASS" if not fails else "FAIL")
        return 1 if fails else 0
    return 0


if __name__ == "__main__":
    sys.exit(main())
