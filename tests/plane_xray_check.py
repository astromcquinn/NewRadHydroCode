#!/usr/bin/env python3
"""X-ray preheating in the plane-parallel mode before the turn-on, and
plane_xray_substeps (docs/plane_parallel.md, docs/xray_heating.md).

usage: plane_xray_check.py RADHYDRO WORKDIR

A 16^3 Zel'dovich box of 0.25 Mpc/h from z = 20, plane_parallel = 1 with
2^3 domains from z_re = 9 (tests/plane_full_check.py's set-up), and
xray_heating = 1 with the production abundance matching
(data/uv_lf/am_params_st_fof.txt) at log L_X/SFR = 40.5 (a few hundred K
by z_re; the C1 1024^3 run uses 39.9, ~100 K at z = 7), run to z = 8.5:
  * the X-rays heat the gas before z_re (mass-weighted T > 30 K at z_re);
  * plane_xray_substeps = 10 takes 10 RT substeps per hydro step before z_re
    (stats nstep_rt), the default the light crossings (0: off, 50-1200 here);
  * the two agree on the mass-weighted T and x_HII of the fields at the
    z_re landing step to 3 per cent, and on T_V, T_M and x_HII,V at the end
    to 2 per cent. The error is first order in the substep length, ~0.35 per
    cent per Myr: with this box's 30-100 Myr hydro steps 1, 3, 10 and 30
    substeps give -16.6, -5.6, -1.6 and -0.4 per cent at z_re
    (docs/plane_parallel.md).
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, HERE)
from plane_full_check import FAILS, check, make_problem, run, steps  # noqa: E402


def am_params():
    d = {}
    for line in open(os.path.join(ROOT, "data", "uv_lf", "am_params_st_fof.txt")):
        line = line.split("#")[0]
        if "=" in line:
            k, v = line.split("=", 1)
            d[k.strip()] = v.strip()
    return d


def main():
    exe, work = sys.argv[1], sys.argv[2]
    os.makedirs(work, exist_ok=True)
    N, z_re = 16, 9.0
    base = {"plane_ndom": 2, "plane_faces": "+x +y", "plane_gamma12": 3, "plane_z_re": z_re,
            "c_reduced_schedule": "0.1 0.99 0.01", "cool_suppress_overdensity": 0}
    p0 = os.path.join(work, "prob_nox")
    make_problem(p0, N, 0.25, 20, "15 12 10 9 8.5", False, base)
    px = os.path.join(work, "prob_x")
    make_problem(px, N, 0.25, 20, "15 12 10 9 8.5", False,
                 {**base, **am_params(), "xray_heating": 1, "xray_log_lx_sfr": 40.5})
    nox = steps(run(exe, p0, os.path.join(work, "nox"), []))
    full = steps(run(exe, px, os.path.join(work, "full"), []))
    one = steps(run(exe, px, os.path.join(work, "ten"), ["plane_xray_substeps=10"]))

    def z(st):
        return np.array([s["redshift"][0] for s in st])

    zf, zo = z(full), z(one)
    nro = np.array([s["nstep_rt"][0] for s in one])
    nrf = np.array([s["nstep_rt"][0] for s in full])
    lo = [i for i, zz in enumerate(zo) if abs(zz - z_re) < 1e-9]
    lf = [i for i, zz in enumerate(zf) if abs(zz - z_re) < 1e-9]
    check(len(lo) == 1 and len(lf) == 1, "both runs land on z_re")
    io, i_f = lo[0], lf[0]
    check(bool(np.all(nro[1:io + 1] == 10)), "10 RT substeps per hydro step before z_re (the landing step included)")
    check(bool(np.all(nrf[1:i_f + 1] > 10)), f"the default takes the light crossings before z_re ({int(nrf[1:i_f + 1].min())}+)")

    # the mass-weighted T and x_HI of the fields at the z_re landing step (stats
    # T_V is dominated by a few near-empty cells, T_M is not written before z_re)
    def fields_mw(out, i):
        d = os.path.join(out, f"step_{i:04d}")
        T = np.fromfile(os.path.join(d, "temp.bin"), dtype=np.float32).astype(np.float64)
        r = np.fromfile(os.path.join(d, "rhog.bin"), dtype=np.float32).astype(np.float64)
        f = np.fromfile(os.path.join(d, "fHI.bin"), dtype=np.float32).astype(np.float64)
        return np.sum(r * T) / r.sum(), np.sum(r * f) / r.sum()
    Tf, xf = fields_mw(os.path.join(work, "full"), i_f)
    To, xo = fields_mw(os.path.join(work, "ten"), io)
    check(Tf > 2 * 15, f"X-rays heat before z_re: mass-weighted T(z_re) = {Tf:.1f} K")
    check(abs(To / Tf - 1) < 0.03, f"mass-weighted T at z_re: {To:.2f} vs {Tf:.2f} K ({To / Tf - 1:+.4f}, < 3 per cent)")
    check(abs((1 - xo) / (1 - xf) - 1) < 0.03, f"mass-weighted x_HII at z_re: {1 - xo:.3e} vs {1 - xf:.3e} (< 3 per cent)")
    for key in ("T_V", "T_M", "xHII_V"):
        r = abs(one[-1][key][0] / full[-1][key][0] - 1)
        check(r < 0.02, f"{key} at the end (z = {zo[-1]:.3f} and {zf[-1]:.3f}) to {r:.4f} < 0.02")

    if FAILS:
        print(f"{len(FAILS)} FAILED")
        sys.exit(1)
    print("all passed")


if __name__ == "__main__":
    main()
