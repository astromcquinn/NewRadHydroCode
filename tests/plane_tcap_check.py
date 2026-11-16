#!/usr/bin/env python3
"""The temperature cap T_max (Config::T_max) in the plane-parallel full loop.

usage: plane_tcap_check.py RADHYDRO WORKDIR

A 16^3 Zel'dovich box of 0.25 Mpc/h from z = 20, plane_parallel = 1 with
2^3 domains from z_re = 9 (tests/plane_full_check.py's set-up), run to
z = 8.5 three times:
  * T_max = 0 (off) and T_max = 1e15 K (never reached) give byte-identical
    stats and fields at every step;
  * T_max = 1.2e4 K binds (the photoheated gas is at ~2e4 K): the hydro
    clamps cells (timing.txt "# cap_T" > 0), the RT temperature field is at
    most T_max (1e-4) at every output after z_re, where the uncapped run
    exceeds 1.3 T_max, and the box is still ionised (x_HII,V > 0.99).
  * T_max in a run whose last zdata1 output is below z = 4 is refused
    (the WHIM; T_max_low_z = 1 overrides).
"""
import filecmp
import glob
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from plane_full_check import FAILS, check, make_problem, run, steps  # noqa: E402


def temps(out):
    return [np.fromfile(p, dtype=np.float32) for p in sorted(glob.glob(os.path.join(out, "step_*", "temp.bin")))]


def main():
    exe, work = sys.argv[1], sys.argv[2]
    os.makedirs(work, exist_ok=True)
    N, z_re, cap = 16, 9.0, 1.2e4
    prob = os.path.join(work, "prob")
    make_problem(prob, N, 0.25, 20, "15 12 10 9 8.5", False,
                 {"plane_ndom": 2, "plane_faces": "+x +y", "plane_gamma12": 3, "plane_z_re": z_re,
                  "c_reduced_schedule": "0.1 0.99 0.01", "cool_suppress_overdensity": 0})
    off = run(exe, prob, os.path.join(work, "off"), [])
    huge = run(exe, prob, os.path.join(work, "huge"), ["T_max=1e15"])
    capd = run(exe, prob, os.path.join(work, "cap"), [f"T_max={cap}"])

    a = sorted(glob.glob(os.path.join(off, "step_*", "*")))
    b = sorted(glob.glob(os.path.join(huge, "step_*", "*")))
    same = len(a) == len(b) > 0 and all(
        os.path.relpath(x, off) == os.path.relpath(y, huge) and filecmp.cmp(x, y, shallow=False) for x, y in zip(a, b))
    check(same, f"T_max = 1e15 byte-identical to T_max = 0 ({len(a)} files)")

    ncap = [int(l.split()[2]) for l in open(os.path.join(capd, "timing.txt")) if l.startswith("# cap_T")]
    check(len(ncap) == 1 and ncap[0] > 0, f"the hydro clamps cells: cap_T = {ncap}")
    st_off, st_cap = steps(off), steps(capd)
    to, tc = temps(off), temps(capd)
    check(len(to) == len(tc) and len(tc) >= 2, f"the same outputs ({len(to)}, {len(tc)})")
    def tmax_after_zre(out, st):
        return max(float(np.fromfile(os.path.join(out, f"step_{i:04d}", "temp.bin"), dtype=np.float32).max())
                   for i, s in enumerate(st) if s["redshift"][0] < z_re - 1e-6 and
                   os.path.isfile(os.path.join(out, f"step_{i:04d}", "temp.bin")))
    tmax_o, tmax_c = tmax_after_zre(off, st_off), tmax_after_zre(capd, st_cap)
    check(tmax_o > 1.3 * cap, f"uncapped max T after z_re {tmax_o:.4g} K > 1.3 T_max")
    check(tmax_c <= cap * (1 + 1e-4), f"capped max T after z_re {tmax_c:.6g} K <= T_max = {cap:g}")
    check(st_cap[-1]["xHII_V"][0] > 0.99, f"ionised at the end: x_HII,V = {st_cap[-1]['xHII_V'][0]:.5f}")

    # a cap in a run that ends below z = 4 (it would remove the WHIM) is
    # refused unless T_max_low_z = 1; the refusal comes before the first step
    low = os.path.join(work, "lowz")
    os.makedirs(low, exist_ok=True)
    cmd = [exe, prob, low, "--quiet", "--set", f"T_max={cap}", "--set", "zdata1=15 12 3.5"]
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=600)
    check(r.returncode != 0 and "WHIM" in r.stdout, "T_max refused in a run ending below z = 4")

    if FAILS:
        print(f"{len(FAILS)} FAILED")
        sys.exit(1)
    print("all passed")


if __name__ == "__main__":
    main()
