#!/usr/bin/env python3
"""photo_cell's photon budget (B17, photon_conserve; docs/rt_multicell.md).

usage: photon_conserve_check.py RT_DRIVER LOGNORMAL_PROBLEM_DIR WORKDIR

Runs rt_driver on lognormal_16src (legacy_rt = legacy_physics = 0,
rt_substep_carry = 1, 20 steps) and compares the ionised hydrogen at the
end, x_HII,M = sum rho x_HII / sum rho, with the counted net ionisations
per H atom of rt_stats.txt (photo-ionisations - recombinations +
collisional ionisations of HI, cumulative over the run). The difference
("leak") is hydrogen ionised by photons the cells did not receive.
Checks (exit 1 on failure):
  * the original limiter (photon_conserve = 0) leaks at k = 1 (> 1e-3:
    the test resolves B17) and more at rt_cells_per_substep = 4;
  * photon_conserve = 1 closes the budget to 3e-4 at k = 1 and k = 4.
"""
import os
import shutil
import subprocess
import sys

import numpy as np

NSTEPS = 20


def kv(fn):
    d = {}
    for line in open(fn):
        if "=" in line:
            k, v = line.split("=", 1)
            try:
                d[k.strip()] = [float(x) for x in v.split()]
            except ValueError:
                pass
    return d


def leak(driver, prob, out, k, pc):
    if os.path.isdir(out):
        shutil.rmtree(out)
    with open(out + ".stdout", "w") as log:
        subprocess.run([driver, prob, out, "--nsteps", str(NSTEPS), "--set", "legacy_rt=0",
                        "--set", "legacy_physics=0", "--set", "rt_substep_carry=1",
                        "--set", f"rt_cells_per_substep={k}", "--set", f"photon_conserve={pc}"],
                       check=True, stdout=log, stderr=subprocess.STDOUT)
    sd = os.path.join(out, f"step_{NSTEPS:04d}")
    fHI = np.fromfile(os.path.join(sd, "fHI.bin"), dtype=np.float32).astype(np.float64)
    rho = np.fromfile(os.path.join(sd, "rhog.bin"), dtype=np.float32).astype(np.float64)
    xM = float(np.sum(rho * (1 - fHI)) / np.sum(rho))
    st = kv(os.path.join(sd, "rt_stats.txt"))
    net = st["nion"][0] - st["nrec"][0] + st["ncol"][0]
    return xM, xM - net


def main():
    drv, prob, work = sys.argv[1], sys.argv[2], sys.argv[3]
    os.makedirs(work, exist_ok=True)
    ok = True
    res = {}
    for k, pc in ((1, 0), (4, 0), (1, 1), (4, 1)):
        xM, lk = leak(drv, prob, os.path.join(work, f"k{k}_pc{pc}"), k, pc)
        res[(k, pc)] = lk
        print(f" k = {k}, photon_conserve = {pc}: x_HII,M = {xM:.5f}, leak = {lk:+.2e} "
              f"({lk / xM:+.2%} of x_HII,M)")
    if not res[(1, 0)] > 1e-3:
        print("  FAIL: the original limiter does not leak here (B17 not resolved)")
        ok = False
    if not res[(4, 0)] > res[(1, 0)]:
        print("  FAIL: the leak does not grow with k")
        ok = False
    for k in (1, 4):
        if abs(res[(k, 1)]) > 3e-4:
            print(f"  FAIL: photon_conserve = 1 leaks {res[(k, 1)]:.2e} at k = {k}")
            ok = False
    print("PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
