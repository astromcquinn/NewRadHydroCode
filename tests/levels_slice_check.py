#!/usr/bin/env python3
"""fields_levels = 1 (radhydro): the plane at the finest resolution.

usage: levels_slice_check.py WORK LAUNCH FLAGS NPFLAG NP DRIVER SRC

Makes the one-source problem of tests/problems_full/make_ifront.py at 16^3,
runs it with two refined levels over fixed nested boxes and fields_levels =
1 under the MPI launcher, and checks the slice of the last step: its size
(16 * 4), the level of every pixel (the boxes, cut by the plane), the
level-0 pixels equal to the base field fHI.bin in that plane, every value
finite and the level-1/2 pixels within [0, 1] for the fractions.
"""
import os
import subprocess
import sys

import numpy as np

N = 16
BOX1 = (4, 4, 4, 11, 11, 11)  # base cells, inclusive (level 1)
BOX2 = (6, 6, 6, 9, 9, 9)     # level 2
INDEX = 7                      # the slice: base cell 7 along z, inside both boxes


def main():
    work, launch, flags, npflag, nps, drv, src = sys.argv[1:8]
    npr = int(nps)
    flags = [f for f in flags.replace(";", " ").split() if f]
    if launch.endswith("srun") and not any("cpus-per-task" in f for f in flags):
        flags.append("--cpus-per-task=" + os.environ.get("OMP_NUM_THREADS", "1"))
    prob = os.path.join(work, "prob")
    out = os.path.join(work, "run")
    subprocess.check_call([sys.executable, os.path.join(src, "tests/problems_full/make_ifront.py"), prob,
                           "--N", str(N), "--box", "2"], stdout=subprocess.DEVNULL)
    box = ",".join(str(v) for v in BOX1 + BOX2)
    cmd = ([launch] + flags + [npflag, str(npr)] if (launch and npr > 1) else []) + [
        drv, prob, out, "--quiet", "--set", "run_steps=6", "--set", "amrex_max_level=2", "--set",
        "amrex_refine_mass=1e9", "--set", "amrex_n_error_buf=0", "--set", "amrex_refine_box=" + box,
        "--set", "amrex_max_grid_size=8", "--set", "fields_levels=1", "--set", f"fields_levels_index={INDEX}"]
    subprocess.check_call(cmd)
    steps = sorted(d for d in os.listdir(out) if d.startswith("step_"))
    d = os.path.join(out, steps[-1])
    t = {}
    for line in open(os.path.join(d, "levels_slice.txt")):
        if " = " in line and not line.startswith("#"):
            k, v = line.split(" = ", 1)
            t[k.strip()] = v.strip()
    bad = []
    n, F = int(t["n"]), int(t["cells_per_base_cell"])
    if (n, F, int(t["finest"]), t["axis"]) != (N * 4, 4, 2, "z"):
        bad.append(f"geometry {t}")
    if int(t["plane"]) != INDEX * 4 + 2:
        bad.append(f"plane {t['plane']} != {INDEX * 4 + 2}")
    a = np.fromfile(os.path.join(d, "levels_slice.bin"), dtype=np.float32).reshape(n, n, int(t["ncomp"]))
    if not np.all(np.isfinite(a)):
        bad.append("non-finite values")
    lev = a[:, :, 5]
    want = np.zeros((n, n))
    c = (np.arange(n) + 0.5) / F  # pixel centres in base cells
    for L, b in ((1, BOX1), (2, BOX2)):
        inx = (c >= b[0]) & (c < b[3] + 1)
        iny = (c >= b[1]) & (c < b[4] + 1)
        want[np.ix_(iny, inx)] = L
    if not np.array_equal(lev, want):
        bad.append(f"level map: {int((lev != want).sum())} pixels differ")
    base = np.fromfile(os.path.join(d, "fHI.bin"), dtype=np.float32).reshape(N, N, N)[INDEX]
    m0 = (lev.reshape(N, F, N, F) == 0).all(axis=(1, 3))
    if not np.array_equal(a[:, :, 1].reshape(N, F, N, F)[:, 0, :, 0][m0], base[m0]):
        bad.append("level-0 pixels differ from fHI.bin")
    fine = lev > 0
    for m, name in ((1, "fHI"), (2, "fHeI"), (3, "fHeII")):
        v = a[:, :, m][fine]
        if v.size == 0 or v.min() < 0 or v.max() > 1:
            bad.append(f"{name} on the refined pixels outside [0, 1]")
    print(f"{steps[-1]}: n {n}, levels {dict(zip(*np.unique(lev, return_counts=True)))}")
    if bad:
        print("FAIL: " + "; ".join(bad))
        sys.exit(1)
    print("PASS")


if __name__ == "__main__":
    main()
