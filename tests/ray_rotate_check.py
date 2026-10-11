#!/usr/bin/env python3
"""Grid-aligned photon deficit B16 and its fix ray_rotate (docs/port_notes.md).

Runs rt_driver on tests/problems/uniform_1src for 15 steps (the deficit of
the original peaks near step 15) with ray_rotate = 0 and 2, and measures the
most negative x_HII of the +-x, +-y axis cells minus the mean of the shell at
the same radius, over steps 5-15 (tools/analysis/axis_deficit.py). Checks:
  * ray_rotate = 0 shows the deficit (< -0.3; it is -0.46 at step 15), so
    the test can see it;
  * ray_rotate = 2 removes it (> -0.05 on every step);
  * the I-front radius (volume of x_HII > 0.5) agrees within 3%, and the
    photons emitted (rt_stats ngam) are identical.

    ray_rotate_check.py RT_DRIVER PROBLEM_DIR OUTDIR
"""
import os, shutil, subprocess, sys
import numpy as np

NSTEPS = 15
N = 32; s0 = 15
i = np.arange(N)
X, Y, Z = np.meshgrid(i, i, i, indexing="ij")
R = np.sqrt((X - s0) ** 2 + (Y - s0) ** 2 + (Z - s0) ** 2)
AX = ((Y == s0) & (Z == s0), (X == s0) & (Z == s0))


def parse(fn):
    d = {}
    for line in open(fn):
        if "=" in line:
            k, v = line.split("=", 1)
            d[k.strip()] = [float(x) for x in v.split()]
    return d


def deficit(x):
    d = []
    for m in AX:
        for rr in range(2, 14):
            on = m & (np.round(R) == rr)
            sh = (np.abs(R - rr) < 0.5) & ~m
            if on.any():
                d.append(x[on].mean() - x[sh].mean())
    return min(d)


def start(driver, prob, out, mode):
    if os.path.isdir(out):
        shutil.rmtree(out)
    os.makedirs(out)
    log = open(out + ".stdout", "w")
    return subprocess.Popen([driver, prob, out, "--nsteps", str(NSTEPS), "--set", f"ray_rotate={mode}"],
                            stdout=log, stderr=subprocess.STDOUT)


def measure(out):
    defs, ngam = [], 0.0
    for n in range(1, NSTEPS + 1):
        sd = os.path.join(out, f"step_{n:04d}")
        x = 1 - np.fromfile(os.path.join(sd, "fHI.bin"), "<f4").reshape(N, N, N, order="F")
        if n >= 5:
            defs.append(deficit(x))
        ngam += parse(os.path.join(sd, "rt_stats.txt"))["ngam"][0]
    Rf = ((x > 0.5).sum() * 3 / (4 * np.pi)) ** (1 / 3)
    return np.array(defs), Rf, ngam


def main():
    driver, prob, odir = sys.argv[1:4]
    procs = {m: start(driver, prob, os.path.join(odir, f"rot{m}"), m) for m in (0, 2)}
    for m, p in procs.items():
        if p.wait() != 0:
            print(f"rt_driver failed for ray_rotate={m}")
            return 1
    d0, R0, n0 = measure(os.path.join(odir, "rot0"))
    d2, R2, n2 = measure(os.path.join(odir, "rot2"))
    print("axis deficit, steps 5-%d:" % NSTEPS)
    print("  ray_rotate=0:", " ".join("%+.3f" % v for v in d0))
    print("  ray_rotate=2:", " ".join("%+.3f" % v for v in d2))
    print("front radius %.3f / %.3f, photons %.6e / %.6e" % (R0, R2, n0, n2))
    ok = True
    if not d0.min() < -0.3:
        print("FAIL: the original deficit is not seen"); ok = False
    if not d2.min() > -0.05:
        print("FAIL: ray_rotate=2 leaves an axis deficit"); ok = False
    if not abs(R2 / R0 - 1) < 0.03:
        print("FAIL: front radius changed"); ok = False
    if n0 != n2:
        print("FAIL: emitted photons differ"); ok = False
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
