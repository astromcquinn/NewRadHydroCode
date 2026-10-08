#!/usr/bin/env python3
"""Quantify deterministic_parallel (colour-major 8^3 domains) vs the legacy
sweep order, with the RT-only driver (legacy_rt = legacy_physics = 0, one
thread; the det result is the same for any thread count).

usage: compare_order.py WORKDIR problem [problem ...]   (tests/problems/<problem>)
"""
import os
import subprocess
import sys

import numpy as np

W = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
S = sys.argv[1]


def parse(fn):
    d = {}
    for line in open(fn):
        if "=" in line:
            k, v = line.split("=", 1)
            d[k.strip()] = [float(x) for x in v.split()]
    return d


def run(p, tag, extra):
    o = os.path.join(S, p, tag)
    os.makedirs(o, exist_ok=True)
    env = dict(os.environ, OMP_NUM_THREADS="1")
    subprocess.run([W + "/build/rt_driver", W + "/tests/problems/" + p, o, "--set", "legacy_rt=0",
                    "--set", "legacy_physics=0"] + extra, check=True, env=env,
                   stdout=subprocess.DEVNULL)
    return o


for p in sys.argv[2:]:
    a = run(p, "legacy_order", [])
    b = run(p, "det", ["--set", "deterministic_parallel=1"])
    print(f"== {p}")
    print(" step  xHII_V(legacy)  rel.diff  xHII_M rel.diff  <T>_M rel.diff  nray rel.diff  "
          "ngam_rays rel.diff  max|dxHII|  rms dxHII  max|dT/T|(ionised)")
    for n in range(1, 11):
        sa, sb = os.path.join(a, f"step_{n:04d}"), os.path.join(b, f"step_{n:04d}")
        fa = np.fromfile(os.path.join(sa, "fHI.bin"), np.float32).astype(float)
        fb = np.fromfile(os.path.join(sb, "fHI.bin"), np.float32).astype(float)
        ta = np.fromfile(os.path.join(sa, "temp.bin"), np.float32).astype(float)
        tb = np.fromfile(os.path.join(sb, "temp.bin"), np.float32).astype(float)
        ra = np.fromfile(os.path.join(sa, "rhog.bin"), np.float32).astype(float)
        xa, xb = 1 - fa, 1 - fb
        xva, xvb = xa.mean(), xb.mean()
        xma, xmb = (ra * xa).sum() / ra.sum(), (ra * xb).sum() / ra.sum()
        tma, tmb = (ra * ta).sum() / ra.sum(), (ra * tb).sum() / ra.sum()
        A, B = parse(os.path.join(sa, "stats.txt")), parse(os.path.join(sb, "stats.txt"))
        ion = xa > 0.5
        dT = np.max(np.abs(tb[ion] / ta[ion] - 1)) if ion.any() else 0
        rel = lambda u, v: (v - u) / u if u else 0.0  # noqa: E731
        print(f" {n:4d} {xva:14.6e} {rel(xva, xvb):+9.2e} {rel(xma, xmb):+9.2e} "
              f"{rel(tma, tmb):+9.2e} {rel(A['nray_active'][0], B['nray_active'][0]):+9.2e} "
              f"{rel(A['ngam_rays_total'][0], B['ngam_rays_total'][0]):+9.2e} "
              f"{np.max(np.abs(xb - xa)):9.2e} {np.sqrt(np.mean((xb - xa) ** 2)):9.2e} {dT:9.2e}")
