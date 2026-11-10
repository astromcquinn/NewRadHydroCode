#!/usr/bin/env python3
"""The on-the-fly diagnostics (diag_dz, igm_diag.hpp) and the movie's |v|
(movie_fields vmag) in a plane-parallel full-loop run.

usage: plane_diag_check.py RADHYDRO WORKDIR

A 32^3 Zel'dovich box of 0.25 Mpc/h from z = 20, 2^3 domains lit from +x
and +y at Gamma_-12 = 3 from z_re = 9 (tests/plane_full_check.py's set-up
with moving gas), outputs at z = 9, 8.5, 8. Checks:
  * diag_dz = 0.2 writes Diag/ek.<ztag> at the first step at or below
    z = 9, 8.8, ..., 8 and at the zdata1 outputs, nothing above z = 9, and
    Diag/tdelta.<ztag> at each of them but the landing step of z_re (no RT
    temperatures yet);
  * at the zdata1 steps E(k) and <v^2> (volume and mass weighted) agree with
    tools/analysis/velocity_power.py on the step's u_gas.bin to 1e-5;
  * the T-Delta histogram holds N^3 cells and a mass of N^3 (sum of Delta)
    to 1e-9, and its mass-weighted mean log T is that of the step's
    temp.bin and u_gas.bin to 1e-6;
  * movie_fields vmag writes slab frames of |v| in km/s (checked for its
    size and range, 0 < max |v| < 1e3 km/s);
  * the diagnostics do not change the run: every step's stats and the last
    step's fields are byte-identical to the run without them.
"""
import glob
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, "tools", "analysis"))
from plane_full_check import FAILS, check, make_problem, run, parse  # noqa: E402
import velocity_power as vp  # noqa: E402


def read_ek(path):
    hdr = open(path).read().splitlines()
    v2 = float(hdr[2].split("=")[1].split()[0])
    v2m = float(hdr[2].split("=")[2])
    d = np.loadtxt(path)
    return d, v2, v2m


def main():
    exe, work = sys.argv[1], sys.argv[2]
    os.makedirs(work, exist_ok=True)
    N = 32
    prob = os.path.join(work, "prob")
    make_problem(prob, N, 0.25, 20, "9 8.5 8", False,
                 {"plane_ndom": 2, "plane_faces": "+x +y", "plane_gamma12": 3, "plane_z_re": 9,
                  "c_reduced_schedule": "0.1 0.99 0.01", "cool_suppress_overdensity": 0,
                  "fields": "fHI,temp,rhog,u_gas", "omegam": 0.305, "omegab": 0.048, "hubble0": 0.68, "ics_check": 0})
    base = run(exe, prob, os.path.join(work, "base"), [])
    diag = run(exe, prob, os.path.join(work, "diag"),
               ["diag_dz=0.2", "movie_dz=0.25", "movie_fields=fHI,vmag", "movie_slab=1"])

    print("outputs")
    eks = sorted(glob.glob(os.path.join(diag, "Diag", "ek.*")))
    zs = sorted((float(open(p).read().splitlines()[1].split()[3]) for p in eks), reverse=True)
    # one output at the first step at or below each of 9, 8.8, ..., 8 (steps
    # of up to ~0.12 in z here) and at every zdata1 step
    targets = [9, 8.8, 8.6, 8.4, 8.2, 8]
    hit = all(any(t - 0.15 <= z <= t + 1e-6 for z in zs) for t in targets)
    check(hit, f"E(k) at the first step below each of {targets}: {[round(z, 3) for z in zs]}")
    check(max(zs) <= 9 + 1e-6, "no diagnostics above z = 9")
    ntd = len(glob.glob(os.path.join(diag, "Diag", "tdelta.*")))
    check(ntd == len(eks) - 1, f"T-Delta at every E(k) output but the z_re landing step ({ntd} of {len(eks)})")

    print("E(k) against velocity_power.py at the zdata1 steps")
    P = vp.params_of(prob)
    om, h0 = float(np.float32(P["omegam"])), float(np.float32(P["hubble0"]))
    shells = vp.Shells(N, 0.25e3, 32)
    worst = 0.0
    nmatch = 0
    for sd in sorted(glob.glob(os.path.join(diag, "step_*"))):
        ug = os.path.join(sd, "u_gas.bin")
        if not os.path.exists(ug):
            continue
        z = parse(os.path.join(sd, "stats.txt"))["redshift"][0]
        tag = vp.ztag(z)
        ek = os.path.join(diag, "Diag", "ek." + tag)
        if not os.path.exists(ek):
            continue
        d, v2, v2m = read_ek(ek)
        E = shells.energy(vp.velocity(ug, N, z, 0.25, "none", om, h0))
        Em = shells.energy(vp.velocity(ug, N, z, 0.25, "mass", om, h0))
        ok = E > 1e-12 * E.max()
        rel = max(np.max(np.abs(d[ok, 4] / E[ok] - 1)), np.max(np.abs(d[ok, 5] / Em[ok] - 1)),
                  abs(v2 / (2 * np.sum(E * shells.dk)) - 1), abs(v2m / (2 * np.sum(Em * shells.dk)) - 1))
        worst = max(worst, rel)
        nmatch += 1

        # the T-Delta histogram against the step's fields
        if not os.path.exists(os.path.join(diag, "Diag", "tdelta." + tag)):
            continue  # the landing step of z_re: no RT temperatures yet
        h = np.loadtxt(os.path.join(diag, "Diag", "tdelta." + tag))
        check(abs(h[:, 2].sum() - N ** 3) < 0.5 and abs(h[:, 3].sum() / N ** 3 - 1) < 1e-9,
              f"z = {z:.3f}: the histogram holds N^3 cells and mass N^3")
        T = np.fromfile(os.path.join(sd, "temp.bin"), dtype=np.float32).astype(np.float64)
        rho = np.fromfile(ug, dtype=np.float32).reshape(-1, 5)[:, 0].astype(np.float64)
        D = rho / rho.mean()
        lt_bins = (h[:, 1] + 0.5) * 0.05
        mean_h = np.sum(h[:, 3] * lt_bins) / np.sum(h[:, 3])
        mean_f = np.sum(D * np.log10(T)) / np.sum(D)
        check(abs(mean_h - mean_f) < 0.026, f"z = {z:.3f}: mass-weighted <log T> {mean_h:.4f} (bins) vs {mean_f:.4f} (fields), within half a bin")
    check(nmatch == 3, f"three zdata1 steps compared ({nmatch})")
    check(worst < 1e-5, f"E(k), E_mass(k), <v^2> agree with velocity_power.py to {worst:.2e} < 1e-5")

    print("movie |v|")
    fr = sorted(glob.glob(os.path.join(diag, "movie", "frame_*_vmag.slab")))
    check(len(fr) > 0, f"{len(fr)} vmag frames")
    if fr:
        v = np.fromfile(fr[-1], dtype=np.float32)
        check(v.size == 3 * N * N and 0 < v.max() < 1e3 and v.min() >= 0, f"vmag slab: {v.size} cells, max {v.max():.2f} km/s")

    print("no effect on the run")
    sb = sorted(glob.glob(os.path.join(base, "step_*", "stats.txt")))
    sdg = sorted(glob.glob(os.path.join(diag, "step_*", "stats.txt")))
    same = len(sb) == len(sdg) and all(open(a).read() == open(b).read() for a, b in zip(sb, sdg))
    check(same, f"every step's stats.txt is identical ({len(sb)} steps)")
    last_b, last_d = os.path.dirname(sb[-1]), os.path.dirname(sdg[-1])
    for f in ("fHI.bin", "temp.bin", "rhog.bin", "u_gas.bin"):
        a, b = os.path.join(last_b, f), os.path.join(last_d, f)
        check(open(a, "rb").read() == open(b, "rb").read(), f"last step {f} byte-identical")

    if FAILS:
        print(f"{len(FAILS)} FAILED")
        sys.exit(1)
    print("all passed")


if __name__ == "__main__":
    main()
