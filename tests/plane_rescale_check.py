#!/usr/bin/env python3
"""plane_packet_rescale (docs/plane_parallel.md "Packet rescaling"): the
sweep's packets in flight rescaled to the current RT substep.

usage: plane_rescale_check.py RADHYDRO WORKDIR [--amrex EXE --launch LAUNCHER FLAGS NPFLAG NP]

The problem of tests/plane_full_check.py (a uniform 16^3 box of 0.25 Mpc/h,
2^3 domains lit from +x and +y at Gamma_-12 = 3 from z_re = 9, the speed of
light 0.1 c -> 0.01 c at x_HII,V = 0.99) with plane_transport = sweep, the
hydro step capped at 0.3 Myr (so the RT substep grows ~10x at the switch
and every substep is short against the domain crossing) and a movie frame
every dz = 5e-4 (frames.txt: the volume-weighted Gamma_HI of the last RT
substep). Checks:
  * plane_packet_rescale = 0 (the default): the packets launched at 0.1 c
    hold a tenth of the photons of a 0.01 c substep, so Gamma_HI,V in the
    first frame after the switch is below half its value just before it
    (the bug the flag fixes; documents the default);
  * plane_packet_rescale = 1: Gamma_HI,V in every frame from the switch to
    the end is within 1 per cent of its value in the last frame before the
    switch (the box is ionised and optically thin; the source flux is
    fixed);
  * the photon bookkeeping closes to 1e-3 in every RT step with the
    rescaling (its photons count as emitted);
  * the rescaled run does not depend on the thread count (1 and 4 threads:
    identical stats.txt and frames.txt);
  * --amrex: the backend run (rt_backend = hydro_backend = pm_backend =
    amrex) with the rescaling has the same Gamma_HI,V history to 1e-4.
"""
import argparse
import filecmp
import glob
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import plane_full_check as P  # noqa: E402

check = P.check


def frames(out):
    """(z, xHII_V, GammaHI_V) per distinct movie frame after the turn-on"""
    rows, seen = [], set()
    for line in open(os.path.join(out, "movie", "frames.txt")):
        if line.startswith("#"):
            continue
        f = line.split()
        z = float(f[3])
        if z in seen or z > 9 - 1e-9:  # the RT frames (the turn-on at z_re = 9)
            continue
        seen.add(z)
        rows.append((z, float(f[4]), float(f[6])))
    return np.array(rows)


def switch_index(fr, out):
    """the first frame after the c_reduced switch (stats.txt: the first step at
    0.01 c starts at the previous step's redshift)"""
    st = P.steps(out)
    i = next(n for n, s in enumerate(st) if s["c_reduced"][0] < 0.1 - 1e-12)
    z_sw = st[i - 1]["redshift"][0]
    return int(np.argmax(fr[:, 0] < z_sw - 1e-6))  # frames.txt: 8 decimals


def bookkeeping(out):
    worst, prev = 0.0, 0.0
    for s in P.steps(out):
        if s["plane_emitted"][0] > 0:
            res = (s["plane_emitted"][0] - s["plane_absorbed"][0] - s["plane_escaped"][0] -
                   (s["ngam_rays_total"][0] - prev)) / s["plane_emitted"][0]
            worst = max(worst, abs(res))
        prev = s["ngam_rays_total"][0]
    return worst


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("exe")
    ap.add_argument("work")
    ap.add_argument("--amrex", default=None)
    ap.add_argument("--launch", nargs=4, metavar=("LAUNCHER", "FLAGS", "NPFLAG", "NP"), default=None)
    a = ap.parse_args()
    os.makedirs(a.work, exist_ok=True)
    launch = []
    if a.launch:
        launch = [a.launch[0]] + [f for f in a.launch[1].split(";") if f] + [a.launch[2], a.launch[3]]

    prob = os.path.join(a.work, "prob")
    P.make_problem(prob, 16, 0.25, 20, "15 12 10 9 8.95", True,
                   {"plane_ndom": 2, "plane_faces": "+x +y", "plane_gamma12": 3, "plane_z_re": 9,
                    "c_reduced_schedule": "0.1 0.99 0.01", "cool_suppress_overdensity": 0,
                    "plane_transport": "sweep", "dt_max_myr": 0.3,
                    "movie_dz": 0.0005, "movie_fields": "fHI", "movie_slab": 1})

    print("plane_packet_rescale = 0 (default)")
    o0 = P.run(a.exe, prob, os.path.join(a.work, "r0"), [])
    r0 = frames(o0)
    i0 = switch_index(r0, o0)
    g_before, g_after = r0[i0 - 1, 2], r0[i0, 2]
    check(g_after < 0.5 * g_before,
          f"Gamma_HI,V drops {g_before:.3e} -> {g_after:.3e} at the switch (z = {r0[i0, 0]:.4f}) without the rescaling")

    print("plane_packet_rescale = 1")
    os.environ["OMP_NUM_THREADS"] = "4"
    o1 = P.run(a.exe, prob, os.path.join(a.work, "r1"), ["plane_packet_rescale=1"])
    r1 = frames(o1)
    i1 = switch_index(r1, o1)
    check(i1 == i0, f"the switch at the same frame ({i1} vs {i0})")
    g_ref = r1[i1 - 1, 2]
    dev = float(np.max(np.abs(r1[i1:, 2] / g_ref - 1)))
    check(dev < 0.01, f"Gamma_HI,V within {dev:.2e} (< 1e-2) of {g_ref:.3e} from the switch to z = {r1[-1, 0]:.3f}")
    worst = bookkeeping(o1)
    check(worst < 1e-3, f"photon bookkeeping closes to {worst:.1e}")
    os.environ["OMP_NUM_THREADS"] = "1"
    o1t = P.run(a.exe, prob, os.path.join(a.work, "r1_t1"), ["plane_packet_rescale=1"])
    same = filecmp.cmp(os.path.join(o1, "movie", "frames.txt"), os.path.join(o1t, "movie", "frames.txt"), shallow=False)
    for f in sorted(glob.glob(os.path.join(o1, "step_*", "stats.txt"))):
        g = os.path.join(o1t, os.path.relpath(f, o1))
        same = same and os.path.exists(g) and filecmp.cmp(f, g, shallow=False)
    check(same, "1 and 4 threads: identical stats.txt and frames.txt")

    if a.amrex:
        print("the AMReX backend")
        am = frames(P.run(a.amrex, prob, os.path.join(a.work, "amrex"), [
            "plane_packet_rescale=1", "rt_backend=amrex", "hydro_backend=amrex", "pm_backend=amrex",
            "amrex_state=1", "amrex_max_grid_size=8", "amrex_hydro_max_grid_size=0",
            "deterministic_parallel=0"], launch))
        n = min(len(am), len(r1))
        ok = n == len(r1) and np.allclose(am[:n, 0], r1[:n, 0], rtol=0, atol=1e-9)
        d = float(np.max(np.abs(am[:n, 2] / r1[:n, 2] - 1))) if n else 1.0
        check(ok and d < 1e-4, f"Gamma_HI,V history of the backend within {d:.1e} (< 1e-4) of the serial run ({n} frames)")

    if P.FAILS:
        sys.exit(f"FAIL: {len(P.FAILS)} check(s)")
    print("all checks passed")


if __name__ == "__main__":
    main()
