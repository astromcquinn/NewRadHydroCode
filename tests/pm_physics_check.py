#!/usr/bin/env python3
"""Physics checks of PM runs (apps/pm_driver or reference/harness_pm).

usage:
  pm_physics_check.py pancake  RUN PROBLEM [--tol-x T] [--tol-v T]
  pm_physics_check.py growth   RUN PROBLEM [--tol T] [--kbins 2] [--a-max A]
  pm_physics_check.py momentum RUN PROBLEM [--tol T]
  pm_physics_check.py same-dm  RUN RUN2

RUN is an output directory (step_NNNN/, Pk/), PROBLEM the problem directory
(params.txt, ics_info.txt). Exit status 0 = pass, 1 = fail.

pancake   The 1D Zel'dovich solution x = q - (D(a)/D_c) sin(kq)/k is exact
          before shell crossing (for a plane wave the force is linear in the
          displacement). D is the growing mode of the gravity-only run: the
          gas density is never advected without hydro, so only the DM
          clusters, D'' = (2/3)(1-f_b) a D (tools/ics/pmcosmo.py). Compares
          the particles' positions and velocities along the wave with the
          analytic ones at the end of every step; errors are relative to the
          analytic amplitude. Also checks that the transverse coordinates do
          not move.
growth    Ratio of the DM power spectrum (calc_power_spectrum, Pk/ files) in
          the lowest k bins to the first output, against (D(a)/D(a_1))^2.
momentum  Total momentum of DM + gas is conserved by the PM scheme (same
          CIC kernel for assignment and interpolation, symmetric Green's
          function, antisymmetric gradient). Checks |P(t)-P(0)| against the
          sum of |p| over all particles and cells.
"""

import argparse
import glob
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "tools", "ics"))
import pmcosmo as pc  # noqa: E402


def kv(path):
    d = {}
    for line in open(path):
        line = line.split("#", 1)[0]
        if "=" in line:
            k, v = line.split("=", 1)
            d[k.strip()] = v.strip()
    return d


def steps(run):
    out = sorted(glob.glob(os.path.join(run, "step_[0-9][0-9][0-9][0-9]")))
    if not out:
        sys.exit(f"no steps in {run}")
    return out


def load_step(d, ngas):
    st = {k: float(v) for k, v in kv(os.path.join(d, "stats.txt")).items()}
    xv = np.fromfile(os.path.join(d, "xv_dm.bin"), "<f4").reshape(6, -1, order="F").astype(float)
    u = np.fromfile(os.path.join(d, "u_gas.bin"), "<f4").reshape(5, -1, order="F").astype(float)
    v = np.fromfile(os.path.join(d, "v_gas.bin"), "<f4").reshape(3, -1, order="F").astype(float)
    return st, xv, u, v


def problem(pdir):
    p = kv(os.path.join(pdir, "params.txt"))
    info = kv(os.path.join(pdir, "ics_info.txt"))
    ngas = int(p.get("Ngas", 32))
    ndm = int(p.get("Ndm", ngas // 2))
    return p, info, ngas, ndm


def check_pancake(a):
    p, info, ngas, ndm = problem(a.problem)
    ax = "xyz".index(info["axis"])
    acr = float(info["a_cross"])
    cluster = 1.0 if info.get("growth", "frozen-gas") == "full" else 1.0 - pc.FB
    g = pc.Growth(cluster)
    k = 2 * np.pi / ndm
    q = np.arange(ndm) + 0.5
    Q = [c.ravel(order="F") for c in np.meshgrid(q, q, q, indexing="ij")]
    qa = Q[ax]
    sq = np.sin(k * qa) / k
    other = [c for c in range(3) if c != ax]
    w = dict(px=0.0, pv=0.0, mx=0.0, mv=0.0, perp=0.0)
    A = V = None
    print(f"pancake along {'xyz'[ax]}, a_cross = {acr}")
    print("  errors relative to the analytic amplitude: per particle (max) and of the\n"
          "  fitted sin(kq) mode; 'lf' = mode error of the continuum force law\n"
          "  integrated with the code's own DKD steps (time-integration error only)")
    print(" step  scalefactor   D/Dc  part_x    part_v    mode_x    mode_v    lf_x      lf_v    transverse")
    for d in steps(a.run):
        st, xv, _, _ = load_step(d, ngas)
        sa = st["scalefactor"]
        amp = g(sa) / g(acr)
        vamp = g.dtau(sa) / g(acr)
        dx = (xv[ax] - (qa - amp * sq) + ndm / 2) % ndm - ndm / 2
        dv = xv[3 + ax] + vamp * sq
        px = np.abs(dx).max() / (amp / k)
        pv = np.abs(dv).max() / (vamp / k)
        # mode amplitudes: x - q = -A sin(kq)/k, v = -V sin(kq)/k
        disp = (xv[ax] - qa + ndm / 2) % ndm - ndm / 2
        Am = -(disp * sq).sum() / (sq * sq).sum()
        Vm = -(xv[3 + ax] * sq).sum() / (sq * sq).sum()
        mx, mv = abs(Am / amp - 1), abs(Vm / vamp - 1)
        # continuum force law D'' = (2/3) c a D with the code's DKD steps
        # (dt_sim, scalefactor at mid-step)
        if A is None:
            A, V = amp, vamp
        else:
            dt, am = st["dt_sim"], sa - st["da_sim"] / 2
            A = A + V * dt / 2
            V = V + (2.0 / 3.0) * cluster * am * A * dt
            A = A + V * dt / 2
        lx, lv = abs(A / amp - 1), abs(V / vamp - 1)
        perp = max(max(np.abs((xv[c] - Q[c] + ndm / 2) % ndm - ndm / 2).max(),
                       np.abs(xv[3 + c]).max()) for c in other)
        print(f" {int(st['istep']):4d}  {sa:.6f}  {amp:6.4f}  {px:8.2e}  {pv:8.2e}  {mx:8.2e}  "
              f"{mv:8.2e}  {lx:8.2e}  {lv:8.2e}  {perp:7.1e}")
        for key, val in (("px", px), ("pv", pv), ("mx", mx), ("mv", mv), ("perp", perp)):
            w[key] = max(w[key], val)
    ok = (w["px"] <= a.tol_x and w["pv"] <= a.tol_v and w["mx"] <= a.tol_mode_x
          and w["mv"] <= a.tol_mode_v and w["perp"] == 0.0)
    print(f"max relative error per particle: x {w['px']:.3e} (tol {a.tol_x}), "
          f"v {w['pv']:.3e} (tol {a.tol_v})")
    print(f"max relative error of the mode : x {w['mx']:.3e} (tol {a.tol_mode_x}), "
          f"v {w['mv']:.3e} (tol {a.tol_mode_v})")
    print(f"transverse motion {w['perp']:.1e} (must be 0)")
    print("RESULT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def check_growth(a):
    p, info, ngas, ndm = problem(a.problem)
    cluster = 1.0 if info.get("growth", "full") == "full" else 1.0 - pc.FB
    g = pc.Growth(cluster)
    rows = []
    for d in steps(a.run):
        st = {k: float(v) for k, v in kv(os.path.join(d, "stats.txt")).items()}
        if st["istep"] == 0:
            continue
        z = st["redshift"]
        tag = ("z=0%6.4f" % z) if z < 10 else ("z=%7.4f" % z) if z < 100 else ("z=%7.3f" % z)
        f = os.path.join(a.run, "Pk", "pk." + tag)
        if os.path.isfile(f):
            rows.append((st["scalefactor"], np.loadtxt(f)))
    if len(rows) < 2:
        print("not enough Pk files")
        return 1
    a1, pk1 = rows[0]
    print(f"DM power growth relative to a = {a1:.5f} (bins k = "
          + ", ".join(f"{x:.3g}" for x in pk1[:a.kbins, 0]) + " h/Mpc)")
    print(" scalefactor  (D/D1)^2   " + "  ".join(f"ratio_k{i + 1}" for i in range(a.kbins))
          + "   max|err|")
    worst = 0.0
    for sa, pk in rows:
        expect = (g(sa) / g(a1)) ** 2
        r = pk[:a.kbins, 1] / pk1[:a.kbins, 1]
        err = np.abs(r / expect - 1).max()
        tested = sa <= a.a_max
        if tested:
            worst = max(worst, err)
        print(f" {sa:.6f}  {expect:9.4f}  " + "  ".join(f"{x:9.4f}" for x in r)
              + f"  {err:9.3e}" + ("" if tested else "  (a > a_max: not tested)"))
    ok = worst <= a.tol
    print(f"max relative deviation from linear growth for a <= {a.a_max}: {worst:.3e} "
          f"(tol {a.tol})")
    print("RESULT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def check_momentum(a):
    p, info, ngas, ndm = problem(a.problem)
    ndg = ngas / ndm                     # < 1 with pm_refine and Ndm > Ngas
    mp = (1 - pc.FB) * ndg ** 3          # DM particle mass (gas cell units)
    P0 = None
    worst = 0.0
    print(" step  P_dm(x,y,z)                          P_gas(x,y,z)                         |dP|/sum|p|")
    for d in steps(a.run):
        st, xv, u, v = load_step(d, ngas)
        pdm = mp * ndg * xv[3:6]                 # gas grid velocity units
        pgas = u[1:4] + u[0] * v                 # MF: rho*(dv + v_grid); non-MF: v = 0
        P = pdm.sum(axis=1) + pgas.sum(axis=1)
        scale = np.abs(pdm).sum() + np.abs(pgas).sum()
        if P0 is None:
            P0 = P
        err = np.abs(P - P0).max() / scale
        worst = max(worst, err)
        print(f" {int(st['istep']):4d}  " + " ".join(f"{x:11.3e}" for x in pdm.sum(axis=1))
              + "  " + " ".join(f"{x:11.3e}" for x in pgas.sum(axis=1)) + f"  {err:9.2e}")
    ok = worst <= a.tol
    print(f"max |P - P0| / sum|p| = {worst:.3e} (tol {a.tol})")
    print("RESULT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def check_same_dm(a):
    """The DM of two runs is bit-identical at every step (a pm_refine > 1
    run and its pm_refine = 1 twin with the same mesh and particles)."""
    s1, s2 = steps(a.run), steps(a.run2)
    ok = len(s1) == len(s2)
    for d1, d2 in zip(s1, s2):
        b1 = open(os.path.join(d1, "xv_dm.bin"), "rb").read()
        b2 = open(os.path.join(d2, "xv_dm.bin"), "rb").read()
        same = b1 == b2
        ok = ok and same
        print(f" {os.path.basename(d1)}: xv_dm {'identical' if same else 'DIFFERENT'}")
    print("RESULT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("pancake")
    s.add_argument("run")
    s.add_argument("problem")
    s.add_argument("--tol-x", type=float, default=0.04)
    s.add_argument("--tol-v", type=float, default=0.07)
    s.add_argument("--tol-mode-x", type=float, default=0.015)
    s.add_argument("--tol-mode-v", type=float, default=0.02)
    s = sub.add_parser("growth")
    s.add_argument("run")
    s.add_argument("problem")
    s.add_argument("--tol", type=float, default=0.03)
    s.add_argument("--kbins", type=int, default=2)
    s.add_argument("--a-max", type=float, default=0.02,
                   help="test only outputs with a <= a_max (linear regime)")
    s = sub.add_parser("momentum")
    s.add_argument("run")
    s.add_argument("problem")
    s.add_argument("--tol", type=float, default=1e-5)
    s = sub.add_parser("same-dm")
    s.add_argument("run")
    s.add_argument("run2")
    a = ap.parse_args()
    return {"pancake": check_pancake, "growth": check_growth,
            "momentum": check_momentum, "same-dm": check_same_dm}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
