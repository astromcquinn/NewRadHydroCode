#!/usr/bin/env python3
"""Physics sanity checks of hydro runs (C++ hydro_driver or Fortran harness).

usage: hydro_physics_check.py RUN_ROOT [--strict]

RUN_ROOT holds one output directory per problem of tests/problems_hydro
(sod_x, sod_y, sod_z, sound_x, random3d, random3d_freeze); missing problems
are skipped.  Checks (numpy only):

  sod_*    density / velocity / pressure profiles against the exact Riemann
           solution (gamma = real(4) 5/3, both periodic interfaces, advected
           by the bulk velocity) at the final time t = sum(dt_sim); L1 errors
           and shock position.  sod_y and sod_z profiles against sod_x.
  sound_x  phase speed of the k = 1 density mode (expected c0 + vbx) and its
           amplitude decay.
  all      conservation of mass (exact up to float32 round-off), lab-frame
           momentum and energy, and of the passive w_gas totals.

The thresholds are deliberately loose (the scheme is 2nd order, N = 32, with
an 8-cell Gaussian moving frame); exit status 1 if any is violated.
"""
import math
import os
import sys

import numpy as np

GAMMA = float(np.float32(5.0) / np.float32(3.0))
DGAS_AVG = float(np.float32(0.045)) / float(np.float32(0.27))


def read(d, name, nc, N):
    a = np.fromfile(os.path.join(d, name + ".bin"), dtype="<f4")
    return a.reshape((nc, N, N, N), order="F").astype(np.float64)


def stats(d):
    out = {}
    for line in open(os.path.join(d, "stats.txt")):
        if "=" in line:
            k, v = line.split("=", 1)
            vals = [float(x.replace("D", "E")) for x in v.split()]
            out[k.strip()] = vals if len(vals) > 1 else vals[0]
    return out


def steps(run):
    return sorted(int(n[5:]) for n in os.listdir(run) if n.startswith("step_"))


def primitive(u, v):
    rho = u[0]
    dv = u[1:4] / rho
    vel = dv + v
    eth = u[4] - 0.5 * rho * (dv ** 2).sum(axis=0)
    p = (GAMMA - 1) * eth
    return rho, vel, p


# ---------------------------------------------------------------- Riemann
def exact_riemann(rl, ul, pl, rr, ur, pr, xi, g=GAMMA):
    """Exact solution of the Riemann problem at similarity coordinates xi."""
    cl, cr = math.sqrt(g * pl / rl), math.sqrt(g * pr / rr)

    def f(p, rk, pk, ck):
        if p > pk:
            A, B = 2 / ((g + 1) * rk), (g - 1) / (g + 1) * pk
            return (p - pk) * math.sqrt(A / (p + B))
        return 2 * ck / (g - 1) * ((p / pk) ** ((g - 1) / (2 * g)) - 1)

    lo, hi = 1e-12 * min(pl, pr), 10 * max(pl, pr)
    for _ in range(200):  # bisection on f_L + f_R + du = 0
        pm = 0.5 * (lo + hi)
        val = f(pm, rl, pl, cl) + f(pm, rr, pr, cr) + (ur - ul)
        if val > 0:
            hi = pm
        else:
            lo = pm
    ps = 0.5 * (lo + hi)
    us = 0.5 * (ul + ur) + 0.5 * (f(ps, rr, pr, cr) - f(ps, rl, pl, cl))
    out = np.zeros((3, xi.size))
    for n, s in enumerate(xi):
        if s < us:  # left of contact
            if ps > pl:  # shock
                sl = ul - cl * math.sqrt((g + 1) / (2 * g) * ps / pl + (g - 1) / (2 * g))
                if s < sl:
                    out[:, n] = (rl, ul, pl)
                else:
                    r = rl * (ps / pl + (g - 1) / (g + 1)) / ((g - 1) / (g + 1) * ps / pl + 1)
                    out[:, n] = (r, us, ps)
            else:  # rarefaction
                csl = cl * (ps / pl) ** ((g - 1) / (2 * g))
                if s < ul - cl:
                    out[:, n] = (rl, ul, pl)
                elif s > us - csl:
                    out[:, n] = (rl * (ps / pl) ** (1 / g), us, ps)
                else:
                    uu = 2 / (g + 1) * (cl + (g - 1) / 2 * ul + s)
                    c = 2 / (g + 1) * (cl + (g - 1) / 2 * (ul - s))
                    r = rl * (c / cl) ** (2 / (g - 1))
                    out[:, n] = (r, uu, pl * (c / cl) ** (2 * g / (g - 1)))
        else:  # right of contact
            if ps > pr:
                sr = ur + cr * math.sqrt((g + 1) / (2 * g) * ps / pr + (g - 1) / (2 * g))
                if s > sr:
                    out[:, n] = (rr, ur, pr)
                else:
                    r = rr * (ps / pr + (g - 1) / (g + 1)) / ((g - 1) / (g + 1) * ps / pr + 1)
                    out[:, n] = (r, us, ps)
            else:
                csr = cr * (ps / pr) ** ((g - 1) / (2 * g))
                if s > ur + cr:
                    out[:, n] = (rr, ur, pr)
                elif s < us + csr:
                    out[:, n] = (rr * (ps / pr) ** (1 / g), us, ps)
                else:
                    uu = 2 / (g + 1) * (-cr + (g - 1) / 2 * ur + s)
                    c = 2 / (g + 1) * (cr - (g - 1) / 2 * (ur - s))
                    r = rr * (c / cr) ** (2 / (g - 1))
                    out[:, n] = (r, uu, pr * (c / cr) ** (2 * g / (g - 1)))
    shock = ur + cr * math.sqrt((g + 1) / (2 * g) * ps / pr + (g - 1) / (2 * g))
    return out, shock


def sod_profile(run, axis, N):
    """Final profile along the tube axis (averaged over the transverse plane)."""
    last = steps(run)[-1]
    d = os.path.join(run, f"step_{last:04d}")
    u, v = read(d, "u_gas", 5, N), read(d, "v_gas", 3, N)
    rho, vel, p = primitive(u, v)
    mean_axes = tuple(a for a in range(3) if a != axis)
    prof = np.array([rho.mean(axis=mean_axes), vel[axis].mean(axis=mean_axes),
                     p.mean(axis=mean_axes)])
    spread = max(float(np.ptp(rho, axis=mean_axes).max()), 0.0)
    t = sum(stats(os.path.join(run, f"step_{s:04d}"))["dt_sim"] for s in steps(run) if s > 0)
    return prof, t, spread


def check_sod(root, fails):
    N = 32
    c_L = 2.0
    rL, rR = DGAS_AVG, 0.125 * DGAS_AVG
    pL = rL * c_L ** 2 / GAMMA
    pR = 0.1 * pL
    bulk_x = np.array([0.25, -0.15, 0.1])
    profs = {}
    for axis, name in enumerate(("sod_x", "sod_y", "sod_z")):
        run = os.path.join(root, name)
        if not os.path.isdir(run):
            continue
        prof, t, spread = sod_profile(run, axis, N)
        vb = np.roll(bulk_x, axis)[axis]
        x = np.arange(N) + 0.5
        # interface at x = N/2: L | R ; interface at x = 0 (= N): R | L
        ex = np.empty((3, N))
        sol1, shock = exact_riemann(rL, 0, pL, rR, 0, pR, (x - N / 2 - vb * t) / t)
        xm = np.where(x < N / 2, x, x - N)        # coordinate relative to face 0
        sol2, _ = exact_riemann(rR, 0, pR, rL, 0, pL, (xm - vb * t) / t)
        right_half = (x >= N / 4) & (x < 3 * N / 4)
        ex[:, right_half] = sol1[:, right_half]
        ex[:, ~right_half] = sol2[:, ~right_half]
        ex[1] += vb
        l1 = [float(np.abs(prof[q] - ex[q]).mean() / np.abs(ex[q]).mean())
              if q != 1 else float(np.abs(prof[1] - ex[1]).mean() / (c_L))
              for q in range(3)]
        # shock position (right-moving shock of the face N/2): where the
        # velocity crosses the mean of post-shock and pre-shock velocities,
        # searched in cells N/2+1 .. 3N/4, linear interpolation between centres
        v_post = float(ex[1][(x > N / 2) & (x < N / 2 + (shock + vb) * t - 0.5)].max(initial=vb))
        vmid = 0.5 * (v_post + vb)
        xs_num = float("nan")
        for i in range(N // 2 + 1, 3 * N // 4):
            if prof[1][i] < vmid <= prof[1][i - 1]:
                f = (prof[1][i - 1] - vmid) / (prof[1][i - 1] - prof[1][i])
                xs_num = x[i - 1] + f
                break
        xs_ex = N / 2 + (shock + vb) * t
        waves_meet = (shock + vb) * t > N / 4 + 1
        print(f"{name}: t = {t:.4f} code units, shock speed {shock:.4f} (+bulk {vb:+.2f}), "
              f"transverse spread {spread:.2e}")
        print(f"   L1 rel: rho {l1[0]:.4f}  v/c_L {l1[1]:.4f}  p {l1[2]:.4f}; "
              f"shock at {xs_num:.2f} (exact {xs_ex:.2f}){'  [waves interacting]' if waves_meet else ''}")
        if spread > 1e-6 * DGAS_AVG:
            fails.append(f"{name}: not uniform across the tube ({spread:.3e})")
        if not (l1[0] < 0.12 and l1[1] < 0.10 and l1[2] < 0.12):
            fails.append(f"{name}: L1 errors too large {l1}")
        if not abs(xs_num - xs_ex) <= 1.5:
            fails.append(f"{name}: shock position {xs_num} vs exact {xs_ex:.2f}")
        profs[name] = (prof, t)
    # The Courant step depends on the orientation (cmax_gas includes
    # |vel - vg|, and vg is smoothed only along the last sweep direction), so
    # tubes are compared with each other only when they reached the same time.
    if "sod_x" in profs:
        for name in ("sod_y", "sod_z"):
            if name in profs:
                (pa, ta), (pb, tb) = profs["sod_x"], profs[name]
                dmax = float(np.abs(pb[0] - pa[0]).max() / rL)
                if abs(ta - tb) <= 1e-9 * ta:
                    print(f"{name} vs sod_x (same t): max |drho|/rho_L = {dmax:.3e}")
                    if dmax > 1e-3:
                        fails.append(f"{name}: differs from sod_x by {dmax:.3e}")
                else:
                    print(f"{name} vs sod_x: different final times ({tb:.4f} vs {ta:.4f}), "
                          f"max |drho|/rho_L = {dmax:.3e} (not checked)")


def check_sound(root, fails):
    run = os.path.join(root, "sound_x")
    if not os.path.isdir(run):
        return
    N, A, c0, vbx = 32, 1e-3, 2.0, 0.2
    k = 2 * np.pi / N
    x = np.arange(N) + 0.5
    ts, ph, amp = [], [], []
    t = 0.0
    for s in steps(run):
        d = os.path.join(run, f"step_{s:04d}")
        if s > 0:
            t += stats(d)["dt_sim"]
        rho = read(d, "u_gas", 5, N)[0].mean(axis=(1, 2))
        drho = rho / rho.mean() - 1
        F = (drho * np.exp(-1j * k * x)).sum()
        ts.append(t)
        ph.append(np.angle(F))
        amp.append(2 * abs(F) / N)
    ph = np.unwrap(np.array(ph))
    ts = np.array(ts)
    slope = np.polyfit(ts, ph, 1)[0]
    c_meas = -slope / k
    rel = (c_meas - (c0 + vbx)) / (c0 + vbx)
    decay = amp[-1] / amp[0]
    print(f"sound_x: t = {ts[-1]:.3f}, phase speed {c_meas:.5f} (expected c0+vbx = "
          f"{c0 + vbx:.3f}, rel err {rel:+.2e}), amplitude {amp[0]:.3e} -> {amp[-1]:.3e} "
          f"(ratio {decay:.4f})")
    if abs(rel) > 0.02:
        fails.append(f"sound_x: phase speed off by {rel:.3e}")
    if not (0.8 < decay <= 1.0001):
        fails.append(f"sound_x: amplitude ratio {decay}")


def check_conservation(root, fails):
    print("conservation (relative change, last step vs step 0):")
    for name in sorted(os.listdir(root)):
        run = os.path.join(root, name)
        if not os.path.isdir(run) or not steps(run):
            continue
        st = [stats(os.path.join(run, f"step_{s:04d}")) for s in steps(run)]
        m0, m1 = st[0]["mass"], st[-1]["mass"]
        e0, e1 = st[0]["energy"], st[-1]["energy"]
        p0, p1 = np.array(st[0]["momentum"]), np.array(st[-1]["momentum"])
        w0, w1 = np.array(st[0]["w_total"]), np.array(st[-1]["w_total"])
        dm = (m1 - m0) / m0
        de = (e1 - e0) / e0
        # momentum relative to the total |momentum| scale (sum |rho v|)
        dp = np.abs(p1 - p0).max() / max(np.abs(p0).max(), 1e-300)
        dw = np.where(w0 != 0, np.abs(w1 - w0) / np.where(w0 != 0, np.abs(w0), 1), np.abs(w1))
        print(f"  {name:16s} mass {dm:+.2e}  momentum {dp:.2e}  energy {de:+.2e}  "
              f"w {dw.max():.2e}")
        # mass: conservative fluxes, only float32 round-off (and the
        # Dgas_min clamps of smooth_cells/summarize_gas, none here)
        if abs(dm) > 1e-6:
            fails.append(f"{name}: mass not conserved ({dm:.3e})")
        # lab-frame momentum and energy are NOT conserved exactly by the MF
        # scheme (frame-velocity source terms, the T floor of check_temp,
        # smooth_cells); require only that the drift stays small
        tol = 1e-6 if name.startswith("sound") else 1e-2
        if dp > tol or abs(de) > tol:
            fails.append(f"{name}: momentum/energy drift {dp:.3e} / {de:.3e} > {tol}")


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    root = sys.argv[1]
    fails = []
    check_sod(root, fails)
    check_sound(root, fails)
    check_conservation(root, fails)
    if fails:
        print("FAILED:")
        for f in fails:
            print("  " + f)
        return 1
    print("PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
