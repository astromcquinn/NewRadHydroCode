#!/usr/bin/env python3
"""The plane-parallel RT-domain mode in the full loop (apps/radhydro,
docs/plane_parallel.md): the turn-on at z_re, the adaptive speed of light,
the two step caps, and (optionally) the AMReX backend against the serial RT.

usage: plane_full_check.py RADHYDRO WORKDIR [--amrex EXE --launch LAUNCHER FLAGS NPFLAG NP]
                           [--dm20] [--set key=value ...]

Default: a uniform 16^3 box of 0.25 Mpc/h (uniform gas, the dark matter at
rest on its lattice, as tests/problems_full/make_ifront.py), hydro and
gravity from z = 20, plane_parallel = 1 with 2^3 domains lit from +x and +y
at Gamma_-12 = 3 from z_re = 9, and the DM20 speed-of-light schedule
c_reduced_schedule = "0.1 0.99 0.01" (0.1 c until x_HI,V < 1 per cent,
then 0.01 c). Checks:
  * the step before z_re ends exactly at z_re (to 1e-9), no photons are
    emitted before it, photons are emitted in every step after it;
  * the schedule: c_reduced = 0.1 while the previous step's x_HII,V < 0.99
    and 0.01 after, exactly one switch, at which the RT substeps per step
    drop by more than 5x;
  * the box is ionised (x_HII,V > 0.99) at the end and heated (T_M in
    1e4 .. 4e4 K) after the first RT step;
  * dt_max_myr = 2: every step from z_re on is at most 2 Myr (the steps
    before are the normal ones); rt_substeps_max = 50: every RT step holds
    at most 50 substeps; both capped runs end within 2 per cent of the
    uncapped x_HII,V and T_M (the same photons, redistributed in time);
  * plane_turnon_cfl (the default 1): the step that starts at z_re is at
    most plane_turnon_courant times the CFL step with the sound speed of
    fully ionised gas at plane_turnon_temp added, dt_sim (cmax_gas +
    c_hot) = plane_turnon_courant in code units,
    c_hot = sqrt(gamma k T (2X + 3Y/4) / m_H) converted here from first
    principles (vel_unit = x_unit / time_unit), i.e. dt_yr <=
    dx_phys / (v_max + c_hot) at z_re up to the growth of time_unit ~ a^2
    over the step; the steps before z_re are those of the uncapped run,
    and with plane_turnon_cfl = 0 the turn-on step is the old one. The
    check runs with plane_turnon_temp = 4e4 and plane_turnon_courant = 0.5
    (36 Myr) to z = 8.5, against the same run with plane_turnon_cfl = 0
    (the ra_max step, 97 Myr, in this uniform, cold 16^3 box at z = 9), and
    the default run's turn-on step is plane_turnon_courant = 0.05 (the
    default, 2026-10-28) times the CFL step with c_hot(2e4 K);
  * --amrex: the backend run (rt_backend = hydro_backend = pm_backend =
    amrex, amrex_state = 1) agrees with the serial run to round-off:
    x_HII,V and T_M to 1e-4 and fHI cell by cell to 1e-3 at the last step.

--dm20 (label nightly): the scaled-down physics sanity check of
docs/plane_parallel.md: 0.256 Mpc/h, 32^3 (DM20: 1024^3), Zel'dovich
initial conditions from z = 100, z_re = 8, Gamma_-12 = 0.3, 2^3 domains,
cooling on everywhere (cool_suppress_overdensity = 0), run to z = 7.
Checks: x_HII,V > 0.99 at z = 7, the mass-weighted temperature of the
ionised gas in 1e4 .. 3e4 K (DM20: about 2e4 K just behind the fronts),
the ionised-gas clumping factor <n^2>/<n>^2 in 1 .. 10 (DM20 Fig. 4:
~ 2-3 at 1 kpc/h; coarser here), and the photon bookkeeping closes to
1e-3 in every RT step.
"""
import argparse
import glob
import os
import shutil
import struct
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
FAILS = []


def check(cond, msg):
    print(("  ok   " if cond else "  FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


def parse(fn):
    d = {}
    for line in open(fn):
        if "=" in line and not line.startswith("#"):
            k, v = line.split("=", 1)
            try:
                d[k.strip()] = [float(x) for x in v.split()]
            except ValueError:
                pass
    return d


def make_problem(out, N, box, z_init, zdata1, uniform, extra):
    os.makedirs(out, exist_ok=True)
    ics = os.path.join(out, "ics")
    cmd = [sys.executable, os.path.join(ROOT, "tools/ics/make_ics.py"), ics, "--kind", "zeldovich",
           "--Ngas", str(N), "--Ndm", str(N if not uniform else N // 2), "--z-init", str(z_init),
           "--box", str(box), "--seed", "7"]
    if uniform:
        cmd += ["--gas", "uniform", "--growth", "frozen-gas"]
    subprocess.check_call(cmd, stdout=subprocess.DEVNULL)
    ndm = N // 2 if uniform else N
    if uniform:
        # the dark matter at rest on its lattice (make_ifront.py): nothing moves
        q = np.arange(ndm, dtype=np.float32) + np.float32(0.5)
        xv = np.zeros((6, ndm, ndm, ndm), dtype="<f4")
        for c, Q in enumerate(np.meshgrid(q, q, q, indexing="ij")):
            xv[c] = Q
        xv.reshape(6, -1, order="F").ravel(order="F").tofile(os.path.join(ics, "xv_dm.init"))
    with open(os.path.join(out, "params.txt"), "w") as f:
        f.write(f"""# plane-parallel RT-domain mode in the full loop (tests/plane_full_check.py)
N = {N}
Ndm = {ndm}
box = {box}
nray_max_factor = 16
z_init = {z_init}
zdata1 = {zdata1}
ics = ics
legacy_rt = 0
legacy_physics = 0
hydro_scheme = unsplit
deterministic_parallel = 1
rt_substep_carry = 1
redshift_per_distance = 1
plane_parallel = 1
photon_consume_dense = 1
fields = fHI,temp,rhog
fields_at_data1 = 1
restart_every = 0
data_files = gas
max_steps = 5000
""")
        for k, v in extra.items():
            f.write(f"{k} = {v}\n")


def run(exe, prob, out, sets, launch=()):
    if os.path.isdir(out):
        shutil.rmtree(out)
    os.makedirs(out)
    cmd = list(launch) + [exe, prob, out, "--quiet"]
    for kv in sets:
        cmd += ["--set", kv]
    with open(out + ".log", "w") as log:
        r = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
    if r.returncode != 0:
        print(open(out + ".log").read()[-4000:])
        sys.exit(f"FAIL: {' '.join(cmd)}")
    return out


def steps(out):
    return [parse(p) for p in sorted(glob.glob(os.path.join(out, "step_*", "stats.txt")))]


def matched(ref, other, z_re):
    """The largest relative difference of T_M and x_HII,V between two runs at
    the reference run's step redshifts below z_re (the other interpolated in
    z; its last step may end at another redshift)."""
    zr = np.array([s["redshift"][0] for s in ref])
    zo = np.array([s["redshift"][0] for s in other])
    worst = 0.0
    for key in ("T_M", "xHII_V"):
        vr = np.array([s[key][0] for s in ref])
        vo = np.array([s[key][0] for s in other])
        sel = (zr < z_re - 1e-6) & (zr >= zo.min()) & np.isfinite(vr)
        vi = np.interp(zr[sel], zo[::-1], vo[::-1])
        worst = max(worst, float(np.max(np.abs(vi / vr[sel] - 1))))
    return worst


# the code's constants (include/radray/constants.hpp; real(4) literals)
F32 = lambda x: float(np.float32(x))
OM, H0, GAMMA, XHY, YHE = F32(0.27), F32(0.70), F32(5.0 / 3.0), F32(0.76), F32(0.24)
H0_CGS, MPC2CM, K_CGS, MH_CGS, YR2SEC = 3.24086e-18, 3.08560e24, 1.38066e-16, 1.67223e-24, 365 * 86400.0


def c_hot_kms(T):
    """The sound speed of fully ionised H + He (mu = 1/(2X + 3Y/4)) at T [km/s]."""
    return np.sqrt(GAMMA * K_CGS * T * (2 * XHY + 3 * YHE / 4) / MH_CGS) / 1e5


def vel_unit_kms(a, box, N):
    """calc_sim_var's vel_unit = x_unit / time_unit at a [km/s per code unit]."""
    x_unit = a * box / H0 / N * MPC2CM
    time_unit = 2 / (3 * H0_CGS * H0 * np.sqrt(OM)) * a * a
    return x_unit / time_unit / 1e5


T_HOT, COURANT = 4e4, 0.5
TURNON_SETS = [f"plane_turnon_temp={T_HOT:g}", f"plane_turnon_courant={COURANT:g}"]


def turnon_checks(a, prob, base, st, N, box, z_re):
    """plane_turnon_cfl: the turn-on step at plane_turnon_temp = 4e4 K and
    plane_turnon_courant = 0.5 (binds here) against the same run with
    plane_turnon_cfl = 0, and the default run's (Courant 0.05 at 2e4 K)."""
    stop = ["zdata1=15 12 10 9 8.5"]
    hot = steps(run(a.exe, prob, os.path.join(a.work, "turnon_hot"), a.set + stop + TURNON_SETS))
    off = steps(run(a.exe, prob, os.path.join(a.work, "turnon_off"), a.set + stop + TURNON_SETS + ["plane_turnon_cfl=0"]))
    land = lambda ss: [i for i, s in enumerate(ss) if abs(s["redshift"][0] - z_re) < 1e-9]
    lh, lo, lb = land(hot), land(off), land(st)
    check(len(lh) == 1 and len(lo) == 1 and len(lb) == 1 and lh == lo == lb, f"the same landing step on z_re ({lh}, {lo}, {lb})")
    if not (lh and lo and lb) or lh[0] + 1 >= len(hot) or lo[0] + 1 >= len(off):
        return
    il = lh[0]
    same = all(hot[i]["dt_yr"][0] == st[i]["dt_yr"][0] and hot[i]["redshift"][0] == st[i]["redshift"][0] for i in range(il + 1))
    check(same, "the steps up to z_re are those of the default run, bit for bit")
    # the rule in code units: dt_sim (cmax_gas + c_hot) = 1 (cmax_gas of stats.txt is the code
    # value at the end of the landing step, the one the turn-on step uses)
    a_re = 1 / (1 + z_re)
    vu = vel_unit_kms(a_re, box, N)
    ch = c_hot_kms(T_HOT) / vu
    cmax = hot[il]["cmax_gas"][0]
    s1, s0 = hot[il + 1], off[il + 1]
    r = s1["dt_sim"][0] * (cmax + ch)
    check(abs(r / COURANT - 1) < 1e-4, f"turn-on step: dt_sim (cmax_gas + c_hot) = {r:.7f} = plane_turnon_courant {COURANT} to 1e-4 (c_hot = {c_hot_kms(T_HOT):.2f} km/s = {ch:.4f} code units)")
    # the same in physical units: dt <= dx_phys / (v_max + c_hot) at z_re, up to the growth of
    # time_unit ~ a^2 over the step (the code velocities are held fixed through the step)
    dx_km = a_re * box / H0 / N * MPC2CM / 1e5
    vmax = hot[il]["vmax_gas"][0] * vu
    bound = COURANT * dx_km / (vmax + c_hot_kms(T_HOT)) / YR2SEC
    grow = ((a_re + s1["da_sim"][0]) / a_re) ** 2
    check(s1["dt_yr"][0] <= bound * grow, f"turn-on step dt = {s1['dt_yr'][0] / 1e6:.3f} Myr <= {COURANT} dx/(v_max + c_hot) = {bound / 1e6:.3f} Myr (x {grow:.4f}, the expansion over the step)")
    check(s0["dt_yr"][0] > 1.5 * s1["dt_yr"][0],
          f"plane_turnon_cfl = 0: the old turn-on step ({s0['dt_yr'][0] / 1e6:.2f} Myr) against {s1['dt_yr'][0] / 1e6:.2f} Myr")
    # the default run (plane_turnon_temp = 2e4 K, plane_turnon_courant = 0.05): its turn-on step
    # is 0.05 times the CFL step with c_hot(2e4 K)
    ch2 = c_hot_kms(2e4) / vu
    r2 = st[il + 1]["dt_sim"][0] * (st[il]["cmax_gas"][0] + ch2)
    check(abs(r2 / 0.05 - 1) < 1e-4,
          f"the default run's turn-on step: dt_sim (cmax_gas + c_hot(2e4 K)) = {r2:.7f} = 0.05 to 1e-4 "
          f"({st[il + 1]['dt_yr'][0] / 1e6:.2f} Myr)")
    return hot


def last_fields(out, N):
    last = sorted(glob.glob(os.path.join(out, "step_*", "fHI.bin")))[-1]
    d = os.path.dirname(last)
    f = lambda n: np.fromfile(os.path.join(d, n + ".bin"), dtype=np.float32).astype(np.float64).reshape(N, N, N)
    return f("fHI"), f("temp"), f("rhog")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("exe")
    ap.add_argument("work")
    ap.add_argument("--amrex", default=None, help="the radhydro of an AMReX build")
    ap.add_argument("--launch", nargs=4, metavar=("LAUNCHER", "FLAGS", "NPFLAG", "NP"), default=None)
    ap.add_argument("--dm20", action="store_true")
    ap.add_argument("--set", action="append", default=[])
    a = ap.parse_args()
    os.makedirs(a.work, exist_ok=True)
    launch = []
    if a.launch:
        launch = [a.launch[0]] + [f for f in a.launch[1].split(";") if f] + [a.launch[2], a.launch[3]]

    if a.dm20:
        return dm20(a)

    N, z_re = 16, 9.0
    prob = os.path.join(a.work, "prob")
    make_problem(prob, N, 0.25, 20, "15 12 10 9 8.5 8 7.5 7", True,
                 {"plane_ndom": 2, "plane_faces": "+x +y", "plane_gamma12": 3, "plane_z_re": z_re,
                  "c_reduced_schedule": "0.1 0.99 0.01", "cool_suppress_overdensity": 0})
    base = run(a.exe, prob, os.path.join(a.work, "base"), a.set)
    st = steps(base)
    print("turn-on at z_re and the schedule")
    zs = np.array([s["redshift"][0] for s in st])
    em = np.array([s["plane_emitted"][0] for s in st])
    land = np.where(np.abs(zs - z_re) < 1e-9)[0]
    check(len(land) == 1, f"one step ends exactly at z_re = {z_re} (|z - z_re| < 1e-9): steps {land}")
    il = land[0] if len(land) else len(zs)
    check(bool(np.all(em[:il + 1] == 0)), "no photons before the turn-on (the landing step included)")
    check(bool(np.all(em[il + 1:] > 0)), "photons emitted in every step after the turn-on")
    # the schedule (c_schedule_check.py's rule): x of the previous step, 0 before the RT wrote its grid
    cs = [s["c_reduced"][0] for s in st]
    prev_x, want = None, []
    for i, s in enumerate(st):
        if i == 0:
            prev_x = s["xHII_V"][0]
            want.append(0.1)
            continue
        x = 0.0 if prev_x == 1.0 else prev_x
        want.append(0.1 if x < 0.99 else 0.01)
        prev_x = s["xHII_V"][0]
    check(cs[1:] == want[1:], "c_reduced follows the schedule 0.1 -> 0.01 at x_HII,V = 0.99 step by step")
    sw = [i for i in range(1, len(cs)) if cs[i] != cs[i - 1]]
    check(len(sw) == 1, f"exactly one switch (steps {sw})")
    if sw:
        # per Myr: the step before the switch is the short turn-on step (plane_turnon_courant 0.05)
        s_before, s_after = st[sw[0] - 1], st[sw[0]]
        r0 = s_before["nstep_rt"][0] / (s_before["dt_yr"][0] / 1e6)
        r1 = s_after["nstep_rt"][0] / (s_after["dt_yr"][0] / 1e6)
        check(r0 / max(r1, 1e-30) > 5, f"the RT substeps per Myr drop {r0:.2f} -> {r1:.2f} at the switch (> 5x)")
    x_end, T_end = st[-1]["xHII_V"][0], st[-1]["T_M"][0]
    check(x_end > 0.99, f"x_HII,V = {x_end:.4f} > 0.99 at the end")
    T_first = st[il + 1]["T_M"][0] if il + 1 < len(st) else 0
    check(1e4 < T_first < 4e4, f"T_M = {T_first:.0f} K in 1e4 .. 4e4 after the first RT step")
    check(5e3 < T_end < 4e4, f"T_M = {T_end:.0f} K in 5e3 .. 4e4 at the end (z = {zs[-1]:.2f}: expansion and recombination cooling)")
    worst = 0.0
    prev = 0.0
    for s in st:
        if s["plane_emitted"][0] > 0:
            res = (s["plane_emitted"][0] - s["plane_absorbed"][0] - s["plane_escaped"][0] -
                   (s["ngam_rays_total"][0] - prev)) / s["plane_emitted"][0]
            worst = max(worst, abs(res))
        prev = s["ngam_rays_total"][0]
    check(worst < 1e-3, f"photon bookkeeping closes to {worst:.1e}")

    print("the turn-on step (plane_turnon_cfl)")
    hot = turnon_checks(a, prob, base, st, N, 0.25, z_re)

    print("the step caps")
    capd = run(a.exe, prob, os.path.join(a.work, "dt_max"), a.set + ["dt_max_myr=2"])
    sd = steps(capd)
    after = [s for s in sd if s["redshift"][0] < z_re - 1e-9]
    before = [s for s in sd if s["redshift"][0] > z_re + 1e-9]
    check(len(after) > 0 and max(s["dt_yr"][0] for s in after) <= 2e6 * (1 + 1e-9),
          f"dt_max_myr = 2: {len(after)} steps after z_re, the longest {max(s['dt_yr'][0] for s in after) / 1e6:.3f} Myr")
    check(len(before) > 0 and max(s["dt_yr"][0] for s in before) > 2e6, "the steps before z_re are not capped")
    check(matched(st, sd, z_re) < 0.05, f"dt_max_myr: T_M(z) and x_HII,V(z) within 5 per cent of the uncapped run after the turn-on ({matched(st, sd, z_re):.3f})")
    caps = run(a.exe, prob, os.path.join(a.work, "rt_max"), a.set + ["rt_substeps_max=50"])
    ss = steps(caps)
    rt = [s for s in ss if s["plane_emitted"][0] > 0]
    check(len(rt) > 0 and max(s["nstep_rt"][0] for s in rt) <= 50, f"rt_substeps_max = 50: {len(rt)} RT steps, at most {max(s['nstep_rt'][0] for s in rt):.0f} substeps")
    check(sum(1 for s in rt if s["nstep_rt"][0] == 50) > len(rt) // 2, "most RT steps hold exactly 50 substeps")
    check(matched(st, ss, z_re) < 0.05, f"rt_substeps_max: T_M(z) and x_HII,V(z) within 5 per cent of the uncapped run after the turn-on ({matched(st, ss, z_re):.3f})")

    if a.amrex:
        print("the AMReX backend against the serial RT")
        am = run(a.amrex, prob, os.path.join(a.work, "amrex"), a.set + [
            "rt_backend=amrex", "hydro_backend=amrex", "pm_backend=amrex", "amrex_state=1",
            "amrex_max_grid_size=8", "amrex_hydro_max_grid_size=0", "deterministic_parallel=0"], launch)
        sa = steps(am)
        check(len(sa) == len(st), f"the same number of steps ({len(sa)} vs {len(st)})")
        n = min(len(sa), len(st))
        dx = max(abs(sa[i]["xHII_V"][0] - st[i]["xHII_V"][0]) for i in range(n))
        dT = max(abs(sa[i]["T_M"][0] / st[i]["T_M"][0] - 1) for i in range(n)
                 if np.isfinite(st[i]["T_M"][0]))  # T_M is nan before the RT writes its grid
        check(dx < 1e-4 and dT < 1e-4, f"x_HII,V and T_M agree to {dx:.1e}, {dT:.1e} at every step")
        fa, Ta, _ = last_fields(am, N)
        fs, Ts, _ = last_fields(base, N)
        d = np.max(np.abs(fa - fs))
        check(d < 1e-3, f"fHI cell by cell within {d:.1e} at the last step")
        ca = [sa[i]["c_reduced"][0] for i in range(n)]
        check(ca == cs[:n], "the same schedule steps on the backend")
        if hot:
            ah = steps(run(a.amrex, prob, os.path.join(a.work, "amrex_turnon_hot"), a.set + [
                "rt_backend=amrex", "hydro_backend=amrex", "pm_backend=amrex", "amrex_state=1",
                "amrex_max_grid_size=8", "amrex_hydro_max_grid_size=0", "deterministic_parallel=0",
                "zdata1=15 12 10 9 8.5"] + TURNON_SETS, launch))
            # the steps to the turn-on step included bit for bit; the later CFL steps
            # follow the hot gas's cmax_gas, which carries the backend's round-off
            m = min(len(ah), len(hot))
            lt = [i for i, s in enumerate(hot) if abs(s["redshift"][0] - z_re) < 1e-9][0] + 1
            same = m > lt and all(ah[i]["dt_yr"][0] == hot[i]["dt_yr"][0] for i in range(lt + 1))
            ddt = max(abs(ah[i]["dt_yr"][0] / hot[i]["dt_yr"][0] - 1) for i in range(lt + 1, m)) if m > lt + 1 else 0.0
            check(len(ah) == len(hot) and same and ddt < 1e-6,
                  f"the turn-on step on the backend: the serial steps bit for bit to the turn-on step, the later ones to {ddt:.1e} ({len(ah)} vs {len(hot)} steps)")

    if FAILS:
        print(f"\n{len(FAILS)} check(s) failed")
        sys.exit(1)
    print("\nall full-loop plane checks passed")


def dm20(a):
    N = 32
    prob = os.path.join(a.work, "dm20")
    make_problem(prob, N, 0.256, 100, "30 20 15 12 10 9 8.5 8 7.5 7", False,
                 {"plane_ndom": 2, "plane_faces": "+x +y", "plane_gamma12": 0.3, "plane_z_re": 8,
                  "c_reduced_schedule": "0.1 0.99 0.01", "cool_suppress_overdensity": 0})
    out = run(a.exe, prob, os.path.join(a.work, "dm20_run"), a.set)
    st = steps(out)
    print("the scaled-down DM20 problem at z = 7")
    s = st[-1]
    check(abs(s["redshift"][0] - 7) < 0.05, f"the run ended at z = {s['redshift'][0]:.3f}")
    check(s["xHII_V"][0] > 0.99, f"x_HII,V = {s['xHII_V'][0]:.4f} > 0.99")
    fHI, T, rho = last_fields(out, N)
    ion = fHI < 0.5
    TM = np.sum(rho[ion] * T[ion]) / np.sum(rho[ion])
    check(1e4 < TM < 3e4, f"T_M of the ionised gas = {TM:.0f} K in 1e4 .. 3e4 (DM20: ~2e4 K)")
    C = np.mean(rho[ion] ** 2) / np.mean(rho[ion]) ** 2
    check(1 < C < 10, f"clumping factor of the ionised gas <n^2>/<n>^2 = {C:.2f} in 1 .. 10 (DM20 Fig. 4: 2-3 at 1 kpc/h)")
    worst, prev = 0.0, 0.0
    for s in st:
        if s["plane_emitted"][0] > 0:
            res = (s["plane_emitted"][0] - s["plane_absorbed"][0] - s["plane_escaped"][0] -
                   (s["ngam_rays_total"][0] - prev)) / s["plane_emitted"][0]
            worst = max(worst, abs(res))
        prev = s["ngam_rays_total"][0]
    check(worst < 1e-3, f"photon bookkeeping closes to {worst:.1e}")
    print(f"    steps {len(st)}, T_V = {s['T_V'][0]:.0f} K, T_M = {s['T_M'][0]:.0f} K")
    if FAILS:
        print(f"\n{len(FAILS)} check(s) failed")
        sys.exit(1)
    print("\nthe DM20 sanity check passed")


if __name__ == "__main__":
    main()
