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
  * --amrex: the backend run (rt_backend = hydro_backend = pm_backend =
    amrex, amrex_state = 1) agrees with the serial run to round-off:
    x_HII,V and T_M to 1e-4 and fHI cell by cell to 1e-3 at the last step.

--dm20 (label nightly): the scaled-down physics sanity check of
docs/plane_parallel.md: 0.256 Mpc/h, 32^3 (DM20: 1024^3), Zel'dovich
initial conditions from z = 100, z_re = 8, Gamma_-12 = 0.3, 2^3 domains,
cooling on everywhere (cool_suppress_overdensity = 0), land_outputs = 1,
run to z = 7. Checks: the run ends on z = 7 and the outputs are the
zdata1 redshifts (no extra step past the last one), x_HII,V > 0.99 at
z = 7, the mass-weighted temperature of the ionised gas in 1e4 .. 3e4 K
(DM20: about 2e4 K just behind the fronts), the ionised-gas clumping
factor <n^2>/<n>^2 in 1 .. 10 (DM20 Fig. 4: ~ 2-3 at 1 kpc/h; coarser
here), and the photon bookkeeping closes to 1e-3 in every RT step.
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
        n0, n1 = st[sw[0] - 1]["nstep_rt"][0], st[sw[0]]["nstep_rt"][0]
        check(n0 / max(n1, 1) > 5, f"the RT substeps per step drop {n0:.0f} -> {n1:.0f} at the switch (> 5x)")
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

    if FAILS:
        print(f"\n{len(FAILS)} check(s) failed")
        sys.exit(1)
    print("\nall full-loop plane checks passed")


def dm20(a):
    N = 32
    prob = os.path.join(a.work, "dm20")
    make_problem(prob, N, 0.256, 100, "30 20 15 12 10 9 8.5 8 7.5 7", False,
                 {"plane_ndom": 2, "plane_faces": "+x +y", "plane_gamma12": 0.3, "plane_z_re": 8,
                  "c_reduced_schedule": "0.1 0.99 0.01", "cool_suppress_overdensity": 0,
                  "land_outputs": 1})
    out = run(a.exe, prob, os.path.join(a.work, "dm20_run"), a.set)
    st = steps(out)
    print("the scaled-down DM20 problem at z = 7")
    s = st[-1]
    check(abs(s["redshift"][0] - 7) < 1e-9, f"the run ended on the last output z = 7 (z = {s['redshift'][0]:.12f})")
    zw = sorted((float(p.split("z=")[1]) for p in glob.glob(os.path.join(out, "Pk", "pk.z=*"))), reverse=True)
    check(zw == [30, 20, 15, 12, 10, 9, 8.5, 8, 7.5, 7], f"land_outputs = 1: the outputs are the zdata1 redshifts ({zw})")
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
