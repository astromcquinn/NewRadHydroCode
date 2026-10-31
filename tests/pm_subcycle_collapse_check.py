#!/usr/bin/env python3
"""A collapsing overdensity on two refined levels, with and without the
sub-cycled gravity of the refined levels (pm_level_subcycle,
docs/stage_g.md section 6), against the uniform run at the finest
resolution with a short time step.

usage: pm_subcycle_collapse_check.py WORKDIR LAUNCH "LAUNCH_FLAGS" NPFLAG NP RADHYDRO SRC
                                     [--base 16] [--fine 64] [--zc 12] [--zend 7] [--sigma 3]

The particles (fine^3, a lattice in a 4 Mpc/h box from z = 100) carry one
spherical Gaussian Zel'dovich displacement towards a point near the box's
centre, psi = -A (q - c) exp(-|q - c|^2 / 2 sigma^2), whose linear central
overdensity 3A reaches 1.686 at z = ZC; the gas starts uniform (no RT):
  A   gas fine^3, no refinement, the step capped at 1/8 of B1's median base
      step (dt_max_myr): the reference;
  B0  gas base^3 with levels 1 and 2 over the central half (level 2 at the
      fine resolution there), pm_level_solve = 1, pm_level_subcycle = 0
      (the levels kicked once per base step);
  B1  the same with pm_level_subcycle = 1 (level 2 kicked 8 times per base
      step).
At ZEND (after the collapse): the halo's centre (the densest particle's
neighbourhood: the centre of mass within 2 fine cells of A's, iterated), its
mass M200 (the particles within the radius where the enclosed mean density
falls to 200 times the mean) and its spherically averaged density profile
(fine-cell shells 0.5-12) for each run. Passes when B1's centre lies within
0.5 fine cells of A's, its M200 within 10 % of A's, its profile within 25 %
rms (log) of A's, and B1 is not worse than B0 in centre, mass and profile
(each by a margin of 10 % of B0's error or 1e-3).
"""
import argparse
import os
import re
import shutil
import subprocess
import sys

import numpy as np

BASE = ["--quiet", "--set", "legacy_rt=0", "--set", "legacy_physics=0", "--set", "rt_backend=amrex",
        "--set", "hydro_scheme=unsplit", "--set", "hydro_backend=amrex", "--set", "amrex_hydro_max_grid_size=0",
        "--set", "pm_backend=amrex", "--set", "amrex_state=1", "--set", "amrex_max_grid_size=8",
        "--set", "deterministic_parallel=0", "--set", "tvir_min=1e30", "--set", "movie_dz=0"]


def run(cmd, log, env=None):
    with open(log, "w") as f:
        r = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, env=env)
    if r.returncode != 0:
        print(open(log).read()[-3000:])
        sys.exit(f"failed: {' '.join(cmd)}")


def sets(*kv):
    out = []
    for x in kv:
        out += ["--set", x]
    return out


def last_dm(run_dir, npart):
    d = sorted(x for x in os.listdir(os.path.join(run_dir, "DM")) if x.startswith("xv_dm"))[-1]
    return np.fromfile(os.path.join(run_dir, "DM", d), "<f4")[: 6 * npart].reshape(-1, 6).astype(float)


def halo(x, n, c0):
    """centre (iterated centre of mass within 2 cells), M200 (particles) and
    the density profile (mean 1) around it"""
    c = np.array(c0, float)
    for _ in range(20):
        d = x - c
        d -= n * np.round(d / n)
        r = np.sqrt((d ** 2).sum(axis=1))
        sel = r < 2.0
        if sel.sum() == 0:
            break
        c2 = c + d[sel].mean(axis=0)
        if np.abs(c2 - c).max() < 1e-4:
            c = c2
            break
        c = c2
    d = x - c
    d -= n * np.round(d / n)
    r = np.sort(np.sqrt((d ** 2).sum(axis=1)))
    k = np.arange(1, r.size + 1)
    over = k / (4 / 3 * np.pi * np.maximum(r, 1e-6) ** 3)  # mean density within r (particles per cell; mean 1)
    i = np.nonzero(over < 200)[0]
    m200 = int(i[0]) if i.size else r.size
    edges = np.concatenate(([0.5], np.geomspace(1.0, 12.0, 9)))
    prof = np.array([((r >= a) & (r < b)).sum() / (4 / 3 * np.pi * (b ** 3 - a ** 3)) for a, b in zip(edges[:-1], edges[1:])])
    return c, m200, prof


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    pos = ("work", "launch", "flags", "npflag", "np", "exe", "src")
    if len(sys.argv) < 8:
        ap.error("usage: pm_subcycle_collapse_check.py " + " ".join(p.upper() for p in pos) + " [options]")
    ap.add_argument("--base", type=int, default=16)
    ap.add_argument("--fine", type=int, default=64)
    ap.add_argument("--zc", type=float, default=12.0)
    ap.add_argument("--zend", type=float, default=7.0)
    ap.add_argument("--sigma", type=float, default=3.0)
    ap.add_argument("--keep", action="store_true", help="reuse existing runs")
    a = ap.parse_args(sys.argv[8:])
    for name, val in zip(pos, sys.argv[1:8]):
        setattr(a, name, val)
    flags = [f for f in a.flags.replace(";", " ").split() if f]
    if a.launch.endswith("srun") and not any("cpus-per-task" in f for f in flags):
        flags.append("--cpus-per-task=" + os.environ.get("OMP_NUM_THREADS", "1"))
    mpi = [a.launch] + flags + [a.npflag, a.np]
    work = a.work
    if os.path.isdir(work) and not a.keep:
        shutil.rmtree(work)
    os.makedirs(work, exist_ok=True)
    nb, nf = a.base, a.fine
    nlev = (nf // nb).bit_length() - 1
    c0 = np.array([nf / 2 + 0.3, nf / 2 - 0.3, nf / 2 + 0.1])  # the collapse centre (fine cells = DM cells)

    smoke = open(os.path.join(a.src, "tests", "problems_full", "smoke", "params.txt")).read()
    for n in (nb, nf):
        p = os.path.join(work, f"p{n}")
        if a.keep and os.path.isfile(os.path.join(p, "params.txt")):
            continue
        ics = os.path.join(p, "ics")
        # a pancake's ICs for the lattice, the gas and the growth factors (v = vamp/amp psi)
        run([sys.executable, os.path.join(a.src, "tools", "ics", "make_ics.py"), ics, "--kind", "pancake",
             "--Ngas", str(n), "--Ndm", str(nf), "--z-init", "100", "--box", "4", "--gas", "uniform",
             "--a-cross", repr(1 / (1 + a.zc))], os.path.join(work, f"ics{n}.log"))
        info = dict(re.findall(r"(?m)^(\w+) = (\S+)$", open(os.path.join(ics, "ics_info.txt")).read()))
        amp, vamp = float(info["D_over_Dcross_init"]), float(info["vamp"])
        q = np.arange(nf) + 0.5
        Q = np.stack(np.meshgrid(q, q, q, indexing="ij")).reshape(3, -1, order="F")
        d = Q - c0[:, None]
        d -= nf * np.round(d / nf)
        A = 1.686 / 3 * amp  # the central linear overdensity 3A = 1.686 D/D(a_c)
        psi = -A * d * np.exp(-(d ** 2).sum(axis=0) / (2 * a.sigma ** 2))
        xv = np.zeros((6, nf ** 3))
        xv[:3] = np.mod(Q + psi, nf)
        xv[3:] = vamp / amp * psi
        xv.astype("<f4").ravel(order="F").tofile(os.path.join(ics, "xv_dm.init"))
        txt = re.sub(r"(?m)^N = .*$", f"N = {n}", smoke)
        txt = re.sub(r"(?m)^Ndm = .*$", f"Ndm = {nf}", txt)
        txt = re.sub(r"(?m)^box = .*$", "box = 4", txt)
        txt = re.sub(r"(?m)^zdata1 = .*$", f"zdata1 = {a.zend:g}", txt)
        open(os.path.join(p, "params.txt"), "w").write(txt)
    q4 = nb // 4
    region = f"{q4},{q4},{q4},{3 * q4 - 1},{3 * q4 - 1},{3 * q4 - 1}"
    levels = sets(f"amrex_max_level={nlev}", "amrex_refine_mass=1e9", f"amrex_refine_box={region}",
                  "amrex_n_error_buf=0", "pm_refine=1", "pm_level_solve=1")
    runs = {"B0": (f"p{nb}", levels + sets("pm_level_subcycle=0")),
            "B1": (f"p{nb}", levels + sets("pm_level_subcycle=1"))}
    for name, (p, extra) in runs.items():
        o = os.path.join(work, name)
        if not (a.keep and os.path.isdir(os.path.join(o, "DM"))):
            run(mpi + [a.exe, os.path.join(work, p), o] + BASE + sets("run_steps=2000", "step_timers=1") + extra,
                os.path.join(work, name + ".log"))
    # A's step: 1/8 of B1's median base step (level 2's substep)
    dts = [float(m.group(1)) for m in re.finditer(r"(?m)^dt_yr\s+(\S+)", "".join(
        open(os.path.join(work, "B1", s, "stats.txt")).read() for s in sorted(os.listdir(os.path.join(work, "B1")))
        if s.startswith("step_") and os.path.isfile(os.path.join(work, "B1", s, "stats.txt"))))]
    dts = [x for x in dts if x > 0]
    cap = float(np.median(dts)) / 8 / 1e6
    print(f"B1: {len(dts)} base steps, median {np.median(dts) / 1e6:.3f} Myr; A's step capped at {cap:.4f} Myr")
    o = os.path.join(work, "A")
    if not (a.keep and os.path.isdir(os.path.join(o, "DM"))):
        run(mpi + [a.exe, os.path.join(work, f"p{nf}"), o] + BASE +
            sets("run_steps=20000", "pm_refine=1", "amrex_max_level=0", f"dt_max_myr={cap!r}"),
            os.path.join(work, "A.log"))
    npart = nf ** 3
    res = {}
    for name in ("A", "B0", "B1"):
        x = last_dm(os.path.join(work, name), npart)[:, :3]
        res[name] = halo(x, nf, c0 if name == "A" else res["A"][0])
    cA, mA, pA = res["A"]
    print(f"A: centre ({cA[0]:.3f}, {cA[1]:.3f}, {cA[2]:.3f}), M200 {mA} particles; profile " +
          " ".join(f"{v:.3g}" for v in pA))
    err = {}
    for name in ("B0", "B1"):
        c, m, p = res[name]
        dc = c - cA
        dc -= nf * np.round(dc / nf)
        ok = (p > 0) & (pA > 0)
        lp = float(np.sqrt(np.mean(np.log(p[ok] / pA[ok]) ** 2)))
        err[name] = (float(np.sqrt((dc ** 2).sum())), abs(m / mA - 1), lp)
        print(f"{name}: centre offset {err[name][0]:.3f} fine cells, M200 {m} ({m / mA - 1:+.3f}), profile rms log "
              f"difference {lp:.3f}; profile " + " ".join(f"{v:.3g}" for v in p))
    fails = []
    e0, e1 = err["B0"], err["B1"]
    for i, (lab, lim) in enumerate((("centre", 0.5), ("M200", 0.10), ("profile", 0.25))):
        if not e1[i] < lim:
            fails.append(f"B1 {lab} {e1[i]:.3f} >= {lim:g}")
        if not e1[i] <= 1.1 * e0[i] + 1e-3:
            fails.append(f"B1 {lab} {e1[i]:.3f} worse than B0's {e0[i]:.3f}")
    if fails:
        print("FAIL:\n  " + "\n  ".join(fails))
        sys.exit(1)
    print("PASS")


if __name__ == "__main__":
    main()
