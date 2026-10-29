#!/usr/bin/env python3
"""The level Poisson solve (pm_level_solve, Stage G, docs/stage_g.md).

usage: pm_level_solve_check.py WORKDIR LAUNCH "LAUNCH_FLAGS" NPFLAG RADHYDRO SRC

1. Equivalence (the test of pm_gas_levels_check.py with the level solve).
   Uniform gas and the same 16^3 dark matter (make_ics.py, seed 7, 4 Mpc/h,
   z = 100 to 5) on 2 ranks, no RT:
     B   gas 32^3, pm_refine = 1 (the 32^3 mesh): the reference;
     A1  gas 16^3, pm_refine = 2, level 1 over the whole box,
         pm_gas_levels = 1 (the uniform 32^3 mesh fed by the levels);
     G1  gas 16^3, pm_refine = 1, level 1 over the whole box,
         pm_level_solve = 1 (level 0 on the 16^3 mesh by the FFT, level 1
         on its own 32^3 mesh by MLMG, periodic since it covers the box).
   G1's level 1 has B's resolution, mesh and particles, so its composite gas
   (level 0) must match B averaged onto 16^3, and its particles B's, to the
   time stepping of the hierarchy: gas density and momentum rms below 1e-2
   and 2e-2, particle positions below 1e-3 cells rms, and no worse than
   1.5x A1 (which has the same mesh everywhere).
2. Momentum. RADRAY_PM_MOMENTUM_CHECK=1 on the smoke problem with two levels
   and pm_level_solve = 1: over the whole box (every level's solve is
   periodic: the DM + gas momentum change of every PM step below 1e-5 of
   the summed gas |change|, as with the uniform mesh) for pm_refine = 1 and
   2, and over half the box (the Dirichlet level solves break the pairwise
   antisymmetry at the level boundaries: below 2e-2, reported).
"""
import glob
import os
import re
import shutil
import subprocess
import sys

import numpy as np

BASE = ["--quiet", "--set", "legacy_rt=0", "--set", "legacy_physics=0", "--set", "rt_backend=amrex",
        "--set", "hydro_scheme=unsplit", "--set", "hydro_backend=amrex", "--set", "amrex_hydro_max_grid_size=0",
        "--set", "pm_backend=amrex", "--set", "amrex_state=1", "--set", "amrex_max_grid_size=8"]


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


def gas(run_dir, n):
    f = sorted(glob.glob(os.path.join(run_dir, "Gas", "u_gas*")))[-1]
    return np.fromfile(f, "<f4")[: 5 * n ** 3].reshape(5, n, n, n, order="F").astype(float)


def xv(run_dir, npart):
    f = sorted(glob.glob(os.path.join(run_dir, "DM", "xv_dm*")))[-1]
    return np.fromfile(f, "<f4")[: 6 * npart].reshape(-1, 6)


def coarse(a):
    n = a.shape[0] // 2
    return a.reshape(n, 2, n, 2, n, 2).mean(axis=(1, 3, 5))


def main():
    work, launch, flags, npflag, exe, src = sys.argv[1:7]
    flags = [f for f in flags.replace(";", " ").split() if f]
    if launch.endswith("srun") and not any("cpus-per-task" in f for f in flags):
        flags.append("--cpus-per-task=" + os.environ.get("OMP_NUM_THREADS", "1"))
    mpi = [launch] + flags + [npflag, "2"]
    if os.path.isdir(work):
        shutil.rmtree(work)
    os.makedirs(work)
    fails = []

    # 1. equivalence
    smoke = open(os.path.join(src, "tests", "problems_full", "smoke", "params.txt")).read()
    for n in (16, 32):
        p = os.path.join(work, f"p{n}")
        run([sys.executable, os.path.join(src, "tools", "ics", "make_ics.py"), os.path.join(p, "ics"),
             "--kind", "zeldovich", "--Ngas", str(n), "--Ndm", "16", "--z-init", "100", "--box", "4",
             "--seed", "7", "--gas", "uniform"], os.path.join(work, f"ics{n}.log"))
        txt = re.sub(r"(?m)^N = .*$", f"N = {n}", smoke)
        txt = re.sub(r"(?m)^Ndm = .*$", "Ndm = 16", txt)
        txt = re.sub(r"(?m)^zdata1 = .*$", "zdata1 = 5", txt)
        open(os.path.join(p, "params.txt"), "w").write(txt)
    whole = sets("amrex_max_level=1", "amrex_refine_mass=1e9", "amrex_refine_box=0,0,0,15,15,15")
    runs = {"B": ("p32", sets("pm_refine=1")),
            "A1": ("p16", sets("pm_refine=2", "pm_gas_levels=1") + whole),
            "G1": ("p16", sets("pm_refine=1", "pm_level_solve=1") + whole)}
    for name, (p, extra) in runs.items():
        run(mpi + [exe, os.path.join(work, p), os.path.join(work, name)] + BASE + sets("run_steps=60") + extra,
            os.path.join(work, name + ".log"))
    b = gas(os.path.join(work, "B"), 32)
    rb = coarse(b[0])
    pb = [coarse(b[1 + m]) / 2 for m in range(3)]  # B's momenta in 16^3-cell velocity units
    xb = xv(os.path.join(work, "B"), 16 ** 3)
    err = {}
    for name in ("A1", "G1"):
        a = gas(os.path.join(work, name), 16)
        e_rho = float(np.sqrt(np.mean((a[0] / rb - 1) ** 2)))
        e_p = float(np.sqrt(sum(np.mean((a[1 + m] - pb[m]) ** 2) for m in range(3)) /
                            sum(np.mean(pb[m] ** 2) for m in range(3))))
        dx = xv(os.path.join(work, name), 16 ** 3)[:, :3].astype(float) - xb[:, :3]
        dx -= 16 * np.round(dx / 16)
        e_x = float(np.sqrt(np.mean(dx ** 2)))
        err[name] = (e_rho, e_p, e_x)
        print(f"{name}: against B (32^3), gas rho rms {e_rho:.3e}, gas momentum rms {e_p:.3e}, "
              f"particle positions rms {e_x:.3e} cells")
    eg, ea = err["G1"], err["A1"]
    for i, (lab, lim) in enumerate((("gas rho", 1e-2), ("gas momentum", 2e-2), ("particles", 1e-3))):
        if not eg[i] < lim:
            fails.append(f"G1 {lab} {eg[i]:.3e} >= {lim:g}")
        if not eg[i] < 1.5 * ea[i] + 1e-5:
            fails.append(f"G1 {lab} {eg[i]:.3e} worse than 1.5x A1's {ea[i]:.3e}")

    # 2. momentum conservation with two levels
    prob = os.path.join(src, "tests", "problems_full", "smoke")
    env = dict(os.environ, RADRAY_PM_MOMENTUM_CHECK="1")
    for r, box, lim in ((1, "0,0,0,15,15,15", 1e-5), (2, "0,0,0,15,15,15", 1e-5), (1, "0,0,0,7,7,7", 2e-2)):
        tag = "whole" if box.endswith("15") else "half"
        o = os.path.join(work, f"mom_r{r}_{tag}")
        log = o + ".log"
        run(mpi + [exe, prob, o] + BASE + sets("run_steps=6", f"pm_refine={r}", "amrex_max_level=2",
                                               "amrex_refine_mass=1e9", f"amrex_refine_box={box}",
                                               "pm_level_solve=1"), log, env)
        vals = [(int(m.group(1)), float(m.group(2))) for m in
                re.finditer(r"PM momentum check: levels (\d+), .*? = (\S+) ", open(log).read())]
        two = [v for lev, v in vals if lev == 2]
        worst = max(two) if two else float("nan")
        print(f"momentum, pm_refine = {r}, levels over the {tag} box: {len(two)} PM steps with two levels, "
              f"worst {worst:.3e}")
        if not two:
            fails.append(f"pm_refine = {r}, {tag}: no PM step with two levels")
        elif not worst < lim:
            fails.append(f"pm_refine = {r}, {tag}: momentum change {worst:.3e} >= {lim:g}")
    if fails:
        print("FAIL:\n  " + "\n  ".join(fails))
        sys.exit(1)
    print("PASS")


if __name__ == "__main__":
    main()
