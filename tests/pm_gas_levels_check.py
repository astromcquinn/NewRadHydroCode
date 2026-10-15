#!/usr/bin/env python3
"""The PM on the refined gas levels (pm_gas_levels, docs/stage_f.md section 4).

usage: pm_gas_levels_check.py WORKDIR LAUNCH "LAUNCH_FLAGS" NPFLAG RADHYDRO SRC

1. Equivalence. Uniform gas and the same 16^3 dark matter (make_ics.py,
   seed 7, 4 Mpc/h, z = 100 to 5) in four runs on 2 ranks, no RT:
     B   gas 32^3, pm_refine = 1 (the 32^3 mesh): the reference;
     C   gas 16^3, pm_refine = 2 (the same mesh), no refinement;
     A0  C with level 1 over the whole box, pm_gas_levels = 0 (the level-0
         coupling: the fine cells get their parent's kick);
     A1  the same with pm_gas_levels = 1.
   A1's level 1 has B's resolution, mesh and particles, so its composite
   gas (level 0) must match B averaged onto 16^3, and its particles B's,
   to the time-stepping differences of the hierarchy: gas density and
   momentum rms below 1e-2 and 2e-2, particle positions below 1e-3 cells
   rms, each at least 10x better than A0.
2. Momentum. RADRAY_PM_MOMENTUM_CHECK=1 on the smoke problem with two levels
   over the whole box and pm_refine = 1, 2, 4 (levels finer than, mixed
   with, and coarser than the mesh): the DM + gas momentum change of every
   PM step below 1e-5 of the summed gas |change|.
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
    runs = {"B": ("p32", sets("pm_refine=1")), "C": ("p16", sets("pm_refine=2")),
            "A0": ("p16", sets("pm_refine=2", "pm_gas_levels=0") + whole),
            "A1": ("p16", sets("pm_refine=2", "pm_gas_levels=1") + whole)}
    for name, (p, extra) in runs.items():
        run(mpi + [exe, os.path.join(work, p), os.path.join(work, name)] + BASE + sets("run_steps=60") + extra,
            os.path.join(work, name + ".log"))
    b = gas(os.path.join(work, "B"), 32)
    rb = coarse(b[0])
    pb = [coarse(b[1 + m]) / 2 for m in range(3)]  # B's momenta in 16^3-cell velocity units
    xb = xv(os.path.join(work, "B"), 16 ** 3)
    err = {}
    for name in ("C", "A0", "A1"):
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
    e1, e0 = err["A1"], err["A0"]
    for i, (lab, lim) in enumerate((("gas rho", 1e-2), ("gas momentum", 2e-2), ("particles", 1e-3))):
        if not e1[i] < lim:
            fails.append(f"A1 {lab} {e1[i]:.3e} >= {lim:g}")
        if not e1[i] * 10 < e0[i]:
            fails.append(f"A1 {lab} {e1[i]:.3e} not 10x below A0's {e0[i]:.3e}")

    # 2. momentum conservation with two levels
    prob = os.path.join(src, "tests", "problems_full", "smoke")
    env = dict(os.environ, RADRAY_PM_MOMENTUM_CHECK="1")
    for r in (1, 2, 4):
        o = os.path.join(work, f"mom_r{r}")
        log = o + ".log"
        run(mpi + [exe, prob, o] + BASE + sets("run_steps=6", f"pm_refine={r}", "amrex_max_level=2",
                                               "amrex_refine_mass=1e9", "amrex_refine_box=0,0,0,15,15,15"),
            log, env)
        vals = [(int(m.group(1)), float(m.group(2))) for m in
                re.finditer(r"PM momentum check: levels (\d+), .*? = (\S+) ", open(log).read())]
        two = [v for lev, v in vals if lev == 2]
        worst = max(two) if two else float("nan")
        print(f"momentum, pm_refine = {r}: {len(two)} PM steps with two levels, worst {worst:.3e}")
        if not two:
            fails.append(f"pm_refine = {r}: no PM step with two levels")
        elif not worst < 1e-5:
            fails.append(f"pm_refine = {r}: momentum change {worst:.3e} >= 1e-5")
    if fails:
        print("FAIL:\n  " + "\n  ".join(fails))
        sys.exit(1)
    print("PASS")


if __name__ == "__main__":
    main()
