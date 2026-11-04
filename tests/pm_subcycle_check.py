#!/usr/bin/env python3
"""Momentum with the sub-cycled gravity of the refined levels
(pm_level_subcycle, docs/stage_g.md section 6).

usage: pm_subcycle_check.py WORKDIR LAUNCH "LAUNCH_FLAGS" NPFLAG RADHYDRO SRC

RADRAY_PM_MOMENTUM_CHECK=1 on the smoke problem (2 ranks, 6 steps) with two
refined levels and pm_level_solve = 1, with pm_level_subcycle = 0 (the
momentum of each PM step, docs/stage_g.md section 3) and 1 (the momentum
all the kicks of a base step add: level 0's at the step's middle and every
refined level's in its substeps, over the gas cells not covered by a finer
level and the particles), against the summed |dp| of the gas kicks:
  * levels 1 and 2 over the whole box (every solve periodic, every particle
    on level 2: the kicks are pairwise antisymmetric): below 1e-5;
  * over half the box (the Dirichlet boundary, and with sub-cycling the
    levels' mass and kicks at different times): below 2e-2, and reported
    next to pm_level_subcycle = 0's.
Then the corrector diagnostic (RADRAY_PM_CORRECTOR_CHECK=1, docs/stage_g.md
section 6.7) on the half-box run: it reports, and the outputs are
byte-identical to the run without it.
"""
import os
import re
import shutil
import subprocess
import sys

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


def main():
    work, launch, flags, npflag, exe, src = sys.argv[1:7]
    flags = [f for f in flags.replace(";", " ").split() if f]
    if launch.endswith("srun") and not any("cpus-per-task" in f for f in flags):
        flags.append("--cpus-per-task=" + os.environ.get("OMP_NUM_THREADS", "1"))
    mpi = [launch] + flags + [npflag, "2"]
    if os.path.isdir(work):
        shutil.rmtree(work)
    os.makedirs(work)
    prob = os.path.join(src, "tests", "problems_full", "smoke")
    env = dict(os.environ, RADRAY_PM_MOMENTUM_CHECK="1")
    fails = []
    for box, lim in (("0,0,0,15,15,15", 1e-5), ("0,0,0,7,7,7", 2e-2)):
        tag = "whole" if box.endswith("15") else "half"
        worst = {}
        for sub in (0, 1):
            o = os.path.join(work, f"mom_{tag}_sub{sub}")
            log = o + ".log"
            run(mpi + [exe, prob, o] + BASE + sets("run_steps=6", "pm_refine=1", "amrex_max_level=2",
                                                   "amrex_refine_mass=1e9", f"amrex_refine_box={box}",
                                                   "pm_level_solve=1", f"pm_level_subcycle={sub}"), log, env)
            txt = open(log).read()
            if sub:
                pat = r"PM momentum check \(sub-cycled\): levels (\d+), \|sum dP of the kicks\| / sum\|dp_gas\| = (\S+),"
            else:
                pat = r"PM momentum check: levels (\d+), .*? = (\S+) "
            two = [float(m.group(2)) for m in re.finditer(pat, txt) if int(m.group(1)) == 2]
            worst[sub] = max(two) if two else float("nan")
            print(f"momentum, levels over the {tag} box, pm_level_subcycle = {sub}: {len(two)} steps with two "
                  f"levels, worst {worst[sub]:.3e}")
            if not two:
                fails.append(f"{tag}, pm_level_subcycle = {sub}: no step with two levels")
        if not worst[1] < lim:
            fails.append(f"{tag}: sub-cycled momentum change {worst[1]:.3e} >= {lim:g}")
    # the corrector diagnostic (RADRAY_PM_CORRECTOR_CHECK=1) reports and
    # changes nothing: the half-box sub-cycled run with and without it
    outs = {}
    for chk in (0, 1):
        o = os.path.join(work, f"corr{chk}")
        e = dict(os.environ)
        e.pop("RADRAY_PM_CORRECTOR_CHECK", None)
        if chk:
            e["RADRAY_PM_CORRECTOR_CHECK"] = "1"
        run(mpi + [exe, prob, o] + BASE + sets("run_steps=6", "pm_refine=1", "amrex_max_level=2", "amrex_refine_mass=1e9",
                                               "amrex_refine_box=0,0,0,7,7,7", "pm_level_solve=1",
                                               "pm_level_subcycle=1", "fields_every=1"), o + ".log", e)
        outs[chk] = o
    rep = [float(m.group(1)) for m in
           re.finditer(r"PM corrector check: level \d+ .*? force near the children (\S+),", open(outs[1] + ".log").read())]
    print(f"corrector check: {len(rep)} reports, force near the children at most "
          f"{max(rep) if rep else float('nan'):.3e} of max |F|")
    if not rep or not all(0 <= v < 1 for v in rep):
        fails.append("corrector check: no or bad reports")
    ndiff = 0
    for root, _, names in os.walk(outs[0]):
        for nm in names:
            if nm == "timing.txt":
                continue
            a = os.path.join(root, nm)
            b = os.path.join(outs[1], os.path.relpath(a, outs[0]))
            if not os.path.isfile(b) or open(a, "rb").read() != open(b, "rb").read():
                ndiff += 1
    print(f"corrector check: {ndiff} files differ from the run without it")
    if ndiff:
        fails.append(f"corrector check changed {ndiff} output files")
    if fails:
        print("FAIL:\n  " + "\n  ".join(fails))
        sys.exit(1)
    print("PASS")


if __name__ == "__main__":
    main()
