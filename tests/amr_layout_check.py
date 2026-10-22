#!/usr/bin/env python3
"""Two refined levels are independent of the rank count and of the fine
levels' box layout (docs/stage_f.md section 3): AMReX chops the fine grids
by the rank count, and the RT hierarchy once rounded each fine box's partial
sum of the deposits it hands to the parent, so 1, 2 and 3 ranks differed
from about step 8 (relative 1e-10).

usage: amr_layout_check.py WORK LAUNCH FLAGS NPFLAG DRIVER SRC

Makes the 8-source problem of tests/problems_full/make_ifront.py at 16^3,
runs it with two nested refined levels on 1 rank (fine boxes of 16 and of
8 cells) and on 2 and 3 ranks, and requires every output file (fields,
stats, snapshots) to be byte-identical to the first run.
"""
import filecmp
import os
import subprocess
import sys

N = 16
BOXES = "0,2,2,7,9,9,1,3,3,6,8,8"  # level 1 and level 2 (base cells, inclusive)


def main():
    work, launch, flags, npflag, drv, src = sys.argv[1:7]
    flags = [f for f in flags.replace(";", " ").split() if f]
    if launch.endswith("srun") and not any("cpus-per-task" in f for f in flags):
        flags.append("--cpus-per-task=" + os.environ.get("OMP_NUM_THREADS", "1"))
    prob = os.path.join(work, "prob")
    subprocess.check_call([sys.executable, os.path.join(src, "tests/problems_full/make_ifront.py"), prob,
                           "--N", str(N), "--box", "2", "--sources", "8src"], stdout=subprocess.DEVNULL)
    common = ["--quiet", "--set", "run_steps=8", "--set", "amrex_max_level=2", "--set", "amrex_refine_mass=1e9",
              "--set", "amrex_n_error_buf=0", "--set", "amrex_refine_box=" + BOXES,
              "--set", "amrex_max_grid_size=8", "--set", "fields_every=1"]
    cases = [("r1_fmg16", 1, ["--set", "amrex_fine_max_grid_size=16"]),
             ("r1_fmg8", 1, ["--set", "amrex_fine_max_grid_size=8"]),
             ("r2", 2, []), ("r3", 3, [])]
    outs = []
    for name, npr, extra in cases:
        out = os.path.join(work, name)
        cmd = ([launch] + flags + [npflag, str(npr)] if (launch and npr > 1) else []) + [drv, prob, out] + common + extra
        subprocess.check_call(cmd)
        outs.append((name, out))
    ref_name, ref = outs[0]
    files = []
    for root, _, names in os.walk(ref):
        for n in names:
            if n.endswith((".bin", ".txt")) or n.startswith(("u_gas", "rt_grid")):
                rel = os.path.relpath(os.path.join(root, n), ref)
                if not rel.startswith("Log") and n != "timing.txt":
                    files.append(rel)
    bad = 0
    for name, out in outs[1:]:
        diff = [f for f in files if not (os.path.exists(os.path.join(out, f)) and
                                         filecmp.cmp(os.path.join(ref, f), os.path.join(out, f), shallow=False))]
        print(f"{name} against {ref_name}: {len(files)} files, {len(diff)} differ" + (f" (first: {diff[0]})" if diff else ""))
        bad += len(diff)
    if not files:
        print("FAIL: no output files")
        sys.exit(1)
    print("PASS" if not bad else "FAIL")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
