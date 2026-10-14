#!/usr/bin/env python3
"""A uniform box with one source at its centre, for the RT with refinement
(docs/stage_f.md section 3): the I-front of a single source in uniform gas
crossing a refinement boundary.

usage: make_ifront.py <outdir> [--N 32] [--box 4] [--L 8e27]

Writes <outdir>/{params.txt, ics/, Halo/}: uniform gas (tools/ics/make_ics.py
--gas uniform --growth frozen-gas, the dark matter at rest on the grid), the
source as a legacy halo catalogue (halo.z=..., sfr.txt, as make_problems.py's
smoke problem) at the box centre, RT from z = 30 with the fixed physics and
the AMReX backends with the gas state on the levels. Refinement is set by the
test: amrex_max_level and amrex_refine_box. Analysis: tools/analysis/ifront_full.py.
"""
import argparse
import os
import shutil
import struct
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))


def tag(z):
    return f"z=0{z:6.4f}" if z < 10 else f"z={z:7.4f}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir")
    ap.add_argument("--N", type=int, default=32)
    ap.add_argument("--box", type=float, default=4.0)
    ap.add_argument("--L", type=float, default=8e27, help="source luminosity record field (photons/s)")
    a = ap.parse_args()
    out = a.outdir
    os.makedirs(out, exist_ok=True)
    subprocess.check_call([sys.executable, os.path.join(ROOT, "tools/ics/make_ics.py"),
                           os.path.join(out, "ics"), "--kind", "zeldovich", "--gas", "uniform",
                           "--growth", "frozen-gas", "--Ngas", str(a.N), "--Ndm", str(a.N // 2),
                           "--z-init", "100", "--box", str(a.box), "--seed", "7"])
    with open(os.path.join(out, "params.txt"), "w") as f:
        f.write(f"""# One source in uniform gas (make_ifront.py): the RT with refinement
N = {a.N}
Ndm = {a.N // 2}
box = {a.box}
nray_max_factor = 16
z_init = 100
zdata1 = 30 28 26 24 22 20
ics = ics
source_mode = legacy
dir_halo = Halo
fof_min_npart = 8
legacy_rt = 0
legacy_physics = 0
c_reduced = 0.1
hydro_scheme = unsplit
rt_backend = amrex
hydro_backend = amrex
pm_backend = amrex
amrex_state = 1
amrex_max_grid_size = 16
amrex_hydro_max_grid_size = 0
deterministic_parallel = 0
fields = fHI,temp
fields_every = 1
restart_every = 0
""")
    sh = os.path.join(out, "Halo")
    shutil.rmtree(sh, ignore_errors=True)
    os.makedirs(sh)
    c = a.box / 2
    # source_type: x(3),v(3) real(4), N int(4), M,R,T,L,sfr,mstar,fgas real(4)
    rec = struct.pack("<6fi7f", c, c, c, 0, 0, 0, 40, 1e10, 0.02, 3e4, a.L, 0.5, 0, 0)
    for zc in (40.0, 10.0):
        with open(os.path.join(sh, "halo." + tag(zc)), "wb") as f:
            f.write(struct.pack("<q", 1) + rec)
    with open(os.path.join(sh, "sfr.txt"), "w") as f:
        for zc, sfr in ((40.0, 0.01), (10.0, 0.01)):
            f.write(f"{zc:10.4f}" + " 0" * 8 + f" {sfr:.6e}\n")
    print("wrote", out)


if __name__ == "__main__":
    main()
