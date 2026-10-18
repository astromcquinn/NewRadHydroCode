#!/usr/bin/env python3
"""A uniform box with one source at its centre, or the 8 sources of
tests/problems/uniform_8src, for the RT with refinement (docs/stage_f.md
section 3): the I-fronts in uniform gas crossing a refinement boundary.

usage: make_ifront.py <outdir> [--N 32] [--box 4] [--L 8e27] [--sources center|8src|<sources.txt>]

--sources 8src (or a sources.txt of the RT validation problems, `i j k Luv sfr`
with 1-based cell indices on the grid of the params.txt next to it): every
source at the centre of its cell on that grid, i.e. at the same place in the
box for any N, with the luminosity L * Luv/max(Luv).

Writes <outdir>/{params.txt, ics/, Halo/}: uniform gas (tools/ics/make_ics.py
--gas uniform --growth frozen-gas) and the dark matter at rest on its lattice
(written over make_ics.py's Zel'dovich particles), so nothing moves under
gravity; the sources as a legacy halo catalogue (halo.z=..., sfr.txt, as
make_problems.py's smoke problem) and their positions in box units
(sources_pos.txt: x y z L); RT from z = 30 with the fixed physics and the
AMReX backends with the gas state on the levels. Refinement is set by the
test: amrex_max_level and amrex_refine_box. Analysis:
tools/analysis/ifront_full.py (one source), ifront_multi.py (several).
"""
import argparse
import os
import shutil
import struct
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))


def tag(z):
    return f"z=0{z:6.4f}" if z < 10 else f"z={z:7.4f}"


def read_sources(spec, box):
    """[(x, y, z, L/Lmax)] in Mpc/h"""
    if spec == "center":
        return [(box / 2, box / 2, box / 2, 1.0)]
    path = os.path.join(ROOT, "tests/problems/uniform_8src/sources.txt") if spec == "8src" else spec
    n = None
    with open(os.path.join(os.path.dirname(os.path.abspath(path)), "params.txt")) as f:
        for line in f:
            k, _, v = line.partition("=")
            if k.strip() == "N":
                n = int(v)
    if n is None:
        raise SystemExit(f"{path}: no N in the params.txt next to it")
    rows = [line.split() for line in open(path) if line.strip()]
    lmax = max(float(r[3]) for r in rows)
    return [tuple((int(r[m]) - 0.5) * box / n for m in range(3)) + (float(r[3]) / lmax,) for r in rows]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir")
    ap.add_argument("--N", type=int, default=32)
    ap.add_argument("--box", type=float, default=4.0)
    ap.add_argument("--L", type=float, default=8e27, help="source luminosity record field (photons/s)")
    ap.add_argument("--sources", default="center", help="center (default), 8src, or a sources.txt")
    a = ap.parse_args()
    src = read_sources(a.sources, a.box)
    out = a.outdir
    os.makedirs(out, exist_ok=True)
    subprocess.check_call([sys.executable, os.path.join(ROOT, "tools/ics/make_ics.py"),
                           os.path.join(out, "ics"), "--kind", "zeldovich", "--gas", "uniform",
                           "--growth", "frozen-gas", "--Ngas", str(a.N), "--Ndm", str(a.N // 2),
                           "--z-init", "35", "--box", str(a.box), "--seed", "7"])
    # the dark matter at rest on its lattice (make_ics.py's order, x fastest):
    # make_ics.py's Zel'dovich displacements would grow structure in the gas
    # under gravity, which the refined runs (PM on level 0 only until Stage
    # G) and the uniform 2N run would follow differently
    ndm = a.N // 2
    q = np.arange(ndm, dtype=np.float32) + np.float32(0.5)
    xv = np.zeros((6, ndm, ndm, ndm), dtype="<f4")
    for c, Q in enumerate(np.meshgrid(q, q, q, indexing="ij")):
        xv[c] = Q
    xv.reshape(6, -1, order="F").ravel(order="F").tofile(os.path.join(out, "ics", "xv_dm.init"))
    with open(os.path.join(out, "params.txt"), "w") as f:
        f.write(f"""# {len(src)} source(s) in uniform gas (make_ifront.py --sources {a.sources}): the RT with refinement
N = {a.N}
Ndm = {a.N // 2}
box = {a.box}
nray_max_factor = 16
z_init = 35
zdata1 = 30 28 26 24 22 20 18 16 14 12
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
    # source_type: x(3),v(3) real(4), N int(4), M,R,T,L,sfr,mstar,fgas real(4)
    recs = b"".join(struct.pack("<6fi7f", x, y, z, 0, 0, 0, 40, 1e10, 0.02, 3e4, a.L * w, 0.5, 0, 0)
                    for x, y, z, w in src)
    for zc in (40.0, 10.0):
        with open(os.path.join(sh, "halo." + tag(zc)), "wb") as f:
            f.write(struct.pack("<q", len(src)) + recs)
    with open(os.path.join(out, "sources_pos.txt"), "w") as f:
        f.write("# x y z (box units) L (record field)\n")
        for x, y, z, w in src:
            f.write(f"{x / a.box:.8f} {y / a.box:.8f} {z / a.box:.8f} {a.L * w:.6e}\n")
    with open(os.path.join(sh, "sfr.txt"), "w") as f:
        for zc, sfr in ((40.0, 0.01), (10.0, 0.01)):
            f.write(f"{zc:10.4f}" + " 0" * 8 + f" {sfr:.6e}\n")
    print("wrote", out)


if __name__ == "__main__":
    main()
