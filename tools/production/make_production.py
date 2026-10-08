#!/usr/bin/env python3
"""Set up the production reionization run (docs/production_run.md).

    make_production.py OUTDIR --N 128 [--box 10] [--z-init 100] [--seed 7]
                       [--c-reduced 0.01] [--nray-max-factor 16] [--threads 4]

Writes OUTDIR/params.txt and the Zel'dovich ICs (tools/ics/make_ics.py,
Ngas = Ndm = N) in OUTDIR/ics. Run with

    OMP_NUM_THREADS=4 build/radhydro OUTDIR OUTDIR/run --quiet

and restart after an interruption with

    build/radhydro OUTDIR OUTDIR/run --quiet --restart OUTDIR/run/restart/checkpoint.bin

The configuration: all fixes on (legacy_rt = legacy_physics = 0), the
deterministic parallel ray tracing, reduced speed of light, exact total RT
substeps (rt_substep_carry), on-the-fly FoF sources, RT from z = 15 (the
first zdata1 output), outputs every unit in z to z = 5.5, movie frames of
x_HI and T every movie_dz, a checkpoint after every hydro step.
"""
import argparse
import os
import subprocess
import sys

here = os.path.dirname(os.path.abspath(__file__))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir")
    ap.add_argument("--N", type=int, required=True)
    ap.add_argument("--box", type=float, default=10.0)
    ap.add_argument("--z-init", type=float, default=100.0)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--c-reduced", type=float, default=0.01)
    ap.add_argument("--nray-max-factor", type=int, default=16)
    ap.add_argument("--movie-dz", type=float, default=0.05)
    ap.add_argument("--ra-max", type=float, default=0.0,
                    help="max da/a per hydro step as ra_max/(1-ra_max); 0: original 0.1")
    ap.add_argument("--zdata1", default="15 14 13 12 11 10 9 8 7.5 7 6.5 6 5.5")
    ap.add_argument("--no-ics", action="store_true", help="params.txt only")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    ics = os.path.join(a.outdir, "ics")
    if not a.no_ics and not os.path.exists(os.path.join(ics, "u_gas.init")):
        subprocess.run([sys.executable, os.path.join(here, "..", "ics", "make_ics.py"), ics,
                        "--kind", "zeldovich", "--Ngas", str(a.N), "--Ndm", str(a.N),
                        "--z-init", str(a.z_init), "--box", str(a.box), "--seed", str(a.seed)],
                       check=True)
    lines = [
        "# Production reionization run (tools/production/make_production.py)",
        f"N = {a.N}",
        f"Ndm = {a.N}",
        f"box = {a.box:g}",
        f"nray_max_factor = {a.nray_max_factor}",
        f"z_init = {a.z_init:g}",
        f"zdata1 = {a.zdata1}",
        "ics = ics",
        "source_mode = onthefly",
        "fof_min_npart = 20",
        "legacy_rt = 0",
        "legacy_physics = 0",
        "deterministic_parallel = 1",
        f"c_reduced = {a.c_reduced:g}",
        "rt_substep_carry = 1",
        f"movie_dz = {a.movie_dz:g}",
        "movie_fields = fHI,temp",
        "restart_every = 1",
        "fields_at_data1 = 1",
        "fields = fHI,fHeI,fHeII,temp,rhog,gamma",
        "data_files = norays",
    ]
    if a.ra_max > 0:
        lines.append(f"ra_max = {a.ra_max:g}")
    with open(os.path.join(a.outdir, "params.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print("wrote", os.path.join(a.outdir, "params.txt"))


if __name__ == "__main__":
    main()
