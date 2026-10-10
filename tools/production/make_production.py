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
substeps (rt_substep_carry), the source-ray frame rotated every RT substep
and source (ray_rotate = 3, the B16 fix; --ray-rotate), on-the-fly FoF sources, RT from z = 15 (the
first zdata1 output), outputs every unit in z to z = 5 (Lyman-alpha forest comparison), movie frames of
x_HI, T and the gas density every movie_dz (tools/movies/make_movie.py,
docs/movies.md), a checkpoint every --restart-every hydro steps (default 1).
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
    ap.add_argument("--movie-fields", default="fHI,temp,rhog",
                    help="rt_grid fields per movie frame (tools/movies/make_movie.py uses "
                         "fHI, temp, rhog); 4N^3 B each per frame")
    ap.add_argument("--ra-max", type=float, default=0.0,
                    help="max da/a per hydro step as ra_max/(1-ra_max); 0: original 0.1")
    ap.add_argument("--zdata1", default="15 14 13 12 11 10 9 8 7.5 7 6.5 6 5.5 5")
    ap.add_argument("--no-ics", action="store_true", help="params.txt only")
    ap.add_argument("--am-params", default="",
                    help="am_params.txt from apps/am_fit: abundance-matched sources "
                         "(sfr_model = abundance, docs/abundance_matching.md)")
    ap.add_argument("--xi-ion", type=float, default=25.38,
                    help="log10 xi_ion with --am-params (production choice 25.38)")
    ap.add_argument("--f-esc", type=float, default=0.0,
                    help="escape fraction with --am-params (calibrated to tau_es); with "
                         "--f-esc-beta, its value at z = --f-esc-z0")
    ap.add_argument("--f-esc-beta", type=float, default=0.0,
                    help="f_esc(z) = f_esc ((1+z)/(1+f_esc_z0))^beta, capped at 1 "
                         "(tools/analysis/emissivity.py); 0: constant f_esc")
    ap.add_argument("--f-esc-z0", type=float, default=6.0)
    ap.add_argument("--restart-every", type=int, default=1)
    ap.add_argument("--ray-rotate", type=int, default=3, choices=(0, 1, 2, 3),
                    help="HEALPix frame rotation of emitted rays (B16 fix, docs/port_notes.md); "
                         "3: every RT substep and source (production), 0: the original")
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
        f"ray_rotate = {a.ray_rotate}",
        f"movie_dz = {a.movie_dz:g}",
        f"movie_fields = {a.movie_fields}",
        f"restart_every = {a.restart_every}",
        "fields_at_data1 = 1",
        "fields = fHI,fHeI,fHeII,temp,rhog,gamma",
        "data_files = norays",
    ]
    if a.ra_max > 0:
        lines.append(f"ra_max = {a.ra_max:g}")
    if a.am_params:
        if not a.f_esc > 0:
            sys.exit("--am-params needs --f-esc")
        lines.append("# Abundance-matched sources (" + a.am_params + ")")
        for l in open(a.am_params):
            l = l.split("#")[0].strip()
            if l:
                lines.append(l)
        lines.append(f"xi_ion = {a.xi_ion:g}")
        lines.append(f"f_esc = {a.f_esc:g}")
        if a.f_esc_beta != 0:
            lines.append(f"f_esc_beta = {a.f_esc_beta:g}")
            lines.append(f"f_esc_z0 = {a.f_esc_z0:g}")
    with open(os.path.join(a.outdir, "params.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print("wrote", os.path.join(a.outdir, "params.txt"))


if __name__ == "__main__":
    main()
