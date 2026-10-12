#!/usr/bin/env python3
"""Set up the production reionization run (docs/production_run.md).

    make_production.py OUTDIR --N 128 [--box 10] [--z-init 100] [--seed 7]
                       [--c-reduced 0.01] [--nray-max-factor 16] [--threads 4]
                       [--hydro-scheme unsplit]

Writes OUTDIR/params.txt and the Zel'dovich ICs (tools/ics/make_ics.py,
Ngas = N, Ndm = pm_refine*N) in OUTDIR/ics. Run with

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
                    help="escape fraction with --am-params (calibrated to tau_es)")
    ap.add_argument("--restart-every", type=int, default=1)
    ap.add_argument("--ray-rotate", type=int, default=3, choices=(0, 1, 2, 3),
                    help="HEALPix frame rotation of emitted rays (B16 fix, docs/port_notes.md); "
                         "3: every RT substep and source (production), 0: the original")
    ap.add_argument("--hydro-scheme", default="", choices=("", "legacy", "unsplit"),
                    help="hydro_scheme line (docs/hydro_unsplit.md): unsplit for production (it "
                         "fixes the cold-gas overheating of the legacy scheme); empty: not written, "
                         "i.e. the code default, legacy")
    ap.add_argument("--redshift-per-distance", type=int, default=0, choices=(0, 1),
                    help="redshift the radiation per distance travelled instead of per time "
                         "(fscale^s, docs/production_run.md); 0: the original (default, "
                         "no line written)")
    ap.add_argument("--pm-refine", type=int, default=1,
                    help="PM mesh refinement (docs/pm_refine.md): writes pm_refine and "
                         "Ndm = pm_refine*N, and makes ICs with pm_refine*N particles per "
                         "dimension whose large scales are those of the N^3 ICs of the same "
                         "seed (--phases-from-ndm N); 1: the original, Ndm = N. The abundance matching (any box, "
                         "any pm_refine): --am-params data/uv_lf/am_params_st_fof.txt "
                         "(docs/pm_refine.md section 6)")
    ap.add_argument("--fof-dt-myr", type=float, default=0.0,
                    help="run the on-the-fly FoF only every this many Myr, reusing its halos "
                         "in between (docs/pm_refine.md, FoF cadence); 0: every hydro step "
                         "(the original, no line written)")
    ap.add_argument("--data-files", default="norays", choices=("all", "norays", "dm"),
                    help="data_files line (dm: only the particles at the zdata1 outputs)")
    a = ap.parse_args()
    ndm = a.pm_refine * a.N
    if a.nray_max_factor * a.N ** 3 >= 2 ** 31:
        sys.exit(f"--nray-max-factor {a.nray_max_factor} at N = {a.N}: Nray_max >= 2^31 overflows the "
                 f"int32 ray links (rays.hpp); use at most {(2 ** 31 - 1) // a.N ** 3}")
    os.makedirs(a.outdir, exist_ok=True)
    ics = os.path.join(a.outdir, "ics")
    if not a.no_ics and not os.path.exists(os.path.join(ics, "u_gas.init")):
        cmd = [sys.executable, os.path.join(here, "..", "ics", "make_ics.py"), ics,
               "--kind", "zeldovich", "--Ngas", str(a.N), "--Ndm", str(ndm),
               "--z-init", str(a.z_init), "--box", str(a.box), "--seed", str(a.seed)]
        if a.pm_refine > 1:
            cmd += ["--phases-from-ndm", str(a.N)]
        subprocess.run(cmd, check=True)
    lines = [
        "# Production reionization run (tools/production/make_production.py)",
        f"N = {a.N}",
        f"Ndm = {ndm}",
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
    ]
    if a.pm_refine > 1:
        lines.append(f"pm_refine = {a.pm_refine}")
    if a.fof_dt_myr > 0:
        lines.append(f"fof_dt_myr = {a.fof_dt_myr:g}")
    if a.hydro_scheme:
        lines.append(f"hydro_scheme = {a.hydro_scheme}")
    lines += [
        f"movie_dz = {a.movie_dz:g}",
        f"movie_fields = {a.movie_fields}",
        f"restart_every = {a.restart_every}",
        "fields_at_data1 = 1",
        "fields = fHI,fHeI,fHeII,temp,rhog,gamma",
        f"data_files = {a.data_files}",
    ]
    if a.ra_max > 0:
        lines.append(f"ra_max = {a.ra_max:g}")
    if a.redshift_per_distance:
        lines.append("redshift_per_distance = 1")
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
    with open(os.path.join(a.outdir, "params.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print("wrote", os.path.join(a.outdir, "params.txt"))


if __name__ == "__main__":
    main()
