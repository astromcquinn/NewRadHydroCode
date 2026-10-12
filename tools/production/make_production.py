#!/usr/bin/env python3
"""Set up the production reionization run (docs/production_run.md).

    make_production.py OUTDIR --N 128 [--box 10] [--z-init 100] [--seed 7]
                       [--c-reduced 0.01] [--nray-max-factor 16] [--threads 4]
                       [--hydro-scheme unsplit] [--redshift-per-distance 1]

Writes OUTDIR/params.txt and the Zel'dovich ICs (tools/ics/make_ics.py,
Ngas = Ndm = N) in OUTDIR/ics. Run with

    OMP_NUM_THREADS=4 build/radhydro OUTDIR OUTDIR/run --quiet

and restart after an interruption with

    build/radhydro OUTDIR OUTDIR/run --quiet --restart OUTDIR/run/restart/checkpoint.bin

The configuration: all fixes on (legacy_rt = legacy_physics = 0), the
deterministic parallel ray tracing, reduced speed of light, exact total RT
substeps (rt_substep_carry), the source-ray frame rotated every RT substep
and source (ray_rotate = 3, the B16 fix; --ray-rotate), the unsplit hydro with dual energy
(hydro_scheme = unsplit, docs/hydro_unsplit.md; --hydro-scheme), radiation redshifted per distance
travelled (redshift_per_distance = 1; --redshift-per-distance), on-the-fly FoF sources, RT from z = 15 (the
first zdata1 output), outputs every unit in z to z = 5 (Lyman-alpha forest comparison), movie frames of
x_HI, T and the gas density every movie_dz (tools/movies/make_movie.py,
docs/movies.md), a checkpoint every --restart-every hydro steps (default 1).

Movie disk: a full-cube frame takes 4N^3 B per field (512 MB at 512^3).
--movie-slab t writes only the three slabs of t cells through the box
centre instead (movie_slab = t: 3 x t x N^2 x 4 B per field, 3.1 MB at
512^3 with t = 1), which is all that make_movie.py, stills.py and reion.py need
(docs/movies.md, "Slab frames"). Default since 2026-10-12 (user): t = 1, the
three central planes (the movies show one plane); 0: full cubes.
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
    ap.add_argument("--movie-slab", type=int, default=1,
                    help="movie frames hold only slabs of this many cells normal to x, y and z "
                         "through the box centre (movie_slab; default 1: the three central "
                         "planes; more cells allow thin projections); 0: full cubes")
    ap.add_argument("--ra-max", type=float, default=0.0,
                    help="max da/a per hydro step as ra_max/(1-ra_max); 0: original 0.1")
    ap.add_argument("--zdata1", default="15 14 13 12 11 10 9 8 7.5 7 6.5 6 5.5 5",
                    help="snapshot redshifts; the RT (sources) starts at the first one, so start "
                         "at >= 15")
    ap.add_argument("--fields", default="fHI,fHeI,fHeII,temp,rhog,gamma",
                    help="step_NNNN field files at each zdata1 output; gammaHI writes only the "
                         "HI photoionization rate (1 float per cell, gamma has 3)")
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
    ap.add_argument("--recomb-case", default="A", choices=("A", "B"),
                    help="recombination case (docs/recomb_radiation.md); production default "
                         "since 2026-10-11: A with the recombination radiation transported "
                         "(--recomb-rad 1); B: on the spot, the original")
    ap.add_argument("--recomb-rad", type=int, default=1, choices=(0, 1),
                    help="transport the ground-state recombination photons (case A only)")
    ap.add_argument("--recomb-rad-block", type=int, default=4,
                    help="rays: recombination emission gathered over b^3 RT cells (4: 79%% "
                         "launched at 78 kpc/h cells, docs/recomb_radiation.md section 4.1)")
    ap.add_argument("--recomb-rad-cell-onspot", type=int, default=1, choices=(0, 1),
                    help="rays: optically thick cells keep their recombination photons on the "
                         "spot (recomb_rad_cell_onspot; default 1, docs/rt_methods_512.md)")
    ap.add_argument("--data-files", default="norays", choices=("all", "norays", "dm", "gas"),
                    help="files of the zdata1 snapshots besides the fields: norays (DM, Gas, "
                         "RT/rt_grid), gas (only Gas/u_gas, enough for igm_stats.py without "
                         "the source-start mean free path)")
    ap.add_argument("--restart-keep", type=int, default=0,
                    help="keep the newest n numbered checkpoints (restart_keep); 0: one "
                         "checkpoint.bin, rewritten (no line written)")
    ap.add_argument("--ray-rotate", type=int, default=3, choices=(0, 1, 2, 3),
                    help="HEALPix frame rotation of emitted rays (B16 fix, docs/port_notes.md); "
                         "3: every RT substep and source (production), 0: the original")
    ap.add_argument("--hydro-scheme", default="unsplit", choices=("legacy", "unsplit"),
                    help="hydro_scheme (docs/hydro_unsplit.md): unsplit for production (it fixes "
                         "the cold-gas overheating of the legacy scheme; user decision "
                         "2026-10-11); legacy: the original, the code default")
    ap.add_argument("--redshift-per-distance", type=int, default=1, choices=(0, 1),
                    help="1: redshift the radiation per distance travelled instead of per time "
                         "(fscale^s, docs/production_run.md; production, user decision "
                         "2026-10-11); 0: the original, per time")
    ap.add_argument("--rt-cells-per-substep", type=int, default=1,
                    help="RT cells a ray crosses per RT substep (rt_cells_per_substep, "
                         "docs/rt_multicell.md); 1: the original (default, no line written)")
    ap.add_argument("--ray-deposit", default="share", choices=("share", "exposure"),
                    help="how rays deposit photons along their path (ray_deposit, "
                         "docs/rt_multicell.md); share: the original (default, no line)")
    ap.add_argument("--photon-conserve", type=int, default=0, choices=(0, 1),
                    help="iterate photo_cell's photon limit to convergence (photon_conserve, "
                         "the B17 fix, docs/rt_multicell.md); 0: the original (default, no line)")
    a = ap.parse_args()
    if a.nray_max_factor * a.N ** 3 >= 2 ** 31:
        sys.exit(f"--nray-max-factor {a.nray_max_factor} at N = {a.N}: Nray_max >= 2^31 overflows the "
                 f"int32 ray links (rays.hpp); use at most {(2 ** 31 - 1) // a.N ** 3}")
    if float(a.zdata1.split()[0]) < 15:
        print(f"warning: --zdata1 starts at {a.zdata1.split()[0]}: the RT starts at the first "
              f"zdata1 redshift, so sources before it do not shine; start the list at >= 15")
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
    ]
    lines.append(f"hydro_scheme = {a.hydro_scheme}")
    lines += [
        f"movie_dz = {a.movie_dz:g}",
        f"movie_fields = {a.movie_fields}",
    ] + ([f"movie_slab = {a.movie_slab}"] if a.movie_slab > 0 else []) + [
        f"restart_every = {a.restart_every}",
    ] + ([f"restart_keep = {a.restart_keep}"] if a.restart_keep > 0 else []) + [
        "fields_at_data1 = 1",
        f"fields = {a.fields}",
        f"data_files = {a.data_files}",
        f"recomb_case = {a.recomb_case}",
    ]
    if a.recomb_case == "A" and a.recomb_rad:
        lines += ["recomb_rad = 1", f"recomb_rad_block = {a.recomb_rad_block}"]
        if a.recomb_rad_cell_onspot:
            lines.append("recomb_rad_cell_onspot = 1")
    elif a.recomb_rad and a.recomb_case == "B":
        sys.exit("--recomb-rad 1 needs --recomb-case A (case B already absorbs them on the spot)")
    if a.ra_max > 0:
        lines.append(f"ra_max = {a.ra_max:g}")
    lines.append(f"redshift_per_distance = {a.redshift_per_distance}")
    if a.photon_conserve:
        lines.append("photon_conserve = 1")
    if a.ray_deposit != "share":
        lines.append(f"ray_deposit = {a.ray_deposit}")
    if a.rt_cells_per_substep != 1:
        lines.append(f"rt_cells_per_substep = {a.rt_cells_per_substep}")
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
