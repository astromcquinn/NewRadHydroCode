#!/usr/bin/env python3
"""Set up the production reionization run (docs/production_run.md).

    make_production.py OUTDIR --N 128 [--box 10] [--z-init 100] [--seed 7]
                       [--c-schedule 0.1,0.25,0.2,0.5,0.3,0.75,0.5 | --c-schedule none --c-reduced 0.2]
                       [--nray-max-factor 16] [--threads 4]
                       [--hydro-scheme unsplit] [--redshift-per-distance 1]
                       [--omegam X --omegab X --hubble0 X --sigma8 X --ns X]
                       [--T-min K --z-floor Z]

Writes OUTDIR/params.txt and the Zel'dovich ICs (tools/ics/make_ics.py,
Ngas = N, Ndm = pm_refine*N or --ndm) in OUTDIR/ics. Run with

    OMP_NUM_THREADS=4 build/radhydro OUTDIR OUTDIR/run --quiet

and restart after an interruption with

    build/radhydro OUTDIR OUTDIR/run --quiet --restart OUTDIR/run/restart/checkpoint.bin

The configuration: all fixes on (legacy_rt = legacy_physics = 0, which includes
B17, photon_conserve, and P11, photon_consume_dense: dense gas consumes photons), the
deterministic parallel ray tracing, reduced speed of light, exact total RT
substeps (rt_substep_carry), the source-ray frame rotated every RT substep
and source (ray_rotate = 3, the B16 fix; --ray-rotate), the unsplit hydro with dual energy
(hydro_scheme = unsplit, docs/hydro_unsplit.md; --hydro-scheme) with dual_e_self = 1 (U7; --dual-e-self), radiation redshifted per distance
travelled (redshift_per_distance = 1; --redshift-per-distance), on-the-fly FoF sources, RT from z = 15 (the
first zdata1 output), outputs every unit in z to z = 5 (Lyman-alpha forest comparison), movie frames of
x_HI, T and the gas density every movie_dz (tools/movies/make_movie.py,
docs/movies.md), a checkpoint every --restart-every hydro steps (default 10: a 512^3 refined
checkpoint is ~240 GB and ~9 min to write; user decision 2026-10-18).

--backend amrex writes the AMReX backends with the state on the levels and
load balancing into params.txt (rt/hydro/pm_backend = amrex, amrex_state =
1, amrex_balance = 1, 32^3 boxes, deterministic_parallel = 0;
docs/stage_e_state.md): the MPI configuration for runs over several nodes,
which tools/hyak/reion_production.slurm launches with srun when it sees
rt_backend = amrex in params.txt (the restarts then use the same
parameters). It needs a build with AMReX (tools/hyak/build_amrex.sh).

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
    # the runtime cosmology (docs/production_run.md section 10): passed to
    # make_ics.py and written to params.txt; default the original's
    for k in ("omegam", "omegab", "hubble0", "sigma8", "ns"):
        ap.add_argument("--" + k, type=float, default=None,
                        help=f"runtime cosmology {k} (default: the original's)")
    ap.add_argument("--T-min", type=float, default=0.0,
                    help="the preheating floor in K at z < --z-floor, a stand-in for X-ray pre-heating "
                         "(T_min, docs/production_run.md section 10; C26: 100 or 1000 K; default 0 = off)")
    ap.add_argument("--z-floor", type=float, default=15.0,
                    help="the redshift the preheating floor switches on at (z_floor; default 15, "
                         "user decision 2026-10-29)")
    ap.add_argument("--c-reduced", type=float, default=None,
                    help="constant reduced speed of light; given alone it turns the default "
                         "schedule off (written as c_reduced; default 0.01 with --c-schedule none)")
    ap.add_argument("--c-schedule", default=None,
                    help="c_reduced_schedule c0,x1,c1,x2,c2,...: c_reduced by the box's x_HII,V "
                         "(production default, user decision 2026-10-18: 0.1c below x = 0.25, "
                         "0.2c to x = 0.5, 0.3c to 0.75, 0.5c after; docs/production_run.md "
                         "section 1); "
                         "'none' for the constant --c-reduced")
    ap.add_argument("--c-schedule-interp", type=int, default=0, choices=(0, 1),
                    help="1: c_reduced_schedule_interp = 1, the schedule continuous in x_HII,V "
                         "(linear between the plateau midpoints; docs/production_run.md section 1)")
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
                         "HI photoionization rate (1 float per cell, gamma has 3); gammaHI_avg its "
                         "average over the hydro step's RT substeps (serial RT only)")
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
    ap.add_argument("--restart-every", type=int, default=10)
    ap.add_argument("--recomb-case", default="B", choices=("A", "B"),
                    help="recombination case (docs/recomb_radiation.md). Default B: on the spot, "
                         "the original (user decision 2026-10-12, so that the production "
                         "configuration runs on rt_backend = amrex, which has no recomb_rad yet); "
                         "A with --recomb-rad 1 transports the recombination radiation (the "
                         "default from 2026-10-11 to 2026-10-12, docs/rt_methods_512.md)")
    ap.add_argument("--recomb-rad", type=int, default=0, choices=(0, 1),
                    help="transport the ground-state recombination photons (case A only)")
    ap.add_argument("--recomb-rad-block", type=int, default=4,
                    help="rays: recombination emission gathered over b^3 RT cells (4: 79%% "
                         "launched at 78 kpc/h cells, docs/recomb_radiation.md section 4.1)")
    ap.add_argument("--recomb-rad-cell-onspot", type=float, default=1.0,
                    help="rays: cells with optical depth > tau at the HI edge keep their "
                         "recombination photons on the spot (recomb_rad_cell_onspot; default 1, "
                         "docs/rt_methods_512.md); 0: off")
    ap.add_argument("--moment-inject", type=int, default=1, choices=(0, 1),
                    help="M1/OTVET (rt_method set by tools/hyak/rt_methods_runs.sh; ignored by "
                         "the rays): 1 = source and recombination photons enter the moments "
                         "after the transport, so they are absorbed in their own cell first "
                         "(moment_inject; default 1, docs/rt_moments.md); 0: before, the original "
                         "order, which leaves dense source cells neutral")
    ap.add_argument("--data-files", default="nodm", choices=("all", "norays", "nodm", "dm", "gas"),
                    help="files of the zdata1 snapshots besides the fields: nodm (default, user "
                         "decision 2026-10-18: Gas and RT/rt_grid, no particles; they stay in "
                         "the checkpoints), norays (DM, Gas, RT/rt_grid), gas (only Gas/u_gas, enough for igm_stats.py without "
                         "the source-start mean free path), dm (only the particles: 26 GB "
                         "per output at 1024^3, docs/pm_refine.md)")
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
    ap.add_argument("--thermal-floor", type=float, default=None,
                    help="thermal_floor in K (docs/thermal_history.md section 4): the hydro floor and, "
                         "with legacy_physics = 0, the clamp of the RT physics floor below 1 K; "
                         "unset: the code default (1 K), or 0.1 K with --plane-parallel (user "
                         "decision 2026-10-29)")
    ap.add_argument("--dual-e-self", type=int, default=1, choices=(0, 1),
                    help="dual_e_self with the unsplit hydro (docs/hydro_unsplit.md section 7, "
                         "U7): 1 for production (rho e evolved as itself where the pressure comes "
                         "from E - KE, no over-cooled or numerically heated IGM; user decision "
                         "2026-10-27); 0: the scheme before the key")
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
    ap.add_argument("--photon-conserve", type=int, default=-1, choices=(-1, 0, 1),
                    help="iterate photo_cell's photon limit to convergence (photon_conserve, "
                         "the B17 fix, docs/rt_multicell.md); -1 (default, no line written): "
                         "follow legacy_physics, i.e. on in production (legacy_physics = 0)")
    ap.add_argument("--photon-consume-dense", type=int, default=-1, choices=(-1, 0, 1),
                    help="every photoionisation consumes its photon, also in cells denser than "
                         "rho_source = 50 times the mean (photon_consume_dense, P11, docs/port_notes.md "
                         "5b); -1 (default, no line written): follow legacy_physics, i.e. on in "
                         "production (legacy_physics = 0, user decision 2026-10-28); 0: the original "
                         "rho_source factor")
    ap.add_argument("--cool-suppress-overdensity", type=float, default=200.0,
                    help="overdensity above which cooling and SN heating are suppressed "
                         "(cool_suppress_overdensity, docs/stage_f.md); 200: the original (default, "
                         "no line written); 0: off (dense gas cools)")
    ap.add_argument("--sn-heating", type=int, default=0, choices=(0, 1),
                    help="photo_cell's SN heating above rho_source = 50 times the mean in "
                         "photoheated gas (sn_heating); 0 (default, user decision 2026-10-30: "
                         "no SN heating in any run): off; 1: the original")
    ap.add_argument("--amrex-max-level", type=int, default=0,
                    help="refinement levels of the gas hierarchy (amrex_max_level, docs/stage_f.md; "
                         "--backend amrex only); 0: none (default, no line written)")
    ap.add_argument("--amrex-refine-mass", type=float, default=8.0,
                    help="a level-l cell is refined above amrex_refine_mass * 8^l times the mean gas "
                         "density (default 8, no line written)")
    ap.add_argument("--sf-overdensity", type=float, default=0.0,
                    help="the star-formation sink: gas above this overdensity forms star particles "
                         "(sf_overdensity, docs/stage_f.md section 2; --backend amrex only); 0: off "
                         "(default, no line written); the user's production choice is a few thousand")
    ap.add_argument("--sf-efficiency", type=float, default=1.0,
                    help="fraction of a sink cell's mass turned into stars per free-fall time "
                         "(sf_efficiency, default 1, no line written)")
    ap.add_argument("--pm-refine", type=int, default=1,
                    help="PM mesh refinement (docs/pm_refine.md): writes pm_refine and "
                         "Ndm = pm_refine*N, and makes ICs with pm_refine*N particles per "
                         "dimension whose large scales are those of the N^3 ICs of the same "
                         "seed (--phases-from-ndm N); 1: the original, Ndm = N. The abundance matching (any box, "
                         "any pm_refine): --am-params data/uv_lf/am_params_st_fof.txt "
                         "(docs/pm_refine.md section 6)")
    ap.add_argument("--ndm", type=int, default=0,
                    help="dark-matter particles per dimension (default 0: pm_refine*N); must divide "
                         "the mesh pm_refine*N, e.g. --pm-refine 4 --ndm 1024 at N = 512 (the 2048^3 "
                         "mesh with 1024^3 particles; ICs with --phases-from-ndm N when ndm > N)")
    ap.add_argument("--pm-gas-levels", type=int, default=1, choices=(0, 1),
                    help="with --amrex-max-level > 0: the refined gas levels enter the PM deposit and "
                         "kick (pm_gas_levels, docs/stage_f.md section 4; default 1, no line written); "
                         "0: the level-0 coupling")
    ap.add_argument("--backend", default="serial", choices=("serial", "amrex"),
                    help="amrex: the AMReX RT, hydro and PM backends with the state on the levels "
                         "and load balancing (rt/hydro/pm_backend = amrex, amrex_state = 1, "
                         "amrex_balance = 1, 32^3 boxes; docs/stage_e_state.md), for MPI runs "
                         "over several nodes with tools/hyak/reion_production.slurm; needs "
                         "--hydro-scheme unsplit and no --recomb-rad. Default: the serial code")
    ap.add_argument("--fof-dt-myr", type=float, default=0.0,
                    help="run the on-the-fly FoF only every this many Myr, reusing its halos "
                         "in between (docs/pm_refine.md, FoF cadence); 0: every hydro step "
                         "(the original, no line written)")
    ap.add_argument("--plane-parallel", action="store_true",
                    help="the plane-parallel RT-domain mode (docs/plane_parallel.md): the faces of "
                         "cubic RT domains as the only sources from --plane-z-re (no halo finder, "
                         "no abundance matching); implies cooling on everywhere unless "
                         "--cool-suppress-overdensity is given")
    ap.add_argument("--plane-transport", choices=("rays", "sweep"), default="sweep",
                    help="the packets in flight as rays of the ray pool, or as a field per lit face "
                         "shifted along its axis (default sweep: no ray pool, a fifth of the memory; the "
                         "same photons to round-off)")
    ap.add_argument("--plane-ndom", type=int, default=8, help="RT domains per dimension (default 8)")
    ap.add_argument("--plane-faces", default="+x +y",
                    help='the emitting faces of every domain (default "+x +y": two orthogonal '
                         'directions, one face each, as DM20 and C26); "+x" with --plane-ndom 1 is one '
                         'front across the whole box')
    ap.add_argument("--plane-gamma12", type=float, default=0.3,
                    help="Gamma_HI in the source cells in 1e-12 s^-1 (default 0.3, C26)")
    ap.add_argument("--plane-z-re", type=float, default=7.0, help="the turn-on redshift (default 7)")
    ap.add_argument("--plane-alpha", type=float, default=1.5, help="J_nu ~ nu^-alpha (default 1.5)")
    ap.add_argument("--dt-max-myr", type=float, default=0.0,
                    help="cap on the hydro step in Myr from the turn-on (0: off; the T4 test)")
    ap.add_argument("--rt-substeps-max", type=int, default=0,
                    help="cap on the RT substeps per hydro step from the turn-on (0: off; the T4 test)")
    a = ap.parse_args()
    if a.plane_parallel:
        # DM20's adaptive speed of light unless given: 0.1 c until x_HI,V < 1 per cent, then 0.01 c
        if a.c_schedule is None and a.c_reduced is None:
            a.c_schedule = "0.1,0.99,0.01"
        if a.N % a.plane_ndom:
            sys.exit(f"--plane-ndom {a.plane_ndom} must divide N = {a.N}")
        if not any(x.startswith("--cool-suppress-overdensity") for x in sys.argv):
            a.cool_suppress_overdensity = 0.0  # cooling on everywhere (user decision 2026-10-26)
        if a.am_params:
            sys.exit("--plane-parallel uses no halo sources: drop --am-params")
        if a.thermal_floor is None:
            a.thermal_floor = 0.1  # the unheated gas ahead of the fronts is below 1 K at z < 8 (user decision 2026-10-29)
        a.zdata1 = a.zdata1  # the outputs; the RT starts at --plane-z-re, not at the first zdata1
    # the speed of light: the default schedule, unless a constant --c-reduced
    # is given alone (then that constant)
    if a.c_schedule is None:
        a.c_schedule = "none" if a.c_reduced is not None else "0.1,0.25,0.2,0.5,0.3,0.75,0.5"
    if a.c_reduced is None:
        a.c_reduced = 0.01 if a.c_schedule in ("", "none") else float(a.c_schedule.split(",")[0])
    ndm = a.ndm if a.ndm > 0 else a.pm_refine * a.N
    if (a.pm_refine * a.N) % ndm:
        sys.exit(f"--ndm {ndm} must divide the PM mesh pm_refine*N = {a.pm_refine * a.N}")
    if a.backend == "amrex":
        if a.hydro_scheme != "unsplit":
            sys.exit("--backend amrex needs --hydro-scheme unsplit")
        if a.recomb_rad:
            sys.exit("--backend amrex has no recombination radiation transport (--recomb-rad 0)")
        if a.N % 32:
            sys.exit("--backend amrex: N must be a multiple of the 32^3 boxes")
    if a.nray_max_factor * a.N ** 3 >= 2 ** 31:
        sys.exit(f"--nray-max-factor {a.nray_max_factor} at N = {a.N}: Nray_max >= 2^31 overflows the "
                 f"int32 ray links (rays.hpp); use at most {(2 ** 31 - 1) // a.N ** 3}")
    if float(a.zdata1.split()[0]) < 15 and not a.plane_parallel:
        print(f"warning: --zdata1 starts at {a.zdata1.split()[0]}: the RT starts at the first "
              f"zdata1 redshift, so sources before it do not shine; start the list at >= 15")
    os.makedirs(a.outdir, exist_ok=True)
    ics = os.path.join(a.outdir, "ics")
    if not a.no_ics and not os.path.exists(os.path.join(ics, "u_gas.init")):
        cmd = [sys.executable, os.path.join(here, "..", "ics", "make_ics.py"), ics,
               "--kind", "zeldovich", "--Ngas", str(a.N), "--Ndm", str(ndm),
               "--z-init", str(a.z_init), "--box", str(a.box), "--seed", str(a.seed)]
        if ndm > a.N:
            cmd += ["--phases-from-ndm", str(a.N)]
        for k in ("omegam", "omegab", "hubble0", "sigma8", "ns"):
            if getattr(a, k) is not None:
                cmd += ["--" + k, repr(getattr(a, k))]
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
        # the serial code's deterministic parallel ray tracing; the backend
        # is deterministic by construction and refuses the flag
        f"deterministic_parallel = {0 if a.backend == 'amrex' else 1}",
        f"c_reduced = {a.c_reduced:g}",
        *([] if a.c_schedule in ("", "none") else
          [f"c_reduced_schedule = {' '.join(x.strip() for x in a.c_schedule.split(','))}"]
          + ([] if not a.c_schedule_interp else ["c_reduced_schedule_interp = 1"])),
        "rt_substep_carry = 1",
        f"ray_rotate = {a.ray_rotate}",
    ]
    for k in ("omegam", "omegab", "hubble0", "sigma8", "ns"):
        if getattr(a, k) is not None:
            lines.append(f"{k} = {getattr(a, k)!r}")
    if a.T_min > 0:
        lines += ["# the preheating floor (docs/production_run.md section 10): a plain switch-on",
                  f"T_min = {a.T_min:g}", f"z_floor = {a.z_floor:g}"]
    if a.pm_refine > 1:
        lines.append(f"pm_refine = {a.pm_refine}")
    if a.fof_dt_myr > 0:
        lines.append(f"fof_dt_myr = {a.fof_dt_myr:g}")
    lines.append(f"hydro_scheme = {a.hydro_scheme}")
    if a.hydro_scheme == "unsplit":
        lines.append(f"dual_e_self = {a.dual_e_self}")
    if a.backend == "amrex":
        lines += [
            "# The AMReX backends, the state on the levels, load balancing "
            "(docs/stage_d.md, stage_e_hydro.md, stage_e_pm.md, stage_e_state.md);",
            "# reion_production.slurm launches MPI ranks when it sees rt_backend = amrex",
            "rt_backend = amrex",
            "hydro_backend = amrex",
            "pm_backend = amrex",
            "amrex_max_grid_size = 32",
            "amrex_hydro_max_grid_size = 0",
            "amrex_state = 1",
            "amrex_balance = 1",
        ]
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
        if a.recomb_rad_cell_onspot > 0:
            lines.append(f"recomb_rad_cell_onspot = {a.recomb_rad_cell_onspot:g}")
    elif a.recomb_rad and a.recomb_case == "B":
        sys.exit("--recomb-rad 1 needs --recomb-case A (case B already absorbs them on the spot)")
    if a.moment_inject:
        lines.append("moment_inject = 1")
    if a.ra_max > 0:
        lines.append(f"ra_max = {a.ra_max:g}")
    lines.append(f"redshift_per_distance = {a.redshift_per_distance}")
    if a.photon_conserve >= 0:
        lines.append(f"photon_conserve = {a.photon_conserve}")
    if a.photon_consume_dense >= 0:
        lines.append(f"photon_consume_dense = {a.photon_consume_dense}")
    if a.ray_deposit != "share":
        lines.append(f"ray_deposit = {a.ray_deposit}")
    if a.cool_suppress_overdensity != 200.0:
        lines.append(f"cool_suppress_overdensity = {a.cool_suppress_overdensity:g}")
    if a.sn_heating != 1:
        lines.append(f"sn_heating = {a.sn_heating}")
    if a.amrex_max_level > 0 or a.sf_overdensity > 0:
        if a.backend != "amrex":
            sys.exit("--amrex-max-level and --sf-overdensity need --backend amrex")
    if a.amrex_max_level > 0:
        lines.append(f"amrex_max_level = {a.amrex_max_level}")
        if a.amrex_refine_mass != 8.0:
            lines.append(f"amrex_refine_mass = {a.amrex_refine_mass:g}")
        if a.pm_gas_levels != 1:
            lines.append(f"pm_gas_levels = {a.pm_gas_levels}")
    if a.sf_overdensity > 0:
        lines.append(f"sf_overdensity = {a.sf_overdensity:g}")
        if a.sf_efficiency != 1.0:
            lines.append(f"sf_efficiency = {a.sf_efficiency:g}")
    if a.rt_cells_per_substep != 1:
        lines.append(f"rt_cells_per_substep = {a.rt_cells_per_substep}")
    if a.thermal_floor is not None:
        lines += ["# the temperature floor (docs/thermal_history.md section 4): the hydro floor and,",
                  "# with legacy_physics = 0, the clamp of the RT physics floor below 1 K",
                  f"thermal_floor = {a.thermal_floor:g}"]
    if a.plane_parallel:
        lines += ["# Plane-parallel RT-domain mode (docs/plane_parallel.md): the domain faces are the",
                  "# only sources from plane_z_re; the halo finder and the sources are not used",
                  "plane_parallel = 1",
                  f"plane_transport = {a.plane_transport}",
                  f"plane_ndom = {a.plane_ndom}",
                  f"plane_faces = {a.plane_faces}",
                  f"plane_gamma12 = {a.plane_gamma12:g}",
                  f"plane_z_re = {a.plane_z_re:g}",
                  "# the hydro steps end exactly on the zdata1 outputs (the same times for every step cap)",
                  "land_outputs = 1"]
        if a.plane_alpha != 1.5:
            lines.append(f"plane_alpha = {a.plane_alpha:g}")
    if a.plane_parallel or a.dt_max_myr > 0 or a.rt_substeps_max > 0:
        lines.append("max_steps = 100000")  # the step caps (and short post-ionisation CFL steps) need more than the original's 1000
    if a.dt_max_myr > 0:
        lines.append(f"dt_max_myr = {a.dt_max_myr:g}")
    if a.rt_substeps_max > 0:
        lines.append(f"rt_substeps_max = {a.rt_substeps_max}")
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
