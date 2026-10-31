#!/usr/bin/env python3
"""T1 of the zoom ladder (docs/zoom_turbulence.md section 5; docs/stage_g.md):
a coarse base with two refined levels over the central eighth, gravity
with the level solve, against the uniform run at the finest resolution.

usage: pm_zoom_t1_check.py WORKDIR LAUNCH "LAUNCH_FLAGS" NPFLAG NP RADHYDRO FIND_HALOS SRC
                           [--base 32] [--fine 128] [--box 10] [--zend 5] [--max-steps 400] [--subcycle]

The same Zel'dovich particles (fine^3, seed 7, BOX Mpc/h, z = 100 to ZEND)
and uniform gas in three runs (no RT):
  A  gas fine^3, no refinement (the reference: the fine^3 mesh everywhere);
  B  gas base^3, levels 1 and 2 over the central eighth (level 2 at the
     fine resolution there), pm_level_solve = 1 with pm_refine = 1: the
     base mesh outside the region by the FFT, the levels by MLMG;
  C  the same hierarchy with the uniform fine^3 PM mesh (pm_gas_levels = 1,
     pm_refine = fine/base): the Stage F coupling, the comparison;
  D  (--subcycle) B with pm_level_subcycle = 1: the refined levels' gravity
     in their own substeps (docs/stage_g.md section 6).
Inside the region (particles by index, the same ICs): the rms and maximum
displacement of B's and C's particles from A's (fine cells; also by
distance from the region's centre), the CIC density on the fine grid of
the region (rms difference, the fraction of cells above overdensities 1,
10 and 100), and the FoF halos (find_halos, b = 0.2, >= NMIN particles)
whose centres lie in the region: A's halos with >= 4 NMIN particles
matched to the nearest in B (and C) within 2 fine cells, the matched
fraction, the mass ratios and the rms log-mass difference.

What to expect: C has the fine mesh everywhere, so its particles are A's
to 1e-2 cells; B's base is 4x coarser outside the region, so the matter
that flows into the region carries the base's errors (0.5-0.9 fine cells
rms outside, 0.4 at the region's edge, 0.14 in its inner half at z = 5
in the 10 Mpc/h test), which is the zoom's resolution, not the solver's
(the solver's own test is the refined-everywhere one of
pm_level_solve_check.py). Passes when B's particles in the inner half of
the region (Chebyshev distance below 3/4 of the region's half-width)
stay within 0.3 fine cells rms of A's, 90 % of A's halos with >= 4 NMIN
particles in the region have a B halo within 2 fine cells, 75 % of those
within 40 % in mass, and B's density PDF fractions agree with A's to 25 %
where A has >= 100 cells above the threshold.
"""
import argparse
import glob
import os
import re
import shutil
import subprocess
import sys

import numpy as np

BASE = ["--quiet", "--set", "legacy_rt=0", "--set", "legacy_physics=0", "--set", "rt_backend=amrex",
        "--set", "hydro_scheme=unsplit", "--set", "hydro_backend=amrex", "--set", "amrex_hydro_max_grid_size=0",
        "--set", "pm_backend=amrex", "--set", "amrex_state=1", "--set", "amrex_max_grid_size=8",
        "--set", "deterministic_parallel=0", "--set", "tvir_min=1e30", "--set", "movie_dz=0"]
NMIN = 32


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


def last_dm(run_dir, npart):
    f = sorted(glob.glob(os.path.join(run_dir, "DM", "xv_dm*")))[-1]
    return f, np.fromfile(f, "<f4")[: 6 * npart].reshape(-1, 6).astype(float)


def read_catalog(path):
    raw = open(path, "rb").read()
    n = int(np.frombuffer(raw[:8], "<i8")[0])
    rec = np.dtype([("x", "<f4", 3), ("v", "<f4", 3), ("N", "<i4"), ("M", "<f4"), ("R", "<f4"), ("T", "<f4"),
                    ("L", "<f4"), ("sfr", "<f4"), ("mstar", "<f4"), ("fgas", "<f4")])
    return np.frombuffer(raw[8:8 + 56 * n], rec)


def cic(x, n, lo, hi):
    """CIC density of the particles with x in [lo, hi) (fine-grid units) on the (hi-lo)^3 grid, mean 1"""
    inside = np.all((x >= lo) & (x < hi), axis=1)
    p = x[inside] - lo
    m = hi - lo
    rho = np.zeros((m + 1, m + 1, m + 1))
    q = p - 0.5
    i0 = np.floor(q).astype(int)
    d = q - i0
    for dx in (0, 1):
        for dy in (0, 1):
            for dz in (0, 1):
                w = (d[:, 0] if dx else 1 - d[:, 0]) * (d[:, 1] if dy else 1 - d[:, 1]) * (d[:, 2] if dz else 1 - d[:, 2])
                np.add.at(rho, (np.clip(i0[:, 0] + dx, 0, m), np.clip(i0[:, 1] + dy, 0, m), np.clip(i0[:, 2] + dz, 0, m)), w)
    rho = rho[:m, :m, :m]
    return rho / (n ** 3 / float(n ** 3))  # particles per cell: one per fine cell on average


def main():
    # the eight positionals first (the launcher flags start with "--", which
    # argparse would take for options), then the options
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    pos = ("work", "launch", "flags", "npflag", "np", "exe", "find_halos", "src")
    if len(sys.argv) < 9:
        ap.error("usage: pm_zoom_t1_check.py " + " ".join(pos.upper() for pos in pos) + " [options]")
    ap.add_argument("--base", type=int, default=32)
    ap.add_argument("--fine", type=int, default=128)
    ap.add_argument("--box", type=float, default=10.0)
    ap.add_argument("--zend", type=float, default=5.0)
    ap.add_argument("--max-steps", type=int, default=400)
    ap.add_argument("--keep", action="store_true", help="reuse existing runs")
    ap.add_argument("--subcycle", action="store_true",
                    help="also D: B with pm_level_subcycle = 1 (docs/stage_g.md section 6), held to B's criteria")
    a = ap.parse_args(sys.argv[9:])
    for name, val in zip(pos, sys.argv[1:9]):
        setattr(a, name, val)
    flags = [f for f in a.flags.replace(";", " ").split() if f]
    if a.launch.endswith("srun") and not any("cpus-per-task" in f for f in flags):
        flags.append("--cpus-per-task=" + os.environ.get("OMP_NUM_THREADS", "1"))
    mpi = [a.launch] + flags + [a.npflag, a.np]
    work = a.work
    if os.path.isdir(work) and not a.keep:
        shutil.rmtree(work)
    os.makedirs(work, exist_ok=True)
    nb, nf = a.base, a.fine
    if nf % nb or (nf // nb) & (nf // nb - 1):
        sys.exit("fine/base must be a power of 2")
    nlev = (nf // nb).bit_length() - 1
    fails = []

    smoke = open(os.path.join(a.src, "tests", "problems_full", "smoke", "params.txt")).read()
    for n in (nb, nf):
        p = os.path.join(work, f"p{n}")
        if a.keep and os.path.isfile(os.path.join(p, "params.txt")):
            continue
        run([sys.executable, os.path.join(a.src, "tools", "ics", "make_ics.py"), os.path.join(p, "ics"),
             "--kind", "zeldovich", "--Ngas", str(n), "--Ndm", str(nf), "--z-init", "100", "--box", str(a.box),
             "--seed", "7", "--gas", "uniform"], os.path.join(work, f"ics{n}.log"))
        txt = re.sub(r"(?m)^N = .*$", f"N = {n}", smoke)
        txt = re.sub(r"(?m)^Ndm = .*$", f"Ndm = {nf}", txt)
        txt = re.sub(r"(?m)^box = .*$", f"box = {a.box:g}", txt)
        txt = re.sub(r"(?m)^zdata1 = .*$", f"zdata1 = {a.zend:g}", txt)
        open(os.path.join(p, "params.txt"), "w").write(txt)
    q = nb // 4
    region = f"{q},{q},{q},{3 * q - 1},{3 * q - 1},{3 * q - 1}"
    levels = sets(f"amrex_max_level={nlev}", "amrex_refine_mass=1e9", f"amrex_refine_box={region}",
                  "amrex_n_error_buf=0", "amrex_fine_max_grid_size=32")
    runs = {"A": (f"p{nf}", sets("pm_refine=1", "amrex_max_level=0")),
            "B": (f"p{nb}", sets("pm_refine=1", "pm_level_solve=1") + levels),
            "C": (f"p{nb}", sets(f"pm_refine={nf // nb}", "pm_gas_levels=1", "pm_level_solve=0") + levels)}
    if a.subcycle:
        runs["D"] = (f"p{nb}", sets("pm_refine=1", "pm_level_solve=1", "pm_level_subcycle=1") + levels)
    cmpn = [n for n in ("B", "C", "D") if n in runs]
    for name, (p, extra) in runs.items():
        o = os.path.join(work, name)
        if a.keep and glob.glob(os.path.join(o, "DM", "xv_dm*")):
            continue
        run(mpi + [a.exe, os.path.join(work, p), o] + BASE + sets(f"run_steps={a.max_steps}", "step_timers=1") + extra,
            os.path.join(work, name + ".log"))
    npart = nf ** 3
    lo, hi = q * (nf // nb), 3 * q * (nf // nb)  # the region in fine cells
    xv = {}
    files = {}
    for name in runs:
        files[name], xv[name] = last_dm(os.path.join(work, name), npart)
    xa = xv["A"][:, :3]
    inside = np.all((xa >= lo) & (xa < hi), axis=1)
    half = (hi - lo) / 2
    rc = np.abs(xa - (lo + half)).max(axis=1)  # Chebyshev distance from the region's centre
    inner = rc < 0.75 * half
    print(f"region: fine cells [{lo}, {hi})^3, {inside.sum()} of {npart} particles of A inside, {inner.sum()} in its inner half")
    stats = {}
    for name in cmpn:
        dx = xv[name][:, :3] - xa
        dx -= nf * np.round(dx / nf)
        d = np.sqrt((dx ** 2).sum(axis=1))
        dv = np.sqrt(((xv[name][:, 3:] - xv["A"][:, 3:]) ** 2).sum(axis=1))
        vr = np.sqrt((xv["A"][:, 3:] ** 2).sum(axis=1).mean())
        stats[name] = (float(np.sqrt((d[inside] ** 2).mean())), float(d[inside].max()),
                       float(np.sqrt((d[~inside] ** 2).mean())), float(np.sqrt((dv[inside] ** 2).mean()) / vr),
                       float(np.sqrt((d[inner] ** 2).mean())))
        print(f"{name} vs A: particles inside the region: rms {stats[name][0]:.3e} (inner half {stats[name][4]:.3e}), "
              f"max {stats[name][1]:.3e} fine cells; outside: rms {stats[name][2]:.3e}; velocities inside: rms "
              f"{stats[name][3]:.3e} of the rms velocity")
        shells = [(0, 0.25), (0.25, 0.5), (0.5, 0.75), (0.75, 1), (1, 1.25), (1.25, 2)]
        print("   rms displacement by Chebyshev distance from the region's centre (in half-widths): " +
              ", ".join(f"{a:g}-{b:g}: {np.sqrt((d[(rc >= a * half) & (rc < b * half)] ** 2).mean()):.3f}"
                        for a, b in shells))
    dbc = xv["B"][:, :3] - xv["C"][:, :3]
    dbc -= nf * np.round(dbc / nf)
    print(f"B vs C: particles inside: rms {np.sqrt((dbc[inside] ** 2).sum(axis=1).mean()):.3e} fine cells")
    # the density PDF of the region
    pdf = {}
    for name in runs:
        rho = cic(xv[name][:, :3], nf, lo, hi)
        fr = [float((rho > t).mean()) for t in (2, 11, 101)]
        pdf[name] = (rho, fr)
        print(f"{name}: region density: cells above overdensity 1 / 10 / 100: {fr[0]:.4f} / {fr[1]:.4f} / {fr[2]:.5f}")
    for name in cmpn:
        d = pdf[name][0] - pdf["A"][0]
        print(f"{name} vs A: region density rms difference / rms {np.sqrt((d ** 2).mean()) / np.sqrt((pdf['A'][0] ** 2).mean()):.3e}")
    # the halos
    z = float(re.search(r"z=(\d+\.\d+)", os.path.basename(files["A"])).group(1))
    cats = {}
    for name in runs:
        out = os.path.join(work, f"halos_{name}.bin")
        run([a.find_halos, files[name], str(z), out, "--b", "0.2", "--nmin", str(NMIN), "--box", str(a.box),
             "--ngas", str(nf)], os.path.join(work, f"halos_{name}.log"))
        c = read_catalog(out)
        x = c["x"] / a.box * nf  # fine cells
        keep = np.all((x >= lo) & (x < hi), axis=1)
        cats[name] = (x[keep], c["N"][keep].astype(float))
        print(f"{name}: {len(c)} halos with >= {NMIN} particles, {keep.sum()} inside the region, "
              f"{(c['N'][keep] >= 4 * NMIN).sum()} with >= {4 * NMIN}")
    xa_h, na_h = cats["A"]
    big = na_h >= 4 * NMIN
    match = {}
    for name in cmpn:
        xb, nbh = cats[name]
        npos, nmass, dm2, ratios = 0, 0, [], []
        for i in np.nonzero(big)[0]:
            d = xb - xa_h[i]
            d -= nf * np.round(d / nf)
            dd = np.sqrt((d ** 2).sum(axis=1))
            j = int(np.argmin(dd)) if len(dd) else -1
            if j >= 0 and dd[j] < 2.0:
                npos += 1
                ratios.append(nbh[j] / na_h[i])
                dm2.append(np.log10(nbh[j] / na_h[i]) ** 2)
                if abs(np.log10(nbh[j] / na_h[i])) < np.log10(1.4):
                    nmass += 1
        nb_ = max(1, int(big.sum()))
        match[name] = (npos / nb_, nmass / nb_)
        print(f"{name} vs A: halos with >= {4 * NMIN} particles inside the region: {npos} of {big.sum()} have a {name} halo "
              f"within 2 fine cells ({npos / nb_:.3f}), {nmass} of those within 40 % in mass ({nmass / nb_:.3f}); "
              f"rms dlog10 M of the matched {np.sqrt(np.mean(dm2)) if dm2 else float('nan'):.3e}; mass ratios "
              + " ".join(f"{x:.2f}" for x in ratios))
    # the criteria (B, and D with --subcycle)
    ncell = float((hi - lo) ** 3)
    for nm in [n for n in ("B", "D") if n in runs]:
        if not stats[nm][4] < 0.3:
            fails.append(f"{nm} particles in the inner half of the region: rms {stats[nm][4]:.3e} >= 0.3 fine cells")
        if big.sum() >= 5:
            if not match[nm][0] >= 0.9:
                fails.append(f"{nm} halos found within 2 cells {match[nm][0]:.3f} < 0.9")
            if not match[nm][1] >= 0.75:
                fails.append(f"{nm} halos within 40 % in mass {match[nm][1]:.3f} < 0.75")
        for t, (fa, fb) in enumerate(zip(pdf["A"][1], pdf[nm][1])):
            if fa * ncell >= 100 and not abs(fb / fa - 1) < 0.25:
                fails.append(f"{nm} density PDF fraction {t} {fb:.4g} vs A {fa:.4g}")
    if fails:
        print("FAIL:\n  " + "\n  ".join(fails))
        sys.exit(1)
    print("PASS")


if __name__ == "__main__":
    main()
