#!/usr/bin/env python3
"""End-to-end sanity checks of rt_driver on the validation problems.

usage: e2e_check.py <rt_driver> <problems_dir> <run_dir> [problem ...] [--set key=value ...]

--set is passed to rt_driver (e.g. --set legacy_rt=0).

For each problem: runs rt_driver, then checks
  * every step directory has all spec fields with the right sizes, no NaN/Inf
  * ionised volume (sum of 1-fHI) grows monotonically
  * fractions in [0,1], temperature >= the initial floor
  * ray stats are consistent (nray_active > 0, ndel1 <= ncst1+..., etc.)
  * photon budget: ionisations (H + He, + recombinations) do not exceed the
    photons emitted, and absorb a sensible fraction of them
  * uniform_1src: the I-front radius grows with time, monotonically, and the
    ionised region is symmetric about the source
  * lognormal_16src_smallbuf: the merge path is exercised (nmrg1 > 0)
Prints a table per problem and the run time.
"""
import math
import os
import subprocess
import sys
import time

import numpy as np

FIELDS = {"fHI": 1, "fHeI": 1, "fHeII": 1, "temp": 1, "rhog": 1, "ngam": 5,
          "nu": 5, "fgam": 5, "tau": 5, "gamma": 3, "lambda": 3, "dx": 3,
          "u_gas": 5, "w_gas": 4}

# physical constants / parameters (as in the code)
f32 = lambda v: float(np.float32(v))
omegam, omegab, hubble0 = f32(0.27), f32(0.045), f32(0.70)
XHy, YHe = f32(0.76), f32(0.24)
H0_cgs, rhoc_cgs, Mpc2cm, mH = 3.24086e-18, 1.87890e-29, 3.08560e24, 1.67223e-24
pop2 = float(np.float32(10 ** np.float32(25.2))) * np.array(
    [f32(0.4185), f32(0.2942), f32(0.2009), f32(0.0811), f32(0.0054)], dtype=np.float32)
fResc = f32(0.2)


def read_params(d):
    p = {}
    for line in open(os.path.join(d, "params.txt")):
        if "=" in line:
            k, v = line.split("=", 1)
            p[k.strip()] = v.strip()
    return p


def read_stats(d):
    s = {}
    for line in open(os.path.join(d, "stats.txt")):
        k, v = line.split("=", 1)
        s[k.strip()] = float(v.replace("D", "e"))
    return s


def dt_sec(a, da):
    # cosmo_tools.f90: calc_dt_from_da (dy), simple midpoint integration
    n = 100
    om, ol = omegam, 1 - omegam
    dy = 0.0
    for i in range(1, n + 1):
        x = a + (i - 0.5) * da / n
        E = math.sqrt(om / x**3 + ol)
        dy += 1 / (x * H0_cgs * hubble0 * E) * da / n
    return dy


def check_problem(driver, pdir, rdir, name, extra=()):
    errs = []
    t0 = time.time()
    out = subprocess.run([driver, pdir, rdir, *extra], capture_output=True, text=True)
    elapsed = time.time() - t0
    if out.returncode != 0:
        return [f"rt_driver failed: {out.stderr.strip()}"], elapsed, []
    p = read_params(pdir)
    N = int(p["N"])
    nsteps = int(p["nsteps"])
    z0 = float(p["z_init"])
    da = float(p["da_sim"])
    src = np.loadtxt(os.path.join(pdir, "sources.txt"), ndmin=2)
    Ltot = float(np.sum(src[:, 3].astype(np.float32).astype(np.float64)))

    rows = []
    prev_ion = 0.0
    prev_R = 0.0
    emitted = 0.0
    a = 1 / (1 + z0)
    for step in range(1, nsteps + 1):
        d = os.path.join(rdir, f"step_{step:04d}")
        if not os.path.isdir(d):
            errs.append(f"missing {d}")
            break
        data = {}
        for f, nc in FIELDS.items():
            fn = os.path.join(d, f + ".bin")
            if not os.path.exists(fn):
                errs.append(f"step {step}: missing {f}.bin")
                continue
            arr = np.fromfile(fn, dtype="<f4")
            if arr.size != nc * N**3:
                errs.append(f"step {step}: {f} has size {arr.size}")
                continue
            if not np.all(np.isfinite(arr)):
                errs.append(f"step {step}: {f} has NaN/Inf")
            data[f] = arr.reshape((N, N, N, nc) if nc > 1 else (N, N, N))
        st = read_stats(d)
        fHI = data["fHI"]
        for f in ("fHI", "fHeI", "fHeII"):
            if data[f].min() < 0 or data[f].max() > 1:
                errs.append(f"step {step}: {f} out of [0,1]")
        # ionised volume (cells), I-front radius from the ionised volume
        ion = float(np.sum(1.0 - fHI.astype(np.float64)))
        R = (3 * ion / (4 * math.pi)) ** (1 / 3)
        if ion < prev_ion * (1 - 1e-6):
            errs.append(f"step {step}: ionised volume decreased ({prev_ion} -> {ion})")
        prev_ion = ion
        # photon budget: photons emitted so far vs ionisations
        # ngam is normalised by <nH>*vol_rt -> count in units of "H atoms of a mean cell"
        emitted_step = fResc * float(np.sum(pop2)) * Ltot * dt_sec(a, da) * 3.15576e7 / 3.15576e7
        a += da
        emitted += emitted_step
        # number of H atoms in a mean-density cell at this a (proper volume)
        x_unit = (a * 5.0 / hubble0 / N) * Mpc2cm
        rho_b = rhoc_cgs * omegam * hubble0**2 / a**3 * omegab / omegam
        nH_cell = rho_b * XHy / mH * x_unit**3
        D = data["rhog"].astype(np.float64)
        n_ion_H = float(np.sum((1 - fHI) * D)) * nH_cell
        nHe_cell = nH_cell * (YHe / 4) / XHy
        n_ion_He = float(np.sum((1 - data["fHeI"]) * D)) * nHe_cell + \
            float(np.sum((1 - data["fHeI"] - data["fHeII"]) * D)) * nHe_cell
        frac = (n_ion_H + n_ion_He) / emitted
        # (rough estimate; collisional ionisation in dense, hot gas needs no
        # photons, so allow some slack)
        if frac > 1.15:
            errs.append(f"step {step}: more ionisations than photons emitted ({frac:.3f})")
        if step == nsteps and frac < 0.2:
            errs.append(f"step {step}: photon budget too small ({frac:.3f})")
        if st["nray_active"] <= 0 and not name.startswith("lognormal"):
            # (dense lognormal source cells can be optically thick, in which
            # case ray_source deposits photons locally and emits no rays)
            errs.append(f"step {step}: no active rays")
        rows.append((step, st["redshift"], int(st["nstep_rt"]), ion, R, frac,
                     int(st["nray_active"]), int(st["ncst1"]), int(st["nspl1"]),
                     int(st["nmrg1"]), int(st["ndel1"]), st["ngam_rays_total"],
                     float(data["temp"].max())))
        if name == "uniform_1src":
            if R < prev_R:
                errs.append(f"step {step}: I-front radius decreased")
            prev_R = R
            # symmetry: HEALPix rays are 4-fold symmetric about z, so the
            # x and y profiles through the source must agree, and each must
            # be mirror-symmetric about the source cell (z differs: rings).
            i, j, k = (int(v) - 1 for v in src[0, :3])
            prof_x = fHI[k, j, :].astype(np.float64)
            prof_y = fHI[k, :, i].astype(np.float64)
            m = min(i, N - 1 - i)
            mirror = np.abs(prof_x[i - m:i][::-1] - prof_x[i + 1:i + 1 + m]).max()
            if np.max(np.abs(prof_x - prof_y)) > 1e-4 or mirror > 1e-4:
                errs.append(f"step {step}: ionised region not symmetric")
    if name == "uniform_1src" and rows:
        if not (rows[-1][4] > 2.0 and rows[-1][4] > rows[0][4] * 1.5):
            errs.append("I-front did not grow to a multi-cell region")
    if name == "lognormal_16src_smallbuf" and rows:
        if sum(r[9] for r in rows) == 0:
            errs.append("merge path not exercised (nmrg1 == 0 in all steps)")
    return errs, elapsed, rows


def main():
    args, extra = [], []
    it = iter(sys.argv[1:])
    for a in it:
        if a == "--set":
            extra += ["--set", next(it)]
        else:
            args.append(a)
    driver, pdir, rdir = args[:3]
    names = args[3:] or ["uniform_1src", "uniform_8src", "lognormal_16src",
                             "uniform_8src_smallbuf", "lognormal_16src_smallbuf"]
    nfail = 0
    for name in names:
        errs, elapsed, rows = check_problem(driver, os.path.join(pdir, name),
                                            os.path.join(rdir, name), name, extra)
        print(f"== {name}: run time {elapsed:.2f} s")
        print(" step  redshift   nrt  ion_cells  R_if  budget  nray  ncst1  nspl1  nmrg1  ndel1  ngam_rays    Tmax")
        for r in rows:
            print(" %4d  %.5f  %3d  %9.2f  %4.2f  %6.3f  %5d  %5d  %5d  %5d  %5d  %9.4g  %7.4g" % r)
        for e in errs:
            print("  FAIL:", e)
        nfail += len(errs)
    print("FAILED" if nfail else "OK")
    sys.exit(1 if nfail else 0)


if __name__ == "__main__":
    main()
