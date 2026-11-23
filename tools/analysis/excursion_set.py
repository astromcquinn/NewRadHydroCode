#!/usr/bin/env python3
"""The 21cmFAST excursion-set criterion on the gas density of the 512^3
RT-methods suite (B1 of docs/taueff_21cmfast.md).

    excursion_set.py run RUN --steps 61,69,74,78,84 --out DIR [--rmax 15]
                     [--rmax-mfp IGM.json] [--mmin 1.28e8] [--threads N]
    excursion_set.py run-mass RUN --steps ... --out DIR [--alpha-star A] [--pop2-start P]
    excursion_set.py standalone RUN --steps ... --out DIR [--hii-dim 512] [--dim 1024]
    excursion_set.py figure DIR RUN:label [RUN:label ...] --outdir FIGDIR [--es NAME,...]
    excursion_set.py --selftest

`run-mass` (B1m) uses 21cmFAST's mass-dependent sources (E-INTEGRAL) with
the suite's abundance-matched slope, N_ion ~ M^a(z), and solves POP2_ION for
x_V. `run` needs py21cmfast (tested with 4.3.1, the conda env of
docs/taueff_21cmfast.md); `figure` and the self-test need numpy and scipy only.

The criterion, as 21cmFAST 4.3.1 applies it with SOURCE_MODEL = CONST-ION-EFF
(src/py21cmfast/src/IonisationBox.c: setup_radii, copy_filter_transform,
calculate_conditional_nion, find_ionised_regions; hmf.c:
nion_conditional_acg, fcoll_conditional_eps; filtering.c: real_tophat_filter;
Constants.c: delta_c_sph = 1.686, l_factor = 0.620350491):

* the density contrast delta of each cell is smoothed with a real-space top
  hat of radius R (in k space, W = 3 (sin kR - kR cos kR)/(kR)^3), for R from
  R_max = min(R_BUBBLE_MAX, l_factor BOX_LEN) down to R_min =
  max(R_BUBBLE_MIN, l_factor dx) in factors of DELTA_R_HII_FACTOR = 1.1 (the
  smallest R uses the unsmoothed delta), delta_R clipped to > -1;
* f_coll(R) = erfc[(delta_c - delta_R)/D(z) / sqrt(2 (sigma^2(M_min) -
  sigma^2(M_R)))] (extended Press-Schechter, sigma at z = 0, M_R the mass
  of the sphere R at the mean matter density);
* n_ion(R) = zeta f_coll(R) times the factor that sets its box mean to the
  unconditional zeta f_coll(> M_min) of the HMF (ST by default; fix_mean);
* a cell is ionised if n_ion(R) > 1 on any R ("central pixel",
  IONISE_ENTIRE_SPHERE = False); otherwise x_HII = n_ion(R_min), clipped
  to [0, 1].

`run` puts the suite's gas density (rhog / <rhog> - 1 of the reference run's
step) into a 21cmFAST PerturbedField (the FAQ route "Can I instantiate my own
OutputStruct objects?": new() and set(); the density is then marked computed),
passes a dummy InitialConditions (not read by CONST-ION-EFF), and calls
compute_ionization_field. zeta is chosen so that the volume-weighted x_HII
equals the reference run's at that step (a numpy mirror of the criterion,
excursion(), gives the starting value; the secant method on 21cmFAST itself
finishes). M_min is the suite's smallest source halo by default (20 particles
of 6.4e6 Msun), and R_BUBBLE_MIN is raised to the radius of M_min (21cmFAST's
EPS f_coll is undefined where sigma(M_R) > sigma(M_min)). One variant per
R_max: the 21cmFAST default 15 Mpc and, with --rmax-mfp, the reference run's
lambda_mfp (ionised starts, comoving) at each z.

Output: DIR/es_<variant>/{params.txt, run/step_NNNN/{fHI.bin, stats.txt}}
(the ionised fraction as a neutral-fraction cube, readable like a run) and
DIR/es_<variant>/zeta.txt.

`figure` compares every run (the first is the reference) and every es_*
variant of DIR at the reference's steps: Delta^2 of x_HII, the ratio to the
reference, r(k) with it, the fraction of cells whose x_HII differs from the
reference's by more than 0.5, and the x_HII maps (the slab of
rt_suite_figures.py). `--es rmax15` shows only that variant (labelled
"21cmFAST"); B2 only if `standalone` is listed.
"""
import argparse
import glob
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

DELTA_C = 1.686           # Constants.c delta_c_sph
L_FACTOR = 0.620350491    # Constants.c l_factor, (4 pi / 3)^(-1/3)
FRACT_FLOAT_ERR = 1e-7    # Constants.h


# ---------------------------------------------------------------------------
# The criterion in numpy (mirror of IonisationBox.c, for zeta, the cross-check
# and the self-test)
# ---------------------------------------------------------------------------

def radii_21cmfast(rmin, rmax, factor=1.1):
    """IonisationBox.c setup_radii: R_i = rmin factor^i, the last one = rmax."""
    n = int(math.log(rmax / rmin) / math.log(factor) + 1)
    out = []
    for i in range(n):
        R = rmin * factor ** i
        if R > rmax - FRACT_FLOAT_ERR:
            out.append(rmax)
            break
        out.append(R)
    return np.array(out)


def tophat_k(N, box, R):
    """filtering.c real_tophat_filter on the rfftn grid (box and R in the
    same length unit)."""
    kf = 2 * np.pi / box
    n = np.fft.fftfreq(N, 1.0 / N)
    nz = np.arange(N // 2 + 1)
    k = kf * np.sqrt(n[:, None, None] ** 2 + n[None, :, None] ** 2 + nz[None, None, :] ** 2)
    kR = k * R
    with np.errstate(divide="ignore", invalid="ignore"):
        W = 3.0 * (np.sin(kR) - kR * np.cos(kR)) / kR ** 3
    return np.where(kR < 1e-4, 1 - kR * kR / 10, W)


def fcoll_eps(delta_R, growth, sig_min, sig_R):
    """hmf.c fcoll_conditional_eps."""
    from scipy.special import erfc
    if sig_R >= sig_min:
        raise ValueError("sigma(M_R) >= sigma(M_min): R below the radius of M_min")
    s = math.sqrt(sig_min ** 2 - sig_R ** 2)
    return erfc((DELTA_C - delta_R) / growth / (math.sqrt(2) * s))


def excursion(delta, box, radii, nion, unconditional=None, workers=1):
    """Max over R of the (mean-fixed) n_ion per unit zeta, and its value at
    the smallest R. delta (N^3), box and radii in one length unit; nion(dR,
    i) the n_ion per unit zeta of the smoothed field dR at radius radii[i];
    unconditional: the box mean every R is rescaled to (None: no rescaling,
    fix_mean off). A cell is ionised for zeta * gmax > 1, otherwise x_HII =
    clip(zeta * g0, 0, 1) (find_ionised_regions)."""
    import scipy.fft as sfft
    N = delta.shape[0]
    dk = sfft.rfftn(delta.astype(np.float64), workers=workers)
    gmax = np.zeros(delta.shape, np.float32)
    g0 = None
    fixes = []
    for i in range(len(radii) - 1, -1, -1):  # largest R first, as 21cmFAST
        if i > 0:
            dR = sfft.irfftn(dk * tophat_k(N, box, radii[i]), s=delta.shape, workers=workers)
            dR = np.clip(dR, -1.0, 1e6)
        else:
            dR = delta.astype(np.float64)
        dR = np.maximum(dR, -1.0 + FRACT_FLOAT_ERR)
        n = nion(dR, i)
        fix = unconditional / n.mean() if unconditional is not None else 1.0
        fixes.append(fix)
        np.maximum(gmax, (fix * n).astype(np.float32), out=gmax)
        if i == 0:
            g0 = (fix * n).astype(np.float32)
        del dR, n
    return gmax, g0, fixes[::-1]


def xhii_of_zeta(gmax, g0, zeta):
    return np.where(zeta * gmax > 1, np.float32(1), np.clip(zeta * g0, 0, 1).astype(np.float32))


def zeta_for_xv(gmax, g0, target, lo=1e-3, hi=1e6):
    """Bisection in log zeta on the volume-weighted x_HII."""
    for _ in range(80):
        mid = math.sqrt(lo * hi)
        if xhii_of_zeta(gmax, g0, mid).mean() < target:
            lo = mid
        else:
            hi = mid
        if hi / lo < 1 + 1e-6:
            break
    return math.sqrt(lo * hi)


# ---------------------------------------------------------------------------
# Runs
# ---------------------------------------------------------------------------

def read_stats(d):
    out = {}
    for line in open(os.path.join(d, "stats.txt")):
        if " = " in line:
            k, v = line.split(" = ", 1)
            try:
                out[k.strip()] = float(v)
            except ValueError:
                out[k.strip()] = v.strip()
    return out


def run_dir(path):
    path = os.path.abspath(path.rstrip("/"))
    return os.path.join(path, "run") if os.path.isdir(os.path.join(path, "run")) else path


def box_of(rd):
    p = os.path.join(os.path.dirname(rd), "params.txt")
    for line in open(p):
        k, _, v = line.partition("=")
        if k.strip() == "box":
            return float(v.split("#")[0])
    return None


def load_cube(path):
    a = np.fromfile(path, dtype=np.float32)
    N = round(a.size ** (1 / 3))
    return a.reshape(N, N, N)


def snapshots(rd):
    """{step: (z, dir)} of the steps with fHI.bin."""
    out = {}
    for d in sorted(glob.glob(os.path.join(rd, "step_[0-9][0-9][0-9][0-9]"))):
        if os.path.exists(os.path.join(d, "fHI.bin")):
            out[int(d[-4:])] = (read_stats(d)["redshift"], d)
    return out


# ---------------------------------------------------------------------------
# 21cmFAST
# ---------------------------------------------------------------------------

def p21_inputs(N, box_mpc, zeta, rmax, rmin, log10_mmin, threads, cosmo, dim=None, seed=1):
    import py21cmfast as p21
    inp = p21.InputParameters.from_template("const-zeta", random_seed=seed)
    return inp.evolve_input_structs(
        HII_DIM=N, DIM=dim or N, BOX_LEN=box_mpc, N_THREADS=threads,
        OMm=cosmo["omegam"], OMb=cosmo["omegab"], hlittle=cosmo["hubble0"],
        SIGMA_8=cosmo["sigma8"], POWER_INDEX=cosmo["ns"],
        HII_EFF_FACTOR=zeta, R_BUBBLE_MAX=rmax, R_BUBBLE_MIN=rmin,
        M_MIN_in_Mass=True, M_TURN_STELLAR_FEEDBACK=log10_mmin,
        HII_FILTER="spherical-tophat", IONISE_ENTIRE_SPHERE=False)


def p21_ionize(delta, z, inp):
    """compute_ionization_field on our density; returns the IonizedBox."""
    import py21cmfast as p21
    pf = p21.PerturbedField.new(inputs=inp, redshift=z)
    pf.set("density", np.ascontiguousarray(delta, dtype=np.float32))
    pf.set("velocity_z", np.zeros(delta.shape, np.float32))
    for k in ("density", "velocity_z"):  # set() leaves them "initialized"
        setattr(pf, k, pf.arrays[k].computed())
    ic = p21.InitialConditions.new(inputs=inp, dummy=True)
    return p21.compute_ionization_field(perturbed_field=pf, initial_conditions=ic, inputs=inp)


def p21_sigmas(inp, z, mmin, radii):
    """sigma(M_min), sigma(M_R) at z = 0 and D(z) from the 21cmFAST backend."""
    from py21cmfast.wrapper import cfuncs
    D = cfuncs.get_growth_factor(inputs=inp, redshift=z)  # also initialises the C state
    rho = cfuncs.get_condition_mass(inp, 1.0) / (4.0 / 3.0 * np.pi)
    m = np.concatenate([[mmin], 4.0 / 3.0 * np.pi * radii ** 3 * rho])
    sig, _ = cfuncs.evaluate_sigma(inputs=inp, masses=m)
    return float(D), float(sig[0]), sig[1:], rho


COSMO_DEFAULT = {"omegam": 0.27, "omegab": 0.045, "hubble0": 0.70, "sigma8": 0.80, "ns": 0.96}


def read_cosmo(top):
    """params.txt omegam, omegab, hubble0, sigma8, ns (the code's defaults
    otherwise, include/radray/cosmo_params.hpp)."""
    c = dict(COSMO_DEFAULT)
    for line in open(os.path.join(top, "params.txt")):
        k, _, v = line.partition("=")
        if k.strip() in c:
            c[k.strip()] = float(v.split("#")[0])
    return c


def cmd_run(a):
    import time
    rd = run_dir(a.run)
    top = os.path.dirname(rd)
    cosmo = read_cosmo(top)
    h = cosmo["hubble0"]
    box_h = box_of(rd)
    box = box_h / h                      # Mpc (21cmFAST lengths have no h)
    snaps = snapshots(rd)
    steps = [int(s) for s in a.steps.split(",")]
    mfp = {}
    if a.rmax_mfp:
        import json
        d = json.load(open(a.rmax_mfp))
        for s in d["snapshots"]:
            m = s["mfp"].get("ionized", s["mfp"]["random"])
            mfp[s["step"]] = m["lambda_fit"] * (1 + s["z"])  # comoving Mpc
    variants = [("rmax%g" % a.rmax, lambda st: a.rmax)]
    if mfp:
        variants.append(("rmaxmfp", lambda st: mfp[st]))
    threads = a.threads or int(os.environ.get("SLURM_CPUS_PER_TASK", "1"))
    log10_mmin = math.log10(a.mmin)
    for st in steps:
        z, sd = snaps[st]
        rho = load_cube(os.path.join(sd, "rhog.bin")).astype(np.float64)
        delta = (rho / rho.mean() - 1.0).astype(np.float32)
        del rho
        N = delta.shape[0]
        target = float(1.0 - load_cube(os.path.join(sd, "fHI.bin")).astype(np.float64).mean())
        for vname, rmax_of in variants:
            t0 = time.time()
            rmax = min(rmax_of(st), L_FACTOR * box)
            # R_min: the cell (l_factor dx) or the radius of M_min, whichever is larger
            inp0 = p21_inputs(N, box, 30.0, rmax, L_FACTOR * box / N, log10_mmin, threads, cosmo)
            from py21cmfast.wrapper import cfuncs
            cfuncs.get_growth_factor(inputs=inp0, redshift=z)
            rho_m = cfuncs.get_condition_mass(inp0, 1.0) / (4.0 / 3.0 * np.pi)
            r_mmin = (a.mmin / (4.0 / 3.0 * np.pi * rho_m)) ** (1 / 3)
            rmin = max(L_FACTOR * box / N, 1.01 * r_mmin)
            radii = radii_21cmfast(rmin, rmax)
            inp = p21_inputs(N, box, 30.0, rmax, rmin, log10_mmin, threads, cosmo)
            D, smin, sR, _ = p21_sigmas(inp, z, a.mmin, radii)
            ib = p21_ionize(delta, z, inp)
            unc = float(ib.nion_unconditional_acg) / 30.0   # per unit zeta
            # numpy mirror: the starting zeta and, later, the cross-check
            gmax, g0, fixes = excursion(
                delta, box, radii, lambda dR, i: fcoll_eps(dR, D, smin, sR[i]), unc, workers=threads)
            zeta_np = zeta_for_xv(gmax, g0, target)
            # secant in log zeta on 21cmFAST itself
            hist = []

            def xv_of(zeta):
                ibz = p21_ionize(delta, z, p21_inputs(N, box, zeta, rmax, rmin, log10_mmin,
                                                      threads, cosmo))
                x = float(1.0 - np.asarray(ibz.neutral_fraction, np.float64).mean())
                hist.append((zeta, x))
                print(f"  step {st} z {z:.3f} {vname}: zeta {zeta:.5g} -> x_V {x:.5f} "
                      f"(target {target:.5f})", flush=True)
                return x, ibz

            z1 = zeta_np
            x1, ib1 = xv_of(z1)
            z0, x0 = z1 * (1.05 if x1 < target else 0.95), None
            x0, ib0 = xv_of(z0)
            best = min(((x1, z1, ib1), (x0, z0, ib0)), key=lambda t: abs(t[0] - target))
            for _ in range(8):
                if abs(best[0] - target) < a.tol or x1 == x0:
                    break
                lz = math.log(z1) + (target - x1) * (math.log(z1) - math.log(z0)) / (x1 - x0)
                z0, x0 = z1, x1
                z1 = math.exp(lz)
                x1, ib1 = xv_of(z1)
                if abs(x1 - target) < abs(best[0] - target):
                    best = (x1, z1, ib1)
            xv, zeta, ibz = best
            xhii = (1.0 - np.asarray(ibz.neutral_fraction, np.float32))
            xnp = xhii_of_zeta(gmax, g0, zeta)
            mism = float(np.mean((xhii > 0.999) != (xnp > 0.999)))
            dmean = float(np.abs(xhii - xnp).mean())
            outd = os.path.join(a.out, "es_" + vname)
            sdir = os.path.join(outd, "run", f"step_{st:04d}")
            os.makedirs(sdir, exist_ok=True)
            with open(os.path.join(outd, "params.txt"), "w") as f:
                f.write(f"# excursion_set.py run (21cmFAST criterion on {top})\nN = {N}\nbox = {box_h}\n")
            (1.0 - xhii).astype(np.float32).tofile(os.path.join(sdir, "fHI.bin"))
            with open(os.path.join(sdir, "stats.txt"), "w") as f:
                f.write(f"redshift = {z}\nxHII_V = {xv}\nzeta = {zeta}\ntarget_xHII_V = {target}\n"
                        f"R_max_Mpc = {rmax}\nR_min_Mpc = {rmin}\nn_radii = {len(radii)}\n"
                        f"M_min_Msun = {a.mmin}\nsigma_Mmin = {smin}\ngrowth = {D}\n"
                        f"nion_unconditional_per_zeta = {unc}\nzeta_numpy = {zeta_np}\n"
                        f"numpy_mismatch_ionised = {mism}\nnumpy_mean_abs_dx = {dmean}\n"
                        f"source_step = {sd}\n")
            zt = os.path.join(outd, "zeta.txt")
            new = not os.path.exists(zt)
            with open(zt, "a") as f:
                if new:
                    f.write("# step z target_xV xV zeta zeta_numpy R_max[Mpc] R_min[Mpc] n_radii "
                            "M_min[Msun] mismatch_numpy mean|dx|_numpy fix(R_min) fix(R_max) seconds\n")
                f.write(f"{st} {z:.4f} {target:.5f} {xv:.5f} {zeta:.5g} {zeta_np:.5g} {rmax:.4g} "
                        f"{rmin:.4g} {len(radii)} {a.mmin:.4g} {mism:.3e} {dmean:.3e} "
                        f"{fixes[0]:.4g} {fixes[-1]:.4g} {time.time() - t0:.0f}\n")
            print(f"step {st} z {z:.3f} {vname}: zeta {zeta:.4g} (numpy {zeta_np:.4g}), x_V {xv:.4f} "
                  f"target {target:.4f}; numpy vs 21cmFAST: {mism:.2e} of the cells differ in "
                  f"'ionised', <|dx|> {dmean:.2e}", flush=True)
            del gmax, g0, xhii, xnp, ib, ibz


def am_slope(top, z):
    """The low-mass slope a(z) of the suite's abundance-matched L_UV(M)
    (docs/abundance_matching.md: L = 2 L_p / ((M/M_p)^-a + (M/M_p)^-b),
    a = a0 + a1 x, x = clamp(z, z_lo, z_hi) - z0), from the run's params.txt."""
    p = {}
    for line in open(os.path.join(top, "params.txt")):
        if "=" in line and not line.lstrip().startswith("#"):
            k, v = [t.strip() for t in line.split("=", 1)]
            p[k] = v
    x = min(max(z, float(p["am_zlo"])), float(p["am_zhi"])) - float(p["am_z0"])
    return float(p["am_a0"]) + float(p["am_a1"]) * x


def p21_inputs_mass(N, box_mpc, pop2, rmax, rmin, log10_mturn, alpha_star, threads, cosmo,
                    alpha_esc=0.0):
    """21cmFAST's mass-dependent source model (SOURCE_MODEL = E-INTEGRAL):
    N_ion per halo ~ M f_*(M) f_esc(M) POP2_ION, f_* = F_STAR10 (M/1e10)^ALPHA_STAR,
    f_esc = F_ESC10 (M/1e10)^ALPHA_ESC, duty cycle exp(-M_turn/M), M_turn =
    max(M_TURN_STELLAR_FEEDBACK, the atomic-cooling mass); POP2_ION is the free
    normalisation (it scales N_ion linearly, like zeta). Otherwise as p21_inputs."""
    import py21cmfast as p21
    inp = p21.InputParameters.from_template("const-zeta", random_seed=1)
    return inp.evolve_input_structs(
        HII_DIM=N, DIM=N, BOX_LEN=box_mpc, N_THREADS=threads,
        OMm=cosmo["omegam"], OMb=cosmo["omegab"], hlittle=cosmo["hubble0"],
        SIGMA_8=cosmo["sigma8"], POWER_INDEX=cosmo["ns"],
        SOURCE_MODEL="E-INTEGRAL", ALPHA_STAR_ACG=alpha_star, ALPHA_ESC=alpha_esc,
        POP2_ION=pop2, R_BUBBLE_MAX=rmax, R_BUBBLE_MIN=rmin,
        M_MIN_in_Mass=True, M_TURN_STELLAR_FEEDBACK=log10_mturn,
        HII_FILTER="spherical-tophat", IONISE_ENTIRE_SPHERE=False)


def cmd_run_mass(a):
    """B1 with the mass-dependent source model: ALPHA_STAR = a(z) - 1 from the
    suite's abundance matching (N_ion ~ M^a below the turnover M_p ~ 4e11),
    ALPHA_ESC = 0 (the suite's f_esc depends on z only), M_turn = M_min;
    POP2_ION solved so that x_V = the rays'. Writes OUT/es_<name>/ as `run`."""
    import time
    rd = run_dir(a.run)
    top = os.path.dirname(rd)
    cosmo = read_cosmo(top)
    box_h = box_of(rd)
    box = box_h / cosmo["hubble0"]
    snaps = snapshots(rd)
    threads = a.threads or int(os.environ.get("SLURM_CPUS_PER_TASK", "1"))
    log10_mt = math.log10(a.mmin)
    outd = os.path.join(a.out, "es_" + a.name)
    for st in [int(s) for s in a.steps.split(",")]:
        t0 = time.time()
        z, sd = snaps[st]
        rho = load_cube(os.path.join(sd, "rhog.bin")).astype(np.float64)
        delta = (rho / rho.mean() - 1.0).astype(np.float32)
        del rho
        N = delta.shape[0]
        target = float(1.0 - load_cube(os.path.join(sd, "fHI.bin")).astype(np.float64).mean())
        alpha = (am_slope(top, z) - 1.0) if a.alpha_star is None else a.alpha_star
        rmax = min(a.rmax, L_FACTOR * box)
        rmin = a.rmin                       # B1's R_min (radius of 1.28e8 Msun x 1.01)

        def xv_of(pop2):
            inp = p21_inputs_mass(N, box, pop2, rmax, rmin, log10_mt, alpha, threads, cosmo)
            ib = p21_ionize(delta, z, inp)
            x = float(1.0 - np.asarray(ib.neutral_fraction, np.float64).mean())
            print(f"  step {st} z {z:.3f}: POP2_ION {pop2:.5g} -> x_V {x:.5f} (target {target:.5f})",
                  flush=True)
            return x, ib

        xv, pop2, ib = solve_zeta(xv_of, target, a.pop2_start, a.tol, factor=1.5, maxcall=a.maxcall)
        xhii = (1.0 - np.asarray(ib.neutral_fraction, np.float32))
        sdir = os.path.join(outd, "run", f"step_{st:04d}")
        os.makedirs(sdir, exist_ok=True)
        with open(os.path.join(outd, "params.txt"), "w") as f:
            f.write(f"# excursion_set.py run-mass (21cmFAST E-INTEGRAL on {top})\nN = {N}\nbox = {box_h}\n")
        (1.0 - xhii).astype(np.float32).tofile(os.path.join(sdir, "fHI.bin"))
        with open(os.path.join(sdir, "stats.txt"), "w") as f:
            f.write(f"redshift = {z}\nxHII_V = {xv}\npop2_ion = {pop2}\ntarget_xHII_V = {target}\n"
                    f"alpha_star = {alpha}\nalpha_esc = 0\nM_turn_Msun = {a.mmin}\n"
                    f"R_max_Mpc = {rmax}\nR_min_Mpc = {rmin}\nsource_step = {sd}\n")
        with open(os.path.join(outd, f"pop2_step{st}.txt"), "w") as f:
            f.write("# step z target_xV xV POP2_ION alpha_star seconds\n"
                    f"{st} {z:.4f} {target:.5f} {xv:.5f} {pop2:.5g} {alpha:.4f} {time.time() - t0:.0f}\n")
        print(f"step {st} z {z:.3f}: POP2_ION {pop2:.4g}, alpha_star {alpha:.3f}, x_V {xv:.4f} "
              f"target {target:.4f}", flush=True)


def solve_zeta(xv_of, target, zeta0, tol, factor=1.3, maxcall=16):
    """zeta with x_V(zeta) = target (x_V rises with zeta): expand from
    zeta0 by `factor` until the target is bracketed, then secant steps in
    ln zeta kept inside the bracket (bisection when a step leaves it).
    xv_of(zeta) returns (x_V, result). Returns (x_V, zeta, result) of the
    closest call."""
    calls = []

    def call(z):
        x, r = xv_of(z)
        calls.append((x, z, r))
        return x

    lo = hi = None                      # (ln zeta, x) below / above the target
    lz = math.log(zeta0)
    x = call(math.exp(lz))
    while len(calls) < maxcall:
        if abs(x - target) < tol:
            break
        if x < target:
            lo = (lz, x)
        else:
            hi = (lz, x)
        if lo is None:
            lz = hi[0] - math.log(factor)
        elif hi is None:
            lz = lo[0] + math.log(factor)
        else:
            # secant between the bracket ends; bisection if degenerate or outside
            if hi[1] > lo[1]:
                lz = lo[0] + (target - lo[1]) * (hi[0] - lo[0]) / (hi[1] - lo[1])
            if not (min(lo[0], hi[0]) < lz < max(lo[0], hi[0])) or hi[1] <= lo[1] or \
                    abs(lz - lo[0]) < 1e-3 * abs(hi[0] - lo[0]) or abs(hi[0] - lz) < 1e-3 * abs(hi[0] - lo[0]):
                lz = 0.5 * (lo[0] + hi[0])
            if abs(hi[0] - lo[0]) < 1e-5:
                break
        x = call(math.exp(lz))
    return min(calls, key=lambda t: abs(t[0] - target))


def cmd_standalone(a):
    """B2: a 21cmFAST run of its own (its ICs, 2LPT density at each z), the
    suite's cosmology, box and resolution; zeta at each z matched to the
    reference run's x_V there (coeval boxes: no evolution is needed without
    spin temperature and recombinations)."""
    import time
    import py21cmfast as p21
    p21.config["direc"] = os.path.join(os.path.abspath(a.out), "cache")  # never the home default
    rd = run_dir(a.run)
    top = os.path.dirname(rd)
    cosmo = read_cosmo(top)
    box_h = box_of(rd)
    box = box_h / cosmo["hubble0"]
    snaps = snapshots(rd)
    steps = [int(s) for s in a.steps.split(",")]
    threads = a.threads or int(os.environ.get("SLURM_CPUS_PER_TASK", "1"))
    N = a.hii_dim
    log10_mmin = math.log10(a.mmin)
    rmax = min(a.rmax, L_FACTOR * box)
    inp0 = p21_inputs(N, box, 30.0, rmax, L_FACTOR * box / N, log10_mmin, threads, cosmo,
                      dim=a.dim, seed=a.seed)
    from py21cmfast.wrapper import cfuncs
    cfuncs.get_growth_factor(inputs=inp0, redshift=7.0)
    rho_m = cfuncs.get_condition_mass(inp0, 1.0) / (4.0 / 3.0 * np.pi)
    rmin = max(L_FACTOR * box / N, 1.01 * (a.mmin / (4.0 / 3.0 * np.pi * rho_m)) ** (1 / 3))

    def inputs(zeta):
        return p21_inputs(N, box, zeta, rmax, rmin, log10_mmin, threads, cosmo, dim=a.dim, seed=a.seed)

    t0 = time.time()
    zeta_prev = None
    ic = p21.compute_initial_conditions(inputs=inputs(30.0))
    print(f"ICs {a.dim}^3 -> {N}^3, seed {a.seed}: {time.time() - t0:.0f} s", flush=True)
    outd = os.path.join(a.out, "es_standalone")
    os.makedirs(outd, exist_ok=True)
    with open(os.path.join(outd, "params.txt"), "w") as f:
        f.write(f"# excursion_set.py standalone: 21cmFAST {p21.__version__}, own ICs (seed {a.seed}, "
                f"DIM {a.dim}), 2LPT\nN = {N}\nbox = {box_h}\n")
    for st in steps:
        z, sd = snaps[st]
        target = float(1.0 - load_cube(os.path.join(sd, "fHI.bin")).astype(np.float64).mean())
        t1 = time.time()
        pf = p21.perturb_field(redshift=z, initial_conditions=ic, inputs=inputs(30.0))
        hist = []

        def xv_of(zeta):
            ib = p21.compute_ionization_field(perturbed_field=pf, initial_conditions=ic,
                                              inputs=inputs(zeta))
            x = float(1.0 - np.asarray(ib.neutral_fraction, np.float64).mean())
            hist.append((zeta, x))
            print(f"  step {st} z {z:.3f}: zeta {zeta:.5g} -> x_V {x:.5f} (target {target:.5f})",
                  flush=True)
            return x, ib

        xv, zeta, ib = solve_zeta(xv_of, target, zeta_prev or a.zeta0, a.tol)
        zeta_prev = zeta
        sdir = os.path.join(outd, "run", f"step_{st:04d}")
        os.makedirs(sdir, exist_ok=True)
        np.asarray(ib.neutral_fraction, np.float32).tofile(os.path.join(sdir, "fHI.bin"))
        (1.0 + np.asarray(pf.density, np.float32)).tofile(os.path.join(sdir, "rhog.bin"))
        with open(os.path.join(sdir, "stats.txt"), "w") as f:
            f.write(f"redshift = {z}\nxHII_V = {xv}\nzeta = {zeta}\ntarget_xHII_V = {target}\n"
                    f"R_max_Mpc = {rmax}\nR_min_Mpc = {rmin}\nM_min_Msun = {a.mmin}\n"
                    f"seed = {a.seed}\nDIM = {a.dim}\n")
        zt = os.path.join(outd, "zeta.txt")
        new = not os.path.exists(zt)
        with open(zt, "a") as f:
            if new:
                f.write("# step z target_xV xV zeta R_max[Mpc] R_min[Mpc] M_min[Msun] seconds\n")
            f.write(f"{st} {z:.4f} {target:.5f} {xv:.5f} {zeta:.5g} {rmax:.4g} {rmin:.4g} "
                    f"{a.mmin:.4g} {time.time() - t1:.0f}\n")
        print(f"step {st} z {z:.3f}: zeta {zeta:.4g}, x_V {xv:.4f} (target {target:.4f})", flush=True)
        del pf, ib


# ---------------------------------------------------------------------------
# Figure 8
# ---------------------------------------------------------------------------

def cmd_figure(a):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import gamma_power as gp
    from rt_suite_figures import COLORS, STYLES

    specs = []
    for r in a.runs:
        p, _, lab = r.partition(":")
        specs.append((run_dir(p), lab or os.path.basename(p.rstrip("/"))))
    pick = [x for x in a.es.split(",") if x] if a.es else None     # --es: the variants to show
    es = sorted(d for d in glob.glob(os.path.join(a.esdir, "es_*")) if not d.endswith("_standalone")
                and (pick is None or os.path.basename(d)[3:] in pick))
    if pick is not None and len(es) == 1 and "standalone" not in pick:
        es_specs = [(os.path.join(es[0], "run"), "21cmFAST")]
    else:
        es_specs = [(os.path.join(d, "run"), "21cmFAST ES" + ("" if i == 0 else f" ({os.path.basename(d)[3:]})"))
                    for i, d in enumerate(es)]
    # B2, the standalone 21cmFAST run: another realisation, so only its
    # Delta^2 and the ratio are shown (no r(k), no cell differences, no map)
    sa = os.path.join(a.esdir, "es_standalone")
    if os.path.isdir(os.path.join(sa, "run")) and (pick is None or "standalone" in pick):
        es_specs.append((os.path.join(sa, "run"), "21cmFAST own ICs (B2)"))
    all_specs = specs + es_specs
    nrun, nes = len(specs), len(es_specs)
    box = box_of(specs[0][0])
    ref = snapshots(specs[0][0])
    es_steps = set()
    for d, _ in es_specs:
        es_steps |= set(snapshots(d))
    steps = sorted(s for s in ref if s in es_steps) if a.steps == "" else \
        [int(s) for s in a.steps.split(",")]
    os.makedirs(a.outdir, exist_ok=True)
    colors = COLORS[:nrun] + ["#6250d6", "#9a8ad0", "#3a2a96"][:nes]
    styles = STYLES[:nrun] + ["-", ":", "--"][:nes]
    stat_only = [l.endswith("(B2)") for _, l in all_specs]
    for i, so in enumerate(stat_only):
        if so:
            colors[i], styles[i] = "#888888", "-."
    lw = [1.6] + [1.1] * (nrun - 1) + [1.8] * nes
    binner = None
    rows = {}     # step -> per spec (P_auto, P_cross, frac, slab, z, xV)
    for st in steps:
        z0, d0 = ref[st]
        res = []
        fr = None
        xr = None
        for i, (rd, lab) in enumerate(all_specs):
            sn = snapshots(rd)
            if i >= nrun:
                cand = {s: v for s, v in sn.items() if s == st}
            else:
                cand = sn
            if not cand:
                res.append(None)
                continue
            s_ = min(cand, key=lambda s: abs(cand[s][0] - z0))
            zz, d = cand[s_]
            if abs(zz - z0) > a.tol:
                res.append(None)
                continue
            x = 1.0 - load_cube(os.path.join(d, "fHI.bin")).astype(np.float64)
            if binner is None:
                binner = gp.Binner(x.shape[0], box, a.nbins)
            N = x.shape[0]
            lo, hi = N // 2 - a.slab // 2, N // 2 - a.slab // 2 + a.slab
            slab = x[lo:hi].mean(axis=0)
            fx = np.fft.rfftn(x)
            if i == 0:
                fr, xr = fx, x
            pa = binner.power(fx)
            pc = binner.power(fx, fr)
            frac = float(np.mean(np.abs(x - xr) > 0.5))
            res.append((pa, pc, frac, slab, zz, float(x.mean())))
            del x, fx
        rows[st] = res
        print(f"step {st} z {z0:.3f} done", flush=True)
    k = binner.k
    labels = [l for _, l in all_specs]
    # ---- table
    with open(os.path.join(a.outdir, "es_compare.txt"), "w") as f:
        f.write(f"# excursion_set.py figure: x_HII of every run and of the 21cmFAST excursion set "
                f"(B1, docs/taueff_21cmfast.md); box {box} Mpc/h\n")
        f.write("# fraction of cells with |x_HII - x_HII,ref| > 0.5 (ref = " + labels[0] + ")\n")
        f.write("# step z_ref " + " ".join(l.replace(" ", "_") for l in labels) + "\n")
        for st in steps:
            f.write(f"{st} {ref[st][0]:.4f} " + " ".join(
                "nan" if r is None else f"{r[2]:.5f}" for r in rows[st]) + "\n")
        f.write("# volume-weighted x_HII\n")
        for st in steps:
            f.write(f"{st} {ref[st][0]:.4f} " + " ".join(
                "nan" if r is None else f"{r[5]:.5f}" for r in rows[st]) + "\n")
        for st in steps:
            f.write(f"\n# step {st} z {ref[st][0]:.4f}: k [h/Mpc] then per run Delta2, P/P_ref, r\n")
            f.write("# runs: " + " ".join(l.replace(" ", "_") for l in labels) + "\n")
            pr = rows[st][0][0]
            cols = [k]
            for r in rows[st]:
                if r is None:
                    cols += [np.full_like(k, np.nan)] * 3
                    continue
                cols += [k ** 3 * r[0] / (2 * np.pi ** 2), r[0] / np.maximum(pr, 1e-300),
                         r[1] / np.sqrt(np.maximum(r[0] * pr, 1e-300))]
            for row in np.array(cols).T:
                f.write(" ".join(f"{v:.6g}" for v in row) + "\n")
    # ---- spectra
    plt.rcParams.update({"font.size": 9, "axes.titlesize": 9, "axes.labelsize": 9,
                         "xtick.labelsize": 8, "ytick.labelsize": 8, "figure.dpi": a.dpi})
    nc = len(steps)
    fig, axes = plt.subplots(3, nc + 1, figsize=(2.6 * (nc + 1) + 0.4, 6.6), squeeze=False)
    for ic, st in enumerate(steps):
        pr = rows[st][0][0]
        for ir, r in enumerate(rows[st]):
            if r is None:
                continue
            kw = dict(color=colors[ir], ls=styles[ir], lw=lw[ir], label=labels[ir])
            axes[0, ic].loglog(k, k ** 3 * r[0] / (2 * np.pi ** 2), **kw)
            if ir > 0:
                axes[1, ic].semilogx(k, r[0] / pr, **kw)
                if not stat_only[ir]:
                    axes[2, ic].semilogx(k, r[1] / np.sqrt(r[0] * pr), **kw)
        axes[0, ic].set_title(f"z = {ref[st][0]:.2f}, $x_V$ = {rows[st][0][5]:.2f}")
        axes[1, ic].axhline(1, color="0.6", lw=0.8)
        axes[1, ic].set_yscale("log")
        axes[1, ic].set_ylim(0.1, 10)
        axes[2, ic].set_ylim(0, 1.05)
        axes[2, ic].set_xlabel(r"$k$ [$h$ Mpc$^{-1}$]")
        for r_ in range(3):
            axes[r_, ic].sharex(axes[0, 0])
    axes[0, 0].set_ylabel(r"$\Delta^2_{x_{\rm HII}}(k)$")
    axes[1, 0].set_ylabel(r"$P / P_{\rm rays}$")
    axes[2, 0].set_ylabel(r"$r(k)$ with the rays")
    # last column: the fraction of cells differing by > 0.5, and the legend
    axf = axes[0, nc]
    zs = np.array([ref[st][0] for st in steps])
    for ir in range(1, len(labels)):
        if stat_only[ir]:
            continue
        v = np.array([np.nan if rows[st][ir] is None else rows[st][ir][2] for st in steps])
        axf.plot(zs, v, color=colors[ir], ls=styles[ir], lw=lw[ir], marker="o", ms=2.5)
    axf.set_xlabel("z")
    axf.set_ylabel(r"fraction $|x_{\rm HII} - x_{\rm HII,rays}| > 0.5$")
    axes[1, nc].axis("off")
    axes[2, nc].axis("off")
    h_, l_ = axes[0, 0].get_legend_handles_labels()
    axes[1, nc].legend(h_, l_, frameon=False, fontsize=7, loc="upper left")
    for ax in axes.ravel():
        ax.grid(True, which="major", color="0.9", lw=0.6)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(a.outdir, "es_pk_xHII.png"))
    plt.close(fig)
    # ---- maps
    mcols = [i for i in range(len(all_specs)) if not stat_only[i]]
    ncol = len(mcols)
    w = 1.9
    fig, axes = plt.subplots(nc, ncol, figsize=(w * ncol + 0.6, w * nc + 0.3), squeeze=False)
    im = None
    for irow, st in enumerate(steps):
        for icol, ir in enumerate(mcols):
            r = rows[st][ir]
            ax = axes[irow, icol]
            ax.set_xticks([])
            ax.set_yticks([])
            if r is None:
                ax.text(0.5, 0.5, "-", ha="center", va="center", transform=ax.transAxes)
                continue
            im = ax.imshow(r[3].T, origin="lower", extent=(0, box, 0, box), vmin=0, vmax=1,
                           cmap="cividis", interpolation="nearest")
            if irow == 0:
                ax.set_title(labels[ir], fontsize=8)
            if icol == 0:
                ax.set_ylabel(f"z = {ref[st][0]:.2f}, $x_V$ = {r[5]:.2f}")
    fig.subplots_adjust(left=0.04, right=0.92, bottom=0.02, top=0.95, wspace=0.03, hspace=0.03)
    cax = fig.add_axes([0.93, 0.1, 0.012, 0.8])
    fig.colorbar(im, cax=cax, label=r"$x_{\rm HII}$")
    fig.savefig(os.path.join(a.outdir, "es_maps_xHII.png"))
    plt.close(fig)
    print("wrote", a.outdir)


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------

def _overlap(R0, R, d):
    """Volume of the intersection of spheres of radii R0, R at distance d."""
    if d >= R0 + R:
        return 0.0
    if d <= abs(R0 - R):
        return 4 / 3 * math.pi * min(R0, R) ** 3
    return (math.pi * (R0 + R - d) ** 2 *
            (d * d + 2 * d * (R0 + R) - 3 * (R0 - R) ** 2) / (12 * d))


def selftest():
    """A single spherical overdensity delta0 of radius R0 (cells) in a
    background that keeps the mean at zero. Without the mean fix the
    criterion is local in delta_R(d), the top-hat average around a cell at
    distance d from the centre, known analytically from the overlap of two
    spheres: the cell is ionised if zeta f_coll(delta_R(d), R) > 1 on any
    R. The ionised region of excursion() must be a ball of that radius
    (within one cell), larger than R0, and contain no cell beyond it."""
    N, R0, d0 = 64, 5.0, 2.0
    zeta, D, smin = 60.0, 0.3, 2.2
    radii = radii_21cmfast(1.2, 20.0)
    sR = 1.6 * (radii / radii[0]) ** -0.45      # a falling sigma(R) < sigma(M_min)
    c = N / 2
    x = np.arange(N) + 0.5 - c
    r = np.sqrt(x[:, None, None] ** 2 + x[None, :, None] ** 2 + x[None, None, :] ** 2)
    inside = r < R0
    db = -d0 * inside.sum() / (N ** 3 - inside.sum())
    delta = np.where(inside, d0, db).astype(np.float32)
    v0 = inside.sum()          # the grid sphere's volume, for the analytic overlap
    R0e = (3 * v0 / (4 * math.pi)) ** (1 / 3)
    gmax, g0, _ = excursion(delta, float(N), radii,
                            lambda dR, i: fcoll_eps(dR, D, smin, sR[i]), None)
    xh = xhii_of_zeta(gmax, g0, zeta)
    ion = xh >= 1

    def ionised_at(d):
        for i, R in enumerate(radii):
            VR = 4 / 3 * math.pi * R ** 3
            ov = _overlap(R0e, R, d)
            dR = (d0 * ov + db * (VR - ov)) / VR
            if zeta * float(fcoll_eps(np.array([dR]), D, smin, sR[i])[0]) > 1:
                return True
        return False

    ds = np.linspace(0, N / 2, 2001)
    flags = np.array([ionised_at(d) for d in ds])
    r_pred = float(ds[np.nonzero(flags)[0].max()]) if flags.any() else 0.0
    contiguous = flags.any() and flags[:np.nonzero(flags)[0].max() + 1].all()
    r_grid = (3 * ion.sum() / (4 * math.pi)) ** (1 / 3)
    r_far = float(r[ion].max()) if ion.any() else 0.0
    ok = (contiguous and r_pred > R0 + 2 and abs(r_grid - r_pred) < 1.0
          and r_far < r_pred + 1.0 and not ion[r > r_pred + 1.0].any())
    print(f"sphere delta0 {d0} R0 {R0} cells: criterion radius {r_pred:.2f} (contiguous {contiguous}), "
          f"grid ionised ball radius {r_grid:.2f}, farthest ionised cell {r_far:.2f}: "
          f"{'ok' if ok else 'FAILED'}")
    # zeta bisection reproduces a target
    zt = zeta_for_xv(gmax, g0, float(xh.mean()))
    ok2 = abs(xhii_of_zeta(gmax, g0, zt).mean() - xh.mean()) < 1e-4
    print(f"zeta for the same x_V: {zt:.4g} (used {zeta}): {'ok' if ok2 else 'FAILED'}")
    ok = ok and ok2
    try:
        import py21cmfast  # noqa: F401
    except ImportError:
        print("py21cmfast not importable: the 21cmFAST comparison is skipped")
        return ok
    # 21cmFAST itself on a sphere: physical units, mean fix on
    Nn, boxm, zz = 64, 32.0, 8.0
    dx = boxm / Nn
    rs = 3.0 / dx
    delta = np.where(r < rs, 3.0, 0.0)
    delta = (delta - delta.mean()).astype(np.float32)
    cosmo = dict(COSMO_DEFAULT)
    mmin = 1e9
    inp0 = p21_inputs(Nn, boxm, 30.0, 12.0, L_FACTOR * dx, math.log10(mmin), 1, cosmo)
    from py21cmfast.wrapper import cfuncs
    cfuncs.get_growth_factor(inputs=inp0, redshift=zz)
    rho_m = cfuncs.get_condition_mass(inp0, 1.0) / (4.0 / 3.0 * np.pi)
    rmin = max(L_FACTOR * dx, 1.01 * (mmin / (4 / 3 * np.pi * rho_m)) ** (1 / 3))
    rads = radii_21cmfast(rmin, 12.0)
    zeta3 = 15.0
    inp = p21_inputs(Nn, boxm, zeta3, 12.0, rmin, math.log10(mmin), 1, cosmo)
    Dg, sm, sr, _ = p21_sigmas(inp, zz, mmin, rads)
    ib = p21_ionize(delta, zz, inp)
    unc = float(ib.nion_unconditional_acg) / zeta3
    gm, g0_, _ = excursion(delta, boxm, rads, lambda dR, i: fcoll_eps(dR, Dg, sm, sr[i]), unc)
    xa = 1 - np.asarray(ib.neutral_fraction, np.float64)
    xb = xhii_of_zeta(gm, g0_, zeta3)
    mism = float(np.mean((xa > 0.999) != (xb > 0.999)))
    ra = (3 * (xa > 0.999).sum() / (4 * math.pi)) ** (1 / 3) * dx
    ok3 = mism < 2e-3 and (xa > 0.999).sum() > (r < rs).sum()
    print(f"21cmFAST on a sphere (R = 3 Mpc, delta 3, z = {zz}, zeta {zeta3}): ionised ball "
          f"{ra:.2f} Mpc, the numpy criterion differs in {mism:.1e} of the cells: "
          f"{'ok' if ok3 else 'FAILED'}")
    return ok and ok3


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "--selftest":
        try:
            import scipy  # noqa: F401
        except ImportError:
            print("selftest skipped: no scipy")
            sys.exit(77)
        ok = selftest()
        print("selftest", "passed" if ok else "FAILED")
        sys.exit(0 if ok else 1)
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("run")
    r.add_argument("--steps", required=True)
    r.add_argument("--out", required=True)
    r.add_argument("--rmax", type=float, default=15.0, help="R_BUBBLE_MAX [Mpc] (21cmFAST default 15)")
    r.add_argument("--rmax-mfp", default="", help="igm_stats.py JSON of the run: R_max = lambda_mfp")
    r.add_argument("--mmin", type=float, default=1.28e8, help="M_min [Msun] (20 particles of 6.4e6)")
    r.add_argument("--tol", type=float, default=0.002, help="|x_V - target|")
    r.add_argument("--threads", type=int, default=0)
    m = sub.add_parser("run-mass", help="B1 with 21cmFAST's mass-dependent sources (E-INTEGRAL)")
    m.add_argument("run")
    m.add_argument("--steps", required=True)
    m.add_argument("--out", required=True)
    m.add_argument("--name", default="rmax15_mass")
    m.add_argument("--rmax", type=float, default=15.0)
    m.add_argument("--rmin", type=float, default=0.095, help="R_BUBBLE_MIN [Mpc] (B1's)")
    m.add_argument("--mmin", type=float, default=1.28e8, help="M_TURN_STELLAR_FEEDBACK [Msun]")
    m.add_argument("--alpha-star", type=float, default=None,
                   help="default: a(z) - 1 from the run's abundance matching")
    m.add_argument("--pop2-start", type=float, default=5000.0)
    m.add_argument("--maxcall", type=int, default=12)
    m.add_argument("--tol", type=float, default=0.002)
    m.add_argument("--threads", type=int, default=0)
    b = sub.add_parser("standalone")
    b.add_argument("run", help="the reference run (cosmology, box, the x_V targets)")
    b.add_argument("--steps", required=True)
    b.add_argument("--out", required=True)
    b.add_argument("--hii-dim", type=int, default=512)
    b.add_argument("--dim", type=int, default=1024)
    b.add_argument("--seed", type=int, default=12345)
    b.add_argument("--rmax", type=float, default=15.0)
    b.add_argument("--mmin", type=float, default=1.28e8)
    b.add_argument("--zeta0", type=float, default=30.0)
    b.add_argument("--tol", type=float, default=0.002)
    b.add_argument("--threads", type=int, default=0)
    f = sub.add_parser("figure")
    f.add_argument("esdir")
    f.add_argument("runs", nargs="+")
    f.add_argument("--outdir", required=True)
    f.add_argument("--steps", default="")
    f.add_argument("--tol", type=float, default=0.06)
    f.add_argument("--nbins", type=int, default=24)
    f.add_argument("--slab", type=int, default=8)
    f.add_argument("--dpi", type=int, default=200)
    f.add_argument("--es", default="", help="comma list of the variants to show (es_<name>: rmax15, rmaxmfp; "
                   "standalone for B2); default all")
    a = ap.parse_args()
    if a.cmd == "run":
        cmd_run(a)
    elif a.cmd == "run-mass":
        cmd_run_mass(a)
    elif a.cmd == "standalone":
        cmd_standalone(a)
    else:
        cmd_figure(a)


if __name__ == "__main__":
    main()
