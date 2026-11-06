#!/usr/bin/env python3
"""One self-shielding absorber in a background of sources, with every RT
method (Phase 2 (c) of docs/absorbers_closures.md; the cheap form of Wu,
McQuinn & Eisenstein 2021 section 4): an r^-2 overdensity at the centre of a
periodic box, illuminated by a 4^3 lattice of equal sources (the 8 nearest
at equal distance from the absorber), run with apps/rt_driver (RT only, the
gas static) to a steady state, at several resolutions of the same physical
absorber.

usage: absorber_single.py run --exe build/rt_driver --outdir DIR
                              [--N 64,128] [--methods rays,m1,otvet,cone,split,screen,rays_p110]
                              [--box 5] [--z 6] [--gamma 5e-13] [--delta-c 1000]
                              [--r-c 5] [--delta-min 10] [--nsteps 40] [--da 1e-4]
                              [--c-reduced 0.3] [--threads 16] [--analyse-only]
       absorber_single.py make PROBDIR --N 64 [--box 5] ... [--luv L] [--no-absorber]
       absorber_single.py selftest

The absorber: delta(r) = delta_c for r < r_c, delta_c (r_c/r)^2 outside,
cut where it reaches delta_min (r_t = r_c sqrt(delta_c/delta_min)); r_c in
proper kpc. Its gas starts neutral at 10^4 K; the rest of the box is the
mean density, ionised (x_HI = 1e-4, HeII), 10^4 K. The sources' Luv is
set so that the rays give Gamma_HI = --gamma at the absorber's position in
the same box without the absorber (a calibration run, 5 steps). Every
run uses legacy_rt = legacy_physics = 0, so photon_consume_dense = 1
(P11); the method rays_p110 is the rays with photon_consume_dense = 0, the
rule of the 512^3 suite.

Problem files (docs/validation_spec.md): gas.bin (9, N, N, N) float32
Fortran order, u_gas(1:5) then w_gas(1:4) in code units at z_init
(kernels/rt_cell.hpp gas_state_from_avg: u1 = rho/rho_unit, u5 = e_th/
eden_unit with e_th = n_eff k T/(gamma - 1), w1 = n_HI/H_unit = x_HI X u1,
w2 = x_HeI Y u1, w3 = x_HeII Y u1, w4 = 0), sources.txt (i j k Luv sfr,
1-based), params.txt.

Analysis (the last step and the one before, for the steadiness): spherical
shells around the centre: the mean x_HI, Gamma_HI, T and delta; the
case-B recombination rate and the photoionisations inside r_t (the code's
alpha_B lookup of absorber_budget.py), the shell radius where x_HI
reaches 0.5, the HI column through the centre along x. Writes
DIR/single.txt (tables), DIR/single.json and DIR/absorber_single.png.
"""
import argparse
import json
import math
import os
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from absorber_budget import MH_CGS, YHE, alpha_b  # noqa: E402
from igm_stats import (HUBBLE0, MPC_CM, OMEGAB, OMEGAM, OMEGAR, RHOC_CGS, XHY,  # noqa: E402
                       H0_CGS, K_CGS, SIGMA_912)

GAMMA_IDEAL = 5.0 / 3.0
METHOD_KEYS = {
    "rays": [],
    "rays_p110": ["photon_consume_dense=0"],
    "m1": ["rt_method=m1"],
    "otvet": ["rt_method=otvet"],
    "cone": ["rt_method=otvet", "otvet_closure=cone"],
    "split": ["rt_method=otvet", "otvet_split=1"],
    "screen": ["rt_method=otvet", "otvet_closure=screen"],
}


class Units:
    """cosmology.cpp calc_sim_var at scale factor a for an N^3 box of `box`
    Mpc/h (the constants of constants.hpp)."""

    def __init__(self, N, box, a):
        h = HUBBLE0
        self.a, self.N, self.box = a, N, box
        self.x_unit = a * box / h / N * MPC_CM
        self.rho_unit = RHOC_CGS * OMEGAM * h * h / a ** 3
        self.time_unit = 2 / (3 * H0_CGS * h * math.sqrt(OMEGAM)) * a * a
        self.vel_unit = self.x_unit / self.time_unit
        self.eden_unit = self.rho_unit * self.vel_unit ** 2
        self.rhob_avg = self.rho_unit * OMEGAB / OMEGAM
        self.nH_avg = self.rhob_avg * XHY / MH_CGS
        self.nHe_avg = self.rhob_avg * YHE / (4 * MH_CGS)
        self.H = H0_CGS * h * math.sqrt(OMEGAM / a ** 3 + (1 - OMEGAM) + OMEGAR / a ** 4)
        self.dl_pkpc = self.x_unit / MPC_CM * 1e3


def gas_record(delta, xHI, xHeI, xHeII, T, u):
    """u_gas(1:5), w_gas(1:4) of cells with these properties (arrays)."""
    u1 = delta * (OMEGAB / OMEGAM)
    nH = delta * u.nH_avg
    nHe = delta * u.nHe_avg
    ne = (1 - xHI) * nH + nHe * (xHeII + 2 * (1 - xHeI - xHeII))
    neff = nH + nHe + ne
    eth = neff * K_CGS * T / (GAMMA_IDEAL - 1)
    u5 = eth / u.eden_unit
    w1 = xHI * XHY * u1
    w2 = xHeI * YHE * u1
    w3 = xHeII * YHE * u1
    z = np.zeros_like(u1)
    return np.stack([u1, z, z, z, u5, w1, w2, w3, z], axis=-1)


def absorber_profile(r_pkpc, delta_c, r_c, delta_min):
    d = np.where(r_pkpc < r_c, delta_c, delta_c * (r_c / np.maximum(r_pkpc, 1e-30)) ** 2)
    return np.where(d >= delta_min, d, 1.0)


def make_problem(dirpath, N, box, z, delta_c, r_c, delta_min, luv, absorber=True,
                 nsteps=40, da=1e-4, c_reduced=0.3, nsrc_side=4):
    a = 1.0 / (1 + z)
    u = Units(N, box, a)
    os.makedirs(dirpath, exist_ok=True)
    # cell centres relative to the box centre, in proper kpc (periodic distance)
    # the centre at cell N/2 (a lattice-cell centre of the sources below)
    idx = np.arange(N)
    d1 = idx - N / 2.0
    d1 = np.minimum(np.abs(d1), N - np.abs(d1)) * u.dl_pkpc
    r = np.sqrt(d1[:, None, None] ** 2 + d1[None, :, None] ** 2 + d1[None, None, :] ** 2)
    if absorber:
        delta = absorber_profile(r, delta_c, r_c, delta_min)
    else:
        delta = np.ones((N, N, N))
    inside = delta > 1.0
    xHI = np.where(inside, 1.0, 1e-4)
    xHeI = np.where(inside, 1.0, 0.0)
    xHeII = np.where(inside, 0.0, 1.0)
    T = np.full((N, N, N), 1e4)
    rec = gas_record(delta, xHI, xHeI, xHeII, T, u).astype(np.float32)   # [k, j, i, 9]
    # (9, N, N, N) Fortran order = component fastest, then i, j, k: the
    # [k, j, i, comp] C array written contiguously
    rec.tofile(os.path.join(dirpath, "gas.bin"))
    with open(os.path.join(dirpath, "sources.txt"), "w") as f:
        step = N // nsrc_side
        for k in range(nsrc_side):
            for j in range(nsrc_side):
                for i in range(nsrc_side):
                    # 1-based cell indices at N/8 + m N/4 (the absorber at N/2 is a lattice-cell centre)
                    f.write(f"{step // 2 + i * step + 1} {step // 2 + j * step + 1} "
                            f"{step // 2 + k * step + 1} {luv:.6e} {luv / 2e28:.6e}\n")
    with open(os.path.join(dirpath, "params.txt"), "w") as f:
        f.write(f"N = {N}\nbox = {box}\nz_init = {z}\nda_sim = {da}\nnsteps = {nsteps}\n"
                f"c_reduced = {c_reduced}\nlegacy_rt = 0\nlegacy_physics = 0\n"
                f"deterministic_parallel = 1\nnray_max_factor = 15\n")
    with open(os.path.join(dirpath, "README.txt"), "w") as f:
        f.write(f"absorber_single.py: {'an r^-2 absorber' if absorber else 'no absorber'}, "
                f"delta_c {delta_c}, r_c {r_c} pkpc, cut at delta {delta_min}; {nsrc_side ** 3} "
                f"sources of Luv {luv:.3e}; cell {u.dl_pkpc:.2f} pkpc at z = {z}\n")
    return u, r, delta


def read_step(outdir, step, N):
    d = os.path.join(outdir, f"step_{step:04d}")
    fHI = np.fromfile(os.path.join(d, "fHI.bin"), np.float32).reshape(N, N, N).astype(np.float64)
    temp = np.fromfile(os.path.join(d, "temp.bin"), np.float32).reshape(N, N, N).astype(np.float64)
    rhog = np.fromfile(os.path.join(d, "rhog.bin"), np.float32).reshape(N, N, N).astype(np.float64)
    fHeI = np.fromfile(os.path.join(d, "fHeI.bin"), np.float32).reshape(N, N, N).astype(np.float64)
    fHeII = np.fromfile(os.path.join(d, "fHeII.bin"), np.float32).reshape(N, N, N).astype(np.float64)
    gam = np.fromfile(os.path.join(d, "gamma.bin"), np.float32).reshape(N, N, N, 3)[..., 0].astype(np.float64)
    return fHI, fHeI, fHeII, temp, rhog, gam


def last_step(outdir):
    steps = sorted(int(x[5:]) for x in os.listdir(outdir) if x.startswith("step_"))
    return steps


def analyse_run(outdir, N, u, r, r_t, step=None):
    steps = last_step(outdir)
    if step is None:
        step = steps[-1]
    fHI, fHeI, fHeII, T, rhog, gam = read_step(outdir, step, N)
    V = u.x_unit ** 3
    nH = rhog * u.nH_avg
    nHI = fHI * nH
    nHII = nH - nHI
    ne = nHII + rhog * u.nHe_avg * (fHeII + 2 * (1 - fHeI - fHeII))
    R = alpha_b(T) * ne * nHII * V
    P = gam * nHI * V
    out = {"step": step}
    # shells of one cell width out to 3 r_t
    edges = np.arange(0, 3 * r_t + u.dl_pkpc, u.dl_pkpc)
    ib = np.digitize(r.ravel(), edges) - 1
    prof = {"r": (0.5 * (edges[1:] + edges[:-1])).tolist()}
    for name, arr in (("xHI", fHI), ("gamma", gam), ("T", T), ("delta", rhog)):
        s = np.bincount(ib[ib >= 0], weights=arr.ravel()[ib >= 0], minlength=len(edges) - 1)
        n = np.bincount(ib[ib >= 0], minlength=len(edges) - 1)
        prof[name] = (s / np.maximum(n, 1))[:len(edges) - 1].tolist()
    out["profile"] = prof
    m = r <= r_t
    out["R_in"] = float(R[m].sum())
    out["P_in"] = float(P[m].sum())
    out["R_box"] = float(R.sum())
    out["P_box"] = float(P.sum())
    out["n_in"] = int(m.sum())
    out["xHI_centre"] = float(fHI[N // 2, N // 2, N // 2])
    out["gamma_centre"] = float(gam[N // 2, N // 2, N // 2])
    out["N_HI_x"] = float((nHI[N // 2, N // 2, :] * u.x_unit).sum())
    out["mass_HI_in"] = float(nHI[m].sum() * V)
    # the radius where the shell-mean x_HI reaches 0.5, from the outside in
    x = np.array(prof["xHI"])
    rr = np.array(prof["r"])
    out["r_half"] = float("nan")
    for i in range(len(x) - 1, 0, -1):
        if x[i] < 0.5 <= x[i - 1]:
            f = (0.5 - x[i]) / (x[i - 1] - x[i])
            out["r_half"] = float(rr[i] - f * (rr[i] - rr[i - 1]))
            break
    if np.isnan(out["r_half"]) and x[0] >= 0.5:
        out["r_half"] = float(rr[0])
    # the mean Gamma in the shell just outside the absorber
    m2 = (r > r_t) & (r <= r_t + 2 * u.dl_pkpc)
    out["gamma_rim"] = float(gam[m2].mean()) if m2.any() else float("nan")
    return out


def run_driver(exe, prob, out, keys, threads):
    cmd = [exe, prob, out] + sum((["--set", k] for k in keys), [])
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    with open(out + ".log", "w") as log:
        rc = subprocess.call(cmd, stdout=log, stderr=subprocess.STDOUT, env=env)
    if rc != 0:
        raise SystemExit(f"rt_driver failed ({rc}): {' '.join(cmd)}; see {out}.log")


def cmd_run(a):
    os.makedirs(a.outdir, exist_ok=True)
    Ns = [int(x) for x in a.N.split(",")]
    methods = a.methods.split(",")
    results = {}
    r_t = a.r_c * math.sqrt(a.delta_c / a.delta_min)
    for N in Ns:
        key = f"N{N}"
        results[key] = {}
        u = Units(N, a.box, 1.0 / (1 + a.z))
        base = os.path.join(a.outdir, key)
        cal_prob, cal_out = os.path.join(base, "cal_prob"), os.path.join(base, "cal_rays")
        luv0 = 1e28
        if not a.analyse_only:
            make_problem(cal_prob, N, a.box, a.z, a.delta_c, a.r_c, a.delta_min, luv0,
                         absorber=False, nsteps=a.cal_steps, da=a.da, c_reduced=a.c_reduced)
            run_driver(a.exe, cal_prob, cal_out, [], a.threads)
        # the calibration: Gamma at the centre of the empty box after cal_steps
        st = last_step(cal_out)[-1]
        _, _, _, Tc, rc_, gc = read_step(cal_out, st, N)
        g_centre = float(gc[N // 2, N // 2, N // 2])
        luv = luv0 * a.gamma / g_centre
        results[key]["calibration"] = {"luv0": luv0, "gamma_centre": g_centre, "luv": luv,
                                       "step": st, "rhog_centre": float(rc_[N // 2, N // 2, N // 2]),
                                       "T_centre": float(Tc[N // 2, N // 2, N // 2]),
                                       "dl_pkpc": u.dl_pkpc,
                                       "dt_myr": a.da / (u.a * u.H) / (365.25 * 86400 * 1e6)}
        print(f"{key}: cell {u.dl_pkpc:.2f} pkpc, Gamma(centre) {g_centre:.3e} at Luv {luv0:.1e} "
              f"-> Luv {luv:.3e} for {a.gamma:.1e}", flush=True)
        prob = os.path.join(base, "prob")
        u, r, delta = make_problem(prob, N, a.box, a.z, a.delta_c, a.r_c, a.delta_min, luv,
                                   absorber=True, nsteps=a.nsteps, da=a.da, c_reduced=a.c_reduced)
        results[key]["r_t"] = r_t
        results[key]["n_cells_in"] = int((r <= r_t).sum())
        results[key]["delta_max"] = float(delta.max())
        for mth in methods:
            out = os.path.join(base, mth)
            if not a.analyse_only:
                run_driver(a.exe, prob, out, METHOD_KEYS[mth], a.threads)
            steps = last_step(out)
            res = analyse_run(out, N, u, r, r_t)
            prev = analyse_run(out, N, u, r, r_t, step=steps[-2]) if len(steps) > 1 else None
            res["R_in_prev"] = prev["R_in"] if prev else float("nan")
            results[key][mth] = res
            print(f"  {mth}: R_in {res['R_in']:.3e} (prev {res['R_in_prev']:.3e}), x_HI(centre) "
                  f"{res['xHI_centre']:.3e}, r_half {res['r_half']:.1f} pkpc, N_HI {res['N_HI_x']:.2e}",
                  flush=True)
    with open(os.path.join(a.outdir, "single.json"), "w") as f:
        json.dump({"args": vars(a), "results": results}, f, indent=1)
    write_report(results, methods, a)
    make_figure(results, methods, a)


def write_report(results, methods, a):
    L = [f"# absorber_single.py: one r^-2 absorber (delta_c {a.delta_c}, r_c {a.r_c} pkpc, cut at "
         f"delta {a.delta_min}, r_t {a.r_c * math.sqrt(a.delta_c / a.delta_min):.0f} pkpc) in a box of "
         f"{a.box} Mpc/h at z = {a.z}, 64 sources set for Gamma = {a.gamma:.1e} at the centre; "
         f"c_reduced {a.c_reduced}, {a.nsteps} steps of da {a.da}\n"]
    for key, R in results.items():
        c = R["calibration"]
        L.append(f"## {key}: cell {c['dl_pkpc']:.2f} pkpc, step {c['dt_myr']:.2f} Myr, "
                 f"{R['n_cells_in']} cells inside r_t, delta_max {R['delta_max']:.0f}, "
                 f"Luv {c['luv']:.3e} (Gamma {c['gamma_centre']:.3e} at {c['luv0']:.0e})\n")
        L.append("| method | R(<r_t) | / rays | steady (R_prev/R) | P(<r_t) | x_HI(centre) | "
                 "r_half [pkpc] | N_HI(x) [cm^-2] | Gamma(centre) | Gamma(rim) | M_HI(<r_t)/rays |")
        L.append("|---|---|---|---|---|---|---|---|---|---|---|")
        ref = R.get("rays")
        for m in methods:
            r = R.get(m)
            if not r:
                continue
            L.append(f"| {m} | {r['R_in']:.3e} | {r['R_in'] / ref['R_in'] if ref else float('nan'):.3f} | "
                     f"{r['R_in_prev'] / r['R_in']:.3f} | {r['P_in']:.3e} | {r['xHI_centre']:.3e} | "
                     f"{r['r_half']:.1f} | {r['N_HI_x']:.2e} | {r['gamma_centre']:.2e} | "
                     f"{r['gamma_rim']:.2e} | {r['mass_HI_in'] / ref['mass_HI_in'] if ref else float('nan'):.3f} |")
        L.append("")
    with open(os.path.join(a.outdir, "single.txt"), "w") as f:
        f.write("\n".join(L) + "\n")


def make_figure(results, methods, a):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    colors = {"rays": "black", "rays_p110": "0.5", "m1": "#2d7bd6", "otvet": "#d9541e",
              "cone": "#2a9d6b", "split": "#8e5bc2", "screen": "#c99a1b"}
    keys = list(results)
    fig, axes = plt.subplots(3, len(keys), figsize=(4.2 * len(keys), 9), squeeze=False)
    for j, key in enumerate(keys):
        R = results[key]
        for m in methods:
            r = R.get(m)
            if not r:
                continue
            p = r["profile"]
            rr = np.array(p["r"])
            axes[0, j].plot(rr, p["xHI"], color=colors.get(m, None), lw=1.5, label=m)
            axes[1, j].plot(rr, p["gamma"], color=colors.get(m, None), lw=1.5)
            axes[2, j].plot(rr, p["T"], color=colors.get(m, None), lw=1.5)
        axes[0, j].set_title(f"{key}: cell {R['calibration']['dl_pkpc']:.1f} pkpc")
        axes[0, j].set_yscale("log")
        axes[0, j].set_ylim(1e-6, 1.5)
        axes[1, j].set_yscale("log")
        axes[2, j].set_yscale("log")
        axes[2, j].set_xlabel("r [pkpc]")
        for ax in axes[:, j]:
            ax.axvline(R["r_t"], color="0.7", lw=0.8)
    axes[0, 0].set_ylabel("shell-mean $x_{\\rm HI}$")
    axes[1, 0].set_ylabel("shell-mean $\\Gamma_{\\rm HI}$ [s$^{-1}$]")
    axes[2, 0].set_ylabel("shell-mean T [K]")
    axes[0, 0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(a.outdir, "absorber_single.png"), dpi=130)


def selftest():
    """The gas record round trip (gas_state_from_avg in reverse) and the
    profile: no driver needed."""
    u = Units(32, 5.0, 1 / 7.0)
    delta = np.array([1.0, 100.0])
    rec = gas_record(delta, np.array([1e-4, 1.0]), np.array([0.0, 1.0]), np.array([1.0, 0.0]),
                     np.array([1e4, 1e4]), u)
    errs = []
    # rho back: u1 rho_unit / rhob_avg = delta
    got = rec[:, 0] * u.rho_unit / u.rhob_avg
    if not np.allclose(got, delta, rtol=1e-12):
        errs.append(f"delta round trip {got}")
    # T back: (gamma-1) e_th / (k n_eff)
    nH = delta * u.nH_avg
    nHe = delta * u.nHe_avg
    xHI = np.array([1e-4, 1.0])
    ne = (1 - xHI) * nH + nHe * np.array([1.0, 0.0])
    neff = nH + nHe + ne
    T = (GAMMA_IDEAL - 1) * rec[:, 4] * u.eden_unit / (K_CGS * neff)
    if not np.allclose(T, 1e4, rtol=1e-10):
        errs.append(f"T round trip {T}")
    # n_HI back: w1 H_unit = w1 rho_unit/mH
    nHI = rec[:, 5] * u.rho_unit / MH_CGS
    if not np.allclose(nHI, xHI * nH, rtol=1e-10):
        errs.append(f"nHI round trip {nHI} vs {xHI * nH}")
    # the profile
    r = np.array([0.0, 5.0, 10.0, 50.0, 60.0])
    d = absorber_profile(r, 1000.0, 5.0, 10.0)
    if not np.allclose(d, [1000, 1000, 250, 10, 1]):
        errs.append(f"profile {d}")
    # a small problem is written and read back
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        uu, rr, dd = make_problem(tmp, 16, 5.0, 6.0, 1000.0, 5.0, 10.0, 1e28)
        g = np.fromfile(os.path.join(tmp, "gas.bin"), np.float32).reshape(16, 16, 16, 9)
        if abs(float(g[8, 8, 8, 0]) - 1000 * OMEGAB / OMEGAM) > 1e-3:
            errs.append(f"centre density {g[8, 8, 8, 0]}")
        if abs(float(g[0, 0, 0, 0]) - OMEGAB / OMEGAM) > 1e-6:
            errs.append(f"corner density {g[0, 0, 0, 0]}")
        src = open(os.path.join(tmp, "sources.txt")).read().split("\n")
        if len([s for s in src if s.strip()]) != 64:
            errs.append("64 sources expected")
    if errs:
        print("absorber_single selftest FAILED:\n  " + "\n  ".join(errs))
        return 1
    print("absorber_single selftest OK")
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd")
    r = sub.add_parser("run")
    r.add_argument("--exe", required=True)
    r.add_argument("--outdir", required=True)
    r.add_argument("--N", default="64,128")
    r.add_argument("--methods", default="rays,m1,otvet,cone,split,screen,rays_p110")
    r.add_argument("--box", type=float, default=5.0)
    r.add_argument("--z", type=float, default=6.0)
    r.add_argument("--gamma", type=float, default=5e-13)
    r.add_argument("--delta-c", type=float, default=1000.0)
    r.add_argument("--r-c", type=float, default=5.0, help="core radius [proper kpc]")
    r.add_argument("--delta-min", type=float, default=10.0)
    r.add_argument("--nsteps", type=int, default=40)
    r.add_argument("--cal-steps", type=int, default=5)
    r.add_argument("--da", type=float, default=1e-4)
    r.add_argument("--c-reduced", type=float, default=0.3)
    r.add_argument("--threads", type=int, default=16)
    r.add_argument("--analyse-only", action="store_true")
    m = sub.add_parser("make")
    m.add_argument("probdir")
    m.add_argument("--N", type=int, default=64)
    m.add_argument("--box", type=float, default=5.0)
    m.add_argument("--z", type=float, default=6.0)
    m.add_argument("--delta-c", type=float, default=1000.0)
    m.add_argument("--r-c", type=float, default=5.0)
    m.add_argument("--delta-min", type=float, default=10.0)
    m.add_argument("--luv", type=float, default=1e28)
    m.add_argument("--no-absorber", action="store_true")
    m.add_argument("--nsteps", type=int, default=40)
    m.add_argument("--da", type=float, default=1e-4)
    m.add_argument("--c-reduced", type=float, default=0.3)
    sub.add_parser("selftest")
    a = ap.parse_args()
    if a.cmd == "selftest":
        return selftest()
    if a.cmd == "make":
        make_problem(a.probdir, a.N, a.box, a.z, a.delta_c, a.r_c, a.delta_min, a.luv,
                     absorber=not a.no_absorber, nsteps=a.nsteps, da=a.da, c_reduced=a.c_reduced)
        return 0
    if a.cmd == "run":
        cmd_run(a)
        return 0
    ap.print_help()
    return 2


if __name__ == "__main__":
    sys.exit(main())
