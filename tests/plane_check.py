#!/usr/bin/env python3
"""Physics checks of the plane-parallel RT-domain mode (docs/plane_parallel.md)
on the RT-only driver.

usage: plane_check.py RT_DRIVER UNIFORM_PROBLEM_DIR WORKDIR [--only NAME ...]
                      [--set key=value ...] [--threads n]

Every check builds its own problem from the uniform gas row of
UNIFORM_PROBLEM_DIR (tests/problems/uniform_1src: mean density, neutral, at
temp_min), scaled in density, ionisation and temperature, and runs rt_driver
with plane_parallel = 1 (legacy_rt = legacy_physics = 0). --set is passed
on to every run (e.g. --set rt_backend=amrex with --launch for the backend).

Checks (each exits 1 on failure):
  gamma     Gamma_HI normalisation: in an ionised, optically thin box the
            HI photoionisation rate of the source cells (and of every lit
            cell) is plane_gamma12 * 1e-12 s^-1 to 1e-3, for 1, 2 and 4 RT
            cells per substep; the bins of plane_spectrum.txt are the
            equal-photon bins of J_nu ~ nu^-1.5 between 1 and 4 Ry; the
            optical depths of the grid are n_HI sigma_HI(nu_n) len.
  slab      Attenuation: photons through a neutral slab of 16 cells with
            tau_cell(bin 1) = 0.05 are transmitted as prod_j (1 - tau_nj)
            per bin (the first-order-per-cell absorption of photo_cell, to
            2e-3) and as exp(-tau_n) to 6 per cent, harder bins less
            attenuated; the gas outside the slab is transparent.
  front     One plane-parallel front crossing the whole box (plane_ndom =
            1, one face): in neutral uniform gas at the mean density the
            ionised column x_I(t) (cells) follows the photon budget
            dx/dt = [F - alpha_B(T) n_e n_HII x] / (n_H + n_He) to 5 per
            cent once x_I > 3 cells (c_reduced = 1, the front at ~0.03 c);
            (at 100x the mean density with t_rec ~ 24 Myr over 2.3 t_rec: the
            H column equals the code's photoionisations + collisional
            ionisations - recombinations, and the effective recombination
            coefficient is case B at the interior temperature to a factor 1.5);
            with c_reduced = 0.1 and 100x the flux the front is limited by
            the reduced speed of light and moves at c_hat (to 10 per cent).
  conserve  Photon conservation (the front and slab runs): per hydro step
            emitted = absorbed + escaped + (in flight at the end - in
            flight at the start), to 1e-3 of the emitted photons.
  domains   plane_ndom = 2, faces +x +y: the 8 domains ionise identically
            (their x_HII,V agree to 1e-5 at every step); no ray crosses a
            domain edge: in the thin ionised box the source cells of the
            second domain along each axis get exactly one ray's photons
            (Gamma_HI = the target, not twice it) and no cell exceeds it.
  threads   deterministic_parallel = 1: the front run with 1 and 2 threads
            is byte-identical (fields and stats, the plane counters included).
  sweep     plane_transport = sweep against rays: the thin front (one face)
            and the 2^3-domain problem (two faces) with 1 and 2 cells per
            substep give the same fHI, temp and photon fields (to 3e-5,
            float round-off: one packet per cell and direction either
            way) and the same photon bookkeeping to 2e-4.
"""
import argparse
import os
import shutil
import subprocess
import sys

import numpy as np

# physical_constants.f90 / cosmo_parameters.f90 (include/radray/constants.hpp)
C_CGS = 2.99792e10
K_CGS = 1.38066e-16
MH = 1.67223e-24
MHE = 4 * MH
H0_CGS = 3.24086e-18
RHOC = 1.87890e-29
MPC = 3.08560e24
HI_EV = float(np.float32(13.60))
OM, OB, H0, XHY, YHE = float(np.float32(0.27)), float(np.float32(0.045)), float(np.float32(0.70)), \
    float(np.float32(0.76)), float(np.float32(0.24))
GAM = float(np.float32(5.0) / np.float32(3.0))
YR = 365 * 86400.0

EXTRA = []
LAUNCH = []
THREADS = 1


def units(z, box, N):
    """calc_sim_var (src/cosmology.cpp) for Ngas = Nrt = N."""
    a = 1 / (1 + z)
    x_unit = a * box / H0 / N * MPC
    rho_unit = RHOC * OM * H0 ** 2 / a ** 3
    time_unit = 2 / (3 * H0_CGS * H0 * np.sqrt(OM)) * a ** 2
    vel_unit = x_unit / time_unit
    eden_unit = rho_unit * vel_unit ** 2
    return dict(a=a, x_unit=x_unit, len_rt=x_unit, rho_unit=rho_unit, eden_unit=eden_unit,
                nH_avg=rho_unit * OB / OM * XHY / MH, nHe_avg=rho_unit * OB / OM * YHE / MHE,
                H_unit=rho_unit / MH, He_unit=rho_unit / MHE)


def gas_row(u, delta, fHI, fHeI, fHeII, T):
    """One gas.bin row (u(5), w(4)): overdensity delta, ionisation fractions,
    temperature T (the energy from the number of particles the code counts)."""
    rho = delta * OB / OM
    nH = rho * u["rho_unit"] * XHY / MH
    nHe = rho * u["rho_unit"] * YHE / MHE
    ne = nH * (1 - fHI) + nHe * (fHeII + 2 * (1 - fHeI - fHeII))
    neff = nH + nHe + ne
    eth = neff * K_CGS * T / (GAM - 1) / u["eden_unit"]
    return [rho, 0.0, 0.0, 0.0, eth, rho * XHY * fHI, rho * YHE * fHeI, rho * YHE * fHeII, 0.0]


def da_for_substeps(z, box, N, nsub, c_reduced):
    """da_sim for which calc_timestep gives nsub RT substeps (one cell each)."""
    a = 1 / (1 + z)
    dl_rt = box / H0 * MPC / N
    # dl_c = 2c/(H0 h sqrt(om)) (sqrt(a + da) - sqrt(a)) = nsub dl_rt / c_reduced
    x = nsub * dl_rt / c_reduced * (H0_CGS * H0 * np.sqrt(OM)) / (2 * C_CGS)
    return (np.sqrt(a) + x) ** 2 - a


def verner_HI(nu):
    nu = np.asarray(nu, dtype=np.float64)
    s = lambda v: float(np.float32(v))
    x = nu / s(0.4298)
    sig = 5.475e-14 * (x - 1) ** 2 / x ** s(4.0185) / (1 + np.sqrt(x / s(32.88))) ** s(2.963)
    return np.where(nu < HI_EV, 0.0, sig)


def write_problem(dst, N, box, z, rows, da, nsteps, extra_params):
    """rows: an (N,N,N,9) array of gas rows, [k,j,i]."""
    os.makedirs(dst, exist_ok=True)
    np.asarray(rows, dtype=np.float32).ravel().tofile(os.path.join(dst, "gas.bin"))
    with open(os.path.join(dst, "params.txt"), "w") as f:
        f.write(f"N = {N}\nbox = {float(box)!r}\nz_init = {z}\nda_sim = {float(da)!r}\nnsteps = {nsteps}\n"
                "legacy_rt = 0\nlegacy_physics = 0\nplane_parallel = 1\nnray_max_factor = 16\n")
        for k, v in extra_params.items():
            f.write(f"{k} = {v}\n")
    open(os.path.join(dst, "sources.txt"), "w").close()  # no point sources


def run(driver, prob, out, sets=(), threads=None):
    if os.path.isdir(out):
        shutil.rmtree(out)
    os.makedirs(out)
    cmd = LAUNCH + [driver, prob, out] + EXTRA
    for kv in sets:
        cmd += ["--set", kv]
    cmd += ["--set", f"threads={threads if threads else THREADS}"]
    with open(out + ".stdout", "w") as log:
        r = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
    if r.returncode != 0:
        print(open(out + ".stdout").read()[-4000:])
        sys.exit(f"FAIL: {' '.join(cmd)}")
    return out


def parse(fn):
    d = {}
    for line in open(fn):
        if "=" in line and not line.startswith("#"):
            k, v = line.split("=", 1)
            try:
                d[k.strip()] = [float(x) for x in v.split()]
            except ValueError:
                pass
    return d


def field(sd, name, N, nc=1):
    a = np.fromfile(os.path.join(sd, name + ".bin"), dtype=np.float32).astype(np.float64)
    return a.reshape(N, N, N, nc) if nc > 1 else a.reshape(N, N, N)


FAILS = []


def check(cond, msg):
    print(("  ok   " if cond else "  FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


# ---------------------------------------------------------------------------
def check_gamma(driver, work, row0):
    print("gamma: the Gamma_HI normalisation, the bins, the optical depths")
    N, box, z, g12 = 16, 0.1, 7.0, 0.3
    u = units(z, box, N)
    fHI = 1e-4
    rows = np.tile(gas_row(u, 1.0, fHI, 1e-4, 1.0, 1e4), (N, N, N, 1))
    spec = None
    for k in (1, 2, 4):
        prob = os.path.join(work, f"gamma_k{k}")
        write_problem(prob, N, box, z, rows, da_for_substeps(z, box, N, 40 * k, 1.0), 1,
                      {"plane_ndom": 1, "plane_faces": "+x", "plane_gamma12": g12, "rt_cells_per_substep": k})
        out = run(driver, prob, prob + "_run")
        gam = field(os.path.join(out, "step_0001"), "gamma", N, 3)[..., 0]
        target = g12 * 1e-12
        src = gam[:, :, 0] / target - 1
        allc = gam / target - 1
        check(np.max(np.abs(src)) < 1e-3, f"k = {k}: source cells Gamma_HI/target - 1 = {np.max(np.abs(src)):.2e} < 1e-3")
        # the lit cells: the target minus the attenuation across the box
        # (tau_cell ~ 1e-4 at fHI = 1e-4: ~2e-3 over 16 cells)
        check(np.max(np.abs(allc)) < 5e-3, f"k = {k}: every cell Gamma_HI/target - 1 = {np.max(np.abs(allc)):.2e} < 5e-3")
        if spec is None:
            spec = parse(os.path.join(out, "plane_spectrum.txt"))
            tau = field(os.path.join(out, "step_0001"), "tau", N, 5)[0, 0, 5]
            fHI_f = field(os.path.join(out, "step_0001"), "fHI", N)[0, 0, 5]
            fHeI_f = field(os.path.join(out, "step_0001"), "fHeI", N)[0, 0, 5]
    # the bins: equal photon number of nu^-1.5 between 1 and 4 Ry
    alpha = spec["alpha"][0]
    edges = np.array([spec[f"bin {n}"][0] for n in range(5)] + [spec["bin 4"][1]])
    want = HI_EV * (1 - np.arange(6) / 5 * (1 - 4.0 ** -alpha)) ** (-1 / alpha)
    check(np.allclose(edges, want, rtol=1e-10), f"bin edges {edges.round(3)} eV = the equal-photon edges")
    nu = np.array([spec[f"bin {n}"][2] for n in range(5)])
    frac = np.array([spec[f"bin {n}"][3] for n in range(5)])
    check(np.allclose(frac, 0.2, rtol=1e-12), "photon fractions 1/5")
    inside = (nu > edges[:-1]) & (nu < edges[1:])
    check(bool(np.all(inside)), f"effective frequencies {nu.round(2)} eV inside their bins")
    sig = verner_HI(nu)
    check(np.allclose(sig, [spec[f"bin {n}"][4] for n in range(5)], rtol=1e-6), "sigma_HI(nu_n) as written")
    # the optical depths of the grid: (n_HI sigma_HI(nu_n) + n_HeI sigma_HeI(nu_n)) len
    sig_HeI = np.array([spec[f"bin {n}"][5] for n in range(5)])
    tau_want = (u["nH_avg"] * fHI_f * sig + u["nHe_avg"] * fHeI_f * sig_HeI) * u["len_rt"]
    check(np.allclose(tau, tau_want, rtol=1e-3), f"tau of the groups = (n_HI sigma_HI + n_HeI sigma_HeI) len ({np.max(np.abs(tau / tau_want - 1)):.1e})")
    return spec


def check_slab(driver, work, spec):
    print("slab: attenuation through a neutral slab")
    N, box, z, g12 = 32, 0.005, 7.0, 3e-3  # a small box: Gamma t < 1e-3 over the run
    u = units(z, box, N)
    nu = np.array([spec[f"bin {n}"][2] for n in range(5)])
    sig = verner_HI(nu)
    tau_cell = 0.05
    delta = tau_cell / (u["nH_avg"] * sig[0] * u["len_rt"])
    rows = np.tile(gas_row(u, 1e-3, 1.0, 1.0, 0.0, 100.0), (N, N, N, 1))
    rows[:, :, 8:24, :] = gas_row(u, delta, 1.0, 1.0, 0.0, 100.0)
    prob = os.path.join(work, "slab")
    write_problem(prob, N, box, z, rows, da_for_substeps(z, box, N, 48, 1.0), 1,
                  {"plane_ndom": 1, "plane_faces": "+x", "plane_gamma12": g12})
    out = run(driver, prob, prob + "_run")
    sd = os.path.join(out, "step_0001")
    ngam = field(sd, "ngam", N, 5)
    fgam = field(sd, "fgam", N, 5)
    tau = field(sd, "tau", N, 5)
    dep = np.where(fgam > 0, ngam / np.maximum(fgam, 1e-30), 0.0)  # deposited = transmitted / fgam
    d0 = dep[:, :, 0, :].mean(axis=(0, 1))
    d1 = dep[:, :, 24, :].mean(axis=(0, 1))
    T = d1 / d0
    tau_slab = tau[0, 0, 8:24, :]
    prod = np.prod(1 - np.minimum(tau_slab, 1), axis=0)
    expo = np.exp(-tau_slab.sum(axis=0))
    check(np.max(np.abs(T / prod - 1)) < 2e-3, f"transmission per bin {T.round(4)} = prod(1 - tau_j) {prod.round(4)} to 2e-3 ({np.max(np.abs(T / prod - 1)):.1e})")
    check(np.max(np.abs(T / expo - 1)) < 0.06, f"transmission = exp(-tau) {expo.round(4)} to 6 per cent ({np.max(np.abs(T / expo - 1)):.3f})")
    check(bool(np.all(np.diff(T) > 0)), "harder bins less attenuated")
    before = dep[:, :, 7, :].mean(axis=(0, 1)) / d0
    check(np.max(np.abs(before - 1)) < 1e-3, f"the gas before the slab is transparent ({np.max(np.abs(before - 1)):.1e})")
    check(abs(tau_slab[0, 0] / tau_cell - 1) < 1e-3, f"tau_cell(bin 1) = {tau_slab[0, 0]:.4f} as built")
    return out


def columns(sd, N):
    """Ionised columns in cells along x: H, and He (singly) in H-atom units."""
    fHI = field(sd, "fHI", N)
    fHeI = field(sd, "fHeI", N)
    return np.sum(1 - fHI) / N ** 2, np.sum(1 - fHeI) / N ** 2 * YHE / MHE * MH / XHY


def budget_run(driver, work, name, spec, N, box, z, g12, delta, nsub, nsteps, c_reduced, tol):
    """The ionised column (H + He in H-atom units) against the photon
    budget: (photons that left the rays)/n_H minus the recombinations,
    integrated with alpha_B at the measured temperature of the ionised gas."""
    u = units(z, box, N)
    rows = np.tile(gas_row(u, delta, 1.0, 1.0, 0.0, 100.0), (N, N, N, 1))
    prob = os.path.join(work, name)
    write_problem(prob, N, box, z, rows, da_for_substeps(z, box, N, nsub, c_reduced), nsteps,
                  {"plane_ndom": 1, "plane_faces": "+x", "plane_gamma12": g12, "c_reduced": c_reduced,
                   "photon_consume_dense": 1})  # every ionisation consumes its photon (the dense box)
    out = run(driver, prob, prob + "_run")
    F = parse(os.path.join(out, "plane_spectrum.txt"))["flux"][0]  # this run's flux (Gamma_-12 = g12)
    nH, nHe = delta * u["nH_avg"], delta * u["nHe_avg"]
    t, col, colH, inflight, Tion, bal, rec = [0.0], [0.0], [0.0], [0.0], [], [0.0], [0.0]
    sec = 0.0
    for n in range(1, nsteps + 1):
        sd = os.path.join(out, f"step_{n:04d}")
        st = parse(os.path.join(sd, "rt_stats.txt"))
        s2 = parse(os.path.join(sd, "stats.txt"))
        sec += st["dtsec_rt"][0] * s2["nstep_rt"][0]
        t.append(sec)
        xH, xHe = columns(sd, N)
        col.append(xH + xHe)
        colH.append(xH)
        # the code's running totals of H photoionisations, collisional
        # ionisations and recombinations (rt_stats nion, ncol, nrec: per H
        # atom at the mean density, averaged over the box; cumulative over the
        # run), as columns in cells of the local n_H = delta <n_H>: N / delta
        f = N / delta
        bal.append((st["nion"][0] + st["ncol"][0] - st["nrec"][0]) * f)
        rec.append(st["nrec"][0] * f)
        # photons in flight: per <nH> per cell volume, over N^2 face cells -> a column in cells of n_H
        inflight.append(s2["ngam_rays_total"][0] / N ** 2 / delta)
        fHI = field(sd, "fHI", N)
        temp = field(sd, "temp", N)
        ion = fHI < 0.01  # the ionised interior, without the front cells
        Tion.append(temp[ion].mean() if ion.any() else 1e4)
    t, col, inflight = np.array(t), np.array(col), np.array(inflight)
    Tm = float(np.mean(Tion))
    alpha = 2.59e-13 * (Tm / 1e4) ** -0.7  # case B, Hui & Gnedin 1997 form
    ne = nH + nHe
    x, pred = 0.0, [0.0]
    for i in range(1, len(t)):
        dt = (t[i] - t[i - 1]) / 400
        for _ in range(400):
            x += (F - alpha * ne * nH * x * u["len_rt"]) / nH / u["len_rt"] * dt
        pred.append(x - inflight[i])
    pred = np.array(pred)
    colH, bal, rec = np.array(colH), np.array(bal), np.array(rec)
    sel = col > 3
    err = np.abs(col[sel] / pred[sel] - 1)
    lS = F / (alpha * ne * nH) / u["len_rt"]
    print(f"    {name}: t [Myr] {np.round(t[1:] / 3.156e13, 2)}\n    column {np.round(col[1:], 2)}\n    budget {np.round(pred[1:], 2)}"
          f"  (T_ion = {Tm:.0f} K, Stromgren length {lS:.1f} cells, t_rec = {1 / (alpha * ne) / 3.156e13:.1f} Myr)")
    # the H column against the code's running totals: photoionisations +
    # collisional ionisations - recombinations (the units of the statistics)
    errb = np.abs(colH[sel] / bal[sel] - 1)
    check(sel.sum() >= 3 and np.max(errb) < 0.01, f"{name}: the H column = (nion + ncol - nrec)/n_H of rt_stats to 1 per cent ({np.max(errb) if sel.any() else -1:.4f})")
    if rec[-1] > 0.5:
        # the effective recombination coefficient of the run, alpha_eff =
        # recombinations / int n_e n_H x dt with the measured column, against
        # case B at the mean interior temperature (which varies along the
        # column and in time; the fit is accurate to ~10 per cent)
        integ = np.sum(0.5 * (col[1:] + col[:-1]) * u["len_rt"] * nH * ne * np.diff(t))
        alpha_eff = rec[-1] * nH * u["len_rt"] / integ
        check(1 / 1.5 < alpha_eff / alpha < 1.5, f"{name}: alpha_eff = {alpha_eff:.2e} cm^3/s within a factor 1.5 of case B at T_ion ({alpha:.2e}); the photon budget with the fit and without collisional ionisations is off by {np.max(err):.3f}")
    else:
        check(sel.sum() >= 3 and np.max(err) < tol, f"{name}: the ionised column follows the photon budget to {tol * 100:.0f} per cent ({np.max(err) if sel.any() else -1:.3f})")
    check(col[-1] < N - 2, f"{name}: the front has not reached the far face")
    return out


def check_front(driver, work, spec):
    print("front: one front across the whole box against the photon budget")
    N, box, z, g12 = 32, 0.25, 7.0, 0.3
    # thin: recombinations negligible (t_rec ~ 2 Gyr), the front at ~0.03 c
    out = budget_run(driver, work, "front", spec, N, box, z, g12, 1.0, 32, 8, 1.0, 0.03)
    # dense: 100x the mean density and 10x the flux, c_hat = 0.1 c (the front
    # at ~0.004 c): recombinations matter (t_rec ~ 20 Myr, l_S ~ 16 cells)
    budget_run(driver, work, "front_dense", spec, N, box, z, 10 * g12, 100.0, 100, 12, 0.1, 0.15)
    # light-limited front: c_reduced = 0.1 and 100x the flux
    u = units(z, box, N)
    rows = np.tile(gas_row(u, 1.0, 1.0, 1.0, 0.0, 100.0), (N, N, N, 1))
    prob2 = os.path.join(work, "front_rsla")
    write_problem(prob2, N, box, z, rows, da_for_substeps(z, box, N, 24, 0.1), 1,
                  {"plane_ndom": 1, "plane_faces": "+x", "plane_gamma12": 100 * g12, "c_reduced": 0.1})
    out2 = run(driver, prob2, prob2 + "_run")
    sd = os.path.join(out2, "step_0001")
    nstep = parse(os.path.join(sd, "stats.txt"))["nstep_rt"][0]
    xI2 = columns(sd, N)[0]
    # rays move one cell per substep: a light-limited front sits at nstep cells
    check(abs(xI2 / nstep - 1) < 0.1, f"light-limited front: x_I = {xI2:.2f} cells after {nstep:.0f} substeps at c_hat (to 10 per cent)")
    F = spec["flux"][0]
    thin = F * 100 / (u["nH_avg"] + u["nHe_avg"]) * (24 * u["len_rt"] / (0.1 * C_CGS)) / u["len_rt"]
    check(thin > 3 * nstep, f"(the thin front would be at {thin:.0f} cells: limited by c_hat)")
    return out


def check_conserve(runs):
    print("conserve: emitted = absorbed + escaped + d(in flight)")
    for out in runs:
        steps = sorted(d for d in os.listdir(out) if d.startswith("step_"))
        prev = 0.0
        worst = 0.0
        for d in steps:
            st = parse(os.path.join(out, d, "stats.txt"))
            inflight = st["ngam_rays_total"][0]
            em, ab, es = st["plane_emitted"][0], st["plane_absorbed"][0], st["plane_escaped"][0]
            res = (em - ab - es - (inflight - prev)) / em
            worst = max(worst, abs(res))
            prev = inflight
        check(worst < 1e-3, f"{os.path.basename(out)}: |residual| / emitted <= {worst:.1e} over {len(steps)} steps")


def check_domains(driver, work):
    print("domains: plane_ndom = 2, faces +x +y")
    N, box, z, g12 = 32, 0.25, 7.0, 0.3
    L = N // 2
    u = units(z, box, N)
    rows = np.tile(gas_row(u, 1.0, 1.0, 1.0, 0.0, 100.0), (N, N, N, 1))
    prob = os.path.join(work, "domains")
    write_problem(prob, N, box, z, rows, da_for_substeps(z, box, N, 16, 1.0), 6,
                  {"plane_ndom": 2, "plane_faces": "+x +y", "plane_gamma12": g12, "c_reduced": 1})
    out = run(driver, prob, prob + "_run")
    worst = 0.0
    for n in range(1, 7):
        fHI = field(os.path.join(out, f"step_{n:04d}"), "fHI", N)
        xs = []
        for a in range(2):
            for b in range(2):
                for c in range(2):
                    xs.append(np.mean(1 - fHI[c * L:(c + 1) * L, b * L:(b + 1) * L, a * L:(a + 1) * L]))
        xs = np.array(xs)
        if xs.mean() > 0:
            worst = max(worst, np.max(np.abs(xs / xs.mean() - 1)))
    check(worst < 1e-5 and xs.mean() > 0.05, f"the 8 domains ionise identically (x_HII,V spread {worst:.1e}, mean {xs.mean():.3f})")
    # no ray crosses a domain edge: the thin ionised box
    N2 = 16
    u2 = units(z, 0.1, N2)
    rows2 = np.tile(gas_row(u2, 1.0, 1e-4, 1e-4, 1.0, 1e4), (N2, N2, N2, 1))
    prob2 = os.path.join(work, "domains_thin")
    write_problem(prob2, N2, 0.1, z, rows2, da_for_substeps(z, 0.1, N2, 40, 1.0), 1,
                  {"plane_ndom": 2, "plane_faces": "+x +y", "plane_gamma12": g12})
    out2 = run(driver, prob2, prob2 + "_run")
    gam = field(os.path.join(out2, "step_0001"), "gamma", N2, 3)[..., 0] / (g12 * 1e-12)
    # every cell is lit by one +x and one +y ray per substep: 2 targets; the
    # source planes of the second domain (i or j = N2/2) would get 3 if the
    # first domain's rays were not deleted at its far face
    err = np.max(np.abs(gam / 2.0 - 1))
    check(err < 2e-3, f"thin box: every cell gets both directions' flux once, Gamma_HI = 2 targets ({err:.1e}); the second domain's source planes (i, j = {N2 // 2}) are not lit twice")
    return out


def check_threads(driver, work, spec):
    print("threads: deterministic_parallel = 1 with 1 and 2 threads")
    N, box, z, g12 = 32, 0.25, 7.0, 0.3
    u = units(z, box, N)
    rows = np.tile(gas_row(u, 1.0, 1.0, 1.0, 0.0, 100.0), (N, N, N, 1))
    prob = os.path.join(work, "threads")
    write_problem(prob, N, box, z, rows, da_for_substeps(z, box, N, 16, 1.0), 3,
                  {"plane_ndom": 2, "plane_faces": "+x +y", "plane_gamma12": g12, "deterministic_parallel": 1})
    outs = [run(driver, prob, prob + f"_t{t}", threads=t) for t in (1, 2)]
    same = True
    for d in sorted(os.listdir(outs[0])):
        if not d.startswith("step_"):
            continue
        for f in sorted(os.listdir(os.path.join(outs[0], d))):
            a = open(os.path.join(outs[0], d, f), "rb").read()
            b = open(os.path.join(outs[1], d, f), "rb").read()
            if a != b:
                same = False
                print(f"    differs: {d}/{f}")
    check(same, "1 and 2 threads byte-identical")


def check_sweep(driver, work):
    print("sweep: plane_transport = sweep against rays")
    N, box, z, g12 = 32, 0.25, 7.0, 0.3
    u = units(z, box, N)
    rows = np.tile(gas_row(u, 1.0, 1.0, 1.0, 0.0, 100.0), (N, N, N, 1))
    cases = [("sweep_front", 1, "+x", 1), ("sweep_dom", 2, "+x +y", 1), ("sweep_dom_k2", 2, "+x +y", 2)]
    for name, ndom, faces, k in cases:
        prob = os.path.join(work, name)
        write_problem(prob, N, box, z, rows, da_for_substeps(z, box, N, 16 * k, 1.0), 4,
                      {"plane_ndom": ndom, "plane_faces": faces, "plane_gamma12": g12, "c_reduced": 1,
                       "rt_cells_per_substep": k})
        outs = [run(driver, prob, prob + "_" + t, sets=[f"plane_transport={t}"]) for t in ("rays", "sweep")]
        worst = {"fHI": 0.0, "temp": 0.0, "ngam": 0.0}
        wst = 0.0
        for n in range(1, 5):
            sds = [os.path.join(o, f"step_{n:04d}") for o in outs]
            for fld, nc in (("fHI", 1), ("temp", 1), ("ngam", 5)):
                a, b = (field(sd, fld, N, nc) for sd in sds)
                sc = np.maximum(np.abs(a), np.abs(b)).max()
                worst[fld] = max(worst[fld], float(np.max(np.abs(a - b)) / sc) if sc > 0 else 0.0)
            sa, sb = (parse(os.path.join(sd, "stats.txt")) for sd in sds)
            for key in ("plane_emitted", "plane_absorbed", "plane_escaped", "ngam_rays_total"):
                ref = max(abs(sa[key][0]), abs(sb[key][0]), 1e-300)
                wst = max(wst, abs(sa[key][0] - sb[key][0]) / ref if ref > 1e-12 else 0.0)
        # float round-off through the front's nonlinear response: the two
        # transports sum the same deposits in the same order per cell, but
        # the packets' positions and radii are rebuilt from cell indices
        check(max(worst.values()) < 3e-5, f"{name}: fHI, temp, ngam agree to {max(worst.values()):.1e} ({worst})")
        check(wst < 2e-4, f"{name}: the photon bookkeeping agrees to {wst:.1e}")


def main():
    global EXTRA, LAUNCH, THREADS
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("driver")
    ap.add_argument("problem")
    ap.add_argument("work")
    ap.add_argument("--only", nargs="*", default=None)
    ap.add_argument("--set", action="append", default=[])
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--launch", nargs=4, metavar=("LAUNCHER", "FLAGS", "NPFLAG", "NP"), default=None)
    a = ap.parse_args()
    EXTRA = [x for kv in a.set for x in ("--set", kv)]
    THREADS = a.threads
    if a.launch:
        LAUNCH = [a.launch[0]] + [f for f in a.launch[1].split(";") if f] + [a.launch[2], a.launch[3]]
    os.makedirs(a.work, exist_ok=True)
    g = np.fromfile(os.path.join(a.problem, "gas.bin"), dtype=np.float32).reshape(-1, 9)
    if not np.all(g == g[0]):
        sys.exit("the problem gas is not uniform: " + a.problem)
    only = set(a.only) if a.only else None
    want = lambda n: only is None or n in only
    spec = check_gamma(a.driver, a.work, g[0]) if want("gamma") or want("slab") or want("front") or want("threads") else None
    runs = []
    if want("slab"):
        runs.append(check_slab(a.driver, a.work, spec))
    if want("front"):
        runs.append(check_front(a.driver, a.work, spec))
    if want("conserve"):
        check_conserve(runs)
    if want("domains"):
        check_domains(a.driver, a.work)
    if want("threads"):
        check_threads(a.driver, a.work, spec)
    if want("sweep"):
        check_sweep(a.driver, a.work)
    if FAILS:
        print(f"\n{len(FAILS)} check(s) failed")
        sys.exit(1)
    print("\nall plane-parallel checks passed")


if __name__ == "__main__":
    main()
