"""Cosmological problems of the hydro suite: Zel'dovich pancake with gas and
dark matter, and the adiabatic cosmological box (docs/hydro_tests.md).

The cosmology is the code's (cosmo_parameters.f90: Omega_m = 0.27,
Omega_L = 0.73, Omega_b = 0.045, h = 0.7; compile-time constants of the
solver, not changed here). Units are the code's supercomoving units
(tools/ics/pmcosmo.py).
"""
import math
import os
import shutil
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools", "ics"))
import hs_common as hc  # noqa: E402
import pmcosmo as pc  # noqa: E402

# physical_constants.f90 values used by calc_sim_var
K_CGS = 1.38066e-16
MH_CGS = 1.67223e-24
H0_CGS = 3.24086e-18
MPC2CM = 3.08560e24


def units(a, box, N):
    """calc_sim_var: x_unit, time_unit, vel_unit, temp1_unit (neutral gas)"""
    h0 = pc.HUBBLE0
    x_unit = a * box / h0 / N * MPC2CM
    time_unit = 2 / (3 * H0_CGS * h0 * math.sqrt(pc.OMEGAM)) * a * a
    vel_unit = x_unit / time_unit
    temp1_unit = vel_unit ** 2 * (hc.GAMMA - 1) / K_CGS * MH_CGS / (hc.XHY + hc.YHE / 4.0)
    return {"x_unit": x_unit, "time_unit": time_unit, "vel_unit": vel_unit,
            "temp1_unit": temp1_unit}


def tau_of_a(a_grid):
    """code time tau(a) - tau(a_grid[0]) (dtau = da / (a dlna/dtau))"""
    a = np.asarray(a_grid, float)
    f = 1.0 / (a * pc.dlna_dtau(a))
    return np.concatenate(([0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(a))))


def zeldovich_1d(x, A, k, c, iters=60):
    """Lagrangian q of the Eulerian position x for x - c = (q - c) - A sin(k(q-c))/k
    (A < 1: unique). Returns q."""
    q = x.copy()
    for _ in range(iters):
        f = q - A * np.sin(k * (q - c)) / k - x
        q = q - f / (1 - A * np.cos(k * (q - c)))
    return q


def gen_pancake(outdir, N=64, Ndm=None, box=32.0, z_init=50.0, z_c=5.0, T_init=1000.0,
                z_out=(20, 10, 7, 5.5, 5, 4, 3, 2, 1, 0.5, 0), bulk_frac=0.05,
                dlna_max=0.02):
    """Zel'dovich pancake along x (one wavelength = box, comoving box Mpc/h)
    with gas and dark matter, both on the exact growing mode (all matter
    clusters, Growth(1)), shell crossing at z_c, collapse plane at the box
    centre. Gas: exact Zel'dovich density/velocity at the cell centres,
    T = T_init (rho/rho_mean)^(2/3) at z_init (uniform entropy). A bulk
    velocity in y and z (bulk_frac times the peak infall speed at z_c), the
    same for gas and DM, avoids the zero-velocity pencils of quirk H1."""
    Ndm = Ndm or N
    g = pc.Growth(1.0)
    a_i, a_c = 1 / (1 + z_init), 1 / (1 + z_c)
    A_i = g(a_i) / g(a_c)
    Adot = g.dtau(a_i) / g(a_c)            # d(A)/dtau at a_i
    k = 2 * np.pi / Ndm                    # DM-cell units
    c = Ndm / 2.0
    # dark matter: lattice q = n + 1/2
    q = np.arange(Ndm) + 0.5
    Q = np.meshgrid(q, q, q, indexing="ij")
    vmax_c = g.dtau(a_c) / g(a_c) / k      # peak infall speed at z_c (DM cells / time)
    bulk_dm = np.array([0.0, 0.7 * bulk_frac * vmax_c, -bulk_frac * vmax_c])
    xv = np.zeros((6, Ndm, Ndm, Ndm))
    s = np.sin(k * (Q[0] - c)) / k
    xv[0] = np.mod(Q[0] - A_i * s, Ndm)
    xv[1] = Q[1]
    xv[2] = Q[2]
    xv[3] = -Adot * s
    xv[4] = bulk_dm[1]
    xv[5] = bulk_dm[2]
    # gas: cell centres in DM units
    r = N / Ndm
    xg = (np.arange(N) + 0.5) / r
    qg = zeldovich_1d(xg, A_i, k, c)
    fb = pc.FB
    rho1 = fb / (1 - A_i * np.cos(k * (qg - c)))
    v1 = -Adot * np.sin(k * (qg - c)) / k * r            # gas cells / time
    U = units(a_i, box, N)
    T1 = T_init * (rho1 / fb) ** (2.0 / 3.0)
    eth1 = rho1 * T1 / U["temp1_unit"]
    shape = (N, N, N)
    rho = np.broadcast_to(rho1[:, None, None], shape).copy()
    eth = np.broadcast_to(eth1[:, None, None], shape).copy()
    vel = np.zeros((3,) + shape)
    vel[0] = v1[:, None, None]
    vel[1] = bulk_dm[1] * r
    vel[2] = bulk_dm[2] * r
    params = {"N": N, "Ndm": Ndm, "box": box, "mode": "cosmo", "z_init": z_init,
              "z_out": list(z_out), "dlna_max": dlna_max, "pk": 0}
    hc.write_problem(outdir, params, rho, vel, eth,
                     extra_files={"xv_dm.init": xv.reshape(6, -1, order="F")})
    meta = dict(test="pancake", N=N, Ndm=Ndm, box=box, z_init=z_init, z_c=z_c, T_init=T_init,
                z_out=list(z_out), A_i=A_i, k_dm=k, centre_dm=c, bulk_gas=list(bulk_dm * r),
                dlna_max=dlna_max, fb=fb)
    hc.save_json(os.path.join(outdir, "meta.json"), meta)
    return meta


def gen_cosmobox(outdir, N=64, Ndm=None, box=64.0, z_init=50.0, seed=7,
                 z_out=(10, 5, 3, 2, 1, 0.5, 0), dlna_max=None):
    """Adiabatic cosmological box: Zel'dovich ICs from tools/ics/make_ics.py
    (Eisenstein-Hu P(k), sigma8 = 0.8, gas follows the DM), T = T_init
    (1+delta)^(2/3) with T_init the adiabatically decoupled CMB temperature
    T_cmb0 (1+z)^2/(1+z_dec), z_dec = 150. No radiative transfer, no
    cooling. dlna_max = None keeps the original da/a <= 1/9 limit."""
    Ndm = Ndm or N // 2
    os.makedirs(outdir, exist_ok=True)
    tmp = os.path.join(outdir, "ics")
    subprocess.run([sys.executable, os.path.join(HERE, "..", "..", "tools", "ics", "make_ics.py"),
                    tmp, "--kind", "zeldovich", "--Ngas", str(N), "--Ndm", str(Ndm),
                    "--z-init", str(z_init), "--seed", str(seed), "--box", str(box)],
                   check=True, capture_output=True)
    u = np.fromfile(os.path.join(tmp, "u_gas.init"), dtype="<f4").reshape(
        (5, N, N, N), order="F").astype(np.float64)
    a_i = 1 / (1 + z_init)
    T_init = pc.TCMB0 * (1 + z_init) ** 2 / 151.0
    U = units(a_i, box, N)
    rho = u[0]
    eth = rho * T_init * (rho / pc.FB) ** (2.0 / 3.0) / U["temp1_unit"]
    u[4] = u[4] + eth
    hc.write_f32(os.path.join(outdir, "u_gas.bin"), u)
    shutil.copy(os.path.join(tmp, "xv_dm.init"), os.path.join(outdir, "xv_dm.init"))
    params = {"N": N, "Ndm": Ndm, "box": box, "mode": "cosmo", "z_init": z_init,
              "z_out": " ".join(repr(float(z)) for z in z_out), "pk": 1}
    if dlna_max is not None:
        params["dlna_max"] = dlna_max
    with open(os.path.join(outdir, "params.txt"), "w") as f:
        for k_, v in params.items():
            f.write(f"{k_} = {v}\n")
    meta = dict(test="cosmobox", N=N, Ndm=Ndm, box=box, z_init=z_init, seed=seed,
                z_out=list(z_out), T_init=T_init, dlna_max=dlna_max)
    hc.save_json(os.path.join(outdir, "meta.json"), meta)
    return meta
