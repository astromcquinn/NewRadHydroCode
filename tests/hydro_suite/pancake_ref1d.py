#!/usr/bin/env python3
"""High-resolution 1-D reference solution of the Zel'dovich pancake with gas
and dark matter (independent of the RadHydro solver).

Lagrangian staggered-mesh hydrodynamics (von Neumann-Richtmyer artificial
viscosity, Caramana-type quadratic + linear terms) for the gas and exact
1-D sheet gravity for gas zones and DM sheets, in the code's supercomoving
variables (no expansion terms for gamma = 5/3; Poisson equation
d2phi/dx2 = (2a/3)(rho - 1), lengths in units of the wavelength, mean total
density 1, gas fraction fb). The pancake is symmetric about x = 1/2 and the
domain [0,1] has reflecting walls at 0 and 1 (= the symmetry planes of the
periodic problem).

usage: pancake_ref1d.py <meta.json of a pancake problem> <out.npz> [--ng 4096]
Writes, at each z of meta["z_out"], the Eulerian gas profile (zone centres,
rho [mean matter density units], v [lambda per code time], e/rho and
T [K]) and the DM sheet positions/velocities.
"""
import argparse
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools", "ics"))
import hs_common as hc  # noqa: E402
import cosmo_gen as cg  # noqa: E402
import pmcosmo as pc  # noqa: E402


def a_of_tau_table(a0, a1, n=200000):
    a = np.exp(np.linspace(np.log(a0), np.log(a1), n))
    return a, cg.tau_of_a(a)


def run(meta, ng=4096, ndm=4096, cq=2.0, cl=0.25, cfl=0.25):
    g = hc.GAMMA
    fb = meta["fb"]
    G = pc.Growth(1.0)
    a_i, a_c = 1 / (1 + meta["z_init"]), 1 / (1 + meta["z_c"])
    A_i = G(a_i) / G(a_c)
    Adot = G.dtau(a_i) / G(a_c)
    k = 2 * np.pi
    # gas zones: Lagrangian boundaries q_b, uniform mass
    qb = np.linspace(0, 1, ng + 1)
    Xb = qb - A_i * np.sin(k * (qb - 0.5)) / k
    Vb = -Adot * np.sin(k * (qb - 0.5)) / k
    m = fb / ng * np.ones(ng)
    rho = m / np.diff(Xb)
    qc = 0.5 * (qb[1:] + qb[:-1])
    U = cg.units(a_i, meta["box"], 1)            # lambda units (N = 1)
    rho_ex = fb / (1 - A_i * np.cos(k * (qc - 0.5)))
    e = meta["T_init"] * (rho_ex / fb) ** (2 / 3) / U["temp1_unit"]   # e per mass
    K = e / rho ** (2 / 3)                       # code-unit entropy (initial)
    K0 = meta["T_init"] / cg.units(a_i, meta["box"], 1)["temp1_unit"] / fb ** (2 / 3)
    # DM sheets
    qd = (np.arange(ndm) + 0.5) / ndm
    Xd = qd - A_i * np.sin(k * (qd - 0.5)) / k
    Vd = -Adot * np.sin(k * (qd - 0.5)) / k
    md = (1 - fb) / ndm

    a_tab, tau_tab = a_of_tau_table(a_i, 1.0)
    z_out = sorted(meta["z_out"], reverse=True)
    tau_out = [float(np.interp(math.log(1 / (1 + z)), np.log(a_tab), tau_tab)) for z in z_out]
    tau_end = tau_out[-1]

    def accel_grav(a, xs, ms):
        """g(x) at the positions xs of all masses (sorted cumulative mass)"""
        order = np.argsort(xs, kind="stable")
        mcum = np.cumsum(ms[order]) - 0.5 * ms[order]
        gx = np.empty_like(xs)
        gx[order] = -(2 * a / 3) * ((mcum - 0.5) - (xs[order] - 0.5))
        return gx

    def forces(a, Xb, Xd, e, rho, qv):
        p = (g - 1) * rho * e
        P = p + qv
        xc = 0.5 * (Xb[1:] + Xb[:-1])
        allx = np.concatenate((xc, Xd))
        allm = np.concatenate((m, np.full(ndm, md)))
        gall = accel_grav(a, allx, allm)
        gz = gall[:ng]
        gd = gall[ng:]
        # boundary acceleration: pressure + gravity (interpolated zone gravity)
        ab = np.zeros(ng + 1)
        mb = 0.5 * (m[1:] + m[:-1])
        ab[1:-1] = -(P[1:] - P[:-1]) / mb + 0.5 * (gz[1:] + gz[:-1])
        return ab, gd, p

    def visc(rho, dVb, cs):
        """artificial viscosity with a shock switch: only in zones whose
        velocity jump is converging and much larger than that of the zones
        four zones away on both sides (a smooth cold infall, where any
        viscosity would destroy the entropy, never triggers it)"""
        du = dVb  # velocity difference across the zone
        ad = np.abs(du)
        n = ad.size
        i = np.arange(n)
        nb = np.maximum(ad[np.clip(i - 4, 0, n - 1)], ad[np.clip(i + 4, 0, n - 1)])
        shock = (du < 0) & (ad > 4 * nb)
        # widen the switch by one zone on each side
        shock = shock | np.roll(shock, 1) | np.roll(shock, -1)
        return np.where(shock & (du < 0), rho * (cq * du * du + cl * cs * ad), 0.0)

    tau = 0.0
    nstep = 0
    outs = []
    io = 0
    qv = np.zeros(ng)
    t0 = time.time()
    a = a_i
    ab, gd, p = forces(a, Xb, Xd, e, rho, qv)
    while io < len(tau_out):
        cs = np.sqrt(g * (g - 1) * np.maximum(e, 0))
        dx = np.diff(Xb)
        dvz = np.abs(np.diff(Vb))
        dt_h = cfl * np.min(dx / (cs + 2 * dvz + 1e-30))
        dt_g = 0.02 / math.sqrt(2 * a / 3 * max(rho.max(), 1.0))
        dt = min(dt_h, dt_g, tau_out[io] - tau)
        # kick-drift-kick
        a_mid = float(np.interp(tau + dt / 2, tau_tab, a_tab))
        Vb_h = Vb + 0.5 * dt * ab
        Vd_h = Vd + 0.5 * dt * gd
        Vb_h[0] = Vb_h[-1] = 0
        Xb_new = Xb + dt * Vb_h
        Xd = Xd + dt * Vd_h
        lo, hi = Xd < 0, Xd > 1                     # reflect at the walls
        Xd = np.where(lo, -Xd, np.where(hi, 2 - Xd, Xd))
        Vd_h = np.where(lo | hi, -Vd_h, Vd_h)
        rho_new = m / np.diff(Xb_new)
        qv = visc(0.5 * (rho + rho_new), np.diff(Vb_h), cs)
        dvol = 1 / rho_new - 1 / rho                # specific volume change
        e = (e - (0.5 * p + qv) * dvol) / (1 + 0.5 * (g - 1) * rho_new * dvol)
        Xb, rho = Xb_new, rho_new
        tau += dt
        a = float(np.interp(tau, tau_tab, a_tab))
        ab, gd, p = forces(a, Xb, Xd, e, rho, qv)
        Vb = Vb_h + 0.5 * dt * ab
        Vb[0] = Vb[-1] = 0
        Vd = Vd_h + 0.5 * dt * gd
        nstep += 1
        if tau >= tau_out[io] - 1e-14:
            Ua = cg.units(a, meta["box"], 1)
            xc = 0.5 * (Xb[1:] + Xb[:-1])
            outs.append({"z": z_out[io], "a": a, "tau": tau, "x": xc.copy(), "rho": rho.copy(),
                         "v": 0.5 * (Vb[1:] + Vb[:-1]), "e": e.copy(),
                         "T": e * Ua["temp1_unit"], "K": e / rho ** (2 / 3) / K0,
                         "xd": Xd.copy(), "vd": Vd.copy(), "nstep": nstep,
                         "wall": time.time() - t0})
            print(f"  ref1d: z = {z_out[io]:5.2f} tau = {tau:.4f} steps = {nstep} "
                  f"rho_max = {rho.max():.3g} T_max = {outs[-1]['T'].max():.3g} K "
                  f"({time.time() - t0:.1f} s)", flush=True)
            io += 1
    return outs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("meta")
    ap.add_argument("out")
    ap.add_argument("--ng", type=int, default=4096)
    args = ap.parse_args()
    meta = hc.load_json(args.meta)
    outs = run(meta, ng=args.ng, ndm=args.ng)
    flat = {}
    for i, o in enumerate(outs):
        for k_, v in o.items():
            flat[f"{i}_{k_}"] = v
    flat["n"] = len(outs)
    np.savez_compressed(args.out, **flat)
    return 0


if __name__ == "__main__":
    sys.exit(main())
