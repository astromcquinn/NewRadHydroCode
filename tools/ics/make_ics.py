#!/usr/bin/env python3
"""Initial conditions for the RadHydro PM (gravity) port.

Writes the two files that io_tools.f90:read_init_conditions reads (the
original expects them in ``Init/``):

  xv_dm.init  float32 xv_dm(6,Npart), Fortran order, stream (no record
              markers): x,y,z in DM grid cells [0,Ndm), vx,vy,vz = dx/dtau in
              DM grid cells per code time unit;
  u_gas.init  float32 u_gas(5,Ngas,Ngas,Ngas): rho (units of the mean total
              matter density, mean Omega_b/Omega_m), rho*v (v in GAS grid
              cells per time unit), and the kinetic energy rho v^2/2.

What the code does with them (read_init_conditions, init_hydro):
  * u_gas(5) += u_gas(1)*Tgas_min, so the file holds the energy WITHOUT the
    thermal part (thermal energy is set to the minimum temperature);
  * w_gas is not read: it is set to fully neutral, zero metallicity
    (w = XHy rho, YHe rho, 0, 0);
  * v_gas (the MF grid velocity) is not read: init_hydro zeroes it, so the
    full gas velocity has to be in u_gas(2:4);
  * NOTE: with Ncpu = 1 the original never reads u_gas.init (the parallel
    read loop only reaches case(2) when Ncpu >= 2); see
    reference/harness_pm/PATCHES.md.

Particles start on a lattice at q = (i+1/2, j+1/2, k+1/2) DM cells. Gas cell
(i,j,k) (1-based) is centred at gas coordinate i-1/2, i.e. DM coordinate
(i-1/2)/N_dm_gas, which is how pm_tools.f90 maps particles onto the gas grid
(x_gas = x_dm*x_dm_gas - 0.5).

Kinds:
  zeldovich  Gaussian random field with the Eisenstein-Hu no-wiggle linear
             P(k) (n_s, sigma8 of cosmo_parameters.f90) at z_init,
             Zel'dovich displacements x = q + psi, v = (dD/dtau / D) psi.
  pancake    one plane wave along --axis: x = q - (D/D_c) sin(k q)/k, with
             k the box fundamental; shell crossing at a = --a-cross.  The 1D
             Zel'dovich solution is exact until shell crossing.

Gas:
  zeldovich  gas density and velocity follow the DM linear field (Eulerian,
             at the gas cell centres);
  uniform    rho = Omega_b/Omega_m, v = 0.

Growth (velocities and the pancake's D(a)):
  full        all matter clusters (D'' = (2/3) a D), the usual choice;
  frozen-gas  only the DM clusters (D'' = (2/3)(1-f_b) a D): the exact
              growing mode of a gravity-only run (no hydro), where the gas
              density never moves.  Use with --gas uniform for analytic
              checks.

More particles than gas cells (Ndm = 2 Ngas with pm_refine = 2,
docs/pm_refine.md): the gas fields keep only the modes |m| < Ngas/2 of the
particles' field (the gas grid cannot hold the rest).

--phases-from-ndm N (zeldovich, N < Ndm): the modes |m| < N/2 of the white
noise are those of an Ndm = N run with the same --seed (the same draws,
rescaled), the rest come from an independent stream. The large scales then
match that run's ICs (e.g. the 512^3 production ICs, seed 7, under a
1024^3-particle pm_refine = 2 run).

Usage:
  make_ics.py OUTDIR --kind zeldovich --Ngas 32 [--Ndm 16] [--z-init 300]
              [--seed 1] [--gas zeldovich|uniform] [--growth full|frozen-gas]
              [--phases-from-ndm N]
  make_ics.py OUTDIR --kind pancake --axis x --a-cross 0.05 ...
Also writes OUTDIR/ics_info.txt with the parameters used.
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pmcosmo as pc  # noqa: E402


def kgrid(n):
    """Integer wavenumbers m (fftfreq*n) along one axis, Nyquist included."""
    return np.fft.fftfreq(n, d=1.0 / n)


def eval_on_lattice(field_k, nsrc, nout, offset):
    """Evaluate the real field with Fourier coefficients field_k (numpy fftn
    convention on an nsrc^3 grid, Nyquist already zeroed) at the points
    x = (n + offset)*nsrc/nout, n = 0..nout-1 (in units of source cells).
    """
    m = kgrid(nsrc)
    if nout < nsrc:
        # the output grid holds only |m| < nout/2 (its Nyquist is dropped too)
        sel = np.nonzero(np.abs(m) < nout / 2)[0]
        field_k = field_k[np.ix_(sel, sel, sel)]
        m = m[sel]
    big = np.zeros((nout, nout, nout), dtype=complex)
    idx = np.where(m >= 0, m, m + nout).astype(int)
    ph = np.exp(2j * np.pi * m * offset / nout)       # shift by offset output cells
    shifted = field_k * ph[:, None, None] * ph[None, :, None] * ph[None, None, :]
    big[np.ix_(idx, idx, idx)] = shifted
    out = np.fft.ifftn(big) * (nout ** 3 / nsrc ** 3)
    return out.real


def zeldovich(args, ndm, ngas, a_init):
    rng = np.random.default_rng(args.seed)
    nc = ndm ** 3
    L = args.box if args.box is not None else pc.BOX
    g_full = pc.Growth(1.0)
    g = pc.Growth(1.0 if args.growth == "full" else 1.0 - pc.FB)

    if args.phases_from_ndm:
        nc0 = args.phases_from_ndm
        if not (0 < nc0 < ndm and nc0 % 2 == 0):
            sys.exit("--phases-from-ndm must be even and below Ndm")
        # the fine modes from an independent stream, the coarse ones from the
        # Ndm = nc0 run's draws (E|w_k|^2 = n^3: rescaled by (ndm/nc0)^1.5)
        wk = np.fft.fftn(np.random.default_rng((args.seed, ndm)).standard_normal((ndm,) * 3))
        wc = np.fft.fftn(rng.standard_normal((nc0,) * 3))
        mc = kgrid(nc0)
        ic = np.nonzero(np.abs(mc) < nc0 // 2)[0]
        mk = mc[ic].astype(int)
        jf = np.where(mk >= 0, mk, mk + ndm)
        wk[np.ix_(jf, jf, jf)] = wc[np.ix_(ic, ic, ic)] * (ndm / nc0) ** 1.5
        del wc
    else:
        white = rng.standard_normal((ndm, ndm, ndm))
        wk = np.fft.fftn(white)
        del white
    m = kgrid(ndm)
    mx, my, mz = np.meshgrid(m, m, m, indexing="ij")
    msq = mx ** 2 + my ** 2 + mz ** 2
    kphys = 2 * np.pi / L * np.sqrt(msq)
    pk = np.zeros_like(kphys)
    nz = msq > 0
    # amplitude at z_init: the usual normalisation with the full growth factor
    D_ratio = g_full(a_init)
    pk[nz] = pc.power_z0(kphys[nz]) * D_ratio ** 2
    delta_k = wk * np.sqrt(nc * pk / L ** 3)
    # zero the Nyquist planes so that every field is real and odd modes vanish
    nyq = (np.abs(mx) == ndm // 2) | (np.abs(my) == ndm // 2) | (np.abs(mz) == ndm // 2)
    delta_k[nyq] = 0
    delta_k[0, 0, 0] = 0
    del wk, kphys, pk, nyq

    # displacement in DM cells: psi_k = i kvec delta_k / k^2 (k per DM cell)
    kg = 2 * np.pi / ndm
    ksq = np.where(nz, msq * kg ** 2, 1.0)
    psi_k = [1j * (kg * mm) * delta_k / ksq for mm in (mx, my, mz)]

    vfac = g.dtau(a_init) / g(a_init)       # v = (dD/dtau)/D * psi
    sig = np.std(np.fft.ifftn(delta_k).real)

    # particles: q = n + 1/2 (DM cells)
    psi_p = [eval_on_lattice(p, ndm, ndm, 0.5) for p in psi_k]
    q = np.arange(ndm) + 0.5
    Q = np.meshgrid(q, q, q, indexing="ij")
    xv = np.empty((6, ndm, ndm, ndm))
    for c in range(3):
        xv[c] = np.mod(Q[c] + psi_p[c], ndm)
        xv[3 + c] = vfac * psi_p[c]
    xv = xv.reshape(6, -1, order="F")

    # gas at cell centres (n + 1/2 gas cells)
    fb = pc.FB
    if args.gas == "uniform":
        rho = np.full((ngas,) * 3, fb)
        vel = [np.zeros((ngas,) * 3)] * 3
    else:
        dg = eval_on_lattice(delta_k, ndm, ngas, 0.5)
        rho = fb * (1 + dg)
        xdg = ngas / ndm
        vel = [xdg * vfac * eval_on_lattice(p, ndm, ngas, 0.5) for p in psi_k]
    info = {"sigma_delta_dm_cell": sig, "D_init": g(a_init), "f_init": g.f(a_init),
            "vfac": vfac, "max_disp_dm_cells": float(np.max(np.abs(np.array(psi_p))))}
    return xv, rho, vel, info


def pancake(args, ndm, ngas, a_init):
    g = pc.Growth(1.0 if args.growth == "full" else 1.0 - pc.FB)
    ax = "xyz".index(args.axis)
    k = 2 * np.pi / ndm                     # box fundamental, per DM cell
    amp = g(a_init) / g(args.a_cross)       # D/D_c
    vamp = g.dtau(a_init) / g(args.a_cross)
    q = np.arange(ndm) + 0.5
    Q = np.meshgrid(q, q, q, indexing="ij")
    xv = np.zeros((6, ndm, ndm, ndm))
    for c in range(3):
        xv[c] = Q[c]
    s = np.sin(k * Q[ax]) / k
    xv[ax] = np.mod(Q[ax] - amp * s, ndm)
    xv[3 + ax] = -vamp * s
    xv = xv.reshape(6, -1, order="F")
    fb = pc.FB
    if args.gas == "uniform":
        rho = np.full((ngas,) * 3, fb)
        vel = [np.zeros((ngas,) * 3)] * 3
    else:
        # linear Eulerian gas: delta = (D/Dc) cos(k x), v along the axis
        xg = (np.arange(ngas) + 0.5) * ndm / ngas      # DM units
        X = np.meshgrid(xg, xg, xg, indexing="ij")[ax]
        rho = fb * (1 + amp * np.cos(k * X))
        vel = [np.zeros((ngas,) * 3) for _ in range(3)]
        vel[ax] = -(ngas / ndm) * vamp * np.sin(k * X) / k
    info = {"D_over_Dcross_init": amp, "a_cross": args.a_cross,
            "k_dm": k, "axis": args.axis, "vamp": vamp}
    return xv, rho, vel, info


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir")
    ap.add_argument("--kind", choices=("zeldovich", "pancake"), default="zeldovich")
    ap.add_argument("--Ngas", type=int, default=32)
    ap.add_argument("--Ndm", type=int, default=None,
                    help="default Ngas/2; a divisor or a multiple of Ngas")
    ap.add_argument("--phases-from-ndm", type=int, default=0,
                    help="zeldovich: take the white-noise modes |m| < N/2 from an Ndm = N run "
                         "with the same seed (see above)")
    ap.add_argument("--z-init", type=float, default=pc.ZINIT)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--gas", choices=("zeldovich", "uniform"), default="zeldovich")
    ap.add_argument("--growth", choices=("full", "frozen-gas"), default="full")
    ap.add_argument("--axis", choices=("x", "y", "z"), default="x")
    ap.add_argument("--a-cross", type=float, default=0.05)
    ap.add_argument("--box", type=float, default=None,
                    help="box size in Mpc/h (default 5, cosmo_parameters.f90); only "
                         "the P(k) normalisation depends on it")
    args = ap.parse_args(argv)

    ngas = args.Ngas
    ndm = args.Ndm or ngas // 2
    if ngas % ndm and ndm % ngas:
        ap.error("Ngas must be a multiple or a divisor of Ndm")
    a_init = 1.0 / (1.0 + args.z_init)

    if args.kind == "zeldovich":
        xv, rho, vel, info = zeldovich(args, ndm, ngas, a_init)
    else:
        xv, rho, vel, info = pancake(args, ndm, ngas, a_init)

    u = np.empty((5, ngas, ngas, ngas))
    u[0] = rho
    for c in range(3):
        u[1 + c] = rho * vel[c]
    u[4] = 0.5 * rho * (vel[0] ** 2 + vel[1] ** 2 + vel[2] ** 2)

    os.makedirs(args.outdir, exist_ok=True)
    np.asarray(xv, dtype="<f4").ravel(order="F").tofile(os.path.join(args.outdir, "xv_dm.init"))
    np.asarray(u, dtype="<f4").ravel(order="F").tofile(os.path.join(args.outdir, "u_gas.init"))
    with open(os.path.join(args.outdir, "ics_info.txt"), "w") as f:
        f.write("# written by tools/ics/make_ics.py\n")
        av = list(argv if argv is not None else sys.argv[1:])
        av = [os.path.basename(os.path.normpath(x)) if x == args.outdir else x for x in av]
        f.write("command = make_ics.py " + " ".join(av) + "\n")
        for k, v in vars(args).items():
            if (k == "box" and v is None) or (k == "phases_from_ndm" and not v):
                continue
            if k == "outdir":
                v = os.path.basename(os.path.normpath(v))
            f.write(f"{k} = {v}\n")
        f.write(f"Ndm = {ndm}\na_init = {a_init!r}\n")
        for k, v in info.items():
            f.write(f"{k} = {v}\n")
    print(f"wrote {args.outdir}: Npart = {ndm**3}, Ngas = {ngas}, " +
          ", ".join(f"{k}={v}" for k, v in info.items()))


if __name__ == "__main__":
    main()
