#!/usr/bin/env python3
"""Reionization history and Thomson optical depth.

    reion.py RUN_DIR [--out tau.txt]

reads RUN_DIR/movie/frames.txt (apps/radhydro movie_dz output: frame, istep,
istep_rt, z, <x_HII>_V, <x_HII>_M) and prints tau_es(z), the CMB
electron-scattering optical depth to redshift z, and the total tau_es.

tau_es(z) = int_0^z sigma_T c n_e(z') dt/dz' dz', with
  n_e = n_H0 (1+z)^3 x_HII,M (1 + y He),  y = Y/(4X) = n_He/n_H,
  He singly ionised with H, doubly ionised (1 + 2y) below z_HeII = 3
  (assumption: HeII reionization is not in the run);
  x_HII = 1 below the last frame's redshift (the run ends ionised; checked:
  a warning is printed if <x_HII>_M < 0.99 there);
  H(z) = H0 sqrt(Om (1+z)^3 + OL), the simulation cosmology
  (cosmo_parameters.f90: Om = 0.27, Ob = 0.045, h = 0.7, X = 0.76).
The mass-weighted x_HII is the right weight: n_e is proportional to the
local gas density times x_HII (clumping is included).

Planck 2018 (Planck Collaboration VI 2020, A&A 641, A6): tau = 0.054 +- 0.007.
"""
import argparse
import os

import numpy as np

SIGMA_T = 6.6524587e-25   # cm^2
C = 2.99792458e10         # cm/s
MPC = 3.0856775814913673e24
MH = 1.6735575e-24        # g, hydrogen atom
G = 6.6743e-8
OM, OB, H0, X, Y = 0.27, 0.045, 0.70, 0.76, 0.24
Z_HEII = 3.0
PLANCK_TAU, PLANCK_SIG = 0.054, 0.007


def _H(z):
    H0s = H0 * 100e5 / MPC
    return H0s * np.sqrt(OM * (1 + z) ** 3 + 1 - OM)


def nH0():
    H0s = H0 * 100e5 / MPC
    rho_c = 3 * H0s ** 2 / (8 * np.pi * G)
    return X * OB * rho_c / MH


def read_frames(run_dir):
    f = np.loadtxt(os.path.join(run_dir, "movie", "frames.txt"), ndmin=2)
    return f[:, 0].astype(int), f[:, 3], f[:, 4], f[:, 5]  # frame, z, xV, xM


def tau_history(z_frames, xM_frames, nz=4000):
    """tau_es(z) at the frame redshifts (and the total). Returns (z, tau)."""
    order = np.argsort(z_frames)
    zf, xf = z_frames[order], xM_frames[order]
    y = Y / (4 * X)
    zmax = zf.max()
    zz = np.linspace(0, zmax, nz)
    x = np.where(zz < zf.min(), 1.0, np.interp(zz, zf, xf))
    he = np.where(zz < Z_HEII, 1 + 2 * y, 1 + y)
    integrand = SIGMA_T * C * nH0() * (1 + zz) ** 3 * x * he / ((1 + zz) * _H(zz))
    tau = np.concatenate([[0], np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(zz))])
    return zz, tau


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("--out")
    a = ap.parse_args()
    fr, z, xV, xM = read_frames(a.run_dir)
    if xM[np.argmin(z)] < 0.99:
        print(f"# WARNING: <x_HII>_M = {xM[np.argmin(z)]:.3f} at the last frame z = {z.min():.2f}; "
              "tau assumes full ionization below it")
    zz, tau = tau_history(z, xM)
    tz = np.interp(z, zz, tau)
    lines = ["# frame z xHII_V xHII_M tau_es(<z)"]
    for i in range(len(z)):
        lines.append(f"{fr[i]} {z[i]:.4f} {xV[i]:.6f} {xM[i]:.6f} {tz[i]:.5f}")
    # Midpoint and end of reionization
    o = np.argsort(z)
    for lev in (0.5, 0.99):
        above = np.where(xV[o][::-1] >= lev)[0]  # descending z
        if len(above):
            zs = z[o][::-1]
            k = above[0]
            zl = zs[k] if k == 0 else np.interp(lev, [xV[o][::-1][k - 1], xV[o][::-1][k]],
                                                [zs[k - 1], zs[k]])
            lines.append(f"# z(x_HII,V = {lev}) = {zl:.3f}")
    lines.append(f"# tau_es total = {tau[-1]:.5f}  (Planck 2018: {PLANCK_TAU} +- {PLANCK_SIG})")
    txt = "\n".join(lines) + "\n"
    if a.out:
        open(a.out, "w").write(txt)
    print("\n".join(l for l in lines if l.startswith("#")))


if __name__ == "__main__":
    main()
