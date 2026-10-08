#!/usr/bin/env python3
"""Write the Phase-1 validation test problems (docs/validation_spec.md).

    python3 reference/harness/make_test_problem.py [--N 32] [--outdir tests/problems]
                                                   [--force] [names...]

Problems (default: all three)
  uniform_1src     uniform gas at mean baryon density, 1 source in cell (N/2,N/2,N/2)
  uniform_8src     uniform gas, 8 sources at fixed fractional positions
  lognormal_16src  seeded log-normal density field, 16 sources in the densest cells
  uniform_8src_smallbuf  as uniform_8src with nray_max_factor = 4 (spec problem 4;
                   NOTE: with N=32 it still never merges, too few rays)
  lognormal_16src_smallbuf  as lognormal_16src with nray_max_factor = 1 (merges from step ~5)

Existing problem directories are left untouched unless --force is given
(another tool may own tests/problems; this generator is a fallback for the
Fortran harness).  numpy only.

Files per problem (formats from the spec)
  params.txt   key = value  (N, z_init, da_sim, nsteps[, nray_max_factor])
  gas.bin      float32 little-endian, Fortran order (9,N,N,N):
               u_gas(1:5) = rho, rho*v (=0), E ;  w_gas(1:4) = rho_HI, rho_HeI, rho_HeII, rho_Z
  sources.txt  "i j k Luv sfr" per line, 1-based RT cell indices

Initial state (as reference/Radhydro io_tools.f90: read_init_conditions)
  rho   = (Omega_b/Omega_m) * delta        (code units: rho_unit = mean matter density)
  v     = 0
  E     = rho * Tgas_min,  Tgas_min = temp_min/temp1_unit  (=> T = temp_min, neutral gas)
  w_gas = (XHy*rho, YHe*rho, 0, 0)          (fully neutral, zero metals)
temp_min and temp1_unit are evaluated at z_init exactly as calc_sim_var
(cosmo_tools.f90) does, including the single-precision literals of
cosmo_parameters.f90 (e.g. omegab = 0.045 is the real(4) 0.045 widened to
real(8)).  At z_init = 10: temp_min = 1.516 K.

Source strengths.  ngam emitted per RT substep (per <nH> per RT-cell volume) is
fResc*sum(pop2_init)*Luv*dtsec_rt/(vol_rt*nH_avg).  With N=32, z=10, da=1e-4
(dtsec_rt ~ 2.1e12 s, nstep_rt = 12) Luv = 2e28 gives ~6e52 ionising photons/s,
enough to ionise a sphere of radius ~4-5 RT cells in 10 hydro steps (~8 Myr).
sfr only needs to exceed sfr_min = 1e-6 (cosmo_tools.f90); it marks source cells.
"""
import argparse
import os
import sys

import numpy as np

f32 = lambda x: float(np.float32(x))  # noqa: E731  (Fortran real(4) literal -> real(8))

# cosmo_parameters.f90 (single-precision literals assigned to real(8) parameters)
BOX = 5.0
OMEGAM = f32(0.27)
OMEGAB = f32(0.045)
HUBBLE0 = f32(0.70)
XHY = f32(0.76)
YHE = f32(0.24)
GAMMA_IDEAL = f32(np.float32(5.0) / np.float32(3.0))
TCMB0 = f32(2.725)
# physical_constants.f90 (double literals)
K_CGS = 1.38066e-16
MH_CGS = 1.67223e-24
H0_CGS = 3.24086e-18
MPC2CM = 3.08560e24


def temp_units(z, n):
    """temp_min and Tgas_min (=temp_min/temp1_unit) as in calc_sim_var."""
    a = 1.0 / (1.0 + z)
    x_unit = (a * BOX / HUBBLE0 / n) * MPC2CM
    time_unit = 2.0 / (3.0 * H0_CGS * HUBBLE0 * np.sqrt(OMEGAM)) * a**2
    vel_unit = x_unit / time_unit
    # eden_unit/rho_unit = vel_unit**2
    temp1_unit = vel_unit**2 * (GAMMA_IDEAL - 1) / K_CGS * MH_CGS / (1 * XHY + 1 * YHE / f32(4.0))
    x = 1e3 * (OMEGAB * HUBBLE0**2) ** f32(0.4)
    if 1 + z > x:
        temp_min = TCMB0 / a
    else:
        temp_min = max(TCMB0 / x / a**2, 1.0)
    return temp_min, temp_min / temp1_unit


def make_gas(delta, z, n):
    temp_min, tgas_min = temp_units(z, n)
    rho = (OMEGAB / OMEGAM * delta).astype(np.float32)
    gas = np.zeros((9, n, n, n), dtype=np.float32, order="F")
    gas[0] = rho
    # u(5) = 0 + u(1)*Tgas_min  (real(4) = real(4)*real(8))
    gas[4] = (rho.astype(np.float64) * tgas_min).astype(np.float32)
    gas[5] = (XHY * rho.astype(np.float64)).astype(np.float32)
    gas[6] = (YHE * rho.astype(np.float64)).astype(np.float32)
    return gas, temp_min


def lognormal_delta(n, seed=20121, sigma=1.3, kcut=0.5):
    """Deterministic log-normal overdensity field with <delta> = 1.

    Gaussian field: white noise filtered with P(k) ~ k^-2 * exp(-(k/kcut)^2)
    (k in units of the Nyquist frequency), scaled to unit variance; then
    delta = exp(sigma*G) / mean(exp(sigma*G)), clipped to >= 1e-3 so that
    rho > Dgas_min = 1e-4 <rhob> (no density floor is triggered).
    """
    rng = np.random.default_rng(seed)
    g = rng.standard_normal((n, n, n))
    gk = np.fft.rfftn(g)
    kx = np.fft.fftfreq(n) * 2
    kz = np.fft.rfftfreq(n) * 2
    k2 = kx[:, None, None] ** 2 + kx[None, :, None] ** 2 + kz[None, None, :] ** 2
    k2[0, 0, 0] = 1.0
    filt = np.sqrt(np.exp(-k2 / kcut**2) / k2)
    filt[0, 0, 0] = 0.0
    g = np.fft.irfftn(gk * filt, s=(n, n, n), axes=(0, 1, 2))
    g /= g.std()
    d = np.exp(sigma * g)
    d /= d.mean()
    d = np.maximum(d, 1e-3)
    # array indices are (i, j, k); written with order='F' so i is fastest
    return d


def write_problem(path, n, gas, sources, z_init, da_sim, nsteps, note, nray_max_factor=None):
    os.makedirs(path, exist_ok=True)
    # Fortran order, little-endian float32: component index fastest
    np.asfortranarray(gas).astype("<f4").ravel(order="F").tofile(os.path.join(path, "gas.bin"))
    with open(os.path.join(path, "sources.txt"), "w") as f:
        for (i, j, k, luv, sfr) in sources:
            f.write(f"{i:d} {j:d} {k:d} {luv:.6e} {sfr:.6e}\n")
    with open(os.path.join(path, "params.txt"), "w") as f:
        f.write(f"N = {n}\n")
        f.write(f"z_init = {z_init}\n")
        f.write(f"da_sim = {da_sim}\n")
        f.write(f"nsteps = {nsteps}\n")
        if nray_max_factor is not None:
            f.write(f"nray_max_factor = {nray_max_factor}\n")
    with open(os.path.join(path, "README.txt"), "w") as f:
        f.write(note)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("names", nargs="*",
                    default=["uniform_1src", "uniform_8src", "lognormal_16src", "uniform_8src_smallbuf",
                             "lognormal_16src_smallbuf"])
    ap.add_argument("--N", type=int, default=32)
    ap.add_argument("--z_init", type=float, default=10.0)
    ap.add_argument("--da_sim", type=float, default=1e-4)
    ap.add_argument("--nsteps", type=int, default=10)
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--outdir", default=os.path.normpath(os.path.join(here, "..", "..", "tests", "problems")))
    ap.add_argument("--force", action="store_true", help="overwrite existing problems")
    args = ap.parse_args()
    n = args.N

    for name in args.names:
        path = os.path.join(args.outdir, name)
        if os.path.exists(os.path.join(path, "gas.bin")) and not args.force:
            print(f"{name}: exists, skipped ({path})")
            continue

        nrf = None
        if name == "uniform_1src":
            delta = np.ones((n, n, n))
            c = n // 2
            sources = [(c, c, c, 2.0e28, 1.0)]
            note = "Uniform neutral gas at mean baryon density, T = temp_min, one source.\n"
        elif name in ("uniform_8src", "uniform_8src_smallbuf"):
            delta = np.ones((n, n, n))
            fr = [(0.25, 0.25, 0.25), (0.75, 0.25, 0.50), (0.25, 0.75, 0.75), (0.75, 0.75, 0.25),
                  (0.50, 0.50, 0.90), (0.10, 0.50, 0.50), (0.60, 0.10, 0.70), (0.40, 0.90, 0.10)]
            sources = [(1 + int(a * n), 1 + int(b * n), 1 + int(c * n), 1.0e28, 0.5)
                       for (a, b, c) in fr]
            note = "Uniform neutral gas at mean baryon density, T = temp_min, 8 sources.\n"
            if name == "uniform_8src_smallbuf":
                nrf = 4
                note += ("nray_max_factor = 4 (Nray_max = 4*N^3): the ngam_ray thresholds become\n"
                         "nonzero once the ray count exceeds j*Nray_max/64, exercising merging.\n")
        elif name in ("lognormal_16src", "lognormal_16src_smallbuf"):
            d = lognormal_delta(n)
            # d[i,j,k] with i = x index (Fortran fastest)
            delta = d
            order = np.argsort(d.ravel(order="F"), kind="stable")[::-1][:16]
            sources = []
            for r, idx in enumerate(order):
                i, j, k = np.unravel_index(idx, (n, n, n), order="F")
                luv = 1.0e28 * (1.0 + 0.1 * r)
                sources.append((int(i) + 1, int(j) + 1, int(k) + 1, luv, luv / 2.0e28))
            note = ("Log-normal density field (seed 20121, sigma_lnrho=1.3, P(k)~k^-2 exp(-(k/0.5k_nyq)^2)),\n"
                    "neutral, T = temp_min, 16 sources in the densest cells (densest first),\n"
                    "Luv = 1e28*(1+0.1*rank), sfr = Luv/2e28.\n"
                    f"delta: min {d.min():.4g} max {d.max():.4g}\n")
            if name == "lognormal_16src_smallbuf":
                nrf = 1
                note += ("nray_max_factor = 1 (Nray_max = N^3): ngam_ray thresholds become nonzero\n"
                         "and rays merge from hydro step ~5 (uniform_8src_smallbuf has too few rays\n"
                         "to merge).\n")
        else:
            print(f"unknown problem {name}", file=sys.stderr)
            sys.exit(1)

        gas, temp_min = make_gas(delta, args.z_init, n)
        write_problem(path, n, gas, sources, args.z_init, args.da_sim, args.nsteps, note, nrf)
        print(f"{name}: written to {path} (N={n}, temp_min={temp_min:.6g} K, {len(sources)} sources)")


if __name__ == "__main__":
    main()
