#!/usr/bin/env python3
"""Problem generators of the hydro accuracy suite (docs/hydro_tests.md).

usage: generate.py <test> <outdir> [key=value ...]
  tests: tube, sedov, sound, pancake, cosmobox, kh

Each generator writes a problem directory for tests/hydro_suite/suite_driver
(params.txt, u_gas.bin, optional v_gas.bin, xv_dm.init) and a meta.json with
everything the analysis needs. Every problem carries a uniform bulk velocity
with three non-zero components where the flow would otherwise have a pencil
of exactly zero velocity (quirk H1 of docs/port_notes_hydro.md: the MF
advect/euler sub-step loop never terminates when the signal speed of a whole
pencil is zero; the C++ throws).
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hs_common as hc  # noqa: E402

# Riemann problems in dimensionless units (rho, u, p), domain [0, 1),
# interfaces at 0.5 (L|R) and 0 (R|L, periodic), and the end time.
TUBES = {
    # Sod (1978)
    "sod": dict(L=(1.0, 0.0, 1.0), R=(0.125, 0.0, 0.1), t_end=0.12),
    # strong shock, pressure ratio 1e5 (Toro test 3 states)
    "strong": dict(L=(1.0, 0.0, 1000.0), R=(1.0, 0.0, 0.01), t_end=0.005),
    # 123 problem (Einfeldt et al. 1991): double rarefaction, near-vacuum
    "123": dict(L=(1.0, -2.0, 0.4), R=(1.0, 2.0, 0.4), t_end=0.07),
}

# direction of the tube normal (grid axes), and the period along the normal
DIRS = {"x": (1, 0, 0), "y": (0, 1, 0), "z": (0, 0, 1), "xy": (1, 1, 0), "xyz": (1, 1, 1)}

DEFAULT_BULK = (0.25, -0.15, 0.1)


def tube_geometry(direction, N):
    """Integer direction m, unit normal n, period along the normal Ln (cells)
    and the cell-centre coordinate s = m.x (x = cell centre, 0-based i+1/2),
    with the sub-cell offset distribution of s (for cell averages)."""
    m = np.array(DIRS[direction], dtype=float)
    n = m / np.linalg.norm(m)
    Ln = N / np.linalg.norm(m)
    c = np.arange(N) + 0.5
    X = np.meshgrid(c, c, c, indexing="ij")
    s = sum(m[a] * X[a] for a in range(3))
    sub = (np.arange(8) + 0.5) / 8 - 0.5
    O = np.meshgrid(sub, sub, sub, indexing="ij")
    off = sum(m[a] * O[a] for a in range(3)).ravel()
    return m, n, Ln, s, off


def tube_scaling(kind):
    """velocity scale U (code velocity per dimensionless velocity): the left
    sound speed is 2 cells per code time unit."""
    rl, ul, pl = TUBES[kind]["L"]
    return 2.0 / math.sqrt(hc.GAMMA * pl / rl)


def tube_states(kind):
    T = TUBES[kind]
    U = tube_scaling(kind)
    conv = lambda st: (st[0] * hc.DGAS_AVG, st[1] * U, st[2] * hc.DGAS_AVG * U * U)  # noqa
    return conv(T["L"]), conv(T["R"]), U


def tube_initial_1d(kind, xi):
    """dimensionless conserved variables at the coordinate xi in [0,1)"""
    T = TUBES[kind]
    xi = np.mod(xi, 1.0)
    left = xi < 0.5
    out = []
    for q in range(3):
        out.append(np.where(left, T["L"][q], T["R"][q]))
    r, u, p = out
    return r, r * u, p / (hc.GAMMA - 1) + 0.5 * r * u * u


def gen_tube(outdir, kind="sod", direction="x", N=32, bulk=None, z_init=1000.0,
             nout=1, t_frac=1.0):
    """Shock tube along `direction`. The dimensionless problem lives on the
    normal coordinate xi = s/N mod 1 (s = m.x), i.e. on a period Ln = N/|m|
    cells; code time t = t_dimless * Ln / U."""
    m, n, Ln, s, off = tube_geometry(direction, N)
    if bulk is None:
        # the default bulk velocity rotates with an axis-aligned tube, so that
        # x, y and z tubes are the same problem up to the sweep order
        bulk = np.roll(DEFAULT_BULK, "xyz".index(direction)) if direction in ("x", "y", "z") \
            else DEFAULT_BULK
    bulk = [float(b) for b in bulk]
    (rL, uL, pL), (rR, uR, pR), U = tube_states(kind)
    t_end = TUBES[kind]["t_end"] * t_frac * Ln / U
    # cell averages of the conserved variables (dimensionless) over the cell
    keys, inv = np.unique(np.mod(s, N), return_inverse=True)
    xs = (keys[:, None] + off[None, :]) / N
    r, mom, en = (a.mean(axis=1) for a in tube_initial_1d(kind, xs.ravel().reshape(xs.shape)))
    r = r[inv].reshape(s.shape)
    mom = mom[inv].reshape(s.shape)
    en = en[inv].reshape(s.shape)
    # to code units
    rho = r * hc.DGAS_AVG
    un = mom / r * U
    etot = en * hc.DGAS_AVG * U * U
    eth = etot - 0.5 * rho * un ** 2
    vel = np.empty((3,) + rho.shape)
    for a in range(3):
        vel[a] = un * n[a] + bulk[a]
    t_out = [t_end * (k + 1) / nout for k in range(nout)]
    params = {"N": N, "mode": "hydro", "z_init": z_init, "t_out": t_out}
    hc.write_problem(outdir, params, rho, vel, eth)
    meta = dict(test="tube", kind=kind, direction=direction, N=N, bulk=list(bulk),
                z_init=z_init, t_out=t_out, U=U, Ln=Ln, normal=list(n), m=list(m),
                L=TUBES[kind]["L"], R=TUBES[kind]["R"])
    hc.save_json(os.path.join(outdir, "meta.json"), meta)
    return meta


# ------------------------------------------------------------------ Sedov
def gen_sedov(outdir, N=64, t_end=10.0, nout=4, R_end_frac=0.35, bulk_frac=(0.11, 0.07, 0.05),
              mach_end=40.0, z_init=1e4, sigma=1.0):
    """Point explosion at the box centre (cell-centre (N/2-1/2)*3, 0-based),
    energy E in a Gaussian of `sigma` cells, uniform background rho0 =
    Dgas_avg. E is chosen so that R(t_end) = R_end_frac*N; the background
    sound speed so that the final shock Mach number is mach_end; the bulk
    velocity is bulk_frac times the final shock speed (H1)."""
    rho0 = hc.DGAS_AVG
    xi0 = hc.sedov_profile()["xi0"]
    R_end = R_end_frac * N
    E = (R_end / xi0) ** 5 * rho0 / t_end ** 2
    Rdot_end = 0.4 * R_end / t_end
    c0 = Rdot_end / mach_end
    p0 = rho0 * c0 ** 2 / hc.GAMMA
    bulk = [b * Rdot_end for b in bulk_frac]
    c = np.arange(N) + 0.5
    X = np.meshgrid(c, c, c, indexing="ij")
    centre = (N // 2) + 0.5   # a cell centre
    r2 = sum((X[a] - centre) ** 2 for a in range(3))
    w = np.exp(-r2 / (2 * sigma ** 2))
    w /= w.sum()
    eth = p0 / (hc.GAMMA - 1) + E * w
    rho = np.full((N, N, N), rho0)
    vel = np.empty((3, N, N, N))
    for a in range(3):
        vel[a] = bulk[a]
    t_out = [t_end * 2.0 ** (k - nout + 1) for k in range(nout)]
    params = {"N": N, "mode": "hydro", "z_init": z_init, "t_out": t_out}
    hc.write_problem(outdir, params, rho, vel, eth)
    meta = dict(test="sedov", N=N, E=E, rho0=rho0, p0=p0, centre=[centre] * 3, bulk=bulk,
                t_out=t_out, sigma=sigma, z_init=z_init, R_end=R_end)
    hc.save_json(os.path.join(outdir, "meta.json"), meta)
    return meta


# ------------------------------------------------------------------ sound
def gen_sound(outdir, N=32, A=1e-4, periods=1.0, nout=4, bulk=(0.2, 0.3, -0.25), c0=2.0,
              direction="x", z_init=1e4):
    """Travelling sound wave, one wavelength per box along `direction`
    (x or xy), exact eigenvector, cell-averaged. Phase speed c0 + vb.n."""
    m, n, Ln, s, off = tube_geometry(direction, N)
    rho0 = hc.DGAS_AVG
    p0 = rho0 * c0 ** 2 / hc.GAMMA
    k = 2 * np.pi
    keys, inv = np.unique(np.mod(s, N), return_inverse=True)
    xs = (keys[:, None] + off[None, :]) / N
    # average of sin over the cell (in xi)
    sn = np.sin(k * xs).mean(axis=1)[inv].reshape(s.shape)
    rho = rho0 * (1 + A * sn)
    un = c0 * A * sn
    p = p0 * (1 + hc.GAMMA * A * sn)
    vel = np.empty((3, N, N, N))
    for a in range(3):
        vel[a] = un * n[a] + bulk[a]
    eth = p / (hc.GAMMA - 1)
    T = Ln / c0
    t_out = [periods * T * (q + 1) / nout for q in range(nout)]
    params = {"N": N, "mode": "hydro", "z_init": z_init, "t_out": t_out}
    hc.write_problem(outdir, params, rho, vel, eth)
    meta = dict(test="sound", N=N, A=A, c0=c0, bulk=list(bulk), direction=direction,
                t_out=t_out, Ln=Ln, m=list(m), normal=list(n), period=T, z_init=z_init)
    hc.save_json(os.path.join(outdir, "meta.json"), meta)
    return meta


# ------------------------------------------------------------------ KH
def gen_kh(outdir, N=64, dv=1.0, a_shear=None, A=1e-3, mode=1, bulk=(0.13, 0.07, 0.11),
           c0=2.0, t_end=None, nout=16, z_init=1e5):
    """Kelvin-Helmholtz: two shear layers vx = +-dv/2 (tanh profile of width
    a_shear cells) at y = N/4 and 3N/4, equal density and pressure, Mach
    dv/c0, seeded by a vy perturbation with x-wavenumber `mode`/N.
    Bulk (0.13, 0.07, 0.11) (H1)."""
    if a_shear is None:
        a_shear = N / 32.0
    rho0 = hc.DGAS_AVG
    p0 = rho0 * c0 ** 2 / hc.GAMMA
    c = np.arange(N) + 0.5
    X, Y, Z = np.meshgrid(c, c, c, indexing="ij")
    vx = 0.5 * dv * (np.tanh((Y - N / 4) / a_shear) - np.tanh((Y - 3 * N / 4) / a_shear) - 1)
    k = 2 * np.pi * mode / N
    env = np.exp(-((Y - N / 4) ** 2) / (2 * (N / 16) ** 2)) + np.exp(-((Y - 3 * N / 4) ** 2) / (2 * (N / 16) ** 2))
    vy = A * dv * np.sin(k * X) * env
    vel = np.stack([vx + bulk[0], vy + bulk[1], np.zeros_like(vx) + bulk[2]])
    rho = np.full((N, N, N), rho0)
    eth = np.full((N, N, N), p0 / (hc.GAMMA - 1))
    if t_end is None:
        t_end = 8.0 / (k * dv / 2)   # ~8 e-foldings of the sharp-interface rate
    t_out = [t_end * (q + 1) / nout for q in range(nout)]
    params = {"N": N, "mode": "hydro", "z_init": z_init, "t_out": t_out}
    hc.write_problem(outdir, params, rho, vel, eth)
    meta = dict(test="kh", N=N, dv=dv, a_shear=a_shear, A=A, mode=mode, k=k, c0=c0,
                bulk=list(bulk), t_out=t_out, z_init=z_init)
    hc.save_json(os.path.join(outdir, "meta.json"), meta)
    return meta


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    test, outdir = sys.argv[1], sys.argv[2]
    kw = {}
    for a in sys.argv[3:]:
        k, v = a.split("=", 1)
        try:
            kw[k] = eval(v, {})  # numbers, tuples
        except Exception:
            kw[k] = v
    fn = {"tube": gen_tube, "sedov": gen_sedov, "sound": gen_sound, "kh": gen_kh}
    if test in ("pancake", "cosmobox"):
        import cosmo_gen
        fn = {"pancake": cosmo_gen.gen_pancake, "cosmobox": cosmo_gen.gen_cosmobox}
    print(fn[test](outdir, **kw))
    return 0


if __name__ == "__main__":
    sys.exit(main())
