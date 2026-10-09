"""Shared helpers of the hydro accuracy suite (docs/hydro_tests.md).

Problem files are written in the format of tests/hydro_suite/suite_driver.cpp
(= apps/hydro_driver): float32, Fortran order, component index fastest,
supercomoving code units (lengths in gas cells, velocities in cells per code
time unit, densities in units of the mean matter density, so the mean gas
density is Dgas_avg = Omega_b/Omega_m = 1/6).

numpy only.
"""
import json
import math
import os
import subprocess
import time

import numpy as np

F32 = np.float32
GAMMA = float(F32(5.0) / F32(3.0))     # gamma_ideal = 5./3 (real(4) division)
XHY = float(F32(0.76))
YHE = float(F32(0.24))
OMEGAM = float(F32(0.27))
OMEGAB = float(F32(0.045))
DGAS_AVG = OMEGAB / OMEGAM

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))


# ---------------------------------------------------------------- files
def write_f32(path, a):
    """a: array (nc, N, N, N) indexed [c, i, j, k]; written in Fortran order."""
    np.asarray(a, dtype="<f4").ravel(order="F").tofile(path)


def read_f32(path, nc, N):
    a = np.fromfile(path, dtype="<f4")
    return a.reshape((nc, N, N, N), order="F").astype(np.float64)


def write_problem(d, params, rho, vel, eth, vg=None, extra_files=None):
    """Write u_gas.bin (+v_gas.bin) and params.txt; vel is the lab velocity."""
    os.makedirs(d, exist_ok=True)
    if vg is None:
        vg = np.zeros_like(vel)
    else:
        write_f32(os.path.join(d, "v_gas.bin"), vg)
    dv = vel - vg
    u = np.empty((5,) + rho.shape)
    u[0] = rho
    u[1:4] = rho * dv
    u[4] = eth + 0.5 * rho * (dv ** 2).sum(axis=0)
    write_f32(os.path.join(d, "u_gas.bin"), u)
    with open(os.path.join(d, "params.txt"), "w") as f:
        for k, v in params.items():
            if isinstance(v, (list, tuple, np.ndarray)):
                v = " ".join(repr(float(x)) for x in v)
            f.write(f"{k} = {v}\n")
    for name, arr in (extra_files or {}).items():
        np.asarray(arr, dtype="<f4").ravel(order="F").tofile(os.path.join(d, name))


def read_info(d):
    out = {}
    for line in open(os.path.join(d, "info.txt")):
        if "=" in line:
            k, v = line.split("=", 1)
            vals = [float(x) for x in v.split()]
            out[k.strip()] = vals if len(vals) > 1 else vals[0]
    return out


def outputs(run):
    return sorted(os.path.join(run, n) for n in os.listdir(run) if n.startswith("out_"))


def read_state(d, N):
    """rho, lab velocity (3,...), pressure, info of one output directory.
    With hydro_scheme = unsplit the pressure comes from the internal-energy
    field e_gas.bin (the dual-energy temperature, docs/hydro_unsplit.md),
    which the driver writes only for that scheme."""
    u = read_f32(os.path.join(d, "u_gas.bin"), 5, N)
    v = read_f32(os.path.join(d, "v_gas.bin"), 3, N)
    rho = u[0]
    dv = u[1:4] / rho
    vel = dv + v
    ef = os.path.join(d, "e_gas.bin")
    if os.path.exists(ef):
        eth = read_f32(ef, 1, N)[0]
    else:
        eth = u[4] - 0.5 * rho * (dv ** 2).sum(axis=0)
    return rho, vel, (GAMMA - 1) * eth, read_info(d)


def read_log(run):
    a = np.loadtxt(os.path.join(run, "log.txt"), ndmin=2)
    names = ["istep", "t", "a", "z", "dt", "cmax", "amax", "mass", "px", "py", "pz",
             "energy", "thermal", "wall"]
    return {n: a[:, i] for i, n in enumerate(names)}


def driver_path():
    for cand in (os.environ.get("HYDRO_SUITE_DRIVER"),
                 os.path.join(ROOT, "build", "tests", "hydro_suite", "hydro_suite_driver")):
        if cand and os.path.exists(cand):
            return cand
    raise SystemExit("hydro_suite_driver not found (build it or set HYDRO_SUITE_DRIVER)")


def run_driver(pdir, odir, sets=(), driver=None, log=None):
    """Run the driver; returns the wall time [s]. Raises on failure."""
    driver = driver or driver_path()
    cmd = [driver, pdir, odir, "--quiet"]
    for s in sets:
        cmd += ["--set", s]
    t0 = time.time()
    r = subprocess.run(cmd, capture_output=True, text=True)
    wall = time.time() - t0
    if log:
        log(r.stdout.strip())
    if r.returncode != 0:
        raise RuntimeError(f"driver failed on {pdir}: {r.stderr.strip()}")
    return wall


def _compact(o):
    """round floats to 6 significant digits (results files stay small)"""
    if isinstance(o, dict):
        return {k: _compact(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_compact(v) for v in o]
    if isinstance(o, (float, np.floating)):
        x = float(o)
        return float(f"{x:.6g}") if math.isfinite(x) else x
    if isinstance(o, np.integer):
        return int(o)
    return o


def save_json(path, obj):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        json.dump(_compact(obj), f, separators=(",", ":"), sort_keys=True, default=float)


def load_json(path):
    with open(path) as f:
        return json.load(f)


# ---------------------------------------------------------------- Riemann
def riemann_star(rl, ul, pl, rr, ur, pr, g=GAMMA):
    """p*, u* of the exact Riemann problem (Toro ch. 4), Newton iteration."""
    cl, cr = math.sqrt(g * pl / rl), math.sqrt(g * pr / rr)

    def fk(p, rk, pk, ck):
        if p > pk:
            A, B = 2 / ((g + 1) * rk), (g - 1) / (g + 1) * pk
            q = math.sqrt(A / (p + B))
            return (p - pk) * q, q * (1 - (p - pk) / (2 * (B + p)))
        pr_ = p / pk
        return (2 * ck / (g - 1) * (pr_ ** ((g - 1) / (2 * g)) - 1),
                1 / (rk * ck) * pr_ ** (-(g + 1) / (2 * g)))

    # two-rarefaction guess, then Newton (with bisection safeguard)
    p = ((cl + cr - 0.5 * (g - 1) * (ur - ul)) /
         (cl / pl ** ((g - 1) / (2 * g)) + cr / pr ** ((g - 1) / (2 * g)))) ** (2 * g / (g - 1))
    p = max(p, 1e-14 * (pl + pr))
    for _ in range(200):
        f1, d1 = fk(p, rl, pl, cl)
        f2, d2 = fk(p, rr, pr, cr)
        f = f1 + f2 + ur - ul
        pn = p - f / (d1 + d2)
        if pn <= 0:
            pn = 0.1 * p
        if abs(pn - p) < 1e-14 * (pn + p):
            p = pn
            break
        p = pn
    f1, _ = fk(p, rl, pl, cl)
    f2, _ = fk(p, rr, pr, cr)
    u = 0.5 * (ul + ur) + 0.5 * (f2 - f1)
    return p, u


def riemann_waves(rl, ul, pl, rr, ur, pr, g=GAMMA):
    """Star state and wave speeds: dict with ps, us, rsl, rsr and speeds of
    the left wave (head, tail), contact, right wave (tail, head)."""
    ps, us = riemann_star(rl, ul, pl, rr, ur, pr, g)
    cl, cr = math.sqrt(g * pl / rl), math.sqrt(g * pr / rr)
    gm = (g - 1) / (g + 1)
    w = {"ps": ps, "us": us}
    if ps > pl:
        s = ul - cl * math.sqrt((g + 1) / (2 * g) * ps / pl + (g - 1) / (2 * g))
        w.update(left="shock", lhead=s, ltail=s, rsl=rl * (ps / pl + gm) / (gm * ps / pl + 1))
    else:
        csl = cl * (ps / pl) ** ((g - 1) / (2 * g))
        w.update(left="rarefaction", lhead=ul - cl, ltail=us - csl, rsl=rl * (ps / pl) ** (1 / g))
    if ps > pr:
        s = ur + cr * math.sqrt((g + 1) / (2 * g) * ps / pr + (g - 1) / (2 * g))
        w.update(right="shock", rhead=s, rtail=s, rsr=rr * (ps / pr + gm) / (gm * ps / pr + 1))
    else:
        csr = cr * (ps / pr) ** ((g - 1) / (2 * g))
        w.update(right="rarefaction", rhead=ur + cr, rtail=us + csr, rsr=rr * (ps / pr) ** (1 / g))
    return w


def riemann_sample(rl, ul, pl, rr, ur, pr, xi, g=GAMMA):
    """(rho, u, p) of the exact solution at similarity coordinates xi=x/t."""
    w = riemann_waves(rl, ul, pl, rr, ur, pr, g)
    ps, us = w["ps"], w["us"]
    cl, cr = math.sqrt(g * pl / rl), math.sqrt(g * pr / rr)
    xi = np.asarray(xi, dtype=float)
    r = np.empty_like(xi)
    u = np.empty_like(xi)
    p = np.empty_like(xi)
    for n, s in enumerate(xi.flat):
        if s < us:
            if s < w["lhead"]:
                st = (rl, ul, pl)
            elif s >= w["ltail"]:
                st = (w["rsl"], us, ps)
            else:  # inside the left rarefaction fan
                uu = 2 / (g + 1) * (cl + (g - 1) / 2 * ul + s)
                c = 2 / (g + 1) * (cl + (g - 1) / 2 * (ul - s))
                st = (rl * (c / cl) ** (2 / (g - 1)), uu, pl * (c / cl) ** (2 * g / (g - 1)))
        else:
            if s > w["rhead"]:
                st = (rr, ur, pr)
            elif s <= w["rtail"]:
                st = (w["rsr"], us, ps)
            else:
                uu = 2 / (g + 1) * (-cr + (g - 1) / 2 * ur + s)
                c = 2 / (g + 1) * (cr - (g - 1) / 2 * (ur - s))
                st = (rr * (c / cr) ** (2 / (g - 1)), uu, pr * (c / cr) ** (2 * g / (g - 1)))
        r.flat[n], u.flat[n], p.flat[n] = st
    return r, u, p


def cell_average(fun, x_lo, x_hi, nsub=64):
    """Average of fun(x) (returns tuple of arrays) over [x_lo, x_hi] per cell."""
    s = (np.arange(nsub) + 0.5) / nsub
    xs = x_lo[:, None] + (x_hi - x_lo)[:, None] * s[None, :]
    vals = fun(xs.ravel())
    return tuple(v.reshape(xs.shape).mean(axis=1) for v in vals)


# ---------------------------------------------------------------- Sedov
def sedov_profile(g=GAMMA, n=6000):
    """Self-similar Sedov-Taylor solution (spherical, uniform medium), in the
    closed parametric form of Kamm & Timmes (2007, LA-UR-07-2849) for j = 3,
    w = 0 (equivalent to Landau & Lifshitz §106).

    Returns dict with lam = r/R (increasing, 0..1) and the profiles
    normalised to the post-shock values: rho/rho2, v/v2, p/p2 with
    rho2 = rho0 (g+1)/(g-1), v2 = 2/(g+1) Rdot, p2 = 2/(g+1) rho0 Rdot^2;
    xi0 such that R = xi0 (E t^2/rho0)^(1/5) (from the energy integral) and
    mass_check = M(<R)/(4/3 pi R^3 rho0) (should be 1).
    """
    j = 3.0
    a = (j + 2) * (g + 1) / 4
    b = (g + 1) / (g - 1)
    c = (j + 2) * g / 2
    d = (j + 2) * (g + 1) / ((j + 2) * (g + 1) - 2 * (2 + j * (g - 1)))
    e = (2 + j * (g - 1)) / 2
    al0 = 2 / (j + 2)
    al2 = -(g - 1) / (2 * (g - 1) + j)
    al1 = (j + 2) * g / (2 + j * (g - 1)) * (2 * j * (2 - g) / (g * (j + 2) ** 2) - al2)
    al3 = j / (2 * (g - 1) + j)
    al4 = (j + 2) * j / (j * (2 - g)) * al1
    al5 = -2 * j / (j * (2 - g))
    V0, V2 = 2 / ((j + 2) * g), 4 / ((j + 2) * (g + 1))
    w = np.unique(np.concatenate((np.logspace(-15, 0, n), 1 - np.logspace(-15, 0, n),
                                  np.linspace(0, 1, n))))
    w = w[(w > 0) & (w <= 1)]
    V = V0 + (V2 - V0) * w
    x1 = a * V
    x2 = b * (c * V - 1)
    x3 = d * (1 - e * V)
    x4 = b * (1 - c * V / g)
    lam = x1 ** (-al0) * x2 ** (-al2) * x3 ** (-al1)
    vel = x1 * lam
    rho = x2 ** al3 * x3 ** al4 * x4 ** al5
    pres = x1 ** (al0 * j) * x3 ** (al4 - 2 * al1) * x4 ** (1 + al5)
    order = np.argsort(lam)
    lam, rho, vel, pres = lam[order], rho[order], vel[order], pres[order]
    lam = np.concatenate(([0.0], lam))
    rho = np.concatenate(([0.0], rho))
    vel = np.concatenate(([0.0], vel))
    pres = np.concatenate(([pres[0]], pres))
    # energy and mass integrals for R = 1, t = 1, rho0 = 1 (Rdot = 2/5)
    Rdot = 0.4
    dens = rho * (g + 1) / (g - 1)
    v = vel * 2 / (g + 1) * Rdot
    p = pres * 2 / (g + 1) * Rdot ** 2
    en = 0.5 * dens * v ** 2 + p / (g - 1)
    E = np.trapezoid(4 * np.pi * lam ** 2 * en, lam)
    M = np.trapezoid(4 * np.pi * lam ** 2 * dens, lam) / (4 * np.pi / 3)
    return {"lam": lam, "rho": rho, "v": vel, "p": pres, "xi0": E ** (-0.2), "mass_check": M}


def sedov_radius(E, rho0, t, g=GAMMA):
    return sedov_profile(g)["xi0"] * (E * t ** 2 / rho0) ** 0.2
