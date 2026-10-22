#!/usr/bin/env python3
"""Run-level tests of the runtime cosmology and the preheating floor
(params.txt omegam, omegab, hubble0, sigma8, ns, T_min, z_floor;
docs/production_run.md section 10).

usage: cosmology_check.py MODE WORKDIR [args]

  ics MAKE_ICS                   make_ics.py: the original's values written
                                 explicitly give byte-identical ICs; a Planck
                                 cosmology gives other ICs and records itself
                                 in ics_info.txt.
  rt RT_DRIVER PROBLEM           rt_driver default against Planck: the
                                 redshift schedule (fixed da) is the same to
                                 round-off,
                                 the clock, the RT substep count and the units
                                 in rt_stats.txt follow independent closed
                                 forms of each cosmology, the fields differ.
  full RADHYDRO PROBLEM MAKE_ICS [radhydro args]
                                 radhydro with Planck ICs: refused without
                                 the keys or with another z_init (the
                                 ics_info.txt check), runs with them; the
                                 first step's redshift is the original's and
                                 the year at z_init is Planck's.
  floor RT_DRIVER [--amrex LAUNCHER FLAGS NPFLAG NP]
                                 a cooling test: gas of a range of densities,
                                 temperatures and ionisation states without
                                 sources. T_min = 0, and a floor at redshifts
                                 the run never reaches, are bit-identical to
                                 the default; T_min = 1000 K holds every cell
                                 at >= 1000 K, stops the cells that cool
                                 through 1000 K there, lifts the cold neutral
                                 cells the RT physics otherwise skips, and
                                 leaves the cells that stay warmer bit-identical;
                                 the same with rt_method = m1. --amrex: the
                                 floored run on rt_backend = amrex (NP ranks,
                                 FLAGS ;-separated) is bit-identical to the
                                 serial one (legacy_rt = 0; no sources, so
                                 no ray order enters).

Exit status 0 = pass.
"""
import glob
import math
import os
import shutil
import subprocess
import sys

import numpy as np

f32 = lambda x: float(np.float32(x))  # noqa: E731

# the original (cosmo_parameters.f90) and Planck 2018 (TT,TE,EE+lowE+
# lensing+BAO, from memory: docs/production_run.md section 10; only a
# cosmology different from the original matters here)
ORIG = {"omegam": 0.27, "omegab": 0.045, "hubble0": 0.70, "sigma8": 0.80, "ns": 0.96}
PLANCK = {"omegam": 0.3111, "omegab": 0.04897, "hubble0": 0.6766, "sigma8": 0.8102,
          "ns": 0.9665}

# physical_constants.f90
C_CGS, H0_CGS, MPC2CM, RHOC_CGS, MH_CGS = 2.99792e10, 3.24086e-18, 3.08560e24, 1.87890e-29, 1.67223e-24
SEC2YR = 1.0 / (365 * 86400)
XHY, TCMB0 = f32(0.76), f32(2.725)

FAIL = []


def check(cond, msg):
    if not cond:
        FAIL.append(msg)
        print("FAIL:", msg)


def close(a, b, rtol, msg):
    ok = abs(a - b) <= rtol * max(abs(a), abs(b))
    check(ok, f"{msg}: {a!r} vs {b!r} (rtol {rtol})")


def sets(c):
    return [x for k, v in c.items() for x in ("--set", f"{k}={v}")]


def run(cmd, log, expect_ok=True):
    with open(log, "w") as f:
        r = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
    text = open(log).read()
    if expect_ok and r.returncode != 0:
        print(text[-3000:])
        sys.exit(f"command failed: {' '.join(cmd)}")
    return r.returncode, text


def kv(path):
    d = {}
    for line in open(path):
        if "=" in line:
            k, v = (t.strip() for t in line.split("=", 1))
            d[k] = v
    return d


class Cosmo:
    """Independent closed forms of the code's cosmology (float32 inputs)."""

    def __init__(self, c):
        self.om, self.ob, self.h = f32(c["omegam"]), f32(c["omegab"]), f32(c["hubble0"])
        self.ol = 1 - self.om
        self.orad = f32(4.48e-7) * (1 + f32(0.69)) * TCMB0 ** 4 / self.h ** 2
        self.x = 1e3 * (self.ob * self.h ** 2) ** f32(0.4)

    def E(self, a):
        return math.sqrt(self.om / a ** 3 + self.ol + self.orad / a ** 4)

    def year0(self, a):  # init_variables: matter + Lambda, no radiation
        return (2 / (3 * H0_CGS * self.h * math.sqrt(self.ol)) *
                math.asinh(math.sqrt(self.ol / self.om * a ** 3)) * SEC2YR)

    def dyr(self, a1, a2, n=2000):  # Simpson of da/(a H)
        s = 0.0
        for i in range(n + 1):
            a = a1 + (a2 - a1) * i / n
            w = 1 if i in (0, n) else (4 if i % 2 else 2)
            s += w / (a * self.E(a))
        return s * (a2 - a1) / n / 3 / (H0_CGS * self.h) * SEC2YR

    def nH_avg(self, a):
        return RHOC_CGS * self.ob * self.h ** 2 / a ** 3 * XHY / MH_CGS

    def light_cells(self, a1, da, box, N):  # calc_timestep: dl_c/dl_rt
        dl_c = 2 * C_CGS / (H0_CGS * self.h * math.sqrt(self.om)) * (math.sqrt(a1 + da) - math.sqrt(a1))
        return dl_c / (box / self.h * MPC2CM / N)

    def temp_min(self, a):
        return TCMB0 if 1 / a > self.x else max(TCMB0 / self.x / a ** 2, 1.0)


# ---------------------------------------------------------------------------
def mode_ics(work, make_ics):
    base = ["--kind", "zeldovich", "--Ngas", "16", "--Ndm", "8", "--z-init", "100", "--box", "4",
            "--seed", "7"]
    dirs = {}
    for name, extra in (("default", []),
                        ("explicit", [x for k, v in ORIG.items() for x in ("--" + k, str(v))]),
                        ("planck", [x for k, v in PLANCK.items() for x in ("--" + k, str(v))])):
        d = os.path.join(work, name)
        run([sys.executable, make_ics, d] + base + extra, d + ".log")
        dirs[name] = d
    for f in ("xv_dm.init", "u_gas.init"):
        a = open(os.path.join(dirs["default"], f), "rb").read()
        check(a == open(os.path.join(dirs["explicit"], f), "rb").read(),
              f"ics: {f} with the original's values written explicitly differs from the default")
        check(a != open(os.path.join(dirs["planck"], f), "rb").read(), f"ics: Planck {f} equals the default")
    for name, c in (("default", ORIG), ("planck", PLANCK)):
        info = kv(os.path.join(dirs[name], "ics_info.txt"))
        for k, v in c.items():
            check(k in info and float(info[k]) == f32(v), f"ics: {name} ics_info.txt {k} = {info.get(k)}")
    # the white-noise draws are the seed's: the Planck density field is the
    # same random field with another P(k), so the two are strongly correlated
    ug = [np.fromfile(os.path.join(dirs[n], "u_gas.init"), "<f4").reshape(-1, 5)[:, 0]
          for n in ("default", "planck")]
    r = np.corrcoef(ug[0], ug[1])[0, 1]
    check(r > 0.95, f"ics: default and Planck gas densities correlate only {r:.3f}")


# ---------------------------------------------------------------------------
def mode_rt(work, drv, prob):
    p = kv(os.path.join(prob, "params.txt"))
    N, z0, da, nsteps = int(p["N"]), float(p["z_init"]), float(p["da_sim"]), int(p["nsteps"])
    box = 5.0  # rt_driver: the original's box
    out = {}
    for name, c in (("default", None), ("planck", PLANCK)):
        o = os.path.join(work, name)
        shutil.rmtree(o, ignore_errors=True)
        run([drv, prob, o] + (sets(c) if c else []), o + ".log")
        out[name] = o
    steps = sorted(os.path.basename(d) for d in glob.glob(os.path.join(out["default"], "step_*")))
    check(len(steps) == nsteps, f"rt: {len(steps)} step directories, expected {nsteps}")
    for name, c in (("default", ORIG), ("planck", PLANCK)):
        cm = Cosmo(c)
        a = 1 / (1 + z0)
        year = cm.year0(a)
        for st in steps:
            s = kv(os.path.join(out[name], st, "stats.txt"))
            r = kv(os.path.join(out[name], st, "rt_stats.txt"))
            ns = int(s["nstep_rt"])
            check(ns == max(round(cm.light_cells(a, da, box, N)), 1),
                  f"rt {name} {st}: nstep_rt {ns} vs {cm.light_cells(a, da, box, N):.3f} light cells")
            year += cm.dyr(a, a + da)
            a += da
            close(float(r["year"]), year, 1e-6, f"rt {name} {st}: year")
            close(float(r["nH_avg"]), cm.nH_avg(a), 1e-12, f"rt {name} {st}: nH_avg")
            close(float(r["vol_rt"]), (a * box / cm.h / N * MPC2CM) ** 3, 1e-12, f"rt {name} {st}: vol_rt")
            close(float(r["temp_min"]), cm.temp_min(a), 1e-12, f"rt {name} {st}: temp_min")
    # the redshift schedule depends on da only: equal to round-off (the
    # original advances the clock by da/nstep_rt per RT substep, and
    # nstep_rt depends on the cosmology)
    for st in steps:
        for k in ("redshift", "scalefactor"):
            a0 = float(kv(os.path.join(out["default"], st, "stats.txt"))[k])
            a1 = float(kv(os.path.join(out["planck"], st, "stats.txt"))[k])
            close(a0, a1, 1e-13, f"rt {st}: {k} (default vs Planck)")
    last = steps[-1]
    for f in ("fHI", "temp"):
        x0 = np.fromfile(os.path.join(out["default"], last, f + ".bin"), "<f4")
        x1 = np.fromfile(os.path.join(out["planck"], last, f + ".bin"), "<f4")
        check(not np.array_equal(x0, x1), f"rt: the Planck {f} field equals the default")


# ---------------------------------------------------------------------------
def mode_full(work, rh, prob, make_ics, *extra_args):
    p = kv(os.path.join(prob, "params.txt"))
    z0 = float(p["z_init"])
    ics = os.path.join(work, "ics_planck")
    info = kv(os.path.join(prob, p["ics"], "ics_info.txt"))
    cmd = [sys.executable, make_ics, ics, "--kind", "zeldovich", "--Ngas", info["Ngas"],
           "--Ndm", info["Ndm"], "--z-init", info["z_init"], "--box", info["box"],
           "--seed", info["seed"]]
    run(cmd + [x for k, v in PLANCK.items() for x in ("--" + k, str(v))], ics + ".log")
    common = ["--quiet", "--set", f"ics={ics}"] + list(extra_args)
    # refused: the ICs are Planck's, the run the original's
    rc, text = run([rh, prob, os.path.join(work, "nokeys")] + common, os.path.join(work, "nokeys.log"),
                   expect_ok=False)
    check(rc != 0 and "ics_info.txt" in text and "omegam" in text,
          "full: a run without the cosmology keys on Planck ICs was not refused")
    # refused: another z_init
    rc, text = run([rh, prob, os.path.join(work, "zinit")] + common + sets(PLANCK) +
                   ["--set", f"z_init={z0 + 1:g}"], os.path.join(work, "zinit.log"), expect_ok=False)
    check(rc != 0 and "z_init" in text, "full: a run with another z_init than the ICs was not refused")
    # runs: Planck with its keys, and the default for comparison (no --quiet:
    # the calc_timestep summary has the year and the redshifts of each step)
    texts = {}
    for name, extra in (("planck", ["--set", f"ics={ics}"] + sets(PLANCK)), ("default", [])):
        _, texts[name] = run([rh, prob, os.path.join(work, name)] + extra + list(extra_args),
                             os.path.join(work, name + ".log"))

    def first_step(text):
        lines = text.splitlines()
        i = next(i for i, l in enumerate(lines) if l.strip().startswith("Step number"))
        red = [float(x) for x in lines[i + 1].split(":")[1].split()]
        yr = [float(x) for x in lines[i + 3].split(":")[1].split()]
        return red, yr

    r0, y0 = first_step(texts["default"])
    r1, y1 = first_step(texts["planck"])
    check(r0 == r1, f"full: the first step's redshifts differ: {r0} vs {r1} (da = a ra_max/(1-ra_max))")
    a = 1 / (1 + z0)
    for name, c, y in (("default", ORIG, y0), ("planck", PLANCK, y1)):
        close(y[0], Cosmo(c).year0(a), 2e-7, f"full {name}: year at z_init")  # printed as real(4)
    check(y0[0] != y1[0], "full: the year at z_init is the same in both cosmologies")
    check("cosmology omegam 0.3111" in texts["planck"], "full: the Planck run does not report its cosmology")
    check("cosmology" not in texts["default"], "full: the default run reports a cosmology line")


# ---------------------------------------------------------------------------
def write_cooling_problem(d, N=16, seed=3):
    """No sources; overdensities 0.3-100, T 100 K - 30000 K, fHI 0.02-1."""
    os.makedirs(d, exist_ok=True)
    rng = np.random.default_rng(seed)
    n = N ** 3
    c = Cosmo(ORIG)
    z0 = 10.0
    a = 1 / (1 + z0)
    fb = c.ob / c.om
    delta = 10 ** rng.uniform(np.log10(0.3), 2, n)
    T = 10 ** rng.uniform(2, np.log10(3e4), n)
    fHI = np.where(rng.random(n) < 0.4, 1.0, 10 ** rng.uniform(np.log10(0.02), 0, n))
    rho = fb * delta                                      # units of the mean matter density
    # thermal energy per volume in grid units: n_eff k T/(gamma-1)/eden_unit
    x_unit = a * 5.0 / c.h / N * MPC2CM
    t_unit = 2 / (3 * H0_CGS * c.h * math.sqrt(c.om)) * a ** 2
    rho_unit = RHOC_CGS * c.om * c.h ** 2 / a ** 3
    eden = rho_unit * (x_unit / t_unit) ** 2
    rhog = rho * rho_unit
    nH, nHe = rhog * XHY / MH_CGS, rhog * f32(0.24) / (4 * MH_CGS)
    ne = nH * (1 - fHI)
    eth = (nH + nHe + ne) * 1.38066e-16 * T / (f32(5.0 / 3.0) - 1) / eden
    g = np.zeros((n, 9), "<f4")
    g[:, 0] = rho
    g[:, 4] = eth
    g[:, 5] = XHY * rho * fHI                             # w = XHy rho fHI, YHe rho (neutral He)
    g[:, 6] = f32(0.24) * rho
    g.tofile(os.path.join(d, "gas.bin"))
    open(os.path.join(d, "sources.txt"), "w").write("# no sources: a cooling test\n")
    open(os.path.join(d, "params.txt"), "w").write(
        f"N = {N}\nz_init = {z0}\nda_sim = 0.004\nnsteps = 8\nnray_max_factor = 1\n")
    return z0, 1 / (a + 8 * 0.004) - 1


def mode_floor(work, drv, *amrex):
    prob = os.path.join(work, "problem")
    z0, z1 = write_cooling_problem(prob)
    T_MIN = 1000.0
    runs = {"default": [], "zero": ["--set", "T_min=0", "--set", "z_floor=20"],
            "late": ["--set", f"T_min={T_MIN:g}", "--set", f"z_floor={z1 - 0.5:g}"],
            "floor": ["--set", f"T_min={T_MIN:g}", "--set", f"z_floor={z0 + 5:g}"]}
    out = {}
    for name, extra in runs.items():
        o = os.path.join(work, name)
        shutil.rmtree(o, ignore_errors=True)
        run([drv, prob, o] + extra, o + ".log")
        out[name] = o
    steps = sorted(os.path.basename(d) for d in glob.glob(os.path.join(out["default"], "step_*")))
    files = sorted(os.path.relpath(f, out["default"])
                   for f in glob.glob(os.path.join(out["default"], "step_*", "*")))
    # T_min = 0 and a floor below the run's redshifts: bit-identical
    for name in ("zero", "late"):
        bad = [f for f in files
               if open(os.path.join(out["default"], f), "rb").read() != open(os.path.join(out[name], f), "rb").read()]
        check(not bad, f"floor: {name} differs from the default in {bad[:5]}")
    check("T_min" in open(out["floor"] + ".log").read(), "floor: rt_driver does not report the floor")

    def field(name, st, f, nc=1):
        return np.fromfile(os.path.join(out[name], st, f + ".bin"), "<f4").reshape(-1, nc)

    T0 = np.stack([field("default", st, "temp")[:, 0] for st in steps])   # (steps, cells)
    T1 = np.stack([field("floor", st, "temp")[:, 0] for st in steps])
    fHI0 = field("default", steps[0], "fHI")[:, 0]
    # The floor holds where photo_cell applies it, once per RT substep; the
    # outputs come after the substep's clock update, over which the grid-unit
    # energy cools adiabatically (T ~ a^-2): by 2 da/(a nstep_rt) at most.
    st0 = [kv(os.path.join(out["floor"], st, "stats.txt")) for st in steps]
    tol = max(2 * 0.004 / (float(x["scalefactor"]) * int(x["nstep_rt"])) for x in st0) * 1.05 + 1e-5
    check(T1.min() >= T_MIN * (1 - tol), f"floor: a temperature below T_min (1 - {tol:.2e}): {T1.min()}")
    # the cooling test: cells above T_min at the first output that cool below it
    cooled = (T0[0] > 1.05 * T_MIN) & (T0[-1] < 0.9 * T_MIN)
    check(cooled.sum() >= 20, f"floor: only {cooled.sum()} cells cool through T_min without the floor")
    check(np.all(np.abs(T1[-1, cooled] / T_MIN - 1) <= tol), "floor: the cooled cells do not end at T_min")
    # cold neutral cells: the RT physics skips them (no heating or cooling;
    # without the floor they only follow the expansion, T ~ a^-2), with the
    # floor they are lifted to T_min
    cold = (T0 < 0.9 * T_MIN).all(axis=0) & (fHI0 >= 0.99)
    check(cold.sum() >= 20, f"floor: only {cold.sum()} cold neutral cells")
    check(np.all(np.abs(T1[:, cold] / T_MIN - 1) <= tol),
          "floor: the cold neutral cells are not at T_min at every step")
    a = np.array([float(x["scalefactor"]) for x in st0])
    adiabatic = T0[0, cold][None, :] * (a[0] / a[:, None]) ** 2
    check(np.allclose(T0[:, cold], adiabatic, rtol=1e-5),
          "floor: without the floor the cold neutral cells do not just expand (T ~ a^-2)")
    # cells that stay warmer than T_min are bit-identical (no sources: the
    # cells are independent)
    warm = (T0 > 1.001 * T_MIN).all(axis=0)
    check(warm.sum() >= 100, f"floor: only {warm.sum()} warm cells")
    for st in steps:
        for f, nc in (("temp", 1), ("fHI", 1), ("fHeI", 1), ("fHeII", 1), ("u_gas", 5), ("w_gas", 4)):
            a0, a1 = field("default", st, f, nc), field("floor", st, f, nc)
            if f in ("u_gas", "w_gas"):
                a0, a1 = a0.reshape(-1, nc), a1.reshape(-1, nc)
            check(np.array_equal(a0[warm], a1[warm]), f"floor {st}: {f} differs in cells warmer than T_min")
    # the moment RT (rt_method = m1) goes through the same cell physics
    for name, extra in (("m1", []), ("m1_floor", runs["floor"])):
        o = os.path.join(work, name)
        shutil.rmtree(o, ignore_errors=True)
        run([drv, prob, o, "--set", "rt_method=m1"] + extra, o + ".log")
        out[name] = o
    Tm0 = np.stack([field("m1", st, "temp")[:, 0] for st in steps])
    Tm1 = np.stack([field("m1_floor", st, "temp")[:, 0] for st in steps])
    check(Tm1.min() >= T_MIN * (1 - tol), f"floor m1: a temperature below T_min: {Tm1.min()}")
    wm = (Tm0 > 1.001 * T_MIN).all(axis=0)
    check(wm.sum() >= 100 and np.array_equal(Tm0[:, wm], Tm1[:, wm]),
          "floor m1: the cells warmer than T_min differ")
    check(((Tm0[0] > 1.05 * T_MIN) & (Tm0[-1] < 0.9 * T_MIN)).sum() >= 20, "floor m1: no cooling test")
    # the AMReX backend against the serial RT, both legacy_rt = 0, floored
    if amrex:
        if amrex[0] != "--amrex" or len(amrex) != 5:
            sys.exit("floor: --amrex LAUNCHER FLAGS NPFLAG NP")
        launcher, flags, npflag, np_ = amrex[1:]
        launch = [launcher] + [f for f in flags.split(";") if f] + [npflag, np_]
        for name, extra in (("serial0", []), ("amrex", ["--set", "rt_backend=amrex",
                                                        "--set", "amrex_max_grid_size=8"])):
            o = os.path.join(work, name)
            shutil.rmtree(o, ignore_errors=True)
            run(launch + [drv, prob, o, "--set", "legacy_rt=0"] + runs["floor"] + extra, o + ".log")
            out[name] = o
        for st in steps:
            for f in ("temp", "fHI", "fHeI", "fHeII", "u_gas", "w_gas"):
                a0 = open(os.path.join(out["serial0"], st, f + ".bin"), "rb").read()
                a1 = open(os.path.join(out["amrex"], st, f + ".bin"), "rb").read()
                check(a0 == a1, f"floor amrex {st}: {f} differs from the serial RT")
        Ta = field("amrex", steps[-1], "temp")[:, 0]
        check(Ta.min() >= T_MIN * (1 - tol), f"floor amrex: a temperature below T_min: {Ta.min()}")
    print(f"floor: {cooled.sum()} cells cooled through {T_MIN:g} K without the floor and stop at it "
          f"with it; {cold.sum()} cold neutral cells lifted; {warm.sum()} warm cells bit-identical")


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    mode, work = sys.argv[1], sys.argv[2]
    os.makedirs(work, exist_ok=True)
    a = sys.argv[3:]
    {"ics": lambda: mode_ics(work, *a), "rt": lambda: mode_rt(work, *a),
     "full": lambda: mode_full(work, *a), "floor": lambda: mode_floor(work, *a)}[mode]()
    if FAIL:
        sys.exit(f"{len(FAIL)} failure(s)")
    print("OK")


if __name__ == "__main__":
    main()
