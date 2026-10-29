#!/usr/bin/env python3
"""Runner of the hydro accuracy suite (docs/hydro_tests.md).

usage: run_suite.py <test> [--tier quick|full] [--scheme legacy|unsplit]
                    [--amr L] [--refine-mass M] [--rundir DIR] [--check]
                    [--jobs n] [--driver PATH] [--results DIR] [--reuse] [--clean]

  test   tubes | galilean | sedov | sound | pancake | cosmobox | kh | all
  tier   quick: the small cases used by CTest (label hydro_suite);
         full:  everything (label nightly, and the numbers in the docs)
  scheme legacy (default): the original moving-frame TVD scheme;
         unsplit: the unsplit MUSCL-Hancock + HLLC scheme with dual energy
         (docs/hydro_unsplit.md; the driver gets --set hydro_scheme=unsplit)
  --amr L   run on the AMReX hydro backend with L refinement levels
            (hydro_backend=amrex, amrex_max_level=L, amrex_refine_mass=M,
            docs/stage_f.md; needs a driver built with RADRAY_AMREX and the
            unsplit scheme); the results carry the suffix _amrL and the
            analysis reads the level-0 (averaged-down) fields
  --check   apply the pass/fail tolerances (exit status 1 on failure); they
            are per scheme (tolerances.py)
  --reuse   do not re-run cases whose output already exists
  --clean   remove the problem and run directories of each test once it is
            analysed (the full tier writes ~7 GB per scheme; CTest nightly)
  --set K=V an extra driver key for every run (repeatable; e.g.
            dual_e_self=1). The results keep their names: use --results
Results go to <results>/<test>_<tier>[_unsplit].json (default
tests/hydro_suite/results for the full tier, <rundir>/results for quick);
runs to <rundir>/<case>. Figures: analyze.py figures [... scheme].
"""
import argparse
import json
import os
import shutil
import sys
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hs_common as hc  # noqa: E402
import generate as gen  # noqa: E402
import analyze as an  # noqa: E402
import tolerances as tol  # noqa: E402


class Ctx:
    def __init__(self, args):
        self.args = args
        self.scheme = args.scheme
        self.rundir = os.path.abspath(args.rundir)
        self.driver = args.driver or hc.driver_path()
        self.fails = []
        self.made = []
        os.makedirs(self.rundir, exist_ok=True)

    def sets(self, extra=None):
        """driver --set list: the scheme (legacy runs stay exactly as before)"""
        out = list(extra or ())
        if self.scheme != "legacy":
            out.append(f"hydro_scheme={self.scheme}")
        if self.args.amr > 0:
            out += ["hydro_backend=amrex", f"amrex_max_level={self.args.amr}",
                    f"amrex_refine_mass={self.args.refine_mass:g}"]
        return out + list(self.args.set)

    def case(self, name, genfun, **kw):
        """generate (always, cheap) and return (problem dir, run dir)"""
        pdir = os.path.join(self.rundir, "problems", name)
        odir = os.path.join(self.rundir, name)
        meta = genfun(pdir, **kw)
        self.made += [pdir, odir]
        return pdir, odir, meta

    def clean(self):
        """--clean: remove the problem and run directories made so far"""
        for d in self.made:
            shutil.rmtree(d, ignore_errors=True)
        self.made = []

    def run_many(self, cases, sets=None):
        """cases: list of (name, pdir, odir); runs in parallel; returns walls"""
        def one(c):
            name, pdir, odir = c
            done = os.path.exists(os.path.join(odir, "timing.txt"))
            if self.args.reuse and done:
                return name, read_wall(odir)
            print(f"  running {name} ...", flush=True)
            w = hc.run_driver(pdir, odir, self.sets(sets), driver=self.driver)
            print(f"  {name}: {w:.1f} s", flush=True)
            return name, w
        with ThreadPoolExecutor(max_workers=self.args.jobs) as ex:
            return dict(ex.map(one, cases))

    def check(self, ok, msg):
        if not ok:
            self.fails.append(msg)
            print("  FAIL: " + msg)


def read_wall(odir):
    for line in open(os.path.join(odir, "timing.txt")):
        if line.startswith("wall ="):
            return float(line.split("=")[1])
    return float("nan")


# ----------------------------------------------------------------------
def test_tubes(ctx, tier):
    Ns = (32, 64) if tier == "quick" else (32, 64, 128)
    specs = []
    for kind in ("sod", "strong", "123"):
        dirs = ("x", "y", "z", "xy", "xyz") if kind == "sod" else ("x",)
        for d in dirs:
            for N in Ns:
                if tier == "quick" and kind == "sod" and d in ("xy", "xyz") and N > 32:
                    continue
                specs.append((kind, d, N))
    cases, metas = [], {}
    for kind, d, N in specs:
        name = f"tube_{kind}_{d}_N{N}"
        pdir, odir, meta = ctx.case(name, gen.gen_tube, kind=kind, direction=d, N=N)
        cases.append((name, pdir, odir))
        metas[name] = meta
    walls = ctx.run_many(cases)
    res = {}
    for name, pdir, odir in cases:
        r = an.analyze_tube(odir, metas[name])
        r["wall"] = walls[name]
        r["kind"], r["direction"] = metas[name]["kind"], metas[name]["direction"]
        res[name] = r
        print(f"  {name:22s} t={r['t']:.3f} steps={r['nsteps']:3d} L1 rho {r['l1']['rho']:.4f} "
              f"u {r['l1']['u']:.4f} p {r['l1']['p']:.4f}  "
              + " ".join(f"{k}:{v['err_cells']:+.2f}" for k, v in r["positions"].items())
              + f"  min rho {r['min_rho']:.3g} p {r['min_p']:.3g}  {r['wall']:.1f}s")
    # convergence orders
    orders = {}
    for kind in ("sod", "strong", "123"):
        for d in ("x", "y", "z", "xy", "xyz"):
            keys = [f"tube_{kind}_{d}_N{N}" for N in Ns if f"tube_{kind}_{d}_N{N}" in res]
            if len(keys) >= 2:
                NN = [res[k]["N"] for k in keys]
                orders[f"{kind}_{d}"] = {q: an.fit_order(NN, [res[k]["l1"][q] for k in keys])
                                         for q in ("rho", "u", "p")}
    # orientation: x vs y vs z (same N): H2
    orient = {}
    for N in Ns:
        a = res.get(f"tube_sod_x_N{N}")
        for d in ("y", "z", "xy", "xyz"):
            b = res.get(f"tube_sod_{d}_N{N}")
            if a and b and abs(a["t"] - b["t"]) < 1e-9 * a["t"] and d in ("y", "z"):
                orient[f"{d}_vs_x_N{N}"] = an.compare_profiles(a, b)
    out = {"cases": res, "orders": orders, "orientation": orient}
    tol.check_tubes(ctx, out, tier)
    return out


def test_galilean(ctx, tier):
    """Sod along x with different bulk velocities; in the comoving frame
    the result should not depend on the bulk velocity."""
    N = 32 if tier == "quick" else 64
    meta0 = gen.gen_tube(os.path.join(ctx.rundir, "problems", "_tmp"), kind="sod",
                         direction="x", N=N)
    t_end = meta0["t_out"][-1]
    # shifts of the profile by the bulk motion: 1/64, 1/4, 1/2 and 1 box
    shifts = (N / 64, N / 4, N / 2, float(N))
    cases, metas = [], {}
    for sh in shifts:
        for vt in (0.1, 2.0):
            vb = (sh / t_end, -1.5 * vt, vt)
            name = f"gal_N{N}_shift{sh:g}_vt{vt:g}"
            pdir, odir, meta = ctx.case(name, gen.gen_tube, kind="sod", direction="x", N=N,
                                        bulk=vb)
            cases.append((name, pdir, odir))
            metas[name] = meta
    walls = ctx.run_many(cases)
    res = {}
    for name, pdir, odir in cases:
        r = an.analyze_tube(odir, metas[name])
        r["wall"] = walls[name]
        r["bulk"] = metas[name]["bulk"]
        res[name] = r
    ref = res[cases[0][0]]
    for name in res:
        res[name]["diff_vs_ref"] = an.compare_profiles(ref, res[name])
        r = res[name]
        print(f"  {name:28s} vb=({r['bulk'][0]:6.3f},{r['bulk'][1]:5.2f},{r['bulk'][2]:4.2f}) "
              f"steps {r['nsteps']:3d} L1 rho {r['l1']['rho']:.4f}  max|drho| vs ref "
              f"{r['diff_vs_ref']['rho']:.2e}  {r['wall']:.1f}s")
    # H1 demonstration: zero bulk velocity along the tube (the legacy scheme
    # stops with the H1 error; the lab-frame unsplit scheme runs it, and its
    # result is compared with the moving reference)
    name = f"gal_N{N}_zero_bulk_x"
    pdir, odir, meta = ctx.case(name, gen.gen_tube, kind="sod", direction="x", N=N,
                                bulk=(0.0, -0.15, 0.1))
    zero = None
    try:
        hc.run_driver(pdir, odir, ctx.sets(), driver=ctx.driver)
        h1 = "ran" + (" (unexpected)" if ctx.scheme == "legacy" else "")
        zero = an.analyze_tube(odir, meta)
        zero["bulk"] = meta["bulk"]
        zero["diff_vs_ref"] = an.compare_profiles(ref, zero)
        print(f"  {name:28s} vb=(0,-0.15,0.1) steps {zero['nsteps']:3d} L1 rho "
              f"{zero['l1']['rho']:.4f}  max|drho| vs ref {zero['diff_vs_ref']['rho']:.2e}")
    except RuntimeError as e:
        h1 = str(e).split("error:")[-1].strip()
    print(f"  zero bulk vx: {h1}")
    out = {"cases": res, "N": N, "t_end": t_end, "h1_zero_bulk": h1, "zero_bulk": zero}
    tol.check_galilean(ctx, out, tier)
    return out


# ----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("test")
    ap.add_argument("--tier", default="quick", choices=("quick", "full"))
    ap.add_argument("--scheme", default="legacy", choices=("legacy", "unsplit"))
    ap.add_argument("--amr", type=int, default=0)
    ap.add_argument("--refine-mass", type=float, default=8.0)
    ap.add_argument("--rundir", default=None)
    ap.add_argument("--results", default=None)
    ap.add_argument("--driver", default=None)
    ap.add_argument("--jobs", type=int, default=3)
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--reuse", action="store_true")
    ap.add_argument("--clean", action="store_true")
    ap.add_argument("--set", action="append", default=[], metavar="KEY=VALUE")
    args = ap.parse_args()
    suffix = "" if args.scheme == "legacy" else "_" + args.scheme
    if args.amr > 0:
        if args.scheme != "unsplit":
            ap.error("--amr needs --scheme unsplit")
        suffix += f"_amr{args.amr}"
    if args.rundir is None:
        args.rundir = os.path.join(hc.ROOT, "build", "hydro_suite_runs", args.tier + suffix)
    if args.results is None:
        args.results = (os.path.join(hc.HERE, "results") if args.tier == "full"
                        else os.path.join(args.rundir, "results"))
    ctx = Ctx(args)
    tests = {"tubes": test_tubes, "galilean": test_galilean}
    try:
        import tests_extra
        tests.update(tests_extra.TESTS)
    except ImportError:
        pass
    names = list(tests) if args.test == "all" else args.test.split(",")
    for name in names:
        print(f"== {name} ({args.tier}, {args.scheme})", flush=True)
        t0 = time.time()
        try:
            res = tests[name](ctx, args.tier)
        finally:
            if args.clean:
                ctx.clean()
        res["suite_wall"] = time.time() - t0
        res["scheme"] = args.scheme
        res["amr"] = args.amr
        res["refine_mass"] = args.refine_mass
        hc.save_json(os.path.join(args.results, f"{name}_{args.tier}{suffix}.json"), res)
        print(f"== {name}: {res['suite_wall']:.1f} s", flush=True)
    if ctx.fails:
        print("FAILED:\n  " + "\n  ".join(ctx.fails))
        return 1 if args.check else 0
    print("PASS" if args.check else "done")
    return 0


if __name__ == "__main__":
    sys.exit(main())
