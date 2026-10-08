#!/usr/bin/env python3
"""Self-test for tools/compare/compare.py.

Fabricates small (N=8, Nfreq=5) output directories in the validation-spec
format, injects known differences, and checks that compare.py detects them,
locates the worst cell correctly (1-based Fortran indices) and returns the
right exit codes.

Run:  python3 tools/compare/test_compare.py      (or: python3 -m pytest tools/compare)
"""

import io
import json
import os
import shutil
import sys
import tempfile
import unittest
from contextlib import redirect_stdout, redirect_stderr

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import compare  # noqa: E402

N = 8
NFREQ = 5
NSTEPS = 3


def fields_for(step, rng):
    """Return dict name -> array (ncomp, N, N, N) float32 with plausible values."""
    shp = (N, N, N)
    f = {}
    x = np.linspace(0, 1, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    r = np.sqrt((X - .5) ** 2 + (Y - .5) ** 2 + (Z - .5) ** 2)
    ion = np.clip(1.0 - np.exp(-(0.2 * step) / (r + 0.05)), 0, 1)
    f["fHI"] = np.maximum(1.0 - ion, 1e-6)[None]
    f["fHeI"] = np.maximum(1.0 - ion, 1e-6)[None]
    f["fHeII"] = (0.9 * ion)[None]
    f["temp"] = (100.0 + 2e4 * ion)[None]
    f["rhog"] = np.exp(rng.normal(0, 0.5, shp))[None]
    f["ngam"] = rng.random((NFREQ,) + shp) * 1e-3
    f["nu"] = np.broadcast_to(np.array([15.71, 21.06, 29.64, 42.88, 74.12])[:, None, None, None],
                              (NFREQ,) + shp).copy()
    f["fgam"] = rng.random((NFREQ,) + shp)
    f["tau"] = rng.random((NFREQ,) + shp) * 10
    f["gamma"] = rng.random((3,) + shp) * 1e-12
    f["lambda"] = rng.random((3,) + shp) * 1e-23
    f["dx"] = rng.normal(0, 1, (3,) + shp)
    f["u_gas"] = rng.random((5,) + shp)
    f["w_gas"] = rng.random((4,) + shp)
    return {k: np.asarray(v, dtype=np.float32) for k, v in f.items()}


def stats_for(step):
    return {
        "istep": step, "redshift": 10.0 - 0.01 * step, "scalefactor": 1 / (11.0 - 0.01 * step),
        "nstep_rt": 12, "nray_active": 1000 * step, "ncst1": 900 * step,
        "nmrg1": 50 * step, "nspl1": 40 * step, "ndel1": 3 * step,
        "ngam_rays_total": 0.123456789012 * step,
    }


def write_run(root, mutate=None, mutate_stats=None, ext=""):
    rng = np.random.default_rng(1234)
    for s in range(1, NSTEPS + 1):
        d = os.path.join(root, f"step_{s:04d}")
        os.makedirs(d)
        f = fields_for(s, rng)
        st = stats_for(s)
        if mutate:
            mutate(s, f)
        if mutate_stats:
            mutate_stats(s, st)
        for k, a in f.items():
            # a is (c, i, j, k); write Fortran order (c fastest)
            a.astype("<f4").ravel(order="F").tofile(os.path.join(d, k + ext))
        with open(os.path.join(d, "stats.txt"), "w") as fh:
            for k, v in st.items():
                fh.write(f"{k} = {v!r}\n" if isinstance(v, float) else f"{k} = {v}\n")


def run_compare(*args):
    out, err = io.StringIO(), io.StringIO()
    with redirect_stdout(out), redirect_stderr(err):
        rc = compare.main(list(args))
    return rc, out.getvalue() + err.getvalue()


class CompareTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="cmp_test_")
        self.ref = os.path.join(self.tmp, "ref")
        write_run(self.ref)

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def new(self, **kw):
        d = os.path.join(self.tmp, "new")
        write_run(d, **kw)
        return d

    def report(self, new, *extra):
        js = os.path.join(self.tmp, "rep.json")
        rc, out = run_compare(self.ref, new, "--json", js, *extra)
        with open(js) as f:
            return rc, out, json.load(f)

    # ------------------------------------------------------------------
    def test_identical_passes(self):
        rc, out, rep = self.report(self.new())
        self.assertEqual(rc, 0, out)
        self.assertTrue(rep["pass"])
        self.assertEqual(rep["N"], N)
        self.assertEqual(rep["nfreq"], NFREQ)
        self.assertTrue(all(r["n_diff"] == 0 for r in rep["fields"]))
        self.assertIn("RESULT: PASS", out)

    def test_roundoff_within_tolerance(self):
        def m(s, f):
            f["temp"] = (f["temp"].astype(np.float64) * (1 + 3e-7)).astype(np.float32)
        rc, out, rep = self.report(self.new(mutate=m))
        self.assertEqual(rc, 0, out)
        t = [r for r in rep["fields"] if r["field"] == "temp"]
        self.assertTrue(any(r["n_diff"] > 0 for r in t))
        self.assertTrue(all(r["ok"] for r in t))

    def test_injected_scalar_difference_located(self):
        # 0-based (i,j,k) = (2,5,6) -> 1-based (3,6,7), at step 2 and later
        def m(s, f):
            if s >= 2:
                f["fHI"][0, 2, 5, 6] *= 1.01
        rc, out, rep = self.report(self.new(mutate=m))
        self.assertEqual(rc, 1, out)
        self.assertEqual(rep["first_divergence"], 2)
        r = [r for r in rep["fields"] if r["field"] == "fHI" and r["step"] == 2][0]
        self.assertFalse(r["ok"])
        self.assertEqual(r["n_diff"], 1)
        self.assertEqual(r["n_fail"], 1)
        w = r["worst"]
        self.assertEqual((w["comp"], w["i"], w["j"], w["k"]), (1, 3, 6, 7))
        self.assertAlmostEqual(r["max_rel"], 0.01 / 1.01, places=4)
        r1 = [r for r in rep["fields"] if r["field"] == "fHI" and r["step"] == 1][0]
        self.assertTrue(r1["ok"])
        self.assertIn("FIRST DIVERGENCE beyond tolerance at step 0002", out)

    def test_fortran_flat_index_convention(self):
        # Independently of numpy reshapes: poke the raw file at the Fortran
        # flat offset of (c,i,j,k) = (2,4,1,7) (1-based) in a (3,N,N,N) field.
        new = self.new()
        c, i, j, k = 2, 4, 1, 7
        flat = (c - 1) + 3 * ((i - 1) + N * ((j - 1) + N * (k - 1)))
        p = os.path.join(new, "step_0001", "lambda")
        raw = np.fromfile(p, dtype="<f4")
        raw[flat] = raw[flat] * 3 + 1e-20
        raw.tofile(p)
        rc, out, rep = self.report(new)
        self.assertEqual(rc, 1)
        r = [r for r in rep["fields"] if r["field"] == "lambda" and r["step"] == 1][0]
        w = r["worst"]
        self.assertEqual((w["comp"], w["i"], w["j"], w["k"]), (c, i, j, k))

    def test_injected_multicomponent_difference_located(self):
        # ngam component 4 (1-based) at 0-based (7,0,3) -> (8,1,4)
        def m(s, f):
            if s == 3:
                f["ngam"][3, 7, 0, 3] += 0.5
                f["ngam"][1, 1, 1, 1] *= 1 + 2e-5   # smaller, also failing
        rc, out, rep = self.report(self.new(mutate=m))
        self.assertEqual(rc, 1)
        r = [r for r in rep["fields"] if r["field"] == "ngam" and r["step"] == 3][0]
        self.assertEqual(r["n_diff"], 2)
        self.assertEqual(r["n_fail"], 2)
        w = r["worst"]
        self.assertEqual((w["comp"], w["i"], w["j"], w["k"]), (4, 8, 1, 4))
        self.assertEqual(rep["first_divergence"], 3)

    def test_floor_suppresses_tiny_values(self):
        # a value ~1e-12 of the field max changing by 100% is below the floor
        def m(s, f):
            f["gamma"][0, 0, 0, 0] = 1e-25
        def mref(s, f):
            f["gamma"][0, 0, 0, 0] = 2e-25
        shutil.rmtree(self.ref)
        write_run(self.ref, mutate=mref)
        rc, out, rep = self.report(self.new(mutate=m))
        self.assertEqual(rc, 0, out)
        g = [r for r in rep["fields"] if r["field"] == "gamma"]
        self.assertTrue(all(r["n_diff"] == 1 and r["ok"] for r in g))

    def test_nan_detected(self):
        def m(s, f):
            if s == 1:
                f["tau"][2, 4, 4, 4] = np.nan
        rc, out, rep = self.report(self.new(mutate=m))
        self.assertEqual(rc, 1)
        r = [r for r in rep["fields"] if r["field"] == "tau" and r["step"] == 1][0]
        self.assertEqual(r["n_nan_mismatch"], 1)
        w = r["worst"]
        self.assertEqual((w["comp"], w["i"], w["j"], w["k"]), (3, 5, 5, 5))

    def test_integer_stat_exact(self):
        def ms(s, st):
            if s == 2:
                st["nmrg1"] += 1
        rc, out, rep = self.report(self.new(mutate_stats=ms))
        self.assertEqual(rc, 1)
        bad = [r for r in rep["stats"] if not r["ok"]]
        self.assertEqual([(b["step"], b["key"]) for b in bad], [(2, "nmrg1")])
        self.assertEqual(rep["first_divergence"], 2)
        # a relative tolerance on integer stats can be granted explicitly
        rc2, out2 = run_compare(self.ref, os.path.join(self.tmp, "new"),
                                "--stats-int-tol", "0.1")
        self.assertEqual(rc2, 0, out2)

    def test_float_stat_tolerance(self):
        def ms(s, st):
            st["ngam_rays_total"] *= (1 + 1e-9)
        rc, out, rep = self.report(self.new(mutate_stats=ms))
        self.assertEqual(rc, 0, out)
        def ms2(s, st):
            if s == 3:
                st["ngam_rays_total"] *= 1.001
        shutil.rmtree(os.path.join(self.tmp, "new"))
        rc, out, rep = self.report(self.new(mutate_stats=ms2))
        self.assertEqual(rc, 1)

    def test_field_rtol_override(self):
        def m(s, f):
            f["temp"][0, 1, 2, 3] *= 1.001
        new = self.new(mutate=m)
        rc, out = run_compare(self.ref, new)
        self.assertEqual(rc, 1)
        rc, out = run_compare(self.ref, new, "--field-rtol", "temp=1e-2")
        self.assertEqual(rc, 0, out)

    def test_histories(self):
        def m(s, f):
            f["fHI"][:] = f["fHI"] * 0.5
        rc, out, rep = self.report(self.new(mutate=m))
        self.assertEqual(rc, 1)
        for s in map(str, range(1, NSTEPS + 1)):
            hr, hn = rep["history_ref"][s], rep["history_new"][s]
            # xHII = 1 - fHI ; new fHI halved
            self.assertGreater(hn["xHII_V"], hr["xHII_V"])
            self.assertIn("xHII_M", hr)
            self.assertIn("T_V", hr)
            self.assertIn("ngam_rays_total", hr)

    def test_missing_step_and_field(self):
        new = self.new()
        shutil.rmtree(os.path.join(new, "step_0003"))
        rc, out = run_compare(self.ref, new)
        self.assertEqual(rc, 2, out)
        self.assertIn("steps only in REF: [3]", out)
        new2 = os.path.join(self.tmp, "new2")
        write_run(new2)
        os.remove(os.path.join(new2, "step_0001", "lambda"))
        rc, out = run_compare(self.ref, new2)
        self.assertEqual(rc, 2, out)
        self.assertIn("MISSING", out)

    def test_extension_and_explicit_N(self):
        new = os.path.join(self.tmp, "newbin")
        write_run(new, ext=".bin")
        rc, out = run_compare(self.ref, new, "--N", str(N), "--nfreq", str(NFREQ))
        self.assertEqual(rc, 0, out)

    def test_plots(self):
        try:
            import matplotlib  # noqa: F401
        except ImportError:
            self.skipTest("matplotlib not installed")
        def m(s, f):
            f["temp"][0, 4, 4, 4] *= 2
        pd = os.path.join(self.tmp, "plots")
        rc, out = run_compare(self.ref, self.new(mutate=m), "--plots", pd, "--plot-steps", "all")
        self.assertEqual(rc, 1)
        for s in range(1, NSTEPS + 1):
            self.assertTrue(os.path.isfile(os.path.join(pd, f"slice_step_{s:04d}.png")))
        self.assertTrue(os.path.isfile(os.path.join(pd, "history.png")))


NPART = 5 ** 3


def write_pm_run(root, mutate=None):
    """A PM-format run (pm_driver / harness_pm): no scalar RT fields, N is
    inferred from u_gas; xv_dm(6,Npart) particles."""
    rng = np.random.default_rng(99)
    for s in range(0, 3):
        d = os.path.join(root, f"step_{s:04d}")
        os.makedirs(d)
        f = {"u_gas": rng.random((5, N, N, N)), "v_gas": rng.normal(0, 1, (3, N, N, N)),
             "phi": rng.normal(0, 1, (1, N, N, N)), "xv_dm": rng.random((6, NPART)) * N}
        f = {k: np.asarray(v, dtype=np.float32) for k, v in f.items()}
        if mutate:
            mutate(s, f)
        for k, a in f.items():
            a.astype("<f4").ravel(order="F").tofile(os.path.join(d, k + ".bin"))
        with open(os.path.join(d, "stats.txt"), "w") as fh:
            fh.write(f"istep = {s}\nnstep_rt = 7\namax_dm = {0.1 * s!r}\n")


class ComparePMTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="cmp_pm_test_")
        self.ref = os.path.join(self.tmp, "ref")
        write_pm_run(self.ref)

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def report(self, mutate=None):
        new = os.path.join(self.tmp, "new")
        write_pm_run(new, mutate)
        js = os.path.join(self.tmp, "rep.json")
        rc, out = run_compare(self.ref, new, "--rtol", "0", "--atol-frac", "0", "--json", js)
        with open(js) as f:
            return rc, out, json.load(f)

    def test_identical_strict(self):
        rc, out, rep = self.report()
        self.assertEqual(rc, 0, out)
        self.assertEqual(rep["N"], N)
        names = {r["field"] for r in rep["fields"]}
        self.assertEqual(names, {"u_gas", "v_gas", "phi", "xv_dm"})

    def test_particle_difference_located(self):
        # one ulp in vy (comp 5) of particle ip = 17 (1-based) at step 1
        def m(s, f):
            if s == 1:
                f["xv_dm"][4, 16] = np.nextafter(f["xv_dm"][4, 16], np.float32(1e9))
        rc, out, rep = self.report(m)
        self.assertEqual(rc, 1, out)
        r = [r for r in rep["fields"] if r["field"] == "xv_dm" and r["step"] == 1][0]
        self.assertEqual(r["n_diff"], 1)
        self.assertEqual((r["worst"]["comp"], r["worst"]["ip"]), (5, 17))
        self.assertEqual(rep["first_divergence"], 1)
        self.assertIn("ip=17", out)

    def test_v_gas_difference_located(self):
        def m(s, f):
            if s == 2:
                f["v_gas"][2, 1, 2, 3] += 1.0
        rc, out, rep = self.report(m)
        self.assertEqual(rc, 1)
        r = [r for r in rep["fields"] if r["field"] == "v_gas" and r["step"] == 2][0]
        w = r["worst"]
        self.assertEqual((w["comp"], w["i"], w["j"], w["k"]), (3, 2, 3, 4))


def write_hydro_run(root, dmom=0.0, dfield=False):
    """Hydro-only layout: u_gas, w_gas, v_gas and multi-valued stats."""
    rng = np.random.default_rng(99)
    for s in range(0, NSTEPS + 1):
        d = os.path.join(root, f"step_{s:04d}")
        os.makedirs(d)
        for name, nc in (("u_gas", 5), ("w_gas", 4), ("v_gas", 3)):
            a = rng.random((nc, N, N, N)).astype(np.float32)
            if dfield and name == "v_gas" and s == NSTEPS:
                a[2, 1, 2, 3] = np.nextafter(a[2, 1, 2, 3], np.float32(2))
            a.astype("<f4").ravel(order="F").tofile(os.path.join(d, name + ".bin"))
        with open(os.path.join(d, "stats.txt"), "w") as fh:
            fh.write(f"istep = {s}\nmass =  3.0720000950098038E+003\n")
            fh.write(f"momentum =  1.5E+002 {-2.0 + dmom * s!r} 3.25D0\n")


class HydroLayoutTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="cmp_hydro_")
        self.ref = os.path.join(self.tmp, "ref")
        write_hydro_run(self.ref)

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def test_hydro_identical(self):
        new = os.path.join(self.tmp, "new")
        write_hydro_run(new)
        rc, out = run_compare(self.ref, new, "--rtol", "0", "--atol-frac", "0")
        self.assertEqual(rc, 0, out)
        self.assertIn(f"N = {N}", out)
        self.assertIn("mass", out)

    def test_hydro_multivalued_stat_and_ulp(self):
        new = os.path.join(self.tmp, "new")
        write_hydro_run(new, dmom=1e-12)
        rc, out = run_compare(self.ref, new, "--rtol", "0", "--atol-frac", "0", "-q")
        self.assertEqual(rc, 1, out)
        self.assertIn("stat momentum", out)
        new2 = os.path.join(self.tmp, "new2")
        write_hydro_run(new2, dfield=True)
        rc, out = run_compare(self.ref, new2, "--rtol", "0", "--atol-frac", "0", "-q")
        self.assertEqual(rc, 1, out)
        self.assertIn("v_gas", out)
        self.assertIn("(c=3,2,3,4)", out)


if __name__ == "__main__":
    unittest.main(verbosity=2)
