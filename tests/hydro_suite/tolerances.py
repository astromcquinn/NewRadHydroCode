"""Pass/fail tolerances of the hydro accuracy suite.

Set from the measured accuracy (docs/hydro_tests.md, results/*.json) plus a
margin (typically 25-50 % on errors, or a fixed band around a measured
value), so that a change of the solver's accuracy is caught while round-off
level differences (thread count, compiler) are not. Each bound quotes the
value it was set from. Known deficiencies of the scheme are checked as
regression guards (the measured value must not get worse), not as
physics passes.
"""
import math


def _get(d, *keys):
    for k in keys:
        if d is None:
            return None
        d = d.get(k)
    return d


def check_tubes(ctx, out, tier):
    cs = out["cases"]
    # L1 bounds (rho, u, p) per (kind, N): measured (x/y/z/xy/xyz max) * ~1.3
    L1 = {("sod", 32): (0.052, 0.059, 0.050),    # measured <= 0.040, 0.045, 0.039
          ("sod", 64): (0.036, 0.039, 0.031),    # measured <= 0.028, 0.030, 0.024
          ("sod", 128): (0.019, 0.020, 0.016),   # measured <= 0.015, 0.016, 0.012
          ("strong", 32): (0.30, 0.050, 0.047),  # measured 0.232, 0.037, 0.036
          ("strong", 64): (0.23, 0.031, 0.028),  # measured 0.177, 0.023, 0.022
          ("strong", 128): (0.15, 0.016, 0.016),  # measured 0.116, 0.012, 0.012
          ("123", 32): (0.092, 0.22, 0.33),      # measured 0.071, 0.169, 0.252
          ("123", 64): (0.060, 0.10, 0.25),      # measured 0.046, 0.077, 0.196
          ("123", 128): (0.039, 0.061, 0.20)}    # measured 0.030, 0.047, 0.155
    for name, r in cs.items():
        key = (r["kind"], r["N"])
        if key in L1:
            for q, lim in zip(("rho", "u", "p"), L1[key]):
                ctx.check(r["l1"][q] <= lim, f"{name}: L1 {q} {r['l1'][q]:.4f} > {lim}")
        # transverse uniformity (exact symmetry of the problem): round-off only
        ctx.check(r["spread"] < 2e-5, f"{name}: transverse spread {r['spread']:.2e}")
        # mass: float32 round-off only
        ctx.check(abs(r["dmass"]) < 1e-6, f"{name}: mass drift {r['dmass']:.2e}")
        # positivity (123 and strong shock)
        ctx.check(r["min_rho"] > 0 and r["min_p"] > 0, f"{name}: non-positive rho or p")
        # lab-frame energy (H9): measured |dE| <= 3.2e-3 (123), 1.8e-3 (strong), 3.4e-4 (Sod)
        lim = {"sod": 6e-4, "strong": 4e-3, "123": 1.5e-2}[r["kind"]]
        ctx.check(abs(r["denergy"]) < lim, f"{name}: energy drift {r['denergy']:.2e} > {lim}")
        if r["kind"] == "sod" and r["N"] >= 64:
            # measured |err| <= 0.45 (shock), 0.45 (contact) cells
            for w in ("shock", "contact"):
                e = r["positions"][w]["err_cells"]
                ctx.check(math.isfinite(e) and abs(e) < 0.8,
                          f"{name}: {w} position error {e:+.2f} cells")
            # 10-90 % contact width measured 5.0-5.8 cells
            cw = r.get("contact_width")
            if cw is not None and math.isfinite(cw):
                ctx.check(cw < 7.0, f"{name}: contact width {cw:.2f} cells")
    for key, o in out["orders"].items():
        if key.startswith("sod"):
            # measured 0.62-0.85 (32-64-128), 0.55-0.75 (32-64)
            ctx.check(o["rho"] > 0.45, f"{key}: density L1 order {o['rho']:.2f} < 0.45")


def check_galilean(ctx, out, tier):
    for name, r in out["cases"].items():
        d = r["diff_vs_ref"]
        # measured max |drho|/rho_L <= LIM (see docs)
        lim = 0.06 if out["N"] <= 32 else 0.06
        ctx.check(d["rho_l1"] < 0.012, f"{name}: L1 difference to the reference "
                  f"{d['rho_l1']:.3e}")
        ctx.check(d["rho"] < lim, f"{name}: max difference to the reference {d['rho']:.3e}")
    ctx.check("H1" in out["h1_zero_bulk"],
              f"H1: zero bulk velocity did not raise the expected error: {out['h1_zero_bulk']}")


def check_sedov(ctx, out, tier):
    for name, r in out["cases"].items():
        N = r["N"]
        o = r["outputs"][-1]
        # shock radius: measured +10.7 % (N=32), +X % (64), +Y % (128)
        lim = {32: 0.16, 64: 0.10, 128: 0.07}[N]
        ctx.check(abs(o["R_err_rel"]) < lim, f"{name}: R error {o['R_err_rel']:+.3f}")
        lim = {32: 0.03, 64: 0.02, 128: 0.015}[N]
        ctx.check(o["aniso_rms"] < lim, f"{name}: anisotropy rms {o['aniso_rms']:.3f}")
        ctx.check(abs(o["dmass"]) < 1e-6, f"{name}: mass drift {o['dmass']:.2e}")
        ctx.check(abs(o["denergy"]) < 0.02, f"{name}: energy drift {o['denergy']:.2e}")
        ctx.check(abs(r["R_slope"] - 0.4) < 0.08, f"{name}: R ~ t^{r['R_slope']:.3f}")


def check_sound(ctx, out, tier):
    cs = out["cases"]
    for name, r in cs.items():
        ctx.check(abs(r["denergy"]) < 1e-6, f"{name}: energy drift {r['denergy']:.2e}")
    lim = {16: (0.14, 6e-3, 0.19), 32: (0.038, 1.6e-3, 0.031), 64: (0.012, 6e-4, 6e-3),
           128: (4e-3, 2e-4, 2e-3)}
    for N, (l1, sp, dmp) in lim.items():
        r = cs.get(f"sound_x_N{N}")
        if not r:
            continue
        ctx.check(r["l1_end"] < l1, f"sound N={N}: L1/A {r['l1_end']:.3e} > {l1}")
        ctx.check(abs(r["speed_err"]) < sp, f"sound N={N}: speed error {r['speed_err']:.2e}")
        ctx.check(r["damping_per_period"] < dmp,
                  f"sound N={N}: damping {r['damping_per_period']:.2e} > {dmp}")
    ctx.check(out["orders_fine"]["l1"] > 1.5, f"sound: L1 order {out['orders_fine']['l1']:.2f}")


def check_pancake(ctx, out, tier):
    for name, r in out["cases"].items():
        for o in r["outputs"]:
            zd = o.get("zeldovich")
            if zd and o["z"] >= 10 - 1e-6:
                # density/velocity follow Zel'dovich: measured L1 rho 2.6e-2,
                # v 7e-3 at N=32, z=10 (smaller at higher N)
                ctx.check(zd["l1_rho"] < 0.035, f"{name} z={o['z']}: Zel'dovich L1 rho "
                          f"{zd['l1_rho']:.3e}")
                ctx.check(zd["l1_v"] < 0.01, f"{name} z={o['z']}: Zel'dovich L1 v "
                          f"{zd['l1_v']:.3e}")
            ctx.check(o["spread"] < 1e-3, f"{name} z={o['z']}: not plane-parallel "
                      f"({o['spread']:.2e})")


def check_cosmobox(ctx, out, tier):
    for name, r in out["cases"].items():
        for o in r["outputs"]:
            if "pk_gas_over_dm_lowk" in o:
                g = o["pk_gas_over_dm_lowk"]
                ctx.check(0.8 < g[0] < 1.15, f"{name} z={o['z']}: P_gas/P_dm(k_min) {g[0]:.3f}")
                lin = o["pk_dm_over_linear_lowk"]
                if o["z"] >= 3 - 1e-6:
                    ctx.check(0.8 < lin[0] < 1.3,
                              f"{name} z={o['z']}: P_dm/linear(k_min) {lin[0]:.3f}")
            ctx.check(abs(o["dmass"]) < 1e-5, f"{name}: mass drift {o['dmass']:.2e}")


def check_kh(ctx, out, tier):
    for name, r in out["cases"].items():
        if r["N"] >= 64:
            rel = r["rate"] / r["rate_linear"] - 1
            ctx.check(math.isfinite(rel) and abs(rel) < 0.3,
                      f"{name}: growth rate {r['rate']:.4f} vs {r['rate_linear']:.4f}")
