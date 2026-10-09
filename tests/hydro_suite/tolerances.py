"""Pass/fail tolerances of the hydro accuracy suite.

Set from the measured accuracy (docs/hydro_tests.md, results/*.json) plus a
margin (typically 25-50 % on errors, or a fixed band around a measured
value), so that a change of the solver's accuracy is caught while round-off
level differences (thread count, compiler) are not. Each bound quotes the
value it was set from. Known deficiencies of the scheme are checked as
regression guards (the measured value must not get worse), not as
physics passes.

Two schemes (ctx.scheme): "legacy", the original moving-frame TVD scheme,
and "unsplit", the unsplit MUSCL-Hancock + HLLC scheme with dual energy
(docs/hydro_unsplit.md), each with its own limits from its own measured
values (results/*_full.json and results/*_full_unsplit.json).
"""
import math


def _get(d, *keys):
    for k in keys:
        if d is None:
            return None
        d = d.get(k)
    return d


def check_tubes(ctx, out, tier):
    if ctx.scheme == "unsplit":
        return check_tubes_unsplit(ctx, out, tier)
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
        # (Sod, 123: <= 2e-7 measured). The strong shock breaks the planar
        # symmetry: 5.3e-5 (N=64), 3.1e-2 (N=128) -- regression guard only
        lim = 0.05 if r["kind"] == "strong" else 2e-5
        ctx.check(r["spread"] < lim, f"{name}: transverse spread {r['spread']:.2e}")
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
    if ctx.scheme == "unsplit":
        return check_galilean_unsplit(ctx, out, tier)
    # the profile is advected through the grid by vb*t cells (1/64 .. 1 box);
    # advection diffusion makes it depend on vb (measured at N = 64: L1
    # difference <= 0.022, max <= 0.047 for a full box)
    for name, r in out["cases"].items():
        d = r["diff_vs_ref"]
        ctx.check(d["rho_l1"] < 0.03, f"{name}: L1 difference to the reference "
                  f"{d['rho_l1']:.3e}")
        ctx.check(d["rho"] < 0.065, f"{name}: max difference to the reference {d['rho']:.3e}")
    # the transverse bulk velocity (0.1 vs 2 cells/time): measured <= 8.8e-4
    names = list(out["cases"])
    for a, b in zip(names[::2], names[1::2]):
        da = out["cases"][a]["profile"]["rho"]
        db = out["cases"][b]["profile"]["rho"]
        diff = max(abs(x - y) for x, y in zip(da, db))
        ctx.check(diff < 2e-3, f"{a} vs {b}: transverse bulk changes rho by {diff:.2e}")
    ctx.check("H1" in out["h1_zero_bulk"],
              f"H1: zero bulk velocity did not raise the expected error: {out['h1_zero_bulk']}")


def check_sedov(ctx, out, tier):
    if ctx.scheme == "unsplit":
        return check_sedov_unsplit(ctx, out, tier)
    for name, r in out["cases"].items():
        N = r["N"]
        o = r["outputs"][-1]
        # shock radius: measured +10.7 % (N=32), +4.1 % (64), +0.9 % (128)
        lim = {32: 0.15, 64: 0.065, 128: 0.03}[N]
        ctx.check(abs(o["R_err_rel"]) < lim, f"{name}: R error {o['R_err_rel']:+.3f}")
        # rms of R over 400 directions: measured 1.85 %, 0.98 %, 0.50 %
        lim = {32: 0.025, 64: 0.0135, 128: 0.007}[N]
        ctx.check(o["aniso_rms"] < lim, f"{name}: anisotropy rms {o['aniso_rms']:.3f}")
        ctx.check(abs(o["dmass"]) < 1e-6, f"{name}: mass drift {o['dmass']:.2e}")
        # lab-frame energy is lost (H9; measured -0.86 %, -2.4 %, -4.9 % at
        # N = 32, 64, 128): regression guard
        lim = {32: 0.013, 64: 0.035, 128: 0.065}[N]
        ctx.check(abs(o["denergy"]) < lim, f"{name}: energy drift {o['denergy']:.2e}")
        # L1 (rho/rho2): measured 0.087, 0.051, 0.020
        lim = {32: 0.115, 64: 0.068, 128: 0.028}[N]
        ctx.check(o["l1"]["rho"] < lim, f"{name}: L1 rho {o['l1']['rho']:.3f}")
        # R ~ t^alpha over t = 1.25..10: measured 0.341, 0.364, 0.381
        ctx.check(abs(r["R_slope"] - 0.4) < 0.075, f"{name}: R ~ t^{r['R_slope']:.3f}")


def check_sound(ctx, out, tier):
    if ctx.scheme == "unsplit":
        return check_sound_unsplit(ctx, out, tier)
    cs = out["cases"]
    for name, r in cs.items():
        ctx.check(abs(r["denergy"]) < 1e-6, f"{name}: energy drift {r['denergy']:.2e}")
    # (L1/A, |phase speed error|/c, amplitude loss per period); measured
    # (0.108, 5.6e-3, 0.146), (0.031, 1.7e-3, 0.024), (7.4e-3, 2.4e-4, 1.8e-3),
    # (3.2e-3, 7.9e-5, 2.9e-4) for N = 16, 32, 64, 128 (A = 1e-3)
    lim = {16: (0.14, 7.5e-3, 0.19), 32: (0.04, 2.3e-3, 0.031), 64: (0.0098, 3.5e-4, 2.4e-3),
           128: (4.2e-3, 1.2e-4, 4e-4)}
    for N, (l1, sp, dmp) in lim.items():
        r = cs.get(f"sound_x_N{N}")
        if not r:
            continue
        ctx.check(r["l1_end"] < l1, f"sound N={N}: L1/A {r['l1_end']:.3e} > {l1}")
        ctx.check(abs(r["speed_err"]) < sp, f"sound N={N}: speed error {r['speed_err']:.2e}")
        ctx.check(r["damping_per_period"] < dmp,
                  f"sound N={N}: damping {r['damping_per_period']:.2e} > {dmp}")
    # error of the k = 1 Fourier mode: measured order 2.5 (16..128), 2.2 (16..64)
    ctx.check(out["orders"]["k1"] > 1.8, f"sound: k=1 error order {out['orders']['k1']:.2f}")


def check_pancake(ctx, out, tier):
    if ctx.scheme == "unsplit":
        return check_pancake_unsplit(ctx, out, tier)
    for name, r in out["cases"].items():
        for o in r["outputs"]:
            zd = o.get("zeldovich")
            if zd and o["z"] >= 10 - 1e-6:
                # density/velocity follow Zel'dovich: measured L1 rho 2.6e-2,
                # v 7e-3 at N=32, z=10 (smaller at higher N)
                ctx.check(zd["l1_rho"] < 0.035, f"{name} z={o['z']}: Zel'dovich L1 rho "
                          f"{zd['l1_rho']:.3e}")
                ctx.check(zd["l1_v"] < 0.013, f"{name} z={o['z']}: Zel'dovich L1 v "
                          f"{zd['l1_v']:.3e}")
            # accretion shock vs the 1-D Lagrangian reference (z <= 1, N >= 64):
            # measured |error| <= 0.93 cells (N = 128, z = 0.5)
            if r["N"] >= 64 and o["z"] <= 1 + 1e-6 and "shock_err_cells" in o:
                e = o["shock_err_cells"]
                ctx.check(math.isfinite(e) and abs(e) < 1.3,
                          f"{name} z={o['z']}: shock position error {e:+.2f} cells")
            # transverse symmetry breaking: measured <= 2.2e-2 of rho_max
            ctx.check(o["spread"] < 0.04, f"{name} z={o['z']}: not plane-parallel "
                      f"({o['spread']:.2e})")


def check_cosmobox(ctx, out, tier):
    if ctx.scheme == "unsplit":
        return check_cosmobox_unsplit(ctx, out, tier)
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
    if ctx.scheme == "unsplit":
        return check_kh_unsplit(ctx, out, tier)
    for name, r in out["cases"].items():
        amps = r["amp"]
        ctx.check(all(math.isfinite(a) for a in amps) and max(amps) < 0.5,
                  f"{name}: KH mode amplitude not finite/bounded")
        if r["N"] >= 128:
            # layer resolved by 4 cells: grows at ~0.6 of the incompressible
            # linear rate (Mach 0.5 shear and numerical viscosity, see docs)
            rel = r["rate"] / r["rate_linear"]
            ctx.check(math.isfinite(rel) and 0.4 < rel < 1.1,
                      f"{name}: growth rate {r['rate']:.4f} vs {r['rate_linear']:.4f}")


# ======================================================================
# The unsplit scheme (docs/hydro_unsplit.md). PROVISIONAL limits until the
# full tier has been measured; every value is replaced by measured * margin.
# ======================================================================
def check_tubes_unsplit(ctx, out, tier):
    cs = out["cases"]
    L1 = {("sod", 32): (0.040, 0.045, 0.039),
          ("sod", 64): (0.028, 0.030, 0.024),
          ("sod", 128): (0.015, 0.016, 0.012),
          ("strong", 32): (0.30, 0.050, 0.047),
          ("strong", 64): (0.23, 0.031, 0.028),
          ("strong", 128): (0.15, 0.016, 0.016),
          ("123", 32): (0.092, 0.22, 0.33),
          ("123", 64): (0.060, 0.10, 0.25),
          ("123", 128): (0.039, 0.061, 0.20)}
    for name, r in cs.items():
        key = (r["kind"], r["N"])
        if key in L1:
            for q, lim in zip(("rho", "u", "p"), L1[key]):
                ctx.check(r["l1"][q] <= lim, f"{name}: L1 {q} {r['l1'][q]:.4f} > {lim}")
        ctx.check(r["spread"] < 1e-5, f"{name}: transverse spread {r['spread']:.2e}")
        ctx.check(abs(r["dmass"]) < 1e-6, f"{name}: mass drift {r['dmass']:.2e}")
        ctx.check(r["min_rho"] > 0 and r["min_p"] > 0, f"{name}: non-positive rho or p")
        # conservative scheme: energy to round-off apart from the floors
        ctx.check(abs(r["denergy"]) < 1e-4, f"{name}: energy drift {r['denergy']:.2e}")
        if r["kind"] == "sod" and r["N"] >= 64:
            for w in ("shock", "contact"):
                e = r["positions"][w]["err_cells"]
                ctx.check(math.isfinite(e) and abs(e) < 0.8,
                          f"{name}: {w} position error {e:+.2f} cells")
            cw = r.get("contact_width")
            if cw is not None and math.isfinite(cw):
                ctx.check(cw < 5.0, f"{name}: contact width {cw:.2f} cells")
    for key, o in out["orders"].items():
        if key.startswith("sod"):
            ctx.check(o["rho"] > 0.6, f"{key}: density L1 order {o['rho']:.2f} < 0.6")


def check_galilean_unsplit(ctx, out, tier):
    for name, r in out["cases"].items():
        d = r["diff_vs_ref"]
        ctx.check(d["rho_l1"] < 0.05, f"{name}: L1 difference to the reference "
                  f"{d['rho_l1']:.3e}")
        ctx.check(d["rho"] < 0.1, f"{name}: max difference to the reference {d['rho']:.3e}")
    names = list(out["cases"])
    for a, b in zip(names[::2], names[1::2]):
        da = out["cases"][a]["profile"]["rho"]
        db = out["cases"][b]["profile"]["rho"]
        diff = max(abs(x - y) for x, y in zip(da, db))
        # the CFL substep count depends on the transverse speed (U3):
        # measured 6.3e-3 at N = 32
        ctx.check(diff < 1e-2, f"{a} vs {b}: transverse bulk changes rho by {diff:.2e}")
    # a lab-frame scheme has no H1: the zero-bulk tube runs
    ctx.check(out["h1_zero_bulk"].startswith("ran"),
              f"zero bulk velocity did not run: {out['h1_zero_bulk']}")
    z = out.get("zero_bulk")
    if z:
        ctx.check(z["diff_vs_ref"]["rho_l1"] < 0.05,
                  f"zero bulk: L1 difference to the reference {z['diff_vs_ref']['rho_l1']:.3e}")


def check_sedov_unsplit(ctx, out, tier):
    for name, r in out["cases"].items():
        N = r["N"]
        o = r["outputs"][-1]
        lim = {32: 0.15, 64: 0.065, 128: 0.03}[N]
        ctx.check(abs(o["R_err_rel"]) < lim, f"{name}: R error {o['R_err_rel']:+.3f}")
        lim = {32: 0.04, 64: 0.0135, 128: 0.007}[N]
        ctx.check(o["aniso_rms"] < lim, f"{name}: anisotropy rms {o['aniso_rms']:.3f}")
        ctx.check(abs(o["dmass"]) < 1e-6, f"{name}: mass drift {o['dmass']:.2e}")
        ctx.check(abs(o["denergy"]) < 1e-3, f"{name}: energy drift {o['denergy']:.2e}")
        lim = {32: 0.115, 64: 0.068, 128: 0.028}[N]
        ctx.check(o["l1"]["rho"] < lim, f"{name}: L1 rho {o['l1']['rho']:.3f}")
        ctx.check(abs(r["R_slope"] - 0.4) < 0.075, f"{name}: R ~ t^{r['R_slope']:.3f}")


def check_sound_unsplit(ctx, out, tier):
    cs = out["cases"]
    for name, r in cs.items():
        ctx.check(abs(r["denergy"]) < 1e-6, f"{name}: energy drift {r['denergy']:.2e}")
    lim = {16: (0.14, 7.5e-3, 0.19), 32: (0.04, 2.3e-3, 0.031), 64: (0.0098, 3.5e-4, 2.4e-3),
           128: (4.2e-3, 1.2e-4, 4e-4)}
    for N, (l1, sp, dmp) in lim.items():
        r = cs.get(f"sound_x_N{N}")
        if not r:
            continue
        ctx.check(r["l1_end"] < l1, f"sound N={N}: L1/A {r['l1_end']:.3e} > {l1}")
        ctx.check(abs(r["speed_err"]) < sp, f"sound N={N}: speed error {r['speed_err']:.2e}")
        ctx.check(r["damping_per_period"] < dmp,
                  f"sound N={N}: damping {r['damping_per_period']:.2e} > {dmp}")
    ctx.check(out["orders"]["k1"] > 1.8, f"sound: k=1 error order {out['orders']['k1']:.2f}")


def check_pancake_unsplit(ctx, out, tier):
    for name, r in out["cases"].items():
        for o in r["outputs"]:
            zd = o.get("zeldovich")
            if zd and o["z"] >= 10 - 1e-6:
                ctx.check(zd["l1_rho"] < 0.035, f"{name} z={o['z']}: Zel'dovich L1 rho "
                          f"{zd['l1_rho']:.3e}")
                ctx.check(zd["l1_v"] < 0.013, f"{name} z={o['z']}: Zel'dovich L1 v "
                          f"{zd['l1_v']:.3e}")
            if zd and r["N"] >= 64 and o["z"] >= 7 - 1e-6:
                # the point of the dual-energy formulation: the cold infall
                # stays adiabatic (median T/T_adiabatic within a factor of 2)
                ctx.check(0.5 < zd["med_T"] < 2.0,
                          f"{name} z={o['z']}: median T/T_adiabatic {zd['med_T']:.3g}")
            if r["N"] >= 64 and o["z"] <= 1 + 1e-6 and "shock_err_cells" in o:
                e = o["shock_err_cells"]
                ctx.check(math.isfinite(e) and abs(e) < 1.3,
                          f"{name} z={o['z']}: shock position error {e:+.2f} cells")
            ctx.check(o["spread"] < 0.04, f"{name} z={o['z']}: not plane-parallel "
                      f"({o['spread']:.2e})")


def check_cosmobox_unsplit(ctx, out, tier):
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
            if o["z"] >= 10 - 1e-6:
                # unshocked gas near its adiabat (or the 1 K floor), not 1e4 K
                ctx.check(o["T_med"] < 100.0,
                          f"{name} z={o['z']}: median T {o['T_med']:.3g} K (adiabatic "
                          f"{o['T_adiabatic_mean']:.3g} K)")


def check_kh_unsplit(ctx, out, tier):
    for name, r in out["cases"].items():
        amps = r["amp"]
        ctx.check(all(math.isfinite(a) for a in amps) and max(amps) < 0.5,
                  f"{name}: KH mode amplitude not finite/bounded")
        if r["N"] >= 128:
            rel = r["rate"] / r["rate_linear"]
            ctx.check(math.isfinite(rel) and 0.4 < rel < 1.1,
                      f"{name}: growth rate {r['rate']:.4f} vs {r['rate_linear']:.4f}")
