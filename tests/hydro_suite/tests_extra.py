"""Sedov, sound, pancake, cosmological box and Kelvin-Helmholtz sections of
the hydro suite runner (run_suite.py)."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hs_common as hc  # noqa: E402
import generate as gen  # noqa: E402
import analyze as an  # noqa: E402
import cosmo_gen as cg  # noqa: E402
import tolerances as tol  # noqa: E402


def strip_profiles(res, keep):
    """drop the bulky profile arrays except for the cases in `keep`"""
    return res


# ----------------------------------------------------------------------
def test_sedov(ctx, tier):
    Ns = (32, 64) if tier == "quick" else (32, 64, 128)
    cases, metas = [], {}
    for N in Ns:
        name = f"sedov_N{N}"
        pdir, odir, meta = ctx.case(name, gen.gen_sedov, N=N)
        cases.append((name, pdir, odir))
        metas[name] = meta
    walls = ctx.run_many(cases)
    res = {}
    for name, pdir, odir in cases:
        r = an.analyze_sedov(odir, metas[name])
        r["wall"] = walls[name]
        res[name] = r
        for o in r["outputs"]:
            print(f"  {name:12s} t={o['t']:6.3f} steps={o['nsteps']:4d} R={o['R_mean_dirs']:7.3f} "
                  f"(exact {o['R_exact']:7.3f}, {100 * o['R_err_rel']:+5.1f}%) aniso rms "
                  f"{100 * o['aniso_rms']:.2f}% axes/diag {100 * o['axes_over_diag']:+.2f}% "
                  f"peak {o['peak_rho_over_rho2']:.2f} L1 rho {o['l1']['rho']:.3f} "
                  f"p {o['l1']['p']:.3f} v {o['l1']['v']:.3f} dE {o['denergy']:+.2e}")
        print(f"  {name}: R ~ t^{r.get('R_slope', float('nan')):.3f} (0.4), {r['wall']:.1f} s")
    orders = {}
    keys = [f"sedov_N{N}" for N in Ns]
    for q in ("rho", "p", "v"):
        orders["l1_" + q] = an.fit_order(Ns, [res[k]["outputs"][-1]["l1"][q] for k in keys])
    orders["R_err"] = an.fit_order(Ns, [abs(res[k]["outputs"][-1]["R_err_rel"]) for k in keys])
    orders["aniso"] = an.fit_order(Ns, [res[k]["outputs"][-1]["aniso_rms"] for k in keys])
    print("  orders:", {k: round(v, 2) for k, v in orders.items()})
    out = {"cases": res, "orders": orders}
    tol.check_sedov(ctx, out, tier)
    return out


# ----------------------------------------------------------------------
def test_sound(ctx, tier):
    specs = [("x", N) for N in ((16, 32, 64) if tier == "quick" else (16, 32, 64, 128))]
    specs += [("xy", N) for N in ((32,) if tier == "quick" else (32, 64))]
    cases, metas = [], {}
    for d, N in specs:
        name = f"sound_{d}_N{N}"
        pdir, odir, meta = ctx.case(name, gen.gen_sound, N=N, direction=d, nout=8)
        cases.append((name, pdir, odir))
        metas[name] = meta
    walls = ctx.run_many(cases)
    res = {}
    for name, pdir, odir in cases:
        r = an.analyze_sound(odir, metas[name])
        r["wall"] = walls[name]
        res[name] = r
        print(f"  {name:16s} steps {r['nsteps']:4d} phase speed err {r['speed_err']:+.2e} "
              f"damping/period {r['damping_per_period']:.2e} L1/A {r['l1_end']:.3e} "
              f"k=1 error/A {r['k1_err_end']:.3e} "
              f"{r['wall']:.1f}s")
    Ns = [N for d, N in specs if d == "x"]
    keys = [f"sound_x_N{N}" for N in Ns]
    orders = {"l1": an.fit_order(Ns, [res[k]["l1_end"] for k in keys]),
              "k1": an.fit_order(Ns, [res[k]["k1_err_end"] for k in keys]),
              "speed": an.fit_order(Ns, [abs(res[k]["speed_err"]) for k in keys]),
              "damping": an.fit_order(Ns, [res[k]["damping_per_period"] for k in keys])}
    # orders from the two finest grids only (asymptotic regime)
    orders_fine = {"l1": an.fit_order(Ns[-2:], [res[k]["l1_end"] for k in keys[-2:]]),
                   "k1": an.fit_order(Ns[-2:], [res[k]["k1_err_end"] for k in keys[-2:]]),
                   "speed": an.fit_order(Ns[-2:], [abs(res[k]["speed_err"]) for k in keys[-2:]]),
                   "damping": an.fit_order(Ns[-2:], [res[k]["damping_per_period"]
                                                     for k in keys[-2:]])}
    print("  orders (all):", {k: round(v, 2) for k, v in orders.items()},
          " (finest pair):", {k: round(v, 2) for k, v in orders_fine.items()})
    out = {"cases": res, "orders": orders, "orders_fine": orders_fine}
    tol.check_sound(ctx, out, tier)
    return out


# ----------------------------------------------------------------------
def pancake_reference(ctx, meta, ng):
    path = os.path.join(ctx.rundir, f"pancake_ref1d_ng{ng}_T{meta['T_init']:g}.npz")
    if not os.path.exists(path):
        import pancake_ref1d
        outs = pancake_ref1d.run(meta, ng=ng, ndm=ng)
        flat = {f"{i}_{k}": v for i, o in enumerate(outs) for k, v in o.items()}
        flat["n"] = len(outs)
        np.savez_compressed(path, **flat)
    return np.load(path)


def test_pancake(ctx, tier):
    if tier == "quick":
        specs = [dict(N=32, T_init=1000.0, z_out=(10, 5.5, 3))]
        ng = 512
    else:
        specs = [dict(N=64, T_init=1000.0), dict(N=128, T_init=1000.0),
                 dict(N=64, T_init=1e4), dict(N=32, T_init=1000.0)]
        ng = 2048
    cases, metas = [], {}
    for s in specs:
        name = f"pancake_N{s['N']}_T{s['T_init']:g}"
        pdir, odir, meta = ctx.case(name, cg.gen_pancake, **s)
        cases.append((name, pdir, odir))
        metas[name] = meta
    walls = ctx.run_many(cases)
    res = {}
    refs = {}
    for name, pdir, odir in cases:
        meta = metas[name]
        key = meta["T_init"]
        if key not in refs:
            refs[key] = pancake_reference(ctx, meta, ng)
        r = an.analyze_pancake(odir, meta, refs[key])
        r["wall"] = walls[name]
        r["T_init"] = meta["T_init"]
        res[name] = r
        for o in r["outputs"]:
            zd = o.get("zeldovich")
            s = (f"  {name:20s} z={o['z']:5.2f} steps={o['nsteps']:4d} rho_max={o['rho_max']:7.2f} "
                 f"T_max={o['T_max']:.3g}")
            if zd:
                s += (f" | Zel'dovich L1 rho {zd['l1_rho']:.2e} v {zd['l1_v']:.2e} "
                      f"T {zd['l1_T']:.2e} (max {zd['max_T']:.2e})")
            else:
                s += (f" | shock {np.nanmean(o['shock_halfwidth']) * r['N']:.2f} cells "
                      f"(ref {np.nanmean(o['ref']['shock_halfwidth']) * r['N']:.2f})"
                      if "ref" in o else "")
                if "K_away_median_dev" in o:
                    s += f" K-1 away: median {o['K_away_median_dev']:.2e} max {o['K_away_max_dev']:.2e}"
            print(s)
        print(f"  {name}: {r['wall']:.1f} s")
    out = {"cases": res, "ref_ng": ng}
    tol.check_pancake(ctx, out, tier)
    return out


# ----------------------------------------------------------------------
def analyze_cosmobox(run, meta):
    import pmcosmo as pc
    N = meta["N"]
    fb = cg.pc.FB
    res = {"N": N, "outputs": []}
    G = pc.Growth(1.0)
    pk0 = None
    for d in hc.outputs(run):
        rho, vel, p, info = hc.read_state(d, N)
        z = info["redshift"]
        T = p / ((hc.GAMMA - 1) * rho) * info["temp1_unit"]
        delta = rho / fb
        o = {"z": z, "nsteps": int(info["istep"]), "wall": info["wall"],
             "T_mw": float(np.sum(rho * T) / np.sum(rho)), "T_vw": float(np.mean(T)),
             "T_med": float(np.median(T)), "delta_max": float(delta.max()),
             "frac_T_floor": float(np.mean(T < 1.05)),
             "frac_T_gt_1e5": float(np.sum(rho[T > 1e5]) / np.sum(rho))}
        # T-rho relation: median log T in bins of log delta
        ld = np.log10(delta).ravel()
        lt = np.log10(np.maximum(T, 1e-3)).ravel()
        bins = np.linspace(-1.5, 3, 19)
        idx = np.digitize(ld, bins) - 1
        med, p16, p84, cen = [], [], [], []
        for b in range(len(bins) - 1):
            s = lt[idx == b]
            if s.size >= 10:
                cen.append(0.5 * (bins[b] + bins[b + 1]))
                med.append(float(np.median(s)))
                p16.append(float(np.percentile(s, 16)))
                p84.append(float(np.percentile(s, 84)))
        o["trho"] = {"logd": cen, "logT_med": med, "logT_16": p16, "logT_84": p84}
        # adiabatic expectation for unshocked gas: T = T_i (a_i/a)^2 (delta)^(2/3)
        o["T_adiabatic_mean"] = meta["T_init"] * ((1 + z) / (1 + meta["z_init"])) ** 2
        if z < 0.01:
            H, xe, ye = np.histogram2d(ld, lt, bins=(np.linspace(-1.5, 3, 46),
                                                     np.linspace(-1, 8.5, 58)))
            o["trho_hist"] = {"H": H.astype(int).tolist(), "x": xe.tolist(), "y": ye.tolist()}
        pkf = os.path.join(d, "pk.txt")
        if os.path.exists(pkf):
            pk = np.loadtxt(pkf)
            o["pk"] = {"k": pk[:, 0].tolist(), "dm": pk[:, 3].tolist(), "dm_raw": pk[:, 1].tolist(),
                       "gas": pk[:, 5].tolist()}
            if pk0 is None:
                pk0 = (pk, info["scalefactor"])
            else:
                a0 = pk0[1]
                grow = (G(info["scalefactor"]) / G(a0)) ** 2
                o["pk_dm_over_linear_lowk"] = (pk[:3, 3] / (pk0[0][:3, 3] * grow)).tolist()
                o["pk_gas_over_dm_lowk"] = (pk[:3, 5] / pk[:3, 3]).tolist()
        lg = hc.read_log(run)
        o["dmass"] = float(info["mass"] / lg["mass"][0] - 1) if len(lg["mass"]) else 0.0
        res["outputs"].append(o)
    return res


def test_cosmobox(ctx, tier):
    if tier == "quick":
        specs = [dict(N=32, box=64.0, z_out=(10, 3))]
    else:
        specs = [dict(N=64, box=64.0), dict(N=64, box=64.0, dlna_max=0.02),
                 dict(N=32, box=64.0)]
    cases, metas = [], {}
    for s in specs:
        name = f"cosmobox_N{s['N']}" + (f"_dlna{s['dlna_max']:g}" if s.get("dlna_max") else "")
        pdir, odir, meta = ctx.case(name, cg.gen_cosmobox, **s)
        cases.append((name, pdir, odir))
        metas[name] = meta
    walls = ctx.run_many(cases)
    res = {}
    for name, pdir, odir in cases:
        r = analyze_cosmobox(odir, metas[name])
        r["wall"] = walls[name]
        r["dlna_max"] = metas[name]["dlna_max"]
        res[name] = r
        for o in r["outputs"]:
            print(f"  {name:22s} z={o['z']:6.2f} steps={o['nsteps']:4d} T_mw={o['T_mw']:.3g} K "
                  f"T_vw={o['T_vw']:.3g} delta_max={o['delta_max']:.3g} f(T>1e5)="
                  f"{o['frac_T_gt_1e5']:.3f} floor {o['frac_T_floor']:.3f} "
                  + (f"Pdm/lin {np.round(o['pk_dm_over_linear_lowk'], 2)} "
                     f"Pgas/Pdm {np.round(o['pk_gas_over_dm_lowk'], 2)}"
                     if "pk_gas_over_dm_lowk" in o else ""))
        print(f"  {name}: {r['wall']:.1f} s")
    out = {"cases": res}
    tol.check_cosmobox(ctx, out, tier)
    return out


# ----------------------------------------------------------------------
def analyze_kh(run, meta):
    N = meta["N"]
    k = meta["k"]
    c = np.arange(N) + 0.5
    ts, amps = [], []
    for d in hc.outputs(run):
        rho, vel, p, info = hc.read_state(d, N)
        vy = vel[1].mean(axis=2) - meta["bulk"][1]     # (x, y)
        # mode amplitude of vy(x) at the interfaces (|y - N/4| < N/8)
        sel = np.abs(c - N / 4) < N / 8
        F = np.abs((vy[:, sel] * np.exp(-1j * k * c)[:, None]).sum(axis=0)) * 2 / N
        ts.append(info["t"])
        amps.append(float(F.max()))
    ts, amps = np.array(ts), np.array(amps)
    # growth rate: fit log amplitude where it is between 3x initial and 0.1 dv
    sel = (ts > 0.25 * ts[-1]) & (amps < 0.1 * meta["dv"])
    rate = float(np.polyfit(ts[sel], np.log(amps[sel]), 1)[0]) if sel.sum() >= 3 else float("nan")
    return {"N": N, "t": ts.tolist(), "amp": amps.tolist(), "rate": rate}


def kh_linear_rate(k, dv, a, c0):
    """growth rate of the tanh shear layer vx = (dv/2) tanh(y/a), uniform
    density, from the incompressible Rayleigh equation solved numerically
    (shooting); compressibility at Mach dv/c0 = 0.5 lowers it slightly
    (the Mach number correction is not included)."""
    from numpy.linalg import eigvals
    # finite-difference eigenproblem of the Rayleigh equation
    # (U - c)(phi'' - k^2 phi) - U'' phi = 0 on y in [-L, L], phi(+-L) = 0
    L = 12 * a
    n = 600
    y = np.linspace(-L, L, n + 2)[1:-1]
    h = y[1] - y[0]
    U = 0.5 * dv * np.tanh(y / a)
    Upp = -dv / a ** 2 * np.tanh(y / a) / np.cosh(y / a) ** 2
    D2 = (np.diag(np.full(n - 1, 1.0), -1) - 2 * np.eye(n) + np.diag(np.full(n - 1, 1.0), 1)) / h ** 2
    Lap = D2 - k ** 2 * np.eye(n)
    A = np.diag(U) @ Lap - np.diag(Upp)
    cs = eigvals(np.linalg.solve(Lap, A))
    return float(k * np.max(cs.imag))


def test_kh(ctx, tier):
    Ns = (32,) if tier == "quick" else (32, 64, 128)
    cases, metas = [], {}
    for N in Ns:
        name = f"kh_N{N}"
        # same physical problem at every N: two wavelengths per box, layer
        # width N/32 cells (k a = 0.39), Mach 0.5 shear, run for ~6 e-folds
        # of the linear growth rate
        k = 2 * np.pi * 2 / N
        t_end = 6.0 / kh_linear_rate(k, 1.0, N / 32.0, 2.0)
        pdir, odir, meta = ctx.case(name, gen.gen_kh, N=N, mode=2, a_shear=N / 32.0,
                                    A=1e-2, nout=24, t_end=t_end)
        cases.append((name, pdir, odir))
        metas[name] = meta
    walls = ctx.run_many(cases)
    res = {}
    for name, pdir, odir in cases:
        m = metas[name]
        r = analyze_kh(odir, m)
        r["rate_linear"] = kh_linear_rate(m["k"], m["dv"], m["a_shear"], m["c0"])
        r["rate_sharp"] = m["k"] * m["dv"] / 2
        r["wall"] = walls[name]
        res[name] = r
        print(f"  {name}: growth rate {r['rate']:.4f}, linear theory (tanh layer, incompressible) "
              f"{r['rate_linear']:.4f}, sharp interface {r['rate_sharp']:.4f}  {r['wall']:.1f}s")
    out = {"cases": res}
    tol.check_kh(ctx, out, tier)
    return out


TESTS = {"sedov": test_sedov, "sound": test_sound, "pancake": test_pancake,
         "cosmobox": test_cosmobox, "kh": test_kh}
