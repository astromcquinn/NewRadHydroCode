#!/usr/bin/env python3
"""Analysis of the hydro accuracy suite runs (docs/hydro_tests.md).

Library of per-test analysis functions used by run_suite.py, plus the
figure makers (`python3 analyze.py figures <results_dir>` re-draws all
figures from the stored result JSON files).
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hs_common as hc  # noqa: E402
import generate as gen  # noqa: E402


# ======================================================================
# Shock tubes
# ======================================================================
def tube_profile(state, meta):
    """Average the 3-D state over the cells of equal normal coordinate.
    Returns xi (cell centre, dimensionless), rho, u_n (lab), p, and the
    maximal spread within a group (relative)."""
    N = meta["N"]
    m, n, Ln, s, off = gen.tube_geometry(meta["direction"], N)
    rho, vel, p, info = state
    un = sum(vel[a] * n[a] for a in range(3))
    keys, inv = np.unique(np.mod(s, N), return_inverse=True)
    inv = inv.ravel()
    cnt = np.bincount(inv)
    prof = []
    spread = 0.0
    for q in (rho, un, p):
        qf = q.ravel()
        mean = np.bincount(inv, weights=qf) / cnt
        mx = np.full(keys.size, -np.inf)
        mn = np.full(keys.size, np.inf)
        np.maximum.at(mx, inv, qf)
        np.minimum.at(mn, inv, qf)
        spread = max(spread, float(np.max((mx - mn) / (np.abs(mean).max() + 1e-300))))
        prof.append(mean)
    return keys / N, prof[0], prof[1], prof[2], spread, off / N


def tube_exact(meta, t, xi_c, off):
    """Exact solution (dimensionless rho, u_n, p) averaged over the cells
    with centres xi_c (comoving frame shift removed by the caller)."""
    L, R = meta["L"], meta["R"]
    U, Ln = meta["U"], meta["Ln"]
    td = t * U / Ln   # dimensionless time

    def fun(x):
        x = np.mod(x, 1.0)
        out = [np.empty_like(x) for _ in range(3)]
        # interface 0.5 (L|R) owns [0.25, 0.75), interface 0 (R|L) the rest
        a = (x >= 0.25) & (x < 0.75)
        r1 = hc.riemann_sample(*L, *R, (x[a] - 0.5) / td)
        xb = np.where(x[~a] >= 0.75, x[~a] - 1.0, x[~a])
        r2 = hc.riemann_sample(*R, *L, xb / td)
        for q in range(3):
            out[q][a] = r1[q]
            out[q][~a] = r2[q]
        return out

    xs = xi_c[:, None] + off[None, :]
    vals = fun(xs.ravel())
    # average the conserved variables, then convert back
    r = vals[0].reshape(xs.shape)
    u = vals[1].reshape(xs.shape)
    p = vals[2].reshape(xs.shape)
    g = hc.GAMMA
    R_ = r.mean(axis=1)
    M_ = (r * u).mean(axis=1)
    E_ = (p / (g - 1) + 0.5 * r * u * u).mean(axis=1)
    U_ = M_ / R_
    P_ = (g - 1) * (E_ - 0.5 * R_ * U_ * U_)
    return R_, U_, P_


def crossing(x, y, level, lo, hi, direction=-1):
    """first position in [lo,hi] where y crosses `level` (linear
    interpolation between samples), scanning in +x; NaN if none."""
    sel = np.where((x >= lo) & (x <= hi))[0]
    for k in range(len(sel) - 1):
        i, j = sel[k], sel[k + 1]
        if (y[i] - level) * (y[j] - level) <= 0 and y[i] != y[j]:
            return x[i] + (level - y[i]) / (y[j] - y[i]) * (x[j] - x[i])
    return float("nan")


def analyze_tube(run, meta, out_index=-1):
    N = meta["N"]
    outs = hc.outputs(run)
    d = outs[out_index]
    state = hc.read_state(d, N)
    info = state[3]
    t = info["t"]
    xi, rho, un, p, spread, off = tube_profile(state, meta)
    U, Ln = meta["U"], meta["Ln"]
    vb_n = float(np.dot(meta["bulk"], meta["normal"]))
    # comoving frame: shift by the bulk motion along the normal
    shift = vb_n * t / Ln
    xi_c = np.mod(xi - shift, 1.0)
    order = np.argsort(xi_c)
    xi_c = xi_c[order]
    r_n = rho[order] / hc.DGAS_AVG
    u_n = (un[order] - vb_n) / U
    p_n = p[order] / (hc.DGAS_AVG * U * U)
    r_e, u_e, p_e = tube_exact(meta, t, xi_c, off)
    L, R = meta["L"], meta["R"]
    cL = math.sqrt(hc.GAMMA * L[2] / L[0])
    l1 = {"rho": float(np.mean(np.abs(r_n - r_e)) / L[0]),
          "u": float(np.mean(np.abs(u_n - u_e)) / cL),
          "p": float(np.mean(np.abs(p_n - p_e)) / L[2])}
    main = (xi_c >= 0.25) & (xi_c < 0.75)     # the L|R problem at 0.5
    l1_main = {"rho": float(np.mean(np.abs(r_n - r_e)[main]) / L[0]),
               "u": float(np.mean(np.abs(u_n - u_e)[main]) / cL),
               "p": float(np.mean(np.abs(p_n - p_e)[main]) / L[2])}
    # wave positions of the L|R problem at 0.5 (in units of the sampling
    # spacing along the normal = 1/N in xi)
    w = hc.riemann_waves(*L, *R)
    td = t * U / Ln
    res = {"t": t, "N": N, "l1": l1, "l1_main": l1_main, "spread": spread,
           "min_rho_main": float(r_n[main].min()), "min_p_main": float(p_n[main].min()), "nsteps": int(info["istep"]),
           "min_rho": float(r_n.min()), "min_p": float(p_n.min()),
           "exact_min_rho": float(r_e.min()), "exact_min_p": float(p_e.min())}
    pos = {}
    if meta["kind"] in ("sod", "strong"):
        rsl, rsr = w["rsl"], w["rsr"]
        lv = {"rarefaction": ((L[0] + rsl) / 2, 0.5 + 0.5 * (w["lhead"] + w["ltail"]) * td),
              "contact": ((rsl + rsr) / 2, 0.5 + w["us"] * td),
              "shock": ((rsr + R[0]) / 2, 0.5 + w["rhead"] * td)}
        lo, hi = 0.5 + w["lhead"] * td - 2.0 / N, 0.5 + w["rhead"] * td + 3.0 / N
        # restrict each search to its own wave region
        bounds = {"rarefaction": (lo, 0.5 + (w["ltail"] + w["us"]) / 2 * td),
                  "contact": (0.5 + (w["ltail"] + w["us"]) / 2 * td,
                              0.5 + (w["us"] + w["rhead"]) / 2 * td),
                  "shock": (0.5 + (w["us"] + w["rhead"]) / 2 * td, hi)}
        for name, (level, x_ex) in lv.items():
            a, b = bounds[name]
            xn = crossing(xi_c, r_n, level, a, b)
            pos[name] = {"num": xn * N, "exact": x_ex * N, "err_cells": (xn - x_ex) * N}
        # widths (10%-90%) of the contact and the shock, in cells along the normal
        def width(a, b, hi_v, lo_v):
            x10 = crossing(xi_c, r_n, lo_v + 0.9 * (hi_v - lo_v), a, b)
            x90 = crossing(xi_c, r_n, lo_v + 0.1 * (hi_v - lo_v), a, b)
            return (x90 - x10) * N
        res["contact_width"] = width(*bounds["contact"], rsl, rsr)
        res["shock_width"] = width(*bounds["shock"], rsr, R[0])
        # overshoots / undershoots near the shock
        res["max_rho_over"] = float(r_n.max() / max(rsr, L[0]) - 1)
    res["positions"] = pos
    res["profile"] = {"xi": xi_c.tolist(), "rho": r_n.tolist(), "u": u_n.tolist(),
                      "p": p_n.tolist(), "rho_ex": r_e.tolist(), "u_ex": u_e.tolist(),
                      "p_ex": p_e.tolist()}
    # conservation (lab frame)
    lg = hc.read_log(run)
    i0 = hc.read_info(outs[0])
    res["dmass"] = float(lg["mass"][-1] / i0["mass"] - 1)
    res["denergy"] = float(lg["energy"][-1] / i0["energy"] - 1)
    p0 = np.array(i0["momentum"])
    res["dmom"] = float(np.abs(np.array([lg["px"][-1], lg["py"][-1], lg["pz"][-1]]) - p0).max()
                        / (i0["mass"] * cL * U))
    return res


def compare_profiles(res_a, res_b):
    """max |difference| of two tube profiles (dimensionless, same grid)"""
    pa, pb = res_a["profile"], res_b["profile"]
    out = {}
    for q in ("rho", "u", "p"):
        a, b = np.array(pa[q]), np.array(pb[q])
        out[q] = float(np.max(np.abs(a - b)))
        out[q + "_l1"] = float(np.mean(np.abs(a - b)))
    return out


def fit_order(Ns, errs):
    Ns, errs = np.asarray(Ns, float), np.asarray(errs, float)
    ok = np.isfinite(errs) & (errs > 0)
    if ok.sum() < 2:
        return float("nan")
    return float(-np.polyfit(np.log(Ns[ok]), np.log(errs[ok]), 1)[0])


# ======================================================================
# Sedov-Taylor
# ======================================================================
def fibonacci_dirs(n):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    th = np.pi * (1 + 5 ** 0.5) * i
    return np.stack([np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)], axis=1)


def trilinear_periodic(f, pts):
    """f[N,N,N] on cell centres (i+1/2); pts (...,3) in cells"""
    N = f.shape[0]
    x = pts - 0.5
    i0 = np.floor(x).astype(int)
    w = x - i0
    out = 0
    for dx in (0, 1):
        for dy in (0, 1):
            for dz in (0, 1):
                ww = ((w[..., 0] if dx else 1 - w[..., 0]) * (w[..., 1] if dy else 1 - w[..., 1])
                      * (w[..., 2] if dz else 1 - w[..., 2]))
                out = out + ww * f[(i0[..., 0] + dx) % N, (i0[..., 1] + dy) % N,
                                   (i0[..., 2] + dz) % N]
    return out


def outer_crossing(x, y, lvl):
    """outermost x where y >= lvl, interpolated to the crossing"""
    k = np.where(y >= lvl)[0]
    if k.size == 0:
        return float("nan")
    k = k.max()
    if k + 1 >= len(x):
        return float("nan")
    return float(x[k] + (y[k] - lvl) / (y[k] - y[k + 1]) * (x[k + 1] - x[k]))


def analyze_sedov(run, meta):
    N = meta["N"]
    rho0, E, p0 = meta["rho0"], meta["E"], meta["p0"]
    g = hc.GAMMA
    sp = hc.sedov_profile()
    xi0 = sp["xi0"]
    outs = hc.outputs(run)
    i0 = hc.read_info(outs[0])
    res = {"N": N, "outputs": []}
    dirs = fibonacci_dirs(400)
    axes = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]], float)
    diag = np.array([[a, b, c] for a in (1, -1) for b in (1, -1) for c in (1, -1)],
                    float) / 3 ** 0.5
    c = np.arange(N) + 0.5
    X = np.meshgrid(c, c, c, indexing="ij")
    for d in outs[1:]:
        rho, vel, p, info = hc.read_state(d, N)
        t = info["t"]
        ctr = np.mod(np.array(meta["centre"]) + np.array(meta["bulk"]) * t, N)
        dx = [np.mod(X[a] - ctr[a] + N / 2, N) - N / 2 for a in range(3)]
        r = np.sqrt(sum(q ** 2 for q in dx))
        vr = sum((vel[a] - meta["bulk"][a]) * dx[a] for a in range(3)) / np.maximum(r, 1e-12)
        R_ex = xi0 * (E * t ** 2 / rho0) ** 0.2
        Rdot = 0.4 * R_ex / t
        nb = N
        edges = np.arange(nb + 1) * 0.5
        idx = np.digitize(r.ravel(), edges) - 1
        ok = idx < nb
        cnt = np.bincount(idx[ok], minlength=nb)

        def binned(q):
            return np.bincount(idx[ok], weights=q.ravel()[ok], minlength=nb) / np.maximum(cnt, 1)
        rc = binned(r)
        good = cnt > 0
        lam = r / R_ex
        rho2 = rho0 * (g + 1) / (g - 1)
        v2 = 2 / (g + 1) * Rdot
        p2 = 2 / (g + 1) * rho0 * Rdot ** 2
        inside = lam < 1
        ex_rho = np.where(inside, np.interp(lam, sp["lam"], sp["rho"]) * rho2, rho0)
        ex_p = np.where(inside, np.interp(lam, sp["lam"], sp["p"]) * p2, p0)
        ex_v = np.where(inside, np.interp(lam, sp["lam"], sp["v"]) * v2, 0.0)
        sel = lam < 1.2
        l1 = {"rho": float(np.abs(rho - ex_rho)[sel].mean() / rho2),
              "p": float(np.abs(p - ex_p)[sel].mean() / p2),
              "v": float(np.abs(vr - ex_v)[sel].mean() / v2)}
        rs = np.arange(0, 0.49 * N, 0.25)

        def ray_radius(D):
            pts = ctr[None, None, :] + D[:, None, :] * rs[None, :, None]
            f = trilinear_periodic(rho, pts)
            return np.array([outer_crossing(rs, row, 0.5 * (rho0 + row.max())) for row in f])
        Rd = ray_radius(dirs)
        Ra = ray_radius(axes)
        Rg = ray_radius(diag)
        prho = binned(rho)
        pr = np.where(good, prho, rho0)
        R_prof = outer_crossing(rc, pr, 0.5 * (rho0 + pr.max()))
        o = {"t": t, "nsteps": int(info["istep"]), "R_exact": R_ex,
             "R_mean_dirs": float(np.nanmean(Rd)), "R_profile": R_prof,
             "R_err_rel": float(np.nanmean(Rd) / R_ex - 1),
             "aniso_rms": float(np.nanstd(Rd) / np.nanmean(Rd)),
             "aniso_maxmin": float((np.nanmax(Rd) - np.nanmin(Rd)) / np.nanmean(Rd)),
             "R_axes": Ra.tolist(), "R_diag": Rg.tolist(),
             "axes_over_diag": float(np.nanmean(Ra) / np.nanmean(Rg) - 1),
             "peak_rho_over_rho2": float(pr.max() / rho2),
             "l1": l1, "wall": info["wall"],
             "denergy": float(info["energy"] / i0["energy"] - 1),
             "dmass": float(info["mass"] / i0["mass"] - 1),
             "profile": {"r": rc[good].tolist(), "rho": prho[good].tolist(),
                         "p": binned(p)[good].tolist(), "v": binned(vr)[good].tolist(),
                         "rho_ex": binned(ex_rho)[good].tolist(),
                         "p_ex": binned(ex_p)[good].tolist(),
                         "v_ex": binned(ex_v)[good].tolist(), "rho2": rho2, "p2": p2,
                         "v2": v2, "R": R_ex},
             "R_dirs": Rd.tolist()}
        res["outputs"].append(o)
    ts = np.array([o["t"] for o in res["outputs"]])
    Rs = np.array([o["R_mean_dirs"] for o in res["outputs"]])
    if len(ts) >= 2:
        res["R_slope"] = float(np.polyfit(np.log(ts), np.log(Rs), 1)[0])
    res["dirs"] = dirs.tolist()
    return res


# ======================================================================
# Zel'dovich pancake
# ======================================================================
def shock_halfwidth(x, K, thr=10.0):
    """distance from the centre (x = 0.5) to the first point, walking
    outwards on each side, where K drops below thr (linear in log K);
    returns (left, right) in units of the wavelength."""
    out = []
    for side in (-1, 1):
        sel = np.where(side * (x - 0.5) >= 0)[0]
        d = np.abs(x[sel] - 0.5)
        o = np.argsort(d)
        d, k = d[o], np.log(np.maximum(K[sel][o], 1e-30))
        if k[0] < math.log(thr):
            out.append(0.0)
            continue
        j = np.where(k < math.log(thr))[0]
        if j.size == 0:
            out.append(float("nan"))
            continue
        j = j[0]
        f = (k[j - 1] - math.log(thr)) / (k[j - 1] - k[j])
        out.append(float(d[j - 1] + f * (d[j] - d[j - 1])))
    return out


def analyze_pancake(run, meta, ref=None):
    import cosmo_gen as cg
    sys.path.insert(0, os.path.join(hc.ROOT, "tools", "ics"))
    import pmcosmo as pc
    N = meta["N"]
    fb = meta["fb"]
    G = pc.Growth(1.0)
    a_i, a_c = 1 / (1 + meta["z_init"]), 1 / (1 + meta["z_c"])
    K0 = meta["T_init"] / cg.units(a_i, meta["box"], N)["temp1_unit"] / fb ** (2 / 3)
    x = (np.arange(N) + 0.5) / N
    outs = hc.outputs(run)
    res = {"N": N, "outputs": []}
    refz = {}
    if ref is not None:
        for i in range(int(ref["n"])):
            refz[round(float(ref[f"{i}_z"]), 6)] = i
    for d in outs[1:]:
        rho3, vel3, p3, info = hc.read_state(d, N)
        z, a = info["redshift"], info["scalefactor"]
        rho = rho3.mean(axis=(1, 2))
        v = vel3[0].mean(axis=(1, 2))
        eps = (p3 / ((hc.GAMMA - 1) * rho3)).mean(axis=(1, 2))
        spread = float(np.ptp(rho3, axis=(1, 2)).max() / rho.max())
        T = eps * info["temp1_unit"]
        K = eps / rho ** (2 / 3) / K0
        o = {"z": z, "a": a, "nsteps": int(info["istep"]), "wall": info["wall"],
             "spread": spread, "rho_max": float(rho.max() / fb), "T_max": float(T.max()),
             "profile": {"x": x.tolist(), "rho": (rho / fb).tolist(), "v": v.tolist(),
                         "T": T.tolist(), "K": K.tolist()}}
        if z > meta["z_c"] + 1e-6:
            A = G(a) / G(a_c)
            q = cg.zeldovich_1d(x.copy(), A, 2 * np.pi, 0.5)
            cq = np.cos(2 * np.pi * (q - 0.5))
            rho_ex = fb / (1 - A * cq)
            v_ex = -G.dtau(a) / G(a_c) * np.sin(2 * np.pi * (q - 0.5)) / (2 * np.pi) * N
            T_ex = K0 * rho_ex ** (2 / 3) * info["temp1_unit"]
            o["zeldovich"] = {"A": float(A),
                              "l1_rho": float(np.mean(np.abs(rho - rho_ex)) / np.mean(rho_ex)),
                              "max_rho": float(np.max(np.abs(rho / rho_ex - 1))),
                              "l1_v": float(np.mean(np.abs(v - v_ex)) / np.abs(v_ex).max()),
                              "max_v": float(np.max(np.abs(v - v_ex)) / np.abs(v_ex).max()),
                              "l1_T": float(np.mean(np.abs(T / T_ex - 1))),
                              "max_T": float(np.max(np.abs(T / T_ex - 1)))}
            o["profile"].update(rho_ex=(rho_ex / fb).tolist(), v_ex=v_ex.tolist(),
                                T_ex=T_ex.tolist())
            o["K_dev_max"] = float(np.max(np.abs(K - 1)))
        # shock (entropy jump) position and entropy away from the shock
        hw = shock_halfwidth(x, K)
        o["shock_halfwidth"] = hw
        xs = np.nanmean(hw)
        away = np.abs(x - 0.5) > xs + 3.0 / N
        Tfloor = 1.0 * 1.0001
        unfloored = T > 1.05 * Tfloor
        if away.any():
            o["K_away_max_dev"] = float(np.max(np.abs(K[away & unfloored] - 1))) \
                if (away & unfloored).any() else float("nan")
            o["K_away_median_dev"] = float(np.median(np.abs(K[away] - 1)))
            o["frac_floor_away"] = float(np.mean(~unfloored[away]))
        if ref is not None and round(z, 6) in refz:
            i = refz[round(z, 6)]
            xr, rr = ref[f"{i}_x"], ref[f"{i}_rho"] / fb
            Kr, Tr, vr = ref[f"{i}_K"], ref[f"{i}_T"], ref[f"{i}_v"] * N
            hw_r = shock_halfwidth(xr, Kr)
            o["ref"] = {"shock_halfwidth": hw_r, "rho_max": float(rr.max()),
                        "T_max": float(Tr.max())}
            # reference averaged over the code's cells (mass-weighted rho)
            edges = np.arange(N + 1) / N
            idx = np.clip(np.digitize(xr, edges) - 1, 0, N - 1)
            dm = np.full(xr.size, fb / xr.size)
            mass = np.bincount(idx, weights=dm, minlength=N)
            rho_c = mass * N / fb
            o["profile"]["rho_ref_cell"] = rho_c.tolist()
            o["profile"]["ref_x"] = xr[::4].tolist()
            o["profile"]["ref_rho"] = rr[::4].tolist()
            o["profile"]["ref_T"] = Tr[::4].tolist()
            o["profile"]["ref_v"] = vr[::4].tolist()
            o["profile"]["ref_K"] = Kr[::4].tolist()
            o["ref"]["l1_rho_cell"] = float(np.mean(np.abs(rho / fb - rho_c)) / np.mean(rho_c))
            o["shock_err_cells"] = float((xs - np.nanmean(hw_r)) * N)
        res["outputs"].append(o)
    return res
