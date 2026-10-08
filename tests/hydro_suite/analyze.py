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
            # the same crossing on the cell-averaged exact profile
            xe = crossing(xi_c, r_e, level, a, b)
            pos[name] = {"num": xn * N, "exact": xe * N, "err_cells": (xn - xe) * N}
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
def shock_halfwidth(x, K, thr=1e3):
    """distance from the centre (x = 1/2) of the accretion shocks: walking
    inwards from the void (x = 0 or 1) on each side, the first point where
    the entropy K/K_init exceeds thr (linear in log K). The centre itself
    holds the cold, adiabatically compressed gas of the first shell
    crossing, so the search must start outside. Units of the wavelength."""
    out = []
    lt = np.log(thr)
    for side in (-1, 1):
        sel = np.where(side * (x - 0.5) >= 0)[0]
        d = np.abs(x[sel] - 0.5)
        o = np.argsort(-d)                     # outside -> inside
        d, k = d[o], np.log(np.maximum(K[sel][o], 1e-30))
        j = np.where(k > lt)[0]
        if j.size == 0:
            out.append(float("nan"))
            continue
        j = j[0]
        if j == 0:
            out.append(float(d[0]))
            continue
        f = (lt - k[j - 1]) / (k[j] - k[j - 1])
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
        # accretion shock position: walking in from the void, the first point
        # where T exceeds 10 % of the profile maximum (the same rule is
        # applied to the reference; an entropy threshold cannot be used for
        # the code, whose pre-shock entropy is spuriously raised, see docs)
        hw = shock_halfwidth(x, T / T.max(), thr=0.1)
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
            hw_r = shock_halfwidth(xr, Tr / Tr.max(), thr=0.1)
            o["ref"] = {"shock_halfwidth": hw_r, "rho_max": float(rr.max()),
                        "T_max": float(Tr.max())}
            # reference averaged over the code's cells (mass-weighted rho)
            edges = np.arange(N + 1) / N
            idx = np.clip(np.digitize(xr, edges) - 1, 0, N - 1)
            dm = np.full(xr.size, fb / xr.size)
            mass = np.bincount(idx, weights=dm, minlength=N)
            rho_c = mass * N / fb
            o["profile"]["rho_ref_cell"] = rho_c.tolist()
            o["profile"]["ref_x"] = xr[::8].tolist()
            o["profile"]["ref_rho"] = rr[::8].tolist()
            o["profile"]["ref_T"] = Tr[::8].tolist()
            o["profile"]["ref_v"] = vr[::8].tolist()
            o["profile"]["ref_K"] = Kr[::8].tolist()
            o["ref"]["l1_rho_cell"] = float(np.mean(np.abs(rho / fb - rho_c)) / np.mean(rho_c))
            if z <= meta["z_c"] + 1e-6:  # an accretion shock exists only after the caustic
                o["shock_err_cells"] = float((xs - np.nanmean(hw_r)) * N)
        res["outputs"].append(o)
    return res


# ======================================================================
# Sound waves
# ======================================================================
def analyze_sound(run, meta):
    """Phase speed, amplitude and L1 error of the travelling sound wave."""
    N = meta["N"]
    m, n, Ln, s, off = gen.tube_geometry(meta["direction"], N)
    keys, inv = np.unique(np.mod(s, N), return_inverse=True)
    inv = inv.ravel()
    cnt = np.bincount(inv)
    xi = keys / N
    A, c0 = meta["A"], meta["c0"]
    rho0 = hc.DGAS_AVG
    vb_n = float(np.dot(meta["bulk"], meta["normal"]))
    ts, ph, amp, l1, k1err = [], [], [], [], []
    # exact (cell-averaged) wave shape
    sub = off / N
    for d in hc.outputs(run):
        rho, vel, p, info = hc.read_state(d, N)
        t = info["t"]
        prof = np.bincount(inv, weights=rho.ravel(), minlength=keys.size) / cnt
        drho = prof / rho0 - 1
        F = np.sum(drho * np.exp(-2j * np.pi * xi)) * 2 / xi.size
        ts.append(t)
        ph.append(np.angle(F))
        amp.append(abs(F))
        shift = (c0 + vb_n) * t / Ln
        ex = A * np.sin(2 * np.pi * (xi[:, None] + sub[None, :] - shift)).mean(axis=1)
        l1.append(float(np.mean(np.abs(drho - ex)) / A))
        Fex = np.sum(ex * np.exp(-2j * np.pi * xi)) * 2 / xi.size
        k1err.append(float(abs(F - Fex) / A))
    ts = np.array(ts)
    ph = np.unwrap(np.array(ph))
    # the k = 1 mode sin(2 pi (xi - s)) has phase -pi/2 - 2 pi s
    speed = -np.polyfit(ts, ph, 1)[0] / (2 * np.pi) * Ln
    amp = np.array(amp)
    # amplitude of the k=1 mode of the cell-averaged exact profile
    a_ex = abs(np.sum(np.sin(2 * np.pi * (xi[:, None] + sub[None, :])).mean(axis=1)
                      * np.exp(-2j * np.pi * xi)) * 2 / xi.size) * A
    T = meta["period"]
    damping = 1 - amp[-1] / a_ex
    per_period = 1 - (amp[-1] / a_ex) ** (T / ts[-1])
    info = hc.read_info(hc.outputs(run)[-1])
    i0 = hc.read_info(hc.outputs(run)[0])
    return {"N": N, "direction": meta["direction"], "t": ts.tolist(), "amp": amp.tolist(),
            "speed": float(speed), "speed_exact": c0 + vb_n,
            "speed_err": float((speed - vb_n) / c0 - 1),
            "damping_per_period": float(per_period), "damping_total": float(damping),
            "l1": l1, "l1_end": l1[-1], "k1_err": k1err, "k1_err_end": k1err[-1],
            "nsteps": int(info["istep"]),
            "wall": info["wall"], "periods": ts[-1] / T,
            "denergy": float(info["energy"] / i0["energy"] - 1)}


# ======================================================================
# Figures
# ======================================================================
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#4a3aa7", "#e87ba4", "#008300", "#e34948"]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#d9d8d4"


def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 8, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                         "xtick.color": INK2, "ytick.color": INK2, "axes.grid": True,
                         "grid.color": GRID, "grid.linewidth": 0.5, "lines.linewidth": 1.2,
                         "legend.frameon": False, "axes.spines.top": False,
                         "axes.spines.right": False, "savefig.dpi": 110})
    return plt


def _save(fig, figdir, name):
    os.makedirs(figdir, exist_ok=True)
    p = os.path.join(figdir, f"hydro_suite_{name}.png")
    fig.savefig(p, bbox_inches="tight")
    print("wrote", p)


def fig_tubes(R, figdir):
    plt = _plt()
    cs = R["cases"]
    fig, ax = plt.subplots(3, 3, figsize=(9.5, 6.5))
    cols = [("sod", "Sod (x)"), ("strong", "strong shock, p ratio 1e5 (x)"),
            ("123", "123 double rarefaction (x)")]
    for j, (kind, title) in enumerate(cols):
        for i, q in enumerate(("rho", "u", "p")):
            a = ax[i, j]
            first = True
            for n, N in enumerate((32, 64, 128)):
                c = cs.get(f"tube_{kind}_x_N{N}")
                if not c:
                    continue
                pr = c["profile"]
                xi = np.array(pr["xi"])
                if first:
                    a.plot(xi, pr[q + "_ex"], color=INK2, lw=1.0, label="exact")
                    first = False
                a.plot(xi, pr[q], ".", ms=2.5 if N == 128 else 3.5, color=C[n], label=f"N={N}")
            a.set_xlim(0, 1)
            if q == "p" and kind == "strong":
                a.set_yscale("log")
            if j == 0:
                a.set_ylabel({"rho": r"$\rho/\rho_L$", "u": r"$u/U$", "p": r"$p/p_L$"}[q])
            if i == 0:
                a.set_title(title, fontsize=8.5)
            if i == 2:
                a.set_xlabel(r"$x/L$ (comoving frame)")
    ax[0, 0].legend(loc="upper right", fontsize=7)
    fig.suptitle("Shock tubes at the final time (dimensionless; two interfaces, at 0 and 1/2)",
                 fontsize=9)
    fig.tight_layout()
    _save(fig, figdir, "tubes")
    plt.close(fig)

    # convergence
    fig, ax = plt.subplots(1, 3, figsize=(9.5, 3.0))
    marks = {"x": "o", "y": "s", "z": "^", "xy": "D", "xyz": "v"}
    for k, (kind, d) in enumerate([("sod", d) for d in ("x", "y", "z", "xy", "xyz")]):
        Ns = [N for N in (32, 64, 128) if f"tube_{kind}_{d}_N{N}" in cs]
        e = [cs[f"tube_{kind}_{d}_N{N}"]["l1"]["rho"] for N in Ns]
        ax[0].loglog(Ns, e, marker=marks[d], color=C[k], label=f"Sod {d}", ms=4)
    for k, kind in enumerate(("strong", "123")):
        Ns = [N for N in (32, 64, 128) if f"tube_{kind}_x_N{N}" in cs]
        for q, ls in (("rho", "-"),):
            e = [cs[f"tube_{kind}_x_N{N}"]["l1"][q] for N in Ns]
            ax[0].loglog(Ns, e, ls, marker="o", color=C[5 + k], label=f"{kind} x", ms=4)
    Nr = np.array([32, 128])
    ax[0].loglog(Nr, 0.04 * (Nr / 32.0) ** -1, ":", color=INK2, lw=0.8)
    ax[0].loglog(Nr, 0.3 * (Nr / 32.0) ** -0.5, ":", color=INK2, lw=0.8)
    ax[0].text(100, 0.04 * (100 / 32) ** -1, r"$N^{-1}$", color=INK2)
    ax[0].text(100, 0.3 * (100 / 32) ** -0.5, r"$N^{-1/2}$", color=INK2)
    ax[0].set_xlabel("N (cells per period)")
    ax[0].set_ylabel(r"L1($\rho$)/$\rho_L$")
    ax[0].set_title("density L1 error", fontsize=8.5)
    ax[0].legend(fontsize=6.5, ncol=2)
    for k, d in enumerate(("x", "y", "z", "xy", "xyz")):
        Ns = [N for N in (64, 128) if f"tube_sod_{d}_N{N}" in cs]
        ax[1].plot(Ns, [cs[f"tube_sod_{d}_N{N}"].get("contact_width", np.nan) for N in Ns],
                   marker=marks[d], color=C[k], label=f"contact {d}")
        ax[1].plot(Ns, [cs[f"tube_sod_{d}_N{N}"].get("shock_width", np.nan) for N in Ns],
                   "--", marker=marks[d], color=C[k], ms=3)
    ax[1].set_xscale("log", base=2)
    ax[1].set_xlabel("N")
    ax[1].set_ylabel("10-90% width [cells along the normal]")
    ax[1].set_title("Sod: contact (solid) and shock (dashed) width", fontsize=8.5)
    ax[1].legend(fontsize=6.5)
    for k, d in enumerate(("x", "y", "z", "xy", "xyz")):
        Ns = [N for N in (32, 64, 128) if f"tube_sod_{d}_N{N}" in cs]
        for w, ls in (("shock", "-"), ("contact", "--")):
            ax[2].plot(Ns, [cs[f"tube_sod_{d}_N{N}"]["positions"][w]["err_cells"] for N in Ns],
                       ls, marker=marks[d], color=C[k], ms=4,
                       label=f"{d}" if w == "shock" else None)
    ax[2].axhline(0, color=INK2, lw=0.6)
    ax[2].set_xscale("log", base=2)
    ax[2].set_xlabel("N")
    ax[2].set_ylabel("position error [cells]")
    ax[2].set_title("Sod: shock (solid) / contact (dashed) position", fontsize=8.5)
    ax[2].legend(fontsize=6.5)
    fig.tight_layout()
    _save(fig, figdir, "tube_convergence")
    plt.close(fig)


def fig_galilean(R, figdir):
    plt = _plt()
    cs = R["cases"]
    fig, ax = plt.subplots(1, 2, figsize=(9.5, 3.0))
    names = sorted(cs, key=lambda n: (cs[n]["bulk"][0], abs(cs[n]["bulk"][2])))
    for k, n in enumerate(names):
        pr = cs[n]["profile"]
        b = cs[n]["bulk"]
        lab = f"vx={b[0]:.2f}, |v_t|={np.hypot(b[1], b[2]):.2f}"
        ax[0].plot(pr["xi"], pr["rho"], "-", color=C[k % 8], lw=0.9, label=lab)
        ax[1].plot(pr["xi"], np.array(pr["rho"]) - np.array(cs[names[0]]["profile"]["rho"]),
                   color=C[k % 8], lw=0.9)
    ax[0].plot(pr["xi"], pr["rho_ex"], color=INK2, lw=0.8, ls=":", label="exact")
    ax[0].set_ylabel(r"$\rho/\rho_L$")
    ax[1].set_ylabel(r"$\rho - \rho_{\rm ref}$")
    for a in ax:
        a.set_xlabel(r"$x/L$ (comoving frame)")
    ax[0].legend(fontsize=6)
    ax[0].set_title(f"Sod x, N={R['N']}, different bulk velocities", fontsize=8.5)
    ax[1].set_title("difference to the smallest bulk velocity", fontsize=8.5)
    fig.tight_layout()
    _save(fig, figdir, "galilean")
    plt.close(fig)


def fig_sedov(R, figdir):
    plt = _plt()
    cs = R["cases"]
    fig, ax = plt.subplots(1, 4, figsize=(11, 2.9))
    for k, N in enumerate((32, 64, 128)):
        c = cs.get(f"sedov_N{N}")
        if not c:
            continue
        o = c["outputs"][-1]
        pr = o["profile"]
        r = np.array(pr["r"]) / pr["R"]
        for j, (q, nrm) in enumerate((("rho", "rho2"), ("p", "p2"), ("v", "v2"))):
            if k == 0:
                sp = hc.sedov_profile()
                lam = np.concatenate((sp["lam"], [1.0, 1.0, 1.4]))
                key = {"rho": "rho", "p": "p", "v": "v"}[q]
                rho0 = 1 / ((hc.GAMMA + 1) / (hc.GAMMA - 1))
                ext = {"rho": [1.0, rho0, rho0], "p": [1.0, 0, 0], "v": [1.0, 0, 0]}[q]
                ax[j].plot(lam, np.concatenate((sp[key], ext)), color=INK2, lw=1.0,
                           label="self-similar")
            ax[j].plot(r, np.array(pr[q]) / pr[nrm], ".", ms=2.5, color=C[k], label=f"N={N}")
            ax[j].set_xlim(0, 1.4)
            ax[j].set_xlabel(r"$r/R_{\rm exact}$")
            ax[j].set_ylabel({"rho": r"$\rho/\rho_2$", "p": r"$p/p_2$", "v": r"$v_r/v_2$"}[q])
    ax[0].legend(fontsize=6.5)
    ax[0].set_title("Sedov-Taylor at the last output", fontsize=8.5)
    # anisotropy: R(direction) vs angle to the nearest axis
    for k, N in enumerate((32, 64, 128)):
        c = cs.get(f"sedov_N{N}")
        if not c:
            continue
        o = c["outputs"][-1]
        dirs = np.abs(np.array(c["dirs"]))
        ang = np.degrees(np.arccos(dirs.max(axis=1)))
        ax[3].plot(ang, np.array(o["R_dirs"]) / o["R_exact"], ".", ms=2, color=C[k],
                   label=f"N={N}")
    ax[3].set_xlabel("angle to the nearest grid axis [deg]")
    ax[3].set_ylabel(r"$R(\hat n)/R_{\rm exact}$")
    ax[3].set_title("shock radius by direction", fontsize=8.5)
    ax[3].legend(fontsize=6.5)
    fig.tight_layout()
    _save(fig, figdir, "sedov")
    plt.close(fig)


def fig_sound(R, figdir):
    plt = _plt()
    cs = R["cases"]
    fig, ax = plt.subplots(1, 3, figsize=(9.5, 2.8))
    for k, d in enumerate(("x", "xy")):
        Ns = [N for N in (16, 32, 64, 128) if f"sound_{d}_N{N}" in cs]
        if not Ns:
            continue
        g = lambda key: [abs(cs[f"sound_{d}_N{N}"][key]) for N in Ns]  # noqa
        ax[0].loglog(Ns, g("k1_err_end"), "o-", color=C[k], label=f"k=1 mode, along {d}", ms=4)
        ax[0].loglog(Ns, g("l1_end"), "s--", color=C[k], label=f"L1, along {d}", ms=3, lw=0.8)
        ax[1].loglog(Ns, g("speed_err"), "o-", color=C[k], ms=4)
        ax[2].loglog(Ns, g("damping_per_period"), "o-", color=C[k], ms=4)
    Nr = np.array([16, 128])
    for a, y0 in zip(ax, (0.1, 4e-3, 0.15)):
        a.loglog(Nr, y0 * (Nr / 16.0) ** -2, ":", color=INK2, lw=0.8)
        a.text(64, y0 * 4 ** -2 * 1.3, r"$N^{-2}$", color=INK2)
        a.set_xlabel("cells per wavelength")
    ax[0].set_ylabel("error / A after one period")
    ax[1].set_ylabel("|phase speed error| / c")
    ax[2].set_ylabel("amplitude loss per period")
    ax[0].legend(fontsize=7)
    fig.suptitle("Linear sound wave (A = 1e-3), one period", fontsize=9)
    fig.tight_layout()
    _save(fig, figdir, "sound")
    plt.close(fig)


def fig_pancake(R, figdir):
    plt = _plt()
    cs = R["cases"]
    fig, ax = plt.subplots(2, 3, figsize=(10, 5.6))
    pre_z, post_z = 7.0, 0.0
    for row, zt in enumerate((pre_z, post_z)):
        refdone = False
        for k, name in enumerate(n for n in ("pancake_N32_T1000", "pancake_N64_T1000",
                                             "pancake_N128_T1000") if n in cs):
            c = cs[name]
            o = [o for o in c["outputs"] if abs(o["z"] - zt) < 1e-6]
            if not o:
                continue
            pr = o[0]["profile"]
            x = np.array(pr["x"]) - 0.5
            if not refdone:
                if row == 0 and "rho_ex" in pr:
                    xs = x
                    ax[row, 0].plot(xs, pr["rho_ex"], color=INK2, lw=1, label="Zel'dovich")
                    ax[row, 1].plot(xs, np.array(pr["v_ex"]) / c["N"], color=INK2, lw=1)
                    ax[row, 2].plot(xs, pr["T_ex"], color=INK2, lw=1)
                if "ref_x" in pr:
                    xr = np.array(pr["ref_x"]) - 0.5
                    ax[row, 0].plot(xr, pr["ref_rho"], color=INK, lw=0.7, ls="--",
                                    label="1-D Lagrangian ref.")
                    ax[row, 1].plot(xr, np.array(pr["ref_v"]) / c["N"], color=INK, lw=0.7,
                                    ls="--")
                    ax[row, 2].plot(xr, pr["ref_T"], color=INK, lw=0.7, ls="--")
                refdone = True
            ax[row, 0].plot(x, pr["rho"], ".", ms=3, color=C[k], label=f"N={c['N']}")
            ax[row, 1].plot(x, np.array(pr["v"]) / c["N"], ".", ms=3, color=C[k])
            ax[row, 2].plot(x, pr["T"], ".", ms=3, color=C[k])
        ax[row, 0].set_yscale("log")
        ax[row, 2].set_yscale("log")
        ax[row, 0].set_ylabel(r"$\rho_{\rm gas}/\bar\rho_{\rm gas}$")
        ax[row, 1].set_ylabel(r"$v_x$ [box / code time]")
        ax[row, 2].set_ylabel("T [K]")
        for a in ax[row]:
            a.set_xlim(-0.5, 0.5)
            a.set_title(f"z = {zt:g}" + (" (before the caustic, z_c = 5)" if row == 0 else
                                         " (after)"), fontsize=8.5)
        ax[row, 0].legend(fontsize=6.5)
    for a in ax[1]:
        a.set_xlabel(r"$x/\lambda$")
    ax[1, 2].set_ylim(1e-1, 1e8)
    ax[1, 0].set_ylim(1e-1, 1e3)
    fig.tight_layout()
    _save(fig, figdir, "pancake")
    plt.close(fig)


def fig_cosmobox(R, figdir):
    plt = _plt()
    cs = R["cases"]
    fig, ax = plt.subplots(1, 3, figsize=(10, 3.0))
    main = cs.get("cosmobox_N64_dlna0.02") or cs.get("cosmobox_N64")
    for k, o in enumerate(o for o in main["outputs"] if "pk" in o and o["z"] in (10.0, 3.0, 0.0)):
        pk = o["pk"]
        kk = np.array(pk["k"])
        kn = np.pi * (main["N"] // 2) / 64.0      # particle Nyquist (Ndm = N/2, box 64)
        sel = kk <= kn
        ax[0].loglog(kk[sel], np.array(pk["dm"])[sel], "-", color=C[k], label=f"DM z={o['z']:g}")
        ax[0].loglog(pk["k"], pk["gas"], "--", color=C[k], label=f"gas z={o['z']:g}")
    ax[0].set_xlabel("k [h/Mpc]")
    ax[0].set_ylabel(r"$\Delta^2(k)$")
    ax[0].legend(fontsize=6)
    ax[0].set_title("power spectra", fontsize=8.5)
    o = main["outputs"][-1]
    h = o["trho_hist"]
    H = np.array(h["H"]).T
    ax[1].pcolormesh(h["x"], h["y"], np.log10(np.maximum(H, 0.5)), cmap="Greys", shading="auto")
    tr = o["trho"]
    ax[1].plot(tr["logd"], tr["logT_med"], color=C[0], label="median")
    ld = np.array(tr["logd"])
    tad = np.log10(np.maximum(o["T_adiabatic_mean"] * 10 ** (2 * ld / 3), 1.0))
    ax[1].plot(ld, tad, "--", color=C[1], label="adiabatic (1 K floor)")
    ax[1].legend(fontsize=6.5, loc="upper left")
    ax[1].set_xlabel(r"log $\delta_{\rm gas}$")
    ax[1].set_ylabel("log T [K]")
    ax[1].set_title(f"T-ρ at z = {o['z']:g} (N=64)", fontsize=8.5)
    for k, (name, c) in enumerate(cs.items()):
        zz = [o["z"] for o in c["outputs"]]
        ax[2].semilogy(zz, [o["T_mw"] for o in c["outputs"]], "o-", ms=3, color=C[k],
                       label=name.replace("cosmobox_", ""))
    ax[2].invert_xaxis()
    ax[2].set_xlabel("z")
    ax[2].set_ylabel("mass-weighted T [K]")
    ax[2].legend(fontsize=6.5)
    ax[2].set_title("mass-weighted temperature", fontsize=8.5)
    fig.tight_layout()
    _save(fig, figdir, "cosmobox")
    plt.close(fig)


def fig_kh(R, figdir):
    plt = _plt()
    cs = R["cases"]
    fig, ax = plt.subplots(1, 1, figsize=(4.5, 3.0))
    for k, name in enumerate(sorted(cs, key=lambda n: cs[n]["N"])):
        c = cs[name]
        t = np.array(c["t"]) * c["rate_linear"]
        ax.semilogy(t, c["amp"], "o-", ms=2.5, color=C[k],
                    label=f"N={c['N']} (layer {c['N'] / 32:g} cells)")
    ax.semilogy(t, 0.005 * np.exp(t - 0.5), ":", color=INK2, lw=1, label="linear theory")
    ax.set_xlabel(r"$t\,\sigma_{\rm lin}$")
    ax.set_ylabel(r"$|\hat v_y(k)|$")
    ax.set_title("Kelvin-Helmholtz: growth of the seeded mode", fontsize=8.5)
    ax.legend(fontsize=7)
    fig.tight_layout()
    _save(fig, figdir, "kh")
    plt.close(fig)


def figures(results_dir, figdir, tier="full"):
    for name, fn in (("tubes", fig_tubes), ("galilean", fig_galilean), ("sedov", fig_sedov),
                     ("sound", fig_sound), ("pancake", fig_pancake),
                     ("cosmobox", fig_cosmobox), ("kh", fig_kh)):
        p = os.path.join(results_dir, f"{name}_{tier}.json")
        if os.path.exists(p):
            fn(hc.load_json(p), figdir)


if __name__ == "__main__":
    if len(sys.argv) >= 2 and sys.argv[1] == "figures":
        rd = sys.argv[2] if len(sys.argv) > 2 else os.path.join(hc.HERE, "results")
        fd = sys.argv[3] if len(sys.argv) > 3 else os.path.join(hc.ROOT, "docs", "figures")
        tier = sys.argv[4] if len(sys.argv) > 4 else "full"
        figures(rd, fd, tier)
    else:
        print("usage: analyze.py figures [results_dir] [figure_dir] [tier]")
