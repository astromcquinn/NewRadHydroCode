#!/usr/bin/env python3
"""Empirical entropy of ray fields (docs/ray_compression.md, step 3).

usage: entropy.py SNAPSHOT_STEM [SNAPSHOT_STEM ...] [--json out.json]

Each continuous field is first quantised at a stated precision (the
precision of the structural encoding, see structural.py), then:

  H(X)        plug-in entropy of the quantised field, bits per ray
  H(X|Y)      plug-in conditional entropy, H(X,Y) - H(Y)
  xent(X|Y)   split-half cross-entropy: p(x|y) is counted on even rays
              (add-1/2 smoothing, backing off to p(x) for unseen y) and
              evaluated on odd rays. This is what a coder with a fitted
              table would pay; the plug-in H(X|Y) under-estimates it when
              Y has many values.

It also measures, for rays that were never merged, how exactly direction,
position and dr follow from (origin, nside, ipix, r).
"""
import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import snap  # noqa: E402

try:
    import healpy as hp
except ImportError:  # pragma: no cover
    hp = None


# ---------------------------------------------------------------------------
# entropy estimators
# ---------------------------------------------------------------------------

def _codes(cols):
    """Map rows of the given int columns to dense integer codes."""
    if len(cols) == 0:
        return None
    m = np.stack([np.asarray(c, dtype=np.int64) for c in cols], axis=1)
    _, inv = np.unique(m, axis=0, return_inverse=True)
    return inv.ravel()


def H(*cols):
    c = _codes(cols)
    if c is None or len(c) == 0:
        return 0.0
    p = np.bincount(c) / len(c)
    p = p[p > 0]
    return float(-(p * np.log2(p)).sum())


def Hc(x, y):
    """H(X|Y) for X, Y lists of int columns."""
    if not y:
        return H(*x)
    return H(*(list(x) + list(y))) - H(*y)


def xent(x, y=(), alpha=0.5):
    """Split-half cross-entropy of X given Y (bits/ray)."""
    xc = _codes(x)
    n = len(xc)
    if n < 4:
        return float("nan")
    K = int(xc.max()) + 1
    tr = np.arange(n) % 2 == 0
    te = ~tr
    # marginal with smoothing over the joint alphabet
    cnt = np.bincount(xc[tr], minlength=K).astype(np.float64)
    pm = (cnt + alpha) / (cnt.sum() + alpha * K)
    if not y:
        return float(-np.log2(pm[xc[te]]).mean())
    yc = _codes(y)
    key = yc * K + xc
    ktr, kcnt = np.unique(key[tr], return_counts=True)
    ytr, ycnt = np.unique(yc[tr], return_counts=True)
    # look up counts for the test rows
    pos = np.searchsorted(ktr, key[te])
    pos = np.minimum(pos, len(ktr) - 1)
    cxy = np.where(ktr[pos] == key[te], kcnt[pos], 0).astype(np.float64)
    posy = np.searchsorted(ytr, yc[te])
    posy = np.minimum(posy, len(ytr) - 1)
    cy = np.where(ytr[posy] == yc[te], ycnt[posy], 0).astype(np.float64)
    # Back off to the marginal: Dirichlet prior with mean p(x) and
    # concentration beta = 1, p = (c(x,y) + beta p(x)) / (c(y) + beta)
    beta = 1.0
    p = (cxy + beta * pm[xc[te]]) / (cy + beta)
    return float(-np.log2(p).mean())


# ---------------------------------------------------------------------------
# quantisers (precision of the structural encoding)
# ---------------------------------------------------------------------------

def q_log2(v, step_bits):
    """round(log2 v * 2^step_bits); 0 -> a reserved code far below."""
    v = np.asarray(v, dtype=np.float64)
    out = np.full(v.shape, -(1 << 30), dtype=np.int64)
    pos = v > 0
    out[pos] = np.round(np.log2(v[pos]) * (1 << step_bits)).astype(np.int64)
    return out


def predicted_geometry(a, d, N):
    """For rays with a known source: the image of the source cell centre
    whose distance best matches r, the predicted pixel from the cell centre,
    and the index of the true pixel among the predicted pixel and its 8
    neighbours (9 = not found)."""
    src = a["src"].astype(np.int64)
    ok = src >= 0
    xs = d["src_x"]
    cc = d["ijk"] + 0.5
    r = d["r"]
    base = cc - xs
    base = (base + N / 2) % N - N / 2
    # candidate images within |n| <= nimg box lengths
    nimg = int(np.ceil(np.nanmax(np.where(ok, r, 0)) / N)) + 1
    best = base.copy()
    bestgap = np.abs(np.sqrt((base ** 2).sum(1)) - r)
    rng = range(-nimg, nimg + 1)
    for i in rng:
        for j in rng:
            for k in rng:
                if i == j == k == 0:
                    continue
                v = base + N * np.array([i, j, k], dtype=np.float64)
                g = np.abs(np.sqrt((v ** 2).sum(1)) - r)
                better = g < bestgap
                best[better] = v[better]
                bestgap[better] = g[better]
    dist = np.sqrt((best ** 2).sum(1))
    nside = a["nside"].astype(np.int64)
    ipix = a["ipix"].astype(np.int64)
    pred = np.zeros(len(a), dtype=np.int64)
    vv = best.copy()
    vv[dist == 0] = [0, 0, 1]
    pred = hp.vec2pix(nside, vv[:, 0], vv[:, 1], vv[:, 2], nest=True)
    sym = np.full(len(a), 9, dtype=np.int64)
    sym[pred == ipix] = 0
    nb = hp.get_all_neighbours(nside, pred, nest=True)  # (8, n)
    for m in range(8):
        hit = (nb[m] == ipix) & (sym == 9)
        sym[hit] = m + 1
    sym[~ok] = -1
    return dist, sym, bestgap


# ---------------------------------------------------------------------------

def analyse(stem):
    a, meta = snap.load(stem)
    d = snap.derive(a, meta)
    N = meta["N"]
    n = len(a)
    out = {"stem": str(stem), "nrays": n, "N": N, "redshift": meta["redshift"],
           "substep": meta["substep"], "xHII_vol": meta["xHII_vol"],
           "nsources": meta["nsources"]}
    nm = d["never_merged"]
    cl = d["clamped"]
    exact = nm & ~cl & (a["src"] >= 0)
    out["frac_never_merged"] = float(nm.mean())
    out["frac_clamped"] = float(cl.mean())
    out["frac_origin_exact"] = float(exact.mean())

    # ---- determinism for never-merged rays ----
    det = {}
    if exact.any():
        det["dev_max_rad"] = float(d["dev_abs"][exact].max())
        det["dev_max_pix"] = float(np.max(np.hypot(d["dev_th"], d["dev_ph"])[exact]))
        rc = np.abs(d["x_resid_cen"][exact]).max(1)
        det["x_resid_cen_p50"] = float(np.percentile(rc, 50))
        det["x_resid_cen_p99"] = float(np.percentile(rc, 99))
        det["x_resid_cen_max"] = float(rc.max())
        det["dr_resid_max"] = float(np.abs(d["dr_resid"][exact]).max())
        det["lside_eq_lside_r"] = float((d["lside"] == d["lside_r"])[exact].mean())
    if (nm & cl).any():
        rc = np.abs(d["x_resid_cen"][nm & cl]).max(1)
        det["clamped_x_resid_p50"] = float(np.nanpercentile(rc, 50))
        det["clamped_x_resid_max"] = float(np.nanmax(rc))
    det["dr_resid_max_all"] = float(np.abs(d["dr_resid"]).max())
    det["nu_init_frac"] = [float(v) for v in d["nu_is_init"].mean(0)]
    out["determinism"] = det

    # ---- quantised fields ----
    offq = {b: np.floor(d["off"] * (1 << b)).astype(np.int64) for b in (10, 16)}
    lside = d["lside"]
    lsr = d["lside_r"]
    ipix = a["ipix"].astype(np.int64)
    r8 = np.round(d["r"] * 8).astype(np.int64)
    r1k = np.round(d["r"] * 1024).astype(np.int64)
    dth = np.round(d["dev_th"] * 128).astype(np.int64)
    dph = np.round(d["dev_ph"] * 128).astype(np.int64)
    ntq = q_log2(d["ntot"], 7)  # 2^-7 octave = 0.54 %
    with np.errstate(divide="ignore", invalid="ignore"):
        shq = [q_log2(d["shape"][:, i], 7) for i in range(5)]  # 2^-7 octave
    nuq = [np.round(np.clip(d["nu_pos"][:, i], 0, 1) * 255).astype(np.int64) for i in range(5)]
    taub = np.round(snap.tau_ctx(a) * 4).astype(np.int64)  # 0.25 dex bins of min tau
    cell = a["cell"].astype(np.int64)
    src = a["src"].astype(np.int64)
    birth = a["birth"].astype(np.int64)
    flag = (~nm).astype(np.int64)

    F = {}

    def put(name, raw_bits, x, ctxs, sel=None):
        s = slice(None) if sel is None else sel
        xs = [c[s] for c in x]
        row = {"raw_bits": raw_bits, "n": int(len(xs[0])), "H": H(*xs)}
        for cname, y in ctxs.items():
            ys = [c[s] for c in y]
            row["H|" + cname] = Hc(xs, ys)
            row["xent|" + cname] = xent(xs, ys)
        row["xent"] = xent(xs)
        F[name] = row

    put("merged_flag", 1, [flag], {"lside": [lside], "tau": [taub]})
    put("off_x(16b)", 16, [offq[16][:, 0]], {"off_y": [offq[16][:, 1]]})
    put("off_x(10b)", 10, [offq[10][:, 0]], {"dx_sign": [np.sign(a["dx"][:, 0]).astype(np.int64)]})
    put("lside", 4, [lside], {"lside(r)": [lsr], "lside(r),merged": [lsr, flag]})
    put("ipix", 30, [ipix], {"lside": [lside]})
    put("r(1/8)", 16, [r8], {"lside": [lside]})
    put("dr(1/8)-pred", 16, [np.round(d["dr_resid"] * 8).astype(np.int64)], {})
    put("dev_th(merged)", 8, [dth], {"lside": [lside]}, sel=~nm)
    put("dev_ph(merged)", 8, [dph], {"dev_th": [dth]}, sel=~nm)
    put("dev_th(never)", 8, [dth], {}, sel=nm)
    put("log2 ntot(1/128)", 32, [ntq], {"lside": [lside], "lside,tau": [lside, taub],
                                        "lside,merged": [lside, flag],
                                        "cell": [cell]})
    for i in range(5):
        put("shape%d(1/128)" % i, 12, [shq[i]], {"tau": [taub], "lside,tau": [lside, taub],
                                                  "cell": [cell]})
    put("shape(joint)", 48, shq[:4], {"cell": [cell], "lside,tau": [lside, taub]})
    for i in range(5):
        put("nu%d(8b)" % i, 8, [nuq[i]], {"merged": [flag]})
    put("nu(joint)", 40, nuq, {"merged": [flag], "birth": [birth]})

    # never-merged, exact origin: origin-based conditional entropies
    if exact.sum() > 10:
        dist, sym, gap = predicted_geometry(a, d, N)
        e = exact
        put("src|cell (never)", 30, [src], {"cell": [cell]}, sel=e)
        dr8 = np.round((d["r"] - dist) * 8).astype(np.int64)
        dr1k = np.round((d["r"] - dist) * 1024).astype(np.int64)
        put("r-|v|(1/8) (never)", 16, [dr8], {"lside": [lside]}, sel=e)
        put("r-|v|(1/1024) (never)", 23, [dr1k], {}, sel=e)
        put("ipix|pred (never)", 30, [sym], {"lside": [lside]}, sel=e)
        out["ipix_pred_hit"] = [float(np.mean(sym[e] == k)) for k in range(10)]
        put("log2 ntot (never)", 32, [ntq], {"src,lside,birth": [src, lside, birth]}, sel=e)
        put("shape(joint) (never)", 48, shq[:4], {"src,birth": [src, birth],
                                                  "src,birth,lside": [src, birth, lside]},
            sel=e)
    out["fields"] = F
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stems", nargs="+")
    ap.add_argument("--json")
    args = ap.parse_args()
    res = []
    for st in args.stems:
        r = analyse(st)
        res.append(r)
        print("\n== %s: %d rays, z=%.3f, xHII=%.3f, never-merged %.3f, clamped %.3f" % (
            st, r["nrays"], r["redshift"], r["xHII_vol"], r["frac_never_merged"],
            r["frac_clamped"]))
        print("  determinism:", json.dumps(r["determinism"]))
        if "ipix_pred_hit" in r:
            print("  ipix pred hits:", ["%.3f" % v for v in r["ipix_pred_hit"]])
        for k, v in r["fields"].items():
            extras = "  ".join("%s=%.2f" % (kk, vv) for kk, vv in v.items()
                               if kk not in ("raw_bits", "n", "H", "xent"))
            print("  %-24s n=%-7d raw=%-3d H=%6.2f xent=%6.2f  %s" % (
                k, v["n"], v["raw_bits"], v["H"], v["xent"], extras))
    if args.json:
        with open(args.json, "w") as f:
            json.dump(res, f, indent=1)


if __name__ == "__main__":
    main()
