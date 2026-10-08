#!/usr/bin/env python3
"""Compare two RT output directories written in the phase-1 validation format.

See docs/validation_spec.md and tools/compare/README.md.

Layout expected in REF_DIR and NEW_DIR:

    step_0001/ fHI fHeI fHeII temp rhog       (N^3        float32, Fortran order)
               ngam nu fgam tau               (Nfreq,N^3  float32)
               gamma lambda dx                (3,N^3      float32)
               u_gas (5,N^3)  w_gas (4,N^3)
               stats.txt                       key = value
    step_0002/ ...

A field file may be called ``<name>``, ``<name>.bin``, ``<name>.dat``,
``<name>.raw`` or ``<name>.f32``.  Any other raw file whose size is a multiple
of 4*N^3 bytes is compared as well (as an "extra" field).

Exit status: 0 = all within tolerance, 1 = tolerance exceeded,
2 = structural problem (missing steps/files, size mismatch, bad arguments).
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import sys
from dataclasses import dataclass, field, asdict
from typing import Dict, List, Optional, Tuple

import numpy as np

# --------------------------------------------------------------------------
# Format description
# --------------------------------------------------------------------------

# name -> number of components ("F" means Nfreq)
KNOWN_FIELDS: Dict[str, object] = {
    "fHI": 1, "fHeI": 1, "fHeII": 1, "temp": 1, "rhog": 1,
    "ngam": "F", "nu": "F", "fgam": "F", "tau": "F",
    "gamma": 3, "lambda": 3, "dx": 3,
    "u_gas": 5, "w_gas": 4,
}
FIELD_ORDER = list(KNOWN_FIELDS)
EXTENSIONS = ("", ".bin", ".dat", ".raw", ".f32")
STATS_FILE = "stats.txt"

# stats keys compared exactly (integers) unless --stats-int-tol is given
STATS_INT_KEYS = ("istep", "nstep_rt", "nray_active",
                  "ncst1", "nmrg1", "nspl1", "ndel1")
# stats keys compared with a relative tolerance
STATS_FLOAT_KEYS = ("redshift", "scalefactor", "ngam_rays_total")

STEP_RE = re.compile(r"^step_(\d{4,})$")


class StructureError(Exception):
    """Inputs cannot be compared (missing files, inconsistent sizes...)."""


# --------------------------------------------------------------------------
# Reading
# --------------------------------------------------------------------------

def list_steps(d: str) -> Dict[int, str]:
    if not os.path.isdir(d):
        raise StructureError(f"not a directory: {d}")
    out = {}
    for name in os.listdir(d):
        m = STEP_RE.match(name)
        if m and os.path.isdir(os.path.join(d, name)):
            out[int(m.group(1))] = os.path.join(d, name)
    return dict(sorted(out.items()))


def find_field_file(stepdir: str, name: str) -> Optional[str]:
    for ext in EXTENSIONS:
        p = os.path.join(stepdir, name + ext)
        if os.path.isfile(p):
            return p
    return None


def field_name_of(path: str) -> str:
    base = os.path.basename(path)
    for ext in EXTENSIONS[1:]:
        if base.endswith(ext):
            return base[: -len(ext)]
    return base


def infer_N(stepdir: str) -> int:
    """Infer N from a scalar field (preferred) or from any raw file."""
    for name in ("fHI", "temp", "rhog", "fHeI", "fHeII"):
        p = find_field_file(stepdir, name)
        if p:
            nbytes = os.path.getsize(p)
            n = round((nbytes / 4) ** (1.0 / 3.0))
            if n ** 3 * 4 == nbytes:
                return n
            raise StructureError(f"{p}: size {nbytes} is not 4*N^3")
    raise StructureError(f"cannot infer N in {stepdir}: no scalar field found "
                         f"(use --N)")


def infer_nfreq(stepdir: str, N: int) -> Optional[int]:
    for name in ("ngam", "nu", "fgam", "tau"):
        p = find_field_file(stepdir, name)
        if p:
            nbytes = os.path.getsize(p)
            if nbytes % (4 * N ** 3):
                raise StructureError(f"{p}: size {nbytes} not a multiple of 4*N^3")
            return nbytes // (4 * N ** 3)
    return None


def ncomp_of(name: str, nfreq: Optional[int]) -> Optional[int]:
    c = KNOWN_FIELDS.get(name)
    if c == "F":
        return nfreq
    return c  # int or None


def read_field(path: str, N: int, ncomp: Optional[int]) -> np.ndarray:
    """Return array of shape (ncomp, N, N, N) indexed [c, i, j, k] (0-based),
    i.e. Fortran (c,i,j,k) layout with c fastest in memory."""
    nbytes = os.path.getsize(path)
    cell = 4 * N ** 3
    if nbytes % cell:
        raise StructureError(f"{path}: size {nbytes} not a multiple of 4*N^3={cell}")
    nc = nbytes // cell
    if ncomp is not None and nc != ncomp:
        raise StructureError(f"{path}: expected {ncomp} components, file has {nc}")
    a = np.fromfile(path, dtype="<f4")
    return a.reshape((nc, N, N, N), order="F")


def read_stats(path: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    if not os.path.isfile(path):
        return out
    with open(path) as f:
        for line in f:
            line = line.split("#", 1)[0].strip()
            if not line or "=" not in line:
                continue
            k, v = line.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def parse_num(s: str) -> float:
    # accept Fortran D exponents
    return float(s.replace("D", "E").replace("d", "e"))


# --------------------------------------------------------------------------
# Comparison
# --------------------------------------------------------------------------

@dataclass
class FieldResult:
    step: int
    field: str
    ncomp: int
    max_abs: float = 0.0
    max_rel: float = 0.0
    mean_abs: float = 0.0
    n_diff: int = 0          # cells (values) that differ at all
    n_fail: int = 0          # values exceeding tolerance
    n_nan_mismatch: int = 0
    worst: Optional[Dict[str, object]] = None   # location of worst rel diff
    ref_max: float = 0.0
    rtol: float = 0.0
    atol: float = 0.0
    ok: bool = True


@dataclass
class StatResult:
    step: int
    key: str
    ref: Optional[str]
    new: Optional[str]
    diff: Optional[float]
    rel: Optional[float]
    ok: bool


def compare_arrays(a: np.ndarray, b: np.ndarray, step: int, name: str,
                   rtol: float, atol_frac: float) -> FieldResult:
    """a = reference, b = new; both (nc, N, N, N) float32."""
    a64 = a.astype(np.float64)
    b64 = b.astype(np.float64)
    nan_a = ~np.isfinite(a64)
    nan_b = ~np.isfinite(b64)
    finite = ~(nan_a | nan_b)
    ref_max = float(np.max(np.abs(a64[~nan_a]))) if np.any(~nan_a) else 0.0
    # Absolute floor: values smaller than atol are considered "zero-like".
    atol = atol_frac * ref_max
    if atol == 0.0:
        atol = np.finfo(np.float32).tiny

    absd = np.where(finite, np.abs(a64 - b64), 0.0)
    denom = np.maximum(np.maximum(np.abs(a64), np.abs(b64)), atol)
    rel = np.where(finite, absd / denom, 0.0)
    # NaN/Inf: mismatch unless both are identical non-finite values
    same_nonfinite = (nan_a & nan_b) & ((np.isnan(a64) & np.isnan(b64)) | (a64 == b64))
    nan_mis = (nan_a | nan_b) & ~same_nonfinite

    fail = (absd > atol + rtol * np.maximum(np.abs(a64), np.abs(b64))) | nan_mis
    differ = (a != b) & ~same_nonfinite

    r = FieldResult(step=step, field=name, ncomp=a.shape[0], rtol=rtol, atol=atol,
                    ref_max=ref_max)
    r.max_abs = float(absd.max()) if absd.size else 0.0
    r.max_rel = float(rel.max()) if rel.size else 0.0
    r.mean_abs = float(absd[finite].mean()) if np.any(finite) else 0.0
    r.n_diff = int(differ.sum())
    r.n_fail = int(fail.sum())
    r.n_nan_mismatch = int(nan_mis.sum())
    r.ok = r.n_fail == 0
    if r.n_diff:
        # worst = NaN mismatch if any, else max relative diff (ties -> max abs)
        if r.n_nan_mismatch:
            flat = int(np.flatnonzero(nan_mis.ravel(order="F"))[0])
        else:
            score = rel + 1e-300 * absd
            flat = int(np.argmax(score.ravel(order="F")))
        c, i, j, k = np.unravel_index(flat, a.shape, order="F")
        r.worst = {
            "comp": int(c) + 1, "i": int(i) + 1, "j": int(j) + 1, "k": int(k) + 1,
            "ref": float(a64[c, i, j, k]), "new": float(b64[c, i, j, k]),
            "abs": float(absd[c, i, j, k]), "rel": float(rel[c, i, j, k]),
        }
    return r


def compare_stats(sa: Dict[str, str], sb: Dict[str, str], step: int,
                  int_tol: float, float_rtol: float) -> List[StatResult]:
    out: List[StatResult] = []
    keys = list(STATS_INT_KEYS) + list(STATS_FLOAT_KEYS)
    keys += sorted((set(sa) | set(sb)) - set(keys))
    for k in keys:
        va, vb = sa.get(k), sb.get(k)
        if va is None and vb is None:
            continue
        if va is None or vb is None:
            out.append(StatResult(step, k, va, vb, None, None, False))
            continue
        try:
            fa, fb = parse_num(va), parse_num(vb)
        except ValueError:
            out.append(StatResult(step, k, va, vb, None, None, va == vb))
            continue
        d = fb - fa
        rel = abs(d) / max(abs(fa), abs(fb), 1e-300) if d != 0 else 0.0
        if k in STATS_INT_KEYS:
            ok = abs(d) <= int_tol * max(abs(fa), 1.0) if int_tol > 0 else d == 0
        else:
            ok = rel <= float_rtol or (math.isnan(fa) and math.isnan(fb))
        out.append(StatResult(step, k, va, vb, d, rel, ok))
    return out


# --------------------------------------------------------------------------
# Global histories
# --------------------------------------------------------------------------

HISTORY_KEYS = ("xHII_V", "xHII_M", "xHeII_V", "xHeII_M", "xHeIII_V",
                "xHeIII_M", "T_V", "T_M", "ngam_grid", "ngam_rays_total",
                "nray_active")


def history_row(stepdir: str, N: int, nfreq: Optional[int]) -> Dict[str, float]:
    h: Dict[str, float] = {}

    def get(name):
        p = find_field_file(stepdir, name)
        return None if p is None else read_field(p, N, ncomp_of(name, nfreq)).astype(np.float64)

    fHI, fHeI, fHeII = get("fHI"), get("fHeI"), get("fHeII")
    temp, rhog, ngam = get("temp"), get("rhog"), get("ngam")
    w = None
    if rhog is not None:
        w = rhog[0] / max(rhog[0].sum(), 1e-300)

    def both(tag, arr):
        if arr is None:
            return
        h[tag + "_V"] = float(arr.mean())
        if w is not None:
            h[tag + "_M"] = float((arr * w).sum())

    both("xHII", None if fHI is None else 1.0 - fHI[0])
    both("xHeII", None if fHeII is None else fHeII[0])
    both("xHeIII", None if (fHeI is None or fHeII is None) else 1.0 - fHeI[0] - fHeII[0])
    both("T", None if temp is None else temp[0])
    if ngam is not None:
        h["ngam_grid"] = float(ngam.sum())
    st = read_stats(os.path.join(stepdir, STATS_FILE))
    for k in ("ngam_rays_total", "nray_active", "redshift"):
        if k in st:
            try:
                h[k] = parse_num(st[k])
            except ValueError:
                pass
    return h


# --------------------------------------------------------------------------
# Plotting (optional)
# --------------------------------------------------------------------------

def make_plots(outdir: str, ref_steps, new_steps, slice_steps, steps, N, nfreq,
               hist_ref, hist_new):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("WARNING: matplotlib not available; skipping plots", file=sys.stderr)
        return
    os.makedirs(outdir, exist_ok=True)
    kmid = N // 2
    for s in slice_steps:
        fig, axes = plt.subplots(2, 3, figsize=(12, 7.5), constrained_layout=True)
        for row, (name, logscale) in enumerate((("fHI", True), ("temp", True))):
            pa = find_field_file(ref_steps[s], name)
            pb = find_field_file(new_steps[s], name)
            if not (pa and pb):
                continue
            a = read_field(pa, N, 1)[0, :, :, kmid].astype(np.float64)
            b = read_field(pb, N, 1)[0, :, :, kmid].astype(np.float64)
            if logscale:
                fa, fb = np.log10(np.maximum(a, 1e-30)), np.log10(np.maximum(b, 1e-30))
                lab = f"log10 {name}"
            else:
                fa, fb, lab = a, b, name
            vmin, vmax = min(fa.min(), fb.min()), max(fa.max(), fb.max())
            for col, (arr, title) in enumerate(((fa, "ref"), (fb, "new"))):
                im = axes[row, col].imshow(arr.T, origin="lower", vmin=vmin, vmax=vmax,
                                           cmap="viridis")
                axes[row, col].set_title(f"{title}: {lab}  (k={kmid + 1})")
                fig.colorbar(im, ax=axes[row, col], shrink=0.8)
            rel = (b - a) / np.maximum(np.maximum(np.abs(a), np.abs(b)), 1e-30)
            lim = max(float(np.abs(rel).max()), 1e-12)
            im = axes[row, 2].imshow(rel.T, origin="lower", cmap="RdBu_r", vmin=-lim, vmax=lim)
            axes[row, 2].set_title(f"(new-ref)/max: {name}  max={lim:.2e}")
            fig.colorbar(im, ax=axes[row, 2], shrink=0.8)
        fig.suptitle(f"step {s:04d} mid-plane slice")
        fn = os.path.join(outdir, f"slice_step_{s:04d}.png")
        fig.savefig(fn, dpi=100)
        plt.close(fig)
    # histories
    keys = [k for k in HISTORY_KEYS if any(k in hist_ref[s] for s in steps)]
    if keys:
        ncol = 3
        nrow = (len(keys) + ncol - 1) // ncol
        fig, axes = plt.subplots(nrow, ncol, figsize=(4 * ncol, 3 * nrow),
                                 constrained_layout=True, squeeze=False)
        for ax, k in zip(axes.ravel(), keys):
            xs = [s for s in steps if k in hist_ref[s] and k in hist_new[s]]
            ax.plot(xs, [hist_ref[s][k] for s in xs], "o-", label="ref")
            ax.plot(xs, [hist_new[s][k] for s in xs], "x--", label="new")
            ax.set_title(k)
            ax.set_xlabel("step")
            ax.legend(fontsize=8)
        for ax in axes.ravel()[len(keys):]:
            ax.axis("off")
        fig.savefig(os.path.join(outdir, "history.png"), dpi=100)
        plt.close(fig)
    print(f"Plots written to {outdir}")


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def parse_overrides(items: List[str]) -> Dict[str, float]:
    out = {}
    for it in items or []:
        if "=" not in it:
            raise SystemExit(f"bad --field-rtol '{it}', expected NAME=VALUE")
        k, v = it.split("=", 1)
        out[k.strip()] = float(v)
    return out


def fmt_loc(w: Optional[Dict[str, object]], ncomp: int) -> str:
    if not w:
        return "-"
    c = f"c={w['comp']}," if ncomp > 1 else ""
    return (f"({c}{w['i']},{w['j']},{w['k']}) ref={w['ref']:.6g} new={w['new']:.6g}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("ref_dir")
    ap.add_argument("new_dir")
    ap.add_argument("--N", type=int, default=None, help="grid size (default: infer)")
    ap.add_argument("--nfreq", type=int, default=None, help="Nfreq (default: infer)")
    ap.add_argument("--rtol", type=float, default=1e-5,
                    help="relative tolerance on fields (default 1e-5)")
    ap.add_argument("--atol-frac", type=float, default=1e-6,
                    help="absolute floor per field = atol_frac * max|ref field| "
                         "(default 1e-6); values below it count as zero")
    ap.add_argument("--field-rtol", action="append", default=[], metavar="NAME=RTOL",
                    help="per-field rtol override (repeatable)")
    ap.add_argument("--stats-int-tol", type=float, default=0.0,
                    help="relative tolerance on integer ray stats (default 0 = exact)")
    ap.add_argument("--stats-rtol", type=float, default=None,
                    help="relative tolerance on float stats (default: --rtol)")
    ap.add_argument("--fields", default=None,
                    help="comma-separated subset of fields to compare")
    ap.add_argument("--ignore", default="",
                    help="comma-separated fields to skip")
    ap.add_argument("--steps", default=None,
                    help="comma-separated step numbers (default: all common)")
    ap.add_argument("--plots", default=None, metavar="OUTDIR",
                    help="write PNG slices/histories (needs matplotlib)")
    ap.add_argument("--plot-steps", default="last",
                    help="'last', 'all', or comma-separated steps for slice plots")
    ap.add_argument("--json", default=None, metavar="FILE",
                    help="write the full report as JSON")
    ap.add_argument("-q", "--quiet", action="store_true",
                    help="only print failing fields/stats and the summary")
    args = ap.parse_args(argv)

    try:
        return run(args)
    except StructureError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2


def run(args) -> int:
    rtol_over = parse_overrides(args.field_rtol)
    stats_rtol = args.rtol if args.stats_rtol is None else args.stats_rtol
    ignore = {s for s in args.ignore.split(",") if s}

    ref_steps = list_steps(args.ref_dir)
    new_steps = list_steps(args.new_dir)
    if not ref_steps:
        raise StructureError(f"no step_NNNN directories in {args.ref_dir}")
    if not new_steps:
        raise StructureError(f"no step_NNNN directories in {args.new_dir}")
    structural = False
    only_ref = sorted(set(ref_steps) - set(new_steps))
    only_new = sorted(set(new_steps) - set(ref_steps))
    if only_ref:
        print(f"WARNING: steps only in REF: {only_ref}")
        structural = True
    if only_new:
        print(f"WARNING: steps only in NEW: {only_new}")
        structural = True
    steps = sorted(set(ref_steps) & set(new_steps))
    if args.steps:
        want = {int(s) for s in args.steps.split(",")}
        steps = [s for s in steps if s in want]
    if not steps:
        raise StructureError("no common steps to compare")

    s0 = steps[0]
    N = args.N or infer_N(ref_steps[s0])
    N_new = args.N or infer_N(new_steps[s0])
    if N != N_new:
        raise StructureError(f"grid size differs: ref N={N}, new N={N_new}")
    nfreq = args.nfreq or infer_nfreq(ref_steps[s0], N)
    print(f"REF  {args.ref_dir}\nNEW  {args.new_dir}")
    print(f"N = {N}, Nfreq = {nfreq}, steps = {steps[0]}..{steps[-1]} ({len(steps)}), "
          f"rtol = {args.rtol:g}, atol_frac = {args.atol_frac:g}, "
          f"int-stat tol = {args.stats_int_tol:g}")

    field_results: List[FieldResult] = []
    stat_results: List[StatResult] = []
    hist_ref: Dict[int, Dict[str, float]] = {}
    hist_new: Dict[int, Dict[str, float]] = {}
    first_fail: Optional[int] = None
    missing: List[str] = []

    for s in steps:
        da, db = ref_steps[s], new_steps[s]
        # field list: known fields in order, then extras found in either dir
        names = [n for n in FIELD_ORDER
                 if find_field_file(da, n) or find_field_file(db, n)]
        extra = set()
        for d in (da, db):
            for fn in os.listdir(d):
                p = os.path.join(d, fn)
                if fn == STATS_FILE or not os.path.isfile(p):
                    continue
                nm = field_name_of(fn)
                if nm in KNOWN_FIELDS:
                    continue
                if os.path.getsize(p) % (4 * N ** 3) == 0 and os.path.getsize(p) > 0:
                    extra.add(nm)
        names += sorted(extra)
        if args.fields:
            sel = [f for f in args.fields.split(",") if f]
            names = [n for n in names if n in sel]
        names = [n for n in names if n not in ignore]

        if not args.quiet:
            print(f"\n=== step {s:04d} ===")
            print(f"  {'field':<8} {'max_abs':>11} {'max_rel':>11} {'mean_abs':>11} "
                  f"{'n_diff':>9} {'n_fail':>8}  worst (1-based c,i,j,k)")
        step_ok = True
        for n in names:
            pa, pb = find_field_file(da, n), find_field_file(db, n)
            if not (pa and pb):
                msg = f"step {s:04d}: field {n} missing in {'REF' if not pa else 'NEW'}"
                print(f"  {n:<8} MISSING in {'REF' if not pa else 'NEW'}")
                missing.append(msg)
                structural = True
                continue
            nc = ncomp_of(n, nfreq)
            a = read_field(pa, N, nc)
            b = read_field(pb, N, nc)
            if a.shape != b.shape:
                raise StructureError(f"step {s}: {n} shape {a.shape} vs {b.shape}")
            r = compare_arrays(a, b, s, n, rtol_over.get(n, args.rtol), args.atol_frac)
            field_results.append(r)
            if not r.ok:
                step_ok = False
            if not args.quiet or not r.ok:
                flag = "" if r.ok else "  <-- FAIL"
                prefix = "" if not args.quiet else f"step {s:04d} "
                print(f"  {prefix}{n:<8} {r.max_abs:11.4e} {r.max_rel:11.4e} "
                      f"{r.mean_abs:11.4e} {r.n_diff:9d} {r.n_fail:8d}  "
                      f"{fmt_loc(r.worst, r.ncomp)}{flag}")

        sa = read_stats(os.path.join(da, STATS_FILE))
        sb = read_stats(os.path.join(db, STATS_FILE))
        if not sa or not sb:
            print(f"  stats.txt missing in {'REF' if not sa else 'NEW'}")
            missing.append(f"step {s:04d}: stats.txt missing")
            structural = True
        srs = compare_stats(sa, sb, s, args.stats_int_tol, stats_rtol)
        stat_results.extend(srs)
        bad = [x for x in srs if not x.ok]
        if bad:
            step_ok = False
        if not args.quiet:
            print("  stats:")
            for x in srs:
                flag = "" if x.ok else "  <-- FAIL"
                d = "" if x.diff is None else f" diff={x.diff:.6g} rel={x.rel:.3e}"
                print(f"    {x.key:<16} ref={x.ref!s:<24} new={x.new!s:<24}{d}{flag}")
        else:
            for x in bad:
                print(f"  step {s:04d} stat {x.key}: ref={x.ref} new={x.new}  <-- FAIL")

        hist_ref[s] = history_row(da, N, nfreq)
        hist_new[s] = history_row(db, N, nfreq)
        if not step_ok and first_fail is None:
            first_fail = s

    # ---- histories
    print("\n=== global histories (ref | new | rel diff) ===")
    keys = [k for k in HISTORY_KEYS if any(k in hist_ref[s] or k in hist_new[s] for s in steps)]
    for k in keys:
        print(f"  {k}")
        for s in steps:
            a, b = hist_ref[s].get(k), hist_new[s].get(k)
            if a is None or b is None:
                print(f"    {s:04d}  ref={a}  new={b}")
                continue
            rel = abs(b - a) / max(abs(a), abs(b), 1e-300) if a != b else 0.0
            print(f"    {s:04d}  {a:16.9e} {b:16.9e} {rel:10.3e}")

    # ---- summary
    nfail_fields = sum(1 for r in field_results if not r.ok)
    nfail_stats = sum(1 for r in stat_results if not r.ok)
    print("\n=== summary ===")
    print(f"  fields compared: {len(field_results)}, failing: {nfail_fields}")
    print(f"  stats compared:  {len(stat_results)}, failing: {nfail_stats}")
    if field_results:
        worst = max(field_results, key=lambda r: r.max_rel)
        print(f"  worst field: step {worst.step:04d} {worst.field} max_rel={worst.max_rel:.3e} "
              f"at {fmt_loc(worst.worst, worst.ncomp)}")
    if first_fail is not None:
        print(f"  FIRST DIVERGENCE beyond tolerance at step {first_fail:04d}")
    if missing:
        print(f"  structural problems: {len(missing)}")
    passed = nfail_fields == 0 and nfail_stats == 0 and not structural
    print("  RESULT:", "PASS" if passed else "FAIL")

    if args.json:
        rep = {
            "ref_dir": args.ref_dir, "new_dir": args.new_dir, "N": N, "nfreq": nfreq,
            "steps": steps, "rtol": args.rtol, "atol_frac": args.atol_frac,
            "fields": [asdict(r) for r in field_results],
            "stats": [asdict(r) for r in stat_results],
            "history_ref": hist_ref, "history_new": hist_new,
            "first_divergence": first_fail, "missing": missing,
            "only_ref": only_ref, "only_new": only_new, "pass": passed,
        }
        with open(args.json, "w") as f:
            json.dump(rep, f, indent=1, default=str)

    if args.plots:
        if args.plot_steps == "all":
            psteps = steps
        elif args.plot_steps == "last":
            psteps = [steps[-1]]
        else:
            psteps = [int(x) for x in args.plot_steps.split(",") if int(x) in steps]
        make_plots(args.plots, ref_steps, new_steps, psteps, steps, N, nfreq,
                   hist_ref, hist_new)

    if structural and nfail_fields == 0 and nfail_stats == 0:
        return 2
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
