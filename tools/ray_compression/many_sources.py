#!/usr/bin/env python3
"""Source-count statistics of ray snapshots (docs/ray_compression.md §5.3).

usage: many_sources.py DUMP_DIR [--codec build/ray_codec] [--json out.json]

For each snapshot in DUMP_DIR: the number of sources, the never-merged
fraction, how many distinct sources the rays of a cell come from, how often
consecutive rays in a cell share a source (what the coder's "same source"
flag sees), the entropy of the source id per ray, and the bytes per ray of
the prototype coder with its raw source id (6-bit total, 10-bit offsets).
"""
import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import snap  # noqa: E402


def H(x):
    _, c = np.unique(x, return_counts=True)
    p = c / c.sum()
    return float(-(p * np.log2(p)).sum())


def analyse(stem, codec, prev):
    a, meta = snap.load(stem)
    d = snap.derive(a, meta)
    n = len(a)
    nm = d["never_merged"]
    src = a["src"].astype(np.int64)
    cell = a["cell"].astype(np.int64)
    out = {"stem": str(stem), "substep": meta["substep"], "redshift": meta["redshift"],
           "xHII_vol": meta["xHII_vol"], "nsources": meta["nsources"], "nrays": n,
           "never_merged": float(nm.mean())}
    known = src >= 0
    out["frac_known_source"] = float(known.mean())
    if known.any():
        # distinct sources per cell, over cells with known-source rays
        key = cell[known] * (1 << 32) + src[known]
        uk = np.unique(key)
        cells_k, per_cell = np.unique(uk >> 32, return_counts=True)
        out["sources_per_cell_mean"] = float(per_cell.mean())
        out["sources_per_cell_p90"] = float(np.percentile(per_cell, 90))
        out["rays_per_cell_mean"] = float(np.bincount(cell[known]).astype(float)[cells_k].mean())
        # consecutive rays in a cell with the same source (snapshot order = chain order)
        same = (src[1:] == src[:-1]) & (cell[1:] == cell[:-1]) & known[1:] & known[:-1]
        out["frac_same_source_as_prev"] = float(same.sum() / max(known.sum() - 1, 1))
        out["H_src_bits"] = H(src[known])
        out["log2_nsources"] = float(np.log2(max(meta["nsources"], 1)))
    if codec:
        r = json.loads(subprocess.run([codec, prev, str(stem), "--reps", "1", "--mbits", "6",
                                       "--offbits", "10"], check=True, capture_output=True,
                                      text=True).stdout)
        out["bytes_entropy"] = r["bytes_entropy"]
        out["bits_src"] = r["bits_by_group"]["src"]
        out["bits_ipix"] = r["bits_by_group"]["ipix"]
        out["class_frac"] = r["class_frac"]
        out["decode_ns"] = r["decode_ns_per_ray"]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dump_dir")
    ap.add_argument("--codec", default=str(Path(__file__).resolve().parents[2] / "build" / "ray_codec"))
    ap.add_argument("--json")
    args = ap.parse_args()
    stems = [str(p)[:-4] for p in sorted(Path(args.dump_dir).glob("rays_s*.npy"))]
    res = []
    print("%8s %6s %6s %8s %8s %6s %8s %8s %7s %7s %7s" % (
        "substep", "z", "xHII", "nsrc", "rays", "never", "src/cell", "same", "H(src)", "B/ray", "bits_src"))
    for i, st in enumerate(stems):
        r = analyse(st, args.codec, stems[i - 1] if i else st)
        res.append(r)
        print("%8d %6.2f %6.3f %8d %8d %6.3f %8.2f %8.3f %7.2f %7.2f %7.2f" % (
            r["substep"], r["redshift"], r["xHII_vol"], r["nsources"], r["nrays"], r["never_merged"],
            r.get("sources_per_cell_mean", float("nan")), r.get("frac_same_source_as_prev", float("nan")),
            r.get("H_src_bits", float("nan")), r.get("bytes_entropy", float("nan")),
            r.get("bits_src", float("nan"))), flush=True)
    if args.json:
        json.dump(res, open(args.json, "w"), indent=1)


if __name__ == "__main__":
    main()
