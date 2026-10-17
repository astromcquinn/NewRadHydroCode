#!/usr/bin/env python3
"""Tabulate the hydro suite's errors with and without refinement.

usage: hydro_suite_amr.py OUT [--tier quick] [--tests tubes,sedov] [--md FILE]

OUT holds <tier>_unsplit/results/<test>_<tier>_unsplit.json (uniform) and
<tier>_unsplit_amrL/results/<test>_<tier>_unsplit_amrL.json (L refinement
levels on the same base grids; tools/hyak/hydro_suite_amr.slurm). The
analysis of a refined run reads its level-0 fields, i.e. the fine solution
averaged down, so a refined run on N^3 should sit between the uniform N^3
and (2^L N)^3 runs. Prints (and with --md writes) one table per test: the
L1 errors of every case for uniform / amr1 / amr2.
"""
import argparse
import json
import os


def load(out, tier, test, sfx):
    p = os.path.join(out, f"{tier}_unsplit{sfx}", "results", f"{test}_{tier}_unsplit{sfx}.json")
    return json.load(open(p)) if os.path.exists(p) else None


def rows_tube(res):
    for name, r in sorted(res["cases"].items()):
        yield name, {q: r["l1"][q] for q in ("rho", "u", "p")}, r.get("wall", float("nan"))


def rows_sedov(res):
    for name, r in sorted(res["cases"].items()):
        o = r["outputs"][-1]
        yield name, {"rho": o["l1"]["rho"], "p": o["l1"]["p"], "v": o["l1"]["v"],
                     "R_err_rel": o["R_err_rel"], "aniso": o["aniso_rms"]}, r.get("wall", float("nan"))


ROWS = {"tubes": rows_tube, "sedov": rows_sedov}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--tier", default="quick")
    ap.add_argument("--tests", default="tubes,sedov")
    ap.add_argument("--md", default=None)
    a = ap.parse_args()
    lines = []
    for test in a.tests.split(","):
        variants = [("uniform", ""), ("amr1", "_amr1"), ("amr2", "_amr2")]
        data = {v: load(a.out, a.tier, test, s) for v, s in variants}
        data = {v: d for v, d in data.items() if d}
        if not data:
            continue
        rows = {v: dict((n, (m, w)) for n, m, w in ROWS.get(test, rows_tube)(d)) for v, d in data.items()}
        names = sorted(set().union(*(set(r) for r in rows.values())))
        metrics = list(next(iter(rows.values()))[names[0]][0].keys())
        head = "| case | " + " | ".join(f"{m} {v}" for m in metrics for v in rows) + " | wall " + "/".join(rows) + " [s] |"
        lines += [f"\n**{test}** ({a.tier}, unsplit; refine_mass of the refined runs in the JSON)\n", head,
                  "|" + "---|" * (head.count("|") - 1)]
        for n in names:
            cells = []
            for m in metrics:
                for v in rows:
                    cells.append(f"{rows[v][n][0][m]:.4f}" if n in rows[v] else "")
            walls = "/".join(f"{rows[v][n][1]:.0f}" if n in rows[v] else "" for v in rows)
            lines.append(f"| {n} | " + " | ".join(cells) + f" | {walls} |")
    txt = "\n".join(lines)
    print(txt)
    if a.md:
        open(a.md, "w").write(txt + "\n")


if __name__ == "__main__":
    main()
