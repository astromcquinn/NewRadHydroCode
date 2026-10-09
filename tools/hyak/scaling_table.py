#!/usr/bin/env python3
"""Scaling table from tools/production/speedup.sh outputs.

    scaling_table.py OUTROOT T1 T2 ...

Reads OUTROOT/t<T>/timing.txt (the step line and the "# RT timers" line)
and prints Markdown tables of the wall time per component and the speedup
relative to the first thread count, and the parallel efficiency.
"""
import os
import re
import sys

STEP_COLS = ["hydro", "pm", "halos", "rt", "summaries", "io", "total"]
RT_KEYS = ["prep", "calc_rt_1", "calc_ray_1", "ray_tracing_1", "radiative_physics",
           "ray_tracing_3", "sync", "post_rt"]
LABEL = {"halos": "FoF + sources", "rt": "RT total", "total": "step total"}


def read(root, t):
    path = os.path.join(root, f"t{t}", "timing.txt")
    out = {}
    try:
        lines = open(path).read().splitlines()
    except OSError:
        return None
    steps = [l.split() for l in lines if l.strip() and not l.startswith("#")]
    if steps:
        # several steps (STEPS > 1): sum them
        for i, k in enumerate(STEP_COLS):
            out[k] = sum(float(r[5 + i]) for r in steps)
        out["nray"] = int(steps[-1][4])
        out["nstep_rt"] = sum(int(r[2]) for r in steps)
        out["z"] = float(steps[-1][1])
    for l in lines:
        if "RT timers" in l:
            m = re.search(r"substeps (\d+)", l)
            out["substeps"] = int(m.group(1))
            for k in RT_KEYS:
                m = re.search(rf"{k} ([0-9.eE+-]+)", l)
                if m:
                    out[k] = float(m.group(1))
        if l.startswith("# total"):
            m = re.search(r"peak_rss_MB ([0-9.]+)", l)
            if m:
                out["rss"] = float(m.group(1))
    return out


def main():
    root = sys.argv[1]
    ts = [int(t) for t in sys.argv[2:]]
    data = {t: read(root, t) for t in ts}
    ts = [t for t in ts if data[t]]
    if not ts:
        print("no timings found")
        return
    b = data[ts[0]]
    rss = ", ".join("%.1f" % (data[t].get("rss", 0) / 1024) for t in ts)
    print(f"step ending at z = {b.get('z', float('nan')):.4f}, {b.get('nstep_rt', 0)} RT substeps, "
          f"{b.get('nray', 0) / 1e6:.2f}M live rays at the end; peak RSS {rss} GB\n")
    rows = [k for k in RT_KEYS + STEP_COLS if k in b]
    print("| component | " + " | ".join(f"{t} thr [s]" for t in ts) + " | "
          + " | ".join(f"speedup {t}" for t in ts[1:]) + " |")
    print("|---" * (1 + 2 * len(ts) - 1) + "|")
    for k in rows:
        v = [data[t].get(k, float("nan")) for t in ts]
        sp = [v[0] / x if x > 0 else float("nan") for x in v[1:]]
        print(f"| {LABEL.get(k, k)} | " + " | ".join(f"{x:.2f}" for x in v) + " | "
              + " | ".join(f"{x:.1f}" for x in sp) + " |")
    print()
    print("| threads | step [s] | speedup | efficiency | RT speedup |")
    print("|---|---|---|---|---|")
    for t in ts:
        s = b["total"] / data[t]["total"]
        r = b["rt"] / data[t]["rt"] if data[t].get("rt") else float("nan")
        print(f"| {t} | {data[t]['total']:.1f} | {s:.1f} | {s * ts[0] / t:.2f} | {r:.1f} |")


if __name__ == "__main__":
    main()
