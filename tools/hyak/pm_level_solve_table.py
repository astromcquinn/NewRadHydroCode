#!/usr/bin/env python3
"""The table of tools/hyak/pm_level_solve_profile.slurm (and of a local run of
the same configurations): per configuration the mean wall-clock per step of
the PM, the hydro and the whole step from <run>/timing.txt, and the PM
backend's cumulative phase timers of the last step (the "step NNNN timers"
line of <run>.log) divided by the step count.

usage: pm_level_solve_table.py OUT [config ...]
"""
import os
import re
import sys


def timing(path):
    """timing.txt: the columns pm, hydro, step (seconds per step), mean over the steps"""
    cols, rows = None, []
    for line in open(path):
        if line.startswith("#"):
            if "step" in line and cols is None and ("pm" in line):
                cols = line[1:].split()
            continue
        v = line.split()
        if cols and len(v) == len(cols):
            rows.append(dict(zip(cols, map(float, v))))
    if not rows:
        return None
    n = len(rows)
    return {k: sum(r[k] for r in rows) / n for k in rows[0]}, n


def main():
    out = sys.argv[1]
    cfgs = sys.argv[2:] or ["level_solve", "gas_levels", "uniform"]
    print(f"{'config':12s} {'steps':>5s} {'step [s]':>9s} {'PM [s]':>8s} {'hydro [s]':>9s}  PM backend phases per step [s]")
    for c in cfgs:
        run = os.path.join(out, c)
        log = run + ".log"
        if not os.path.isfile(log):
            print(f"{c:12s}  (no run)")
            continue
        t = timing(os.path.join(run, "timing.txt")) if os.path.isfile(os.path.join(run, "timing.txt")) else None
        lines = [l for l in open(log) if "timers:" in l and "AMReX PM timers" in l]
        last = lines[-1] if lines else ""
        m = re.search(r"steps (\d+)", last)
        nst = int(m.group(1)) if m else 0
        pm = re.search(r"AMReX PM timers \[s\]: (.*?); steps", last)
        ph = ""
        if pm and nst:
            parts = []
            for key in ("drift", "deposit", "fft", "force", "kick_dm", "kick_gas"):
                mm = re.search(key + r" ([0-9.]+)", pm.group(1))
                if mm:
                    parts.append(f"{key} {float(mm.group(1)) / nst:.3f}")
            ph = " ".join(parts)
            ls = re.search(r"level solve: (.*?); levels (\d+)(.*)", last)
            if ls:
                lp = []
                for key in ("sync", "deposit", "solve", "force", "kick"):
                    mm = re.search(key + r" ([0-9.]+)", ls.group(1))
                    if mm:
                        lp.append(f"{key} {float(mm.group(1)) / nst:.3f}")
                ph += " | levels: " + " ".join(lp) + " " + ls.group(3).strip()
            gl = re.search(r"gas levels: (.*?)$", last.strip())
            if gl:
                lp = []
                for key in ("deposit", "kick"):
                    mm = re.search(key + r" ([0-9.]+)", gl.group(1))
                    if mm:
                        lp.append(f"{key} {float(mm.group(1)) / nst:.3f}")
                ph += " | gas levels: " + " ".join(lp)
        if t:
            tt, n = t
            print(f"{c:12s} {n:5d} {tt.get('step', float('nan')):9.2f} {tt.get('pm', float('nan')):8.3f} "
                  f"{tt.get('hydro', float('nan')):9.2f}  {ph}")
        else:
            print(f"{c:12s} {nst:5d} {'':9s} {'':8s} {'':9s}  {ph}")


if __name__ == "__main__":
    main()
