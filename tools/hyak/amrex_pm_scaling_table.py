#!/usr/bin/env python3
"""Tables of tools/hyak/amrex_pm_scaling.slurm (docs/stage_e_pm.md section 4.1):
per run and per hydro step the PM time of timing.txt and the PM backend's
phase timers (the differences of the cumulative "step NNNN timers:" lines
of radhydro --set step_timers=1), the FoF and P(k) time on rank 0, and the
peak RSS of rank 0 and over the ranks.

usage: amrex_pm_scaling_table.py OUT [run ...]   (default: every run in OUT)
"""
import glob
import os
import re
import sys

PHASES = ("drift", "deposit", "SumBoundary", "fft", "force", "kick_dm", "kick_gas", "gather")


def timers(log):
    """{istep: {phase: cumulative seconds, 'halos':, 'pk':}} from the step_timers lines."""
    out = {}
    for line in open(log):
        m = re.match(r"\s*step\s+(\d+) timers:(.*)", line)
        if not m:
            continue
        d = {}
        for ph in PHASES + ("prep",):
            mm = re.search(r"\b%s ([0-9.]+)" % ph, m.group(2))
            if mm:
                d[ph] = float(mm.group(1))
        d["halos"] = float(re.search(r"halos ([0-9.]+)", m.group(2)).group(1))
        d["pk"] = float(re.search(r"P\(k\) ([0-9.]+)", m.group(2)).group(1))
        out[int(m.group(1))] = d
    return out


def timing(run):
    """{istep: (t_pm, t_halos, t_total)} from timing.txt."""
    out = {}
    for line in open(os.path.join(run, "timing.txt")):
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        out[int(f[0])] = (float(f[6]), float(f[7]), float(f[11]))
    return out


def rss(log):
    r0 = lo = hi = None
    for line in open(log):
        m = re.search(r"peak RSS ([0-9.]+) MB", line)
        if m and "radhydro:" in line:
            r0 = float(m.group(1)) / 1024
        m = re.search(r"peak RSS over (\d+) ranks: min ([0-9.]+) MB, max ([0-9.]+) MB", line)
        if m:
            lo, hi = float(m.group(2)) / 1024, float(m.group(3)) / 1024
    return r0, lo, hi


def main():
    out = sys.argv[1]
    runs = sys.argv[2:] or sorted(d for d in glob.glob(os.path.join(out, "*"))
                                  if os.path.isdir(d) and os.path.exists(d + ".log"))
    print("| run | step | PM [s] | " + " | ".join(PHASES) + " | FoF [s] | P(k) [s] | step [s] |")
    print("|" + "---|" * (len(PHASES) + 6))
    mem = []
    for run in runs:
        name = os.path.basename(run)
        log = run + ".log"
        if not os.path.exists(os.path.join(run, "timing.txt")):
            print(f"| {name} | failed | |")
            continue
        tm = timing(run)
        tr = timers(log)
        prev = None
        for st in sorted(tm):
            t_pm, t_halo, t_tot = tm[st]
            cur = tr.get(st)
            if cur is not None:
                dif = {k: cur.get(k, 0) - (prev or {}).get(k, 0) for k in cur}
                ph = " | ".join(f"{dif.get(p, 0):.2f}" if p in cur else "" for p in PHASES)
                pk = dif["pk"]
                prev = cur
            else:
                ph = " | ".join("" for _ in PHASES)
                pk = float("nan")
            print(f"| {name} | {st} | {t_pm:.2f} | {ph} | {t_halo:.2f} | {pk:.2f} | {t_tot:.2f} |")
        mem.append((name, rss(log)))
    print()
    print("| run | peak RSS rank 0 [GB] | min over ranks [GB] | max over ranks [GB] |")
    print("|---|---|---|---|")
    for name, (r0, lo, hi) in mem:
        f = lambda x: "" if x is None else f"{x:.1f}"
        print(f"| {name} | {f(r0)} | {f(lo)} | {f(hi)} |")


if __name__ == "__main__":
    main()
