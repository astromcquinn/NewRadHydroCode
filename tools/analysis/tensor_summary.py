#!/usr/bin/env python3
"""Collect tensor_compare.py logs of several snapshots into one table and a
figure (docs/otvet_closure.md §3): errors of each closure against the rays'
Eddington tensor as a function of redshift.

usage: tensor_summary.py LOG_OR_JSON... [--closures otvet split ...] [--png FILE]
                         [--txt FILE]

Several analyses of the same snapshot (same z) are merged into one row.

The screened kernels are logged as atten<lambda> with lambda in cells; the
--closures names 'atten' and 'split_atten' pick the one at lambda_eff.
"""
import argparse
import re

import numpy as np


def untruncate(names):
    """The logs truncate the cell- and block-table names to 10 characters, so
    every split_<x> reads 'split_atte' or 'split_cone'; each one follows its
    <x> column."""
    out = list(names)
    for i, n in enumerate(out):
        if n in ("split_atte", "split_cone") and i > 0:
            out[i] = "split_" + out[i - 1]
    return out


def parse(path):
    """A metrics_z*.json (tensor_compare.py, untruncated names) or a log.txt.
    """
    if path.endswith(".json"):
        import json
        with open(path) as f:
            return json.load(f)
    return parse_log(path)


def parse_log(path):
    """One log -> dict(z, xv, lam, cell{regime: {closure: err}},
    A{regime: {closure: A}}, Atrue{regime}, block{regime: {closure: err}},
    blockA{closure}, force{closure: [err/true per k]}, r{closure: [...]}, k)."""
    out = {"cell": {}, "A": {}, "Atrue": {}, "noise": {}, "block": {}, "blockA": {}, "force": {}, "r": {}}
    sec = None
    names = None
    last = None
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("z = "):
                m = re.search(r"z = ([\d.]+).*x_V = ([\d.]+)\s+lambda_eff = ([\d.]+)", line)
                out["z"], out["xv"], out["lam"] = map(float, m.groups())
            elif line.startswith("cell level"):
                sec = "cell"
            elif line.startswith("regime") and sec == "cell":
                names = untruncate(line.split()[4:])
            elif line.startswith("anisotropy A of each"):
                sec = "A"
            elif line.startswith("radiation force"):
                sec = "force"
            elif line.startswith("anisotropy A vs distance"):
                sec = None
            elif line.startswith("block 8"):
                sec = "block"
            elif line.startswith("regime") and sec == "block":
                bnames = untruncate(line.split()[1:])
            elif line.startswith("block anisotropy"):
                for m in re.finditer(r"(\S+) ([\d.]+)", line.split(":", 1)[1].replace("|", " ")):
                    out["blockA"][m.group(1)] = float(m.group(2))
            elif sec == "cell" and names and re.match(r"^\S.*\s[\d.]+\s+[\d.]+\s+[\d.]+", line):
                reg = line[:26].strip()
                v = line[26:].split()
                out["Atrue"][reg] = float(v[1])
                out["noise"][reg] = float(v[2])
                out["cell"][reg] = dict(zip(names, map(float, v[3:])))
            elif sec == "A" and names and line[:26].strip():
                reg = line[:26].strip()
                out["A"][reg] = dict(zip(names, map(float, line[26:].split())))
            elif sec == "force":
                t = line.split()
                if not t:
                    continue
                if t[0] == "k":
                    out["k"] = list(map(float, t[2:]))
                elif t[0] == "r" and last:
                    out["r"][last] = list(map(float, t[1:]))
                else:
                    last = t[0]
                    out["force"][last] = list(map(float, t[1:]))
            elif sec == "block" and line[:26].strip() and not line.startswith("regime"):
                reg = line[:26].strip()
                out["block"][reg] = dict(zip(bnames, map(float, line[26:].split())))
    return out


def merge(parsed):
    """Analyses of the same snapshot (several runs of tensor_compare.py with
    different closures) become one entry; the first one's value wins where
    they overlap. Sorted by decreasing z."""
    by = {}
    for d in parsed:
        k = round(d["z"], 3)
        if k not in by:
            by[k] = d
            continue
        m = by[k]
        for t in ("cell", "A", "block"):
            for reg, row in d[t].items():
                for c, v in row.items():
                    m[t].setdefault(reg, {}).setdefault(c, v)
        for t in ("blockA", "force", "r"):
            for c, v in d[t].items():
                m[t].setdefault(c, v)
        for t in ("Atrue", "noise"):
            for reg, v in d[t].items():
                m[t].setdefault(reg, v)
    return [by[k] for k in sorted(by, reverse=True)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("logs", nargs="+")
    ap.add_argument("--closures", nargs="*", default=["otvet", "split", "atten", "cone2c2t1", "m1"])
    ap.add_argument("--png", default=None)
    ap.add_argument("--txt", default=None)
    a = ap.parse_args()
    runs = merge(parse(p) for p in a.logs)
    lines = []

    def key_of(d, table, name):
        """'atten'/'split_atten' -> the screened kernel at lambda_eff; other
        names as logged (cell tables: the first 10 characters)."""
        if name in ("atten", "split_atten"):
            best = None
            for k in table:
                m = re.fullmatch(name + r"(\d+)", k)
                if m and (best is None or abs(int(m.group(1)) - d["lam"]) < abs(int(best[1]) - d["lam"])):
                    best = (k, m.group(1))
            return best[0] if best else None
        for k in table:
            if k == name or (len(k) == 10 and name.startswith(k)):
                return k
        return None

    def col(d, table, name):
        k = key_of(d, table, name)
        return table[k] if k else np.nan

    hdr = "%6s %5s %6s | %-6s" % ("z", "x_V", "lam", "A_true") + "".join("%16s" % c[:15] for c in a.closures)
    for title, get in [
        ("cell error |h - h_true|, ionised x > 0.9 (S0-weighted)", lambda d, c: col(d, d["cell"].get("ionised x>0.9", {}), c)),
        ("cell error, ionised cells with Gamma < 0.5 <Gamma>", lambda d, c: col(d, d["cell"].get("ion G<0.5<G>", {}), c)),
        ("cell error, shadowed cells N < N_loc/2", lambda d, c: col(d, d["cell"].get("N<N_loc/2 (shadowed)", {}), c)),
        ("block 8^3 error, ionised blocks", lambda d, c: col(d, d["block"].get("ionised blocks x>0.9", {}), c)),
        ("block 8^3 anisotropy (truth in the A_true column)", lambda d, c: col(d, d["blockA"], c)),
    ]:
        lines.append("## " + title)
        lines.append(hdr)
        for d in runs:
            at = d["Atrue"].get("ionised x>0.9", np.nan) if "anisotropy" not in title else d["blockA"].get("truth", np.nan)
            lines.append("%6.2f %5.2f %6.0f | %6.3f" % (d["z"], d["xv"], d["lam"], at)
                         + "".join("%16.3f" % get(d, c) for c in a.closures))
        lines.append("")
    for ik in range(len(runs[0]["k"])):
        k = runs[0]["k"][ik]
        if k > 16:
            continue
        lines.append("## radiation force error power / true power at k = %.2f h/Mpc (r in brackets)" % k)
        lines.append("%6s %5s %6s |" % ("z", "x_V", "lam") + "".join("%15s" % c[:15] for c in a.closures))
        for d in runs:
            row = "%6.2f %5.2f %6.0f |" % (d["z"], d["xv"], d["lam"])
            for c in a.closures:
                key = key_of(d, d["force"], c)
                e = d["force"].get(key, [np.nan] * 8)[ik]
                r = d["r"].get(key, [np.nan] * 8)[ik]
                row += "%15s" % ("%.3f (%.2f)" % (e, r))
            lines.append(row)
        lines.append("")
    txt = "\n".join(lines)
    print(txt)
    if a.txt:
        open(a.txt, "w").write(txt + "\n")
    if a.png:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        z = np.array([d["z"] for d in runs])
        fig, ax = plt.subplots(2, 3, figsize=(15, 8.5))
        ax = ax.ravel()
        panels = [
            ("cell error, ionised", lambda d, c: col(d, d["cell"].get("ionised x>0.9", {}), c)),
            ("cell error, Γ < 0.5⟨Γ⟩", lambda d, c: col(d, d["cell"].get("ion G<0.5<G>", {}), c)),
            ("block 8³ error, ionised", lambda d, c: col(d, d["block"].get("ionised blocks x>0.9", {}), c)),
        ]
        for i, (t, get) in enumerate(panels):
            for c in a.closures:
                ax[i].plot(z, [get(d, c) for d in runs], "o-", label=c)
            ax[i].set_title(t)
            ax[i].set_xlabel("z")
            ax[i].invert_xaxis()
        ax[0].legend(fontsize=8)
        # block anisotropy
        ax[3].plot(z, [d["blockA"].get("truth", np.nan) for d in runs], "k-", lw=3, label="truth")
        for c in a.closures:
            ax[3].plot(z, [col(d, d["blockA"], c) for d in runs], "o-", label=c)
        ax[3].set_title("block 8³ anisotropy A")
        ax[3].invert_xaxis()
        ax[3].set_xlabel("z")
        for j, ik in enumerate([0, 3]):
            for c in a.closures:
                vals = [d["force"].get(key_of(d, d["force"], c), [np.nan] * 8)[ik] for d in runs]
                ax[4 + j].semilogy(z, vals, "o-", label=c)
            ax[4 + j].set_title("force error power / true, k = %.2f h/Mpc" % runs[0]["k"][ik])
            ax[4 + j].invert_xaxis()
            ax[4 + j].set_xlabel("z")
        # x_V on top axis of first panel
        for i in range(6):
            ax[i].grid(alpha=0.3)
        fig.suptitle("a priori closure errors against the rays' Eddington tensor (256³, 20 Mpc/h); "
                     "x_V: " + ", ".join("%.2f" % d["xv"] for d in runs))
        fig.tight_layout()
        fig.savefig(a.png, dpi=110)


if __name__ == "__main__":
    main()
