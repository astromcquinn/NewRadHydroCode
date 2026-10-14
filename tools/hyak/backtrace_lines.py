#!/usr/bin/env python3
"""Resolve AMReX Backtrace.<rank>.<thread> files to source lines.

A hung or cancelled backend run (SIGTERM: scancel, the time limit) writes
one file per rank into the run directory (amrex.handle_sigterm = 1 in
src/amrex/rt_level.cpp) with raw addresses, since addr2line in the signal
handler is too slow for the first rank's MPI_Abort to wait for. This
script runs addr2line (-g1 line tables, tools/hyak/build_amrex.sh) on the
executable and prints each distinct stack once with the ranks that share
it, deepest frame first:

    python3 tools/hyak/backtrace_lines.py build-amrex/radhydro $OUT/amrex_n4/Backtrace.*
"""
import collections
import os
import re
import subprocess
import sys


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    exe = sys.argv[1]
    base = os.path.basename(exe)
    pat = re.compile(r"^\s*(\d+):\s+(\S+)\(([^)]*)\)\s+\[(0x[0-9a-f]+)\]")
    stacks = collections.OrderedDict()
    for f in sys.argv[2:]:
        m = re.search(r"Backtrace\.(?:prog\d+\.)?(\d+)(?:\.(\d+))?$", f)
        rank = f"{m.group(1)}.{m.group(2) or 0}" if m else f
        frames = []
        for line in open(f, errors="replace"):
            m = pat.match(line)
            if not m:
                continue
            path, sym, addr = m.group(2), m.group(3), m.group(4)
            if os.path.basename(path) == base:
                frames.append((addr, None))
            else:
                frames.append((None, f"{os.path.basename(path)}({sym})"))
        addrs = [a for a, _ in frames if a]
        res = {}
        if addrs:
            out = subprocess.run(["addr2line", "-Cpfie", exe] + addrs, stdout=subprocess.PIPE, universal_newlines=True).stdout
            lines = [l for l in out.splitlines() if l.strip()]
            # -i inlines add lines starting with " (inlined by)"; keep the outermost per address
            res_lines = [l for l in lines if not l.startswith(" (inlined")]
            for a, l in zip(addrs, res_lines):
                res[a] = re.sub(r" at .*?/NewRadHydroCode/", " at ", l)
        text = "\n".join(res.get(a, s) if a else s for a, s in frames if (a and a in res) or s)
        if not text:
            text = "(no frames: the file is empty, the rank was killed before it wrote)"
        stacks.setdefault(text, []).append(rank)
    for text, ranks in sorted(stacks.items(), key=lambda kv: -len(kv[1])):
        print(f"=== {len(ranks)} rank(s): {' '.join(ranks)}")
        print(text)
        print()


if __name__ == "__main__":
    main()
