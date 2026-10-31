#!/usr/bin/env python3
"""amrex_pk = level (docs/rank0_streaming.md step 11): the power spectrum of
the data1 outputs on the PM and gas levels.

usage: pk_level_check.py --launch LAUNCHER FLAGS NPFLAG NP RADHYDRO PROBLEM WORKDIR [common args ...]

Three runs of PROBLEM with amrex_state = 1 (the common args):
  A  amrex_pk = host on 1 rank;
  B  amrex_pk = level on 1 rank;
  C  amrex_pk = level on NP ranks.
Checks: every file of B equals A's except Pk/pk.*, whose numbers agree to a
relative 1e-4 of each column's largest value (the host and level FFTs and
DM deposits differ at round-off); every file of C equals B's, Pk/ included
(the level spectrum is rank independent). Exit 0 on success.
"""
import glob
import os
import shutil
import subprocess
import sys

SKIP = ("timing.txt", "ray.log", "rt.log", "stdout.log")


def run(cmd, out):
    if os.path.isdir(out):
        shutil.rmtree(out)
    os.makedirs(out)
    with open(out + ".log", "w") as f:
        r = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
    if r.returncode != 0:
        print(open(out + ".log").read()[-2000:])
        sys.exit("failed: " + " ".join(cmd))


def files_of(d):
    return sorted(os.path.relpath(f, d) for f in glob.glob(os.path.join(d, "**", "*"), recursive=True)
                  if os.path.isfile(f) and os.path.basename(f) not in SKIP
                  and not os.path.relpath(f, d).startswith("Log/"))


def numbers(path):
    return [[float(x) for x in line.split()] for line in open(path) if line.strip()]


def main():
    a = sys.argv[1:]
    if a[0] != "--launch":
        sys.exit(__doc__)
    launcher, flags, npflag, np = a[1:5]
    launch = [launcher] + [f for f in flags.split(";") if f] + [npflag]
    exe, prob, work = a[5:8]
    common = a[8:]
    A, B, C = (os.path.join(work, x) for x in "ABC")
    run(launch + ["1", exe, prob, A] + common + ["--set", "amrex_pk=host"], A)
    run(launch + ["1", exe, prob, B] + common + ["--set", "amrex_pk=level"], B)
    run(launch + [np, exe, prob, C] + common + ["--set", "amrex_pk=level"], C)
    bad = []
    fa = files_of(A)
    npk = 0
    for rel in fa:
        pa, pb = os.path.join(A, rel), os.path.join(B, rel)
        if not os.path.exists(pb):
            bad.append(rel + ": missing in B")
            continue
        if rel.startswith("Pk/"):
            npk += 1
            xa, xb = numbers(pa), numbers(pb)
            if len(xa) != len(xb) or any(len(r) != len(s) for r, s in zip(xa, xb)):
                bad.append(rel + ": shape")
                continue
            ncol = len(xa[0])
            for c in range(ncol):
                scale = max(abs(r[c]) for r in xa) or 1.0
                worst = max(abs(r[c] - s[c]) for r, s in zip(xa, xb)) / scale
                if worst > 1e-4:
                    bad.append(f"{rel}: column {c + 1} differs by {worst:.2e} (relative to its largest value)")
        elif open(pa, "rb").read() != open(pb, "rb").read():
            bad.append(rel + ": A and B differ")
    if npk == 0:
        bad.append("no Pk/ files")
    for rel in files_of(B):
        pc_ = os.path.join(C, rel)
        if not os.path.exists(pc_) or open(pc_, "rb").read() != open(os.path.join(B, rel), "rb").read():
            bad.append(rel + ": B and C (" + np + " ranks) differ")
    for b in bad[:40]:
        print("  FAIL", b)
    if bad:
        sys.exit(f"{len(bad)} problems")
    print(f"{len(fa)} files: A = B except {npk} Pk/ files (within 1e-4), B = C on {np} ranks")


if __name__ == "__main__":
    main()
