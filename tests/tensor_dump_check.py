#!/usr/bin/env python3
"""tensor_dump (docs/otvet_closure.md §1) is a pure diagnostic: run radhydro
with tensor_dump = 0 and 1 and require

  - every output other than the dump byte-identical between the two runs;
  - the dump written at each zdata1 snapshot (tensor.bin, tensor_odd.bin,
    otvet_h.bin, sources.txt), with tr(S2) = S0 and |S1| <= S0 to float
    precision in every cell above float32 underflow, and some radiation.

usage: tensor_dump_check.py WORKDIR RADHYDRO PROBLEM [radhydro args ...]
"""
import glob
import os
import shutil
import subprocess
import sys

import numpy as np

DUMP = ("tensor.bin", "tensor_odd.bin", "otvet_h.bin", "sources.txt")
SKIP = ("timing.txt", "ray.log", "rt.log", "stdout.log")


def run(drv, prob, out, args):
    if os.path.isdir(out):
        shutil.rmtree(out)
    os.makedirs(out)
    with open(out + ".log", "w") as f:
        r = subprocess.run([drv, prob, out] + args, stdout=f, stderr=subprocess.STDOUT)
    if r.returncode != 0:
        print(open(out + ".log").read()[-2000:])
        sys.exit(f"radhydro failed: {out}")


def main():
    work, drv, prob = sys.argv[1], sys.argv[2], sys.argv[3]
    args = sys.argv[4:]
    o0, o1 = os.path.join(work, "dump0"), os.path.join(work, "dump1")
    run(drv, prob, o0, args + ["--set", "tensor_dump=0"])
    run(drv, prob, o1, args + ["--set", "tensor_dump=1"])
    files = sorted(os.path.relpath(f, o1) for f in glob.glob(os.path.join(o1, "**", "*"), recursive=True)
                   if os.path.isfile(f) and os.path.basename(f) not in SKIP)
    bad = []
    dumps = []
    for rel in files:
        if os.path.basename(rel) in DUMP:
            if os.path.basename(rel) == "tensor.bin":
                dumps.append(os.path.dirname(rel))
            continue
        g = os.path.join(o0, rel)
        if not os.path.exists(g) or open(g, "rb").read() != open(os.path.join(o1, rel), "rb").read():
            bad.append(rel)
    extra = [os.path.relpath(f, o0) for f in glob.glob(os.path.join(o0, "**", "*"), recursive=True)
             if os.path.isfile(f) and os.path.basename(f) in DUMP]
    print(f"{len(files)} files, {len(dumps)} dumps: {' '.join(dumps)}")
    for x in bad[:20]:
        print("  DIFF", x)
    for x in extra[:5]:
        print("  dump written with tensor_dump = 0:", x)
    ok = not bad and not extra and dumps
    for d in dumps:
        for f in DUMP:
            if not os.path.exists(os.path.join(o1, d, f)):
                print("  missing", d, f)
                ok = False
        S = np.fromfile(os.path.join(o1, d, "tensor.bin"), dtype=np.float32).reshape(-1, 10).astype(np.float64)
        S0 = S[:, 0]
        # float32 underflow leaves a few denormal cells (S0 ~ 1e-45) whose
        # trace is meaningless; they have no weight in any S0-weighted metric
        lit = S0 > 1e-20 * S0.max()
        tr = S[:, 4] + S[:, 5] + S[:, 6]
        f = np.sqrt((S[:, 1:4] ** 2).sum(1))
        etr = np.abs(tr[lit] / S0[lit] - 1).max() if lit.any() else np.nan
        ef = (f[lit] / S0[lit]).max() if lit.any() else np.nan
        print(f"  {d}: {lit.mean():.3f} of the cells lit, max |tr h - 1| = {etr:.2e}, max |f| = {ef:.6f}")
        if not lit.any() or etr > 1e-5 or ef > 1 + 1e-5 or not np.isfinite(S).all():
            ok = False
    print("OK" if ok else "FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
