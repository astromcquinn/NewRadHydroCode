#!/usr/bin/env python3
"""The nside_merge key (Config::nside_merge, docs/rt_moments.md section 8.6).

usage: nside_merge_check.py WORKDIR DRIVER PROBLEM [driver args ...]

Runs DRIVER PROBLEM <WORKDIR>/<case> [driver args] for: no key, nside_merge = 2
(the default), 4 and 1, and requires
  - nside_merge = 2 byte-identical to no key (every output file but timing and
    logs, as thread_check.py): the default keeps the original;
  - nside_merge = 4 and 1 to change the fields (the key reaches the merge),
    with the volume-mean ionized fraction within 10 % of the default's;
  - nside_merge = 3 to be refused.
Exit 0 if all hold.
"""
import glob
import os
import shutil
import subprocess
import sys

import numpy as np

SKIP = ("timing.txt", "ray.log", "rt.log", "stdout.log")


def run(work, drv, prob, args, name, extra):
    o = os.path.join(work, name)
    if os.path.isdir(o):
        shutil.rmtree(o)
    os.makedirs(o)
    with open(o + ".log", "w") as f:
        r = subprocess.run([drv, prob, o] + args + extra, stdout=f, stderr=subprocess.STDOUT)
    return o, r.returncode


def files(o):
    return sorted(os.path.relpath(f, o) for f in glob.glob(os.path.join(o, "**", "*"), recursive=True)
                  if os.path.isfile(f) and os.path.basename(f) not in SKIP)


def same(a, b):
    fa = files(a)
    return bool(fa) and fa == files(b) and all(
        open(os.path.join(a, r), "rb").read() == open(os.path.join(b, r), "rb").read() for r in fa)


def xmean(o):
    # the last field file of the neutral fraction (fHI) or the ionized fraction
    c = sorted(glob.glob(os.path.join(o, "**", "*fHI*.bin"), recursive=True)) or \
        sorted(glob.glob(os.path.join(o, "**", "*xHII*.bin"), recursive=True))
    if not c:
        return None
    a = np.fromfile(c[-1], dtype=np.float32)
    return float(1 - a.mean()) if "fHI" in os.path.basename(c[-1]) else float(a.mean())


def main():
    work, drv, prob = sys.argv[1], sys.argv[2], sys.argv[3]
    args = sys.argv[4:]
    os.makedirs(work, exist_ok=True)
    out, ok = {}, True
    for name, extra in (("nokey", []), ("n2", ["--set", "nside_merge=2"]),
                        ("n4", ["--set", "nside_merge=4"]), ("n1", ["--set", "nside_merge=1"])):
        out[name], rc = run(work, drv, prob, args, name, extra)
        if rc != 0:
            print(open(out[name] + ".log").read()[-2000:])
            sys.exit(f"driver failed for {name}")
    _, rc = run(work, drv, prob, args, "n3", ["--set", "nside_merge=3"])
    print(f"nside_merge = 3 refused: {rc != 0}")
    ok &= rc != 0
    s2 = same(out["nokey"], out["n2"])
    print(f"nside_merge = 2 byte-identical to no key ({len(files(out['nokey']))} files): {s2}")
    ok &= s2
    x0 = xmean(out["nokey"])
    for n in ("n4", "n1"):
        d = not same(out["nokey"], out[n])
        x = xmean(out[n])
        close = x0 is not None and x is not None and abs(x / x0 - 1) < 0.1
        print(f"{n}: changes the outputs {d}; <x_HII> {x} against {x0}: within 10 % {close}")
        ok &= d and close
    print("PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
