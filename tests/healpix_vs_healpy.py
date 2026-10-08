#!/usr/bin/env python3
"""Compare radray's HEALPix NESTED routines with healpy.

usage: healpix_vs_healpy.py <path/to/healpix_tool> [nvec]

vec2pix: ~1e6 random (non-normalised) vectors spread over nside = 2^0..2^13;
         pixel numbers must agree exactly.
pix2vec: random pixels at every nside; vectors must agree to ~1e-15.
Exit 77 (CTest "skipped") if healpy is not installed.
"""
import os
import subprocess
import sys
import tempfile

import numpy as np

try:
    import healpy as hp
except ImportError:
    print("healpy not available, skipping")
    sys.exit(77)


def main():
    tool = sys.argv[1]
    nvec = int(sys.argv[2]) if len(sys.argv) > 2 else 1_000_000
    rng = np.random.default_rng(12345)
    nsides = [2**l for l in range(14)]
    fails = 0
    with tempfile.TemporaryDirectory() as td:
        # ---- vec2pix ----
        ns = rng.choice(nsides, size=nvec)
        v = rng.standard_normal((nvec, 3)) * rng.uniform(0.1, 10, size=(nvec, 1))
        # add some special directions: axes, near-poles, near z = +-2/3
        sp = []
        for nside in nsides:
            for vec in ([1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0], [0, 0, 1],
                        [0, 0, -1], [1, 1, 0], [1e-12, 0, 1], [0, 1e-9, -1],
                        [np.sqrt(5) / 3, 0, 2 / 3], [0, np.sqrt(5) / 3, -2 / 3]):
                sp.append([nside] + list(vec))
        sp = np.array(sp, dtype=np.float64)
        rec = np.concatenate([np.column_stack([ns.astype(np.float64), v]), sp])
        fin = os.path.join(td, "v.f64")
        fout = os.path.join(td, "p.i32")
        rec.tofile(fin)
        subprocess.check_call([tool, "vec2pix", fin, fout])
        mine = np.fromfile(fout, dtype=np.int32)
        ref = np.empty(len(rec), dtype=np.int64)
        for nside in nsides:
            sel = rec[:, 0] == nside
            ref[sel] = hp.vec2pix(nside, rec[sel, 1], rec[sel, 2], rec[sel, 3], nest=True)
        bad = np.nonzero(mine != ref)[0]
        print(f"vec2pix: {len(rec)} vectors, {len(bad)} mismatches")
        for b in bad[:20]:
            print("  mismatch", rec[b], mine[b], ref[b])
        fails += len(bad)

        # ---- pix2vec ----
        recs = []
        for nside in nsides:
            npix = 12 * nside * nside
            if npix <= 20000:
                p = np.arange(npix)
            else:
                p = rng.integers(0, npix, size=20000)
            recs.append(np.column_stack([np.full(len(p), nside), p]))
        recs = np.concatenate(recs).astype(np.int32)
        fin = os.path.join(td, "p.i32")
        fout = os.path.join(td, "v.f64")
        recs.tofile(fin)
        subprocess.check_call([tool, "pix2vec", fin, fout])
        mine = np.fromfile(fout, dtype=np.float64).reshape(-1, 3)
        ref = np.empty_like(mine)
        for nside in nsides:
            sel = recs[:, 0] == nside
            ref[sel] = np.array(hp.pix2vec(nside, recs[sel, 1], nest=True)).T
        err = np.abs(mine - ref).max()
        print(f"pix2vec: {len(recs)} pixels, max |dv| = {err:.3e}")
        if err > 1e-14:
            fails += 1
    if fails:
        print("FAILED")
        sys.exit(1)
    print("OK")


if __name__ == "__main__":
    main()
