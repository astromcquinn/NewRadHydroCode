#!/usr/bin/env python3
"""Exact linear transfer function of the discrete first-order moment scheme
in the uniform mean free path test (docs/rt_moments.md section 7.1).

usage: moment_transfer_exact.py N frac [nsub] [--code TABLE] > out.txt
lambda = N / frac cells; nsub subcycles per substep (default 4). The table is
on the k shells of moment_transfer.py. With --code TABLE (a table written by
moment_transfer.py for the same N and lambda) it adds the ratio of the code's
T to this solution for every column of TABLE whose name contains m1 or otvet
and inj1, inj0 or inj2 (first match; names containing _o2 are skipped).

One RT substep of the code (tools/analysis/moment_transfer.cpp, the transport
kernels of the solver), in cell units:

  nsub subcycles with C = 1/nsub of
      U' = U - C sum_d [ (f_d(U_R) - f_d(U_L))/2 - (U_R - 2 U + U_L)/2 ],
  U = (N, F), f_d = (F_d, P_jd), P = D N; then (or before, moment_inject = 0)
  the injection N += dS; then the chemistry sees N; then N, F *= 1 - f,
  f = 1 - exp(-1/lambda).

Linearised about the uniform state (D = I/3 + dh, dh the response of the
OTVET tensor to the sources; dh = 0 for M1, whose chi = 1/3 + O(F^2/N^2)),
one subcycle acts on u = (dN, phi = s.F) of a Fourier mode as

  u' = m u + b,   m = [[a, -iC], [-iC |s|^2/3, a]],   b = (0, -iC Nbar w dL),
  a = 1 - C sum_d (1 - cos k_d),  s_d = sin k_d,  w = s_j s_d dh_jd / dL,

so a substep is u -> M u + G b with M = m^nsub, G = 1 + m + ... + m^(nsub-1).
With t = u / (Nbar_seen dL), Nbar_seen = dSbar / f, the steady states are

  inject after:   t = (1-f) M t + (1-f) G bh + f e_N       (Nbar = (1-f) Nbar_seen in the transport)
  inject before:  t = (1-f) M t + G bh + f M e_N           (Nbar = Nbar_seen in the transport)

  first cell:     t = (1-f) M t + (1-f) G bh + f (W, s.m)  (moment_inject = 2: the injection of
                  "inject after" spread by the weights of mk::first_cell_weights, W(k) = sum_o w_o
                  cos k.x_o in N and s.m(k) = -i sum_o (s.m_o) sin k.x_o in phi)

with bh = (0, -iC w), and T = |t_N|. For M1 (w = 0) this is
T_1 = f [(1 - (1-f) M)^-1]_NN and T_0 = f [(1 - (1-f) M)^-1 M]_NN.

The tensor response is that of the code's kernel (otvet_kernel of
rt_moments.cpp: K_jd = d_j d_d / |d|^4 at the nearest image, 4/3 on the
diagonal at d = 0, off-diagonal terms zero on the half-box planes):
dh_jd(k) = [K_jd(k) - (delta_jd / 3) tr K(k)] / tr K(0) dL(k).
"""
import re
import sys

import numpy as np

args = sys.argv[1:]
code = None
if "--code" in args:
    i = args.index("--code")
    code = args[i + 1]
    del args[i:i + 2]
if len(args) < 2:
    sys.exit(__doc__)
N, frac = int(args[0]), int(args[1])
nsub = int(args[2]) if len(args) > 2 else 4
lam = N / frac
f = 1 - np.exp(-1 / lam)
C = 1.0 / nsub

k1 = 2 * np.pi * np.fft.fftfreq(N)
kz, ky, kx = np.meshgrid(k1, k1, k1, indexing="ij")  # axes (z, y, x) as moment_transfer.py
kk = np.sqrt(kx**2 + ky**2 + kz**2)
sx, sy, sz = np.sin(kx), np.sin(ky), np.sin(kz)
s2 = sx**2 + sy**2 + sz**2
a = 1 - C * ((1 - np.cos(kx)) + (1 - np.cos(ky)) + (1 - np.cos(kz)))


def tensor_response():
    """w(k) = s_j s_d dh_jd / dL from the code's kernel."""
    d1 = np.fft.fftfreq(N, 1.0 / N)  # nearest-image offsets 0, 1, .., -1 (N/2 -> -N/2: even in d)
    z, y, x = np.meshgrid(d1, d1, d1, indexing="ij")
    r2 = x * x + y * y + z * z
    r2[0, 0, 0] = 1.0
    wgt = 1.0 / (r2 * r2)
    tx, ty, tz = 2 * np.abs(x) == N, 2 * np.abs(y) == N, 2 * np.abs(z) == N
    comp = {"xx": wgt * x * x, "yy": wgt * y * y, "zz": wgt * z * z,
            "xy": np.where(tx | ty, 0.0, wgt * x * y), "xz": np.where(tx | tz, 0.0, wgt * x * z),
            "yz": np.where(ty | tz, 0.0, wgt * y * z)}
    for c in ("xx", "yy", "zz"):
        comp[c][0, 0, 0] = 4.0 / 3
    for c in ("xy", "xz", "yz"):
        comp[c][0, 0, 0] = 0.0
    K = {c: np.fft.fftn(v).real for c, v in comp.items()}  # even kernels: real spectra
    tr = K["xx"] + K["yy"] + K["zz"]
    tr0 = tr[0, 0, 0]
    w = (sx * sx * (K["xx"] - tr / 3) + sy * sy * (K["yy"] - tr / 3) + sz * sz * (K["zz"] - tr / 3)
         + 2 * (sx * sy * K["xy"] + sx * sz * K["xz"] + sy * sz * K["yz"])) / tr0
    return w


def first_cell():
    """(w[3,3,3], m[3,3,3,3]) indexed [dz+1, dy+1, dx+1]: mk::first_cell_weights (kernels/moment_cell.hpp),
    the same midpoint rule over 6 nq^2 directions."""
    nq = 512
    h = 2.0 / nq
    g = -1.0 + (np.arange(nq) + 0.5) * h
    u, v = np.meshgrid(g, g, indexing="ij")
    w = np.zeros((3, 3, 3))
    m = np.zeros((3, 3, 3, 3))
    tot = 0.0
    for face in range(6):
        ax, sgn = face // 2, (-1.0 if face % 2 else 1.0)
        p = [None] * 3
        p[ax], p[(ax + 1) % 3], p[(ax + 2) % 3] = np.full_like(u, sgn), u, v
        r = np.sqrt(p[0]**2 + p[1]**2 + p[2]**2)
        wq = h * h / r**3
        d = np.stack([q / r for q in p])                       # (3, nq, nq)
        t = np.where(np.abs(d) > 0.5, 0.5 / np.maximum(np.abs(d), 1e-300), np.inf)
        order = np.argsort(t, axis=0)
        ts = np.take_along_axis(t, order, axis=0)
        off = np.zeros(d.shape, dtype=int)
        t0 = np.zeros_like(u)
        for l in range(4):
            t1 = np.minimum(ts[l], 1.0) if l < 3 else np.ones_like(u)
            ln = np.maximum(t1 - t0, 0.0) * wq
            idx = (off[2] + 1, off[1] + 1, off[0] + 1)
            np.add.at(w, idx, ln)
            for q in range(3):
                np.add.at(m[..., q], idx, ln * d[q])
            if l < 3:
                hit = np.isfinite(ts[l])
                for q in range(3):
                    sel = hit & (order[l] == q)
                    off[q][sel] = np.sign(d[q][sel]).astype(int)
                t0 = np.where(hit, ts[l], 1.0)  # no further face: the path ended at r = 1
        tot += wq.sum()
    return w / tot, m / tot


def first_cell_spectra():
    """W(k) (real) and s.m(k) (as a complex array) of the first-cell injection"""
    w, m = first_cell()
    W = np.zeros_like(a)
    sm = np.zeros_like(a, dtype=complex)
    for dz in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dx in (-1, 0, 1):
                ph = kx * dx + ky * dy + kz * dz
                W += w[dz + 1, dy + 1, dx + 1] * np.cos(ph)
                mo = m[dz + 1, dy + 1, dx + 1]
                sm += -1j * (sx * mo[0] + sy * mo[1] + sz * mo[2]) * np.sin(ph)
    return W, sm


def mat_mul(A, B):
    return (A[0] * B[0] + A[1] * B[2], A[0] * B[1] + A[1] * B[3], A[2] * B[0] + A[3] * B[2], A[2] * B[1] + A[3] * B[3])


def solve(closure):
    one, zero = np.ones_like(a, dtype=complex), np.zeros_like(a, dtype=complex)
    m = (a + 0j, -1j * C * one, -1j * C * s2 / 3, a + 0j)  # (NN, Np, pN, pp)
    M, G = (one, zero, zero, one), (zero, zero, zero, zero)
    for _ in range(nsub):
        G = tuple(g + x for g, x in zip(G, M))  # 1 + m + ... + m^(j-1)
        M = mat_mul(m, M)
    bh = -1j * C * tensor_response() if closure == "otvet" else None
    # A t = rhs, A = 1 - (1-f) M
    A = (1 - (1 - f) * M[0], -(1 - f) * M[1], -(1 - f) * M[2], 1 - (1 - f) * M[3])
    det = A[0] * A[3] - A[1] * A[2]
    out = {}
    W1, sm1 = first_cell_spectra()
    for inj in (1, 0, 2):
        rN = f * one if inj == 1 else f * W1 if inj == 2 else f * M[0]
        rp = zero if inj == 1 else f * sm1 if inj == 2 else f * M[2]
        if bh is not None:
            amp = (1 - f) if inj else 1.0
            rN = rN + amp * G[1] * bh
            rp = rp + amp * G[3] * bh
        out[inj] = np.abs((A[3] * rN - A[1] * rp) / det)
    return out


# the shells of moment_transfer.py
nbin = 14
edges = np.geomspace(2 * np.pi / N * 0.99, np.pi * 1.001, nbin + 1)
idx = np.digitize(kk.ravel(), edges) - 1
ok = (idx >= 0) & (idx < nbin)
cnt = np.bincount(idx[ok], minlength=nbin)


def binned(q):
    return np.bincount(idx[ok], weights=q.ravel()[ok], minlength=nbin) / np.maximum(cnt, 1)


kb = binned(kk)
cols = {}
for cl in ("m1", "otvet"):
    T = solve(cl)
    for inj in (1, 0, 2):
        cols["%s_inj%d" % (cl, inj)] = np.sqrt(binned(T[inj] ** 2))


def code_table(path):
    """the columns of a moment_transfer.py table: {name: T(k)}"""
    lines = open(path).read().split("\n")
    hdr = [ln for ln in lines if ln.startswith("# k*lam")][0]
    names = re.split(r"\s{2,}", hdr[2:].strip())
    rows = [re.sub(r"\([^)]*\)", " ", ln).split() for ln in lines if ln and not ln.startswith("#")]
    arr = np.array(rows, dtype=float)
    return {n: arr[:, i] for i, n in enumerate(names)}


ratio = {}
if code:
    tab = code_table(code)
    if len(tab["k*lam"]) != nbin or np.max(np.abs(tab["k*lam"] / (kb * lam) - 1)) > 0.01:
        sys.exit("%s is not on the shells of N = %d, lambda = %g" % (code, N, lam))
    for c in cols:
        cl, inj = c.split("_")
        for n in tab:
            if cl in n and inj in n and "_o2" not in n and "lin" not in n:
                ratio[c] = tab[n] / cols[c]
                break
names = list(cols) + ["code/" + c for c in ratio]
print("# exact linear transfer function of the discrete scheme (tools/analysis/moment_transfer_exact.py): "
      "N = %d, lambda = %g cells, nsub = %d" % (N, lam, nsub))
print("# k*lam  " + "  ".join("%15s" % c for c in names))
for b in range(nbin):
    print("%7.2f  " % (kb[b] * lam) + "  ".join("%15.4e" % cols[c][b] for c in cols)
          + "".join("  %15.3f" % ratio[c][b] for c in ratio))
for c in ratio:
    print("# code/%s: %.3f .. %.3f" % (c, ratio[c].min(), ratio[c].max()))
