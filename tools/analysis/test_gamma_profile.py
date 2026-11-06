#!/usr/bin/env python3
"""Checks of gamma_profile.py on synthetic cubes (needs numpy, scipy,
matplotlib): the FFT stack against a direct average over the sources, one
streaming source against its exact shell means, and the command line on
two fake runs (the luv, rho and gamma source modes find the same sources).

usage: python3 tools/analysis/test_gamma_profile.py
"""
import os
import subprocess
import sys
import tempfile

import numpy as np
import scipy.fft as sfft

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gamma_profile as gpr  # noqa: E402


def mi_dist(N, c):
    d = [np.minimum(np.abs(np.arange(N) - c[a]), N - np.abs(np.arange(N) - c[a])) for a in range(3)]
    return np.sqrt(d[0][:, None, None] ** 2 + d[1][None, :, None] ** 2 + d[2][None, None, :] ** 2)


def direct(cells, g, N, rmax):
    edges = gpr.EDGES[gpr.EDGES <= rmax + 1e-9]
    out = np.zeros(len(edges) - 1)
    for c in cells:
        r = mi_dist(N, c)
        for s in range(len(edges) - 1):
            m = (r >= edges[s]) & (r < edges[s + 1])
            out[s] += g[m].mean()
    return out / len(cells)


def check_stack():
    rng = np.random.default_rng(3)
    N, rmax = 24, 10
    g = rng.random((N, N, N)).astype(np.float32)
    cells = rng.integers(0, N, size=(15, 3))
    idx, cnt, _, _ = gpr.shell_index(N, rmax)
    ghat = sfft.rfftn(g)
    fft = gpr.stacked(cells, ghat, N, idx, cnt, 1)
    ref = direct(cells, g, N, rmax)
    err = np.max(np.abs(fft - ref) / np.abs(ref))
    assert err < 1e-5, err
    print(f"stack vs direct: max rel. error {err:.1e}")


def check_one_source():
    N, rmax = 32, 14
    c = np.array([5, 17, 30])
    r = mi_dist(N, c)
    g = (1 + 50.0 / (r ** 2 + 0.3)).astype(np.float32)
    idx, cnt, _, edges = gpr.shell_index(N, rmax)
    p = gpr.stacked(c[None, :], sfft.rfftn(g), N, idx, cnt, 1)
    for s in range(len(edges) - 1):
        m = (r >= edges[s]) & (r < edges[s + 1])
        assert abs(p[s] / g[m].mean() - 1) < 1e-5, (s, p[s], g[m].mean())
    print("one source: shell means exact")


def write_run(root, name, g, fhi, rho, luv, z):
    d = os.path.join(root, name, "run", "step_0001")
    os.makedirs(d)
    with open(os.path.join(d, "stats.txt"), "w") as f:
        f.write(f"istep = 1\nredshift = {z}\n")
    g.astype(np.float32).tofile(os.path.join(d, "gammaHI_avg.bin"))
    fhi.astype(np.float32).tofile(os.path.join(d, "fHI.bin"))
    rho.astype(np.float32).tofile(os.path.join(d, "rhog.bin"))
    luv.astype(np.float32).tofile(os.path.join(d, "Luv.bin"))
    with open(os.path.join(root, name, "params.txt"), "w") as f:
        f.write("box = 10\n")
    return os.path.join(root, name, "run")


def check_cli():
    rng = np.random.default_rng(7)
    N = 32
    cells = np.array([[4, 4, 4], [20, 9, 27], [12, 25, 15], [28, 18, 6], [8, 14, 22], [24, 28, 12]])
    L = np.array([1.0, 3.0, 2.0, 6.0, 4.0, 5.0])
    stream = np.ones((N, N, N))
    diffuse = np.ones((N, N, N))
    rho = np.ones((N, N, N)) + 0.1 * rng.random((N, N, N))
    luv = np.zeros((N, N, N))
    for c, l in zip(cells, L):
        r = mi_dist(N, c)
        stream += 20 * l / (r ** 2 + 0.3)
        diffuse += 0.5 * l / (r + 1.0)
        rho[tuple(c)] = 300 + 10 * l
        luv[tuple(c)] = l
    fhi = np.zeros((N, N, N))
    with tempfile.TemporaryDirectory() as root:
        a = write_run(root, "rays", stream, fhi, rho, luv, 5.5)
        b = write_run(root, "m1", diffuse, fhi, rho, luv, 5.51)
        outs = {}
        for mode in ["luv", "rho", "gamma:3"]:
            out = os.path.join(root, "out_" + mode.replace(":", "_"))
            subprocess.run([sys.executable, os.path.join(HERE, "gamma_profile.py"), f"{a}:rays", f"{b}:M1",
                            "--outdir", out, "--z", "5.5", "--sources", mode, "--ngroups", "3", "--rmax", "14",
                            "--workers", "1"], check=True, stdout=subprocess.DEVNULL)
            rows = [l.split() for l in open(os.path.join(out, "gamma_profile.txt")) if not l.startswith("#")]
            assert os.path.exists(os.path.join(out, "gamma_profile.png"))
            outs[mode] = {(r[1], r[2], int(r[4])): float(r[6]) for r in rows}
            assert {r[3] for r in rows if r[2] == "all"} == {"6"}, mode
        for k, v in outs["luv"].items():
            assert abs(outs["rho"][k] / v - 1) < 1e-5, (k, v, outs["rho"][k])
            assert abs(outs["gamma:3"][k] / v - 1) < 1e-5, (k, v, outs["gamma:3"][k])
        top_rays = [outs["luv"][("rays", "top", s)] for s in range(4)]
        top_m1 = [outs["luv"][("M1", "top", s)] for s in range(4)]
        assert top_rays[0] > top_rays[1] > top_rays[2] > top_rays[3], top_rays
        assert top_rays[1] / top_m1[1] > 2, (top_rays, top_m1)
    print("command line: luv, rho and gamma sources agree; the streaming profile falls off")


if __name__ == "__main__":
    check_stack()
    check_one_source()
    check_cli()
    print("test_gamma_profile: all passed")
