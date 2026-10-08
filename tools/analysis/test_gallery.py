#!/usr/bin/env python3
"""Quick-look figures for the validation test problems.

usage: test_gallery.py RT_LEGACY RT_FIXED HYDRO_RUNS PM_RUNS OUTDIR

  RT_LEGACY / RT_FIXED : rt_driver output roots (one subdir per problem in
                         tests/problems) for legacy_rt=legacy_physics=1 and =0
  HYDRO_RUNS           : hydro_driver output root (sod_x, ... subdirs)
  PM_RUNS              : pm_driver output root (pancake_x subdir)
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "tests"))
import hydro_physics_check as hpc  # noqa: E402

N = 32


def field(d, name, nc=1):
    a = np.fromfile(os.path.join(d, name + ".bin"), "<f4")
    return a.reshape((nc, N, N, N) if nc > 1 else (N, N, N), order="F")


def last_step(run):
    return sorted(s for s in os.listdir(run) if s.startswith("step_"))[-1]


def rt_figure(leg, fix, out):
    probs = ["uniform_1src", "uniform_8src", "lognormal_16src"]
    fig, ax = plt.subplots(3, 4, figsize=(15, 10.5), constrained_layout=True)
    for r, p in enumerate(probs):
        dl = os.path.join(leg, p, "step_0010")
        df = os.path.join(fix, p, "step_0010")
        fHI_l, fHI_f = field(dl, "fHI"), field(df, "fHI")
        T_f, rho = field(df, "temp"), field(df, "rhog")
        # slice through the plane holding most ionised cells
        k = int(np.argmax((1 - fHI_f).sum(axis=(0, 1))))
        im = ax[r, 0].imshow(np.log10(rho[:, :, k].T), origin="lower", cmap="viridis")
        fig.colorbar(im, ax=ax[r, 0], shrink=0.8, label="log10 ρ/⟨ρ_b⟩")
        ax[r, 0].set_title(f"{p}: gas density (k={k + 1})")
        im = ax[r, 1].imshow(np.log10(np.maximum(fHI_f[:, :, k].T, 1e-6)), origin="lower",
                             cmap="magma", vmin=-6, vmax=0)
        fig.colorbar(im, ax=ax[r, 1], shrink=0.8, label="log10 x_HI")
        ax[r, 1].set_title("x_HI, step 10 (fixed)")
        im = ax[r, 2].imshow(T_f[:, :, k].T, origin="lower", cmap="inferno")
        fig.colorbar(im, ax=ax[r, 2], shrink=0.8, label="T [K]")
        ax[r, 2].set_title("temperature, step 10 (fixed)")
        # ionised-fraction histories, legacy vs fixed
        for root, lab, ls in ((leg, "legacy", "-"), (fix, "fixed", "--")):
            xs, xm = [], []
            for s in range(1, 11):
                d = os.path.join(root, p, f"step_{s:04d}")
                f, rg = field(d, "fHI"), field(d, "rhog")
                xs.append(1 - f.mean())
                xm.append(1 - (f * rg).sum() / rg.sum())
            ax[r, 3].plot(range(1, 11), xs, ls, c="C0", label=f"volume, {lab}")
            ax[r, 3].plot(range(1, 11), xm, ls, c="C1", label=f"mass, {lab}")
        ax[r, 3].set_xlabel("hydro step")
        ax[r, 3].set_ylabel("⟨x_HII⟩")
        ax[r, 3].legend(fontsize=7)
        ax[r, 3].set_title("ionised fraction history")
    fig.suptitle("RT test problems (N=32, z=10, 10 steps of da=1e-4): "
                 "legacy = original code, fixed = legacy_rt=legacy_physics=0")
    fig.savefig(os.path.join(out, "tests_rt.png"), dpi=80)
    plt.close(fig)


def ifront_figure(leg, fix, out):
    """I-front radius vs time for one source in uniform gas."""
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    for root, lab, ls in ((leg, "legacy", "o-"), (fix, "fixed", "s--")):
        R = []
        for s in range(1, 11):
            f = field(os.path.join(root, "uniform_1src", f"step_{s:04d}"), "fHI")
            R.append((3 * (1 - f).sum() / (4 * np.pi)) ** (1 / 3))
        ax[0].plot(range(1, 11), R, ls, label=lab)
    t = np.arange(1, 11)
    ax[0].plot(t, R[0] * t ** (1 / 3), "k:", label="∝ t^{1/3} (recombination-free)")
    ax[0].set_xlabel("hydro step")
    ax[0].set_ylabel("I-front radius (equivalent sphere) [cells]")
    ax[0].legend()
    ax[0].set_title("uniform_1src: I-front growth")
    d = os.path.join(fix, "uniform_1src", "step_0010")
    f = field(d, "fHI")
    T = field(d, "temp")
    c = np.unravel_index(np.argmin(f), f.shape)
    g = np.indices(f.shape) - np.array(c)[:, None, None, None]
    g = (g + N // 2) % N - N // 2
    r = np.sqrt((g ** 2).sum(axis=0)).ravel()
    ax[1].semilogy(r, f.ravel(), ".", ms=2, label="x_HI")
    ax[1].set_xlabel("distance from source [cells]")
    ax[1].set_ylabel("x_HI")
    a2 = ax[1].twinx()
    a2.plot(r, T.ravel(), ".", ms=2, c="C3")
    a2.set_ylabel("T [K]", color="C3")
    ax[1].set_xlim(0, 10)
    ax[1].set_title("uniform_1src step 10 (fixed): radial profiles")
    fig.savefig(os.path.join(out, "tests_rt_ifront.png"), dpi=80)
    plt.close(fig)


def sod_figure(hroot, out):
    fig, ax = plt.subplots(3, 3, figsize=(13, 9), constrained_layout=True, sharex=True)
    c_L = 2.0
    rL, rR = hpc.DGAS_AVG, 0.125 * hpc.DGAS_AVG
    pL = rL * c_L ** 2 / hpc.GAMMA
    pR = 0.1 * pL
    bulk = np.array([0.25, -0.15, 0.1])
    x = np.arange(N) + 0.5
    xf = np.linspace(0, N, 1000)
    for axis, name in enumerate(("sod_x", "sod_y", "sod_z")):
        run = os.path.join(hroot, name)
        prof, t, _ = hpc.sod_profile(run, axis, N)
        vb = np.roll(bulk, axis)[axis]
        ex = np.empty((3, xf.size))
        s1, _ = hpc.exact_riemann(rL, 0, pL, rR, 0, pR, (xf - N / 2 - vb * t) / t)
        xm = np.where(xf < N / 2, xf, xf - N)
        s2, _ = hpc.exact_riemann(rR, 0, pR, rL, 0, pL, (xm - vb * t) / t)
        half = (xf >= N / 4) & (xf < 3 * N / 4)
        ex[:, half] = s1[:, half]
        ex[:, ~half] = s2[:, ~half]
        ex[1] += vb
        for q, lab in enumerate(("density", "velocity", "pressure")):
            ax[q, axis].plot(xf, ex[q], "k-", lw=1, label="exact Riemann")
            ax[q, axis].plot(x, prof[q], "o", ms=3, c="C0", label="code")
            ax[q, axis].set_ylabel(lab + " (code units)")
            if q == 0:
                ax[q, axis].set_title(f"{name} (t = {t:.3f}, bulk v = {vb:+.2f})")
        ax[2, axis].set_xlabel("cell along tube axis")
    ax[0, 0].legend()
    fig.suptitle("Sod shock tubes, N=32 (two interfaces, periodic box). "
                 "z differs from x/y: the moving-grid velocity is smoothed only along the "
                 "current sweep (quirk H2)")
    fig.savefig(os.path.join(out, "tests_hydro_sod.png"), dpi=80)
    plt.close(fig)


def pancake_figure(proot, out):
    import pm_physics_check as pmc
    run = os.path.join(proot, "pancake_x")
    prob = os.path.join(ROOT, "tests", "problems_pm", "pancake_x")
    p, info, ngas, ndm = pmc.problem(prob)
    acr = float(info["a_cross"])
    g = pmc.pc.Growth(1.0 if info.get("growth", "frozen-gas") == "full" else 1.0 - pmc.pc.FB)
    k = 2 * np.pi / ndm
    steps = pmc.steps(run)
    pick = [steps[0], steps[len(steps) // 2], steps[-1]]
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5), constrained_layout=True)
    q = np.arange(ndm) + 0.5
    qa = np.meshgrid(q, q, q, indexing="ij")[0].ravel(order="F")
    qs = np.linspace(0, ndm, 400)
    for i, d in enumerate(pick):
        st, xv, _, _ = pmc.load_step(d, ngas)
        a = st["scalefactor"]
        amp = g(a) / g(acr)
        vamp = g.dtau(a) / g(acr)
        sel = slice(0, ndm)  # one row of particles along x
        ax[0].plot(xv[0][sel], xv[3][sel], "o", ms=4, c=f"C{i}", label=f"a = {a:.4f} (D/D_cross = {amp:.2f})")
        ax[0].plot((qs - amp * np.sin(k * qs) / k) % ndm, -vamp * np.sin(k * qs) / k, "-", c=f"C{i}", lw=1)
        disp = (xv[0] - qa + ndm / 2) % ndm - ndm / 2
        ax[1].plot(qa[sel], disp[sel], "o", ms=4, c=f"C{i}")
        ax[1].plot(qs, -amp * np.sin(k * qs) / k, "-", c=f"C{i}", lw=1)
    ax[0].set_xlabel("x [DM cells]")
    ax[0].set_ylabel("v_x (code units)")
    ax[0].set_title("Zel'dovich pancake: phase space (points: code, lines: analytic)")
    ax[0].legend(fontsize=8)
    ax[1].set_xlabel("Lagrangian q [DM cells]")
    ax[1].set_ylabel("displacement x − q")
    ax[1].set_title("displacement vs analytic Zel'dovich")
    fig.savefig(os.path.join(out, "tests_pm_pancake.png"), dpi=80)
    plt.close(fig)


def main():
    leg, fix, hroot, proot, out = sys.argv[1:6]
    os.makedirs(out, exist_ok=True)
    rt_figure(leg, fix, out)
    ifront_figure(leg, fix, out)
    sod_figure(hroot, out)
    pancake_figure(proot, out)


if __name__ == "__main__":
    main()
