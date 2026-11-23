"""Slices of x_HI, T, Delta_gas and |v| from the field outputs of a run
(step_NNNN/{fHI,temp,rhog,u_gas}.bin), one frame per output, and an mp4 of
them. For runs without movie frames (movie_dz = 0): the frames are only the
zdata1 outputs. usage: snapshot_slices.py RUN OUTDIR [--zmax 7.05] [--slice K]
"""
import argparse, glob, os, subprocess, sys
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(__file__))
from velocity_power import params_of, OM, H0, H0_CGS, MPC

ap = argparse.ArgumentParser()
ap.add_argument("run"); ap.add_argument("outdir")
ap.add_argument("--zmax", type=float, default=7.05)
ap.add_argument("--slice", type=int, default=None)
a = ap.parse_args()
P = params_of(a.run); N = int(P["N"]); box = float(P.get("box", 5))
om = float(np.float32(P.get("omegam", OM))); h0 = float(np.float32(P.get("hubble0", H0)))
k0 = N // 2 if a.slice is None else a.slice
os.makedirs(a.outdir, exist_ok=True)

def plane(path, ncomp=1):
    m = np.memmap(path, dtype=np.float32, mode="r")
    return np.array(m[k0 * N * N * ncomp:(k0 + 1) * N * N * ncomp]).reshape(N, N, ncomp) if ncomp > 1 \
        else np.array(m[k0 * N * N:(k0 + 1) * N * N]).reshape(N, N)

steps = []
for sd in sorted(glob.glob(os.path.join(a.run, "step_*"))):
    if not os.path.exists(os.path.join(sd, "u_gas.bin")): continue
    z = float([l.split("=")[1] for l in open(os.path.join(sd, "stats.txt")) if l.startswith("redshift")][0])
    if z <= a.zmax: steps.append((z, sd))
steps.sort(reverse=True)
ext = [0, box * 1e3, 0, box * 1e3]
frames = []
for i, (z, sd) in enumerate(steps):
    fhi, T, rho = plane(sd + "/fHI.bin"), plane(sd + "/temp.bin"), plane(sd + "/rhog.bin")
    u = plane(sd + "/u_gas.bin", 5)
    aa = 1 / (1 + z)
    vel_unit = (aa * box / h0 / N * MPC) / (2 / (3 * H0_CGS * h0 * np.sqrt(om)) * aa ** 2) / 1e5
    v = np.sqrt((u[..., 1:4] ** 2).sum(-1)) / u[..., 0] * vel_unit
    fig, ax = plt.subplots(1, 4, figsize=(16, 4.4))
    panels = [(np.log10(rho / rho.mean()), "log Δ_gas", "viridis", (-1, 2.5)),
              (np.log10(T), "log T [K]", "inferno", (2.5, 5)),
              (np.log10(np.maximum(fhi, 1e-8)), "log x_HI", "Blues", (-6, 0)),
              (v, "|v| [km/s]", "magma", (0, 30))]
    for p, (img, lab, cm, lim) in zip(ax, panels):
        im = p.imshow(img, origin="lower", extent=ext, cmap=cm, vmin=lim[0], vmax=lim[1])
        p.set_title(lab, fontsize=10); p.set_xlabel("kpc/h")
        fig.colorbar(im, ax=p, fraction=0.046, pad=0.03)
    ax[0].set_ylabel("kpc/h")
    fig.suptitle(f"{N}³, slice k = {k0}, z = {z:.2f}", fontsize=11)
    fig.tight_layout()
    f = os.path.join(a.outdir, f"slice_{i:03d}.png"); fig.savefig(f, dpi=110); plt.close(fig)
    frames.append(f); print(f, z, flush=True)
mp4 = os.path.join(a.outdir, "slices.mp4")
subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", "1.5", "-i",
                os.path.join(a.outdir, "slice_%03d.png"), "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2",
                "-pix_fmt", "yuv420p", mp4], check=True)
print("wrote", mp4)
