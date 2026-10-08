#!/usr/bin/env python3
"""Generate the hydro test problems in tests/problems_hydro/ (numpy only).

usage: make_hydro_problems.py [outroot] [--force]

Each problem directory holds
  params.txt   N, z_init, nsteps, timestep (fixed|courant), da_sim
  u_gas.bin    float32 (5,N,N,N): rho, rho*(v-vg) (3), energy in the frame vg
  w_gas.bin    float32 (4,N,N,N): rho_HI, rho_HeI, rho_HeII, rho_Z
  v_gas.bin    float32 (3,N,N,N): grid (frame) velocity vg (optional;
               absent = 0, as init_hydro sets it)
  README.txt   what the problem is
all little-endian, Fortran order (component index fastest), in the code units
of hydro_tools.f90 (supercomoving): lengths in cells, density in units of the
mean MATTER density (mean gas density = Dgas_avg = omegab/omegam = 1/6),
velocities in cells per code time unit, energy density in rho*vel^2.
At z = 10: vel_unit = 134 km/s, 1 code unit of e/rho = 1.767e6 K (neutral).

Problems (see docs/port_notes_hydro.md, section "Test problems"):
  sod_x, sod_y, sod_z  Sod-like shock tube (rho 1 : 0.125, p 1 : 0.1) along
                       x, y or z; the box is periodic, so there are two
                       interfaces (at the cell faces 0 and N/2).  A uniform
                       bulk velocity is added (see below).  Courant step.
  sound_x              small-amplitude (1e-3) travelling sound wave, one
                       wavelength in the box, along x.  Fixed da.
  random3d             3-D log-normal density, compressive random velocity
                       (Mach ~ 20-40), a few near-empty cells below Dgas_min,
                       a few low-density fast cells (above c_limit), random
                       ionisation fractions and metals, non-zero initial
                       grid velocity.  Fixed (large) da.

Bulk velocity: every problem has a uniform velocity with three non-zero
components.  This is REQUIRED, not cosmetic: in the original advect (MF) the
number of sub-steps is nt = ceiling(amax*dt/cfl) with amax = max|vg| along
the pencil; for a pencil whose grid velocity is exactly zero nt = 0 and the
loop `if (it == nt) exit` never terminates (port_notes_hydro.md, H1).
"""

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))

F32 = np.float32
OMEGAM = float(F32(0.27))
OMEGAB = float(F32(0.045))
XHY = float(F32(0.76))
YHE = float(F32(0.24))
GAMMA = float(F32(5.0) / F32(3.0))   # gamma_ideal = 5./3 (real(4) division)
DGAS_AVG = OMEGAB / OMEGAM
DGAS_MIN = DGAS_AVG * 1e-4            # Dgas_avg*10**lgD_min


def write_f32(path, a):
    """a: array (nc, N, N, N) indexed [c, i, j, k]; write Fortran order."""
    np.asarray(a, dtype="<f4").ravel(order="F").tofile(path)


def write_problem(root, name, params, rho, vel, eth, vg=None, w=None, readme="",
                  force=False):
    d = os.path.join(root, name)
    if os.path.isdir(d) and not force:
        print(f"{name}: exists, skipped (use --force)")
        return
    os.makedirs(d, exist_ok=True)
    N = rho.shape[0]
    write_vg = vg is not None
    if vg is None:
        vg = np.zeros((3,) + rho.shape)
    if w is None:
        w = np.zeros((4,) + rho.shape)
        w[0] = XHY * rho
        w[1] = YHE * rho
    dv = vel - vg
    u = np.empty((5,) + rho.shape)
    u[0] = rho
    u[1:4] = rho * dv
    u[4] = eth + 0.5 * rho * (dv ** 2).sum(axis=0)
    write_f32(os.path.join(d, "u_gas.bin"), u)
    write_f32(os.path.join(d, "w_gas.bin"), w)
    if write_vg:
        write_f32(os.path.join(d, "v_gas.bin"), vg)
    with open(os.path.join(d, "params.txt"), "w") as f:
        f.write(f"# hydro test problem '{name}' (reference/harness_hydro/make_hydro_problems.py)\n")
        f.write(f"N = {N}\n")
        for k, v in params.items():
            f.write(f"{k} = {v}\n")
    with open(os.path.join(d, "README.txt"), "w") as f:
        f.write(readme.strip() + "\n")
    print(f"{name}: written to {d}")


def sod(axis, N=32):
    """Sod-like tube along axis (0,1,2); left state in cells 1..N/2."""
    c_L = 2.0                      # left sound speed [cells / code time]
    rho_L, rho_R = 1.0 * DGAS_AVG, 0.125 * DGAS_AVG
    p_L = rho_L * c_L ** 2 / GAMMA
    p_R = 0.1 * p_L
    x = np.arange(N)
    left = x < N // 2
    prof_rho = np.where(left, rho_L, rho_R)
    prof_p = np.where(left, p_L, p_R)
    shape = [1, 1, 1]
    shape[axis] = N
    rho = np.broadcast_to(prof_rho.reshape(shape), (N, N, N)).copy()
    p = np.broadcast_to(prof_p.reshape(shape), (N, N, N)).copy()
    bulk_x = np.array([0.25, -0.15, 0.1])   # bulk velocity for the x tube
    bulk = np.roll(bulk_x, axis)             # rotate with the problem
    vel = np.zeros((3, N, N, N))
    for c in range(3):
        vel[c] = bulk[c]
    eth = p / (GAMMA - 1)
    readme = f"""
Sod-like shock tube along {'xyz'[axis]} (N={N}, periodic, two interfaces at the
cell faces 0 and {N // 2}). Left (cells 1..{N // 2}): rho = Dgas_avg = 1/6,
p = rho*c_L^2/gamma with c_L = {c_L}; right: rho = 0.125*Dgas_avg, p = 0.1 p_L.
Uniform bulk velocity (vx,vy,vz) = {tuple(bulk)} (needed, see H1 in
docs/port_notes_hydro.md); the exact solution is the Riemann solution
advected with the bulk velocity. v_gas = 0 initially. Neutral w_gas.
Time step: the code's own Courant step (timestep = courant).
"""
    return rho, vel, eth, readme, bulk


def sound(N=32, A=1e-3):
    c0 = 2.0
    rho0 = DGAS_AVG
    p0 = rho0 * c0 ** 2 / GAMMA
    k = 2 * np.pi / N
    xc = np.arange(N) + 0.5          # cell centres, cells 1..N
    s = np.sin(k * xc)
    bulk = np.array([0.2, 0.3, -0.25])
    rho1 = rho0 * (1 + A * s)
    p1 = p0 * (1 + GAMMA * A * s)
    vx1 = bulk[0] + c0 * A * s
    rho = np.broadcast_to(rho1[:, None, None], (N, N, N)).copy()
    p = np.broadcast_to(p1[:, None, None], (N, N, N)).copy()
    vel = np.zeros((3, N, N, N))
    vel[0] = vx1[:, None, None]
    vel[1] = bulk[1]
    vel[2] = bulk[2]
    eth = p / (GAMMA - 1)
    readme = f"""
Travelling sound wave along x (N={N}): rho = rho0 (1 + A sin(kx)),
p = p0 (1 + gamma A sin(kx)), vx = vbx + c0 A sin(kx) with A = {A},
k = 2 pi/N, rho0 = Dgas_avg, c0 = sqrt(gamma p0/rho0) = {c0}.
Uniform bulk velocity (vbx,vby,vbz) = {tuple(bulk)}: the wave moves at
c0 + vbx. v_gas = 0 initially. Fixed da.
"""
    return rho, vel, eth, readme, bulk


def gaussian_field(N, rng, slope):
    """Real Gaussian random field with |delta_k|^2 ~ k^slope, unit variance."""
    k = np.fft.fftfreq(N) * N
    kx, ky, kz = np.meshgrid(k, k, k, indexing="ij")
    kk = np.sqrt(kx ** 2 + ky ** 2 + kz ** 2)
    amp = np.where(kk > 0, np.maximum(kk, 1e-12) ** (slope / 2), 0.0)
    amp[kk > N / 3] = 0.0
    noise = rng.normal(size=(N, N, N)) + 1j * rng.normal(size=(N, N, N))
    f = np.fft.ifftn(noise * amp).real
    f -= f.mean()
    return f / f.std(), (kx, ky, kz, kk, amp)


def random3d(N=32, seed=20121008, bulk=(10.0, -8.0, 6.0), vfast=(70.0, -55.0)):
    rng = np.random.default_rng(seed)
    # log-normal density, sigma_ln = 1.5
    g, _ = gaussian_field(N, rng, -2.0)
    sig = 1.5
    delta = np.exp(sig * g - sig ** 2 / 2)
    delta /= delta.mean()
    rho = DGAS_AVG * delta
    # compressive (potential) velocity field, rms ~ 2.5 per component
    phi, (kx, ky, kz, kk, amp) = gaussian_field(N, rng, -4.0)
    phik = np.fft.fftn(phi)
    vel = np.empty((3, N, N, N))
    for c, kc in enumerate((kx, ky, kz)):
        vel[c] = np.fft.ifftn(1j * kc * 2 * np.pi / N * phik).real
    vel *= 2.5 / vel.std()
    for c in range(3):
        vel[c] += bulk[c]
    # temperature: T = 300 K (rho/<rho>)^(2/3) (adiabat), e/rho in code units
    temp1_unit_z10 = 1.767159e6
    T = 300.0 * delta ** (2.0 / 3.0)
    eth = rho * T / temp1_unit_z10
    # near-empty cells (below Dgas_min -> smooth_cells), 6 of them
    flat = rng.choice(N ** 3, size=12, replace=False)
    empty, fast = flat[:6], flat[6:]
    ii = np.unravel_index(empty, (N, N, N), order="F")
    rho[ii] = 0.5 * DGAS_MIN
    eth[ii] = rho[ii] * 1.0 / temp1_unit_z10
    # low-density (below Dgas_frz = 0.01 Dgas_avg) fast hot cells -> c_limit
    ii = np.unravel_index(fast, (N, N, N), order="F")
    rho[ii] = 0.005 * DGAS_AVG
    vel[0][ii] = vfast[0]
    vel[1][ii] = vfast[1]
    eth[ii] = rho[ii] * 1e8 / temp1_unit_z10
    # random ionisation state and metals
    fHI = rng.uniform(0.0, 1.0, (N, N, N))
    fHeI = rng.uniform(0.0, 1.0, (N, N, N))
    fHeII = (1 - fHeI) * rng.uniform(0.0, 1.0, (N, N, N))
    Z = 0.02 * rng.uniform(0.0, 1.0, (N, N, N)) * (delta > 1)
    w = np.empty((4, N, N, N))
    w[0] = XHY * rho * fHI
    w[1] = YHE * rho * fHeI
    w[2] = YHE * rho * fHeII
    w[3] = Z * rho
    # initial grid velocity: a smoothed copy of the velocity (Gaussian, 4 cells)
    k = np.fft.fftfreq(N) * 2 * np.pi
    kx, ky, kz = np.meshgrid(k, k, k, indexing="ij")
    W = np.exp(-0.5 * (kx ** 2 + ky ** 2 + kz ** 2) * 4.0 ** 2)
    vg = np.empty_like(vel)
    for c in range(3):
        vg[c] = np.fft.ifftn(np.fft.fftn(vel[c]) * W).real
    readme = f"""
3-D random problem (N={N}, seed {seed}): log-normal density (sigma_ln = 1.5,
delta from {delta.min():.3g} to {delta.max():.3g}), compressive velocity
(potential flow, rms 2.5 per component, + bulk {bulk}),
T = 300 K (rho/<rho>)^(2/3) (Mach ~ 20-40), 6 cells at rho = Dgas_min/2
(smooth_cells), 6 low-density (0.005 Dgas_avg < Dgas_frz) cells moving at
{vfast} (x,y) with T ~ 1e8 K (above c_limit -> freezing branches), random
fHI, fHeI, fHeII and metallicity, v_gas = Gaussian-smoothed (4 cells)
velocity. Fixed da (large: several sub-cycles per sweep in early steps).
"""
    return rho, vel, eth, vg, w, readme


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    force = "--force" in sys.argv
    root = args[0] if args else os.path.join(ROOT, "tests", "problems_hydro")
    os.makedirs(root, exist_ok=True)
    for axis in range(3):
        rho, vel, eth, readme, _ = sod(axis)
        write_problem(root, f"sod_{'xyz'[axis]}",
                      {"z_init": 10, "nsteps": 5, "timestep": "courant", "da_sim": 0},
                      rho, vel, eth, readme=readme, force=force)
    rho, vel, eth, readme, _ = sound()
    write_problem(root, "sound_x",
                  {"z_init": 10, "nsteps": 24, "timestep": "fixed", "da_sim": 0.004},
                  rho, vel, eth, readme=readme, force=force)
    # large bulk velocity: the moving frame follows it, so advect sub-cycles
    rho, vel, eth, vg, w, readme = random3d()
    write_problem(root, "random3d",
                  {"z_init": 10, "nsteps": 10, "timestep": "fixed", "da_sim": 0.002},
                  rho, vel, eth, vg=vg, w=w, readme=readme, force=force)
    # small bulk velocity, very fast low-density cells: c_limit / freeze /
    # smooth_cells branches
    rho, vel, eth, vg, w, readme = random3d(bulk=(0.3, -0.2, 0.1), vfast=(100.0, -80.0))
    write_problem(root, "random3d_freeze",
                  {"z_init": 10, "nsteps": 10, "timestep": "fixed", "da_sim": 0.002},
                  rho, vel, eth, vg=vg, w=w, readme=readme, force=force)


if __name__ == "__main__":
    main()
