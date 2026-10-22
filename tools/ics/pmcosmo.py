"""Cosmology helpers shared by make_ics.py and the PM physics checks.

The parameters reproduce reference/Radhydro/cosmo_parameters.f90, including
its single-precision literals: ``omegam = 0.27`` is a real(4) literal stored
in a real(8) parameter, so the code really uses float32(0.27). They are the
defaults of the runtime cosmology (params.txt keys omegam, omegab, hubble0,
sigma8, ns; docs/production_run.md section 10): ``set_cosmology`` replaces
them, rounded to float32 as the code rounds the params.txt values, and
``cosmology()`` returns the values in use.

Units are the code's "supercomoving" grid units (cosmo_tools.f90:
calc_sim_var):

* positions in grid cells (DM grid for particles, gas grid for gas);
* time tau with time_unit = 2/(3 H0 h sqrt(om)) a^2, so that
  da/dtau = 2/(3 sqrt(om)) a^3 E(a) (calc_da_from_dt);
* velocities v = dx/dtau in grid cells per time unit;
* density in units of the mean total matter density.

In these units the PM equations are d2x/dtau2 = -grad(phi),
lap(phi) = (2a/3) (rho - 1) (pm_tools.f90: Aphi = -4 pi a/(6 pi)), and the
linear growth factor obeys D'' = (2/3) c a D, where c is the fraction of the
mean density that clusters (c = 1 if all matter moves; c = 1 - Omega_b/Omega_m
in a gravity-only run where the gas density is never advected).
"""

import numpy as np

f32 = lambda x: float(np.float32(x))  # noqa: E731  real(4) literal -> real(8)

BOX = 5.0                     # Mpc/h (cosmo_parameters.f90: box)
ZINIT = 300.0
OMEGAM = f32(0.27)
OMEGAL = 1.0 - OMEGAM
OMEGAB = f32(0.045)
HUBBLE0 = f32(0.70)
NS = f32(0.96)
SIGMA8 = f32(0.80)
TCMB0 = f32(2.725)
OMEGAR = f32(4.48e-7) * (1 + f32(0.69)) * TCMB0 ** 4 / HUBBLE0 ** 2
FB = OMEGAB / OMEGAM

COSMO_KEYS = ("omegam", "omegab", "hubble0", "sigma8", "ns")


def cosmology():
    """The cosmology in use, as params.txt keys -> values (float32-rounded)."""
    return {"omegam": OMEGAM, "omegab": OMEGAB, "hubble0": HUBBLE0,
            "sigma8": SIGMA8, "ns": NS}


def set_cosmology(omegam=None, omegab=None, hubble0=None, sigma8=None, ns=None):
    """Replace the cosmology of this module (None keeps a value).

    Every value is rounded to float32, as radhydro rounds the params.txt keys
    (the original's real(4) literals), so the ICs and the run agree exactly
    and ics_info.txt can be checked against params.txt. The derived values
    follow cosmo_parameters.f90: OMEGAL = 1 - OMEGAM (flat, OMEGAR not
    subtracted), OMEGAR = 4.48e-7 (1 + 0.69) Tcmb0^4 / h^2.
    """
    global OMEGAM, OMEGAL, OMEGAB, HUBBLE0, NS, SIGMA8, OMEGAR, FB
    if omegam is not None:
        OMEGAM = f32(omegam)
    if omegab is not None:
        OMEGAB = f32(omegab)
    if hubble0 is not None:
        HUBBLE0 = f32(hubble0)
    if sigma8 is not None:
        SIGMA8 = f32(sigma8)
    if ns is not None:
        NS = f32(ns)
    if not (0 < OMEGAM <= 1 and 0 < OMEGAB < OMEGAM and 0 < HUBBLE0 < 2
            and SIGMA8 > 0 and NS > 0):
        raise ValueError(f"bad cosmology {cosmology()}")
    OMEGAL = 1.0 - OMEGAM
    OMEGAR = f32(4.48e-7) * (1 + f32(0.69)) * TCMB0 ** 4 / HUBBLE0 ** 2
    FB = OMEGAB / OMEGAM


def set_cosmology_from_info(path):
    """The cosmology recorded in an ics_info.txt (make_ics.py) or a
    params.txt; keys that are absent keep the original's values."""
    vals = {}
    with open(path) as f:
        for line in f:
            line = line.split("#", 1)[0]
            if "=" not in line:
                continue
            k, v = (t.strip() for t in line.split("=", 1))
            if k in COSMO_KEYS:
                vals[k] = float(v.replace("d", "e").replace("D", "e"))
    set_cosmology(**vals)
    return cosmology()


def E(a):
    """H(a)/H0 as in calc_sim_var/calc_da_from_dt (wde = -1)."""
    return np.sqrt(OMEGAM / a ** 3 + OMEGAL + OMEGAR / a ** 4)


def dlna_dtau(a):
    """d ln a / d tau in code time units."""
    return 2.0 / (3.0 * np.sqrt(OMEGAM)) * a ** 2 * E(a)


def growth_table(cluster=1.0, a_min=1e-8, a_max=1.0, n=40000):
    """Growing mode of D'' = (2/3) cluster a D (tau derivatives).

    Integrated in y = ln a with RK4 from deep in the radiation era, starting
    on the constant mode (D = 1, dD/dy = 0), which becomes the Meszaros
    growing mode; the decaying mode is absent up to O(a_min/a_eq).

    Returns (a, D, dD/dtau), D normalised to D(a_max) = 1 when a_max = 1.
    """
    y = np.linspace(np.log(a_min), np.log(a_max), n)
    h = y[1] - y[0]

    def rhs(yy, s):
        a = np.exp(yy)
        g = dlna_dtau(a)
        # g' = dg/dy (numerical derivative of an analytic function)
        eps = 1e-6
        gy = (dlna_dtau(a * np.exp(eps)) - dlna_dtau(a * np.exp(-eps))) / (2 * eps)
        D, Dy = s
        Dyy = ((2.0 / 3.0) * cluster * a * D - g * gy * Dy) / g ** 2
        return np.array([Dy, Dyy])

    s = np.array([1.0, 0.0])
    out = np.empty((n, 2))
    out[0] = s
    for i in range(n - 1):
        k1 = rhs(y[i], s)
        k2 = rhs(y[i] + h / 2, s + h / 2 * k1)
        k3 = rhs(y[i] + h / 2, s + h / 2 * k2)
        k4 = rhs(y[i] + h, s + h * k3)
        s = s + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        out[i + 1] = s
    a = np.exp(y)
    D = out[:, 0]
    Dtau = out[:, 1] * dlna_dtau(a)
    norm = D[-1]
    return a, D / norm, Dtau / norm


class Growth:
    """Interpolated growth factor D(a) and dD/dtau(a) (log-log in a)."""

    def __init__(self, cluster=1.0):
        self.cluster = cluster
        self.a, self.D, self.Dtau = growth_table(cluster)
        self._la = np.log(self.a)

    def __call__(self, a):
        return np.exp(np.interp(np.log(a), self._la, np.log(self.D)))

    def dtau(self, a):
        return np.interp(np.log(a), self._la, self.Dtau)

    def f(self, a):
        """dlnD/dlna."""
        return self.dtau(a) / self(a) / dlna_dtau(a)


def transfer_eh_nowiggle(k):
    """Eisenstein & Hu (1998) zero-baryon-oscillation transfer function.

    k in h/Mpc. Eqs. (26), (28)-(31) of EH98, with Theta_2.7 = Tcmb0/2.7.
    """
    om, ob, h = OMEGAM, OMEGAB, HUBBLE0
    theta = TCMB0 / 2.7
    omh2 = om * h * h
    obh2 = ob * h * h
    fb = ob / om
    s = 44.5 * np.log(9.83 / omh2) / np.sqrt(1 + 10 * obh2 ** 0.75)   # Mpc
    alpha = 1 - 0.328 * np.log(431 * omh2) * fb + 0.38 * np.log(22.3 * omh2) * fb ** 2
    kh = k * h                                                        # 1/Mpc
    gamma_eff = om * h * (alpha + (1 - alpha) / (1 + (0.43 * kh * s) ** 4))
    q = k * theta ** 2 / gamma_eff
    L0 = np.log(2 * np.e + 1.8 * q)
    C0 = 14.2 + 731.0 / (1 + 62.5 * q)
    return L0 / (L0 + C0 * q * q)


def power_unnormalised(k):
    return k ** NS * transfer_eh_nowiggle(k) ** 2


def sigma_R(R, pk, kmin=1e-5, kmax=1e3, n=200000):
    """rms of the top-hat filtered field (R in Mpc/h, pk(k) in (Mpc/h)^3)."""
    lk = np.linspace(np.log(kmin), np.log(kmax), n)
    k = np.exp(lk)
    x = k * R
    W = 3 * (np.sin(x) - x * np.cos(x)) / x ** 3
    integrand = k ** 3 * pk(k) * W ** 2 / (2 * np.pi ** 2)
    return np.sqrt(np.trapezoid(integrand, lk) if hasattr(np, "trapezoid")
                   else np.trapz(integrand, lk))


def power_z0(k):
    """Linear P(k) at z = 0, normalised to sigma8, (Mpc/h)^3."""
    A = (SIGMA8 / sigma_R(8.0, power_unnormalised)) ** 2
    return A * power_unnormalised(k)
