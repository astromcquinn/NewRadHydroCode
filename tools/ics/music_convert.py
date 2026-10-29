#!/usr/bin/env python3
"""Convert MUSIC's generic HDF5 output (format = generic) into this code's
initial-condition files (docs/zoom_turbulence.md section 4.1, thread D;
tools/ics/music/README.md).

Input: the .hdf5 written by MUSIC, its .conf and the log MUSIC writes next to
the run (<conf>_log.txt); or an .npz packed from them (--pack, the fast-suite
test data: numpy only, no h5py).

Output, in OUTDIR (conventions of tools/ics/make_ics.py; float32, Fortran
order, no record markers):

  xv_dm.init       the base level (MUSIC's levelmin) as today: Ndm^3 = N^3
                   particles, one per base cell, at the cell centre plus the
                   displacement, x in base DM cells [0, N), v = dx/dtau in base
                   DM cells per code time unit; the order is the lattice's
                   (i fastest), i.e. particle id n + 1 sits at cell n.
  u_gas.init       the base level's gas: rho (units of the mean total matter
                   density, mean Omega_b/Omega_m), rho v (v in base gas cells
                   per time unit), rho v^2/2 (no thermal part).
  u_gas.init.L<l>  refined level l = 1 .. (levelmax - levelmin): the same five
                   fields on the level's whole grid (MUSIC's region for that
                   level, padding included); velocities in level-0 cell units,
                   as every level of the AMReX hierarchy keeps them
                   (src/amrex/hydro_level.cpp).
  xv_dm.L<l>       the particles of mass class l: one per level-l cell not
                   covered by level l + 1 (the finest level covering each
                   Lagrangian cell), in the level grid's order (i fastest)
                   with the covered cells skipped; units as xv_dm.init.
  ics_levels.txt   one block per level: the grid (offset in that level's
                   cells from the box origin, taken modulo the box; size),
                   the box of the next level in this level's cells, the
                   particle count and the class mass (code units: mean
                   matter density x base gas-cell volume; and M_sun/h).
                   The base particles under level 1 are superseded by the
                   finer classes; xv_dm.init still holds all N^3 of them so
                   that the code reads it as today.
  ics_info.txt     the parameters, the cosmology, the growth factors used.

Conversion (1LPT; Ndm = Ngas = N at every level):
  * densities: MUSIC writes the overdensity delta (CDM and baryons
    separately); the gas gets rho = f_b (1 + r_D delta_b);
  * displacements: MUSIC's dx, dy, dz are in box units; psi = r_D d N
    (base cells);
  * amplitude: r_D = D(a_init) / D_MUSIC(a_init), the code's growth factor
    (tools/ics/pmcosmo.py, all matter clustering, radiation in H) over
    MUSIC's (its ODE without radiation, Omega_r = 0 in MUSIC's log; computed
    here from the exact integral and checked against D+_start in the log).
    The fields then carry sigma_8 today evolved back with the code's own
    growth, the convention of make_ics.py;
  * velocities: MUSIC writes v = vfact psi_theta (vfact = f a H/h in
    km/s per Mpc/h, so the stored numbers are km/s per box length), where
    psi_theta is the displacement of the velocity field (the density's own
    for transfer functions without velocities, e.g. eisenstein; the CAMB
    v_CDM, v_b columns for camb_file, which already carry the growth rate
    f(k) at zstart). The code's velocities are
        v = F r_D (v_MUSIC / vfact_MUSIC) N,
        F = (dD/dtau)/D (the code's growth rate)        without velocity TFs,
        F = d ln a/d tau (the code's expansion rate)     with velocity TFs,
    so that the velocity-to-displacement ratio is the code's growing mode
    exactly when MUSIC has no velocity information, and CAMB's f(k) with the
    code's H(a) when it has.

Cosmology: until the code reads the cosmology at run time (thread C), the
target cosmology is a preset (--cosmology original | planck_c26) or the five
keys; the converter refuses a MUSIC run whose [cosmology] differs from it
(relative 1e-6, so the original's float literals 0.27f etc. match "0.27"),
and writes it into ics_info.txt for the run to check.

Usage:
  music_convert.py MUSIC.hdf5 OUTDIR --conf RUN.conf [--log RUN.conf_log.txt]
                   [--cosmology planck_c26] [--levels 0,1]
  music_convert.py --pack MUSIC.hdf5 OUT.npz --conf RUN.conf [--float32]
  music_convert.py OUT.npz OUTDIR [--cosmology ...]
"""

import argparse
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pmcosmo as pc  # noqa: E402

NGHOST = 4                              # MUSIC's ghost zones (generic plugin)
RHO_CRIT = 2.77536627e11                # h^2 M_sun / Mpc^3
COMP = ("x", "y", "z")


# ---------------------------------------------------------------- MUSIC files
def parse_conf(text):
    """MUSIC's config syntax: [section], key = value, '#' or ';' comments."""
    conf, sec = {}, None
    for line in text.splitlines():
        line = re.split(r"[#;]", line, maxsplit=1)[0].strip()
        if not line:
            continue
        m = re.match(r"^\[(.+)\]$", line)
        if m:
            sec = m.group(1).strip().lower()
            conf.setdefault(sec, {})
            continue
        if "=" in line and sec is not None:
            k, v = line.split("=", 1)
            conf[sec][k.strip()] = v.strip()
    return conf


def parse_log(text):
    """What the converter needs from MUSIC's log: the shift of the volume
    (MUSIC centres the zoom region unless no_shift = yes), D+_start, Omega_r,
    whether the transfer function has velocities."""
    text = re.sub(r"\x1b\[[0-9;]*m", "", text)
    out = {"shift": (0, 0, 0)}
    m = re.search(r"Domain will be shifted by \(\s*(-?\d+),\s*(-?\d+),\s*(-?\d+)\)", text)
    if m:
        out["shift"] = tuple(int(m.group(i)) for i in (1, 2, 3))
    m = re.search(r"D\+_start = ([0-9.eE+-]+)", text)
    if m:
        out["dplus_start"] = float(m.group(1))
    m = re.search(r"Omega_r\s*=\s*([0-9.eE+-]+)", text)
    if m:
        out["omega_r"] = float(m.group(1))
    m = re.search(r"version (\S+)", text)
    if m:
        out["version"] = m.group(1)
    return out


class MusicData:
    """The hierarchy of a generic MUSIC file: levels levelmin..levelmax, each
    with its absolute offset (own cells, MUSIC's shifted frame) and size; the
    fields read one at a time without the ghost zones."""

    def __init__(self, path, conf_path=None, log_path=None):
        self.path = path
        if path.endswith(".npz"):
            self._z = np.load(path)
            self.conf_text = str(self._z["conf_text"])
            self.log_text = str(self._z["log_text"])
            hdr = {k[7:]: self._z[k] for k in self._z.files if k.startswith("header/")}
            self.levelmin = int(self._z["levelmin"])
            self.levelmax = int(self._z["levelmax"])
            self._h5 = None
        else:
            import h5py   # only for MUSIC's own files
            if conf_path is None:
                sys.exit("--conf is needed with an .hdf5 file")
            self._h5 = h5py.File(path, "r")
            self.conf_text = open(conf_path).read()
            log_path = log_path or conf_path + "_log.txt"
            if not os.path.exists(log_path):
                sys.exit(f"MUSIC's log {log_path} not found (--log)")
            self.log_text = open(log_path).read()
            hdr = {k: self._h5["header"][k][()] for k in self._h5["header"]}
            self.levelmin = int(self._h5["header"].attrs["levelmin"])
            self.levelmax = int(self._h5["header"].attrs["levelmax"])
        self.conf = parse_conf(self.conf_text)
        self.log = parse_log(self.log_text)
        self.boxlength = float(self.conf["setup"]["boxlength"])
        self.zstart = float(self.conf["setup"]["zstart"])
        nl = self.levelmax - self.levelmin + 1
        off = np.array([np.asarray(hdr["grid_off_" + c]).reshape(-1)[:nl] for c in COMP], dtype=np.int64).T
        size = np.array([np.asarray(hdr["grid_len_" + c]).reshape(-1)[:nl] for c in COMP], dtype=np.int64).T
        self.header = hdr
        # grid_off is relative, in the parent's cells; absolute in own cells:
        # abs_l = 2 (abs_{l-1} + off_l)   (MUSIC mesh.hh: offset, offset_abs)
        self.offset_abs, self.size = {}, {}
        prev = np.zeros(3, dtype=np.int64)
        for i, lev in enumerate(range(self.levelmin, self.levelmax + 1)):
            a = np.zeros(3, dtype=np.int64) if i == 0 else 2 * (prev + off[i])
            self.offset_abs[lev], self.size[lev] = a, size[i].copy()
            prev = a
        n0 = 1 << self.levelmin
        if not np.all(self.size[self.levelmin] == n0):
            sys.exit(f"levelmin {self.levelmin} does not cover the box: size {self.size[self.levelmin]}")

    def has(self, lev, name):
        key = f"level_{lev:03d}_{name}"
        return key in (self._z.files if self._h5 is None else self._h5)

    def field(self, lev, name):
        key = f"level_{lev:03d}_{name}"
        if self._h5 is not None:
            a = self._h5[key][NGHOST:-NGHOST, NGHOST:-NGHOST, NGHOST:-NGHOST]
            return np.asarray(a, dtype=np.float64)
        a = self._z[key]
        if a.dtype.kind == "U":             # --pack stored a duplicate as a reference
            a = self._z[str(a)]
        return np.asarray(a, dtype=np.float64)


FIELDS = ["DM_rho", "DM_dx", "DM_dy", "DM_dz", "DM_vx", "DM_vy", "DM_vz",
          "BA_rho", "BA_vx", "BA_vy", "BA_vz"]


def pack(h5path, conf_path, log_path, out, float32=False, mantissa_bits=None):
    """Pack a (small) MUSIC run into an .npz the converter reads with numpy
    only: the conf, the log, the header and the fields of levelmin..levelmax
    without ghost zones (MUSIC's potential, the coarse copies below levelmin,
    and DM_rho when there are baryons are dropped; an array equal to one
    already stored, e.g. the gas velocities of a transfer function without
    velocities, is stored as a reference). float32 and mantissa_bits (the
    float32 mantissa rounded to that many bits, so that zlib compresses it)
    keep test data small."""
    md = MusicData(h5path, conf_path, log_path)
    baryons = md.conf["setup"].get("baryons", "no").lower() in ("yes", "true", "1")
    d = {"conf_text": np.array(md.conf_text), "log_text": np.array(md.log_text),
         "levelmin": np.array(md.levelmin), "levelmax": np.array(md.levelmax)}
    for k, v in md.header.items():
        d["header/" + k] = np.asarray(v)
    stored = []
    for lev in range(md.levelmin, md.levelmax + 1):
        for name in FIELDS:
            if not md.has(lev, name) or (baryons and name == "DM_rho"):
                continue
            key = f"level_{lev:03d}_{name}"
            a = md.field(lev, name)
            if float32 or mantissa_bits:
                a = a.astype(np.float32)
            if mantissa_bits:
                drop = 23 - mantissa_bits
                ai = a.view(np.uint32)
                ai[...] = ((ai.astype(np.uint64) + (1 << (drop - 1))) >> drop << drop).astype(np.uint32)
            same = next((k for k, b in stored if b.shape == a.shape and np.array_equal(a, b)), None)
            d[key] = np.array(same) if same else a
            if not same:
                stored.append((key, a))
    np.savez_compressed(out, **d)
    print(f"packed {h5path} -> {out} ({os.path.getsize(out)} bytes)")


# ------------------------------------------------------------- cosmologies
def music_growth(a, omegam, omegal):
    """MUSIC's linear growth (cosmology_calculator.hh: compute_growth, from
    a = 1e-10 with D = a, Omega_r = 0, w = -1): the growing mode
    D(a) = (5/2) Om E(a) int_0^a da'/(a' E)^3, normalised to D(1) = 1, and
    f = dlnD/dlna. Gauss-Legendre in t with a' = a t^(2/5) (the integrand is
    smooth in t), exact to round-off."""
    t, w = np.polynomial.legendre.leggauss(200)
    t, w = 0.5 * (t + 1), 0.5 * w

    def E(x):
        return np.sqrt(omegam / x ** 3 + omegal + (1 - omegam - omegal) / x ** 2)

    def integral(x):
        ap = x * t ** 0.4
        return np.sum(w * x * 0.4 * t ** (-0.6) / (ap * E(ap)) ** 3)

    def D(x):
        return 2.5 * omegam * E(x) * integral(x)

    Ea = E(a)
    dlnE = (-1.5 * omegam / a ** 3 - (1 - omegam - omegal) / a ** 2) / Ea ** 2
    f = dlnE + 1.0 / (a ** 2 * Ea ** 3 * integral(a))
    return float(D(a) / D(1.0)), float(f), float(100.0 * Ea)   # D/D0, f, H/h in km/s/Mpc


def target_cosmology(args):
    keys = [getattr(args, k) for k in pc.Cosmology.KEYS]
    if any(v is not None for v in keys):
        if not all(v is not None for v in keys):
            sys.exit("give all five of --omegam --omegab --hubble0 --sigma8 --ns, or a --cosmology preset")
        return pc.Cosmology(*keys, name="custom")
    return pc.PRESETS[args.cosmology]


def check_cosmology(md, cosmo):
    """Refuse a MUSIC run made for another cosmology."""
    c = md.conf.get("cosmology", {})
    have = {"omegam": float(c["Omega_m"]), "omegab": float(c["Omega_b"]),
            "hubble0": float(c["H0"]) / 100.0, "sigma8": float(c["sigma_8"]), "ns": float(c["n_s"])}
    bad = [f"{k}: MUSIC {have[k]!r}, target {getattr(cosmo, k)!r}" for k in have
           if abs(have[k] / getattr(cosmo, k) - 1) > 1e-6]
    ol = float(c.get("Omega_L", 1 - have["omegam"]))
    if abs(ol + have["omegam"] - 1) > 1e-6:
        bad.append(f"MUSIC is not flat (Omega_m + Omega_L = {ol + have['omegam']}); the code is")
    if bad:
        sys.exit("MUSIC's cosmology differs from the target (--cosmology / the five keys):\n  " + "\n  ".join(bad))
    if md.log.get("omega_r", 0.0) != 0.0:
        sys.exit(f"MUSIC ran with Omega_r = {md.log['omega_r']}: the converter's MUSIC growth assumes 0")
    return have["omegam"], ol


# --------------------------------------------------------------- conversion
def f32_fortran(path, arr):
    np.asarray(arr, dtype="<f4").ravel(order="F").tofile(path)


def convert(md, outdir, cosmo, levels=None):
    lmin, lmax = md.levelmin, md.levelmax
    N = 1 << lmin
    a = 1.0 / (1.0 + md.zstart)
    om_m, ol_m = check_cosmology(md, cosmo)

    # growth: the code's (all matter clustering) and MUSIC's
    g = pc.Growth(1.0, cosmo)
    D_code = float(g(a))
    vfac_code = float(g.dtau(a) / g(a))
    hub_code = float(pc.dlna_dtau(a, cosmo))
    D_m, f_m, H_m = music_growth(a, om_m, ol_m)
    if "dplus_start" in md.log and abs(md.log["dplus_start"] / D_m - 1) > 2e-5:
        sys.exit(f"MUSIC's D+_start {md.log['dplus_start']} against {D_m:.8g} computed here: "
                 "the growth model differs")
    vfact_m = f_m * a * H_m                    # km/s per Mpc/h (MUSIC: get_vfact)
    r_D = D_code / D_m
    transfer = md.conf["cosmology"].get("transfer", "")
    vel_tf = transfer.startswith("camb")       # MUSIC: tf_withvel_ (camb_file)
    F = hub_code if vel_tf else vfac_code
    baryons = md.conf["setup"].get("baryons", "no").lower() in ("yes", "true", "1")
    fb = cosmo.fb
    shift = np.array(md.log["shift"], dtype=np.int64)
    os.makedirs(outdir, exist_ok=True)
    levels = list(range(lmin, lmax + 1)) if levels is None else [lmin + int(x) for x in levels]

    # MUSIC's velocities against its displacements, every level: without
    # velocity TFs v = vfact psi (a units check; the zooms' multigrid solves,
    # accuracy 1e-5, leave a few 1e-6); with them (CAMB) v/psi is f(k), close
    # to 1, so the two must be tightly correlated. Guards against corrupt
    # velocities: MUSIC 2.0a (967651b) writes garbage or partly corrupted
    # velocities in some runs with levelmin = 4 (tools/ics/music/README.md).
    check = {"music_v_over_psi_rel_err": 0.0, "music_v_psi_corr_min": 1.0}
    for lev in range(lmin, lmax + 1):
        for c in COMP:
            d, v = md.field(lev, "DM_d" + c), md.field(lev, "DM_v" + c)
            corr = float(np.corrcoef(d.ravel(), v.ravel())[0, 1]) if d.std() > 0 and v.std() > 0 else float("nan")
            check["music_v_psi_corr_min"] = min(check["music_v_psi_corr_min"], corr) if corr == corr else corr
            if not vel_tf:
                sel = np.abs(d) > 0.1 * np.abs(d).max()
                err = float(np.max(np.abs(v[sel] / (vfact_m * d[sel]) - 1))) if sel.any() else float("nan")
                check["music_v_over_psi_rel_err"] = max(check["music_v_over_psi_rel_err"], err)
                if not err < 1e-3:
                    sys.exit(f"MUSIC's level {lev} DM v{c}/psi is not vfact = {vfact_m}: {err:.2e} "
                             "(corrupt velocities? see tools/ics/music/README.md)")
            elif not corr > 0.99:
                sys.exit(f"MUSIC's level {lev} DM v{c} and d{c} correlate only {corr:.4f} "
                         "(corrupt velocities? see tools/ics/music/README.md)")

    info_levels = []
    for lev in range(lmin, lmax + 1):
        L = lev - lmin
        sc = 1 << L                             # level-L cells per base cell
        n = md.size[lev]
        # back to the unshifted box: the base arrays are rolled, the refined
        # grids keep their data and move their offset (it may wrap)
        off = md.offset_abs[lev] - shift * sc if L > 0 else np.zeros(3, dtype=np.int64)
        off_mod = np.mod(off, N * sc)

        def get(name):
            x = md.field(lev, name)
            if L == 0 and np.any(shift):
                x = np.roll(x, tuple(int(-s) for s in shift), axis=(0, 1, 2))
            return x

        nxt = None
        if lev < lmax:
            nxt = (md.offset_abs[lev + 1] // 2 - md.offset_abs[lev], md.size[lev + 1] // 2)
        covered = np.zeros(tuple(n), dtype=bool)
        if nxt is not None:
            o, s = nxt
            covered[o[0]:o[0] + s[0], o[1]:o[1] + s[1], o[2]:o[2] + s[2]] = True
        npart = int(n.prod() - covered.sum()) if L > 0 else int(N ** 3)
        mass = (1 - fb) * 8.0 ** (-L)
        mass_msun = RHO_CRIT * cosmo.omegam * (md.boxlength / N) ** 3 * (1 - fb) * 8.0 ** (-L)
        info_levels.append(dict(level=L, music_level=lev, offset=off_mod, size=n, next=nxt,
                                npart=npart, mass=mass, mass_msun=mass_msun))
        if lev not in levels:
            continue

        # particles: lattice (cell centres, base-cell units) + displacement
        q = [(off[c] + np.arange(n[c]) + 0.5) / sc for c in range(3)]
        Q = np.meshgrid(*q, indexing="ij")
        keep = None if L == 0 else ~covered
        cnt = npart
        xv = np.empty((6, cnt), dtype=np.float32)
        for c in range(3):
            psi = r_D * N * get("DM_d" + COMP[c])
            vel = F * r_D * N * get("DM_v" + COMP[c]) / vfact_m
            pos = np.mod(Q[c] + psi, N)
            if keep is None:
                xv[c] = pos.ravel(order="F")
                xv[3 + c] = vel.ravel(order="F")
            else:
                xv[c] = pos.ravel(order="F")[keep.ravel(order="F")]
                xv[3 + c] = vel.ravel(order="F")[keep.ravel(order="F")]
            del psi, vel, pos
        del Q
        f32_fortran(os.path.join(outdir, "xv_dm.init" if L == 0 else f"xv_dm.L{L}"), xv)
        del xv

        # gas on the level's grid
        u = np.empty((5,) + tuple(n), dtype=np.float32)
        if baryons:
            rho = fb * (1 + r_D * get("BA_rho"))
            vn = "BA_v"
        else:
            rho = fb * (1 + r_D * get("DM_rho"))
            vn = "DM_v"
        u[0] = rho
        ek = np.zeros(tuple(n))
        for c in range(3):
            vel = F * r_D * N * get(vn + COMP[c]) / vfact_m
            u[1 + c] = rho * vel
            ek += vel * vel
        u[4] = 0.5 * rho * ek
        del ek, rho
        f32_fortran(os.path.join(outdir, "u_gas.init" if L == 0 else f"u_gas.init.L{L}"), u)
        del u

    with open(os.path.join(outdir, "ics_levels.txt"), "w") as f:
        f.write("# written by tools/ics/music_convert.py; one block per level (0 = the base)\n"
                "# offset: the grid's first cell in this level's cells from the box origin (mod the box)\n"
                "# next: the next level's grid in this level's cells, relative to this grid's first cell\n"
                "# mass: per particle, in mean matter density x base gas-cell volume; mass_msun_h in M_sun/h\n"
                "# level 0's particles under 'next' are superseded by the finer classes\n")
        f.write(f"levels = {lmax - lmin + 1}\nN = {N}\nbox = {md.boxlength!r}\n")
        for d in info_levels:
            L = d["level"]
            f.write(f"\n[level {L}]\nmusic_level = {d['music_level']}\n")
            f.write("offset = " + " ".join(map(str, d["offset"])) + "\n")
            f.write("size = " + " ".join(map(str, d["size"])) + "\n")
            if d["next"] is not None:
                f.write("next_offset = " + " ".join(map(str, d["next"][0])) + "\n")
                f.write("next_size = " + " ".join(map(str, d["next"][1])) + "\n")
            f.write(f"gas_file = {'u_gas.init' if L == 0 else f'u_gas.init.L{L}'}\n")
            f.write(f"dm_file = {'xv_dm.init' if L == 0 else f'xv_dm.L{L}'}\n")
            f.write(f"npart = {d['npart']}\nmass = {d['mass']!r}\nmass_msun_h = {d['mass_msun']:.6e}\n")

    info = {"z_init": md.zstart, "a_init": a, "N": N, "Ngas": N, "Ndm": N, "box": md.boxlength,
            "levelmin": lmin, "levelmax": lmax, "transfer": transfer, "velocity_tf": int(vel_tf),
            "baryons": int(baryons), "shift": " ".join(map(str, shift)),
            "D_code": D_code, "D_music": D_m, "r_D": r_D, "f_code": float(g.f(a)), "f_music": f_m,
            "vfac_code": vfac_code, "dlna_dtau_code": hub_code, "vfact_music_kms_per_mpch": vfact_m,
            "velocity_factor": "vfac_code" if not vel_tf else "dlna_dtau_code"}
    info.update(check)
    with open(os.path.join(outdir, "ics_info.txt"), "w") as f:
        f.write("# written by tools/ics/music_convert.py\n")
        f.write(f"source = {os.path.basename(md.path)}\nmusic_version = {md.log.get('version', '?')}\n")
        f.write(f"cosmology = {cosmo.name}\n")
        for k, v in cosmo.as_dict().items():
            f.write(f"{k} = {v!r}\n")
        for k, v in info.items():
            f.write(f"{k} = {v!r}\n" if isinstance(v, float) else f"{k} = {v}\n")
        f.write("# params.txt for this run: N = Ndm = %d, box = %r, z_init = %r, and the five cosmology keys\n"
                % (N, md.boxlength, md.zstart))
    print(f"wrote {outdir}: N = {N}, levels {lmin}..{lmax}, r_D = {r_D:.6f}, "
          f"velocity factor {'d ln a/dtau' if vel_tf else '(dD/dtau)/D'} = {F:.6g}")
    return info, info_levels


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("music", help="MUSIC's generic .hdf5, or an .npz from --pack")
    ap.add_argument("out", help="OUTDIR (or the .npz with --pack)")
    ap.add_argument("--conf")
    ap.add_argument("--log")
    ap.add_argument("--pack", action="store_true", help="pack the MUSIC run into an .npz (test data)")
    ap.add_argument("--float32", action="store_true", help="--pack: store the fields in float32")
    ap.add_argument("--mantissa-bits", type=int, default=None, help="--pack: round float32 mantissas to this")
    ap.add_argument("--cosmology", choices=sorted(pc.PRESETS), default="planck_c26")
    for k in pc.Cosmology.KEYS:
        ap.add_argument("--" + k, type=float, default=None)
    ap.add_argument("--levels", default=None, help="comma-separated levels (0 = base) to write; default all")
    args = ap.parse_args(argv)
    if args.pack:
        pack(args.music, args.conf, args.log, args.out, args.float32, args.mantissa_bits)
        return
    md = MusicData(args.music, args.conf, args.log)
    convert(md, args.out, target_cosmology(args), None if args.levels is None else args.levels.split(","))


if __name__ == "__main__":
    main()
