#!/usr/bin/env python3
"""Generate the full-loop validation problems (tests/problems_full).

usage: make_problems.py <build_dir> [outdir=tests/problems_full]

<build_dir> must contain the C++ radhydro and find_halos executables.

Writes
  ics/           xv_dm.init, u_gas.init: Zel'dovich ICs, Ngas = 32, Ndm = 16,
                 box = 4 Mpc/h, z_init = 100, seed 7 (tools/ics/make_ics.py)
  onthefly/      params.txt: on-the-fly FoF sources (C++ halo finder)
  legacy/        params.txt + Halo/: legacy catalogs halo.z=XX.XXXX and
                 sfr.txt, made with the C++ halo finder and its legacy writer
                 from DM snapshots of a gravity+hydro run of the same ICs
                 (radhydro with zdata2 outputs, before any RT)

  smoke/         a small, fast problem for the sanitizer and live Fortran
                 tests: Ngas = 16, Ndm = 8, box 4, RT from z ~ 30 with one
                 hand-made legacy source (Halo/, written here with the same
                 binary layout) or on-the-fly (no halos form that early)

Both main problems: RT starts after the first zdata1 output (z <= 8.5, the
original idata1 > 1 gate) and the run stops at the step that crosses
z = 7.5: 24 hydro steps, of which the last one does RT (~1670 RT
substeps, ~70-100k rays, about 4 minutes per run). See docs/full_run.md.
"""
import math
import os
import shutil
import struct
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))

BOX = 4
SEED = 7
Z_INIT = 100
NMIN = 8           # FoF minimum particles (Ndm = 16 resolves only massive halos)
ZDATA1 = "8.5 7.5"
CATALOG_Z = [10.0, 9.0, 8.0, 7.0, 6.0]   # sfr.txt rows / catalog redshifts

COMMON = f"""# Full radhydro loop validation problem (tests/problems_full/make_problems.py)
N = 32
Ndm = 16
box = {BOX}
nray_max_factor = 4
z_init = {Z_INIT}
zdata1 = {ZDATA1}
ics = ../ics
fof_min_npart = {NMIN}
legacy_rt = 1
legacy_physics = 1
"""


def hs_sfr(z, rho_m=0.013, alpha=0.012, beta=0.041):
    """Hernquist & Springel (2003) fit, Msun/yr/Mpc^3 (as SourceModelParams)."""
    om = 0.27
    E = math.sqrt(om * (1 + z) ** 3 + 1 - om)
    chi = E ** (2.0 / 3.0)
    return rho_m * chi ** 2 / (1 + alpha * (chi - 1) ** 3 * math.exp(beta * chi ** 1.75))


def tag(z):
    return f"z=0{z:6.4f}" if z < 10 else f"z={z:7.4f}"


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    build = os.path.abspath(sys.argv[1])
    out = os.path.abspath(sys.argv[2]) if len(sys.argv) > 2 else HERE
    radhydro = os.path.join(build, "radhydro")
    find_halos = os.path.join(build, "find_halos")

    ics = os.path.join(out, "ics")
    subprocess.check_call([sys.executable, os.path.join(ROOT, "tools/ics/make_ics.py"), ics,
                           "--kind", "zeldovich", "--Ngas", "32", "--Ndm", "16",
                           "--z-init", str(Z_INIT), "--box", str(BOX), "--seed", str(SEED)])

    for name, extra in (("onthefly", "source_mode = onthefly\n"),
                        ("legacy", "source_mode = legacy\ndir_halo = Halo\n")):
        d = os.path.join(out, name)
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, "params.txt"), "w") as f:
            f.write(COMMON + extra)

    # DM snapshots at the catalog redshifts (no RT before the last one)
    halo = os.path.join(out, "legacy", "Halo")
    shutil.rmtree(halo, ignore_errors=True)
    os.makedirs(halo)
    with tempfile.TemporaryDirectory() as tmp:
        prob = os.path.join(tmp, "p")
        os.makedirs(prob)
        with open(os.path.join(prob, "params.txt"), "w") as f:
            f.write(COMMON + "ics = %s\nzdata1 = 5\nzdata2 = %s\n"
                    % (ics, " ".join(str(z) for z in CATALOG_Z)))
        run = os.path.join(tmp, "run")
        subprocess.check_call([radhydro, prob, run, "--quiet"], stdout=subprocess.DEVNULL)
        snaps = sorted(os.listdir(os.path.join(run, "DM")), reverse=True)
        snaps = [s for s in snaps if s.startswith("xv_dm.z=")]
        # the first snapshot at or below each catalog redshift
        zs = [(float(s.split("=")[1]), s) for s in snaps]
        rows = []
        for zc in CATALOG_Z:
            cand = [(z, s) for z, s in zs if z <= zc]
            z, s = max(cand)
            res = subprocess.run([find_halos, os.path.join(run, "DM", s), repr(z),
                                  os.path.join(halo, "halo." + tag(zc)), "--box", str(BOX),
                                  "--ngas", "32", "--nmin", str(NMIN)],
                                 check=True, capture_output=True, text=True)
            nh = [l for l in res.stdout.splitlines() if "Halos" in l][0].strip()
            print(f"catalog {tag(zc)} from snapshot z = {z}: {nh}")
            rows.append(zc)
    with open(os.path.join(halo, "sfr.txt"), "w") as f:
        for zc in rows:
            # read(un,*) z,dummy(8),sfr: column 10 = comoving SFR density
            f.write(f"{zc:10.4f}" + " 0" * 8 + f" {hs_sfr(zc):.6e}\n")

    # smoke problem
    sm = os.path.join(out, "smoke")
    subprocess.check_call([sys.executable, os.path.join(ROOT, "tools/ics/make_ics.py"),
                           os.path.join(sm, "ics"), "--kind", "zeldovich", "--Ngas", "16",
                           "--Ndm", "8", "--z-init", str(Z_INIT), "--box", str(BOX),
                           "--seed", str(SEED)])
    with open(os.path.join(sm, "params.txt"), "w") as f:
        f.write(f"""# Small full-loop problem for the sanitizer / live Fortran tests
N = 16
Ndm = 8
box = {BOX}
nray_max_factor = 4
z_init = {Z_INIT}
zdata1 = 30 20
ics = ics
source_mode = legacy
dir_halo = Halo
fof_min_npart = {NMIN}
legacy_rt = 1
legacy_physics = 1
""")
    sh = os.path.join(sm, "Halo")
    shutil.rmtree(sh, ignore_errors=True)
    os.makedirs(sh)
    # source_type: x(3),v(3) real(4), N int(4), M,R,T,L,sfr,mstar,fgas real(4)
    rec = struct.pack("<6fi7f", 2.1, 1.9, 2.05, 0, 0, 0, 40, 1e10, 0.02, 3e4,
                      8.0e27 * 0.5, 0.5, 0, 0)
    for zc in (40.0, 20.0):
        with open(os.path.join(sh, "halo." + tag(zc)), "wb") as f:
            f.write(struct.pack("<q", 1) + rec)
    with open(os.path.join(sh, "sfr.txt"), "w") as f:
        for zc, sfr in ((40.0, 0.01), (20.0, 0.02)):
            f.write(f"{zc:10.4f}" + " 0" * 8 + f" {sfr:.6e}\n")
    print("wrote", out)


if __name__ == "__main__":
    main()
