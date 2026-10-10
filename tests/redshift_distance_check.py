#!/usr/bin/env python3
"""Smoke test of redshift_per_distance (docs/production_run.md).

Runs rt_driver on tests/problems/uniform_1src for a few steps with
da_sim = 5e-4 (about 61 RT cells of light per hydro step, so 6 substeps at
c_reduced = 0.1), with rt_method = rays and m1, and checks:
  * c_reduced = 0.1: redshift_per_distance = 1 runs and changes the results
    (the radiation is redshifted by fscale^s, s ~ 0.1, instead of fscale);
    the driver log and the ray_dump json record the flag;
  * c_reduced = 1: redshift_per_distance = 1 is byte-identical to 0 in every
    output file (the radiation moves at c: per time is per distance).

    redshift_distance_check.py RT_DRIVER PROBLEM_DIR OUTDIR
"""
import filecmp, json, os, shutil, subprocess, sys

NSTEPS = 3


def start(driver, prob, out, method, cr, flag, dump=False):
    if os.path.isdir(out):
        shutil.rmtree(out)
    os.makedirs(out)
    cmd = [driver, prob, out, "--nsteps", str(NSTEPS), "--set", "da_sim=0.0005",
           "--set", f"rt_method={method}", "--set", f"c_reduced={cr}",
           "--set", f"redshift_per_distance={flag}"]
    if dump:
        cmd += ["--set", "ray_dump_dir=" + os.path.join(out, "dump"), "--set", "ray_dump_at=1"]
    return subprocess.Popen(cmd, stdout=open(out + ".stdout", "w"), stderr=subprocess.STDOUT)


def step_files(out):
    fs = []
    for n in range(1, NSTEPS + 1):
        sd = os.path.join(out, f"step_{n:04d}")
        fs += [os.path.join(f"step_{n:04d}", f) for f in sorted(os.listdir(sd))]
    return fs


def same(a, b):
    fa, fb = step_files(a), step_files(b)
    if fa != fb or not fa:
        return False
    return all(filecmp.cmp(os.path.join(a, f), os.path.join(b, f), shallow=False) for f in fa)


def main():
    driver, prob, odir = sys.argv[1:4]
    runs = {}
    for m in ("rays", "m1"):
        for cr in ("0.1", "1"):
            for flag in (0, 1):
                tag = f"{m}_c{cr}_f{flag}"
                runs[tag] = start(driver, prob, os.path.join(odir, tag), m, cr, flag,
                                  dump=(m == "rays" and cr == "0.1" and flag == 1))
    ok = True
    for tag, p in runs.items():
        if p.wait() != 0:
            print(f"FAIL: rt_driver failed for {tag}")
            ok = False
    if not ok:
        return 1
    D = lambda t: os.path.join(odir, t)
    for m in ("rays", "m1"):
        if same(D(f"{m}_c0.1_f0"), D(f"{m}_c0.1_f1")):
            print(f"FAIL: {m}, c_reduced = 0.1: redshift_per_distance = 1 changed nothing"); ok = False
        if not same(D(f"{m}_c1_f0"), D(f"{m}_c1_f1")):
            print(f"FAIL: {m}, c_reduced = 1: redshift_per_distance = 1 is not byte-identical"); ok = False
        if " z_per_dist : 1" not in open(D(f"{m}_c0.1_f1") + ".stdout").read():
            print(f"FAIL: {m}: the driver log does not record the flag"); ok = False
    js = [f for f in os.listdir(os.path.join(D("rays_c0.1_f1"), "dump")) if f.endswith(".json")]
    if not js or json.load(open(os.path.join(D("rays_c0.1_f1"), "dump", js[0]))).get(
            "redshift_per_distance") != 1:
        print("FAIL: the ray_dump json does not record redshift_per_distance = 1"); ok = False
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
