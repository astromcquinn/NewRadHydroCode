#!/usr/bin/env python3
"""Convergence checks of hydro_scheme = unsplit (docs/hydro_unsplit.md §6b).

usage: convergence.py run|analyze|figures [--rundir DIR] [--driver PATH]
                      [--threads n] [--only sedov,box64,box10] [--results DIR]

  run      generate the problems and run every case not yet done; the cases
           run concurrently with their own OMP_NUM_THREADS, at most
           --threads cores in use (results do not depend on thread count)
  analyze  analyse every finished case into <results>/convergence_*.json
  figures  docs/figures/hydro_unsplit_convergence.png

Sets (problems and analysis exactly as in the suite, run_suite.py full tier):
  sedov  Sedov-Taylor blast (generate.gen_sedov, analyze.analyze_sedov) at
         N = 32, 64, 128, 256 for both schemes, plus the step and
         initial-blob variants at N = 128 (courant = 0.5, unsplit_cfl = 0.4,
         sigma = 2 cells)
  box64  the adiabatic cosmological box of the suite (64 Mpc/h, Ndm = N/2,
         z_init = 50) at N = 64 and 128 for a sequence of dlna_max (the
         cap of da/a per hydro step; "orig" = ra_max/(1-ra_max) = 1/9)
  box10  the same at the production box: 10 Mpc/h, N = Ndm = 128,
         z_init = 100 (tools/production/make_production.py), to z = 3
  box10n256  the production box at the production resolution N = Ndm = 256
         (not in the default --only)
"""
import argparse
import json
import os
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import hs_common as hc  # noqa: E402
import generate as gen  # noqa: E402
import analyze as an  # noqa: E402
import cosmo_gen as cg  # noqa: E402
import tests_extra as te  # noqa: E402

DLNA_ORIG = float(np.float32(0.1)) / (1 - float(np.float32(0.1)))


# ---------------------------------------------------------------- cases
def cases(only):
    """list of dicts: name, set, problem (name, generator, kwargs), sets, threads, cost"""
    out = []
    if "sedov" in only:
        thr = {32: 4, 64: 8, 128: 16, 256: 48}
        for scheme in ("legacy", "unsplit"):
            for N in (32, 64, 128, 256):
                out.append(dict(name=f"sedov_{scheme}_N{N}", set="sedov", scheme=scheme, N=N,
                                problem=(f"sedov_N{N}", gen.gen_sedov, dict(N=N)),
                                sets=[], threads=thr[N], cost=N ** 4.6))
        for scheme in ("legacy", "unsplit"):
            out.append(dict(name=f"sedov_{scheme}_N128_courant0.5", set="sedov", scheme=scheme,
                            N=128, problem=("sedov_N128", gen.gen_sedov, dict(N=128)),
                            sets=["courant=0.5"], threads=16, cost=2 * 128 ** 4.6))
        out.append(dict(name="sedov_unsplit_N128_cfl0.4", set="sedov", scheme="unsplit", N=128,
                        problem=("sedov_N128", gen.gen_sedov, dict(N=128)),
                        sets=["unsplit_cfl=0.4"], threads=16, cost=2 * 128 ** 4.6))
        out.append(dict(name="sedov_unsplit_N64_courant0.5", set="sedov", scheme="unsplit", N=64,
                        problem=("sedov_N64", gen.gen_sedov, dict(N=64)),
                        sets=["courant=0.5"], threads=8, cost=2 * 64 ** 4.6))
        for scheme in ("legacy", "unsplit"):
            out.append(dict(name=f"sedov_{scheme}_N128_sigma2", set="sedov", scheme=scheme,
                            N=128, problem=("sedov_N128_sigma2", gen.gen_sedov,
                                            dict(N=128, sigma=2.0)),
                            sets=[], threads=16, cost=128 ** 4.6))
    if "box64" in only:
        for N, dl in ((64, ("orig", 0.05, 0.02, 0.01, 0.005)),
                      (128, ("orig", 0.05, 0.03, 0.02, 0.01, 0.005, 0.0025))):
            prob = (f"box64_N{N}", cg.gen_cosmobox, dict(N=N, box=64.0))
            for d in dl:
                out.append(dict(name=f"box64_unsplit_N{N}_dlna{d}", set="box64",
                                scheme="unsplit", N=N, dlna=d, problem=prob,
                                sets=[] if d == "orig" else [f"dlna_max={d}"],
                                threads=4 if N == 64 else 12,
                                cost=N ** 3 * 200 / (min(DLNA_ORIG if d == "orig" else d, 0.03)
                                                     / 0.03)))
            for d in ("orig", 0.02):
                out.append(dict(name=f"box64_legacy_N{N}_dlna{d}", set="box64", scheme="legacy",
                                N=N, dlna=d, problem=prob,
                                sets=[] if d == "orig" else [f"dlna_max={d}"],
                                threads=4 if N == 64 else 12, cost=N ** 3 * 300))
    if "box10" in only:
        prob = ("box10_N128", cg.gen_cosmobox,
                dict(N=128, Ndm=128, box=10.0, z_init=100.0, z_out=(15, 10, 8, 6, 5, 4, 3)))
        for d in ("orig", 0.05, 0.02, 0.01, 0.005, 0.0025):
            out.append(dict(name=f"box10_unsplit_N128_dlna{d}", set="box10", scheme="unsplit",
                            N=128, dlna=d, problem=prob,
                            sets=[] if d == "orig" else [f"dlna_max={d}"], threads=12,
                            cost=128 ** 3 * 400 / (min(DLNA_ORIG if d == "orig" else d, 0.03)
                                                   / 0.03)))
        out.append(dict(name="box10_unsplit_N128_dlna0.005_courant0.5", set="box10",
                        scheme="unsplit", N=128, dlna=0.005, problem=prob,
                        sets=["dlna_max=0.005", "courant=0.5"], threads=12,
                        cost=128 ** 3 * 400 * 8))
        for d in ("orig", 0.02):
            out.append(dict(name=f"box10_legacy_N128_dlna{d}", set="box10", scheme="legacy",
                            N=128, dlna=d, problem=prob,
                            sets=[] if d == "orig" else [f"dlna_max={d}"], threads=12,
                            cost=128 ** 3 * 400))
    if "box10n256" in only:
        prob = ("box10_N256", cg.gen_cosmobox,
                dict(N=256, Ndm=256, box=10.0, z_init=100.0, z_out=(15, 10, 8, 6, 5, 4, 3)))
        for sc, d in (("unsplit", "orig"), ("unsplit", 0.02), ("unsplit", 0.005),
                      ("legacy", "orig")):
            out.append(dict(name=f"box10_{sc}_N256_dlna{d}", set="box10n256", scheme=sc,
                            N=256, dlna=d, problem=prob,
                            sets=[] if d == "orig" else [f"dlna_max={d}"], threads=24,
                            cost=256 ** 3 * 400 / (min(DLNA_ORIG if d == "orig" else d, 0.03)
                                                   / 0.03)))
    return out


# ---------------------------------------------------------------- run
class Cores:
    def __init__(self, n):
        self.free = n
        self.cv = threading.Condition()

    def take(self, k):
        with self.cv:
            while self.free < k:
                self.cv.wait()
            self.free -= k

    def give(self, k):
        with self.cv:
            self.free += k
            self.cv.notify_all()


def run(args):
    todo = cases(args.only)
    pdirs = {}
    for c in todo:   # generate each problem once (serially: make_ics is not cheap)
        pname, fun, kw = c["problem"]
        pdir = os.path.join(args.rundir, "problems", pname)
        if pname not in pdirs:
            if not os.path.exists(os.path.join(pdir, "meta.json")):
                print(f"generating {pname}", flush=True)
                fun(pdir, **kw)
            pdirs[pname] = pdir
    cores = Cores(args.threads)
    todo.sort(key=lambda c: -c["cost"])

    def one(c):
        odir = os.path.join(args.rundir, c["name"])
        if os.path.exists(os.path.join(odir, "done")):
            return
        k = min(c["threads"], args.threads)
        cores.take(k)
        try:
            sets = list(c["sets"]) + ([f"hydro_scheme={c['scheme']}"]
                                      if c["scheme"] != "legacy" else [])
            cmd = [args.driver, pdirs[c["problem"][0]], odir, "--quiet"]
            for s in sets:
                cmd += ["--set", s]
            env = dict(os.environ, OMP_NUM_THREADS=str(k))
            print(f"start {c['name']} ({k} threads): {' '.join(sets)}", flush=True)
            t0 = time.time()
            r = subprocess.run(cmd, capture_output=True, text=True, env=env)
            w = time.time() - t0
            if r.returncode != 0:
                print(f"FAILED {c['name']}: {r.stderr.strip()}", flush=True)
                return
            with open(os.path.join(odir, "done"), "w") as f:
                f.write(f"wall = {w:.1f}\nthreads = {k}\n")
            print(f"done {c['name']}: {w:.0f} s", flush=True)
        finally:
            cores.give(k)

    with ThreadPoolExecutor(max_workers=len(todo)) as ex:
        list(ex.map(one, todo))


# ---------------------------------------------------------------- analysis
def underdense_metrics(run, meta):
    """per output: median log10 T of the delta < 0.8 cells, and the median
    log10 of T over max(T_adiabatic(delta), 1 K) there (0 = on the adiabat
    or the floor), T_ad = T_init ((1+z)/(1+z_init))^2 delta^(2/3)"""
    N = meta["N"]
    res = []
    for d in hc.outputs(run):
        rho, vel, p, info = hc.read_state(d, N)
        z = info["redshift"]
        T = p / ((hc.GAMMA - 1) * rho) * info["temp1_unit"]
        delta = rho / cg.pc.FB
        sel = delta < 0.8
        if not sel.any():   # initial state of a small box
            sel = delta <= np.percentile(delta, 1)
        Tad = meta["T_init"] * ((1 + z) / (1 + meta["z_init"])) ** 2 * delta[sel] ** (2 / 3)
        lt = np.log10(np.maximum(T[sel], 1e-3))
        ex = lt - np.log10(np.maximum(Tad, 1.0))
        res.append({"z": z, "nsteps": int(info["istep"]), "f_under": float(sel.mean()),
                    "logT_under_med": float(np.median(lt)),
                    "logT_under_p84": float(np.percentile(lt, 84)),
                    "excess_med": float(np.median(ex)), "excess_p84": float(np.percentile(ex, 84)),
                    "excess_mean_lin": float(np.mean(10 ** ex) - 1),
                    "frac_under_floor": float(np.mean(T[sel] < 1.05))})
    return res


def analyze(args):
    res = {}
    for c in cases(args.only):
        odir = os.path.join(args.rundir, c["name"])
        if not os.path.exists(os.path.join(odir, "done")):
            print(f"missing {c['name']}")
            continue
        pdir = os.path.join(args.rundir, "problems", c["problem"][0])
        meta = hc.load_json(os.path.join(pdir, "meta.json"))
        info = dict(open(os.path.join(odir, "done")).read().split("\n")[i].split(" = ")
                    for i in range(2))
        if c["set"] == "sedov":
            r = an.analyze_sedov(odir, meta)
            for o in r["outputs"]:
                for k in ("profile", "R_dirs"):
                    o.pop(k, None)
            o = r["outputs"][-1]
            print(f"{c['name']:36s} steps {o['nsteps']:5d} R/Rex-1 "
                  + " ".join(f"{100 * q['R_err_rel']:+6.2f}%" for q in r["outputs"])
                  + f"  Rprof {100 * (o['R_profile'] / o['R_exact'] - 1):+6.2f}%"
                  f"  aniso {100 * o['aniso_rms']:.2f}% ax/dg {100 * o['axes_over_diag']:+.2f}%"
                  f"  peak {o['peak_rho_over_rho2']:.2f} L1 rho {o['l1']['rho']:.4f}"
                  f" dE {o['denergy']:+.2e}  wall {float(info['wall']):.0f}s", flush=True)
        else:
            r = te.analyze_cosmobox(odir, meta)
            for o in r["outputs"]:
                o.pop("trho_hist", None)
                o.pop("pk", None)
            for o, u in zip(r["outputs"], underdense_metrics(odir, meta)):
                o.update(u)
                # the metric of the §6 table: mean of the median log T of
                # the T-rho bins with log delta < -0.1
                t = o["trho"]
                ub = [y for x, y in zip(t["logd"], t["logT_med"]) if x < -0.1]
                o["logT_under_bins"] = float(np.mean(ub)) if ub else float("nan")
            for o in r["outputs"]:
                print(f"{c['name']:40s} z={o['z']:5.2f} steps {o['nsteps']:5d} "
                      f"T_med {o['T_med']:9.3g} T_mw {o['T_mw']:9.3g} logT(d<0.8) "
                      f"{o['logT_under_med']:5.2f} excess {o['excess_med']:5.2f} "
                      f"(p84 {o['excess_p84']:5.2f}) floor {o['frac_under_floor']:.2f} "
                      f"f(T>1e5) {o['frac_T_gt_1e5']:.3f}"
                      + (f" Pg/Pdm {o['pk_gas_over_dm_lowk'][0]:.3f}"
                         if "pk_gas_over_dm_lowk" in o else ""), flush=True)
        r.update({k: c[k] for k in ("set", "scheme", "N", "sets") if k in c})
        r["dlna"] = c.get("dlna")
        r["wall"] = float(info["wall"])
        r["threads"] = int(info["threads"])
        res[c["name"]] = r
    for s in args.only:
        sub = {k: v for k, v in res.items() if v["set"] == s}
        if sub:
            hc.save_json(os.path.join(args.results, f"convergence_{s}.json"), sub)


# ---------------------------------------------------------------- figure
def figures(args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ld = lambda s: hc.load_json(os.path.join(args.results, f"convergence_{s}.json"))  # noqa
    sed, b64, b10 = ld("sedov"), ld("box64"), ld("box10")
    try:
        b10.update(ld("box10n256"))
    except FileNotFoundError:
        pass
    C = {"legacy": "#b4542c", "unsplit": "#1f5fa8"}
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
    a = ax[0]
    Ns = np.array([32, 64, 128, 256])
    for sc in ("legacy", "unsplit"):
        keys = [f"sedov_{sc}_N{N}" for N in Ns if f"sedov_{sc}_N{N}" in sed]
        NN = np.array([sed[k]["N"] for k in keys])
        err = np.array([sed[k]["outputs"][-1]["R_err_rel"] for k in keys])
        dE = np.array([sed[k]["outputs"][-1]["denergy"] for k in keys])
        a.plot(NN, 100 * err, "o-", color=C[sc], label=f"{sc}")
        # the same with the radius rescaled to the energy actually left
        # (R ~ E^(1/5)): separates the position error from the energy leak
        a.plot(NN, 100 * ((1 + err) / (1 + dE) ** 0.2 - 1), "s:", color=C[sc], ms=4,
               label=f"{sc}, rescaled to E(t)")
    a.plot(Ns, 100 * 0.66 / (0.35 * Ns), "--", color="0.5", lw=0.8, label="0.66 cells")
    a.axhline(0, color="0.6", lw=0.8)
    a.axhline(1, color="0.6", lw=0.6, ls=":")
    a.set_xscale("log", base=2)
    a.set_xticks(Ns)
    a.set_xticklabels([str(n) for n in Ns])
    a.set_xlabel("N (R = 0.35 N at t = 10)")
    a.set_ylabel("R/R_exact - 1 at t = 10 [%]")
    a.set_title("Sedov shock radius (400-ray midpoint)")
    a.legend(fontsize=7)

    def dl(v):
        return DLNA_ORIG if v["dlna"] == "orig" else float(v["dlna"])

    cmap = plt.get_cmap("viridis")
    for a, data, title in ((ax[1], b64, "adiabatic box, 64 Mpc/h, N = 128 (dashed: 64)"),
                           (ax[2], b10, "production box, 10 Mpc/h, N = Ndm = 128 (dashed: 256)")):
        ks = [k for k, v in data.items() if v["scheme"] == "unsplit" and
              all(s.startswith("dlna") for s in v["sets"])]
        dls = sorted({dl(data[k]) for k in ks}, reverse=True)
        col = {d: cmap(i / max(len(dls) - 1, 1) * 0.9) for i, d in enumerate(dls)}
        for k in sorted(ks, key=lambda k: (-data[k]["N"], -dl(data[k]))):
            v = data[k]
            o = v["outputs"][1:]
            z = [q["z"] for q in o]
            y = [q["logT_under_bins"] for q in o]
            d = dl(v)
            lab = ("1/9 (orig)" if v["dlna"] == "orig" else f"{d:g}") if v["N"] == 128 else None
            a.plot([1 + q for q in z], y, "o-" if v["N"] == 128 else "s--", color=col[d],
                   ms=3, lw=1.2 if v["N"] == 128 else 0.8, label=lab)
            if v["N"] == 256 and v["dlna"] == "orig":
                a.plot([], [], "s--", color="k", ms=3, lw=0.8, label="N = Ndm = 256")
        for k, v in data.items():
            if v["scheme"] == "legacy" and v["dlna"] == "orig" and v["N"] == 128:
                o = v["outputs"][1:]
                a.plot([1 + q["z"] for q in o], [q["logT_under_bins"] for q in o], "-",
                       color=C["legacy"], lw=1, label="legacy, 1/9")
        a.set_xscale("log")
        ticks = [1, 2, 3, 4, 6, 11] if data is b64 else [4, 5, 6, 7, 9, 11, 16]
        a.set_xticks(ticks)
        a.set_xticklabels([str(t) for t in ticks])
        a.minorticks_off()
        a.invert_xaxis()
        a.set_xlabel("1 + z")
        a.set_ylabel("median log T [K] of the δ < 0.8 bins (mean)")
        a.set_title(title, fontsize=9)
        a.axhline(0, color="0.6", lw=0.8)
        a.legend(fontsize=7, title="max da/a", title_fontsize=7)
    fig.tight_layout()
    out = os.path.join(hc.ROOT, "docs", "figures", "hydro_unsplit_convergence.png")
    fig.savefig(out, dpi=110)
    print("wrote", out)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=("run", "analyze", "figures"))
    ap.add_argument("--rundir", default=os.path.join(hc.ROOT, "build", "hydro_conv"))
    ap.add_argument("--driver", default=None)
    ap.add_argument("--threads", type=int, default=os.cpu_count())
    ap.add_argument("--only", default="sedov,box64,box10")
    ap.add_argument("--results", default=os.path.join(HERE, "results"))
    args = ap.parse_args()
    args.only = args.only.split(",")
    args.rundir = os.path.abspath(args.rundir)
    if args.cmd == "run":
        args.driver = args.driver or hc.driver_path()
        run(args)
    elif args.cmd == "analyze":
        analyze(args)
    else:
        figures(args)


if __name__ == "__main__":
    main()
