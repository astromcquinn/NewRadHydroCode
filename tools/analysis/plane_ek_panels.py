# E(k) small multiples by redshift with the compensated k^(5/3) E(k) row: plane_ek_panels.py EK_JSON OUT_PNG
# (EK_JSON from velocity_power.py --json; docs/plane_parallel.md, the C1 ladder)
import json, sys, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
d = json.load(open(sys.argv[1])); out = sys.argv[2]
cols = {"N64": "#86b6ef", "N128": "#3987e5", "N256": "#1c5cab", "N512": "#0d366b"}
zs = [6.6, 6.2, 6.0, 5.5, 5.0]
ink, muted, grid = "#2b2b2b", "#6b6b6b", "#e6e6e6"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": muted, "axes.labelcolor": ink,
                     "xtick.color": muted, "ytick.color": muted})
fig, ax = plt.subplots(2, len(zs), figsize=(13, 5.6), sharex=True, sharey="row",
                       gridspec_kw={"height_ratios": [3, 2]})
for j, z in enumerate(zs):
    a, b = ax[0, j], ax[1, j]
    for lab, c in cols.items():
        r = [r for r in d if r["label"] == lab and abs(r["z"] - z) < 0.01]
        if not r: continue
        r = r[0]; k = np.array(r["k"], float); E = np.array(r["E"], float); m = E > 0
        n = lab[1:] + "³"
        a.loglog(k[m], E[m], color=c, lw=1.6, label=n)
        b.semilogx(k[m], E[m] * k[m] ** (5 / 3), color=c, lw=1.6)
    kk = np.array([0.25, 1.0])
    a.loglog(kk, 25 * (kk / 0.25) ** (-5 / 3), color=muted, ls="--", lw=1.2)
    if j == 1: a.text(1.1, 2.0, "k$^{-5/3}$", color=muted, fontsize=8)
    for p in (a, b):
        p.axvspan(0.25, 1.0, color="#f3f3f3", zorder=0)
        p.grid(True, which="major", color=grid, lw=0.6); p.set_axisbelow(True)
        for s in ("top", "right"): p.spines[s].set_visible(False)
    a.set_title(f"z = {z}", color=ink)
    b.set_xlabel("k  [h kpc$^{-1}$]")
ax[0, 0].set_ylabel("E(k)  [(km/s)$^2$ kpc/h]")
ax[1, 0].set_ylabel("k$^{5/3}$ E(k)  (flat = Kolmogorov)")
ax[0, 0].set_ylim(1e-3, 3e3); ax[1, 0].set_ylim(0, 12)
ax[0, -1].legend(frameon=False, title="grid", fontsize=8, title_fontsize=8, loc="lower left")
fig.suptitle("Gas velocity power, C1 ladder from C26 ICs (Γ$_{-12}$ = 0.3, z$_{re}$ = 7, 0.25 Mpc/h); shaded: 0.25–1 h/kpc",
             color=ink, fontsize=10)
fig.tight_layout(); fig.savefig(out, dpi=150)
