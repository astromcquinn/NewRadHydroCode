#!/bin/bash
# What refinement costs at 512^3 (docs/stage_f.md section 3): a production
# checkpoint restarted for STEPS hydro steps on the AMReX backends with the
# state on the levels and load balancing, once per entry of LEVELS
# (amrex_max_level; default "0 2"), and the per-step times, the cells per
# level and the peak memory side by side. The refined levels are built at
# the first hydro call by interpolation from the base level, so the probe
# measures the cost of the hierarchy on a realistic density field and ray
# population, not a refined history (the fine gas has had no time to
# collapse further). Default checkpoint: the 512^3 case B run
# reion512_40mpc_f0.23 (40 Mpc/h, legacy_rt = legacy_physics = 0, unsplit,
# pm_refine = 1, docs/stage_e_state.md section 4). pm_refine = 2 adds its own
# known cost on top (docs/stage_e_pm.md section 4.1: PM ~4 s per step on one
# node, FoF and P(k) on rank 0 ~190 + 50 s per step at 1024^3 particles).
#
#   sbatch --nodes=5 --ntasks=120 \
#          --export=ALL,PROB=$RUNS/reion512_40mpc_f0.23 tools/hyak/amr_cost_probe.slurm
#
# Environment:
#   PROB        problem directory (default $RUNS/reion512_40mpc_f0.23)
#   CK          checkpoint (default $PROB/run/restart/checkpoint.bin); it must
#               be from a legacy_rt = 0 run without recomb_rad (the backend)
#   OUT         output root (default $PROB/amr_cost_probe, runs in $OUT/L<l>)
#   STEPS       hydro steps per run (default 2: the first carries the restart
#               and the first regrid)
#   LEVELS      amrex_max_level per run (default "0 2")
#   REFINE_MASS amrex_refine_mass (default 8: a level-l cell refines above
#               8 * 8^l times the mean gas density)
#   RANKS, THREADS  MPI ranks (default the job's --ntasks) x threads (8)
#   COUNT       1 (default): RADRAY_AMR_DEBUG=1 on the refined runs, which
#               prints the cells per level at every regrid (and scans the
#               levels for bad values: a small extra cost in the hydro); 0: off
#   EXTRA       more --set options for every run (e.g. "--set
#               cool_suppress_overdensity=0 --set sf_overdensity=3000")
#   EXE         build with AMReX (default build-amrex/radhydro, from main)
#   HYAK_MODULES, RADRAY_UCX_EAGER  as tools/hyak/amrex_scaling.slurm
# Memory: rank 0 reads the whole checkpoint (~293 GB at 512^3 with its rays,
# docs/stage_e_state.md section 4): keep --mem at the node's memory.
#SBATCH --job-name=amr-cost-probe
#SBATCH --account=astro
#SBATCH --partition=cpu-g2
#SBATCH --nodes=5
#SBATCH --ntasks=120
#SBATCH --cpus-per-task=8
#SBATCH --mem=1400G
#SBATCH --time=10:00:00
#SBATCH --exclusive
#SBATCH --output=%x-%j.out
set -u
cd "${SLURM_SUBMIT_DIR}"
source tools/hyak/env.sh
radray_module_load ${HYAK_MODULES:-ompi/4.1.6-2}
radray_ucx_eager
RUNS="${RUNS:-/gscratch/astro/$USER/radhydro/runs}"
PROB="$(readlink -f "${PROB:-$RUNS/reion512_40mpc_f0.23}")"
CK="$(readlink -f "${CK:-$PROB/run/restart/checkpoint.bin}")"
OUT="${OUT:-$PROB/amr_cost_probe}"
mkdir -p "$OUT" && OUT="$(readlink -f "$OUT")"
STEPS="${STEPS:-2}"
LEVELS="${LEVELS:-0 2}"
REFINE_MASS="${REFINE_MASS:-8}"
RANKS="${RANKS:-${SLURM_NTASKS:-120}}"
THREADS="${THREADS:-8}"
COUNT="${COUNT:-1}"
EXE="$(readlink -f "${EXE:-$PWD/build-amrex/radhydro}")"
[ -x "$EXE" ] || { echo "no executable $EXE (tools/hyak/build_amrex.sh)"; exit 1; }
[ -f "$CK" ] || { echo "no checkpoint $CK"; exit 1; }
echo "host $(hostname), nodes ${SLURM_JOB_NODELIST:-}, job ${SLURM_JOB_ID:-}, $(date)"
echo "EXE $EXE ($(git -C "$SLURM_SUBMIT_DIR" rev-parse --short HEAD 2>/dev/null)), ranks $RANKS x threads $THREADS, checkpoint $CK, steps $STEPS, levels $LEVELS, refine_mass $REFINE_MASS${EXTRA:+, extra: $EXTRA}"
export OMP_NUM_THREADS=$THREADS OMP_PLACES=cores OMP_PROC_BIND=close
ulimit -s unlimited
backend=(--set rt_backend=amrex --set hydro_backend=amrex --set pm_backend=amrex
         --set amrex_max_grid_size=32 --set amrex_hydro_max_grid_size=0 --set amrex_state=1
         --set amrex_balance=1 --set deterministic_parallel=0)
common=(--quiet --restart "$CK" --set run_steps="$STEPS" --set restart_every=0 --set movie_dz=0
        --set fields_every=0 --set step_timers=1)
for L in $LEVELS; do
  run="$OUT/L$L"
  rm -rf "$run" && mkdir -p "$run"
  dbg=()
  [ "$L" -gt 0 ] && [ "$COUNT" = 1 ] && dbg=(RADRAY_AMR_DEBUG=1)
  t0=$(date +%s)
  env "${dbg[@]}" srun --mpi=pmix --ntasks="$RANKS" --cpus-per-task="$THREADS" --cpu-bind=cores \
    --kill-on-bad-exit=1 --chdir="$run" \
    "$EXE" "$PROB" "$run" "${common[@]}" "${backend[@]}" \
    --set amrex_max_level="$L" --set amrex_refine_mass="$REFINE_MASS" ${EXTRA:-} > "$run.log" 2>&1
  echo "L$L: rc=$? in $(( $(date +%s) - t0 )) s"; tail -2 "$run.log"
done

echo "--- per-step timings (istep z nstep_rt Nsource nray t_hydro t_pm t_halos t_rt t_summ t_io t_total):"
for L in $LEVELS; do
  f="$OUT/L$L/timing.txt"
  [ -f "$f" ] && { echo "L$L:"; grep -v '^#' "$f"; }
done
echo "--- cells per level at the last regrid (RADRAY_AMR_DEBUG):"
for L in $LEVELS; do
  [ "$L" -gt 0 ] && grep 'AMR debug: regrid' "$OUT/L$L.log" | tail -1 | sed "s/^/L$L: /"
done
echo "--- timers (last of each kind):"
for L in $LEVELS; do
  for k in 'AMReX RT timers' 'AMReX hydro timers'; do
    grep -h "$k" "$OUT/L$L.log" 2>/dev/null | tail -1 | sed "s/^/L$L: /"
  done
done
# the volume above the refinement thresholds in the base-level gas of the
# unrefined run's last step: the cells the levels need before buffering
for L in $LEVELS; do
  [ "$L" = 0 ] || continue
  rho=$(ls -d "$OUT"/L0/step_*/rhog.bin 2>/dev/null | tail -1)
  [ -n "$rho" ] && python3 - "$rho" "$REFINE_MASS" <<'EOF'
import sys
import numpy as np
r = np.fromfile(sys.argv[1], dtype=np.float32).astype(np.float64)
m, n0 = float(sys.argv[2]), r.size
r /= r.mean()
f1, f2 = np.mean(r > m), np.mean(r > 8 * m)
print(f"--- base-level gas above {m:g} x mean: {f1:.4%} of the volume, above {8 * m:g}: {f2:.4%} "
      f"(the level-0 field: level 2 is underestimated)")
print(f"    unbuffered cells: level 1 {8 * f1 * n0:.3e} ({8 * f1:.3f} x base), "
      f"level 2 >= {64 * f2 * n0:.3e} ({64 * f2:.3f} x base)")
EOF
done
echo "--- ratio of the refined runs to L0 (sum over the steps):"
python3 - "$OUT" $LEVELS <<'EOF'
import sys
from pathlib import Path
import numpy as np
out, levels = Path(sys.argv[1]), sys.argv[2:]
t = {}
for L in levels:
    f = out / f"L{L}" / "timing.txt"
    if f.exists():
        a = np.loadtxt(f, ndmin=2)
        t[L] = a
if "0" in t:
    b = t["0"]
    for L, a in t.items():
        if L == "0":
            continue
        n = min(len(a), len(b))
        # columns: t_hydro 5, t_pm 6, t_halos 7, t_rt 8, t_total 11; the last
        # step only (the first carries the restart and the first regrid)
        for lab, c in (("hydro", 5), ("RT", 8), ("total", 11)):
            print(f"L{L}/L0 {lab}: last step {a[n - 1, c]:.1f} / {b[n - 1, c]:.1f} s = "
                  f"{a[n - 1, c] / max(b[n - 1, c], 1e-9):.2f}x; all steps {a[:n, c].sum() / max(b[:n, c].sum(), 1e-9):.2f}x")
EOF
echo "--- peak RSS of the job steps (sacct, MaxRSS = the largest task):"
sleep 5; sacct -j "$SLURM_JOB_ID" --format=JobID%18,JobName%12,NTasks,MaxRSS,MaxRSSTask,Elapsed -P 2>/dev/null
