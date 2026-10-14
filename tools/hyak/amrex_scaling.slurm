#!/bin/bash
# rt_backend = amrex on one cpu-g2 node: 24 MPI ranks x 8 OpenMP threads
# (one rank per CCD, the best layout of the framework spike,
# docs/framework_decision.md section 3.2) at 256^3, restarted from a
# production checkpoint for STEPS hydro steps, against the serial OpenMP
# code on 192 threads (docs/cluster.md section 6). Outputs in $OUT/{serial,
# amrex}/timing.txt; the per-step RT times are the comparison
# (docs/stage_d.md section 4).
#
#   sbatch --exclusive --export=ALL,PROB=runs/reion256_f0.2,CK=<checkpoint> \
#          tools/hyak/amrex_scaling.slurm
#   # several nodes (24 ranks per node, no serial run):
#   sbatch --exclusive --nodes=2 --ntasks=48 --mem=1400G \
#          --export=ALL,PROB=...,CK=...,SERIAL=0 tools/hyak/amrex_scaling.slurm
#
# Environment:
#   PROB, CK   problem directory and checkpoint (required; the checkpoint
#              must be from a legacy_rt = 0 run without recomb_rad: the
#              backend needs it)
#   OUT        output root (default $PROB/amrex_scaling, with the runs in
#              $OUT/amrex_n<nodes> and $OUT/serial)
#   STEPS      hydro steps per run (default 1)
#   RANKS      MPI ranks (default: the job's --ntasks, 24 per node), THREADS
#              threads per rank (default 8)
#   MGS        amrex_max_grid_size (default 32)
#   EXTRA      more --set options for the backend run, e.g. the state on the
#              levels (docs/stage_e_state.md): "--set hydro_backend=amrex
#              --set pm_backend=amrex --set amrex_hydro_max_grid_size=0
#              --set amrex_state=1" (the serial run gets none of them)
#   EXE        build with the backend (default build-amrex/radhydro)
#   SERIAL     1: also run the serial code on all cores (default on one node)
#   HYAK_MODULES  modules (default "ompi/4.1.6-2")
#   RADRAY_UCX_EAGER  1 (default): UCX eager below 8 MB across nodes (env.sh,
#              RADRAY_UCX_RNDV_THRESH; docs/stage_e_hydro.md section 5,
#              docs/stage_d.md section 4.4); 0: the module's defaults
#
# Memory: every rank holds the full serial state (Stage D: about 15 GB at
# 256^3 and 40-50 GB at 512^3) and rank 0 reads the whole checkpoint, so at
# 512^3 ask for the node's memory (--mem=1400G) and run fewer, wider ranks:
#   sbatch --exclusive --nodes=2 --ntasks=24 --cpus-per-task=16 --mem=1400G \
#          --export=ALL,PROB=...,CK=...,THREADS=16 tools/hyak/amrex_scaling.slurm
# (24 ranks x 8 threads per node overflowed 1.5 TB at 512^3 in the restart.)
# An OOM-killed rank aborts the step (--kill-on-bad-exit); without it the
# other ranks wait in MPI until the time limit. A hung run spins until the
# time limit (give short runs a short --time); the SIGTERM of the limit or
# of scancel leaves Backtrace.<rank>.<thread> files in $OUT/amrex_n<nodes>
# with every rank's call stack (docs/stage_e_state.md section 5).
# Peak RSS per rank: the last line of timing.txt (peak_rss_MB of rank 0 and
# the min/max over the ranks) and sacct's MaxRSS of the job steps, printed
# at the end.
#SBATCH --job-name=amrex-scaling
#SBATCH --account=astro
#SBATCH --partition=cpu-g2
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=8
#SBATCH --mem=512G
#SBATCH --time=06:00:00
#SBATCH --exclusive
#SBATCH --output=%x-%j.out
set -u
cd "${SLURM_SUBMIT_DIR}"
source tools/hyak/env.sh
radray_module_load ${HYAK_MODULES:-ompi/4.1.6-2}
radray_ucx_eager  # UCX eager below 8 MB across nodes (env.sh; RADRAY_UCX_EAGER=0 for the module defaults)
PROB="$(readlink -f "${PROB:?set PROB to the problem directory}")"
CK="$(readlink -f "${CK:?set CK to the checkpoint}")"
OUT="${OUT:-$PROB/amrex_scaling}"
mkdir -p "$OUT" && OUT="$(readlink -f "$OUT")"
STEPS="${STEPS:-1}"
NODES="${SLURM_JOB_NUM_NODES:-1}"
RANKS="${RANKS:-${SLURM_NTASKS:-24}}"
THREADS="${THREADS:-8}"
MGS="${MGS:-32}"
EXE="$(readlink -f "${EXE:-$PWD/build-amrex/radhydro}")"
SERIAL="${SERIAL:-$([ "$NODES" = 1 ] && echo 1 || echo 0)}"
AMREX_OUT="$OUT/amrex_n$NODES"
[ -x "$EXE" ] || { echo "no executable $EXE (tools/hyak/build_amrex.sh)"; exit 1; }
[ -f "$CK" ] || { echo "no checkpoint $CK"; exit 1; }
{
  echo "host $(hostname), nodes ${SLURM_JOB_NODELIST:-}, job ${SLURM_JOB_ID:-}, $(date)"
  echo "EXE $EXE, $NODES nodes, ranks $RANKS x threads $THREADS, max_grid_size $MGS, checkpoint $CK, steps $STEPS${EXTRA:+, extra: $EXTRA}"
  lscpu | grep -E "Model name|^CPU\(s\)|Thread|Core|Socket|NUMA|L3"
} | tee "$OUT/node_n$NODES.txt"
common=(--quiet --restart "$CK" --set run_steps="$STEPS" --set restart_every=0 --set movie_dz=0 --set fields_every=0)

# the backend: one rank per CCD, threads pinned inside it
export OMP_NUM_THREADS=$THREADS OMP_PLACES=cores OMP_PROC_BIND=close
rm -rf "$AMREX_OUT"
mkdir -p "$AMREX_OUT"
# --chdir: a hung run cancelled or timed out leaves Backtrace.<rank>.<thread>
# (where every rank was) in $AMREX_OUT, not in the checkout
srun --mpi=pmix --ntasks="$RANKS" --cpus-per-task="$THREADS" --cpu-bind=cores --kill-on-bad-exit=1 --chdir="$AMREX_OUT" \
  "$EXE" "$PROB" "$AMREX_OUT" "${common[@]}" --set rt_backend=amrex --set deterministic_parallel=0 \
  --set amrex_max_grid_size="$MGS" ${EXTRA:-} > "$AMREX_OUT.log" 2>&1
echo "amrex rc=$?"; tail -3 "$AMREX_OUT.log"

# the serial OpenMP code on the whole node (deterministic_parallel)
if [ "$SERIAL" = 1 ]; then
  export OMP_NUM_THREADS=$((RANKS * THREADS))
  rm -rf "$OUT/serial"
  srun --ntasks=1 --cpus-per-task="$OMP_NUM_THREADS" --cpu-bind=none \
    "$EXE" "$PROB" "$OUT/serial" "${common[@]}" > "$OUT/serial.log" 2>&1
  echo "serial rc=$?"; tail -3 "$OUT/serial.log"
fi
echo "--- per-step timings (istep z nstep_rt Nsource nray t_hydro t_pm t_halos t_rt t_summ t_io t_total):"
for r in "$AMREX_OUT" "$OUT/serial"; do
  [ -f "$r/timing.txt" ] && { echo "$(basename "$r"):"; grep -v '^#' "$r/timing.txt"; grep '^# total' "$r/timing.txt"; }
done
echo "--- RT phase timers:"; grep -h 'RT timers' "$AMREX_OUT.log" "$OUT/serial.log" 2>/dev/null
echo "--- peak RSS of the job steps (sacct, MaxRSS = the largest task):"
sleep 5; sacct -j "$SLURM_JOB_ID" --format=JobID%18,JobName%12,NTasks,MaxRSS,MaxRSSTask,AveRSS,Elapsed -P 2>/dev/null
# the fields of the backend against the serial run (round-off, docs/stage_d.md)
if [ -f "$OUT/serial/timing.txt" ] && [ -f "$AMREX_OUT/timing.txt" ]; then
  echo "--- fields of the last step: backend vs serial (round-off expected, docs/stage_d.md):"
  python3 tools/hyak/amrex_scaling_compare.py "$OUT/serial" "$AMREX_OUT" --rtol 1e-3 | tee "$OUT/compare_n$NODES.txt"
fi
