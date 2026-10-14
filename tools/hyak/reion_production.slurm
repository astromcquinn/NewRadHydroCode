#!/bin/bash
# Production reionization run on one Hyak cpu-g2 node, with automatic restart.
#
# Set up once (an interactive job; the ICs take a few minutes at 256^3):
#   python3 tools/production/make_production.py /gscratch/astro/$USER/radhydro/runs/reion128 --N 128
# Submit from the repository root:
#   sbatch --export=ALL,PROB=/gscratch/astro/$USER/radhydro/runs/reion128 \
#          --cpus-per-task=32 --mem=64G tools/hyak/reion_production.slurm
# A whole node (docs/cluster.md for the thread count):
#   sbatch --exclusive --cpus-per-task=192 --mem=256G --export=ALL,PROB=... \
#          tools/hyak/reion_production.slurm
# The AMReX backends over several nodes (make_production.py --backend amrex,
# which writes rt_backend = amrex into params.txt; docs/stage_e_state.md):
# 24 ranks x 8 threads per node, the module and UCX settings of
# amrex_scaling.slurm, a build from tools/hyak/build_amrex.sh:
#   sbatch --exclusive --nodes=2 --ntasks=48 --cpus-per-task=8 --mem=1400G \
#          --export=ALL,PROB=...,EXE=build-amrex/radhydro,POST=sbatch \
#          tools/hyak/reion_production.slurm
# (RANKS and THREADS override the job's --ntasks and --cpus-per-task; at
# 512^3 ask for the node's memory, rank 0 holds the host mirror and reads
# the checkpoint, docs/stage_e_state.md section 4.)
#
# The job resumes from $RUN/restart/checkpoint.bin when it exists (a requeued
# or resubmitted job continues the run bit-identically; checkpoint writes are
# atomic, docs/production_run.md section 3). With restart_keep_latest = 0 it
# resumes from the newest checkpoint_NNNN.bin. With --requeue it is safe on
# preemptible partitions (ckpt) too. When the run has reached the last zdata1
# redshift it writes $RUN/DONE and later submissions exit at once.
#
# The partition's default time limit is 1 hour: keep --time. Jobs must also
# end before a maintenance reservation (scontrol show reservation), so the
# limit is 4 days; resubmit (or chain with --dependency=afterany:<job>) to
# continue from the last checkpoint.
#
# Environment: PROB (required), RUN (default $PROB/run), EXE (default
# build/radhydro, build-amrex/radhydro for the backends), THREADS (default:
# the job's cpus per task), RANKS (backends: MPI ranks, default the job's
# --ntasks), LAUNCH (a prefix such as "numactl --interleave=all"), EXTRA
# (extra radhydro arguments, e.g. "--set c_reduced=0.1"; they must be the
# same at every restart), HYAK_MODULES (backends: default ompi/4.1.6-2),
# RADRAY_UCX_EAGER (backends: 0 for the module's UCX defaults, env.sh).
# A hung backend run spins until the time limit; its SIGTERM (or scancel's)
# leaves Backtrace.<rank>.<thread> files in $RUN (amrex.handle_sigterm).
# POST: what to do after a successful run or leg (radhydro exit 0):
#   1 (default) run tools/hyak/reion_post.slurm (stills, reion.py,
#     make_movie.py) in this job; 0 nothing; sbatch submit
#     reion_post.slurm as a separate small job (POST_CPUS cores, default 16;
#     POST_MEM, default 64G; log POST_OUT, default $PROB/post-%j.out) and
#     release this node at once.
#SBATCH --job-name=reion
#SBATCH --account=astro
#SBATCH --partition=cpu-g2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=4-00:00:00
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=%x-%j.out
set -u
cd "${SLURM_SUBMIT_DIR}"
source tools/hyak/env.sh
if [ -n "${HYAK_MODULES:-}" ]; then radray_module_load $HYAK_MODULES; fi
PROB="$(readlink -f "${PROB:?set PROB to the problem directory}")"
RUN="${RUN:-$PROB/run}"
mkdir -p "$RUN" && RUN="$(readlink -f "$RUN")"
# the backends: params.txt says rt_backend = amrex (make_production.py
# --backend amrex); MPI ranks under srun, as tools/hyak/amrex_scaling.slurm
BACKEND=0
if grep -Eq '^[[:space:]]*rt_backend[[:space:]]*=[[:space:]]*amrex' "$PROB/params.txt"; then BACKEND=1; fi
if [ "$BACKEND" = 1 ]; then
  EXE="$(readlink -f "${EXE:-build-amrex/radhydro}")"
  radray_module_load ${HYAK_MODULES:-ompi/4.1.6-2}
  radray_ucx_eager
  RANKS="${RANKS:-${SLURM_NTASKS:-1}}"
  export OMP_NUM_THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-8}}"
else
  EXE="${EXE:-build/radhydro}"
  export OMP_NUM_THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-32}}"
fi
export OMP_PLACES=cores OMP_PROC_BIND=close
ulimit -s unlimited
if [ -e "$RUN/DONE" ]; then echo "run already finished"; exit 0; fi
CK="$RUN/restart/checkpoint.bin"
[ -e "$CK" ] || CK="$(ls "$RUN"/restart/checkpoint_*.bin 2>/dev/null | sort | tail -1)"
args=("$PROB" "$RUN" --quiet)
layout="$OMP_NUM_THREADS threads"
[ "$BACKEND" = 1 ] && layout="$RANKS ranks x $OMP_NUM_THREADS threads on ${SLURM_JOB_NUM_NODES:-1} nodes (${SLURM_JOB_NODELIST:-})"
if [ -n "$CK" ] && [ -e "$CK" ]; then
  echo "$(date) restarting from $CK on $(hostname), $layout"
  args+=(--restart "$CK")
else
  echo "$(date) starting $PROB on $(hostname), $layout"
fi
# Free disk check: a checkpoint is rewritten atomically (needs 2x its size)
df -h "$(dirname "$RUN")" 2>/dev/null | tail -1 || true
if [ "$BACKEND" = 1 ]; then
  # one rank per CCD, threads pinned inside it; a rank that dies (OOM) ends
  # the step instead of leaving the others in MPI; the run directory as the
  # working directory for the Backtrace files of a hung run
  /usr/bin/time -v srun --mpi=pmix --ntasks="$RANKS" --cpus-per-task="$OMP_NUM_THREADS" --cpu-bind=cores \
    --kill-on-bad-exit=1 --chdir="$RUN" ${LAUNCH:-} "$EXE" "${args[@]}" ${EXTRA:-} >> "$RUN.log" 2>&1
else
  /usr/bin/time -v ${LAUNCH:-} "$EXE" "${args[@]}" ${EXTRA:-} >> "$RUN.log" 2>&1
fi
rc=$?
echo "$(date) radhydro exit $rc"
grep -E "Maximum resident|Elapsed \(wall" "$RUN.log" | tail -2
if [ $rc -eq 0 ]; then
  touch "$RUN/DONE"
  # Stills, reionization history and tau_es, and the movie
  # (tools/hyak/reion_post.slurm; docs/production_run.md section 8,
  # docs/movies.md)
  case "${POST:-1}" in
    0) echo "POST=0: no post-processing" ;;
    sbatch)
      # a small separate job, so this node is released at once
      if jid=$(sbatch --parsable --job-name="post_${SLURM_JOB_NAME:-reion}" \
                      --cpus-per-task="${POST_CPUS:-16}" --mem="${POST_MEM:-64G}" \
                      --output="${POST_OUT:-$PROB/post-%j.out}" \
                      --export=ALL,RUN="$RUN",JOBS="${POST_CPUS:-16}" tools/hyak/reion_post.slurm); then
        echo "$(date) submitted post-processing job $jid"
      else
        echo "$(date) sbatch of the post-processing failed"
      fi
      ;;
    *) JOBS="${SLURM_CPUS_ON_NODE:-$OMP_NUM_THREADS}" RUN="$RUN" bash tools/hyak/reion_post.slurm ;;
  esac
fi
exit $rc
