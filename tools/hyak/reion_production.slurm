#!/bin/bash
# Production reionization run on one Hyak cpu-g2 node, with automatic restart.
#
# Set up once (an interactive job; the ICs take a few minutes at 256^3):
#   python3 tools/production/make_production.py /gscratch/astro/$USER/runs/reion128 --N 128
# Submit from the repository root:
#   sbatch --export=ALL,PROB=/gscratch/astro/$USER/runs/reion128 \
#          --cpus-per-task=32 --mem=64G tools/hyak/reion_production.slurm
# A whole node (docs/cluster.md for the thread count):
#   sbatch --exclusive --cpus-per-task=192 --mem=256G --export=ALL,PROB=... \
#          tools/hyak/reion_production.slurm
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
# build/radhydro), THREADS (default: the job's cpus), LAUNCH (a prefix such
# as "numactl --interleave=all"), EXTRA (extra radhydro arguments, e.g.
# "--set c_reduced=0.1"; they must be the same at every restart).
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
if [ -n "${HYAK_MODULES:-}" ]; then module load $HYAK_MODULES; fi
PROB="${PROB:?set PROB to the problem directory}"
RUN="${RUN:-$PROB/run}"
EXE="${EXE:-build/radhydro}"
export OMP_NUM_THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-32}}" OMP_PLACES=cores OMP_PROC_BIND=close
ulimit -s unlimited
if [ -e "$RUN/DONE" ]; then echo "run already finished"; exit 0; fi
CK="$RUN/restart/checkpoint.bin"
[ -e "$CK" ] || CK="$(ls "$RUN"/restart/checkpoint_*.bin 2>/dev/null | sort | tail -1)"
args=("$PROB" "$RUN" --quiet)
if [ -n "$CK" ] && [ -e "$CK" ]; then
  echo "$(date) restarting from $CK on $(hostname), $OMP_NUM_THREADS threads"
  args+=(--restart "$CK")
else
  echo "$(date) starting $PROB on $(hostname), $OMP_NUM_THREADS threads"
fi
# Free disk check: a checkpoint is rewritten atomically (needs 2x its size)
df -h "$(dirname "$RUN")" 2>/dev/null | tail -1 || true
/usr/bin/time -v ${LAUNCH:-} "$EXE" "${args[@]}" ${EXTRA:-} >> "$RUN.log" 2>&1
rc=$?
echo "$(date) radhydro exit $rc"
grep -E "Maximum resident|Elapsed \(wall" "$RUN.log" | tail -2
if [ $rc -eq 0 ]; then
  touch "$RUN/DONE"
  # Stills and the ionisation history (tools/movies/stills.py needs the
  # movie fields fHI, temp and rhog)
  if ls "$RUN"/movie/frame_*_rhog.bin >/dev/null 2>&1; then
    python3 tools/movies/stills.py "$RUN" "$RUN/stills" || true
  fi
fi
exit $rc
