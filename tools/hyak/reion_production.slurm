#!/bin/bash
# Production reionization run on one Hyak cpu-g2 node, with automatic restart.
#
# Set up once (login node or an interactive job):
#   python3 tools/production/make_production.py runs/reion$N --N $N \
#       --am-params <am_params.txt> --f-esc <calibrated> --xi-ion 25.38 \
#       --movie-dz 0.025 --movie-fields fHI,temp,rhog --restart-every 1
# Submit:
#   sbatch --export=ALL,PROB=runs/reion$N --cpus-per-task=32 --mem=64G \
#          tools/hyak/reion_production.slurm
#
# The job resumes from $PROB/run/restart/checkpoint.bin when it exists (a
# requeued or resubmitted job continues the run bit-identically; checkpoint
# writes are atomic, docs/production_run.md section 3). With --requeue it is
# safe on preemptible partitions (ckpt) too. When the run has reached the last
# zdata1 redshift it writes $PROB/run/DONE and later submissions exit at once.
#SBATCH --job-name=reion
#SBATCH --account=astro
#SBATCH --partition=cpu-g2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=7-00:00:00
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=%x-%j.out
set -u
cd "${SLURM_SUBMIT_DIR}"
PROB="${PROB:?set PROB to the problem directory}"
RUN="${RUN:-$PROB/run}"
EXE="${EXE:-build/radhydro}"
if [ -n "${HYAK_MODULES:-}" ]; then module load $HYAK_MODULES; fi
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-32}" OMP_PLACES=cores OMP_PROC_BIND=close
if [ -e "$RUN/DONE" ]; then echo "run already finished"; exit 0; fi
CK="$RUN/restart/checkpoint.bin"
args=("$PROB" "$RUN" --quiet)
if [ -e "$CK" ]; then
  echo "$(date) restarting from $CK"
  args+=(--restart "$CK")
else
  echo "$(date) starting $PROB"
fi
# Free disk check: a checkpoint is rewritten atomically (needs 2x its size)
df -h "$RUN" 2>/dev/null | tail -1 || true
"$EXE" "${args[@]}" >> "$PROB/log.txt" 2>&1
rc=$?
echo "$(date) radhydro exit $rc"
if [ $rc -eq 0 ]; then
  touch "$RUN/DONE"
  # Movie and reionization summary (cheap; uses the cores of this job)
  python3 tools/analysis/reion.py "$RUN" --out "$RUN/reion.txt" || true
  python3 tools/movie/make_movie.py "$RUN" --jobs "${OMP_NUM_THREADS}" || true
fi
exit $rc
