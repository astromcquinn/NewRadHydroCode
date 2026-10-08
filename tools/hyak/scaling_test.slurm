#!/bin/bash
# OpenMP scaling of one hydro step (RT-dominated) restarted from a
# production checkpoint, on one cpu-g2 node. Outputs must be byte-identical
# across thread counts (deterministic_parallel); timings in $OUT/t*/timing.txt.
#   sbatch --export=ALL,PROB=runs/reion128,CK=runs/reion128/run/restart/checkpoint.bin \
#          tools/hyak/scaling_test.slurm
#SBATCH --job-name=reion-scaling
#SBATCH --account=astro
#SBATCH --partition=cpu-g2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --mem=200G
#SBATCH --time=04:00:00
#SBATCH --output=%x-%j.out
set -u
cd "${SLURM_SUBMIT_DIR}"
if [ -n "${HYAK_MODULES:-}" ]; then module load $HYAK_MODULES; fi
export OMP_PLACES=cores OMP_PROC_BIND=close
OUT="${OUT:-$PROB/scaling}"
tools/production/speedup.sh "$PROB" "$CK" "$OUT" 1 8 16 32 48 96 192
# Identical outputs across thread counts:
for t in 8 16 32 48 96 192; do
  if diff -rq "$OUT/t1" "$OUT/t$t" -x timing.txt -x '*.log' >/dev/null; then echo "t$t identical"; else echo "t$t DIFFERS"; fi
done
