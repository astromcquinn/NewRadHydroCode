#!/bin/bash
# OpenMP scaling of one hydro step (RT-dominated) restarted from a
# production checkpoint, on a whole cpu-g2 node. Outputs must be
# byte-identical across thread counts (deterministic_parallel); timings in
# $OUT/t*/timing.txt, summary in $OUT/scaling.txt.
#
#   sbatch --exclusive --export=ALL,PROB=runs/reion128,CK=<checkpoint> \
#          tools/hyak/scaling_test.slurm
#
# Environment:
#   PROB, CK   problem directory and checkpoint (required)
#   OUT        output root (default $PROB/scaling)
#   THREADS    thread counts (default "1 8 24 48 96 192"; the first is the
#              baseline of the speedups and of the byte-identity check)
#   LAUNCH     prefix for radhydro, e.g. "numactl --interleave=all"
#   EXE        executable (default build/radhydro)
#   STEPS      hydro steps per run (default 1)
#SBATCH --job-name=reion-scaling
#SBATCH --account=astro
#SBATCH --partition=cpu-g2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --mem=256G
#SBATCH --exclusive
#SBATCH --output=%x-%j.out
set -u
cd "${SLURM_SUBMIT_DIR}"
source tools/hyak/env.sh
if [ -n "${HYAK_MODULES:-}" ]; then radray_module_load $HYAK_MODULES; fi
export OMP_PLACES=cores OMP_PROC_BIND=close
PROB="${PROB:?set PROB to the problem directory}"
CK="${CK:?set CK to the checkpoint}"
OUT="${OUT:-$PROB/scaling}"
THREADS="${THREADS:-1 8 24 48 96 192}"
export EXE="${EXE:-$PWD/build/radhydro}" LAUNCH="${LAUNCH:-}"
mkdir -p "$OUT"
{
  echo "host $(hostname), job ${SLURM_JOB_ID:-}, $(date)"
  echo "EXE $EXE, LAUNCH '${LAUNCH}', OMP_PLACES=$OMP_PLACES OMP_PROC_BIND=$OMP_PROC_BIND"
  echo "checkpoint $CK"
  lscpu | grep -E "Model name|^CPU\(s\)|Thread|Core|Socket|NUMA|L3"
  numactl -H 2>/dev/null | head -4
  uptime
} | tee "$OUT/node.txt"
bash tools/production/speedup.sh "$PROB" "$CK" "$OUT" $THREADS | tee "$OUT/speedup.log"

# Identical outputs across thread counts (all files but the timings/logs):
[ -s "$OUT/t$(echo $THREADS | cut -d" " -f1)/timing.txt" ] || { echo "no timings: speedup.sh failed"; exit 1; }
base="$(echo $THREADS | awk '{print $1}')"
for t in $THREADS; do
  [ "$t" = "$base" ] && continue
  if diff -r "$OUT/t$base" "$OUT/t$t" -x timing.txt -x '*.log' >/dev/null; then
    echo "t$t identical to t$base"
  else
    echo "t$t DIFFERS from t$base"; diff -rq "$OUT/t$base" "$OUT/t$t" -x timing.txt -x '*.log' | head
  fi
done | tee "$OUT/identical.txt"
echo "files compared: $(find "$OUT/t$base" -type f ! -name timing.txt | wc -l)" | tee -a "$OUT/identical.txt"
python3 tools/hyak/scaling_table.py "$OUT" $THREADS | tee "$OUT/scaling.txt"
