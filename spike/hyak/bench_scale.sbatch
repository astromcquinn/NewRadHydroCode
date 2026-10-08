#!/bin/bash
# Multi-node scaling (1, 2, 4 cpu-g2 nodes): strong scaling at 256^3 and
# weak scaling (N^3 per node fixed at ~192^3), MPI x OpenMP with the
# best single-node layout (set LAYOUT, default 24 ranks x 8 threads/node).
#   sbatch --nodes=4 spike/hyak/bench_scale.sbatch        (runs 1, 2, 4)
#SBATCH --job-name=spike-scale
#SBATCH --account=astro
#SBATCH --partition=cpu-g2
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=192
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=08:00:00
#SBATCH --output=spike-scale-%j.log
set -u
source "${SPIKE_REPO:-$SLURM_SUBMIT_DIR}/spike/hyak/env.sh"
RES="$SPIKE_WORK/results"; mkdir -p "$RES"
OUT="$RES/scale_${SLURM_JOB_ID}.csv"
export SER_BIN="$SPIKE_BUILD/serial" AMX_BIN="$SPIKE_BUILD/amrex" PAR_BIN="$SPIKE_BUILD/parthenon"
: "${LAYOUT:=24:8}"
RPN=${LAYOUT%%:*}; TPR=${LAYOUT##*:}
export CODES="amrex parthenon"
mk() { local N=$1 P="$SPIKE_WORK/problems/bench_N$1"
  [ -f "$P/gas.bin" ] || python3 "$SPIKE_REPO/spike/bench/make_bench_problem.py" \
      --N $N --nsrc $((N * N * N / 32768)) --nsteps 4 --nray-max-factor 16 --out "$P"
  echo "$P"; }
for nodes in 1 2 4; do
  [ "$nodes" -le "$SLURM_JOB_NUM_NODES" ] || continue
  export MPIRUN="mpirun --oversubscribe --map-by ppr:$RPN:node:PE=$TPR --bind-to core -H $(scontrol show hostnames | head -n $nodes | sed "s/$/:$RPN/" | paste -sd,)"
  export MPIFLAGS=""
  lay="$((nodes * RPN)):$TPR"
  # strong scaling, 256^3
  "$SPIKE_REPO"/spike/bench/bench.sh "$OUT" "$(mk 256)" "$lay" "32"
  # weak scaling: 192^3, 240^3, 304^3 cells for 1, 2, 4 nodes
  case $nodes in 1) N=192;; 2) N=240;; 4) N=304;; esac
  "$SPIKE_REPO"/spike/bench/bench.sh "$OUT" "$(mk $N)" "$lay" "16" $N 10
done
echo "results in $OUT"
