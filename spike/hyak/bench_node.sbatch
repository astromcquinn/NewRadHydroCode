#!/bin/bash
# Single-node MPI x OpenMP layouts on one cpu-g2 node (192 Genoa cores).
#   source spike/hyak/env.sh && bash spike/hyak/build.sh   (once)
#   sbatch spike/hyak/bench_node.sbatch
# Results: $SPIKE_WORK/results/node_<jobid>.csv (+ the job log).
#SBATCH --job-name=spike-node
#SBATCH --account=astro
#SBATCH --partition=cpu-g2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=08:00:00
#SBATCH --output=spike-node-%j.log
set -u
source "${SPIKE_REPO:-$SLURM_SUBMIT_DIR}/spike/hyak/env.sh"
RES="$SPIKE_WORK/results"; mkdir -p "$RES"
OUT="$RES/node_${SLURM_JOB_ID}.csv"
export SER_BIN="$SPIKE_BUILD/serial" AMX_BIN="$SPIKE_BUILD/amrex" PAR_BIN="$SPIKE_BUILD/parthenon"
# OpenMPI under Slurm: let mpirun place ranks on contiguous core ranges
export MPIRUN="mpirun --oversubscribe"
lscpu | grep -E "Model name|Socket|NUMA node\(s\)|L3"

# Problems: 128^3 and 256^3 log-normal boxes, 512 sources, 4 hydro steps
for N in 128 256; do
  P="$SPIKE_WORK/problems/bench_N${N}"
  [ -f "$P/gas.bin" ] || python3 "$SPIKE_REPO/spike/bench/make_bench_problem.py" \
      --N $N --nsrc 512 --nsteps 4 --nray-max-factor 16 --out "$P"
done

# Layouts ranks:threads with ranks*threads = 192 (plus a few smaller
# counts for the thread-scaling curve). 8 NUMA domains (NPS4 x 2 sockets)
# and 24 CCDs: 8:24 and 24:8 are the NUMA- and CCD-aligned layouts.
LAY192="1:192 2:96 8:24 24:8 48:4 96:2 192:1"
LAYSMALL="1:12 1:48 8:6 24:2"
"$SPIKE_REPO"/spike/bench/bench.sh "$OUT" "$SPIKE_WORK/problems/bench_N128" "$LAYSMALL $LAY192" "32 16" 256 10
CODES="amrex parthenon serial" "$SPIKE_REPO"/spike/bench/bench.sh "$OUT" "$SPIKE_WORK/problems/bench_N256" "$LAY192" "32"
echo "results in $OUT"
