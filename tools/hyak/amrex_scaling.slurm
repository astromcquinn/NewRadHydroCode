#!/bin/bash
# rt_backend = amrex on one cpu-g2 node: 24 MPI ranks x 8 OpenMP threads
# (one rank per CCD, the best layout of the framework spike,
# docs/framework_decision.md section 3.2) at 256^3, restarted from a
# production checkpoint for STEPS hydro steps, against the serial OpenMP
# code on 192 threads (docs/cluster.md section 6). Outputs in $OUT/{serial,
# amrex}/timing.txt; the per-step RT times are the comparison
# (docs/stage_d.md section 4).
#
#   sbatch --exclusive --export=ALL,PROB=runs/reion256_c01,CK=<checkpoint> \
#          tools/hyak/amrex_scaling.slurm
#
# Environment:
#   PROB, CK   problem directory and checkpoint (required; the checkpoint
#              must be from a legacy_rt = 0 run: the backend needs it)
#   OUT        output root (default $PROB/amrex_scaling)
#   STEPS      hydro steps per run (default 1)
#   RANKS      MPI ranks (default 24), THREADS threads per rank (default 8)
#   MGS        amrex_max_grid_size (default 32)
#   EXE        build with the backend (default build-amrex/radhydro)
#   SERIAL     1 (default): also run the serial code on all cores
#   HYAK_MODULES  modules (default "ompi/4.1.6-2")
#SBATCH --job-name=amrex-scaling
#SBATCH --account=astro
#SBATCH --partition=cpu-g2
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=8
#SBATCH --mem=512G
#SBATCH --exclusive
#SBATCH --output=%x-%j.out
set -u
cd "${SLURM_SUBMIT_DIR}"
source tools/hyak/env.sh
module load ${HYAK_MODULES:-ompi/4.1.6-2}
PROB="${PROB:?set PROB to the problem directory}"
CK="${CK:?set CK to the checkpoint}"
OUT="${OUT:-$PROB/amrex_scaling}"
STEPS="${STEPS:-1}"
RANKS="${RANKS:-24}"
THREADS="${THREADS:-8}"
MGS="${MGS:-32}"
EXE="${EXE:-$PWD/build-amrex/radhydro}"
mkdir -p "$OUT"
{
  echo "host $(hostname), job ${SLURM_JOB_ID:-}, $(date)"
  echo "EXE $EXE, ranks $RANKS x threads $THREADS, max_grid_size $MGS, checkpoint $CK, steps $STEPS"
  lscpu | grep -E "Model name|^CPU\(s\)|Thread|Core|Socket|NUMA|L3"
} | tee "$OUT/node.txt"
common=(--quiet --restart "$CK" --set run_steps="$STEPS" --set restart_every=0 --set movie_dz=0 --set fields_every=0)

# the backend: one rank per CCD, threads pinned inside it
export OMP_NUM_THREADS=$THREADS OMP_PLACES=cores OMP_PROC_BIND=close
rm -rf "$OUT/amrex"
srun --ntasks="$RANKS" --cpus-per-task="$THREADS" --cpu-bind=cores \
  "$EXE" "$PROB" "$OUT/amrex" "${common[@]}" --set rt_backend=amrex --set deterministic_parallel=0 \
  --set amrex_max_grid_size="$MGS" > "$OUT/amrex.log" 2>&1
echo "amrex rc=$?"; tail -3 "$OUT/amrex.log"

# the serial OpenMP code on the whole node (deterministic_parallel)
if [ "${SERIAL:-1}" = 1 ]; then
  export OMP_NUM_THREADS=$((RANKS * THREADS))
  rm -rf "$OUT/serial"
  srun --ntasks=1 --cpus-per-task="$OMP_NUM_THREADS" --cpu-bind=none \
    "$EXE" "$PROB" "$OUT/serial" "${common[@]}" > "$OUT/serial.log" 2>&1
  echo "serial rc=$?"; tail -3 "$OUT/serial.log"
fi
echo "--- per-step timings (istep z nstep_rt Nsource nray t_hydro t_pm t_halos t_rt t_summ t_io t_total):"
for r in amrex serial; do [ -f "$OUT/$r/timing.txt" ] && { echo "$r:"; grep -v '^#' "$OUT/$r/timing.txt"; }; done
echo "--- RT phase timers:"; grep -h 'RT timers' "$OUT"/*.log
