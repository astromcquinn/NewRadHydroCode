#!/bin/bash
# hydro_backend = amrex (docs/stage_e_hydro.md) across 1, 2 and 4 cpu-g2
# nodes: hydro_driver on a random3d problem (unsplit scheme), 24 MPI ranks x
# 8 OpenMP threads per node (one rank per CCD), for several box sizes, against
# the serial OpenMP code on one node. This is the measurement of the
# multi-node hydro slowdown of docs/framework_decision.md section 3.3 (the
# spike's AMReX hydro was ~10x slower per cell on 2 and 4 nodes than on one
# node with 16^3 boxes). Outputs must be byte-identical across all layouts.
#
#   python3 tools/hyak/make_hydro_problem.py $RUNS/hydro256 --N 256   # once
#   sbatch --nodes=1 --export=ALL,PROB=$RUNS/hydro256 tools/hyak/amrex_hydro_scaling.slurm
#   sbatch --nodes=2 --export=ALL,PROB=$RUNS/hydro256,SERIAL=0 tools/hyak/amrex_hydro_scaling.slurm
#   sbatch --nodes=4 --export=ALL,PROB=$RUNS/hydro256,SERIAL=0 tools/hyak/amrex_hydro_scaling.slurm
#
# Environment:
#   PROB       problem directory (required; params.txt with hydro_scheme = unsplit)
#   OUT        output root (default $PROB/amrex_hydro_scaling/n<nodes>)
#   STEPS      hydro steps (default 2)
#   RANKS_PER_NODE, THREADS   (default 24 and 8)
#   MGS        box sizes to try (default "16 32 64")
#   EXE        build with the backend (default build-amrex/hydro_driver)
#   SERIAL     1 (default): also run the serial code on all cores of one node
#   HYAK_MODULES  modules (default "ompi/4.1.6-2")
#   RADRAY_UCX_EAGER  1 (default): UCX_RNDV_THRESH=8m UCX_ZCOPY_THRESH=inf
#              (RADRAY_UCX_RNDV_THRESH to change the threshold), without
#              which FillBoundary is 10x slower across nodes
#              (docs/stage_e_hydro.md section 5); 0: the module's defaults
#SBATCH --job-name=amrex-hydro-scaling
#SBATCH --account=astro
#SBATCH --partition=cpu-g2
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=8
#SBATCH --mem=512G
#SBATCH --exclusive
#SBATCH --time=2:00:00
#SBATCH --output=%x-%j.out
set -u
cd "${SLURM_SUBMIT_DIR}"
source tools/hyak/env.sh
radray_module_load ${HYAK_MODULES:-ompi/4.1.6-2}
radray_ucx_eager  # UCX eager below 8 MB across nodes (env.sh; RADRAY_UCX_EAGER=0 for the module defaults)
PROB="${PROB:?set PROB to the problem directory}"
NODES="${SLURM_JOB_NUM_NODES:-1}"
OUT="${OUT:-$PROB/amrex_hydro_scaling/n$NODES}"
STEPS="${STEPS:-2}"
RPN="${RANKS_PER_NODE:-24}"
THREADS="${THREADS:-8}"
MGS="${MGS:-16 32 64}"
EXE="${EXE:-$PWD/build-amrex/hydro_driver}"
SERIAL_EXE="${SERIAL_EXE:-$PWD/build/hydro_driver}"
RANKS=$((RPN * NODES))
mkdir -p "$OUT"
{
  echo "host $(hostname), job ${SLURM_JOB_ID:-}, nodes $NODES, $(date)"
  echo "EXE $EXE, ranks $RANKS ($RPN per node) x threads $THREADS, box sizes $MGS, steps $STEPS"
  lscpu | grep -E "Model name|^CPU\(s\)|Thread|Core|Socket|NUMA|L3"
} | tee "$OUT/node.txt"
common=(--quiet --set nsteps="$STEPS" --set hydro_scheme=unsplit)

export OMP_NUM_THREADS=$THREADS OMP_PLACES=cores OMP_PROC_BIND=close
for m in $MGS; do
  rm -rf "$OUT/mgs$m"
  srun --ntasks="$RANKS" --cpus-per-task="$THREADS" --cpu-bind=cores \
    "$EXE" "$PROB" "$OUT/mgs$m" "${common[@]}" --set hydro_backend=amrex --set amrex_hydro_max_grid_size="$m" \
    > "$OUT/mgs$m.log" 2>&1
  echo "mgs $m rc=$?"; grep -E "hydro_driver: .*steps|AMReX hydro timers" "$OUT/mgs$m.log"
done

if [ "${SERIAL:-1}" = 1 ]; then
  export OMP_NUM_THREADS=$((RPN * THREADS))
  rm -rf "$OUT/serial"
  srun --nodes=1 --ntasks=1 --cpus-per-task="$OMP_NUM_THREADS" --cpu-bind=none \
    "$SERIAL_EXE" "$PROB" "$OUT/serial" "${common[@]}" > "$OUT/serial.log" 2>&1
  echo "serial rc=$?"; grep -E "hydro_driver: .*steps" "$OUT/serial.log"
fi

# byte-identity across box sizes (and against the serial run when present)
first=""
for d in serial $(for m in $MGS; do echo mgs$m; done); do
  [ -d "$OUT/$d/step_0000" ] || continue
  if [ -z "$first" ]; then first=$d; continue; fi
  if diff -rq "$OUT/$first" "$OUT/$d" -x '*.log' -x domain_gas.txt >/dev/null; then
    echo "$d identical to $first"
  else
    echo "$d DIFFERS from $first"; diff -rq "$OUT/$first" "$OUT/$d" -x '*.log' | head -5
  fi
done | tee "$OUT/identical.txt"
