#!/bin/bash
# UCX settings for the AMReX hydro exchange, swept on one fixed node set
# (docs/stage_e_hydro.md section 5): hydro_driver with hydro_backend = amrex
# on a random3d problem, 24 ranks x 8 threads per node, the settings
# interleaved and repeated so that the node pair, the placement and the
# fabric state are the same for every setting (single jobs on different
# node pairs scattered by up to 2x). Prints a table of fill / copyout /
# hydro per run, and checks that every run is byte-identical to the first.
# Two extra runs print UCX's protocol selection (UCX_PROTO_INFO=y) for the
# two thresholds being compared, and one run OpenMPI's pml_monitoring
# message statistics (bytes and messages per peer, rank 0).
#
#   sbatch --nodes=2 --export=ALL,PROB=$RUNS/hydro256 tools/hyak/amrex_hydro_ucx_sweep.slurm
#
# Environment:
#   PROB       problem directory (required; params.txt with hydro_scheme = unsplit)
#   OUT        output root (default $PROB/amrex_hydro_ucx_sweep/n<nodes>)
#   SETTINGS   "label=ENV;label=ENV;..." (ENV: space-separated VAR=value, empty
#              for the module defaults); default: the module defaults and
#              UCX_ZCOPY_THRESH=inf with UCX_RNDV_THRESH 8m, 16m, 32m, inf
#   REPEATS    rounds over the settings (default 3)
#   MGS        amrex_hydro_max_grid_size (default 32), STEPS (default 2)
#   RANKS_PER_NODE, THREADS   (default 24 and 8)
#   PROTO_INFO "label label" settings to run once more with UCX_PROTO_INFO=y
#              (default "8m 32m"; "" for none)
#   EXE        build with the backend (default build-amrex/hydro_driver)
#   HYAK_MODULES  modules (default "ompi/4.1.6-2")
#SBATCH --job-name=amrex-hydro-ucx-sweep
#SBATCH --account=astro
#SBATCH --partition=cpu-g2
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=8
#SBATCH --mem=512G
#SBATCH --exclusive
#SBATCH --time=1:00:00
#SBATCH --output=%x-%j.out
set -u
cd "${SLURM_SUBMIT_DIR}"
source tools/hyak/env.sh
radray_module_load ${HYAK_MODULES:-ompi/4.1.6-2}
# the sweep sets UCX itself: the module defaults are the baseline
PROB="${PROB:?set PROB to the problem directory}"
NODES="${SLURM_JOB_NUM_NODES:-1}"
OUT="${OUT:-$PROB/amrex_hydro_ucx_sweep/n$NODES}"
STEPS="${STEPS:-2}"
RPN="${RANKS_PER_NODE:-24}"
THREADS="${THREADS:-8}"
MGS="${MGS:-32}"
REPEATS="${REPEATS:-3}"
EXE="${EXE:-$PWD/build-amrex/hydro_driver}"
RANKS=$((RPN * NODES))
SETTINGS="${SETTINGS:-default=;8m=UCX_RNDV_THRESH=8m UCX_ZCOPY_THRESH=inf;16m=UCX_RNDV_THRESH=16m UCX_ZCOPY_THRESH=inf;32m=UCX_RNDV_THRESH=32m UCX_ZCOPY_THRESH=inf;inf=UCX_RNDV_THRESH=inf UCX_ZCOPY_THRESH=inf}"
PROTO_INFO="${PROTO_INFO-8m 32m}"
mkdir -p "$OUT"
{
  echo "host $(hostname), job ${SLURM_JOB_ID:-}, nodes $NODES (${SLURM_JOB_NODELIST:-}), $(date)"
  echo "EXE $EXE, ranks $RANKS ($RPN per node) x threads $THREADS, box size $MGS, steps $STEPS, repeats $REPEATS"
  echo "settings: $SETTINGS"
} | tee "$OUT/node.txt"
common=(--quiet --set nsteps="$STEPS" --set hydro_scheme=unsplit --set hydro_backend=amrex --set amrex_hydro_max_grid_size="$MGS")
export OMP_NUM_THREADS=$THREADS OMP_PLACES=cores OMP_PROC_BIND=close

env_of() {  # the ENV of a label
  local IFS=';' s
  for s in $SETTINGS; do
    [ "${s%%=*}" = "$1" ] && { echo "${s#*=}"; return; }
  done
}
labels() { local IFS=';' s; for s in $SETTINGS; do echo "${s%%=*}"; done; }

run_one() {  # label tag extra-env... -> one line of the table
  local label=$1 tag=$2; shift 2
  local dir="$OUT/run_$tag" log="$OUT/$tag.log"
  rm -rf "$dir"
  env $(env_of "$label") "$@" srun --ntasks="$RANKS" --cpus-per-task="$THREADS" --cpu-bind=cores \
    "$EXE" "$PROB" "$dir" "${common[@]}" > "$log" 2>&1
  local rc=$?
  local t; t=$(grep 'AMReX hydro timers' "$log" | tail -1)
  local fill copyout hydro
  fill=$(echo "$t" | sed -E 's/.* fill ([0-9.]+).*/\1/'); copyout=$(echo "$t" | sed -E 's/.* copyout ([0-9.]+).*/\1/')
  hydro=$(grep -o 'hydro [0-9.]* s' "$log" | tail -1 | cut -d' ' -f2)
  local same=ref
  if [ -d "$OUT/run_ref" ]; then
    if diff -rq "$OUT/run_ref" "$dir" -x '*.log' -x domain_gas.txt >/dev/null; then same=identical; else same=DIFFERS; fi
  else
    cp -r "$dir" "$OUT/run_ref"
  fi
  printf "%-10s %-12s rc=%s fill=%-7s copyout=%-7s hydro=%-7s %s\n" "$label" "$tag" "$rc" "$fill" "$copyout" "$hydro" "$same"
  [ "$same" = DIFFERS ] || rm -rf "$dir"
}

echo "--- sweep (label, run, fill / copyout / hydro [s], identity with the first run):"
for r in $(seq 1 "$REPEATS"); do
  for l in $(labels); do run_one "$l" "${l}_$r"; done
done | tee "$OUT/table.txt"

for l in $PROTO_INFO; do
  echo "--- UCX_PROTO_INFO=y with $l: $OUT/proto_$l.log"
  run_one "$l" "proto_$l" UCX_PROTO_INFO=y | tee -a "$OUT/table.txt"
done
echo "--- pml_monitoring (OpenMPI) with 8m: bytes and messages per peer of rank 0, $OUT/monitoring.*"
run_one 8m monitoring OMPI_MCA_pml_monitoring_enable=1 OMPI_MCA_pml_monitoring_enable_output=3 \
  OMPI_MCA_pml_monitoring_filename="$OUT/monitoring" | tee -a "$OUT/table.txt"
ls "$OUT"/monitoring* 2>/dev/null | head -3
f=$(ls "$OUT"/monitoring*0.prof 2>/dev/null | head -1)
[ -n "$f" ] && { echo "rank 0 ($f), first lines:"; head -40 "$f"; }
rm -rf "$OUT/run_ref"
