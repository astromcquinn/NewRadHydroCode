#!/bin/bash
# The I-front of one source in uniform gas with the RT on a refined level
# (docs/stage_f.md section 3): uniform N, refined N (one level over the
# central half of the box, the front crossing its boundary) and uniform 2N,
# then tools/analysis/ifront_full.py --check. One klone node.
#
#   sbatch --export=ALL,BUILD=$PWD/build-sink,OUT=$RUNS/rt_amr_ifront[,N=32][,STEPS=40][,RANKS=8][,HALF=8] \
#          tools/hyak/rt_amr_ifront.slurm
# HALF: the half-width of the refined box around the centre, in base cells
# (8: the central half of a 32^3 box; 4: the front crosses the boundary early)
#SBATCH --job-name=rt-amr-ifront
#SBATCH --account=astro
#SBATCH --partition=cpu-g2
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=4:00:00
set -euo pipefail
if [ -n "${SLURM_SUBMIT_DIR:-}" ] && [ -f "$SLURM_SUBMIT_DIR/tools/hyak/env.sh" ]; then
  src="$SLURM_SUBMIT_DIR"
else
  src="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
fi
source "$src/tools/hyak/env.sh"
radray_module_load ompi/4.1.6-2
BUILD="${BUILD:-$src/build-sink}"
OUT="${OUT:-$RUNS/rt_amr_ifront}"
N="${N:-32}"
STEPS="${STEPS:-40}"
RANKS="${RANKS:-8}"
HALF="${HALF:-8}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
mkdir -p "$OUT"
python3 "$src/tests/problems_full/make_ifront.py" "$OUT/prob$N" --N "$N" > /dev/null
python3 "$src/tests/problems_full/make_ifront.py" "$OUT/prob$((2 * N))" --N "$((2 * N))" > /dev/null
q=$((N / 2 - HALF)); q3=$((N / 2 + HALF - 1))
run() {  # run <name> <problem> [--set ...]
  local name="$1" prob="$2"; shift 2
  rm -rf "$OUT/$name"
  srun --mpi=pmix --ntasks="$RANKS" --cpus-per-task="$OMP_NUM_THREADS" --cpu-bind=cores \
    "$BUILD/radhydro" "$prob" "$OUT/$name" --quiet --set run_steps="$STEPS" "$@" > "$OUT/$name.out" 2>&1
  echo "$name: $(grep -o 'total [0-9.]* s' "$OUT/$name.out" | tail -1)"
}
run uniform$N "$OUT/prob$N"
run amr$N "$OUT/prob$N" --set amrex_max_level=1 --set amrex_refine_mass=1e9 --set amrex_refine_box=$q,$q,$q,$q3,$q3,$q3
run uniform$((2 * N)) "$OUT/prob$((2 * N))"
echo "refined box: base cells $q..$q3 (radius $HALF from the centre)"
python3 "$src/tools/analysis/ifront_full.py" "$OUT/uniform$N" "$OUT/amr$N" --N "$N"
python3 "$src/tools/analysis/ifront_full.py" "$OUT/uniform$((2 * N))" --N "$((2 * N))"
python3 "$src/tools/analysis/ifront_full.py" --N "$N" --check "$OUT/uniform$((2 * N))" "$OUT/uniform$N" "$OUT/amr$N" 0.25
