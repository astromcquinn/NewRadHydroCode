#!/bin/bash
# The I-fronts of the 8 sources of tests/problems/uniform_8src in uniform gas
# with the RT on refined levels (docs/stage_f.md section 3): uniform N,
# refined N (one level over BOX1, where fronts cross the level boundary and
# overlap there) and uniform 2N as the reference; then two levels (BOX1 and
# the nested BOX2) against uniform 4N. tools/analysis/ifront_multi.py
# compares x_HII,V/M(t), the ionised volume around each source and the
# photons in flight. One klone node.
#
#   sbatch --export=ALL,BUILD=$PWD/build-f3,OUT=$RUNS/rt_amr_ifront8[,N=32][,STEPS=40][,RANKS=8][,LEVELS="1 2"] \
#          tools/hyak/rt_amr_ifront8.slurm
# BOX1, BOX2: i0,j0,k0,i1,j1,k1 in base cells (multiples of 4 and 2 for the
# blocking factor 8; amrex_n_error_buf = 0, so a level is exactly its box).
# The defaults (N = 32): level 1 over 0..15 x 4..19 x 4..19 holds sources 1
# and 6 near its faces; their fronts overlap inside it and leave it, the
# fronts of 3, 8 and 2 (through the periodic face) enter it and overlap 6's
# and 1's at its faces. Level 2 over 2..13 x 6..17 x 6..17 holds 1 and 6.
#SBATCH --job-name=rt-amr-ifront8
#SBATCH --account=astro
#SBATCH --partition=cpu-g2
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=2:00:00
set -euo pipefail
if [ -n "${SLURM_SUBMIT_DIR:-}" ] && [ -f "$SLURM_SUBMIT_DIR/tools/hyak/env.sh" ]; then
  src="$SLURM_SUBMIT_DIR"
else
  src="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
fi
source "$src/tools/hyak/env.sh"
radray_module_load ompi/4.1.6-2
BUILD="${BUILD:-$src/build-f3}"
OUT="${OUT:-$RUNS/rt_amr_ifront8}"
N="${N:-32}"
STEPS="${STEPS:-40}"
RANKS="${RANKS:-8}"
LEVELS="${LEVELS:-1 2}"
BOX1="${BOX1:-0,4,4,15,19,19}"
BOX2="${BOX2:-2,6,6,13,17,17}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
mkdir -p "$OUT"
run() {  # run <name> <problem> [--set ...]
  local name="$1" prob="$2"; shift 2
  rm -rf "$OUT/$name"
  srun --mpi=pmix --ntasks="$RANKS" --cpus-per-task="$OMP_NUM_THREADS" --cpu-bind=cores \
    "$BUILD/radhydro" "$prob" "$OUT/$name" --quiet --set run_steps="$STEPS" "$@" > "$OUT/$name.out" 2>&1
  echo "$name: $(grep -o 'total [0-9.]* s' "$OUT/$name.out" | tail -1)"
}
prob() {  # the problem at n^3
  [ -d "$OUT/prob$1" ] || python3 "$src/tests/problems_full/make_ifront.py" "$OUT/prob$1" --N "$1" --sources 8src > /dev/null
}
amr=(--set amrex_refine_mass=1e9 --set amrex_n_error_buf=0)
prob "$N"
run "uniform$N" "$OUT/prob$N"
status=0
for L in $LEVELS; do
  n=$((N << L))
  box="$BOX1"
  [ "$L" -ge 2 ] && box="$BOX1,$BOX2"
  prob "$n"
  run "amr${N}_L$L" "$OUT/prob$N" "${amr[@]}" --set amrex_max_level="$L" --set amrex_refine_box="$box"
  [ -d "$OUT/uniform$n/step_0001" ] || run "uniform$n" "$OUT/prob$n"
  echo "== $L level(s), boxes $box, reference uniform $n^3"
  python3 "$src/tools/analysis/ifront_multi.py" --N "$N" --sources "$OUT/prob$N/sources_pos.txt" \
    "$OUT/uniform$N" "$OUT/amr${N}_L$L" "$OUT/uniform$n" \
    --check "$OUT/uniform$n" "$OUT/uniform$N" "$OUT/amr${N}_L$L" --json "$OUT/ifront8_L$L.json" || status=1
done
exit $status
