#!/bin/bash
# Bitwise validation of the AMReX prototype against the serial reference.
#
#   spike/amrex/validate.sh <serial_build_dir> <amrex_build_dir> <work_dir> [quick|hydro|rt]
#
# hydro: cmp of hydro_amrex vs hydro_serial (n=64, 20 steps, blast and sod)
#        for ranks {1,2,4} x max_grid_size {32,16} x OMP threads {1,2}.
# RT:    tools/compare/compare.py --rtol 0 --atol-frac 0 of rt_amrex vs
#        rt_serial on uniform_8src and lognormal_16src_smallbuf (legacy_rt=0,
#        plus legacy_physics=0 on lognormal) for ranks {1,2,4,8} x
#        max_grid_size {16,8} x OMP threads {1,2} (and 1 rank at mgs 32).
# Prints one PASS/FAIL line per case and a summary; exit status 1 on any FAIL.
set -u
SER=$(realpath "$1"); AMX=$(realpath "$2"); W=$(realpath -m "$3"); MODE=${4:-}
QUICK=""; [ "$MODE" = quick ] && QUICK=1
ROOT=$(cd "$(dirname "$0")/../.." && pwd)
MPIRUN="mpirun --oversubscribe --allow-run-as-root --bind-to none"
ulimit -s unlimited 2>/dev/null || true
mkdir -p "$W"
npass=0; nfail=0
report() {  # name status
  if [ "$2" = 0 ]; then echo "PASS  $1"; npass=$((npass+1)); else echo "FAIL  $1"; nfail=$((nfail+1)); fi
}

# ---------------- hydro ----------------
[ "$MODE" = rt ] || for ic in blast sod; do
  OMP_NUM_THREADS=1 "$SER/hydro_serial" --n 64 --steps 20 --ic $ic --out "$W/hs_$ic.bin" > /dev/null
  for np in 1 2 4; do for mgs in 32 16; do for th in 1 2; do
    [ -n "$QUICK" ] && [ $th = 2 ] && [ $np != 1 ] && continue
    o="$W/ha_${ic}_${np}_${mgs}_${th}.bin"
    OMP_NUM_THREADS=$th $MPIRUN -np $np -x OMP_NUM_THREADS "$AMX/hydro_amrex" --n 64 --steps 20 --ic $ic \
      --max_grid_size $mgs --out "$o" > /dev/null 2>&1
    cmp -s "$W/hs_$ic.bin" "$o"; report "hydro ic=$ic ranks=$np mgs=$mgs threads=$th (cmp)" $?
  done; done; done
done

# ---------------- RT ----------------
rt_case() {  # tag problem sets...
  local tag=$1 prob=$2; shift 2
  local sets=()
  for kv in "$@"; do sets+=(--set "$kv"); done
  local ref="$W/rs_$tag"
  rm -rf "$ref"
  OMP_NUM_THREADS=1 "$SER/rt_serial" "$ROOT/tests/problems/$prob" "$ref" "${sets[@]}" > "$ref.log" 2>&1
  local cfgs="1:32:1"
  for np in 1 2 4 8; do for mgs in 16 8; do for th in 1 2; do cfgs="$cfgs $np:$mgs:$th"; done; done; done
  for c in $cfgs; do
    IFS=: read np mgs th <<< "$c"
    [ -n "$QUICK" ] && [ $th = 2 ] && [ $np != 1 ] && continue
    local out="$W/ra_${tag}_${np}_${mgs}_${th}"
    rm -rf "$out"
    OMP_NUM_THREADS=$th $MPIRUN -np $np -x OMP_NUM_THREADS "$AMX/rt_amrex" "$ROOT/tests/problems/$prob" "$out" \
      "${sets[@]}" --max_grid_size $mgs > "$out.log" 2>&1
    python3 "$ROOT/tools/compare/compare.py" "$ref" "$out" --rtol 0 --atol-frac 0 > "$out.cmp" 2>&1 \
      && grep -q "RESULT: PASS" "$out.cmp"
    report "rt $tag ranks=$np mgs=$mgs threads=$th (compare.py --rtol 0 --atol-frac 0)" $?
  done
}
if [ "$MODE" != hydro ]; then
rt_case u8 uniform_8src legacy_rt=0
rt_case ln16 lognormal_16src_smallbuf legacy_rt=0
rt_case ln16p0 lognormal_16src_smallbuf legacy_rt=0 legacy_physics=0
fi

echo "SUMMARY: $npass PASS, $nfail FAIL"
[ $nfail = 0 ]
