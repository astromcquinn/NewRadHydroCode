#!/bin/bash
# Benchmark matrix for the framework spike (local or inside a Slurm job).
#
#   SER_BIN=... AMX_BIN=... PAR_BIN=... spike/bench/bench.sh <out.csv> <rt_problem_dir> \
#       "<ranks:threads ...>" "<block sizes>" [hydro_n] [hydro_steps]
#
# For each layout ranks:threads and block size it runs rt_amrex and
# rt_parthenon (rt_serial for ranks = 1 only, any thread count), then the
# hydro kernels on an n^3 blast problem, and appends one CSV line per run:
#   kernel,code,ranks,threads,block,wall_s,work,ns_per_unit_core
# where work = ray-substeps (RT) or cell-updates (hydro) and
# ns_per_unit_core = wall*ranks*threads/work in ns.
# Environment: MPIRUN (default "mpirun"), MPIFLAGS (default binds ranks to
# contiguous core ranges with OpenMPI: --map-by slot:PE=<threads> --bind-to core),
# CODES (default "serial amrex parthenon"), REPEAT (default 1, best of).
set -u
OUT=$1; PROB=$2; LAYOUTS=$3; BLOCKS=$4; HN=${5:-0}; HS=${6:-10}
ROOT=$(cd "$(dirname "$0")/../.." && pwd)
MPIRUN=${MPIRUN:-mpirun}
CODES=${CODES:-serial amrex parthenon}
REPEAT=${REPEAT:-1}
TMP=$(mktemp -d)
[ -f "$OUT" ] || echo "kernel,code,ranks,threads,block,wall_s,work,ns_per_unit_core" > "$OUT"

mpiflags() {  # threads
  if [ -n "${MPIFLAGS+x}" ]; then echo "$MPIFLAGS"; else echo "--map-by slot:PE=$1 --bind-to core"; fi
}
best_of() {  # command... ; prints "wall work"
  local best="" work=0
  for r in $(seq "$REPEAT"); do
    "$@" > "$TMP/log" 2>&1 || { echo "FAILED: $*" >&2; tail -5 "$TMP/log" >&2; echo "nan 0"; return; }
    local w; w=$(grep -o "wall [0-9.]* s" "$TMP/log" | tail -1 | awk '{print $2}')
    work=$(grep -o "ray-substeps [0-9]*" "$TMP/log" | tail -1 | awk '{print $2}')
    if [ -z "$best" ] || awk "BEGIN{exit !($w < $best)}"; then best=$w; fi
  done
  echo "$best ${work:-0}"
}
emit() {  # kernel code ranks threads block wall work
  local ns
  ns=$(awk -v w="$6" -v r="$3" -v t="$4" -v n="$7" 'BEGIN{ if (n>0 && w==w) printf "%.1f", 1e9*w*r*t/n; else print "nan"}')
  echo "$1,$2,$3,$4,$5,$6,$7,$ns" >> "$OUT"
  echo "$1,$2,$3,$4,$5,$6,$7,$ns"
}

for lay in $LAYOUTS; do
  np=${lay%%:*}; nt=${lay##*:}
  export OMP_NUM_THREADS=$nt
  for code in $CODES; do
    case $code in
      serial)
        [ "$np" = 1 ] || continue
        read -r w n < <(best_of "$SER_BIN/rt_serial" "$PROB" "$TMP/o" --no-output)
        emit rt serial 1 "$nt" 0 "$w" "$n" ;;
      amrex)
        for b in $BLOCKS; do
          read -r w n < <(best_of $MPIRUN $(mpiflags "$nt") -np "$np" -x OMP_NUM_THREADS \
            "$AMX_BIN/rt_amrex" "$PROB" "$TMP/o" --no-output --max_grid_size "$b")
          emit rt amrex "$np" "$nt" "$b" "$w" "$n"
        done ;;
      parthenon)
        for b in $BLOCKS; do
          read -r w n < <(best_of $MPIRUN $(mpiflags "$nt") -np "$np" -x OMP_NUM_THREADS \
            "$PAR_BIN/rt_parthenon" -i "$ROOT/spike/parthenon/rt.in" "$PROB" "$TMP/o" --no-output --mb "$b")
          emit rt parthenon "$np" "$nt" "$b" "$w" "$n"
        done ;;
    esac
  done
  if [ "$HN" -gt 0 ]; then
    cells=$((HN * HN * HN * HS))
    for code in $CODES; do
      case $code in
        serial)
          [ "$np" = 1 ] || continue
          read -r w n < <(best_of "$SER_BIN/hydro_serial" --n "$HN" --steps "$HS")
          emit hydro serial 1 "$nt" 0 "$w" "$cells" ;;
        amrex)
          for b in $BLOCKS; do
            read -r w n < <(best_of $MPIRUN $(mpiflags "$nt") -np "$np" -x OMP_NUM_THREADS \
              "$AMX_BIN/hydro_amrex" --n "$HN" --steps "$HS" --max_grid_size "$b")
            emit hydro amrex "$np" "$nt" "$b" "$w" "$cells"
          done ;;
        parthenon)
          for b in $BLOCKS; do
            read -r w n < <(best_of $MPIRUN $(mpiflags "$nt") -np "$np" -x OMP_NUM_THREADS \
              "$PAR_BIN/hydro_parthenon" -i "$ROOT/spike/parthenon/hydro.in" --n "$HN" --steps "$HS" --mb "$b")
            emit hydro parthenon "$np" "$nt" "$b" "$w" "$cells"
          done ;;
      esac
    done
  fi
done
rm -rf "$TMP"
