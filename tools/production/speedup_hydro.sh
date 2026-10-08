#!/bin/bash
# Hydro/PM/summaries speedup: the first <nsteps> steps of a problem (before
# RT) with each thread count. speedup_hydro.sh <problem_dir> <outroot> <nsteps> [threads...]
set -u
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exe="$here/../../build/radhydro"
prob="$1"; root="$2"; n="$3"; shift 3
threads=("$@")
[ ${#threads[@]} -eq 0 ] && threads=(1 2 4)
for t in "${threads[@]}"; do
  o="$root/t$t"
  mkdir -p "$o"
  find "$o" -mindepth 1 -delete
  OMP_NUM_THREADS=$t "$exe" "$prob" "$o" --quiet --set max_steps=$n --set movie_dz=0 \
      --set restart_every=0 > "$o.log" 2>&1
  # sum of the per-step component times (steps 2..n: step 1 includes start-up)
  awk -v t=$t '!/^#/ && $1 >= 2 {h+=$6; p+=$7; f+=$8; s+=$10; k++}
       END {printf "threads %d: %d steps, per step: hydro %.3f s, pm %.3f s, halos %.3f s, summaries %.3f s\n", t, k, h/k, p/k, f/k, s/k}' "$o/timing.txt"
done
