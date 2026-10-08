#!/bin/bash
# Speedup of one hydro step restarted from a checkpoint, for 1 and 4 threads.
#   speedup.sh <problem_dir> <checkpoint> <outroot> [threads...] [-- --set key=value ...]
# Writes <outroot>/t<N>/timing.txt (per-component times of the step, RT phase
# timers) and prints them; the outputs of the runs must be identical.
set -u
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exe="$here/../../build/radhydro"
prob="$1"; ck="$2"; root="$3"; shift 3
threads=()
while [ $# -gt 0 ] && [ "$1" != "--" ]; do threads+=("$1"); shift; done
[ $# -gt 0 ] && shift
[ ${#threads[@]} -eq 0 ] && threads=(1 4)
for t in "${threads[@]}"; do
  o="$root/t$t"
  mkdir -p "$o"
  find "$o" -mindepth 1 -delete
  OMP_NUM_THREADS=$t "$exe" "$prob" "$o" --quiet --restart "$ck" --set movie_dz=0 \
      --set restart_every=0 "$@" > "$o.log" 2>&1
  echo "threads $t (rc $?)"
  grep -v "^#" "$o/timing.txt" | tail -1
  grep "RT timers\|peak" "$o.log"
done
