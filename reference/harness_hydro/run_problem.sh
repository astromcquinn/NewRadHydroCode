#!/usr/bin/env bash
# Run the Fortran hydro harness on one test problem.
#
#   run_problem.sh <problem_dir> <outdir> [nsteps] [exe]
#
# Builds build/N<N>/hydro_harness.x for the problem's N if no exe is given.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
prob="$(cd "$1" && pwd)"
out="$2"
nsteps="${3:-}"
N=$(awk -F= '$1 ~ /^ *N *$/ {gsub(/ /,"",$2); print $2}' "$prob/params.txt")
exe="${4:-}"
if [ -z "$exe" ]; then
  make -s -C "$here" N="$N" >/dev/null
  exe="$here/$(make -s -C "$here" N="$N" print-exe)"
fi
export OMP_NUM_THREADS=1
mkdir -p "$out"
status=0
"$exe" "$prob" "$out" $nsteps > "$out/stdout.log" 2>&1 || status=$?
echo "exit status $status; stdout in $out/stdout.log"
exit $status
