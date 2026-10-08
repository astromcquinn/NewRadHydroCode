#!/usr/bin/env bash
# Run the PM problems with the Fortran harness.
#
#   run_all.sh [outroot] [debug]
#
# For each tests/problems_pm/<p> (or $PROBLEMS), reads Ngas, Ndm and mf from
# params.txt, builds the matching executable (make NGAS=.. NDM=.. MF=..) and
# runs it single-threaded into <outroot>/<p> (default reference/harness_pm/runs).
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "$here/../.." && pwd)"
outroot="${1:-$here/runs}"
mode="${2:-opt}"
problems="${PROBLEMS:-zeldovich_3d growth_3d pancake_x pancake_z_nomf}"

param() {  # param <file> <key> <default>
    local v
    v=$(sed -n "s/^[[:space:]]*$2[[:space:]]*=[[:space:]]*\([^[:space:]#]*\).*/\1/p" "$1" | head -1)
    echo "${v:-$3}"
}

export OMP_NUM_THREADS=1
for p in $problems; do
    prob="$root/tests/problems_pm/$p"
    ng=$(param "$prob/params.txt" Ngas 32)
    nd=$(param "$prob/params.txt" Ndm $((ng / 2)))
    mf=$(param "$prob/params.txt" mf 1)
    make -s -C "$here" MODE="$mode" NGAS="$ng" NDM="$nd" MF="$mf" >/dev/null
    exe="$here/$(make -s -C "$here" MODE="$mode" NGAS="$ng" NDM="$nd" MF="$mf" print-exe)"
    out="$outroot/$p"
    rm -rf "$out"
    mkdir -p "$out"
    start=$(date +%s.%N)
    status=0
    "$exe" "$prob" "$out" > "$out/stdout.log" 2>&1 || status=$?
    end=$(date +%s.%N)
    rm -f "$out/Init"
    printf '%s: exit %d, %.2f s (%s)\n' "$p" "$status" "$(echo "$end - $start" | bc)" "$exe"
    [ "$status" -eq 0 ] || exit "$status"
done
