#!/usr/bin/env bash
# Run all validation problems with the harness.
#
#   run_all.sh [outroot] [debug]
#
# For each tests/problems/<p>, reads N and nray_max_factor from params.txt,
# builds the matching executable (make N=.. NRF=.. [MODE=debug]) and runs it
# into <outroot>/<p> (default reference/harness/runs/).  The D2D/ slice dumps
# written by write_data_3 are removed afterwards.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "$here/../.." && pwd)"
outroot="${1:-$here/runs}"
mode="${2:-opt}"
problems="${PROBLEMS:-uniform_1src uniform_8src lognormal_16src uniform_8src_smallbuf lognormal_16src_smallbuf}"

python3 "$here/make_test_problem.py"   # writes tests/problems/* only if missing

param() {  # param <file> <key> <default>
    local v
    v=$(sed -n "s/^[[:space:]]*$2[[:space:]]*=[[:space:]]*\([^[:space:]#]*\).*/\1/p" "$1" | head -1)
    echo "${v:-$3}"
}

for p in $problems; do
    prob="$root/tests/problems/$p"
    n=$(param "$prob/params.txt" N 32)
    nrf=$(param "$prob/params.txt" nray_max_factor 64)
    make -s -C "$here" MODE="$mode" N="$n" NRF="$nrf" >/dev/null
    exe="$here/$(make -s -C "$here" MODE="$mode" N="$n" NRF="$nrf" print-exe)"
    out="$outroot/$p"
    rm -rf "$out"
    start=$(date +%s.%N)
    "$here/run_problem.sh" "$prob" "$out" "" "$exe"
    end=$(date +%s.%N)
    rm -rf "$out/D2D"
    printf '%s: %.1f s (%s)\n' "$p" "$(echo "$end - $start" | bc)" "$exe"
done
du -sh "$outroot"
