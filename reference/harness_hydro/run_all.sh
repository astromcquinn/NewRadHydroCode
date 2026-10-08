#!/usr/bin/env bash
# Run the Fortran hydro harness on every problem in tests/problems_hydro.
#
#   run_all.sh [outroot=reference/harness_hydro/runs] [opt|debug]
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "$here/../.." && pwd)"
outroot="${1:-$here/runs}"
mode="${2:-opt}"
python3 "$here/make_hydro_problems.py" "$root/tests/problems_hydro" >/dev/null
for p in "$root"/tests/problems_hydro/*/; do
  name=$(basename "$p")
  N=$(awk -F= '$1 ~ /^ *N *$/ {gsub(/ /,"",$2); print $2}' "$p/params.txt")
  make -s -C "$here" MODE="$mode" N="$N" >/dev/null
  exe="$here/$(make -s -C "$here" MODE="$mode" N="$N" print-exe)"
  rm -rf "${outroot:?}/$name"
  printf '%-10s ' "$name"
  "$here/run_problem.sh" "$p" "$outroot/$name" "" "$exe"
done
