#!/usr/bin/env bash
# Run pm_driver on every PM problem and compare with Fortran harness runs.
#
#   tests/run_pm_all.sh <pm_driver> <fortran_runs_root> <outroot>
#
# e.g.  reference/fortran/harness_pm/run_all.sh            # Fortran -> reference/fortran/harness_pm/runs
#       tests/run_pm_all.sh build/pm_driver reference/fortran/harness_pm/runs build/pm_runs
# Every .bin and Pk file must be byte-identical and every stats value equal.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
drv="$1"; ref="$2"; out="$3"
problems="${PROBLEMS:-zeldovich_3d growth_3d pancake_x pancake_z_nomf}"
status=0
for p in $problems; do
    rm -rf "$out/$p"
    start=$(date +%s.%N)
    "$drv" "$here/problems_pm/$p" "$out/$p" --quiet > "$out.$p.log"
    end=$(date +%s.%N)
    printf '%-16s C++ %.2f s  ' "$p" "$(echo "$end - $start" | bc)"
    python3 "$here/compare_reference.py" "$ref/$p" "$out/$p" | tail -1 || status=1
    if [ -d "$ref/$p/Pk" ]; then
        diff -rq "$ref/$p/Pk" "$out/$p/Pk" || status=1
    fi
done
exit $status
