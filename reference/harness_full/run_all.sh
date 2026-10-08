#!/usr/bin/env bash
# Run the full Fortran harness on every validation case and (optionally)
# store the compact reference outputs.
#
#   run_all.sh [runs=reference/harness_full/runs] [--store] [--jobs n]
#
# Cases (tests/problems_full, see docs/full_run.md):
#   <problem>_rt1 : legacy_rt = 1, legacy_physics = 1 (the original)
#   <problem>_rt0 : legacy_rt = 0, legacy_physics = 0 (all fixes)
# for problem = smoke, legacy, onthefly (CASES overrides the list). The C++
# build directory (for the on-the-fly halo finder) is RADRAY_BUILD (default
# <repo>/build).
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "$here/../.." && pwd)"
runs="$here/runs"
store=0
jobs=4
while [ $# -gt 0 ]; do
    case "$1" in
        --store) store=1;;
        --jobs) jobs="$2"; shift;;
        *) runs="$1";;
    esac
    shift
done
cases="${CASES:-smoke_rt1 smoke_rt0 legacy_rt1 legacy_rt0 onthefly_rt1 onthefly_rt0}"
mkdir -p "$runs"

# build every executable first (sequentially; the runs below share them)
for c in $cases; do
    p="${c%_rt*}"
    "$here/run_problem.sh" "$root/tests/problems_full/$p" "$runs/.build_$p" \
        --set zdata1=1000 >/dev/null 2>&1 || true
    rm -rf "$runs/.build_$p"
done

run_case() {
    local c="$1" p="${1%_rt*}" f
    case "$c" in
        *_rt1) f="--set legacy_rt=1 --set legacy_physics=1";;
        *_rt0) f="--set legacy_rt=0 --set legacy_physics=0";;
    esac
    local t0
    t0=$(date +%s)
    # shellcheck disable=SC2086
    "$here/run_problem.sh" "$root/tests/problems_full/$p" "$runs/$c" $f >/dev/null
    echo "$c: $(( $(date +%s) - t0 )) s"
}
export -f run_case
export here root runs
printf '%s\n' $cases | xargs -P "$jobs" -I{} bash -c 'run_case {}'

if [ "$store" = 1 ]; then
    for c in $cases; do
        f=""
        case "$c" in legacy_rt1|smoke_*) f="--fields";; esac
        python3 "$here/store_outputs.py" "$runs/$c" "$here/outputs/$c" $f
    done
fi
