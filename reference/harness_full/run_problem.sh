#!/usr/bin/env bash
# Run the full Fortran radhydro harness on one problem.
#
#   run_problem.sh <problem_dir> <outdir> [--debug] [--set key=value ...]
#
# Builds the executable for the problem's N, Ndm, nray_max_factor and box
# (params.txt plus --set overrides) if needed. The on-the-fly source mode
# links the C++ halo finder from RADRAY_BUILD (default <repo>/build), which
# must have been built first (cmake --build build --target radray_halos).
#
# ray_tracing_1 keeps a ~15 MB automatic array on the stack, so the stack
# limit is raised (as for the RT harness). The D2D/ slices written every RT
# substep (write_data_3) are deleted after the run.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "$here/../.." && pwd)"
prob="$(cd "$1" && pwd)"
out="$2"
shift 2
mode=opt
if [ "${1:-}" = "--debug" ]; then mode=debug; shift; fi

param() {  # param <key> <default>: params.txt value, overridden by --set
    local v
    v=$(sed -n "s/^[[:space:]]*$1[[:space:]]*=[[:space:]]*\([^[:space:]#]*\).*/\1/p" \
        "$prob/params.txt" | tail -1)
    local a
    for a in "$@"; do :; done
    local args=("${SETS[@]:-}")
    for a in "${args[@]}"; do
        case "$a" in "$1="*) v="${a#*=}";; esac
    done
    echo "${v:-$2}"
}
SETS=()
args=("$@")
for ((i = 0; i < ${#args[@]}; i++)); do
    if [ "${args[$i]}" = "--set" ]; then SETS+=("${args[$((i + 1))]}"); fi
done

N=$(param N 32)
NDM=$(param Ndm $((N / 2)))
NRF=$(param nray_max_factor 64)
BOX=$(param box 5)
RADRAY_BUILD="${RADRAY_BUILD:-$root/build}"
mk=(make -s -C "$here" MODE="$mode" N="$N" NDM="$NDM" NRF="$NRF" BOX="$BOX"
    RADRAY_BUILD="$RADRAY_BUILD")
"${mk[@]}" >/dev/null
exe="$here/$("${mk[@]}" print-exe)"

ulimit -s unlimited 2>/dev/null || ulimit -s 65536
export OMP_NUM_THREADS=1
export OMP_STACKSIZE=64M

rm -rf "$out"
mkdir -p "$out"
status=0
"$exe" "$prob" "$out" "$@" > "$out/stdout.log" 2>&1 || status=$?
rm -rf "$out/D2D"
echo "exit status $status; stdout in $out/stdout.log"
exit $status
