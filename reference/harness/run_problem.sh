#!/usr/bin/env bash
# Run the Fortran RT harness on one test problem.
#
#   run_problem.sh <problem_dir> <outdir> [nsteps] [exe]
#
# ray_tracing_1 keeps a ~15 MB automatic array (merge1) on the stack, so the
# stack limit is raised first (the original ifort runs needed the same).
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
prob="$1"
out="$2"
nsteps="${3:-}"
exe="${4:-$here/build/rt_harness.x}"

ulimit -s unlimited 2>/dev/null || ulimit -s 65536
export OMP_NUM_THREADS=1
export OMP_STACKSIZE=64M

mkdir -p "$out"
status=0
"$exe" "$prob" "$out" $nsteps > "$out/stdout.log" 2>&1 || status=$?
# write_data_3 dumps 2D slices (data_2d_1/2) every RT substep into D2D/;
# they are not part of the validation outputs.
echo "exit status $status; stdout in $out/stdout.log"
exit $status
