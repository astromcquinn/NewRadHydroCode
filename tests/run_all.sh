#!/usr/bin/env bash
# Run the C++ rt_driver on every problem in tests/problems and optionally
# compare with Fortran reference outputs (e.g. from
# reference/harness/run_all.sh <refroot>).
#
#   tests/run_all.sh <build_dir> <outroot> [<refroot>]
#
# With <refroot>: tools/compare/compare.py (fields + stats.txt) and
# tests/compare_reference.py (byte-identical fields, exact stats.txt and
# rt_stats.txt) for every problem.
set -uo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
build="$1"
out="$2"
ref="${3:-}"
problems="${PROBLEMS:-uniform_1src uniform_8src lognormal_16src uniform_8src_smallbuf lognormal_16src_smallbuf}"
mkdir -p "$out"
status=0
for name in $problems; do
  p="$here/problems/$name"
  t0=$(date +%s.%N)
  "$build/rt_driver" "$p" "$out/$name" > "$out/$name.stdout" || { echo "$name: rt_driver failed"; status=1; continue; }
  t1=$(date +%s.%N)
  printf '=== %-26s run time %6.2f s\n' "$name" "$(echo "$t1 - $t0" | bc)"
  if [ -n "$ref" ]; then
    python3 "$here/../tools/compare/compare.py" -q "$ref/$name" "$out/$name" | tail -6 || status=1
    python3 "$here/compare_reference.py" "$ref/$name" "$out/$name" || status=1
  fi
done
exit $status
