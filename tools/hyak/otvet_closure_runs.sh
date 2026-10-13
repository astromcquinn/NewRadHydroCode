#!/bin/bash
# OTVET with the alternative Eddington tensors (docs/otvet_closure.md) on an
# existing ray run's box: for each closure TAG makes the sibling
# <base>_otvet_<TAG> with the base's params.txt plus "rt_method = otvet" and
# the closure's parameters, and a symlink to the same ics/, so the runs
# differ from <base> and <base>_otvet only in the RT. Then submits
# tools/hyak/reion_production.slurm for each (whole node). The base run is
# not touched.
#
# usage (from the repository root, after tools/hyak/build_hyak.sh):
#   tools/hyak/otvet_closure_runs.sh BASE [TAG ...]
# TAGs (default: split screen cone):
#   split        otvet_split = 1
#   screen       otvet_closure = screen
#   screensplit  otvet_closure = screen, otvet_split = 1
#   cone         otvet_closure = cone (nside 4, coarse 4, r_c 8, budget,
#                blend 4: the defaults)
#   conenb       otvet_closure = cone, otvet_cone_blend = 0
# Environment: DRYRUN=1 (make the directories, print the sbatch commands),
# SBATCH_ARGS (default: a whole node, 2 days), THREADS (default 192),
# NO_TENSOR_DUMP=1 (append tensor_dump = 0, which overrides the base's: the
# dumps are 14 GB per snapshot at 512^3), PARAMS_EXTRA (more "key = value"
# lines, separated by ';'; not EXTRA, which reion_production.slurm passes
# to radhydro as arguments).
set -eu
THREADS="${THREADS:-192}"
SBATCH_ARGS="${SBATCH_ARGS:---exclusive --cpus-per-task=192 --mem=400G --time=2-00:00:00}"
[ $# -ge 1 ] || { echo "usage: $0 BASE [TAG ...]"; exit 2; }
base="${1%/}"
shift
[ $# -gt 0 ] || set -- split screen cone
[ -x build/radhydro ] || { echo "build/radhydro missing: run tools/hyak/build_hyak.sh first"; exit 1; }
[ -f "$base/params.txt" ] || { echo "$base: no params.txt"; exit 1; }
if grep -Eq '^[[:space:]]*rt_method' "$base/params.txt"; then
  echo "$base already sets rt_method: give the ray run"; exit 1
fi
ics="$(cd "$base" && readlink -f ics)"
for tag in "$@"; do
  case "$tag" in
    split) extra="otvet_split = 1" ;;
    screen) extra="otvet_closure = screen" ;;
    screensplit) extra="otvet_closure = screen;otvet_split = 1" ;;
    cone) extra="otvet_closure = cone" ;;
    conenb) extra="otvet_closure = cone;otvet_cone_blend = 0" ;;
    *) echo "unknown closure tag $tag"; exit 2 ;;
  esac
  [ "${NO_TENSOR_DUMP:-0}" = 1 ] && extra="$extra;tensor_dump = 0"
  [ -n "${PARAMS_EXTRA:-}" ] && extra="$extra;$PARAMS_EXTRA"
  dst="${base}_otvet_${tag}"
  mkdir -p "$dst"
  if [ ! -f "$dst/params.txt" ]; then
    {
      echo "# $(basename "$base") with rt_method = otvet, closure $tag (tools/hyak/otvet_closure_runs.sh)"
      grep -v '^[[:space:]]*ics[[:space:]]*=' "$base/params.txt" || true
      echo "ics = $ics"
      echo "rt_method = otvet"
      echo "$extra" | tr ';' '\n' | sed 's/^[[:space:]]*//'
    } > "$dst/params.txt"
  fi
  [ -e "$dst/ics" ] || ln -s "$ics" "$dst/ics"
  cmd=(sbatch $SBATCH_ARGS --job-name="otvcl_$tag" --output="$dst/slurm-%j.out"
       --export=ALL,EXTRA=,PROB="$dst",THREADS="$THREADS" tools/hyak/reion_production.slurm)
  if [ "${DRYRUN:-0}" = 1 ]; then echo "${cmd[*]}"; else "${cmd[@]}"; fi
done
