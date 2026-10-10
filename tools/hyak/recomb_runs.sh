#!/bin/bash
# Case A and case A + recombination radiation vs the original case B on the
# production boxes (docs/recomb_radiation.md section 4). For each base
# problem directory (an existing case-B ray run made by
# tools/production/make_production.py, e.g. reion256_c01), makes the
# siblings
#   <base>_caseA      recomb_case = A
#   <base>_caseA_rr   recomb_case = A, recomb_rad = 1 (+ RR_EXTRA lines)
# with the same params.txt otherwise and a symlink to the same ics/, and
# submits tools/hyak/reion_production.slurm for each (whole node). The base
# runs are not touched or resubmitted.
#
# usage (from the repository root, on a klone login node, after
# tools/hyak/build_hyak.sh has built build/radhydro from this branch):
#   tools/hyak/recomb_runs.sh [BASE ...]
# default BASE: $RUNS/reion256_c01 with RUNS=/gscratch/astro/$USER/radhydro/runs.
# Environment: VARIANTS (default "caseA caseA_rr"), RR_EXTRA (extra
# params.txt lines for caseA_rr, e.g. "recomb_rad_block = 4"), DRYRUN=1
# (make the directories and print the sbatch commands only), SBATCH_ARGS
# (default: a whole node for 2 days), THREADS (default 192).
set -eu
RUNS="${RUNS:-/gscratch/astro/$USER/radhydro/runs}"
VARIANTS="${VARIANTS:-caseA caseA_rr}"
THREADS="${THREADS:-192}"
SBATCH_ARGS="${SBATCH_ARGS:---exclusive --cpus-per-task=192 --mem=256G --time=2-00:00:00}"
[ $# -gt 0 ] || set -- "$RUNS/reion256_c01"
[ -x build/radhydro ] || { echo "build/radhydro missing: run tools/hyak/build_hyak.sh first"; exit 1; }
for base in "$@"; do
  base="${base%/}"
  [ -f "$base/params.txt" ] || { echo "skip $base: no params.txt"; continue; }
  if grep -Eq '^[[:space:]]*(recomb_|rt_method)' "$base/params.txt"; then
    echo "skip $base: already sets recomb_* or rt_method"; continue
  fi
  ics="$(cd "$base" && readlink -f ics)"
  for v in $VARIANTS; do
    dst="${base}_${v}"
    mkdir -p "$dst"
    if [ ! -f "$dst/params.txt" ]; then
      {
        echo "# $(basename "$base") with $v (tools/hyak/recomb_runs.sh)"
        grep -v '^[[:space:]]*ics[[:space:]]*=' "$base/params.txt"
        echo "ics = $ics"
        case "$v" in
          caseA) echo "recomb_case = A" ;;
          caseA_rr)
            echo "recomb_case = A"
            echo "recomb_rad = 1"
            [ -n "${RR_EXTRA:-}" ] && printf '%s\n' "$RR_EXTRA"
            ;;
          *) echo "unknown variant $v" >&2; exit 1 ;;
        esac
      } > "$dst/params.txt"
    fi
    [ -e "$dst/ics" ] || ln -s "$ics" "$dst/ics"
    cmd=(sbatch $SBATCH_ARGS --job-name="reion_$v" --output="$dst/slurm-%j.out"
         --export=ALL,PROB="$dst",THREADS="$THREADS" tools/hyak/reion_production.slurm)
    if [ "${DRYRUN:-0}" = 1 ]; then echo "${cmd[*]}"; else "${cmd[@]}"; fi
  done
done
