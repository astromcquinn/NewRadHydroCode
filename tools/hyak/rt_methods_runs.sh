#!/bin/bash
# Rays vs moment RT (M1, OTVET) on the production boxes (docs/rt_moments.md
# section 6). For each base problem directory (an existing ray run made by
# tools/production/make_production.py, e.g. reion256_c01), makes the siblings
# <base>_m1 and <base>_otvet: the same params.txt plus "rt_method = ..." and
# a symlink to the same ics/, so the three runs differ only in the RT. Then
# submits tools/hyak/reion_production.slurm for each new run (whole node).
# The base ray runs are not touched or resubmitted.
#
# usage (from the repository root, on a klone login node, after
# tools/hyak/build_hyak.sh has built build/radhydro from this branch):
#   tools/hyak/rt_methods_runs.sh [BASE ...]
# default BASE: $RUNS/reion256_c01 $RUNS/reion256 with
# RUNS=/gscratch/astro/$USER/runs.
# Environment: METHODS (default "m1 otvet"), DRYRUN=1 (make the directories
# and print the sbatch commands only), SBATCH_ARGS (default: a whole node,
# --exclusive -c 192 --mem=256G --time=2-00:00:00), THREADS (default 192).
set -eu
RUNS="${RUNS:-/gscratch/astro/$USER/runs}"
METHODS="${METHODS:-m1 otvet}"
THREADS="${THREADS:-192}"
SBATCH_ARGS="${SBATCH_ARGS:---exclusive --cpus-per-task=192 --mem=256G --time=2-00:00:00}"
[ $# -gt 0 ] || set -- "$RUNS/reion256_c01" "$RUNS/reion256"
[ -x build/radhydro ] || { echo "build/radhydro missing: run tools/hyak/build_hyak.sh first"; exit 1; }
for base in "$@"; do
  base="${base%/}"
  [ -f "$base/params.txt" ] || { echo "skip $base: no params.txt"; continue; }
  if grep -Eq '^[[:space:]]*rt_method' "$base/params.txt"; then
    echo "skip $base: already sets rt_method"; continue
  fi
  ics="$(cd "$base" && readlink -f ics)"
  for m in $METHODS; do
    dst="${base}_${m}"
    mkdir -p "$dst"
    if [ ! -f "$dst/params.txt" ]; then
      {
        echo "# $(basename "$base") with rt_method = $m (tools/hyak/rt_methods_runs.sh)"
        grep -v '^[[:space:]]*ics[[:space:]]*=' "$base/params.txt"
        echo "ics = $ics"
        echo "rt_method = $m"
      } > "$dst/params.txt"
    fi
    [ -e "$dst/ics" ] || ln -s "$ics" "$dst/ics"
    cmd=(sbatch $SBATCH_ARGS --job-name="reion_$m" --output="$dst/slurm-%j.out"
         --export=ALL,PROB="$dst",THREADS="$THREADS" tools/hyak/reion_production.slurm)
    if [ "${DRYRUN:-0}" = 1 ]; then echo "${cmd[*]}"; else "${cmd[@]}"; fi
  done
done
