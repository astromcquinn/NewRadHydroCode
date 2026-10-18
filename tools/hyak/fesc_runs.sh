#!/bin/bash
# Escape fraction falling with redshift (docs/fesc_evolution.md) on an
# existing ray run's box: for each F6:BETA makes the sibling
# <DEST>/<PREFIX>_b<BETA> with the base's params.txt, its f_esc replaced by
#   f_esc = F6, f_esc_beta = BETA, f_esc_z0 = 6
# (f_esc(z) = min(f_esc_max, F6 ((1+z)/7)^BETA), src/halo_sources.cpp
# f_esc_at), tensor_dump dropped, and a symlink to the same ics/, so the
# runs differ from the base only in the sources. Then submits
# tools/hyak/reion_production.slurm for each (whole node). The base run is
# not touched.
#
# usage (from the repository root, after tools/hyak/build_hyak.sh):
#   tools/hyak/fesc_runs.sh BASE PREFIX F6:BETA [F6:BETA ...]
# e.g. tools/hyak/fesc_runs.sh $RUNS/otvcl256 $RUNS/fesc256 0.15:2 0.12:3
# A run is named PREFIX_b<BETA>; give PREFIX_b<BETA>_f<F6> by setting
# NAME_F6=1 when two runs share a slope.
# Environment: DRYRUN=1 (make the directories, print the sbatch commands),
# SBATCH_ARGS (default: a whole node, 2 days), THREADS (default 192),
# PARAMS_EXTRA (more "key = value" lines, separated by ';').
set -eu
THREADS="${THREADS:-192}"
SBATCH_ARGS="${SBATCH_ARGS:---exclusive --cpus-per-task=192 --mem=200G --time=2-00:00:00}"
[ $# -ge 3 ] || { echo "usage: $0 BASE PREFIX F6:BETA [F6:BETA ...]"; exit 2; }
base="${1%/}"
prefix="${2%/}"
shift 2
[ -x build/radhydro ] || { echo "build/radhydro missing: run tools/hyak/build_hyak.sh first"; exit 1; }
[ -f "$base/params.txt" ] || { echo "$base: no params.txt"; exit 1; }
ics="$(cd "$base" && readlink -f ics)"
for spec in "$@"; do
  f6="${spec%%:*}"
  beta="${spec#*:}"
  dst="${prefix}_b${beta}"
  [ "${NAME_F6:-0}" = 1 ] && dst="${dst}_f${f6}"
  mkdir -p "$dst"
  if [ ! -f "$dst/params.txt" ]; then
    {
      echo "# $(basename "$base") with f_esc(z) = $f6 ((1+z)/7)^$beta (tools/hyak/fesc_runs.sh)"
      grep -Ev '^[[:space:]]*(ics|f_esc|f_esc_beta|f_esc_z0|tensor_dump)[[:space:]]*=' "$base/params.txt" || true
      echo "ics = $ics"
      echo "f_esc = $f6"
      echo "f_esc_beta = $beta"
      echo "f_esc_z0 = 6"
      [ -n "${PARAMS_EXTRA:-}" ] && echo "$PARAMS_EXTRA" | tr ';' '\n' | sed 's/^[[:space:]]*//'
    } > "$dst/params.txt"
  fi
  [ -e "$dst/ics" ] || ln -s "$ics" "$dst/ics"
  cmd="sbatch $SBATCH_ARGS --job-name=fesc_b${beta} --output=$dst/slurm-%j.out --export=ALL,PROB=$dst,THREADS=$THREADS tools/hyak/reion_production.slurm"
  if [ "${DRYRUN:-0}" = 1 ]; then echo "$cmd"; else $cmd; fi
done
