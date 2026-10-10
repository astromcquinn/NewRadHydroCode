#!/bin/bash
# Analysis and side-by-side movie of case B vs case A vs case A +
# recombination radiation (tools/hyak/recomb_runs.sh,
# docs/recomb_radiation.md section 4). Submit from the repository root,
# after (or with a dependency on) the runs:
#   sbatch --dependency=afterany:<jobs> tools/hyak/recomb_post.slurm
# Environment: RUNS (default /gscratch/astro/$USER/radhydro/runs), BASE
# (default reion256_c01: compares $BASE, ${BASE}_caseA, ${BASE}_caseA_rr),
# RR (default caseA_rr: the recombination-radiation sibling), OUT (default
# $RUNS/recomb_$BASE).
# Writes $OUT/compare.txt, history.png, maps_z*.png, compare.mp4 (+ stills).
#SBATCH --job-name=recomb_post
#SBATCH --account=astro
#SBATCH --partition=cpu-g2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=%x-%j.out
set -u
cd "${SLURM_SUBMIT_DIR}"
source tools/hyak/env.sh
RUNS="${RUNS:-/gscratch/astro/$USER/radhydro/runs}"
BASE="${BASE:-reion256_c01}"
RR="${RR:-caseA_rr}"
OUT="${OUT:-$RUNS/recomb_$BASE}"
b="$RUNS/$BASE"
mkdir -p "$OUT"
python3 tools/analysis/rt_methods_runs.py "$b/run:caseB" "${b}_caseA/run:caseA" \
    "${b}_${RR}/run:caseA+RR" --outdir "$OUT"
python3 tools/movies/compare_movie.py "$OUT/compare.mp4" "$b/run:case B" \
    "${b}_caseA/run:case A" "${b}_${RR}/run:case A + RR" --jobs "${SLURM_CPUS_PER_TASK:-8}"
