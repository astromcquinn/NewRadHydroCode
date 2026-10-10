#!/bin/bash
# Analysis and side-by-side movies of the rays vs moment-RT runs
# (tools/hyak/rt_methods_runs.sh, docs/rt_moments.md section 6). Submit from
# the repository root, after (or with a dependency on) the four runs:
#   sbatch --dependency=afterany:<jobs> tools/hyak/rt_methods_post.slurm
# Environment: RUNS (default /gscratch/astro/$USER/radhydro/runs).
# Writes $RUNS/rt_methods_{c01,c001}/ (compare.txt, history.png, maps_z*.png)
# and $RUNS/rt_methods_{c01,c001}/compare.mp4 (+ compare_stills/).
#SBATCH --job-name=rt_methods_post
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
for pair in "reion256_c01 c01" "reion256 c001"; do
  set -- $pair
  base="$RUNS/$1" out="$RUNS/rt_methods_$2"
  mkdir -p "$out"
  python3 tools/analysis/rt_methods_runs.py "$base/run:rays" "${base}_m1/run:M1" \
      "${base}_otvet/run:OTVET" --outdir "$out"
  python3 tools/movies/compare_movie.py "$out/compare.mp4" "$base/run:rays" \
      "${base}_m1/run:M1" "${base}_otvet/run:OTVET" --jobs "${SLURM_CPUS_PER_TASK:-8}"
done
