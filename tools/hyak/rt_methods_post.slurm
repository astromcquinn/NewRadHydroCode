#!/bin/bash
# Analysis and side-by-side movies of the rays vs moment-RT runs
# (tools/hyak/rt_methods_runs.sh, docs/rt_moments.md section 6,
# docs/rt_methods_512.md). Submit from the repository root, after (or with a
# dependency on) the runs:
#   sbatch --dependency=afterany:<jobs> tools/hyak/rt_methods_post.slurm
# Environment:
#   RUNS  default /gscratch/astro/$USER/radhydro/runs
#   SETS  "base:outname ..." (default "reion256_c01:rt_methods_c01
#         reion256:rt_methods_c001"): compares $RUNS/<base> (rays) with
#         <base>_m1 and <base>_otvet into $RUNS/<outname>/
#   EXTRA "run:label ..." further runs added to every set's comparison and
#         movie (e.g. "reion512_rays_caseB:rays_B")
#   IGM=1 also igm_stats.py on every run (z <= 7.1 snapshots) and
#         igm_plot.py against data/igm/igm_z5_6.txt into $RUNS/<outname>/igm/
#   MOVIE=0 skip the movie
# Writes compare.txt, history.png, maps_z*.png, maps_gamma_z*.png,
# gamma_pdf.png, compare.mp4 + compare_stills/.
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
SETS="${SETS:-reion256_c01:rt_methods_c01 reion256:rt_methods_c001}"
JOBS="${SLURM_CPUS_PER_TASK:-8}"
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
for pair in $SETS; do
  base="$RUNS/${pair%%:*}" out="$RUNS/${pair#*:}"
  mkdir -p "$out"
  specs=("$base/run:rays" "${base}_m1/run:M1" "${base}_otvet/run:OTVET")
  igm=("${pair%%:*}:rays" "${pair%%:*}_m1:M1" "${pair%%:*}_otvet:OTVET")
  for e in ${EXTRA:-}; do
    specs+=("$RUNS/${e%%:*}/run:${e#*:}")
    igm+=("$e")
  done
  python3 tools/analysis/rt_methods_runs.py "${specs[@]}" --outdir "$out"
  if [ "${MOVIE:-1}" != 0 ]; then
    python3 tools/movies/compare_movie.py "$out/compare.mp4" "${specs[@]}" --jobs "$JOBS"
  fi
  if [ "${IGM:-0}" = 1 ]; then
    mkdir -p "$out/igm"
    js=()
    for r in "${igm[@]}"; do
      name="${r%%:*}"
      OMP_NUM_THREADS=1 python3 tools/analysis/igm_stats.py "$RUNS/$name" --nlos 8192 --nskew 8192 \
          --jobs 4 --out "$out/igm/igm_$name.json"
      js+=("$out/igm/igm_$name.json:${r#*:}")
    done
    python3 tools/analysis/igm_plot.py "${js[@]}" --out "$out/igm/igm.png" --table "$out/igm/igm.txt"
  fi
done
