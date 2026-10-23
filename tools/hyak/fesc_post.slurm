#!/bin/bash
# Analysis of a falling-f_esc 512^3 run against the constant-f_esc rays
# (docs/fesc_evolution.md section 4). Only reads the run directories.
# Submit from the repository root:
#   sbatch tools/hyak/fesc_post.slurm
# Environment:
#   RUNS  default /gscratch/astro/$USER/radhydro/runs
#   NEW   the new run (default fesc512_b2.5), label NEWLAB (default fesc)
#   REF   the reference rays run (default reion512_rays; its extra field
#         cubes in <REF>/run_cubes are used for P_Gamma when present),
#         label REFLAB (default rays)
#   SPREAD the rays' own spread for P_Gamma (default reion512_rays_caseB;
#         empty = none)
#   REFTITLE, NEWTITLE the movie's column titles (default "rays, f_esc = 0.23"
#         and "rays, f_esc(z)")
#   For the refined 512^3 run against its control (docs/stage_f.md section 4):
#   NEW=fesc512_b2.5_amr NEWLAB=amr REF=fesc512_b2.5_uniform REFLAB=uniform
#   SPREAD= REFTITLE="uniform 512^3" NEWTITLE="512^3 + 2 levels" IGM=0
#   OUT=$RUNS/fesc_evolution/512_amr sbatch tools/hyak/fesc_post.slurm
#   OUT   default $RUNS/fesc_evolution/512
#   IGM=0, PK=0, MOVIE=0, MAPS=0 skip a part
# Writes OUT/igm/ (igm_stats.py json, igm_plot.py figure and table),
# OUT/gamma_pk*/ (gamma_power.py, also capped at 3 <Gamma>), OUT/compare.mp4
# (compare_movie.py), OUT/{compare.txt,history.png,maps_*.png,gamma_pdf.png}
# (rt_methods_runs.py).
#SBATCH --job-name=fesc_post
#SBATCH --account=astro
#SBATCH --partition=cpu-g2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --output=/gscratch/astro/mcquinn/radhydro/logs/%x-%j.out
set -u
cd "${SLURM_SUBMIT_DIR}"
source tools/hyak/env.sh
RUNS="${RUNS:-/gscratch/astro/$USER/radhydro/runs}"
NEW="${NEW:-fesc512_b2.5}" NEWLAB="${NEWLAB:-fesc}"
REF="${REF:-reion512_rays}" REFLAB="${REFLAB:-rays}" SPREAD="${SPREAD-reion512_rays_caseB}"
REFTITLE="${REFTITLE:-rays, f_esc = 0.23}" NEWTITLE="${NEWTITLE:-rays, f_esc(z)}"
OUT="${OUT:-$RUNS/fesc_evolution/512}"
JOBS="${SLURM_CPUS_PER_TASK:-8}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p "$OUT"
ref_run="$RUNS/$REF/run"
[ -d "$RUNS/$REF/run_cubes" ] && ref_run="$ref_run+$RUNS/$REF/run_cubes"
spread=()
if [ -n "$SPREAD" ]; then
  spr_run="$RUNS/$SPREAD/run"
  [ -d "$RUNS/$SPREAD/run_cubes" ] && spr_run="$spr_run+$RUNS/$SPREAD/run_cubes"
  spread=("$spr_run:rays_B")
fi

if [ "${IGM:-1}" = 1 ]; then
  mkdir -p "$OUT/igm"
  python3 tools/analysis/igm_stats.py "$RUNS/$NEW" --nlos 8192 --nskew 8192 --jobs 4 \
      --out "$OUT/igm/igm_$NEW.json"
fi
if [ "${PK:-1}" = 1 ]; then
  for clip in 0 3; do
    c=(); [ "$clip" != 0 ] && c=(--clip "$clip")
    python3 tools/analysis/gamma_power.py "$ref_run:$REFLAB" "$RUNS/$NEW/run:$NEWLAB" "${spread[@]}" \
        --outdir "$OUT/gamma_pk" --box 40 --zmin 4.9 --zmax 5.6 --tol 0.05 --workers 4 "${c[@]}"
  done
fi
if [ "${MAPS:-1}" = 1 ]; then
  python3 tools/analysis/rt_methods_runs.py "$RUNS/$REF/run:$REFLAB" "$RUNS/$NEW/run:$NEWLAB" --outdir "$OUT"
fi
if [ "${MOVIE:-1}" = 1 ]; then
  python3 tools/movies/compare_movie.py "$OUT/compare.mp4" "$RUNS/$REF/run:$REFTITLE" \
      "$RUNS/$NEW/run:$NEWTITLE" --jobs "$JOBS"
fi
