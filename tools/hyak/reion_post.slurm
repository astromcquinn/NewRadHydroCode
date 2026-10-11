#!/bin/bash
# Post-processing of a production run (tools/hyak/reion_production.slurm):
# stills (tools/movies/stills.py), the reionization history and tau_es
# (tools/analysis/reion.py) and the movie (tools/movies/make_movie.py).
#
# reion_production.slurm runs it after every successful leg: in its own job
# with POST=1 (the default), or submits it as this small separate job with
# POST=sbatch, so the big node is released at once. By hand, from the
# repository root:
#   sbatch --export=ALL,RUN=/gscratch/astro/$USER/radhydro/runs/reion256/run \
#          tools/hyak/reion_post.slurm
# or, inside an allocation, RUN=... bash tools/hyak/reion_post.slurm
#
# Environment: RUN (required: the run directory, <problem>/run), JOBS (movie
# render processes, default: the job's cpus, else 8), MOVIE_ARGS (extra
# make_movie.py arguments). Works with full-cube and slab movie frames
# (movie_slab, docs/movies.md).
#SBATCH --job-name=reion_post
#SBATCH --account=astro
#SBATCH --partition=cpu-g2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=%x-%j.out
set -u
cd "${SLURM_SUBMIT_DIR:-.}"
source tools/hyak/env.sh
if [ -n "${HYAK_MODULES:-}" ]; then module load $HYAK_MODULES; fi
RUN="${RUN:?set RUN to the run directory}"
JOBS="${JOBS:-${SLURM_CPUS_PER_TASK:-8}}"
# one thread per process: the renderers are separate processes, and the
# production job exports OMP_NUM_THREADS = its cores
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
echo "$(date) post-processing $RUN on $(hostname), $JOBS render processes"
t0=$SECONDS
# Stills and the ionisation history (tools/movies/stills.py needs the
# movie fields fHI, temp and rhog)
if ls "$RUN"/movie/frame_*_rhog.bin >/dev/null 2>&1 || ls "$RUN"/movie/frame_*_rhog.slab >/dev/null 2>&1; then
  python3 tools/movies/stills.py "$RUN" "$RUN/stills" || true
  echo "  stills.py done ($((SECONDS - t0)) s)"
fi
# Reionization history and tau_es (docs/production_run.md section 8), and
# the movie (docs/movies.md)
python3 tools/analysis/reion.py "$RUN" --out "$RUN/reion.txt" || true
echo "  reion.py done ($((SECONDS - t0)) s)"
python3 tools/movies/make_movie.py "$RUN" --jobs "$JOBS" ${MOVIE_ARGS:-} || true
echo "$(date) post-processing done ($((SECONDS - t0)) s)"
