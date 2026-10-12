#!/bin/bash
# Extra field cubes at every hydro step of an existing run, written into a
# separate output directory (the run itself is not touched): restart from one
# of its checkpoints (bit-identical continuation) or rerun from the start,
# with fields_every = 1, no movie frames and no checkpoints. Used for the
# Gamma_HI power spectra after overlap (docs/rt_methods_512.md).
# Submit from the repository root on a whole node:
#   sbatch --exclusive --cpus-per-task=192 --mem=300G --time=6:00:00 \
#     --export=ALL,PROB=$RUNS/reion512_rays_m1,CK=$RUNS/reion512_rays_m1/run/restart/checkpoint_0048.bin \
#     tools/hyak/field_cubes.slurm
# Environment: PROB (required), CK (checkpoint; empty: from the start), OUT
# (default $PROB/run_cubes), FIELDS (default fHI,rhog,gammaHI; export it before sbatch,
# since --export=...,FIELDS=a,b splits at the comma), THREADS, EXE.
#SBATCH --job-name=field_cubes
#SBATCH --account=astro
#SBATCH --partition=cpu-g2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --mem=300G
#SBATCH --time=6:00:00
#SBATCH --output=%x-%j.out
set -u
cd "${SLURM_SUBMIT_DIR}"
source tools/hyak/env.sh
PROB="${PROB:?set PROB}"
OUT="${OUT:-$PROB/run_cubes}"
EXE="${EXE:-build/radhydro}"
export OMP_NUM_THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-192}}" OMP_PLACES=cores OMP_PROC_BIND=close
ulimit -s unlimited
args=("$PROB" "$OUT" --quiet --set fields_every=1 --set "fields=${FIELDS:-fHI,rhog,gammaHI}"
      --set movie_dz=0 --set restart_every=0)
[ -n "${CK:-}" ] && args+=(--restart "$CK")
echo "$(date) $EXE ${args[*]} on $(hostname), $OMP_NUM_THREADS threads"
/usr/bin/time -v "$EXE" "${args[@]}" >> "$OUT.log" 2>&1
rc=$?
echo "$(date) radhydro exit $rc"
grep -E "Maximum resident|Elapsed \(wall" "$OUT.log" | tail -2
exit $rc
