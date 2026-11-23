#!/bin/bash
# The photon sink by density of the uniform 512^3 RT-methods suite
# (tools/analysis/absorber_budget.py, docs/absorbers_closures.md): the
# rays (the reference) against the moment methods at the snapshots nearest
# the redshifts of ZLIST. Only reads the run directories. Submit from the
# repository root:
#   sbatch tools/hyak/absorber_budget.slurm
# Environment: RUNS (default /gscratch/astro/$USER/radhydro/runs), RUNLIST
# (run names under $RUNS with a :label each; the first is the reference),
# OUT (default $RUNS/absorbers_wu21), ZLIST, ABS (the reference's absorber
# redshifts), FIELD (gammaHI_avg), JOBS. A 512^3 snapshot takes ~2 min and
# ~2 GB per job (the absorber step ~5 GB).
#SBATCH --job-name=absorber_budget
#SBATCH --account=astro
#SBATCH --partition=cpu-g2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=03:00:00
#SBATCH --output=/gscratch/astro/mcquinn/radhydro/logs/%x-%j.out
set -u
cd "${SLURM_SUBMIT_DIR}"
source tools/hyak/env.sh
RUNS="${RUNS:-/gscratch/astro/$USER/radhydro/runs}"
RUNLIST="${RUNLIST:-reion512_rays:rays reion512_rays_m1:M1 reion512_rays_otvet:OTVET reion512_rays_otvet_cone:cone reion512_rays_otvet_split:split reion512_rays_otvet_screen:screen}"
OUT="${OUT:-$RUNS/absorbers_wu21}"
ZLIST="${ZLIST:-5.98,5.49,5.25,4.99}"
ABS="${ABS:-5.49,4.99}"
FIELD="${FIELD:-gammaHI_avg}"
JOBS="${JOBS:-8}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p "$OUT"
specs=()
for r in $RUNLIST; do
    specs+=("$RUNS/${r%%:*}:${r#*:}")
done
echo "host $(hostname -s), $(date)"
/usr/bin/time -v python3 tools/analysis/absorber_budget.py "${specs[@]}" --outdir "$OUT" \
    --zlist "$ZLIST" --absorbers "$ABS" --field "$FIELD" --jobs "$JOBS"
echo "done $(date)"
