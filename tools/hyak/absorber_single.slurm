#!/bin/bash
# Phase 2 (c) of docs/absorbers_closures.md: one r^-2 absorber in a lattice
# of sources with every RT method (tools/analysis/absorber_single.py), at
# several resolutions, with apps/rt_driver. Submit from the repository root:
#   sbatch tools/hyak/absorber_single.slurm
#   NLIST=16 NSTEPS=3 METHODS=rays,m1 OUT=$RUNS/absorbers_wu21/single_smoke sbatch ... # a smoke test
# Environment: RUNS, EXE (default $RUNS/absorbers_wu21_bin/rt_driver, a build
# with P11), OUT (default $RUNS/absorbers_wu21/single), NLIST (default
# 32,64,128: 32, 16 and 8 pkpc cells in the 5 Mpc/h box at z = 6), METHODS,
# NSTEPS (40), CALSTEPS (0 = NSTEPS), THREADS (the job's cpus), EXTRA (more
# absorber_single.py run arguments).
#SBATCH --job-name=absorber_single
#SBATCH --account=astro
#SBATCH --partition=cpu-g2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=08:00:00
#SBATCH --output=/gscratch/astro/mcquinn/radhydro/logs/%x-%j.out
set -u
cd "${SLURM_SUBMIT_DIR}"
source tools/hyak/env.sh
RUNS="${RUNS:-/gscratch/astro/$USER/radhydro/runs}"
EXE="${EXE:-$RUNS/absorbers_wu21_bin/rt_driver}"
OUT="${OUT:-$RUNS/absorbers_wu21/single}"
NLIST="${NLIST:-32,64,128}"
METHODS="${METHODS:-rays,m1,otvet,cone,split,screen,rays_p110}"
NSTEPS="${NSTEPS:-40}"
CALSTEPS="${CALSTEPS:-0}"
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-16}}"
export OMP_PLACES=cores OMP_PROC_BIND=close
mkdir -p "$OUT"
echo "host $(hostname -s), $(date), exe $EXE"
/usr/bin/time -v python3 tools/analysis/absorber_single.py run --exe "$EXE" --outdir "$OUT" \
    --N "$NLIST" --methods "$METHODS" --nsteps "$NSTEPS" --cal-steps "$CALSTEPS" \
    --threads "$THREADS" ${EXTRA:-}
echo "done $(date)"
