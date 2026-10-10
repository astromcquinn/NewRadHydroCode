#!/bin/bash
# IGM diagnostics (tools/analysis/igm_stats.py, igm_plot.py;
# docs/igm_diagnostics.md) of finished radhydro runs, compared with
# data/igm/igm_z5_6.txt. Only reads the run directories. Submit from the
# repository root:
#   RUNLIST="reion256_f0.2:f0.2 reion256_prod:f0.28" sbatch tools/hyak/igm_diag.slurm
# Environment: RUNS (default /gscratch/astro/$USER/radhydro/runs), RUNLIST
# (run names under $RUNS, each with an optional :label), OUT (default
# $RUNS/igm_diag), NLOS, NSKEW (sightlines per start set, skewers per axis;
# default 8192 each), CONV=1 also runs the convergence checks (sightline
# step, number of sightlines and skewers) on the first run's last snapshot.
# At N = 512 a snapshot takes ~1.5 min plus one scan of the 19 GB rt_grid and
# ~6 GB of resident memory (mostly mapped file pages): use JOBS=4 with 64 GB.
#SBATCH --job-name=igm_diag
#SBATCH --account=astro
#SBATCH --partition=cpu-g2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=/gscratch/astro/mcquinn/radhydro/logs/%x-%j.out
set -u
cd "${SLURM_SUBMIT_DIR}"
source tools/hyak/env.sh
RUNS="${RUNS:-/gscratch/astro/$USER/radhydro/runs}"
RUNLIST="${RUNLIST:-reion256_f0.2:f0.2 reion256_prod:f0.28 reion256_f0.1:f0.1}"
OUT="${OUT:-$RUNS/igm_diag}"
NLOS="${NLOS:-8192}"
NSKEW="${NSKEW:-8192}"
JOBS="${JOBS:-8}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p "$OUT"
specs=()
for r in $RUNLIST; do
    name="${r%%:*}"
    label="${r#*:}"
    python3 tools/analysis/igm_stats.py "$RUNS/$name" --nlos "$NLOS" --nskew "$NSKEW" \
        --jobs "$JOBS" --out "$OUT/igm_$name.json"
    specs+=("$OUT/igm_$name.json:$label")
done
python3 tools/analysis/igm_plot.py "${specs[@]}" --out "$OUT/igm.png"

if [ "${CONV:-0}" = 1 ]; then
    first="${RUNLIST%% *}"
    first="${first%%:*}"
    mkdir -p "$OUT/conv"
    step=$(ls -d "$RUNS/$first"/run/step_* | while read d; do
        [ -f "$d/gamma.bin" ] && basename "$d"; done | tail -1)
    step=$((10#${step#step_}))
    c() { python3 tools/analysis/igm_stats.py "$RUNS/$first" --steps "$step" "$@"; }
    c --nlos 4096 --nskew 4096 --out "$OUT/conv/base.json" &
    c --nlos 4096 --nskew 4096 --ds 0.25 --out "$OUT/conv/ds025.json" &
    c --nlos 32768 --nskew 4096 --out "$OUT/conv/nlos32k.json" &
    c --nlos 4096 --nskew 65536 --out "$OUT/conv/allskew.json" &
    c --nlos 4096 --nskew 4096 --seed 2 --out "$OUT/conv/seed2.json" &
    c --nlos 4096 --nskew 4096 --rmax 80 --out "$OUT/conv/rmax80.json" &
    c --nlos 4096 --nskew 4096 --lya-profile point --out "$OUT/conv/point.json" &
    wait
fi
