#!/bin/bash
# tau_eff distributions and the 21cmFAST comparison of the 512^3 suite
# (docs/taueff_21cmfast.md). One phase per job:
#
#   PHASE=taueff     sbatch -c 24 --mem=140G tools/hyak/taueff_21cmfast.slurm  # Part A: igm_stats.py --taueff, Figure 7
#   PHASE=slabs      sbatch -c 28 --mem=120G tools/hyak/taueff_21cmfast.slurm  # Part A4: movie slabs + calibration
#   PHASE=es         sbatch -c 10 --mem=40G  tools/hyak/taueff_21cmfast.slurm  # B1: the criterion on the suite's density
#   PHASE=standalone sbatch -c 32 --mem=200G tools/hyak/taueff_21cmfast.slurm  # B2: 21cmFAST with its own ICs
#   PHASE=figures    sbatch -c 4  --mem=60G  tools/hyak/taueff_21cmfast.slurm  # Figure 8 (needs es, standalone)
#
# Run from the checkout. The 21cmFAST phases use the conda env
# $RADRAY_SW/conda-envs/21cmfast (python 3.12, gsl, fftw, pip 21cmFAST==4.3.1).
#SBATCH -J taueff21
#SBATCH -A astro
#SBATCH -p cpu-g2
#SBATCH --time=8:00:00
#SBATCH -o taueff21_%j.log
set -eo pipefail
RUNS=${RUNS:-/gscratch/astro/$USER/radhydro/runs}
OUT=${OUT:-$RUNS/rt_closure_paper}
NAMES="rays rays_caseA rays_m1 rays_otvet rays_otvet_split rays_otvet_screen rays_otvet_cone"
LABELS=(rays rays_A+RR M1 OTVET split screen cone)
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1} OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
use_21cmfast() {
  source /mmfs1/sw/miniforge3/25.9.1-0/etc/profile.d/conda.sh
  conda activate "${RADRAY_SW:-/gscratch/astro/$USER/sw}/conda-envs/21cmfast"
}
case "$PHASE" in
  taueff)
    source tools/hyak/env.sh
    export OMP_NUM_THREADS=1
    mkdir -p "$OUT/taueff"
    for r in $NAMES; do
      python3 tools/analysis/igm_stats.py "$RUNS/reion512_$r" --zmin 4.9 --zmax 6.3 --nskew 8192 \
        --taueff 50 --lya-only --jobs 3 --seed 1 --out "$OUT/taueff/igm_taueff_$r.json" &
    done
    wait
    L=""; i=0
    for r in $NAMES; do L="$L $OUT/taueff/igm_taueff_$r.json:${LABELS[$i]}"; i=$((i + 1)); done
    python3 tools/analysis/taueff_plot.py $L --outdir "$OUT/taueff/fig"
    python3 tools/analysis/taueff_plot.py $L --outdir "$OUT/taueff/fig_box" --seg box
    ;;
  slabs)
    source tools/hyak/env.sh
    export OMP_NUM_THREADS=1
    i=0
    for r in $NAMES; do
      python3 tools/analysis/taueff_slabs.py "$RUNS/reion512_$r:${LABELS[$i]}" --out "$OUT/taueff/slabs" \
        --jobs 4 --calib &
      i=$((i + 1))
    done
    wait
    L=""; i=0     # raw tau_eff against z and lambda_mfp (A.6)
    for r in $NAMES; do L="$L ${LABELS[$i]}:$RUNS/rt_methods_512/igm/igm_reion512_$r.json"; i=$((i + 1)); done
    python3 tools/analysis/taueff_z.py --slabs "$OUT/taueff/slabs" --outdir "$OUT/taueff/fig_z" $L
    ;;
  es)
    use_21cmfast
    python tools/analysis/excursion_set.py run "$RUNS/reion512_rays" --steps 61,69,74,78,84 \
      --out "$OUT/es" --rmax 15 --rmax-mfp "$RUNS/rt_methods_512/igm/igm_reion512_rays.json"
    ;;
  standalone)
    use_21cmfast
    python tools/analysis/excursion_set.py standalone "$RUNS/reion512_rays" --steps 61,69,74,78,84 \
      --out "$OUT/es" --hii-dim 512 --dim 1024 --seed 12345
    ;;
  figures)
    source tools/hyak/env.sh
    L=""; i=0
    for r in $NAMES; do L="$L $RUNS/reion512_$r:${LABELS[$i]}"; i=$((i + 1)); done
    python3 tools/analysis/excursion_set.py figure "$OUT/es" $L --outdir "$OUT/es/fig"
    ;;
  *)
    echo "PHASE must be taueff, slabs, es, standalone or figures" >&2
    exit 1
    ;;
esac
