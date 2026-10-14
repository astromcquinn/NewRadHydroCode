#!/bin/bash
# The hydro accuracy suite on the AMReX backend with refinement (Stage F,
# docs/stage_f.md section 2) against the uniform unsplit runs, on one klone
# node. Builds nothing: BUILD must hold a RADRAY_AMREX build with
# tests/hydro_suite/hydro_suite_driver (tools/hyak/build_amrex.sh).
#
#   sbatch --export=ALL,BUILD=$PWD/build-e4,OUT=$RUNS/hydro_suite_amr \
#          [,TESTS=tubes,sedov][,TIER=quick][,JOBS=8] tools/hyak/hydro_suite_amr.slurm
#
# Runs, for each test: the uniform unsplit runs, then --amr 1 and --amr 2
# with REFINE_MASS (default 0.2: everything but the undisturbed state of a
# tube is refined; the Sedov shell at 4x is), results in OUT/<tier>_unsplit*/
# results/*.json; tools/analysis/hydro_suite_amr.py tabulates them.
#SBATCH --job-name=hydro-suite-amr
#SBATCH --account=astro
#SBATCH --partition=cpu-g2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=4:00:00
set -euo pipefail
# sbatch copies the script to the spool: the repository is the submit directory
if [ -n "${SLURM_SUBMIT_DIR:-}" ] && [ -f "$SLURM_SUBMIT_DIR/tools/hyak/env.sh" ]; then
  src="$SLURM_SUBMIT_DIR"
else
  src="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
fi
source "$src/tools/hyak/env.sh"
BUILD="${BUILD:-$src/build-e4}"
OUT="${OUT:-$RUNS/hydro_suite_amr}"
TESTS="${TESTS:-tubes,sedov}"
TIER="${TIER:-quick}"
JOBS="${JOBS:-${SLURM_CPUS_PER_TASK:-4}}"
REFINE_MASS="${REFINE_MASS:-0.2}"
drv="$BUILD/tests/hydro_suite/hydro_suite_driver"
[ -x "$drv" ] || { echo "no $drv"; exit 1; }
export OMP_NUM_THREADS=1
mkdir -p "$OUT"
run() {  # run <suffix> [extra run_suite args]
  local sfx="$1"; shift
  python3 "$src/tests/hydro_suite/run_suite.py" "$TESTS" --tier "$TIER" --scheme unsplit \
      --driver "$drv" --rundir "$OUT/${TIER}_unsplit$sfx" --results "$OUT/${TIER}_unsplit$sfx/results" \
      --jobs "$JOBS" --reuse "$@"
}
run ""
run "_amr1" --amr 1 --refine-mass "$REFINE_MASS"
run "_amr2" --amr 2 --refine-mass "$REFINE_MASS"
python3 "$src/tools/analysis/hydro_suite_amr.py" "$OUT" --tier "$TIER" --tests "$TESTS"
