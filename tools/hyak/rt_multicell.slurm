#!/bin/bash
# Rays crossing k RT cells per substep (rt_cells_per_substep,
# docs/rt_multicell.md): the step-4 comparison of the production
# configuration at 128^3 on one klone cpu-g2 node, in case the 64^3 result
# is marginal. Runs k = 1, 2, 4 (KS) from the same ICs one after the other
# on the whole node, then tools/analysis/rt_multicell_compare.py.
#
# Submit from the repository root (build first: tools/hyak/build_hyak.sh):
#   sbatch --export=ALL,ROOT=/gscratch/astro/$USER/radhydro/runs/multicell128 \
#          tools/hyak/rt_multicell.slurm
# The scheme recommended by docs/rt_multicell.md (exposure deposit and the
# converged photon limit) in a second directory, with the same ICs (the more
# useful of the two if only one is run; leave out ICS= to make new ICs):
#   sbatch --export=ALL,ROOT=$RUNS/multicell128_px,ICS=$RUNS/multicell128/prob_k1/ics,\
#          EXTRA_PROD="--photon-conserve 1 --ray-deposit exposure" tools/hyak/rt_multicell.slurm
# The step-5 repeat at full c (fast fronts in a coarse step; about 10x the
# cost of 0.1c, so stop at z = 9.5, where the 64^3 run had x_HII,V ~ 0.8):
#   sbatch --export=ALL,ROOT=$RUNS/multicell128_c1,C=1,KS="1 4",ZDATA1="15 14 13 12 11 10 9.5",\
#          ICS=$RUNS/multicell128/prob_k1/ics,EXTRA_PROD="--photon-conserve 1 --ray-deposit exposure" \
#          tools/hyak/rt_multicell.slurm
#
# Environment: ROOT (required; created), N (default 128), C (c_reduced,
# default 0.1), KS (default "1 2 4"), THREADS (default: the job's cpus),
# ICS (an existing ics/ directory of the same N to link instead of making
# new ones, e.g. $RUNS/reion128/ics), EXTRA_PROD (more make_production.py
# arguments, e.g. "--am-params data/uv_lf/am_params_128_10mpc.txt
# --f-esc 0.2"), ZDATA1 (the output redshifts, default make_production.py's
# 15 ... 5; the run stops at the last), EXE (default build/radhydro).
#
# Layout: $ROOT/prob_k<k>/params.txt (make_production.py --N $N
# --c-reduced $C --hydro-scheme unsplit --movie-dz 0.05 --movie-slab 8
# --movie-fields fHI --rt-cells-per-substep <k>, the same ICs) and the runs
# $ROOT/cmp/base (k = 1), $ROOT/cmp/k2, $ROOT/cmp/k4. Each run restarts from
# its checkpoint if the job is resubmitted. Results: $ROOT/cmp.json,
# $ROOT/cmp_history.png, $ROOT/cmp_slices.png and the table in this job's
# log. Expected cost (docs/rt_multicell.md, the 64^3 timings scaled by 8 in
# cells and ~2 in substeps): k = 1 about 30-40 min on 192 threads, k = 2 and
# 4 less.
#SBATCH --job-name=multicell
#SBATCH --account=astro
#SBATCH --partition=cpu-g2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH --cpus-per-task=192
#SBATCH --mem=256G
#SBATCH --time=1-00:00:00
#SBATCH --open-mode=append
#SBATCH --output=%x-%j.out
set -u
cd "${SLURM_SUBMIT_DIR}"
source tools/hyak/env.sh
if [ -n "${HYAK_MODULES:-}" ]; then module load $HYAK_MODULES; fi
ROOT="${ROOT:?set ROOT to the output directory}"
N="${N:-128}"
C="${C:-0.1}"
KS="${KS:-1 2 4}"
EXE="${EXE:-build/radhydro}"
export OMP_NUM_THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-192}}" OMP_PLACES=cores OMP_PROC_BIND=close
ulimit -s unlimited
mkdir -p "$ROOT/cmp"
echo "$(date) rt_multicell: N = $N, c_reduced = $C, k = $KS, $OMP_NUM_THREADS threads on $(hostname)"

first=""
for k in $KS; do
  P="$ROOT/prob_k$k"
  if [ ! -f "$P/params.txt" ]; then
    mkdir -p "$P"
    noics=""
    if [ -n "$first" ] || [ -n "${ICS:-}" ]; then noics="--no-ics"; fi
    python3 tools/production/make_production.py "$P" --N "$N" --c-reduced "$C" \
        --hydro-scheme unsplit --movie-dz 0.05 --movie-slab 8 --movie-fields fHI \
        --restart-every 4 --rt-cells-per-substep "$k" $noics ${EXTRA_PROD:-} \
        ${ZDATA1:+--zdata1 "$ZDATA1"} || exit 1
    if [ -n "${ICS:-}" ]; then
      ln -sfn "$ICS" "$P/ics"
    elif [ -n "$first" ]; then
      ln -sfn "$first/ics" "$P/ics"
    fi
  fi
  [ -n "$first" ] || first="$P"
  if [ "$k" = 1 ]; then R="$ROOT/cmp/base"; else R="$ROOT/cmp/k$k"; fi
  if [ -e "$R/DONE" ]; then echo "k = $k already done"; continue; fi
  args=("$P" "$R" --quiet)
  CK="$R/restart/checkpoint.bin"
  if [ -e "$CK" ]; then args+=(--restart "$CK"); fi
  echo "$(date) k = $k -> $R"
  /usr/bin/time -v "$EXE" "${args[@]}" >> "$R.log" 2>&1
  rc=$?
  echo "$(date) k = $k: radhydro exit $rc"
  grep -E "Maximum resident|Elapsed \(wall" "$R.log" | tail -2
  if [ $rc -ne 0 ]; then exit $rc; fi
  touch "$R/DONE"
done

python3 tools/analysis/rt_multicell_compare.py "$ROOT/cmp" --N "$N" --json "$ROOT/cmp.json" \
    --fig "$ROOT/cmp"
for d in "$ROOT"/cmp/*/; do
  echo "== $d"
  python3 tools/analysis/reion.py "$d" | tail -3
done
echo "$(date) done"
