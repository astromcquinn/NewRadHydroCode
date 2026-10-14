#!/bin/bash
# The PM backend (pm_backend = amrex, docs/stage_e_pm.md section 4.1) at
# pm_refine = 2 on klone: STEPS hydro steps restarted from a checkpoint with
# the three backends and the state on the levels (rt_backend = hydro_backend
# = pm_backend = amrex, amrex_state = 1), 24 ranks x 8 threads per cpu-g2
# node, for every (amrex_pm_fft, amrex_max_grid_size) of CONFIGS; on one node
# also the serial code on all 192 cores. Every backend run is compared
# bitwise (tests/amrex_pm_check.py bitwise: every file of every step) with
# the reference run of its FFT (the two FFTs differ at round-off by design,
# docs/stage_e_pm.md section 1.2), and every output file (steps, DM, Pk, halos) is
# hashed into <run>/sha256.txt; the big files of the other runs are then
# removed (KEEP=1 keeps them). The per-step timers are the "step NNNN
# timers:" lines of <run>.log (radhydro --set step_timers=1, cumulative).
#
#   PROB=$RUNS/hmf512_40mpc_r2 CK=$PROB/run/restart/checkpoint.bin
#   sbatch --nodes=1 --ntasks=24 --export=ALL,PROB=$PROB,CK=$CK,OUT=... tools/hyak/amrex_pm_scaling.slurm
#   sbatch --nodes=2 --ntasks=48 --dependency=afterok:<n1 job> --export=... tools/hyak/amrex_pm_scaling.slurm
#
# Environment:
#   PROB, CK   problem directory and checkpoint (required)
#   OUT        output root (default $PROB/amrex_pm_scaling); runs in
#              $OUT/n<nodes>_<fft>_mgs<mgs>, the serial run in $OUT/serial
#   STEPS      hydro steps (default 3)
#   CONFIGS    "fft:mgs ..." (default "pencil:32 amrex:32 pencil:64 amrex:64")
#   REF        the reference run (default $OUT/n1_<fft>_mgs32 for each FFT,
#              which the one-node job writes first and keeps)
#   ZDATA1     the output redshifts of the restart (default "15 5.5 5.35 4":
#              one zdata1 output inside the 3 steps, for the P(k) time and
#              the DM file, and the run continues past the original end)
#   FIELDS     fields written every step (default u_gas,e_gas,rhog: the
#              gas the PM kicks, and the density)
#   RANKS, THREADS  (default the job's --ntasks, 8)
#   SERIAL     1: also the serial code (default on one node)
#   EXE        build with the backends (default build-amrex/radhydro)
#   KEEP       1: keep every run's output files
#SBATCH --job-name=amrex-pm-scaling
#SBATCH --account=astro
#SBATCH --partition=cpu-g2
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=8
#SBATCH --mem=700G
#SBATCH --time=06:00:00
#SBATCH --exclusive
#SBATCH --output=%x-%j.out
set -u
cd "${SLURM_SUBMIT_DIR}"
source tools/hyak/env.sh
radray_module_load ${HYAK_MODULES:-ompi/4.1.6-2}
radray_ucx_eager
PROB="${PROB:?set PROB}"
CK="${CK:?set CK}"
OUT="${OUT:-$PROB/amrex_pm_scaling}"
STEPS="${STEPS:-3}"
NODES="${SLURM_JOB_NUM_NODES:-1}"
RANKS="${RANKS:-${SLURM_NTASKS:-24}}"
THREADS="${THREADS:-8}"
CONFIGS="${CONFIGS:-pencil:32 amrex:32 pencil:64 amrex:64}"
REF_SET="${REF:-}"
ZDATA1="${ZDATA1:-15 5.5 5.35 4}"
FIELDS="${FIELDS:-u_gas,e_gas,rhog}"
EXE="${EXE:-$PWD/build-amrex/radhydro}"
SERIAL="${SERIAL:-$([ "$NODES" = 1 ] && echo 1 || echo 0)}"
KEEP="${KEEP:-0}"
[ -x "$EXE" ] || { echo "no executable $EXE"; exit 1; }
[ -f "$CK" ] || { echo "no checkpoint $CK"; exit 1; }
mkdir -p "$OUT"
echo "host $(hostname), nodes ${SLURM_JOB_NODELIST:-}, job ${SLURM_JOB_ID:-}, $(date)"
echo "EXE $EXE ($(git rev-parse --short HEAD)), $NODES nodes, $RANKS ranks x $THREADS threads, configs $CONFIGS, steps $STEPS"

common=(--quiet --restart "$CK" --set run_steps="$STEPS" --set restart_every=0 --set movie_dz=0
        --set fields_every=1 --set fields="$FIELDS" --set zdata1="$ZDATA1" --set step_timers=1)

# sha256 of every output file of a run (relative paths, sorted)
hash_run() {
  ( cd "$1" && find step_* DM Pk halos -type f 2>/dev/null | sort | xargs -r sha256sum ) > "$1/sha256.txt"
}
# remove the big files of a run (fields, particles), keeping stats, P(k), catalogues and hashes
trim_run() {
  [ "$KEEP" = 1 ] && return
  rm -f "$1"/step_*/*.bin "$1"/DM/*
}

summary="$OUT/summary_n$NODES.txt"
: > "$summary"
export OMP_NUM_THREADS=$THREADS OMP_PLACES=cores OMP_PROC_BIND=close
for cfg in $CONFIGS; do
  fft=${cfg%%:*}; mgs=${cfg##*:}
  run="$OUT/n${NODES}_${fft}_mgs${mgs}"
  REF="${REF_SET:-$OUT/n1_${fft}_mgs32}"
  rm -rf "$run"
  t0=$(date +%s)
  srun --ntasks="$RANKS" --cpus-per-task="$THREADS" --cpu-bind=cores --kill-on-bad-exit=1 \
    "$EXE" "$PROB" "$run" "${common[@]}" --set write_catalogs="$run/halos" \
    --set rt_backend=amrex --set hydro_backend=amrex --set pm_backend=amrex \
    --set amrex_hydro_max_grid_size=0 --set amrex_state=1 --set deterministic_parallel=0 \
    --set amrex_max_grid_size="$mgs" --set amrex_pm_fft="$fft" > "$run.log" 2>&1
  rc=$?
  echo "$run rc=$rc wall $(( $(date +%s) - t0 )) s" | tee -a "$summary"
  tail -4 "$run.log"
  [ $rc = 0 ] || continue
  hash_run "$run"
  if [ "$run" != "$REF" ]; then
    if [ -d "$REF" ]; then
      echo "bitwise $run vs $REF:" | tee -a "$summary"
      python3 tests/amrex_pm_check.py bitwise "$REF" "$run" 2>&1 | tee -a "$summary"
      if cmp -s <(grep -v ' step_' "$REF/sha256.txt") <(grep -v ' step_' "$run/sha256.txt"); then
        echo "  DM, Pk, halos: IDENTICAL" | tee -a "$summary"
      else
        echo "  DM, Pk, halos: DIFFERENT" | tee -a "$summary"
        diff <(grep -v ' step_' "$REF/sha256.txt") <(grep -v ' step_' "$run/sha256.txt") | head | tee -a "$summary"
      fi
    else
      echo "no reference $REF: hashes only" | tee -a "$summary"
    fi
    trim_run "$run"
  fi
done

if [ "$SERIAL" = 1 ]; then
  export OMP_NUM_THREADS=$((RANKS * THREADS))
  run="$OUT/serial"
  rm -rf "$run"
  t0=$(date +%s)
  srun --ntasks=1 --cpus-per-task="$OMP_NUM_THREADS" --cpu-bind=none \
    "$EXE" "$PROB" "$run" "${common[@]}" --set write_catalogs="$run/halos" > "$run.log" 2>&1
  rc=$?
  echo "$run rc=$rc wall $(( $(date +%s) - t0 )) s" | tee -a "$summary"
  tail -4 "$run.log"
  [ $rc = 0 ] && hash_run "$run"
fi
echo "--- peak RSS of the job steps (sacct, MaxRSS = the largest task):"
sleep 5; sacct -j "$SLURM_JOB_ID" --format=JobID%18,JobName%12,NTasks,MaxRSS,MaxRSSTask,AveRSS,Elapsed -P 2>/dev/null
