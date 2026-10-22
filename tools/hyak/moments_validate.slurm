#!/bin/bash
# The moment RT on the AMReX level against a serial production run
# (docs/stage_e_moments.md section 7, phase 5): restart the serial run's
# checkpoint of step STEP on the level (all three backends, the state on
# z-slabs), run NSTEPS hydro steps, and compare every file the level run
# writes for those steps (step_NNNN/, movie frames) byte for byte with the
# serial run's own. Prints the level's timers and the peak RSS over the
# ranks (also in OUT/timing.txt).
#
#   sbatch --nodes=2 --ntasks=32 --cpus-per-task=12 --mem=1400G --exclusive \
#          --export=ALL,SRC=/gscratch/astro/$USER/radhydro/runs/otvcl256_otvet_cone,STEP=40,NSTEPS=1,OUT=... \
#          tools/hyak/moments_validate.slurm
#
# SRC     the serial run (SRC/params.txt, SRC/run/restart/checkpoint_STEP.bin,
#         SRC/run/step_NNNN)
# STEP    the checkpoint's step; NSTEPS the steps to run (default 1)
# OUT     the level run's directory (default SRC/level_<job id>)
# EXE     radhydro of an AMReX build (default build-amrex/radhydro)
# PLANES  amrex_slab_planes (default 0: about N / ranks)
# STATS   amrex_moment_stats: exact (default; the bitwise comparison) | fixed
# EXTRA   further --set arguments
# COMPARE 0: no comparison (timing runs past the serial run's end)
# MODE    host (default): the RT and the hydro on the level, the host state
#         on every rank (amrex_state = 0) and the serial PM, so that the run
#         is bitwise the serial production run's continuation (the AMReX PM
#         agrees with the serial PM to round-off only); state: the
#         production configuration (amrex_state = 1, the three backends; the
#         state on the slabs), for the timings
#
#SBATCH --job-name=moments-validate
#SBATCH --account=astro
#SBATCH --partition=cpu-g2
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --time=04:00:00
#SBATCH --output=%x-%j.out
set -u
cd "${SLURM_SUBMIT_DIR}"
source tools/hyak/env.sh
radray_module_load ${HYAK_MODULES:-ompi/4.1.6-2}
radray_ucx_eager
SRC="$(readlink -f "${SRC:?set SRC to the serial run}")"
STEP="${STEP:?set STEP to the checkpoint step}"
NSTEPS="${NSTEPS:-1}"
OUT="${OUT:-$SRC/level_${SLURM_JOB_ID}}"
EXE="$(readlink -f "${EXE:-build-amrex/radhydro}")"
RANKS="${RANKS:-${SLURM_NTASKS:-1}}"
export OMP_NUM_THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-8}}"
export OMP_PLACES=cores OMP_PROC_BIND=close
ulimit -s unlimited
CK=$(printf "%s/run/restart/checkpoint_%04d.bin" "$SRC" "$STEP")
[ -e "$CK" ] || { echo "no checkpoint $CK"; exit 1; }
mkdir -p "$OUT"
case "${MODE:-host}" in
  host) backends="--set pm_backend=serial --set amrex_state=0" ;;
  state) backends="--set amrex_hydro_max_grid_size=0 --set pm_backend=amrex --set amrex_state=1" ;;
  *) echo "MODE must be host or state"; exit 1 ;;
esac
echo "$(date) MODE ${MODE:-host}: $SRC step $STEP + $NSTEPS on $RANKS ranks x $OMP_NUM_THREADS threads, ${SLURM_JOB_NUM_NODES:-1} nodes (${SLURM_JOB_NODELIST:-})"
/usr/bin/time -v srun --mpi=pmix --ntasks="$RANKS" --cpus-per-task="$OMP_NUM_THREADS" --cpu-bind=cores --kill-on-bad-exit=1 \
  --chdir="$OUT" "$EXE" "$SRC" "$OUT" --quiet --restart "$CK" \
  --set run_steps="$NSTEPS" --set restart_every=0 \
  --set rt_backend=amrex --set hydro_backend=amrex $backends --set amrex_slab=1 --set amrex_slab_planes="${PLANES:-0}" \
  --set amrex_moment_stats="${STATS:-exact}" ${EXTRA:-} > "$OUT.log" 2>&1
rc=$?
echo "$(date) radhydro exit $rc"
grep -E "radhydro:|peak RSS|AMReX moment timers|AMReX PM timers|AMReX hydro timers|Maximum resident|Elapsed \(wall" "$OUT.log"
[ $rc -eq 0 ] || exit $rc
[ "${COMPARE:-1}" = 0 ] && exit 0
# every file of the new steps against the serial run's
n=0
bad=0
for s in $(seq $((STEP + 1)) $((STEP + NSTEPS))); do
  d=$(printf "step_%04d" "$s")
  for f in "$OUT/$d"/*; do
    [ -e "$f" ] || continue
    g="$SRC/run/$d/$(basename "$f")"
    [ -e "$g" ] || { echo "only in the level run: $d/$(basename "$f")"; continue; }
    n=$((n + 1))
    cmp -s "$f" "$g" || { bad=$((bad + 1)); echo "DIFFERS: $d/$(basename "$f")"; }
  done
done
for f in "$OUT"/movie/frame_*; do
  [ -e "$f" ] || continue
  g="$SRC/run/movie/$(basename "$f")"
  [ -e "$g" ] || continue
  n=$((n + 1))
  cmp -s "$f" "$g" || { bad=$((bad + 1)); echo "DIFFERS: movie/$(basename "$f")"; }
done
echo "compared $n files with the serial run: $bad differ"
[ $bad -eq 0 ]
