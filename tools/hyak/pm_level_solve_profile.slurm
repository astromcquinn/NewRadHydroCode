#!/bin/bash
# The level Poisson solve (pm_level_solve = 1, Stage G, docs/stage_g.md) at
# 128^3 with two refined levels on klone: STEPS hydro steps of a gravity +
# hydro problem (no sources) from Zel'dovich ICs, the three backends with
# the state on the levels, for the three gravity couplings of the refined
# levels:
#   level_solve   pm_level_solve = 1 (pm_refine = 1: level 0 on the 128^3
#                 mesh by the FFT, level 1 and 2 on their own meshes by MLMG)
#   gas_levels    pm_gas_levels = 1, pm_refine = 4 (the uniform 512^3 mesh fed
#                 by the levels: Stage F's coupling, the comparison)
#   uniform       no refinement, pm_refine = 1 (the base cost)
# The per-step timers ("step NNNN timers:" of <run>.log: the PM backend's
# phases, the MLMG iterations and residual per level) and the wall-clock
# per step (timing.txt) are the result; tools/hyak/pm_level_solve_table.py
# prints the table. One node by default (24 ranks x 8 threads); more nodes
# with --nodes/--ntasks.
#
#   sbatch tools/hyak/pm_level_solve_profile.slurm
#   sbatch --nodes=2 --ntasks=48 tools/hyak/pm_level_solve_profile.slurm
#
# Environment:
#   OUT       output root (default $RUNS/pm_level_solve_128, RUNS from env.sh)
#   N         base grid (default 128); NDM particles per dimension (default N)
#   BOX       Mpc/h (default 10); ZINIT 100; ZEND 6 (zdata1)
#   REGION    the refined region in level-0 cells, i0,j0,k0,i1,j1,k1
#             (default the central eighth: 32,32,32,95,95,95 for N = 128:
#             level 1 128^3 and level 2 256^3 cells, 17 M cells on level 2,
#             tens of GB for the hydro hierarchy; a 15 GB VM needs the
#             central quarter width, 48,48,48,79,79,79)
#   STEPS     hydro steps (default 10)
#   CONFIGS   "level_solve gas_levels uniform" (default all three)
#   MLMG_ITERS, MLMG_TOL, MLMG_BOTTOM   pm_mlmg_iters (20), pm_mlmg_tol (1e-10),
#             pm_mlmg_bottom_sweeps (8; a rectangular level gets at least its bottom width squared)
#   EXTRA     extra radhydro arguments for every configuration, e.g.
#             EXTRA="--set pm_mlmg_verbose=1" (MLMG's residuals and its own
#             Solve / Iter / Bottom timers per level solve, in <config>.log)
#   RANKS, THREADS  (default the job's --ntasks, 8)
#   EXE       build with the backends (default build-amrex/radhydro, built
#             against the amrex-26.10-mlmg install of build_amrex.sh)
#SBATCH --job-name=pm-level-solve
#SBATCH --account=astro
#SBATCH --partition=cpu-g2
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=8
#SBATCH --mem=300G
#SBATCH --time=04:00:00
#SBATCH --exclusive
#SBATCH --output=%x-%j.out
set -u
cd "${SLURM_SUBMIT_DIR}"
source tools/hyak/env.sh
radray_module_load ${HYAK_MODULES:-ompi/4.1.6-2}
radray_ucx_eager
RUNS="${RUNS:-/gscratch/astro/$USER/radhydro/runs}"
OUT="${OUT:-$RUNS/pm_level_solve_128}"
N="${N:-128}"
NDM="${NDM:-$N}"
BOX="${BOX:-10}"
ZINIT="${ZINIT:-100}"
ZEND="${ZEND:-6}"
q=$((N / 4))
REGION="${REGION:-$q,$q,$q,$((3 * q - 1)),$((3 * q - 1)),$((3 * q - 1))}"
STEPS="${STEPS:-10}"
CONFIGS="${CONFIGS:-level_solve gas_levels uniform}"
MLMG_ITERS="${MLMG_ITERS:-20}"
MLMG_TOL="${MLMG_TOL:-1e-10}"
MLMG_BOTTOM="${MLMG_BOTTOM:-8}"
RANKS="${RANKS:-${SLURM_NTASKS:-24}}"
THREADS="${THREADS:-8}"
EXE="${EXE:-$PWD/build-amrex/radhydro}"
[ -x "$EXE" ] || { echo "no executable $EXE"; exit 1; }
mkdir -p "$OUT"
echo "host $(hostname), nodes ${SLURM_JOB_NODELIST:-}, job ${SLURM_JOB_ID:-}, $(date)"
echo "EXE $EXE ($(git rev-parse --short HEAD)), $RANKS ranks x $THREADS threads, N $N, NDM $NDM, region $REGION, steps $STEPS, configs $CONFIGS"

# the problem: Zel'dovich ICs, the gas following the DM, no sources
prob="$OUT/prob"
if [ ! -f "$prob/ics/xv_dm.init" ]; then
  python3 tools/ics/make_ics.py "$prob/ics" --kind zeldovich --Ngas "$N" --Ndm "$NDM" --z-init "$ZINIT" \
    --box "$BOX" --seed 7 --gas zeldovich
fi
cat > "$prob/params.txt" <<PARAMS
# pm_level_solve profile (tools/hyak/pm_level_solve_profile.slurm)
N = $N
Ndm = $NDM
box = $BOX
nray_max_factor = 4
z_init = $ZINIT
zdata1 = $ZEND
ics = ics
source_mode = legacy
dir_halo = Halo
fof_min_npart = 32
legacy_rt = 0
legacy_physics = 0
PARAMS

common=(--quiet --set run_steps="$STEPS" --set restart_every=0 --set movie_dz=0 --set fields_every=0
        --set step_timers=1 --set hydro_scheme=unsplit --set rt_backend=amrex --set hydro_backend=amrex
        --set pm_backend=amrex --set amrex_state=1 --set amrex_max_grid_size=32 --set amrex_hydro_max_grid_size=0
        --set amrex_fine_max_grid_size=64 --set deterministic_parallel=0 --set tvir_min=1e30)
levels=(--set amrex_max_level=2 --set amrex_refine_mass=1e9 --set amrex_refine_box="$REGION" --set amrex_n_error_buf=0)

for cfg in $CONFIGS; do
  case "$cfg" in
    level_solve) extra=("${levels[@]}" --set pm_level_solve=1 --set pm_refine=1 --set pm_mlmg_iters="$MLMG_ITERS" --set pm_mlmg_tol="$MLMG_TOL" --set pm_mlmg_bottom_sweeps="$MLMG_BOTTOM") ;;
    gas_levels)  extra=("${levels[@]}" --set pm_level_solve=0 --set pm_gas_levels=1 --set pm_refine=4) ;;
    uniform)     extra=(--set amrex_max_level=0 --set pm_refine=1) ;;
    *) echo "unknown config $cfg"; exit 1 ;;
  esac
  run="$OUT/$cfg"
  rm -rf "$run"
  echo "== $cfg: $(date)"
  export OMP_NUM_THREADS="$THREADS"
  # shellcheck disable=SC2086
  srun --mpi=pmix --ntasks="$RANKS" --cpus-per-task="$THREADS" "$EXE" "$prob" "$run" "${common[@]}" "${extra[@]}" ${EXTRA:-} > "$run.log" 2>&1
  echo "exit $?"
  grep "step .* timers" "$run.log" | tail -1
done
python3 tools/hyak/pm_level_solve_table.py "$OUT" $CONFIGS || true
echo "done $(date)"
