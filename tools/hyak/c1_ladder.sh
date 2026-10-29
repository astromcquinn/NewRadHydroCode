#!/bin/bash
# The C1 ladder of docs/zoom_turbulence.md section 3b on klone: the uniform
# small boxes of C26 (Cain et al. 2026) with the plane-parallel RT-domain
# mode (docs/plane_parallel.md), 0.25 Mpc/h at 64^3, 128^3, 256^3 and 512^3,
# z_re = 7, Gamma_-12 = 0.3, 8^3 RT domains lit from +x and +y, DM20's speed
# of light schedule (0.1 c until x_HI,V < 1 per cent, then 0.01 c), cooling on
# everywhere, no preheating unless PREHEAT is set (the pre-heating floor in K
# from z = 15, C26's 100 or 1000); Zel'dovich initial conditions of the existing
# generator from z = 100 (until the MUSIC/CAMB ICs of WP1 land), Ndm = Ngas.
# Also the T4 coupling-interval test at 128^3 (tools/analysis/t4_coupling.py).
#
# usage (from the repository root, after `source tools/hyak/env.sh` and a
# build in build/ (serial) and build-amrex/ (the backend, for 512^3)):
#   tools/hyak/c1_ladder.sh make            # the problems under $RUNS
#   tools/hyak/c1_ladder.sh run [N ...]     # sbatch one job per N (default: all four)
#   tools/hyak/c1_ladder.sh post [N ...]    # E(k), the front profiles, the figures
#   tools/hyak/c1_ladder.sh t4 [make|run|post]   # the T4 variants at $T4_N (128)
# Environment: RUNS (default /gscratch/astro/$USER/radhydro/runs/c1), SEED (7),
# GAMMA12 (0.3), ZRE (7), NDOM (8), PREHEAT (0), EXE (build/radhydro), EXE_AMREX
# (build-amrex/radhydro), THREADS (192), RANKS_512 (12: MPI ranks of the
# 512^3 job, 16 threads each), SBATCH_ARGS, DRYRUN=1 (print the commands).
#
# Outputs per run (RUNS/c1_N<N>/run): the gas state u_gas, fHI, temp and
# rhog at z = 7 (the turn-on), 6.9, 6.8, 6.6, 6.4, 6.2, 6, 5.5, 5, 4.5 and 4
# (C26: E(k) at z = 7.2 ... 4), stats.txt per step with the plane-parallel
# photon bookkeeping. post writes RUNS/post/<N>/velocity_power.png (the
# E(k) ladder for C26 Fig. 3 bottom) and plane_front.png.
set -euo pipefail
here=$(cd "$(dirname "$0")/../.." && pwd)
RUNS=${RUNS:-/gscratch/astro/$USER/radhydro/runs/c1}
SEED=${SEED:-7}
GAMMA12=${GAMMA12:-0.3}
ZRE=${ZRE:-7}
NDOM=${NDOM:-8}
PREHEAT=${PREHEAT:-0}      # the pre-heating floor in K from z = 15 (C26: 100 or 1000; 0 = none)
EXE=${EXE:-$here/build/radhydro}
EXE_AMREX=${EXE_AMREX:-$here/build-amrex/radhydro}
THREADS=${THREADS:-192}
RANKS_512=${RANKS_512:-12}
SBATCH_ARGS=${SBATCH_ARGS:---exclusive --cpus-per-task=192 --mem=256G --time=2-00:00:00}
T4_N=${T4_N:-128}
ZDATA="7 6.9 6.8 6.6 6.4 6.2 6 5.5 5 4.5 4"
phase=${1:-}
shift || true
sizes=${*:-64 128 256 512}

run_cmd() { if [ "${DRYRUN:-0}" = 1 ]; then echo "+ $*"; else "$@"; fi; }

make_problem() {
  local N=$1 dir=$RUNS/c1_N$N backend=serial
  [ "$N" -ge 512 ] && backend=amrex
  mkdir -p "$RUNS"
  run_cmd python3 "$here/tools/production/make_production.py" "$dir" \
    --N "$N" --ndm "$N" --box 0.25 --z-init 100 --seed "$SEED" \
    --plane-parallel --plane-ndom "$NDOM" --plane-z-re "$ZRE" --plane-gamma12 "$GAMMA12" \
    --zdata1 "$ZDATA" --fields fHI,temp,rhog,u_gas --data-files gas --movie-dz 0 \
    --restart-every 20 --restart-keep 1 --nray-max-factor 8 --backend $backend \
    --preheat-floor "$PREHEAT"
}

case "$phase" in
  make)
    for N in $sizes; do make_problem "$N"; done ;;
  run)
    for N in $sizes; do
      dir=$RUNS/c1_N$N
      [ -f "$dir/params.txt" ] || make_problem "$N"
      if [ "$N" -ge 512 ]; then
        run_cmd sbatch $SBATCH_ARGS --job-name=c1_N$N \
          --export=ALL,PROB=$dir,EXE=$EXE_AMREX,RANKS=$RANKS_512,THREADS=$((192 / RANKS_512)),POST=0 \
          "$here/tools/hyak/reion_production.slurm"
      else
        run_cmd sbatch $SBATCH_ARGS --job-name=c1_N$N \
          --export=ALL,PROB=$dir,EXE=$EXE,THREADS=$THREADS,POST=0 \
          "$here/tools/hyak/reion_production.slurm"
      fi
    done ;;
  post)
    mkdir -p "$RUNS/post"
    specs=()
    for N in $sizes; do specs+=("$RUNS/c1_N$N/run:N$N"); done
    run_cmd python3 "$here/tools/analysis/velocity_power.py" "${specs[@]}" --outdir "$RUNS/post" \
      --z 6.9 6.6 6.2 5.5 5 4 --json "$RUNS/post/ek.json"
    run_cmd python3 "$here/tools/analysis/plane_front.py" "${specs[@]}" --outdir "$RUNS/post" \
      --json "$RUNS/post/front.json" ;;
  t4)
    sub=${1:-run}
    root=$RUNS/t4_N$T4_N
    case "$sub" in
      make) run_cmd python3 "$here/tools/analysis/t4_coupling.py" make "$root" --N "$T4_N" --seed "$SEED" --gamma12 "$GAMMA12" --z-re "$ZRE" --ndom "$NDOM" ;;
      run)
        [ -f "$root/prob/params.txt" ] || run_cmd python3 "$here/tools/analysis/t4_coupling.py" make "$root" --N "$T4_N" --seed "$SEED" --gamma12 "$GAMMA12" --z-re "$ZRE" --ndom "$NDOM"
        for v in base dt100 dt10 dt1 dt0.1 sub; do
          run_cmd sbatch $SBATCH_ARGS --job-name=t4_$v --wrap \
            "source $here/tools/hyak/env.sh && python3 $here/tools/analysis/t4_coupling.py run $root --exe $EXE --threads $THREADS --variants $v"
        done ;;
      post) run_cmd python3 "$here/tools/analysis/t4_coupling.py" analyse "$root" ;;
      *) echo "t4: make | run | post" >&2; exit 2 ;;
    esac ;;
  *)
    echo "usage: $0 make|run|post [N ...] | t4 make|run|post" >&2
    exit 2 ;;
esac
