#!/bin/bash
# The C1 ladder of docs/zoom_turbulence.md section 3b on klone: the uniform
# small boxes of C26 (Cain et al. 2026) with the plane-parallel RT-domain
# mode (docs/plane_parallel.md), 0.25 Mpc/h at 64^3, 128^3, 256^3 and 512^3,
# z_re = 7, Gamma_-12 = 0.3, 8^3 RT domains lit from +x and +y, DM20's speed
# of light schedule (0.1 c until x_HI,V < 1 per cent, then 0.01 c), cooling on
# everywhere, no preheating unless PREHEAT is set (the pre-heating floor in K
# from z = 15, C26's 100 or 1000), Ndm = Ngas. Also the T4 coupling-interval
# test at 128^3 (tools/analysis/t4_coupling.py; Zel'dovich ICs as before).
#
# Initial conditions (ICS):
#   music (default): C26's set-up. MUSIC in uniform mode (tools/ics/music/
#     c2_uniform.conf at level log2(N), seed[10] = MUSIC_SEED, so every rung
#     and the 128^3 linear-growth box of tools/hyak/linear_growth.slurm share
#     their modes) with CAMB's separate CDM and baryon transfer functions at
#     z_init = 300, 1LPT, the C26 cosmology (Omega_m 0.305, Omega_b 0.048,
#     h 0.68, sigma_8 0.82, n_s 0.9667), converted by tools/ics/music_convert.py;
#     THERMAL = recfast (default): thermal_init = recfast with the RECFAST
#     table of that cosmology (recfast_table, make_recfast_table.py); floor:
#     the original's floor history. Linear growth from these ICs is checked in
#     docs/thermal_history.md section 6.
#   zeldovich: the generator of make_ics.py from z = 100 in the original
#     cosmology with the floor history (the ladder before 2026-10-30, and the
#     local 64^3 rung of docs/STATUS.md).
# thermal_init = recfast acts on the host gas grid only and is refused with
# the backend's amrex_state = 1, so with recfast every rung runs the serial
# build on the whole node (no ray pool with the default sweep transport; about
# 100 GB at 512^3); with THERMAL = floor the 512^3 rung uses the backend as
# before. BACKEND_512 = serial|amrex overrides.
#
# usage (from the repository root, after `source tools/hyak/env.sh` and a
# build in build/ (serial) and build-amrex/ (the backend, for 512^3 with
# THERMAL = floor)):
#   tools/hyak/c1_ladder.sh make [N ...]    # the problems under $RUNS (MUSIC and CAMB
#                                           # run here: use a compute node, or MAKE_SBATCH=1)
#   tools/hyak/c1_ladder.sh run [N ...]     # sbatch one job per N (default: all four)
#   tools/hyak/c1_ladder.sh post [N ...]    # E(k), the front profiles, the figures
#   tools/hyak/c1_ladder.sh t4 [make|run|post]   # the T4 variants at $T4_N (128)
# Environment: RUNS (default /gscratch/astro/$USER/radhydro/runs/c1), ICS
# (music), THERMAL (recfast with music, floor with zeldovich), ZINIT (300 with
# music, 100 with zeldovich), MUSIC_SEED (2026), SEED (7: the Zel'dovich
# generator's), the cosmology of the MUSIC ICs OMEGAM (0.305), OMEGAB (0.048),
# HUBBLE0 (0.68), SIGMA8 (0.82), NS (0.9667), MUSIC ($RADRAY_SW/music/MUSIC),
# GAMMA12 (0.3), ZRE (7), NDOM (8), PREHEAT (0), EXE (build/radhydro), EXE_AMREX
# (build-amrex/radhydro), THREADS (192), RANKS_512 (12: MPI ranks of a 512^3
# backend job, 16 threads each), BACKEND_512, SBATCH_ARGS, MAKE_SBATCH=1 (make
# submits itself as a job), DRYRUN=1 (print the commands).
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
ICS=${ICS:-music}
case "$ICS" in
  music) THERMAL=${THERMAL:-recfast}; ZINIT=${ZINIT:-300} ;;
  zeldovich) THERMAL=${THERMAL:-floor}; ZINIT=${ZINIT:-100} ;;
  *) echo "ICS must be music or zeldovich" >&2; exit 2 ;;
esac
case "$THERMAL" in recfast|floor) ;; *) echo "THERMAL must be recfast or floor" >&2; exit 2 ;; esac
[ "$ICS" = zeldovich ] && [ "$THERMAL" = recfast ] && { echo "THERMAL = recfast needs ICS = music (the C26 cosmology's table)" >&2; exit 2; }
MUSIC_SEED=${MUSIC_SEED:-2026}
OMEGAM=${OMEGAM:-0.305} OMEGAB=${OMEGAB:-0.048} HUBBLE0=${HUBBLE0:-0.68} SIGMA8=${SIGMA8:-0.82} NS=${NS:-0.9667}
MUSIC=${MUSIC:-${RADRAY_SW:-}/music/MUSIC}
if [ "$THERMAL" = recfast ]; then BACKEND_512=${BACKEND_512:-serial}; else BACKEND_512=${BACKEND_512:-amrex}; fi
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

backend_of() { if [ "$1" -ge 512 ]; then echo "$BACKEND_512"; else echo serial; fi; }

make_problem() {
  local N=$1 dir=$RUNS/c1_N$N backend
  backend=$(backend_of "$N")
  mkdir -p "$RUNS"
  local common=(--N "$N" --ndm "$N" --box 0.25 --z-init "$ZINIT" --seed "$SEED"
    --plane-parallel --plane-ndom "$NDOM" --plane-z-re "$ZRE" --plane-gamma12 "$GAMMA12"
    --zdata1 "$ZDATA" --fields fHI,temp,rhog,u_gas --data-files gas --movie-dz 0
    --restart-every 20 --restart-keep 1 --nray-max-factor 8 --backend "$backend"
    --T-min "$PREHEAT" --z-floor 15)
  if [ "$ICS" = zeldovich ]; then
    run_cmd python3 "$here/tools/production/make_production.py" "$dir" "${common[@]}"
    return
  fi
  # MUSIC + CAMB (as tools/hyak/linear_growth.slurm): the transfer functions
  # at z_init, c2_uniform.conf at level log2(N), the code's ICs from MUSIC's
  # output, the RECFAST table of the cosmology
  local level cos=(--omegam "$OMEGAM" --omegab "$OMEGAB" --hubble0 "$HUBBLE0")
  level=$(python3 -c "import math; print(int(round(math.log2($N))))")
  local md=$RUNS/music tf=$RUNS/music/camb_z$ZINIT.dat conf=$RUNS/music/c1_$N.conf
  run_cmd mkdir -p "$md"
  [ -f "$tf" ] || run_cmd python3 "$here/tools/ics/music/camb_transfer.py" "$tf" --zstart "$ZINIT" \
    --omega-m "$OMEGAM" --omega-b "$OMEGAB" --h "$HUBBLE0" --ns "$NS" --sigma8 "$SIGMA8" --kmax 5000
  if [ ! -f "$md/c1_$N.hdf5" ]; then
    local sedargs=(-e "s/^levelmin *=.*/levelmin    = $level/" -e "s/^levelmin_TF *=.*/levelmin_TF = $level/"
      -e "s/^levelmax *=.*/levelmax    = $level/" -e "s/^boxlength *=.*/boxlength   = 0.25/"
      -e "s/^zstart *=.*/zstart      = $ZINIT/" -e "s/^seed\[10\] *=.*/seed[10] = $MUSIC_SEED/"
      -e "s/^Omega_m *=.*/Omega_m       = $OMEGAM/"
      -e "s/^Omega_L *=.*/Omega_L       = $(python3 -c "print(round(1 - $OMEGAM, 10))")/"
      -e "s/^Omega_b *=.*/Omega_b       = $OMEGAB/"
      -e "s/^H0 *=.*/H0            = $(python3 -c "print(round(100 * $HUBBLE0, 10))")/"
      -e "s/^sigma_8 *=.*/sigma_8       = $SIGMA8/" -e "s/^n_s *=.*/n_s           = $NS/"
      -e "s|^transfer_file *=.*|transfer_file = $tf|" -e "s/^filename *=.*/filename = c1_$N.hdf5/")
    if [ "${DRYRUN:-0}" = 1 ]; then
      echo "+ sed ${sedargs[*]} $here/tools/ics/music/c2_uniform.conf > $conf"
      echo "+ (cd $md && $MUSIC c1_$N.conf)"
    else
      sed "${sedargs[@]}" "$here/tools/ics/music/c2_uniform.conf" > "$conf"
      (cd "$md" && "$MUSIC" "c1_$N.conf" > "c1_$N.music.log" 2>&1)
    fi
  fi
  [ -f "$dir/params.txt" ] || run_cmd python3 "$here/tools/production/make_production.py" "$dir" "${common[@]}" \
    "${cos[@]}" --sigma8 "$SIGMA8" --ns "$NS" --no-ics
  [ -f "$dir/ics/xv_dm.init" ] || run_cmd python3 "$here/tools/ics/music_convert.py" "$md/c1_$N.hdf5" "$dir/ics" \
    --conf "$conf" "${cos[@]}" --sigma8 "$SIGMA8" --ns "$NS"
  if [ "$THERMAL" = recfast ] && [ ! -f "$dir/recfast.txt" ]; then
    run_cmd python3 "$here/tools/thermal/make_recfast_table.py" --out "$dir/recfast.txt" "${cos[@]}"
    if [ "${DRYRUN:-0}" = 1 ]; then
      echo "+ append thermal_init = recfast, recfast_table = $dir/recfast.txt to $dir/params.txt"
    else
      {
        echo "# the RECFAST thermal history of this cosmology from z_init (docs/thermal_history.md section 4)"
        echo "thermal_init = recfast"
        echo "recfast_table = $dir/recfast.txt"
      } >> "$dir/params.txt"
    fi
  fi
}

case "$phase" in
  make)
    if [ "${MAKE_SBATCH:-0}" = 1 ]; then
      run_cmd sbatch --account=astro --partition=cpu-g2 --cpus-per-task=32 --mem=200G --time=06:00:00 \
        --job-name=c1_make --export=ALL,MAKE_SBATCH=0 \
        --wrap "cd $here && source tools/hyak/env.sh && tools/hyak/c1_ladder.sh make $sizes"
    else
      for N in $sizes; do make_problem "$N"; done
    fi ;;
  run)
    for N in $sizes; do
      dir=$RUNS/c1_N$N
      if [ "$ICS" = music ] && [ ! -f "$dir/ics/xv_dm.init" ] && [ "${DRYRUN:-0}" != 1 ]; then
        echo "c1_N$N: no MUSIC ICs yet; run '$0 make $N' first (MUSIC and CAMB, on a compute node or with MAKE_SBATCH=1)" >&2
        exit 2
      fi
      [ -f "$dir/params.txt" ] || make_problem "$N"
      if [ "$(backend_of "$N")" = amrex ]; then
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
