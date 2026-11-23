#!/bin/bash
# Phase 2 (a) of docs/absorbers_closures.md: the suite's configuration at
# 256^3 in 20 Mpc/h (the precheck box: 256^3 gas, 512^3 particles, seed 7,
# pm_refine = 2, the Sheth-Tormen table, f6 = 0.11, beta = 2.5, the
# continuous c-hat schedule, case B, cooling suppression 200, no sink) with
# the rays, M1 and the OTVET cone tensor, each with photon_consume_dense = 1
# (P11, the production default since 2026-10-28) and = 0 (the suite's
# original rule): six runs, one whole cpu-g2 node each (the precheck's rays
# took 1.6-2.9 h; the moment methods about half). Submit from the
# repository root on a login node:
#   tools/hyak/absorber_p11_runs.sh            # make the six directories and submit
#   DRYRUN=1 tools/hyak/absorber_p11_runs.sh   # print the sbatch commands only
# Environment: RUNS (default /gscratch/astro/$USER/radhydro/runs), ICS (the
# precheck ICs, default $RUNS/rt512cal256/ics), EXE (default
# $RUNS/absorbers_wu21_bin/radhydro, a build with P11), METHODS (default
# "rays m1 cone"), P11 (default "1 0"), F6 (0.11), AM_PARAMS, OVERWRITE=1
# (remove existing run directories of the same names), SBATCH_ARGS,
# THREADS. Run names: $RUNS/abs256_<method>_p11<0|1>.
set -eu
RUNS="${RUNS:-/gscratch/astro/$USER/radhydro/runs}"
ICS="${ICS:-$RUNS/rt512cal256/ics}"
EXE="${EXE:-$RUNS/absorbers_wu21_bin/radhydro}"
METHODS="${METHODS:-rays m1 cone}"
P11="${P11:-1 0}"
F6="${F6:-0.11}"
AM_PARAMS="${AM_PARAMS:-data/uv_lf/am_params_st_fof.txt}"
THREADS="${THREADS:-192}"
SBATCH_ARGS="${SBATCH_ARGS:---exclusive --cpus-per-task=192 --mem=300G --time=12:00:00}"
export EXE THREADS
[ -x "$EXE" ] || { echo "$EXE missing"; exit 1; }
[ -f "$ICS/u_gas.init" ] && [ -f "$ICS/xv_dm.init" ] || { echo "ICs missing at $ICS (tools/hyak/rt_suite512.sh precheck makes them)"; exit 1; }
grep -Eq '^Ngas = 256$' "$ICS/ics_info.txt" && grep -Eq '^Ndm = 512$' "$ICS/ics_info.txt" || { echo "$ICS is not the 256^3/512^3 box"; exit 1; }
[ -f "$AM_PARAMS" ] || { echo "$AM_PARAMS missing"; exit 1; }

submit() { if [ "${DRYRUN:-0}" = 1 ]; then echo "DRYRUN: $*"; else "$@"; fi; }

dirs=()
for m in $METHODS; do for p in $P11; do dirs+=("$RUNS/abs256_${m}_p11$p"); done; done
stop=0
for d in "${dirs[@]}"; do
  [ -e "$d" ] || continue
  if [ "${OVERWRITE:-0}" = 1 ]; then echo "removing $d"; [ "${DRYRUN:-0}" = 1 ] || rm -rf "$d"
  else echo "exists: $d"; stop=1; fi
done
[ "$stop" = 0 ] || { echo "set OVERWRITE=1 to remove them"; exit 1; }

for m in $METHODS; do
  for p in $P11; do
    d="$RUNS/abs256_${m}_p11$p"
    # the suite's params.txt (tools/hyak/rt_suite512.sh make_params) at the precheck size, P11 explicit
    python3 tools/production/make_production.py "$d" --N 256 --box 20 --pm-refine 2 --ndm 512 \
      --f-esc "$F6" --movie-dz 0.05 --no-ics --nray-max-factor 15 --am-params "$AM_PARAMS" \
      --f-esc-beta 2.5 --f-esc-z0 6 --xi-ion 25.38 --movie-slab 8 \
      --movie-fields fHI,temp,rhog,gammaHI,gammaHI_avg \
      --fields fHI,fHeI,fHeII,temp,rhog,gammaHI,gammaHI_avg --zdata1 "15 7 6.5 6.25 6 5.75 5.5 5.25 5" \
      --restart-every 10 --restart-keep 1 --data-files gas --cool-suppress-overdensity 200 \
      --recomb-case B --recomb-rad 0 --c-schedule-interp 1 --photon-consume-dense "$p" > /dev/null
    grep -Eq "^photon_consume_dense = $p\$" "$d/params.txt" || { echo "photon_consume_dense = $p missing in $d/params.txt"; exit 1; }
    {
      echo "# Phase 2 (a) of docs/absorbers_closures.md: $m with photon_consume_dense = $p (tools/hyak/absorber_p11_runs.sh)"
      case "$m" in
        rays) ;;
        m1) echo "rt_method = m1" ;;
        otvet) echo "rt_method = otvet" ;;
        cone) echo "rt_method = otvet"; echo "otvet_closure = cone" ;;
        split) echo "rt_method = otvet"; echo "otvet_split = 1" ;;
        screen) echo "rt_method = otvet"; echo "otvet_closure = screen" ;;
        *) echo "unknown method $m" >&2; exit 2 ;;
      esac
    } >> "$d/params.txt"
    ln -sfn "$ICS" "$d/ics"
    submit sbatch $SBATCH_ARGS --job-name="abs256_${m}_p11$p" --output="$d/slurm-%j.out" \
      --export=ALL,PROB="$d",THREADS="$THREADS",EXE="$EXE",POST=sbatch tools/hyak/reion_production.slurm
  done
done
