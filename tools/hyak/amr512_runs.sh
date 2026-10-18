#!/bin/bash
# The first refined production run and its unrefined control (docs/stage_f.md
# sections 3-4): 40 Mpc/h, gas 512^3, 1024^3 particles on the 2048^3 PM mesh
# (pm_refine = 4), 0.1c, the abundance-matched sources with the falling
# escape fraction of the f_esc thread (f_esc(z) = 0.135 ((1+z)/7)^2.5,
# docs/PLAN.md on branch fesc-evolution), on the AMReX backends (case B: the
# backend has no recombination radiation):
#
#   amr512   two levels (amrex_max_level = 2, amrex_refine_mass = 8: a
#            level-l cell refines above 8 * 8^l times the mean gas density),
#            the refined gas in the PM (pm_gas_levels = 1)
#   u512     the same parameters without refinement: the control, so the
#            pair isolates what the refinement does to the IGM (Gamma_HI,
#            lambda_mfp, the history)
#
# Both use the ICs of hmf512_40mpc_r2 (gas 512^3, DM 1024^3 with the large
# scales of the 512^3 production ICs, seed 7; docs/pm_refine.md) through a
# symlink, ray_deposit = exposure (the Stage F decision for rays crossing
# 2^l cells; on both runs), the snapshots of the closure runs, FoF every
# 30 Myr. Then submits tools/hyak/reion_production.slurm for each on 5 nodes
# (24 ranks x 8 threads per node; the second waits for the allocation).
#
# usage (from the repository root on a klone login node, after
# tools/hyak/build_amrex.sh built build-amrex/radhydro from this branch):
#   tools/hyak/amr512_runs.sh
# Environment: RUNS (default /gscratch/astro/$USER/radhydro/runs), ICS
# (default $RUNS/hmf512_40mpc_r2/ics), PREFIX (default $RUNS/fesc512_b2.5:
# the runs are PREFIX_amr and PREFIX_uniform), RUNS_TO_SUBMIT ("amr uniform"),
# F6 (0.135), BETA (2.5), LEVELS (2), REFINE_MASS (8), COOL
# (cool_suppress_overdensity, default 0: dense gas cools, the Stage F
# decision), SF (sf_overdensity, default 3000: "an aggressive sink above a
# few thousand", the Stage F decision; 0 = off), DATA_FILES (default nodm:
# Gas (u, v, w = rho_HI..., e) and RT/rt_grid at the snapshots, no particles;
# norays adds the 26 GB DM/xv_dm per snapshot, gas keeps only Gas/u_gas),
# FOF_DT (30), EXE (build-amrex/radhydro),
# NODES (5), SBATCH_ARGS, DRYRUN=1 (make the directories, print the sbatch
# commands). COOL, SF and everything else apply to both runs.
set -eu
RUNS="${RUNS:-/gscratch/astro/$USER/radhydro/runs}"
ICS="${ICS:-$RUNS/hmf512_40mpc_r2/ics}"
PREFIX="${PREFIX:-$RUNS/fesc512_b2.5}"
RUNS_TO_SUBMIT="${RUNS_TO_SUBMIT:-amr uniform}"
F6="${F6:-0.135}"
BETA="${BETA:-2.5}"
LEVELS="${LEVELS:-2}"
REFINE_MASS="${REFINE_MASS:-8}"
COOL="${COOL:-0}"
SF="${SF:-3000}"
DATA_FILES="${DATA_FILES:-nodm}"
FOF_DT="${FOF_DT:-30}"
EXE="${EXE:-build-amrex/radhydro}"
NODES="${NODES:-5}"
SBATCH_ARGS="${SBATCH_ARGS:---exclusive --nodes=$NODES --ntasks=$((24 * NODES)) --cpus-per-task=8 --mem=1400G --time=4-00:00:00}"
[ -x "$EXE" ] || { echo "$EXE missing: run tools/hyak/build_amrex.sh first"; exit 1; }
# the ICs: 512^3 gas, 1024^3 particles, 40 Mpc/h
info="$ICS/ics_info.txt"
if [ ! -f "$ICS/u_gas.init" ] || [ ! -f "$ICS/xv_dm.init" ] || [ ! -f "$info" ] ||
   ! grep -Eq '^Ngas = 512$' "$info" || ! grep -Eq '^Ndm = 1024$' "$info" || ! grep -Eq '^box = 40(\.0)?$' "$info"; then
  echo "no 512^3 / 1024^3 / 40 Mpc/h ICs at $ICS (ics_info.txt: Ngas = 512, Ndm = 1024, box = 40)."
  echo "Point ICS at them, or make them on a compute node (~300 GB, about an hour):"
  echo "  python3 tools/ics/make_ics.py <dir> --kind zeldovich --Ngas 512 --Ndm 1024 --z-init 100 --box 40 --seed 7 --phases-from-ndm 512"
  exit 1
fi
echo "ICs: $ICS"; grep -E '^(command|Ngas|Ndm|box|seed|z_init) =' "$info"
common=(--N 512 --box 40 --pm-refine 4 --ndm 1024 --no-ics --backend amrex --c-reduced 0.1
        --nray-max-factor 15 --am-params data/uv_lf/am_params_st_fof.txt
        --f-esc "$F6" --f-esc-beta "$BETA" --f-esc-z0 6 --fof-dt-myr "$FOF_DT"
        --zdata1 "15 7 6.5 6.25 6 5.75 5.5 5.25 5" --ray-deposit exposure
        --cool-suppress-overdensity "$COOL" --sf-overdensity "$SF" --data-files "$DATA_FILES")
mp=(python3 tools/production/make_production.py)
for r in amr uniform; do
  d="${PREFIX}_$r"
  if [ ! -f "$d/params.txt" ]; then
    if [ "$r" = amr ]; then
      "${mp[@]}" "$d" "${common[@]}" --amrex-max-level "$LEVELS" --amrex-refine-mass "$REFINE_MASS"
    else
      "${mp[@]}" "$d" "${common[@]}"
    fi
  fi
  [ -e "$d/ics" ] || ln -s "$ICS" "$d/ics"
done
for r in $RUNS_TO_SUBMIT; do
  d="${PREFIX}_$r"
  [ -f "$d/params.txt" ] || { echo "skip $r: no params.txt"; continue; }
  if [ -e "$d/run/DONE" ]; then echo "$r: already done"; continue; fi
  cmd=(sbatch $SBATCH_ARGS --job-name="fesc512_$r" --output="$d/slurm-%j.out"
       --export=ALL,PROB="$d",EXE="$EXE",POST=sbatch tools/hyak/reion_production.slurm)
  if [ "${DRYRUN:-0}" = 1 ]; then echo "${cmd[*]}"; else "${cmd[@]}"; fi
done
echo "progress: tail -f ${PREFIX}_amr/run.log; per step: ${PREFIX}_amr/run/timing.txt"
