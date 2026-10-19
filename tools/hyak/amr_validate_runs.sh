#!/bin/bash
# The refined production configuration at 128^3 against uniform runs
# (docs/stage_f.md section 3, "Validation": "a 128^3 production-
# configuration run with two levels against the uniform 128^3 and 256^3
# runs"). Makes three problem directories under ROOT with
# tools/production/make_production.py (10 Mpc/h, 0.1c, the abundance-matched
# sources of data/uv_lf/am_params_st_fof.txt, the AMReX backends):
#
#   u128    gas 128^3, pm_refine = PM_REFINE (default 4: the 512^3 mesh),
#           256^3 particles
#   amr128  the same ICs and parameters, plus amrex_max_level = LEVELS
#           (default 2) with amrex_refine_mass = REFINE_MASS (default 8:
#           a level-l cell refines above 8 * 8^l times the mean gas density);
#           the refined gas enters the PM on its own cells (pm_gas_levels =
#           1, docs/stage_f.md section 4; PM_GAS_LEVELS=0 for the level-0
#           coupling)
#   u256    gas 256^3, pm_refine = PM_REFINE/2 (the same mesh), the DM ICs
#           of u128 byte for byte (make_ics.py --phases-from-ndm 128), the
#           gas holding the modes the 128^3 grid cannot
#
# so the three runs share the dark matter and the PM mesh and differ only in
# the gas and RT resolution: amr128 should approach u256 where it refines.
# The defaults are the 512^3 target (gas 512^3, two levels, the 2048^3 mesh,
# 1024^3 particles) at a quarter of the size.
# Then submits tools/hyak/reion_production.slurm for each (one whole node
# each, 24 ranks x 8 threads, POST=0). Compare when all three are done:
#
#   python3 tools/analysis/amr_validate_compare.py ROOT
#
# usage (from the repository root on a klone login node, after
# tools/hyak/build_amrex.sh has built build-amrex/radhydro from main):
#   tools/hyak/amr_validate_runs.sh
# Environment: ROOT (default $RUNS/amr_val128), RUNS (default
# /gscratch/astro/$USER/radhydro/runs), LEVELS (2), REFINE_MASS (8),
# PM_REFINE (4; even), PM_GAS_LEVELS (1), F_ESC
# (0.2), COOL (cool_suppress_overdensity, default 200 = the original; 0: dense
# gas cools), SF (sf_overdensity, default 0 = no sink; the user's production
# choice is "a few thousand"), RUNS_TO_SUBMIT (default "u128 amr128 u256"),
# EXE (default build-amrex/radhydro), SBATCH_ARGS (default one whole node for
# 2 days), DRYRUN=1 (make the directories, print the sbatch commands).
# COOL and SF apply to all three runs, so the comparison stays at fixed
# physics. Rerunning the script reuses existing directories and ICs; the
# production script resumes a run from its checkpoint.
set -eu
RUNS="${RUNS:-/gscratch/astro/$USER/radhydro/runs}"
ROOT="${ROOT:-$RUNS/amr_val128}"
LEVELS="${LEVELS:-2}"
REFINE_MASS="${REFINE_MASS:-8}"
PM_REFINE="${PM_REFINE:-4}"
PM_GAS_LEVELS="${PM_GAS_LEVELS:-1}"
[ $((PM_REFINE % 2)) = 0 ] || { echo "PM_REFINE must be even (u256 uses PM_REFINE/2)"; exit 1; }
F_ESC="${F_ESC:-0.2}"
COOL="${COOL:-200}"
SF="${SF:-0}"
EXE="${EXE:-build-amrex/radhydro}"
RUNS_TO_SUBMIT="${RUNS_TO_SUBMIT:-u128 amr128 u256}"
SBATCH_ARGS="${SBATCH_ARGS:---exclusive --nodes=1 --ntasks=24 --cpus-per-task=8 --mem=1400G --time=2-00:00:00}"
[ -x "$EXE" ] || { echo "$EXE missing: run tools/hyak/build_amrex.sh first"; exit 1; }
mkdir -p "$ROOT"
common=(--box 10 --c-reduced 0.1 --am-params data/uv_lf/am_params_st_fof.txt --f-esc "$F_ESC"
        --backend amrex --cool-suppress-overdensity "$COOL" --sf-overdensity "$SF")
mp=(python3 tools/production/make_production.py)

# u128: the ICs (gas 128^3, DM 256^3 with the large scales of the 128^3 ICs)
[ -f "$ROOT/u128/params.txt" ] || "${mp[@]}" "$ROOT/u128" --N 128 --pm-refine "$PM_REFINE" --ndm 256 "${common[@]}"
# amr128: the same ICs
if [ ! -f "$ROOT/amr128/params.txt" ]; then
  "${mp[@]}" "$ROOT/amr128" --N 128 --pm-refine "$PM_REFINE" --ndm 256 --no-ics "${common[@]}" \
    --amrex-max-level "$LEVELS" --amrex-refine-mass "$REFINE_MASS" --pm-gas-levels "$PM_GAS_LEVELS"
fi
[ -e "$ROOT/amr128/ics" ] || ln -s ../u128/ics "$ROOT/amr128/ics"
# u256: the DM of u128 (same seed, Ndm 256, phases from the 128^3 draws)
if [ ! -f "$ROOT/u256/params.txt" ]; then
  "${mp[@]}" "$ROOT/u256" --N 256 --pm-refine $((PM_REFINE / 2)) --ndm 256 --no-ics "${common[@]}"
fi
if [ ! -f "$ROOT/u256/ics/u_gas.init" ]; then
  python3 tools/ics/make_ics.py "$ROOT/u256/ics" --kind zeldovich --Ngas 256 --Ndm 256 \
    --z-init 100 --box 10 --seed 7 --phases-from-ndm 128
fi
if cmp -s "$ROOT/u128/ics/xv_dm.init" "$ROOT/u256/ics/xv_dm.init"; then
  echo "u128 and u256 share the DM ICs (xv_dm.init identical)"
else
  echo "WARNING: xv_dm.init of u128 and u256 differ; the comparison mixes DM realisations"
fi

for r in $RUNS_TO_SUBMIT; do
  d="$ROOT/$r"
  [ -f "$d/params.txt" ] || { echo "skip $r: no params.txt"; continue; }
  if [ -e "$d/run/DONE" ]; then echo "$r: already done"; continue; fi
  cmd=(sbatch $SBATCH_ARGS --job-name="amrval_$r" --output="$d/slurm-%j.out"
       --export=ALL,PROB="$d",EXE="$EXE",POST=0 tools/hyak/reion_production.slurm)
  if [ "${DRYRUN:-0}" = 1 ]; then echo "${cmd[*]}"; else "${cmd[@]}"; fi
done
echo "when the runs are done: python3 tools/analysis/amr_validate_compare.py $ROOT"
