#!/bin/bash
# Linear growth of dark matter and gas from MUSIC ICs with CAMB's separate
# transfer functions (docs/thermal_history.md, "Linear growth from MUSIC
# ICs"): a uniform box, hydro and gravity only, thermal_init = recfast, from
# z_init to the zdata1 outputs; tools/analysis/linear_growth.py compares
# P_cdm(k) and P_b(k) with CAMB's linear theory.
#
# usage (from the repository root, after a serial build in build/):
#   sbatch --export=ALL,ROOT=/gscratch/astro/$USER/radhydro/runs/linear_growth tools/hyak/linear_growth.slurm
# Environment: ROOT (required), N (128), BOX (0.25 Mpc/h), ZINIT (300),
# ZOUT ("100 50 30"), SEED (2026, given at MUSIC level 10 as c2_uniform.conf,
# so the box shares its modes with the C1 ladder), the cosmology OMEGAM
# (0.305), OMEGAB (0.048), HUBBLE0 (0.68), SIGMA8 (0.82), NS (0.9667) (the
# planck_c26 preset of tools/ics/pmcosmo.py, C26), EXE (build/radhydro),
# MUSIC ($RADRAY_SW/music/MUSIC), EXTRA (more radhydro --set arguments).
# Steps already done (their outputs exist) are skipped.
#SBATCH --job-name=linear-growth
#SBATCH --account=astro
#SBATCH --partition=cpu-g2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=%x-%j.out
set -euo pipefail
here=$(pwd)
[ -f "$here/tools/hyak/env.sh" ] || { echo "run from the repository root" >&2; exit 2; }
source "$here/tools/hyak/env.sh" >/dev/null
: "${ROOT:?ROOT is required}"
N=${N:-128}
BOX=${BOX:-0.25}
ZINIT=${ZINIT:-300}
ZOUT=${ZOUT:-100 50 30}
SEED=${SEED:-2026}
OMEGAM=${OMEGAM:-0.305} OMEGAB=${OMEGAB:-0.048} HUBBLE0=${HUBBLE0:-0.68} SIGMA8=${SIGMA8:-0.82} NS=${NS:-0.9667}
EXE=$(readlink -f "${EXE:-build/radhydro}")
MUSIC=${MUSIC:-$RADRAY_SW/music/MUSIC}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_PLACES=cores OMP_PROC_BIND=close
level=$(python3 -c "import math; print(int(round(math.log2($N))))")
cos=(--omegam "$OMEGAM" --omegab "$OMEGAB" --hubble0 "$HUBBLE0")
mkdir -p "$ROOT/music"
echo "$(date) host $(hostname) $(git -C "$here" rev-parse --short HEAD) N=$N box=$BOX z_init=$ZINIT threads=$OMP_NUM_THREADS"

# 1. CAMB transfer functions at z_init (MUSIC's camb_file)
tf=$ROOT/music/camb_z$ZINIT.dat
[ -f "$tf" ] || python3 "$here/tools/ics/music/camb_transfer.py" "$tf" --zstart "$ZINIT" \
  --omega-m "$OMEGAM" --omega-b "$OMEGAB" --h "$HUBBLE0" --ns "$NS" --sigma8 "$SIGMA8" --kmax 5000

# 2. MUSIC: c2_uniform.conf at level log2(N)
conf=$ROOT/music/lg_$N.conf
if [ ! -f "$ROOT/music/lg_$N.hdf5" ]; then
  sed -e "s/^levelmin *=.*/levelmin    = $level/" -e "s/^levelmin_TF *=.*/levelmin_TF = $level/" \
      -e "s/^levelmax *=.*/levelmax    = $level/" -e "s/^boxlength *=.*/boxlength   = $BOX/" \
      -e "s/^zstart *=.*/zstart      = $ZINIT/" -e "s/^seed\[10\] *=.*/seed[10] = $SEED/" \
      -e "s/^Omega_m *=.*/Omega_m       = $OMEGAM/" \
      -e "s/^Omega_L *=.*/Omega_L       = $(python3 -c "print(round(1 - $OMEGAM, 10))")/" \
      -e "s/^Omega_b *=.*/Omega_b       = $OMEGAB/" \
      -e "s/^H0 *=.*/H0            = $(python3 -c "print(round(100 * $HUBBLE0, 10))")/" \
      -e "s/^sigma_8 *=.*/sigma_8       = $SIGMA8/" -e "s/^n_s *=.*/n_s           = $NS/" \
      -e "s|^transfer_file *=.*|transfer_file = $tf|" -e "s/^filename *=.*/filename = lg_$N.hdf5/" \
      "$here/tools/ics/music/c2_uniform.conf" > "$conf"
  (cd "$ROOT/music" && "$MUSIC" "$(basename "$conf")" > music.log 2>&1)
fi

# 3. the code's ICs, the RECFAST table of this cosmology, params.txt
prob=$ROOT/prob
[ -f "$prob/params.txt" ] || python3 "$here/tools/production/make_production.py" "$prob" \
  --N "$N" --ndm "$N" --box "$BOX" --z-init "$ZINIT" "${cos[@]}" --sigma8 "$SIGMA8" --ns "$NS" \
  --plane-parallel --plane-z-re 6 --zdata1 "$ZOUT" --data-files all --movie-dz 0 \
  --restart-every 0 --no-ics
[ -f "$prob/ics/xv_dm.init" ] || python3 "$here/tools/ics/music_convert.py" "$ROOT/music/lg_$N.hdf5" "$prob/ics" \
  --conf "$conf" --omegam "$OMEGAM" --omegab "$OMEGAB" --hubble0 "$HUBBLE0" --sigma8 "$SIGMA8" --ns "$NS"
if [ ! -f "$prob/recfast.txt" ]; then
  python3 "$here/tools/thermal/make_recfast_table.py" --out "$prob/recfast.txt" "${cos[@]}"
  {
    echo "# hydro and gravity only to the last zdata1 output (the plane-parallel RT turns on at"
    echo "# plane_z_re = 6, after the run); the RECFAST thermal history of this cosmology"
    echo "thermal_init = recfast"
    echo "recfast_table = $prob/recfast.txt"
  } >> "$prob/params.txt"
fi

# 4. the run
if [ ! -f "$ROOT/run.done" ]; then
  rm -rf "$ROOT/run"
  /usr/bin/time -v "$EXE" "$prob" "$ROOT/run" ${EXTRA:-} > "$ROOT/run.log" 2>&1
  touch "$ROOT/run.done"
fi

# 5. the comparison with CAMB
python3 "$here/tools/analysis/linear_growth.py" "$ROOT" --z $ZOUT \
  --out "$ROOT/linear_growth.png" --json "$ROOT/linear_growth.json" | tee "$ROOT/linear_growth.txt"
echo "$(date) done"
