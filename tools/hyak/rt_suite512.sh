#!/bin/bash
# The uniform 512^3 rerun of every RT method for the moment-methods paper
# (docs/PLAN.md, "Unified RT-methods suite", decisions of 2026-10-18):
# 40 Mpc/h, 512^3 gas, 1024^3 particles on the 2048^3 PM mesh
# (pm_refine = 4, the uniform-grid refined mesh of docs/pm_refine.md, no
# AMR), the reduced speed of light by ionisation state (c_reduced_schedule,
# the make_production.py default: 0.1c below x_HII,V = 0.25, 0.2c to 0.5,
# 0.3c to 0.75, 0.5c after), f_esc(z) = F6 ((1+z)/7)^BETA with the
# luminosities of the 512^3 calibration (data/uv_lf/am_params_128_10mpc.txt,
# docs/fesc_evolution.md), case B recombination without recombination
# radiation for the base and the moment methods (RECOMB=B, user decision
# 2026-10-18) plus one ray run with case A + transported recombination
# radiation and recomb_rad_cell_onspot = 1 (the caseA sibling; RECOMB=A
# swaps the two), the original cooling suppression (200) and no
# star-formation sink, no dark matter in the outputs (data_files = gas: the
# particles are only in the checkpoints), the substep-averaged Gamma_HI
# (gammaHI_avg) in the snapshots and the slab movie frames. Every run is on
# the serial code, one whole cpu-g2 node each (the moment methods exist only
# there), so all of them get the same physics.
#
# Runs ($RUNS/<name>; the names of the 2026-10-11 set, overwritten):
#   reion512_rays                rays, the base and the reference (case B)
#   reion512_rays_caseA          rays, case A + recombination radiation
#   reion512_rays_m1             M1
#   reion512_rays_otvet          OTVET, thin tensor
#   reion512_rays_otvet_split    OTVET, direct + diffuse split
#   reion512_rays_otvet_screen   OTVET, screened kernel
#   reion512_rays_otvet_cone     OTVET, HI-attenuated cone tensor
#
# usage (repository root, a klone login node, after tools/hyak/build_hyak.sh):
#   tools/hyak/rt_suite512.sh probe     # PROBE_STEPS steps of the base into
#                                       # $BASE/probe: peak memory and PM time
#                                       # at the 2048^3 mesh (/usr/bin/time -v
#                                       # in $BASE/probe-<job>.out)
#   tools/hyak/rt_suite512.sh base      # make and submit the base ray run
#   tools/hyak/rt_suite512.sh siblings  # after the base passed its check
#                                       # (z99, Gamma_HI at z = 5.4): make and
#                                       # submit the six siblings
#   tools/hyak/rt_suite512.sh post      # rt_methods_post.slurm on the set
# OVERWRITE=1 removes an existing directory of a run this phase makes
# (otherwise the script stops and lists them with their sizes). The old
# set's checkpoints are most of its ~1 TB.
#
# Environment: RUNS (default /gscratch/astro/$USER/radhydro/runs), ICS
# (default $RUNS/hmf512_40mpc_r2/ics: 512^3 gas, 1024^3 particles, 40 Mpc/h,
# seed 7 with the 512^3 large-scale phases, docs/pm_refine.md), BASE
# (default $RUNS/reion512_rays), F6 (0.135), BETA (2.5), AM_PARAMS
# (data/uv_lf/am_params_128_10mpc.txt), RECOMB (B, or A), SIBLINGS (default
# "caseA m1 otvet split screen cone"; caseA / caseB is the ray run with the
# other recombination treatment, skipped when it is the base's), RESTART_KEEP
# (1: checkpoint writes are atomic, docs/production_run.md section 3; 2 to
# keep the previous one too), SBATCH_ARGS (default a whole node: --exclusive
# --cpus-per-task=192 --mem=1000G --time=4-00:00:00), THREADS (192), EXE
# (build/radhydro), PROBE_STEPS (3), MAKE_EXTRA (more make_production.py
# arguments), DRYRUN=1 (make the directories, print the sbatch commands).
#
# Disk per run at 512^3 (estimates): snapshots 9 x (7 field cubes of 0.5 GB
# + Gas/u_gas 2.7 GB) = 55 GB; slab movie frames ~400 x 5 fields x 25 MB =
# 50 GB; checkpoints (restart_every = 10, RESTART_KEEP kept): case-B rays
# ~135 GB each, case A + RR rays ~190 GB (1.5x the rays), moments ~45 GB,
# each with the 1024^3 particles (26 GB). With one kept a case-B ray run
# holds ~240 GB while it runs, the case A + RR run ~295 GB, a moment run
# ~150 GB, all ~105 GB once the checkpoints are deleted; the set ~1.3 TB
# while running (1.9 TB with RESTART_KEEP=2), ~0.75 TB kept. Memory:
# the serial PM at the 2048^3 mesh adds ~275 GB of host arrays to the RT's
# 150-250 GB (within the node's 1.5 TB; the probe measures it).
set -eu
RUNS="${RUNS:-/gscratch/astro/$USER/radhydro/runs}"
ICS="${ICS:-$RUNS/hmf512_40mpc_r2/ics}"
BASE="${BASE:-$RUNS/reion512_rays}"
F6="${F6:-0.135}"
BETA="${BETA:-2.5}"
AM_PARAMS="${AM_PARAMS:-data/uv_lf/am_params_128_10mpc.txt}"
RECOMB="${RECOMB:-B}"
case "$RECOMB" in
  A) recomb_args=(--recomb-case A --recomb-rad 1 --recomb-rad-cell-onspot 1) ;;
  B) recomb_args=(--recomb-case B --recomb-rad 0) ;;
  *) echo "RECOMB must be A or B"; exit 2 ;;
esac
RESTART_KEEP="${RESTART_KEEP:-1}"
SIBLINGS="${SIBLINGS:-caseA m1 otvet split screen cone}"
THREADS="${THREADS:-192}"
SBATCH_ARGS="${SBATCH_ARGS:---exclusive --cpus-per-task=192 --mem=1000G --time=4-00:00:00}"
EXE="${EXE:-build/radhydro}"
PROBE_STEPS="${PROBE_STEPS:-3}"
export EXE THREADS   # reach reion_production.slurm through the sibling scripts' --export=ALL
phase="${1:-}"
case "$phase" in probe|base|siblings|post) ;; *) sed -n 2,40p "$0"; exit 2 ;; esac
[ -x "$EXE" ] || { echo "$EXE missing: run tools/hyak/build_hyak.sh first"; exit 1; }
[ -f "$AM_PARAMS" ] || { echo "$AM_PARAMS missing"; exit 1; }

submit() {  # submit CMD...: run or print
  if [ "${DRYRUN:-0}" = 1 ]; then echo "DRYRUN: $*"; else "$@"; fi
}

clear_dirs() {  # clear_dirs DIR...: remove with OVERWRITE=1, else stop
  local d stop=0
  for d in "$@"; do
    [ -e "$d" ] || continue
    if [ "${OVERWRITE:-0}" = 1 ]; then
      echo "removing $d ($(du -sh "$d" 2>/dev/null | cut -f1))"
      [ "${DRYRUN:-0}" = 1 ] || rm -rf "$d"
    else
      echo "exists: $d ($(du -sh "$d" 2>/dev/null | cut -f1))"; stop=1
    fi
  done
  [ "$stop" = 0 ] || { echo "set OVERWRITE=1 to remove the directories above and remake them"; exit 1; }
}

check_ics() {
  local info="$ICS/ics_info.txt"
  if [ ! -f "$ICS/u_gas.init" ] || [ ! -f "$ICS/xv_dm.init" ] || [ ! -f "$info" ] ||
     ! grep -Eq '^Ngas = 512$' "$info" || ! grep -Eq '^Ndm = 1024$' "$info" || ! grep -Eq '^box = 40(\.0)?$' "$info"; then
    echo "no 512^3 gas / 1024^3 particle / 40 Mpc/h ICs at $ICS (ics_info.txt: Ngas = 512, Ndm = 1024, box = 40)."
    echo "Point ICS at them, or make them on a compute node (~300 GB of memory, about an hour):"
    echo "  python3 tools/ics/make_ics.py <dir> --kind zeldovich --Ngas 512 --Ndm 1024 --z-init 100 --box 40 --seed 7 --phases-from-ndm 512"
    exit 1
  fi
}

make_base() {  # the base's params.txt (fresh; an old run of the name needs OVERWRITE=1) and ics link
  check_ics
  local tmp; tmp="$(mktemp -d)"
  # (the defaults: the schedule, unsplit hydro, redshift per distance,
  # case A + RR with onspot = 1, moment_inject = 1, ray_rotate = 3,
  # photon_conserve from legacy_physics = 0)
  python3 tools/production/make_production.py "$tmp" --N 512 --box 40 --pm-refine 4 --no-ics \
    --nray-max-factor 15 --am-params "$AM_PARAMS" --f-esc "$F6" --f-esc-beta "$BETA" --f-esc-z0 6 \
    --xi-ion 25.38 --movie-dz 0.025 --movie-slab 8 --movie-fields fHI,temp,rhog,gammaHI,gammaHI_avg \
    --fields fHI,fHeI,fHeII,temp,rhog,gammaHI,gammaHI_avg --zdata1 "15 7 6.5 6.25 6 5.75 5.5 5.25 5" \
    --restart-every 10 --restart-keep "$RESTART_KEEP" --data-files gas --cool-suppress-overdensity 200 \
    "${recomb_args[@]}" ${MAKE_EXTRA:-} > /dev/null
  # 1024^3 particles on the 2048^3 mesh (make_production.py writes
  # Ndm = pm_refine * N; the code needs only that Ndm divides the mesh)
  sed -i 's/^Ndm = .*/Ndm = 1024/' "$tmp/params.txt"
  grep -Eq '^Ndm = 1024$' "$tmp/params.txt" || { echo "could not set Ndm = 1024"; exit 1; }
  grep -Eq '^pm_refine = 4$' "$tmp/params.txt" || { echo "pm_refine = 4 missing"; exit 1; }
  grep -Eq '^c_reduced_schedule' "$tmp/params.txt" || { echo "c_reduced_schedule missing: is this make_production.py from main after PR #46?"; exit 1; }
  grep -Eq "^recomb_case = $RECOMB\$" "$tmp/params.txt" || { echo "recomb_case = $RECOMB missing in the generated params.txt"; exit 1; }
  if [ -f "$BASE/params.txt" ] && diff -q <(grep -v '^#' "$tmp/params.txt") <(grep -v '^#' "$BASE/params.txt") > /dev/null; then
    : # made by this script already (e.g. by the probe phase): keep it
  else
    clear_dirs "$BASE"   # an old run of this name (or other parameters): OVERWRITE=1 removes it
    mkdir -p "$BASE"
    cp "$tmp/params.txt" "$BASE/params.txt"
  fi
  rm -rf "$tmp"
  [ -e "$BASE/ics" ] || ln -s "$ICS" "$BASE/ics"
  echo "base: $BASE"; grep -E '^(N|Ndm|box|pm_refine|c_reduced|c_reduced_schedule|f_esc|f_esc_beta|recomb_case|recomb_rad|recomb_rad_cell_onspot|data_files|fields|movie_fields|cool_suppress_overdensity|restart_every|restart_keep) ' "$BASE/params.txt" | sed 's/^/  /'
}

case "$phase" in
probe)
  make_base
  clear_dirs "$BASE/probe"
  mkdir -p "$BASE/probe"
  # a short serial run of the base's first steps (no RT yet at z = 100:
  # the hydro, the 1024^3 particles on the 2048^3 mesh, the FoF): peak
  # RSS and the PM time per step (timing.txt, the "Maximum resident set
  # size" line of /usr/bin/time)
  cat > "$BASE/probe/probe.slurm" <<EOF
#!/bin/bash
#SBATCH --job-name=rt512_probe
#SBATCH --account=astro
#SBATCH --partition=cpu-g2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=$BASE/probe-%j.out
set -u
cd "$PWD"
source tools/hyak/env.sh
export OMP_NUM_THREADS=$THREADS OMP_PLACES=cores OMP_PROC_BIND=close
/usr/bin/time -v $EXE "$BASE" "$BASE/probe" --quiet --set run_steps=$PROBE_STEPS --set restart_every=0 --set movie_dz=0
grep -E "^#|^ *[0-9]" "$BASE/probe/timing.txt" | tail -n +1
EOF
  submit sbatch $SBATCH_ARGS --time=06:00:00 "$BASE/probe/probe.slurm"
  ;;
base)
  make_base
  clear_dirs "$BASE/run"   # outputs of an earlier submission of this base
  submit sbatch $SBATCH_ARGS --job-name=reion512_rays --output="$BASE/slurm-%j.out" \
    --export=ALL,PROB="$BASE",THREADS="$THREADS",EXE="$EXE",POST=sbatch tools/hyak/reion_production.slurm
  ;;
siblings)
  [ -f "$BASE/params.txt" ] || { echo "$BASE/params.txt missing: run the base phase first"; exit 1; }
  # the ray run with the other recombination treatment: never the base's own
  if [ "$RECOMB" = B ]; then SIBLINGS="${SIBLINGS//caseB/}"; else SIBLINGS="${SIBLINGS//caseA/}"; fi
  dirs=()
  for s in $SIBLINGS; do
    case "$s" in
      caseA|caseB) dirs+=("${BASE}_$s") ;;
      m1|otvet) dirs+=("${BASE}_$s") ;;
      split|screen|cone) dirs+=("${BASE}_otvet_$s") ;;
      *) echo "unknown sibling $s"; exit 2 ;;
    esac
  done
  clear_dirs "${dirs[@]}"
  ics="$(cd "$BASE" && readlink -f ics)"
  methods=""; closures=""
  for s in $SIBLINGS; do
    case "$s" in
      caseA|caseB)
        d="${BASE}_$s"; mkdir -p "$d"
        {
          if [ "$s" = caseA ]; then
            echo "# $(basename "$BASE") with case A + transported recombination radiation (tools/hyak/rt_suite512.sh)"
          else
            echo "# $(basename "$BASE") with case B and no recombination radiation (tools/hyak/rt_suite512.sh)"
          fi
          grep -Ev '^[[:space:]]*(ics|recomb_case|recomb_rad|recomb_rad_block|recomb_rad_cell_onspot|recomb_rad_lines|recomb_rad_pool)[[:space:]]*=' "$BASE/params.txt"
          echo "ics = $ics"
          if [ "$s" = caseA ]; then
            echo "recomb_case = A"
            echo "recomb_rad = 1"
            echo "recomb_rad_block = 4"
            echo "recomb_rad_cell_onspot = 1"
          else
            echo "recomb_case = B"
            echo "recomb_rad = 0"
          fi
        } > "$d/params.txt"
        ln -s "$ics" "$d/ics"
        submit sbatch $SBATCH_ARGS --job-name="reion512_$s" --output="$d/slurm-%j.out" \
          --export=ALL,PROB="$d",THREADS="$THREADS",EXE="$EXE",POST=sbatch tools/hyak/reion_production.slurm
        ;;
      m1|otvet) methods="$methods $s" ;;
      split|screen|cone) closures="$closures $s" ;;
    esac
  done
  # the moment methods and the closures: the existing sibling scripts (the
  # base's params.txt plus rt_method and the closure keys, the same ics)
  if [ -n "$methods" ]; then
    METHODS="$methods" SBATCH_ARGS="$SBATCH_ARGS" THREADS="$THREADS" DRYRUN="${DRYRUN:-0}" tools/hyak/rt_methods_runs.sh "$BASE"
  fi
  if [ -n "$closures" ]; then
    SBATCH_ARGS="$SBATCH_ARGS" THREADS="$THREADS" DRYRUN="${DRYRUN:-0}" tools/hyak/otvet_closure_runs.sh "$BASE" $closures
  fi
  ;;
post)
  b="$(basename "$BASE")"
  extra=""
  for s in $SIBLINGS; do
    case "$s" in
      caseA) extra="$extra ${b}_caseA:rays_A_RR" ;;
      caseB) extra="$extra ${b}_caseB:rays_B" ;;
      split|screen|cone) extra="$extra ${b}_otvet_$s:$s" ;;
    esac
  done
  submit env RUNS="$RUNS" SETS="$b:rt_methods_512" EXTRA="$extra" IGM=1 sbatch --time=8:00:00 tools/hyak/rt_methods_post.slurm
  ;;
esac
