#!/bin/bash
# Convergence of the Gamma_HI fluctuations of M1 and of the rays
# (docs/rt_moments.md section 8; asked by the user on 2026-11-20 after the
# uniform-opacity test of section 7 showed M1's power above k ~ (2/lambda dx)^1/2
# to be set by the first-order transport).
#
# Two series in the 512^3 suite's configuration (tools/hyak/rt_suite512.sh:
# the same make_production.py arguments) in a 20 Mpc/h box:
#
#  res   the gas / RT grid N = 128, 256, 512 at FIXED dark matter: 512^3
#        particles on the 512^3 PM mesh in every run (pm_refine = 4, 2, 1),
#        ICs with identical particles (make_ics.py --Ndm 512 --phases-from-ndm
#        256: xv_dm.init is byte-identical for every --Ngas), so the halos
#        and with them the sources are the same up to the gas's gravity. For
#        each N the rays, M1 and M1 with moment_inject = 0:
#        $RUNS/rtconv_N<N>_rays, ..._rays_m1, ..._rays_m1_inj0.
#        N = 256 is the suite's cell size (its pre-check configuration).
#        The measurement: P_Gamma(M1)/P_Gamma(rays) at fixed k against N.
#        The uniform-opacity test predicts that M1's ratio falls with N at
#        k >~ 1 h/Mpc after overlap and that the k where M1 overtakes the rays
#        stays at k dx ~ 1.3. (The gas resolution also changes the absorbers
#        and the mean free path: compare at matched x_V and quote lambda.)
#
#  rays  the ray budget at N = 256: nray_max_factor = RAYFACTORS (default
#        "8 32 64 120"; the res series' own run has the suite's 15), every
#        other parameter the same: $RUNS/rtconv_N256_rays_f<F>. Rays fainter
#        than the brightest quarter of the pool are merged into the 48 cones
#        of nside_merge = 2 per cell and only the brightest half may split
#        (calc_ngam_ray), so nray_max_factor sets the rays' angular resolution
#        (Wu, McQuinn & Eisenstein 2021 section 5.2.3: ~l^2 directions
#        converge for k <~ l/lambda). The measurement: P_Gamma against the
#        f120 run, at fixed k. nray_max_factor * N^3 must stay below 2^31, so
#        the 512^3 suite itself cannot go above 15.
#
# usage (repository root, a klone login node, after tools/hyak/build_hyak.sh):
#   tools/hyak/rt_convergence.sh all    # everything in one go: the ICs jobs,
#                                       # every run after its ICs
#                                       # (--dependency=afterok), the analysis
#                                       # after all the runs (afterany)
# or phase by phase:
#   tools/hyak/rt_convergence.sh ics    # the ICs jobs (one per N without ICs)
#   tools/hyak/rt_convergence.sh res    # the rays and M1 at every N of RES
#   tools/hyak/rt_convergence.sh rays   # the ray budgets at N = 256
#   tools/hyak/rt_convergence.sh inj0   # only the M1 moment_inject = 0 runs, beside
#                                       # a res series that is already submitted
#   tools/hyak/rt_convergence.sh post   # gamma_power.py on every set ->
#                                       # $RUNS/rtconv_post/{N<N>,rays}/gamma_pk*.txt
# Environment: RUNS (default /gscratch/astro/$USER/radhydro/runs), RES
# ("128 256 512"), RAYFACTORS ("8 32 64 120"), F6 (0.10), BETA (2.5),
# AM_PARAMS (data/uv_lf/am_params_st_fof.txt; with F6 and the continuous
# c-hat schedule these are the late-end suite's, $RUNS/reion512late, not
# rt_suite512.sh's own defaults: the 128-box table is wrong for this
# particle mass), ICS_<N> (the ICs of N;
# default $RUNS/rtconv_ics_N<N>, and for 256 the suite's pre-check ICs
# $RUNS/rt512cal256/ics when they exist), THREADS (192), EXE
# (build/radhydro), SBATCH_ARGS_<N> (defaults below), MAKE_EXTRA, DRYRUN=1,
# OVERWRITE=1 (as rt_suite512.sh), SCRATCH (/gscratch/scrubbed/$USER/rtconv:
# where run/ and run.log of every run are written; none for $RUNS).
set -eu
RUNS="${RUNS:-/gscratch/astro/$USER/radhydro/runs}"
RES="${RES:-128 256 512}"
RAYFACTORS="${RAYFACTORS:-8 32 64 120}"
F6="${F6:-0.10}"
BETA="${BETA:-2.5}"
AM_PARAMS="${AM_PARAMS:-data/uv_lf/am_params_st_fof.txt}"
THREADS="${THREADS:-192}"
EXE="${EXE:-build/radhydro}"
SBATCH_ARGS_128="${SBATCH_ARGS_128:---exclusive --cpus-per-task=192 --mem=200G --time=12:00:00}"
SBATCH_ARGS_256="${SBATCH_ARGS_256:---exclusive --cpus-per-task=192 --mem=400G --time=2-00:00:00}"
SBATCH_ARGS_512="${SBATCH_ARGS_512:---exclusive --cpus-per-task=192 --mem=800G --time=4-00:00:00}"
# the runs' outputs (run/, run.log) go to SCRATCH/<name>, linked from the
# run directory (user decision 2026-11-20: /gscratch/scrubbed; purged after
# 21 days without access, so analyse and copy the tables in time);
# SCRATCH=none keeps them under $RUNS
SCRATCH="${SCRATCH:-/gscratch/scrubbed/$USER/rtconv}"
NDM=512; BOX=20
export EXE THREADS
phase="${1:-}"
case "$phase" in all|ics|res|rays|inj0|post) ;; *) sed -n 2,55p "$0"; exit 2 ;; esac
declare -A ICSJOB=()   # N -> the job making its ICs (phase all)
RUNJOBS=""             # the run jobs submitted, for the analysis' dependency
[ -x "$EXE" ] || { echo "$EXE missing: run tools/hyak/build_hyak.sh first"; exit 1; }
[ -f "$AM_PARAMS" ] || { echo "$AM_PARAMS missing"; exit 1; }

submit() { if [ "${DRYRUN:-0}" = 1 ]; then echo "DRYRUN: $*"; else "$@"; fi; }
JID=""
submit_job() {  # submit_job SBATCH-ARGS...: sbatch, the job id in JID (DRY<n> with DRYRUN=1)
  if [ "${DRYRUN:-0}" = 1 ]; then
    DRYN=$((${DRYN:-0} + 1)); JID="DRY$DRYN"; echo "DRYRUN [$JID]: sbatch $*"
  else
    JID="$(sbatch --parsable "$@")"; JID="${JID%%;*}"; echo "job $JID: ${*: -1}"
  fi
}

clear_dirs() {
  local d stop=0
  for d in "$@"; do
    [ -e "$d" ] || continue
    if [ "${OVERWRITE:-0}" = 1 ]; then
      echo "removing $d"; [ "${DRYRUN:-0}" = 1 ] || rm -rf "$d"
    else
      echo "exists: $d"; stop=1
    fi
  done
  [ "$stop" = 0 ] || { echo "set OVERWRITE=1 to remove the directories above and remake them"; exit 1; }
}

ics_dir() {  # ics_dir N
  local v="ICS_$1"
  if [ -n "${!v:-}" ]; then echo "${!v}"; return; fi
  if [ "$1" = 256 ] && [ -f "$RUNS/rt512cal256/ics/xv_dm.init" ]; then echo "$RUNS/rt512cal256/ics"; return; fi
  echo "$RUNS/rtconv_ics_N$1"
}

sbatch_args() { local v="SBATCH_ARGS_$1"; echo "${!v}"; }

make_params() {  # make_params DIR N [more make_production.py arguments]: rt_suite512.sh's, in the 20 Mpc/h box
  local dir="$1" n="$2"; shift 2
  [ $((NDM % n)) = 0 ] || { echo "N = $n must divide the $NDM^3 PM mesh"; exit 2; }
  python3 tools/production/make_production.py "$dir" --N "$n" --box "$BOX" --pm-refine $((NDM / n)) --ndm "$NDM" \
    --f-esc "$F6" --movie-dz 0.05 --no-ics \
    --nray-max-factor 15 --am-params "$AM_PARAMS" --f-esc-beta "$BETA" --f-esc-z0 6 \
    --xi-ion 25.38 --movie-slab 8 --movie-fields fHI,temp,rhog,gammaHI,gammaHI_avg \
    --fields fHI,fHeI,fHeII,temp,rhog,gammaHI,gammaHI_avg --zdata1 "15 7 6.5 6.25 6 5.75 5.5 5.25 5" \
    --restart-every 10 --restart-keep 1 --data-files gas --cool-suppress-overdensity 200 \
    --recomb-case B --recomb-rad 0 --c-schedule-interp 1 "$@" ${MAKE_EXTRA:-} > /dev/null
  grep -Eq "^Ndm = $NDM\$" "$dir/params.txt" || { echo "Ndm = $NDM missing in $dir/params.txt"; exit 1; }
}

need_ics() {  # need_ics N: stop unless the ICs of N exist or their job is queued (DRYRUN only warns)
  local d; d="$(ics_dir "$1")"
  [ -z "${ICSJOB[$1]:-}" ] || return 0
  if [ ! -f "$d/xv_dm.init" ] || [ ! -f "$d/u_gas.init" ]; then
    echo "no ICs at $d: run the ics phase first"
    [ "${DRYRUN:-0}" = 1 ] || exit 1
  fi
}

launch() {  # launch DIR N NAME: a whole node on reion_production.slurm, after the ICs job of N if any
  local dep=() sc
  [ -z "${ICSJOB[$2]:-}" ] || dep=(--dependency="afterok:${ICSJOB[$2]}")
  if [ "$SCRATCH" != none ] && [ ! -e "$1/run" ]; then
    sc="$SCRATCH/$(basename "$1")"; mkdir -p "$sc/run"
    ln -s "$sc/run" "$1/run"; ln -sfn "$sc/run.log" "$1/run.log"
  fi
  submit_job $(sbatch_args "$2") "${dep[@]}" --job-name="$3" --output="$1/slurm-%j.out" \
    --export=ALL,PROB="$1",THREADS="$THREADS",EXE="$EXE",POST=sbatch tools/hyak/reion_production.slurm
  RUNJOBS="$RUNJOBS:$JID"
}

make_inj0() {  # make_inj0 N: the M1 sibling of rtconv_N<N>_rays with moment_inject = 0, launched
  local n="$1" d="$RUNS/rtconv_N${1}_rays" m
  # M1 with the sources injected before the transport (moment_inject = 0):
  # without the current substep's photons in the source cells, a white
  # term 1 - exp(-dx/lambda) in T(k) (docs/rt_moments.md section 7)
  m="${d}_m1_inj0"; mkdir -p "$m"
  { echo "# $(basename "$d") with rt_method = m1, moment_inject = 0 (tools/hyak/rt_convergence.sh)"
    grep -Ev '^[[:space:]]*(ics|moment_inject)[[:space:]]*=' "$d/params.txt"
    echo "ics = $(ics_dir "$n")"; echo "rt_method = m1"; echo "moment_inject = 0"; } > "$m/params.txt"
  ln -s "$(ics_dir "$n")" "$m/ics"
  launch "$m" "$n" "rtconv${n}_m1i0"
}

do_ics() {
  for n in $RES; do
    d="$(ics_dir "$n")"
    [ -f "$d/xv_dm.init" ] && { echo "N = $n: ICs at $d"; continue; }
    mkdir -p "$(dirname "$d")"
    job="$RUNS/rtconv_ics_N$n.slurm"
    cat > "$job" <<EOF
#!/bin/bash
#SBATCH --job-name=rtconv_ics$n
#SBATCH --account=astro
#SBATCH --partition=cpu-g2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=02:00:00
#SBATCH --output=$RUNS/rtconv_ics_N$n-%j.out
set -eu
cd "$PWD"
source tools/hyak/env.sh
python3 tools/ics/make_ics.py "$d" --kind zeldovich --Ngas $n --Ndm $NDM --z-init 100 --box $BOX --seed 7 --phases-from-ndm 256
EOF
    submit_job "$job"; ICSJOB[$n]="$JID"
  done
}
do_res() {
  dirs=(); for n in $RES; do dirs+=("$RUNS/rtconv_N${n}_rays" "$RUNS/rtconv_N${n}_rays_m1" "$RUNS/rtconv_N${n}_rays_m1_inj0"); done
  clear_dirs "${dirs[@]}"
  for n in $RES; do
    need_ics "$n"
    d="$RUNS/rtconv_N${n}_rays"
    make_params "$d" "$n"
    ln -s "$(ics_dir "$n")" "$d/ics"
    launch "$d" "$n" "rtconv${n}_rays"
    # the M1 sibling: the same params.txt plus rt_method
    m="${d}_m1"; mkdir -p "$m"
    { echo "# $(basename "$d") with rt_method = m1 (tools/hyak/rt_convergence.sh)"
      grep -v '^[[:space:]]*ics[[:space:]]*=' "$d/params.txt"
      echo "ics = $(ics_dir "$n")"; echo "rt_method = m1"; } > "$m/params.txt"
    ln -s "$(ics_dir "$n")" "$m/ics"
    launch "$m" "$n" "rtconv${n}_m1"
    make_inj0 "$n"
  done
}
do_rays() {
  n=256; need_ics "$n"
  dirs=(); for f in $RAYFACTORS; do dirs+=("$RUNS/rtconv_N${n}_rays_f$f"); done
  clear_dirs "${dirs[@]}"
  for f in $RAYFACTORS; do
    d="$RUNS/rtconv_N${n}_rays_f$f"
    make_params "$d" "$n"
    # the same params.txt with another ray pool (make_production.py's own
    # check: nray_max_factor N^3 < 2^31)
    [ $((f * n * n * n)) -lt 2147483648 ] || { echo "nray_max_factor = $f overflows at N = $n"; exit 2; }
    sed -i -E "s/^nray_max_factor[[:space:]]*=.*/nray_max_factor = $f/" "$d/params.txt"
    grep -Eq "^nray_max_factor = $f\$" "$d/params.txt" || { echo "nray_max_factor not set in $d/params.txt"; exit 1; }
    ln -s "$(ics_dir "$n")" "$d/ics"
    launch "$d" "$n" "rtconv_f$f"
  done
}
do_post() {
  out="$RUNS/rtconv_post"; mkdir -p "$out"
  cmds=()
  for n in $RES; do
    cmds+=("python3 tools/analysis/gamma_power.py $RUNS/rtconv_N${n}_rays/run:rays $RUNS/rtconv_N${n}_rays_m1/run:M1 $RUNS/rtconv_N${n}_rays_m1_inj0/run:M1_inject0 --outdir $out/N$n --box $BOX --field gammaHI_avg --zmin 4.9 --zmax 7.1 --tol 0.05")
  done
  # the ray budgets: the largest first (the reference), then the others and the res series' 15
  ref=""; rest=""
  for f in $(echo $RAYFACTORS | tr ' ' '\n' | sort -nr); do
    if [ -z "$ref" ]; then ref="$RUNS/rtconv_N256_rays_f$f/run:f$f"; else rest="$rest $RUNS/rtconv_N256_rays_f$f/run:f$f"; fi
  done
  cmds+=("python3 tools/analysis/gamma_power.py $ref $RUNS/rtconv_N256_rays/run:f15$rest --outdir $out/rays --box $BOX --field gammaHI_avg --zmin 4.9 --zmax 7.1 --tol 0.05")
  local dep=()
  [ -z "$RUNJOBS" ] || dep=(--dependency="afterany$RUNJOBS")
  for c in "${cmds[@]}"; do
    submit sbatch "${dep[@]}" --account=astro --partition=cpu-g2 --cpus-per-task=16 --mem=128G --time=02:00:00 \
      --job-name=rtconv_post --output="$out/post-%j.out" --wrap="cd $PWD && source tools/hyak/env.sh && $c"
  done
}

case "$phase" in
ics) do_ics ;;
res) do_res ;;
rays) do_rays ;;
inj0)
  # only the moment_inject = 0 siblings, for a res series submitted before
  # they were added (the base's params.txt must exist)
  dirs=(); for n in $RES; do dirs+=("$RUNS/rtconv_N${n}_rays_m1_inj0"); done
  clear_dirs "${dirs[@]}"
  for n in $RES; do
    [ -f "$RUNS/rtconv_N${n}_rays/params.txt" ] || { echo "$RUNS/rtconv_N${n}_rays/params.txt missing: run res first"; exit 1; }
    need_ics "$n"; make_inj0 "$n"
  done
  ;;
post) do_post ;;
all)
  # the 256^3 ray-budget runs need N = 256 ICs even when RES leaves it out
  case " $RES " in *" 256 "*) ;; *) RES="$RES 256" do_ics ;; esac
  do_ics; do_res; do_rays; do_post
  echo "submitted: ICs ${ICSJOB[*]:-(present)}; runs ${RUNJOBS#:}; the tables will be in $RUNS/rtconv_post/{N<N>,rays}/gamma_pk.txt"
  ;;
esac
