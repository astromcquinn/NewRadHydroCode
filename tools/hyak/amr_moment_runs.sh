#!/bin/bash
# The moment RT on the refined levels (docs/moment_amr.md section 4): for the
# rays' refined production run BASE (default $RUNS/fesc512_b2.5_amr: 40 Mpc/h,
# 512^3 + two levels, 1024^3 particles on the 2048^3 mesh, 0.1c, case B,
# f6 = 0.135, the sink, docs/stage_f.md section 3) makes the siblings
# <BASE>_otvet_cone and <BASE>_m1 with the same params.txt, ICs and dark
# matter plus the RT method, so the pair differs in the RT only:
#   rt_method = otvet, otvet_closure = cone (or rt_method = m1)
#   amrex_balance = 0        (the ray-weighted balancing is for the rays)
#   amrex_slab_planes = 4    (the cone tensor's slab shadow of level 0:
#                             a multiple of otvet_cone_coarse)
#   photon_consume_dense = 0 (the base run predates P11; the comparison is
#                             of the methods)
#   restart_every = 10, the suite's fields and movie fields with gammaHI_avg
# and submits tools/hyak/reion_production.slurm for each (5 nodes, 120 ranks
# x 8 threads, as the base run). Submit from the repository root on a login
# node after tools/hyak/build_amrex.sh built build-amrex/radhydro from this
# branch:
#   tools/hyak/amr_moment_runs.sh [cone m1]      # default: cone
#   DRYRUN=1 tools/hyak/amr_moment_runs.sh
# The uniform twin of the pair (the same without refinement, one node):
#   BASE=$RUNS/fesc512_b2.5_uniform NODES=1 SBATCH_ARGS="--exclusive --nodes=1 --ntasks=24 --cpus-per-task=8 --mem=1000G --time=3-00:00:00" tools/hyak/amr_moment_runs.sh cone
# Environment: RUNS, BASE, EXE (default build-amrex/radhydro), NODES (5),
# SBATCH_ARGS, OVERWRITE=1 (remove an existing sibling), EXTRA (more
# "key = value" lines, ';'-separated).
set -eu
RUNS="${RUNS:-/gscratch/astro/$USER/radhydro/runs}"
BASE="${BASE:-$RUNS/fesc512_b2.5_amr}"
EXE="${EXE:-build-amrex/radhydro}"
NODES="${NODES:-5}"
SBATCH_ARGS="${SBATCH_ARGS:---exclusive --nodes=$NODES --ntasks=$((24 * NODES)) --cpus-per-task=8 --mem=1400G --time=4-00:00:00}"
[ -x "$EXE" ] || { echo "$EXE missing: run tools/hyak/build_amrex.sh first"; exit 1; }
[ -f "$BASE/params.txt" ] || { echo "$BASE/params.txt missing"; exit 1; }
# a refined base (the levels: the base level on boxes, the tensor on the slab
# shadow) or a uniform one (the moment level on z-slabs: its twin without
# refinement, e.g. BASE=$RUNS/fesc512_b2.5_uniform NODES=1)
if grep -Eq '^[[:space:]]*amrex_max_level[[:space:]]*=[[:space:]]*[1-9]' "$BASE/params.txt"; then refined=1; else refined=0; fi
ics="$(cd "$BASE" && readlink -f ics)"
submit() { if [ "${DRYRUN:-0}" = 1 ]; then echo "DRYRUN: $*"; else "$@"; fi; }
for m in "${@:-cone}"; do
  case "$m" in
    cone) d="${BASE}_otvet_cone"; keys="rt_method = otvet;otvet_closure = cone" ;;
    m1) d="${BASE}_m1"; keys="rt_method = m1" ;;
    otvet) d="${BASE}_otvet"; keys="rt_method = otvet" ;;
    *) echo "unknown method $m (cone, m1, otvet)"; exit 2 ;;
  esac
  if [ -e "$d" ]; then
    if [ "${OVERWRITE:-0}" = 1 ]; then echo "removing $d"; [ "${DRYRUN:-0}" = 1 ] || rm -rf "$d"; else echo "exists: $d (OVERWRITE=1 to remove)"; exit 1; fi
  fi
  mkdir -p "$d"
  {
    echo "# $(basename "$BASE") with the moment RT on the refined levels: $keys (tools/hyak/amr_moment_runs.sh, docs/moment_amr.md)"
    grep -Ev '^[[:space:]]*(ics|rt_method|otvet_closure|otvet_split|amrex_balance|amrex_slab|amrex_slab_planes|photon_consume_dense|restart_every|fields|movie_fields|ray_deposit)[[:space:]]*=' "$BASE/params.txt"
    echo "ics = $ics"
    echo "$keys" | tr ';' '\n'
    echo "amrex_balance = 0"
    [ "$refined" = 1 ] || echo "amrex_slab = 1"
    echo "amrex_slab_planes = 4"
    echo "photon_consume_dense = 0"
    echo "restart_every = 10"
    echo "fields = fHI,fHeI,fHeII,temp,rhog,gammaHI,gammaHI_avg"
    echo "movie_fields = fHI,temp,rhog,gammaHI,gammaHI_avg"
    [ -n "${EXTRA:-}" ] && echo "$EXTRA" | tr ';' '\n'
  } > "$d/params.txt"
  ln -sfn "$ics" "$d/ics"
  submit sbatch $SBATCH_ARGS --job-name="$(basename "$d")" --output="$d/slurm-%j.out" \
    --export=ALL,PROB="$d",RANKS=$((24 * NODES)),THREADS=8,EXE="$EXE",POST=sbatch tools/hyak/reion_production.slurm
done
