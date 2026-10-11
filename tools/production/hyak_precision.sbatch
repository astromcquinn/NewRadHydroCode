#!/bin/bash
# Lossy ray-precision check on the production configuration (Hyak klone,
# docs/ray_compression.md §6, open item). One job array task per variant;
# each task runs the full radhydro loop at N = 128, 10 Mpc/h, 0.01c,
# nray_max_factor 16, all fixes on, deterministic parallel ray tracing.
#
#   cd NewRadHydroCode
#   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j
#   python3 tools/production/make_production.py RUNS/p128 --N 128 --movie-dz 0
#   sbatch --export=ALL,ROOT=$PWD/RUNS/p128 tools/production/hyak_precision.sbatch
#
# Afterwards, in RUNS/p128:
#   python3 tools/ray_compression/precision_metrics.py RUNS/p128/prec --N 128 --json metrics.json
#   python3 tools/ray_compression/slices.py RUNS/p128/prec p128 --N 128 --out slices.png
# (precision_metrics.py expects RUNS/p128/prec/<problem>/<variant>/step_*; the
# task below lays the outputs out that way, with <problem> = p128.)
#
# Resources: at N = 128 one run needs ~3 GB and ~30 min on 4 threads
# (docs/production_run.md §5); 24 threads per task is a safe over-allocation
# on a 192-core cpu-g2 node. Adjust --account/--partition to your allocation.
#SBATCH --job-name=raycomp-prec
#SBATCH --account=astro
#SBATCH --partition=cpu-g2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --array=0-5
#SBATCH --output=raycomp-prec-%A_%a.out

set -euo pipefail
: "${ROOT:?set ROOT=/path/to/RUNS/p128 (from make_production.py) with --export}"
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Under sbatch the script runs from a spool copy: use the submit directory (the repository root).
exe="${EXE:-${SLURM_SUBMIT_DIR:-$here/../..}/build/radhydro}"

# variant name and the research-only ray_quant settings (ray_quant.hpp)
variants=(
  "base:"
  "d8:ray_quant_bits=8"
  "d4:ray_quant_bits=4"
  "d1:ray_quant_bits=1"
  "r6:ray_quant_bits=6 ray_quant_dither=0"
  "struct_o10_d8:ray_quant_struct=1 ray_quant_offbits=10 ray_quant_bits=8"
)
v="${variants[$SLURM_ARRAY_TASK_ID]}"
name="${v%%:*}"
sets="${v#*:}"
args=()
for kv in $sets; do args+=(--set "$kv"); done

out="$ROOT/prec/p128/$name"
mkdir -p "$out"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores OMP_PROC_BIND=close
echo "variant $name: ${args[*]:-} ($OMP_NUM_THREADS threads)"
# outputs at every zdata1 redshift; the metrics need fHI, temp, rhog
"$exe" "$ROOT" "$out" --quiet --set restart_every=0 --set movie_dz=0 \
    --set "fields=fHI,temp,rhog" --set threads="$OMP_NUM_THREADS" "${args[@]}" \
    > "$out/run.log" 2>&1
touch "$out/done"
echo "done $name"
