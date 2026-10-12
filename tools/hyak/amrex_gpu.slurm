#!/bin/bash
# The AMReX RT backend on one ckpt GPU (A100 or H200) vs one cpu-g2 node
# (docs/stage_d.md section 5; docs/framework_decision.md section 4 holds the
# Parthenon numbers to compare with; PLAN.md: >= 5x on H200 reopens the
# ACCESS question). One hydro step of the 128^3 benchmark from a checkpoint,
# as spike/hyak/bench_node.sbatch. The GPU build is tools/hyak/build_amrex_gpu.sh.
#
#   sbatch --export=ALL,PROB=runs/reion128_c01,CK=<checkpoint> tools/hyak/amrex_gpu.slurm
#
# Environment: PROB, CK (required), OUT (default $PROB/amrex_gpu), STEPS (1),
#   EXE (build-amrex-cuda/radhydro), GPU_PART (ckpt-g2), GPU_TYPE (a100 | h200)
#SBATCH --job-name=amrex-gpu
#SBATCH --account=astro
#SBATCH --partition=ckpt-g2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gpus=1
#SBATCH --mem=128G
#SBATCH --time=04:00:00
#SBATCH --requeue
#SBATCH --output=%x-%j.out
set -u
cd "${SLURM_SUBMIT_DIR}"
source tools/hyak/env.sh
radray_module_load ${HYAK_MODULES:-ompi/4.1.6-2 cuda/12.4.1}
PROB="${PROB:?}"; CK="${CK:?}"
OUT="${OUT:-$PROB/amrex_gpu}"; STEPS="${STEPS:-1}"
EXE="${EXE:-$PWD/build-amrex-cuda/radhydro}"
mkdir -p "$OUT"
{ echo "host $(hostname), job ${SLURM_JOB_ID:-}, $(date)"; nvidia-smi -L; } | tee "$OUT/node.txt"
export OMP_NUM_THREADS=1
rm -rf "$OUT/gpu"
srun --ntasks=1 "$EXE" "$PROB" "$OUT/gpu" --quiet --restart "$CK" --set run_steps="$STEPS" \
  --set restart_every=0 --set movie_dz=0 --set fields_every=0 --set rt_backend=amrex \
  --set deterministic_parallel=0 --set amrex_max_grid_size=64 > "$OUT/gpu.log" 2>&1
echo "gpu rc=$?"; tail -3 "$OUT/gpu.log"
grep -v '^#' "$OUT/gpu/timing.txt"; grep -h 'RT timers' "$OUT/gpu.log"
