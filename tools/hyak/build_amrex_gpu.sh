#!/bin/bash
# CUDA build of the AMReX RT backend on a Hyak ckpt GPU node (A100: sm_80,
# H200/H100: sm_90), docs/stage_d.md section 5. Needs CUDA >= 12.2 (AMReX
# 26.x) and gcc 13. The AMReX CUDA library goes to $RADRAY_SW/amrex-26.10-cuda,
# the code to build-amrex-cuda. The ray kernels compile as CUDA
# (--fmad=false: no fused multiply-adds, as -ffp-contract=off on the host);
# the rate tables are in device memory. Compile-only here; the A100/H200
# runs are tools/hyak/amrex_gpu.slurm.
#
#   salloc -A astro -p ckpt-g2 --gpus=1 -c 16 --mem=64G
#   CUDA_ARCH=80 tools/hyak/build_amrex_gpu.sh
#
# Environment: CUDA_ARCH (80 | 90), HYAK_MODULES (default "ompi/4.1.6-2
# cuda/12.4.1"), CMAKE_EXTRA.
set -euo pipefail
src="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
source "$src/tools/hyak/env.sh"
module load ${HYAK_MODULES:-ompi/4.1.6-2 cuda/12.4.1}
CUDA_ARCH=${CUDA_ARCH:-80}
J=${J:-$(nproc)}
AMREX_DIR="${AMREX_DIR:-$RADRAY_SW/amrex-26.10-cuda}"
if [ ! -f "$AMREX_DIR/lib/cmake/AMReX/AMReXConfig.cmake" ]; then
  mkdir -p "$RADRAY_SW/src"
  [ -d "$RADRAY_SW/src/amrex" ] || git clone --branch 26.10 --depth 1 https://github.com/AMReX-Codes/amrex.git "$RADRAY_SW/src/amrex"
  # CC/CXX = gcc/g++ (not the MPI wrappers): FindMPI then passes the MPI
  # include dirs explicitly, which nvcc needs (spike/hyak/build_gpu.sh)
  CC=gcc CXX=g++ cmake -S "$RADRAY_SW/src/amrex" -B "$RADRAY_SW/src/build-amrex-cuda" -DCMAKE_BUILD_TYPE=Release \
    -DAMReX_GPU_BACKEND=CUDA -DAMReX_CUDA_ARCH="${CUDA_ARCH:0:1}.${CUDA_ARCH:1}" \
    -DAMReX_MPI=ON -DAMReX_OMP=OFF -DAMReX_PARTICLES=ON -DAMReX_SPACEDIM=3 \
    -DAMReX_PARTICLES_PRECISION=SINGLE -DAMReX_FORTRAN=OFF -DAMReX_FASTMATH=OFF \
    -DAMReX_LINEAR_SOLVERS=OFF -DAMReX_PROBINIT=OFF -DAMReX_AMRLEVEL=OFF \
    -DCMAKE_INSTALL_PREFIX="$AMREX_DIR"
  cmake --build "$RADRAY_SW/src/build-amrex-cuda" -j "$J"
  cmake --install "$RADRAY_SW/src/build-amrex-cuda"
fi
CC=gcc CXX=g++ cmake -S "$src" -B "$src/build-amrex-cuda" -DCMAKE_BUILD_TYPE=Release -DRADRAY_AMREX=ON \
  -DAMReX_DIR="$AMREX_DIR/lib/cmake/AMReX" -DCMAKE_CUDA_ARCHITECTURES="$CUDA_ARCH" \
  -DRADRAY_MPIEXEC=srun -DRADRAY_MPIEXEC_FLAGS="" ${CMAKE_EXTRA:-}
cmake --build "$src/build-amrex-cuda" -j "$J" --target rt_driver radhydro
echo "built $src/build-amrex-cuda/{rt_driver,radhydro} for sm_$CUDA_ARCH"
