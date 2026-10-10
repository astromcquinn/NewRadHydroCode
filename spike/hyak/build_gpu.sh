#!/bin/bash
# GPU builds of the prototypes on a Hyak ckpt GPU node (A100: sm_80,
# H200/H100: sm_90). Needs a CUDA module >= 12.2 (AMReX 26.x requirement).
#   salloc -A astro -p ckpt-g2 --gpus=1 -c 16 --mem=64G
#   SPIKE_MODULES="gcc ompi cmake cuda" source spike/hyak/env.sh
#   CUDA_ARCH=80 bash spike/hyak/build_gpu.sh
# Status (2026-10-10): rt_parthenon/hydro_parthenon run on A100 (sm_80) and
# H200 (sm_90) with CUDA 12.4.1; validation and timings in
# docs/framework_decision.md section 4. The AMReX prototype has no CUDA
# target yet (only the AMReX library is built here).
set -euo pipefail
: "${SPIKE_REPO:?source spike/hyak/env.sh first}"
CUDA_ARCH=${CUDA_ARCH:-80}
J=${J:-$(nproc)}
KARCH=$([ "$CUDA_ARCH" = 90 ] && echo HOPPER90 || echo AMPERE80)
cd "$SPIKE_DEPS"
if [ ! -f "$SPIKE_DEPS/inst-amrex-cuda/lib/cmake/AMReX/AMReXConfig.cmake" ]; then
  # CC/CXX=gcc/g++ (not the MPI wrappers the ompi module sets): FindMPI then
  # passes the MPI include dirs explicitly, which nvcc needs (mpi.h not found otherwise)
  CC=gcc CXX=g++ cmake -S amrex -B build-amrex-cuda -DCMAKE_BUILD_TYPE=Release -DAMReX_GPU_BACKEND=CUDA \
    -DAMReX_CUDA_ARCH="${CUDA_ARCH:0:1}.${CUDA_ARCH:1}" -DAMReX_MPI=ON -DAMReX_OMP=OFF \
    -DAMReX_PARTICLES=ON -DAMReX_SPACEDIM=3 -DAMReX_PARTICLES_PRECISION=SINGLE -DAMReX_FORTRAN=OFF \
    -DCMAKE_INSTALL_PREFIX="$SPIKE_DEPS/inst-amrex-cuda"
  cmake --build build-amrex-cuda -j "$J"   # (not `a && b`: set -e ignores a failure of a)
  cmake --install build-amrex-cuda
fi
# The spike's AMReX project has no CUDA target yet (it needs
# enable_language(CUDA) and setup_target_for_cuda_compilation, which could
# not be tested here: AMReX 26.x refuses nvcc < 12.2). Building the AMReX
# CUDA library above is the first check; adding the target is a short
# Stage D task.
cmake -S "$SPIKE_REPO/spike/parthenon" -B "$SPIKE_BUILD/parthenon-cuda" -DCMAKE_BUILD_TYPE=Release \
  -DPARTHENON_SOURCE_DIR="$SPIKE_DEPS/parthenon" -DKokkos_ENABLE_CUDA=ON -DKokkos_ARCH_$KARCH=ON \
  -DCMAKE_CXX_COMPILER="$SPIKE_DEPS/parthenon/external/Kokkos/bin/nvcc_wrapper"
cmake --build "$SPIKE_BUILD/parthenon-cuda" -j "$J"
