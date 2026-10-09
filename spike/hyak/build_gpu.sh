#!/bin/bash
# GPU builds of the prototypes on a Hyak ckpt GPU node (A100: sm_80,
# H200/H100: sm_90). Needs a CUDA module >= 12.2 (AMReX 26.x requirement).
#   salloc -A astro -p ckpt-g2 --gpus=1 -c 16 --mem=64G
#   SPIKE_MODULES="gcc ompi cmake cuda" source spike/hyak/env.sh
#   CUDA_ARCH=80 bash spike/hyak/build_gpu.sh
# Status: the shared kernels compile as CUDA device code (checked with
# nvcc 12.0, sm_80) and rt_parthenon/hydro_parthenon build with Kokkos
# CUDA (compile only, nvcc 12.0); nothing has been RUN on a GPU yet. Known runtime blockers are listed
# in spike/amrex/README.md and spike/parthenon/README.md (host pointers to
# the rate tables in PhysCtx, host-side setup loops, per-block launches).
set -euo pipefail
: "${SPIKE_REPO:?source spike/hyak/env.sh first}"
CUDA_ARCH=${CUDA_ARCH:-80}
J=${J:-$(nproc)}
KARCH=$([ "$CUDA_ARCH" = 90 ] && echo HOPPER90 || echo AMPERE80)
cd "$SPIKE_DEPS"
cmake -S amrex -B build-amrex-cuda -DCMAKE_BUILD_TYPE=Release -DAMReX_GPU_BACKEND=CUDA \
  -DAMReX_CUDA_ARCH="${CUDA_ARCH:0:1}.${CUDA_ARCH:1}" -DAMReX_MPI=ON -DAMReX_OMP=OFF \
  -DAMReX_PARTICLES=ON -DAMReX_SPACEDIM=3 -DAMReX_PARTICLES_PRECISION=SINGLE -DAMReX_FORTRAN=OFF \
  -DCMAKE_INSTALL_PREFIX="$SPIKE_DEPS/inst-amrex-cuda"
cmake --build build-amrex-cuda -j "$J" && cmake --install build-amrex-cuda
# The spike's AMReX project has no CUDA target yet (it needs
# enable_language(CUDA) and setup_target_for_cuda_compilation, which could
# not be tested here: AMReX 26.x refuses nvcc < 12.2). Building the AMReX
# CUDA library above is the first check; adding the target is a short
# Stage D task.
cmake -S "$SPIKE_REPO/spike/parthenon" -B "$SPIKE_BUILD/parthenon-cuda" -DCMAKE_BUILD_TYPE=Release \
  -DPARTHENON_SOURCE_DIR="$SPIKE_DEPS/parthenon" -DKokkos_ENABLE_CUDA=ON -DKokkos_ARCH_$KARCH=ON \
  -DCMAKE_CXX_COMPILER="$SPIKE_DEPS/parthenon/external/Kokkos/bin/nvcc_wrapper"
cmake --build "$SPIKE_BUILD/parthenon-cuda" -j "$J"
