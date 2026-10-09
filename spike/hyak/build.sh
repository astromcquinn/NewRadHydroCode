#!/bin/bash
# Build the spike on Hyak (CPU: MPI + OpenMP, AVX-512 for Zen 4).
#   source spike/hyak/env.sh && bash spike/hyak/build.sh
# Run it on a compute node (salloc -A astro -p cpu-g2 -c 32 --mem 64G) or
# inside the benchmark job: compiling AMReX + Parthenon + Kokkos takes
# ~15 min on 4 cores, ~3 min on 32.
set -euo pipefail
: "${SPIKE_REPO:?source spike/hyak/env.sh first}"
J=${J:-$(nproc)}
ARCHFLAGS=${ARCHFLAGS:--march=znver4}

cd "$SPIKE_DEPS"
[ -d amrex ] || git clone -q --depth 1 --branch "$AMREX_TAG" https://github.com/AMReX-Codes/amrex.git amrex
[ -d parthenon ] || git clone -q --recursive --depth 1 --branch "$PARTHENON_TAG" \
  https://github.com/parthenon-hpc-lab/parthenon.git parthenon

# AMReX: MPI + OpenMP, particles (single precision), linear solvers
if [ ! -f "$SPIKE_DEPS/inst-amrex/lib/cmake/AMReX/AMReXConfig.cmake" ]; then
  cmake -S amrex -B build-amrex -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_FLAGS="$ARCHFLAGS" \
    -DAMReX_MPI=ON -DAMReX_OMP=ON -DAMReX_PARTICLES=ON -DAMReX_SPACEDIM=3 \
    -DAMReX_PARTICLES_PRECISION=SINGLE -DAMReX_FORTRAN=OFF -DAMReX_LINEAR_SOLVERS=ON \
    -DCMAKE_INSTALL_PREFIX="$SPIKE_DEPS/inst-amrex"
  cmake --build build-amrex -j "$J"
  cmake --install build-amrex
fi

# Spike projects. The flags must keep -ffp-contract=off (set by the
# projects); -march=znver4 only changes vector width, not results.
cmake -S "$SPIKE_REPO/spike/serial" -B "$SPIKE_BUILD/serial" -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_FLAGS="$ARCHFLAGS"
cmake --build "$SPIKE_BUILD/serial" -j "$J"

cmake -S "$SPIKE_REPO/spike/amrex" -B "$SPIKE_BUILD/amrex" -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_FLAGS="$ARCHFLAGS" -DAMReX_DIR="$SPIKE_DEPS/inst-amrex/lib/cmake/AMReX"
cmake --build "$SPIKE_BUILD/amrex" -j "$J"

cmake -S "$SPIKE_REPO/spike/parthenon" -B "$SPIKE_BUILD/parthenon" -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_FLAGS="$ARCHFLAGS" -DPARTHENON_SOURCE_DIR="$SPIKE_DEPS/parthenon" \
  -DKokkos_ENABLE_OPENMP=ON -DKokkos_ARCH_ZEN4=ON
cmake --build "$SPIKE_BUILD/parthenon" -j "$J"
echo "build.sh: done; binaries in $SPIKE_BUILD/{serial,amrex,parthenon}"
