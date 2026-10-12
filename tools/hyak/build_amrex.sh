#!/bin/bash
# Build NewRadHydroCode with the AMReX RT backend (rt_backend = amrex,
# docs/stage_d.md) on Hyak klone, and run the `amrex` test label. Run it on
# a compute node:
#
#   salloc -A astro -p cpu-g2 -c 32 --mem=64G --time=3:00:00
#   tools/hyak/build_amrex.sh [build_dir]
#
# or as a batch job from the repository root (sbatch tools/hyak/build_amrex.sh).
# AMReX 26.10 is cloned and built once into $RADRAY_SW/amrex-26.10 (3D, MPI,
# OpenMP, particles in single precision) unless AMREX_DIR points at an
# install. The MPI comes from the klone ompi module (HYAK_MODULES).
#
# Environment overrides:
#   HYAK_MODULES   modules after env.sh (default "ompi/4.1.6-2"; whatever
#                  provides mpicxx/mpirun with gcc 13)
#   AMREX_DIR      an existing AMReX install (lib/cmake/AMReX inside)
#   CMAKE_EXTRA    extra cmake arguments
#   CTEST_LABELS   amrex (default), none, or a label regex
#   RADRAY_MPIEXEC launcher for the tests (default: mpirun; srun on a job)
#
#SBATCH --job-name=radray-amrex-build
#SBATCH --account=astro
#SBATCH --partition=cpu-g2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --output=%x-%j.out
set -euo pipefail
if [ -n "${SLURM_SUBMIT_DIR:-}" ] && [ -f "$SLURM_SUBMIT_DIR/tools/hyak/env.sh" ]; then
  src="$SLURM_SUBMIT_DIR"
else
  src="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
fi
bld="${1:-$src/build-amrex}"
source "$src/tools/hyak/env.sh"
radray_module_load ${HYAK_MODULES:-ompi/4.1.6-2}
ncpu="${SLURM_CPUS_PER_TASK:-8}"
ulimit -s unlimited
echo "host $(hostname), $ncpu cpus, $(date)"
g++ --version | head -1; mpicxx --version | head -1; cmake --version | head -1

AMREX_DIR="${AMREX_DIR:-$RADRAY_SW/amrex-26.10}"
if [ ! -f "$AMREX_DIR/lib/cmake/AMReX/AMReXConfig.cmake" ]; then
  echo "building AMReX 26.10 into $AMREX_DIR"
  mkdir -p "$RADRAY_SW/src"
  [ -d "$RADRAY_SW/src/amrex" ] || git clone --branch 26.10 --depth 1 https://github.com/AMReX-Codes/amrex.git "$RADRAY_SW/src/amrex"
  cmake -S "$RADRAY_SW/src/amrex" -B "$RADRAY_SW/src/build-amrex" -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$AMREX_DIR" -DAMReX_SPACEDIM=3 -DAMReX_MPI=ON -DAMReX_OMP=ON \
    -DAMReX_PARTICLES=ON -DAMReX_PARTICLES_PRECISION=SINGLE -DAMReX_FORTRAN=OFF -DAMReX_FASTMATH=OFF \
    -DAMReX_LINEAR_SOLVERS=OFF -DAMReX_PROBINIT=OFF -DAMReX_AMRLEVEL=OFF
  cmake --build "$RADRAY_SW/src/build-amrex" -j "$ncpu"
  cmake --install "$RADRAY_SW/src/build-amrex"
fi

cmake -S "$src" -B "$bld" -DCMAKE_BUILD_TYPE=Release -DRADRAY_AMREX=ON \
  -DAMReX_DIR="$AMREX_DIR/lib/cmake/AMReX" \
  -DRADRAY_MPIEXEC="${RADRAY_MPIEXEC:-mpirun}" -DRADRAY_MPIEXEC_FLAGS="${RADRAY_MPIEXEC_FLAGS:---bind-to;none}" \
  ${CMAKE_EXTRA:-}
cmake --build "$bld" -j "$ncpu"

labels="${CTEST_LABELS:-amrex}"
case "$labels" in
  none) exit 0 ;;
  *)    ctest --test-dir "$bld" -L "$labels" -LE nightly -j "${CTEST_JOBS:-4}" --output-on-failure ;;
esac
