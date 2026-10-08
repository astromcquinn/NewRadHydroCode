#!/bin/bash
# Build NewRadHydroCode on Hyak klone (Rocky 8, AMD Genoa / Zen 4).
#   tools/hyak/build_hyak.sh [build_dir]
# Needs g++ (C++17), FFTW3 single precision (libfftw3f) and python3 with numpy
# and matplotlib (ffmpeg for the movie). Module names differ between
# clusters and change over time: set HYAK_MODULES to what `module avail`
# shows on klone, e.g.
#   HYAK_MODULES="gcc/13.2.0 fftw/3.3.10" tools/hyak/build_hyak.sh
# If FFTW is not a module, install it in user space (spack, conda
# fftw, or ./configure --enable-float --enable-openmp --prefix=$HOME/opt)
# and pass -DCMAKE_PREFIX_PATH=$HOME/opt via CMAKE_EXTRA.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
src="$here/../.."
bld="${1:-$src/build}"
if [ -n "${HYAK_MODULES:-}" ]; then module load $HYAK_MODULES; fi
# Default flags only: CMakeLists.txt sets -ffp-contract=off, and the
# bit-for-bit validation (ctest) is defined for these flags.
cmake -S "$src" -B "$bld" -DCMAKE_BUILD_TYPE=Release ${CMAKE_EXTRA:-}
cmake --build "$bld" -j "${SLURM_CPUS_PER_TASK:-8}"
ctest --test-dir "$bld" -j "${SLURM_CPUS_PER_TASK:-8}" -LE "long|nightly" | tail -3
