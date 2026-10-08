#!/bin/bash
# One-off software setup on Hyak klone (run inside a Slurm job, not on the
# login node):
#
#   sbatch tools/hyak/setup_sw.sh
#
# (Already done for user mcquinn: /gscratch/astro/mcquinn/sw.)
#
# 1. FFTW 3.3.10, single precision (libfftw3f), SSE2, static, built with gcc 13
#    into $RADRAY_SW/fftw-3.3.10-float. The klone fftw/3.3.10 module is
#    double precision only.
# 2. A conda env with Python, numpy, matplotlib (and healpy, scipy),
#    CMake (the klone cmake modules pull in an older gcc) and ffmpeg in
#    $RADRAY_SW/conda-envs/radray.
#
# The code links only the serial libfftw3f (the FFTs are a small part of the
# PM step), so FFTW is built without --enable-openmp.
#
#SBATCH --job-name=radray-setup
#SBATCH --account=astro
#SBATCH --partition=cpu-g2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=slurm-%x-%j.out
set -euo pipefail

repo="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
source "$repo/tools/hyak/env.sh"
ncpu="${SLURM_CPUS_PER_TASK:-4}"

# --- FFTW (single precision) ------------------------------------------------
if [ ! -f "$RADRAY_FFTW/lib/libfftw3f.a" ]; then
  src="$RADRAY_SW/src"
  mkdir -p "$src"
  cd "$src"
  [ -f fftw-3.3.10.tar.gz ] || curl -fsSLO https://www.fftw.org/fftw-3.3.10.tar.gz
  rm -rf fftw-3.3.10 && tar xzf fftw-3.3.10.tar.gz
  cd fftw-3.3.10
  # SSE2 codelets only. FFTW_ESTIMATE plans differ in the last bits with
  # the SIMD set; only an SSE2-only build reproduces the stored references
  # (AVX2, no-SIMD and Debian-configured builds do not; docs/cluster.md §2).
  ./configure --prefix="$RADRAY_FFTW" --enable-float --enable-sse2 CC=gcc CFLAGS="-O2"
  make -j"$ncpu"
  make install
fi
ls -l "$RADRAY_FFTW/lib/libfftw3f.a"

# --- Python env ---------------------------------------------------------------
if [ ! -x "$RADRAY_CONDA_ENV/bin/cmake" ]; then
  conda create -y -p "$RADRAY_CONDA_ENV" -c conda-forge --override-channels \
        python=3.12 numpy matplotlib scipy healpy cmake ffmpeg
fi
[ -x "$RADRAY_CONDA_ENV/bin/ffmpeg" ] || \
  conda install -y -p "$RADRAY_CONDA_ENV" -c conda-forge --override-channels ffmpeg
source "$repo/tools/hyak/env.sh"   # now activates the env
gcc --version | head -1; cmake --version | head -1; ffmpeg -version | head -1
python -c "import numpy, matplotlib, healpy; print('numpy', numpy.__version__, 'matplotlib', matplotlib.__version__, 'healpy', healpy.__version__)"
