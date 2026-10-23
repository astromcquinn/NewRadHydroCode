#!/bin/bash
# Build MUSIC (cosmo-sims/MUSIC, "MUSIC2": the MUSIC-1 code modernised, with
# zoom initial conditions) on klone, for the zoom ICs of
# docs/zoom_turbulence.md section 4.1 and the C26 boxes of section 3b.
#
#   source tools/hyak/env.sh
#   tools/hyak/setup_music.sh            # ~10 min on a login node
#
# Installs into $RADRAY_SW (default /gscratch/astro/$USER/sw):
#   $RADRAY_SW/src/MUSIC        the sources (git clone)
#   $RADRAY_SW/music/MUSIC      the binary
# and adds to the radray conda env (tools/hyak/setup_sw.sh): gsl, hdf5, the
# double-precision FFTW with threads (MUSIC needs double; the float FFTW of
# env.sh is for the code), pkg-config and camb (for the transfer tables,
# tools/ics/music/camb_transfer.py).
#
# The same cmake recipe builds and runs in the cloud container with the
# system GSL/HDF5/FFTW (2026-10-22); on klone the libraries come from conda,
# so CMAKE_PREFIX_PATH points at the env. ENABLE_PANPHASIA is off (the
# Panphasia generator is not needed; MUSIC's own hierarchical generator
# gives the per-level seeds of the ladder).
set -eu
RADRAY_SW="${RADRAY_SW:-/gscratch/astro/$USER/sw}"
MUSIC_SRC="${MUSIC_SRC:-$RADRAY_SW/src/MUSIC}"
MUSIC_DIR="${MUSIC_DIR:-$RADRAY_SW/music}"
MUSIC_URL="${MUSIC_URL:-https://github.com/cosmo-sims/MUSIC.git}"
ncpu="${NCPU:-8}"

[ -n "${CONDA_PREFIX:-}" ] || { echo "source tools/hyak/env.sh first (the radray conda env must be active)"; exit 1; }
echo "== dependencies into $CONDA_PREFIX"
conda install -y -c conda-forge gsl hdf5 "fftw=*=*openmp*" pkg-config camb >/dev/null || \
  conda install -y -c conda-forge gsl hdf5 fftw pkg-config camb
python -c "import camb; print('camb', camb.__version__)"

echo "== sources into $MUSIC_SRC"
if [ -d "$MUSIC_SRC/.git" ]; then
  git -C "$MUSIC_SRC" pull -q --ff-only
else
  mkdir -p "$(dirname "$MUSIC_SRC")"
  git clone -q "$MUSIC_URL" "$MUSIC_SRC"
fi
git -C "$MUSIC_SRC" log -1 --format='MUSIC at %h (%ad)' --date=short

echo "== build"
mkdir -p "$MUSIC_SRC/build" "$MUSIC_DIR"
cmake -S "$MUSIC_SRC" -B "$MUSIC_SRC/build" -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="$CONDA_PREFIX" -DFFTW3_ROOT="$CONDA_PREFIX" -DHDF5_ROOT="$CONDA_PREFIX" \
  -DENABLE_PANPHASIA=OFF
cmake --build "$MUSIC_SRC/build" -j "$ncpu"
cp "$MUSIC_SRC/build/MUSIC" "$MUSIC_DIR/MUSIC"
echo "== $MUSIC_DIR/MUSIC built; test:"
echo "   cd \$RUNS && mkdir -p ics_test && cd ics_test && cp $PWD/tools/ics/music/test32.conf . &&"
echo "   OMP_NUM_THREADS=8 $MUSIC_DIR/MUSIC test32.conf    # writes test32.hdf5 (generic format)"
