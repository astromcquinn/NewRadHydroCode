# Environment for building and running NewRadHydroCode on UW Hyak klone.
#
#   source tools/hyak/env.sh
#
# Loads gcc 13, activates the Python env (which also provides CMake) and points
# CMake and the Fortran harness Makefiles at the single-precision FFTW built
# by tools/hyak/setup_sw.sh. See docs/cluster.md.

# Lmod's `module` and conda's activation scripts read unset variables, so
# `set -u` is relaxed while this file runs and restored at the end.
case $- in *u*) _radray_u=1; set +u ;; *) _radray_u=0 ;; esac

RADRAY_SW="${RADRAY_SW:-/gscratch/astro/$USER/sw}"
# The same SSE2 build with --enable-openmp (libfftw3f_omp, used only by
# pm_refine > 1, docs/pm_refine.md) when it is installed.
if [ -z "${RADRAY_FFTW:-}" ] && [ -f "$RADRAY_SW/fftw-3.3.10-float-omp/lib/libfftw3f_omp.a" ]; then
  RADRAY_FFTW="$RADRAY_SW/fftw-3.3.10-float-omp"
fi
RADRAY_FFTW="${RADRAY_FFTW:-$RADRAY_SW/fftw-3.3.10-float}"
RADRAY_CONDA_ENV="${RADRAY_CONDA_ENV:-$RADRAY_SW/conda-envs/radray}"
export RADRAY_SW RADRAY_FFTW RADRAY_CONDA_ENV

# gcc 13 (g++, gfortran): the stored references were made with gfortran 13.
# CMake comes from the conda env: every klone cmake module depends_on an
# older gcc (cmake/3.25.1 -> gcc/11.2.0) and would replace gcc 13.
module load gcc/13.2.0

# Keep conda's package cache out of the small home quota.
export CONDA_PKGS_DIRS="$RADRAY_SW/conda-pkgs"
# Miniforge from klone's conda/Miniforge3-25.9.1-0 module, without the
# module: it defines `conda` only in interactive shells and refuses to load
# when an env is already active, so sourcing this file twice would fail.
# Test for the shell function: a conda executable on an exported PATH (sbatch
# --export=ALL from a shell that sourced this file) cannot activate.
declare -F conda >/dev/null || source /mmfs1/sw/miniforge3/25.9.1-0/etc/profile.d/conda.sh
if [ -d "$RADRAY_CONDA_ENV" ]; then
  conda activate "$RADRAY_CONDA_ENV"
fi

# CMake: find_path(fftw3.h) / find_library(fftw3f).
export CMAKE_PREFIX_PATH="$RADRAY_FFTW${CMAKE_PREFIX_PATH:+:$CMAKE_PREFIX_PATH}"
# Harness Makefiles (FFTW_INC ?= /usr/include, FFTW_LIB ?= -lfftw3f).
export FFTW_INC="$RADRAY_FFTW/include"
export FFTW_LIB="-L$RADRAY_FFTW/lib -lfftw3f"

# Load modules from a script that runs under `set -u` (Lmod's init reads
# unset variables such as LD_PRELOAD and the load would fail). klone's
# ompi/4.1.6-2 module also sets CC=mpiCC, the C++ wrapper, which breaks
# CMake's C compiler check when AMReX is built.
radray_module_load() {
  case $- in *u*) _rml_u=1; set +u ;; *) _rml_u=0 ;; esac
  module load "$@"
  _rml_rc=$?
  [ "${CC:-}" = mpiCC ] && export CC=mpicc
  [ "$_rml_u" = 1 ] && set -u
  unset _rml_u
  return $_rml_rc
}

[ "$_radray_u" = 1 ] && set -u
unset _radray_u
