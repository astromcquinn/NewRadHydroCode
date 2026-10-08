# Source this on Hyak klone (login or compute node) before building or
# running the spike:   source spike/hyak/env.sh
#
# Module names differ between Hyak software stacks. Override with e.g.
#   export SPIKE_MODULES="gcc/13.2.0 ompi/4.1.6 cmake/3.28.3"
# before sourcing; check with `module avail gcc ompi cmake cuda`.
# The defaults below load whatever version the stack marks as default.
: "${SPIKE_MODULES:=gcc ompi cmake}"
if command -v module > /dev/null 2>&1; then
  # Lmod's shell functions fail under `set -u` (the sbatch scripts use it)
  case $- in *u*) _spike_u=1; set +u;; *) _spike_u=;; esac
  for m in $SPIKE_MODULES; do
    module load "$m" || echo "env.sh: could not load module $m (set SPIKE_MODULES)" >&2
  done
  [ -n "$_spike_u" ] && set -u; unset _spike_u
fi
# klone's ompi/4.1.6-2 module sets CC=mpiCC (the C++ wrapper), which breaks
# CMake's C compiler check
[ "${CC:-}" = mpiCC ] && export CC=mpicc
# Optional directories to prepend to PATH, e.g. one with a newer cmake and a
# python3 with numpy (klone's cmake modules depend_on gcc/11.2.0 and its
# system python3 is 3.6). On klone (2026-10) we used
#   SPIKE_MODULES="gcc/13.2.0 ompi/4.1.6-2"
#   SPIKE_PATH=/gscratch/astro/$USER/spike/bin   (symlinks to cmake 4.4, python 3.12)
[ -n "${SPIKE_PATH:-}" ] && export PATH="$SPIKE_PATH:$PATH"
for t in mpicxx mpirun cmake g++ python3; do
  command -v $t > /dev/null 2>&1 || echo "env.sh: $t not found in PATH" >&2
done

# Where the repository lives and where dependencies/builds go. The
# dependencies need ~2 GB; /gscratch is the right place on klone.
SPIKE_REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
: "${SPIKE_WORK:=/gscratch/astro/$USER/spike}"
export SPIKE_REPO SPIKE_WORK
export SPIKE_DEPS="$SPIKE_WORK/deps"
export SPIKE_BUILD="$SPIKE_WORK/build"
mkdir -p "$SPIKE_DEPS" "$SPIKE_BUILD"

# Pinned framework versions (the ones validated in spike/*/README.md)
export AMREX_TAG=26.10
export PARTHENON_TAG=v25.12

# cpu-g2: AMD EPYC Genoa, 2 x 96 cores, no SMT. Pin threads to cores.
export OMP_PROC_BIND=close
export OMP_PLACES=cores
