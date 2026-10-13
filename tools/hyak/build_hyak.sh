#!/bin/bash
# Build NewRadHydroCode on Hyak klone (Rocky 8, AMD Genoa / Zen 4) and run
# the fast test suite. Run it on a compute node, not on the login node:
#
#   salloc -A astro -p cpu-g2 -c 32 --mem=64G --time=2:00:00
#   tools/hyak/build_hyak.sh [build_dir]
#
# or as a batch job from the repository root:
#
#   sbatch tools/hyak/build_hyak.sh [build_dir]
#
# The environment (gcc/13.2.0 module, SSE2-only FFTW 3.3.10 float, and a
# conda env with python, numpy, matplotlib, cmake and ffmpeg) comes from
# tools/hyak/env.sh; tools/hyak/setup_sw.sh installs it once. The klone
# fftw module is double precision only. See docs/cluster.md.
#
# Environment overrides:
#   HYAK_MODULES   extra modules to load after env.sh
#   CMAKE_EXTRA    extra cmake arguments, e.g. -DRADRAY_NIGHTLY_TESTS=ON
#   CTEST_LABELS   fast (default: -LE "long|nightly|platform_ref"), none, or a
#                  label regex for ctest -L (long, nightly, ...). The fast
#                  set leaves out the 10 stored RT references that differ
#                  on klone's libm (label platform_ref, docs/cluster.md
#                  section 3); CTEST_LABELS=platform_ref runs just them.
#   CTEST_JOBS     parallel tests (default: the cores of the job)
#
#SBATCH --job-name=radray-build
#SBATCH --account=astro
#SBATCH --partition=cpu-g2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=%x-%j.out
set -euo pipefail
# Under sbatch the script runs from a spool copy: use the submit directory.
if [ -n "${SLURM_SUBMIT_DIR:-}" ] && [ -f "$SLURM_SUBMIT_DIR/tools/hyak/env.sh" ]; then
  src="$SLURM_SUBMIT_DIR"
else
  src="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
fi
bld="${1:-$src/build}"
source "$src/tools/hyak/env.sh"
if [ -n "${HYAK_MODULES:-}" ]; then radray_module_load $HYAK_MODULES; fi
ncpu="${SLURM_CPUS_PER_TASK:-8}"
# The RT harness and the live Fortran comparisons need a large stack.
ulimit -s unlimited

echo "host $(hostname), $ncpu cpus, $(date)"
g++ --version | head -1; gfortran --version | head -1; cmake --version | head -1
python3 -c "import sys, numpy, matplotlib; print('python', sys.version.split()[0], 'numpy', numpy.__version__, 'matplotlib', matplotlib.__version__)"

# Default flags only: CMakeLists.txt sets -ffp-contract=off, and the
# bit-for-bit validation (ctest) is defined for these flags.
cmake -S "$src" -B "$bld" -DCMAKE_BUILD_TYPE=Release ${CMAKE_EXTRA:-}
cmake --build "$bld" -j "$ncpu"

labels="${CTEST_LABELS:-fast}"
case "$labels" in
  none) exit 0 ;;
  fast) sel=(-LE "long|nightly|platform_ref") ;;
  *)    sel=(-L "$labels") ;;
esac
ctest --test-dir "$bld" -j "${CTEST_JOBS:-$ncpu}" "${sel[@]}" --output-on-failure --timeout 7200
echo "done $(date)"
