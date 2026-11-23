#!/bin/bash
# Uniform mean free path test of the rays (tools/analysis/ray_transfer.cpp,
# docs/rt_moments.md section 8.4): T(k) = sqrt(P_N / P_L) of the ray tracer
# for several ray budgets, beside M1 and OTVET with the same sources.
# Run inside an allocation, e.g.
#   sbatch --account=astro --partition=cpu-g2 --cpus-per-task=32 --mem=60G --time=4:00:00 \
#     --export=ALL,N=128,FLIST="8 15 64 120",MOMENTS=1,OUT=$RUNS/ray_transfer tools/hyak/ray_transfer.sh
# env: N       cells per side
#      FLIST   nray_max_factor values (15: the 512^3 suite; 64: the original)
#      FRACS   lambda = N / frac cells (default "4 1")
#      NSRC    source cells (default N^3 / 512: at 192 rays per source cell and
#              substep the pool cannot hold every cell as a source)
#      MOMENTS 1: also M1 and OTVET (moment_inject = 1) with the same sources
#      OUT     output directory
#      BUILD   build directory with libradray.a (default <repo>/build)
# Writes OUT/rays_N<N>_l<frac>_f<F>.{L,N,Navg,log}, OUT/{m1,otvet}_N<N>_l<frac>.{L,N}
# and the table OUT/T_rays_N<N>_l<frac>.txt (moment_transfer.py; the _avg
# columns are the mean over the last 16 substeps).
set -eu
W=${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}
[ -f "$W/tools/analysis/ray_transfer.cpp" ] || { echo "run from the repository root" >&2; exit 1; }
: "${N:?}" "${FLIST:?}" "${OUT:?}"
BUILD=${BUILD:-$W/build}
cd "$W"; unset RADRAY_FFTW; source tools/hyak/env.sh >/dev/null 2>&1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1} OMP_PLACES=cores OMP_PROC_BIND=close
NSRC=${NSRC:-$((N * N * N / 512))}
mkdir -p "$OUT"
tag=${SLURM_JOB_ID:-$$}
CXX="g++ -O2 -fopenmp -ffp-contract=off -Iinclude"
LIBS="$BUILD/libradray.a -L$RADRAY_FFTW/lib -lfftw3f_omp -lfftw3f"
$CXX tools/analysis/ray_transfer.cpp $LIBS -o "$OUT/ray_transfer_$tag"
[ "${MOMENTS:-0}" = 1 ] && $CXX tools/analysis/moment_transfer.cpp $LIBS -o "$OUT/moment_transfer_$tag"

# a uniform problem of N^3 cells: the gas state of tests/problems/uniform_1src
P=$OUT/problem_N$N
if [ ! -f "$P/gas.bin" ]; then
  mkdir -p "$P"
  python3 - "$W/tests/problems/uniform_1src/gas.bin" "$P/gas.bin.tmp" "$N" <<'PY'
import sys
import numpy as np
g = np.fromfile(sys.argv[1], dtype=np.float32).reshape(-1, 9)
assert np.all(g == g[0])
np.tile(g[0], (int(sys.argv[3]) ** 3, 1)).astype(np.float32).tofile(sys.argv[2])
PY
  mv "$P/gas.bin.tmp" "$P/gas.bin"
  printf 'N = %d\nz_init = 10.0\nda_sim = 0.0001\nnsteps = 1\n' "$N" > "$P/params.txt"
fi

for frac in ${FRACS:-4 1}; do
  lam=$((N / frac)); pre=""
  for f in $FLIST; do
    p=$OUT/rays_N${N}_l${frac}_f$f
    [ -f "$p.N" ] || "$OUT/ray_transfer_$tag" "$P" "$p" --lambda $lam --nsrc "$NSRC" --set nray_max_factor="$f" \
        --set legacy_rt=0 --set legacy_physics=0 --set deterministic_parallel=1 --set ray_rotate=3 \
        --set threads="$OMP_NUM_THREADS" > "$p.log" 2>&1
    ln -sf "$(basename "$p").L" "${p}_avg.L"; ln -sf "$(basename "$p").Navg" "${p}_avg.N"
    pre="$pre $p ${p}_avg"
  done
  if [ "${MOMENTS:-0}" = 1 ]; then
    steps=$((20 * lam)); [ $steps -ge 400 ] || steps=400
    for m in m1 otvet; do
      q=$OUT/${m}_N${N}_l${frac}
      [ -f "$q.N" ] || "$OUT/moment_transfer_$tag" "$N" $lam "$NSRC" $m $steps "$q" 1 1 4 2>/dev/null
      pre="$pre $q"
    done
  fi
  python3 tools/analysis/moment_transfer.py "$N" $lam $pre > "$OUT/T_rays_N${N}_l${frac}.txt"
done
rm -f "$OUT/ray_transfer_$tag" "$OUT/moment_transfer_$tag"
