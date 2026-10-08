#!/usr/bin/env bash
# Lossy-precision runs for docs/ray_compression.md (step 6).
#
# usage: run_precision.sh OUTDIR [NPAR]
#
# Runs the RT test problems with the production ray settings
# (legacy_rt = legacy_physics = 0, nray_max_factor = 16) for a baseline and
# for each ray-quantisation variant (include/radray/ray_quant.hpp), then
# keeps only the fields the metrics need. Evaluate with
# precision_metrics.py OUTDIR.
set -euo pipefail
out=${1:?outdir}
npar=${2:-3}
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/../.." && pwd)
drv=${RT_DRIVER:-$root/build/rt_driver}
mkdir -p "$out"

# problem:nsteps
problems="uniform_1src:240 uniform_8src:80 lognormal_16src:100"
# variant name and extra --set arguments
variants=(
  "base:"
  "d1:ray_quant_bits=1"
  "d2:ray_quant_bits=2"
  "d3:ray_quant_bits=3"
  "d4:ray_quant_bits=4"
  "d6:ray_quant_bits=6"
  "d8:ray_quant_bits=8"
  "d10:ray_quant_bits=10"
  "r3:ray_quant_bits=3 ray_quant_dither=0"
  "r6:ray_quant_bits=6 ray_quant_dither=0"
  "struct:ray_quant_struct=1"
  "struct_d4:ray_quant_struct=1 ray_quant_bits=4"
  "nrmf32:nray_max_factor=32"
)

jobs=()
for pn in $problems; do
  p=${pn%%:*}; n=${pn##*:}
  for v in "${variants[@]}"; do
    name=${v%%:*}; sets=${v#*:}
    args=""
    for kv in $sets; do args="$args --set $kv"; done
    jobs+=("$p|$n|$name|$args")
  done
done

run_one() {
  IFS='|' read -r p n name args <<<"$1"
  d="$out/$p/$name"
  if [ -f "$d/done" ]; then return 0; fi
  mkdir -p "$d"
  # shellcheck disable=SC2086
  "$drv" "$root/tests/problems/$p" "$d" --nsteps "$n" --set legacy_rt=0 --set legacy_physics=0 \
      --set nray_max_factor=16 $args > "$d/run.log" 2>&1
  # keep fHI, temp, rhog and the stats
  find "$d" -name '*.bin' ! -name fHI.bin ! -name temp.bin ! -name rhog.bin -delete
  touch "$d/done"
  echo "done $p $name"
}
export -f run_one
export out drv root
printf '%s\n' "${jobs[@]}" | xargs -P "$npar" -I{} bash -c 'run_one "$@"' _ {}
