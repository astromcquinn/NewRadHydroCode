# NewRadHydroCode

A C++ cosmological radiative-transfer (RT) code built on the adaptive
ray-tracing algorithm of Hy Trac's RadHydro (Trac & Cen 2007; Trac, Cen &
Loeb 2008). The roadmap is in [docs/PLAN.md](docs/PLAN.md).

**Status: phase 1, the faithful port.** `radray` is a serial C++17 port of the
RT part of `reference/Radhydro`. That covers ray tracing, HEALPix ray
splitting and merging, photo-ionisation and heating chemistry, and the
cosmology and units the RT uses. On every validation problem it reproduces
the gfortran reference harness (`reference/harness`) **bit for bit**: all
output fields are byte-identical and all ray and RT statistics are identical.
This holds for the original code and with the optional RT and physics fixes
(`legacy_rt`, `legacy_physics`) switched on in both codes.

## Layout

```
include/radray/   headers (namespace radray)
  fortran.hpp     Fortran-semantics helpers (INT/NINT/MOD, real(4) sums, x**n)
  constants.hpp   physical constants, cosmo/RT parameters (exact Fortran values)
  healpix.hpp     HEALPix NESTED pix2vec/vec2pix (port of HEALPix 3.11 Fortran)
  cosmology.hpp   calc_sim_var, calc_dt_from_da, time stepping
  rt_rates.hpp    atomic rates, photoionisation cross sections
  grid.hpp        RT grid / gas grid (structure of arrays)
  rays.hpp        ray storage pool (ray_index, free buffers)
  ray_tracing.hpp ray_tracing_tools.f90
  rt.hpp          rt_tools.f90 (+ domains, sources)
  simulation.hpp  all state, in explicit structs
  io.hpp          problem input, validation output
src/              implementation
apps/rt_driver.cpp   RT-only driver implementing docs/validation_spec.md
tests/            unit tests, end-to-end checks, comparison with the reference
  problems/       validation problems (made by reference/harness/make_test_problem.py)
reference/        original Fortran (read-only) and the gfortran harness
tools/compare/    field and stats comparison tool
tools/plot_cooling.py  cooling-function figure (docs/figures/cooling_rates.png)
docs/             plan, validation spec, design notes, port notes
```

## Build

You need CMake 3.16 or newer and a C++17 compiler (tested with g++ 13.3). It
has no other dependencies.

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release     # -O2
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Other options:

* `-DRADRAY_SANITIZE=ON` builds with `-fsanitize=address,undefined`.
* The build always uses `-ffp-contract=off -fno-fast-math`. Bitwise
  agreement with the Fortran needs plain IEEE arithmetic: no FMA contraction
  and no reassociation. Do not add `-march=native` or `-ffast-math`.

## Run

```
./build/rt_driver tests/problems/uniform_1src out/uniform_1src [--nsteps n] [--set key=value ...]
```

`--set` overrides a `params.txt` entry. Two switches select fixes to the
original code; both default to 1, the original behaviour:

* `--set legacy_rt=0` turns on the RT bookkeeping fixes
  ([docs/port_notes.md](docs/port_notes.md) §5a);
* `--set legacy_physics=0` turns on the atomic-physics fixes: case B
  recombination cooling, HeI collisional excitation, HeII dielectronic
  recombination, the HeIII case B rate and collisional-ionisation cooling
  (§5b, with the cooling function in
  [docs/figures/cooling_rates.png](docs/figures/cooling_rates.png)).

The recommended configuration is `--set legacy_rt=0 --set legacy_physics=0`.

The driver reads `params.txt`, `gas.bin` and `sources.txt`. It then runs
`nsteps` hydro steps. Each step is: clock update, source assignment, then
`radiative_transfer`. After each step it writes
`out/<problem>/step_NNNN/{<field>.bin, stats.txt, rt_stats.txt}`, plus
`ray.log`. The formats are in [docs/validation_spec.md](docs/validation_spec.md).

To run every problem and compare with Fortran reference outputs:

```
reference/harness/run_all.sh /tmp/ref              # Fortran reference (needs gfortran)
tests/run_all.sh build /tmp/cpp /tmp/ref           # C++ runs + comparison
```

Each N = 32 problem takes about 0.5–2 s (10 hydro steps × 12 RT substeps).

## Tests

| test | what |
|---|---|
| `fortran_semantics` | INT/NINT/FLOOR/MOD helpers. Every folded constant is checked against a gfortran dump (`tests/fortran_consts`). |
| `healpix` | pix2vec/vec2pix round trips and nesting, nside 1..8192 |
| `healpix_vs_healpy` | 1e6 random vectors: exact pixel agreement with healpy. pix2vec agrees to 1.1e-15. Skipped if healpy is missing. |
| `rays` | ray_segment, nside/radius relations, chain ops, split and merge conservation |
| `end_to_end` | all problems: I-front growth and symmetry, photon budget, ray stats, merge path |
| `rates` | the `legacy_physics = 0` rate coefficients against independent evaluations of the published fits (Hui & Gnedin 1997, Cen 1992) |
| `end_to_end_legacy0`, `end_to_end_fixed` | as `end_to_end`, with `legacy_rt=0`, and with `legacy_rt=0 legacy_physics=0` |
| `ref_<problem>` | byte comparison with the stored harness outputs `reference/harness/outputs` |
| `ref_legacy0_<problem>`, `ref_fixed_<problem>` | the same for `legacy_rt=0` (`outputs_legacy_rt0`) and for `legacy_rt=0 legacy_physics=0` (`outputs_legacy_rt0_phys0`) |

How the Fortran maps onto the C++, the ambiguous constructs and the bugs
found in the original are described in
[docs/port_notes.md](docs/port_notes.md).
