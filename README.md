# NewRadHydroCode

A C++ cosmological radiative-transfer (RT) code built on the adaptive
ray-tracing algorithm of Hy Trac's RadHydro (Trac & Cen 2007; Trac, Cen &
Loeb 2008). The roadmap is in [docs/PLAN.md](docs/PLAN.md).

**Status (2026-10-09).** The port is complete and bit-identical to the
gfortran references: `apps/radhydro` runs the full `radhydro.f90` loop
(hydro + PM gravity + on-the-fly FoF halos and sources + RT), and the RT,
hydro and PM parts have their own drivers and Fortran harnesses. On top of
the port: optional RT and physics fixes behind `legacy_rt` / `legacy_physics`,
a reduced speed of light, OpenMP with thread-count-independent results,
restarts, an unsplit MUSCL-Hancock + HLLC hydro scheme with dual energy
(`hydro_scheme = unsplit`), a hydro test suite, a ray-compression study and
the AMReX-vs-Parthenon framework spike. What is done, in progress and next
is in [docs/STATUS.md](docs/STATUS.md); the decisions are in
[docs/PLAN.md](docs/PLAN.md) and the long-term plan in
[docs/ROADMAP.md](docs/ROADMAP.md).

**Phase 1, the faithful port.** `radray` is a C++17 port of
`reference/Radhydro`. On every validation problem it reproduces the gfortran
reference harnesses **bit for bit**: all output fields are byte-identical
and all ray and RT statistics are identical. This holds for the original
code and with the optional RT and physics fixes (`legacy_rt`,
`legacy_physics`) switched on in both codes.

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
  hydro.hpp, pm.hpp, fft.hpp   hydro_tools.f90, pm_tools.f90, fft_tools.f90 (FFTW)
  kernels/hydro_mh.hpp         unsplit MUSCL-Hancock + HLLC kernel (hydro_scheme = unsplit)
  kernels/ray_kernels.hpp, rt_cell.hpp, healpix_hd.hpp
                               the ray and RT-cell kernels, shared by the serial port and
                               the AMReX backend (docs/stage_d.md)
  rt_backend.hpp               rt_backend = serial | amrex, ray_encoding = float | rstruct
  halos.hpp, halo_sources.hpp  FoF halo finder, source model, legacy catalogs
  checkpoint.hpp  restart files
  simulation.hpp  all state, in explicit structs
  io.hpp          problem input, validation output
src/              implementation; src/amrex/ the AMReX RT backend (rt_level.cpp, the
                  ray encodings; built with -DRADRAY_AMREX=ON, a stub otherwise)
apps/             rt_driver (RT only), hydro_driver, pm_driver, radhydro (full
                  loop), find_halos (FoF on a snapshot)
tests/            unit tests, end-to-end checks, comparison with the references
  problems*/      validation problems per component; hydro_suite/ (docs/hydro_tests.md)
reference/        original Fortran (read-only) and the gfortran harnesses
                  (harness: RT, harness_hydro, harness_pm, harness_full, harness_common)
tools/            compare/ (field and stats comparison), ics/ (Zel'dovich ICs),
                  analysis/ (figures), production/ (production parameters, timing),
                  movies/ (stills), ray_compression/ (ray dumps, codec, precision),
                  hyak/ (UW Hyak klone environment, build and job scripts; docs/cluster.md)
spike/            framework spike: AMReX and Parthenon prototypes (docs/framework_decision.md);
                  serial/rt_serial is the bitwise reference of the AMReX backend
docs/             plan, status, roadmap, validation spec, port notes, design docs, figures
```

## Build

You need CMake 3.16 or newer, a C++17 compiler (tested with g++ 13.3) and
FFTW3 in single precision (`libfftw3-dev`; the PM gravity uses it, and both
codes use it in place of MKL). The live Fortran comparison tests need
`gfortran`; the Python tools need numpy (and matplotlib for figures).

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release     # -O2
cmake --build build -j
ctest --test-dir build --output-on-failure
```

The AMReX backend of the ray tracing (`rt_backend = amrex`, MPI + OpenMP,
`docs/stage_d.md`) is optional: `-DRADRAY_AMREX=ON` with
`-DAMReX_DIR=<install>/lib/cmake/AMReX` (AMReX 26.10 built 3-D with MPI,
OpenMP and particles in single precision), or `-DRADRAY_FETCH_AMREX=ON` to
download and build it as a subproject; it needs a C++20 compiler and MPI.
Its tests carry the label `amrex` (`ctest -L amrex`; `tools/hyak/build_amrex.sh`
on klone).

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
  [docs/figures/cooling_rates.png](docs/figures/cooling_rates.png)), and
  the converged photon limit of `photo_cell` (B17, `photon_conserve`).

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

Hydro, PM and full-loop tests carry the labels `hydro`, `pm` and `full`;
`threads` (bit-identical outputs for 1, 2 and 4 OpenMP threads), `restart`
(bit-identical restart) and `physics` (`ifront_rsla`, the reduced speed of
light). The four N = 32 full-loop reference tests (`full_ref_legacy_*`,
`full_ref_onthefly_*`, about 5 minutes each) and the ray-pool guard test
carry the label `nightly` and are disabled by default:

```
cmake -DRADRAY_NIGHTLY_TESTS=ON build && ctest --test-dir build -L nightly
```

Every test runs with `OMP_NUM_THREADS=2` unless it sets the thread count
itself.

## Full radiation-hydrodynamics run

```
./build/radhydro tests/problems_full/onthefly out/onthefly [--set key=value ...] [--quiet]
reference/harness_full/run_problem.sh tests/problems_full/onthefly out/f_onthefly   # Fortran
python3 tools/compare/compare.py out/f_onthefly out/onthefly --rtol 0 --atol-frac 0
```

Parameters, outputs, validation and the cost model are in
[docs/full_run.md](docs/full_run.md). OpenMP (`OMP_NUM_THREADS`, results
independent of the thread count) is described in
[docs/parallel.md](docs/parallel.md); the reduced speed of light
(`c_reduced`), the deterministic parallel ray tracing, restarts, movie
frames and the production-run cost in
[docs/production_run.md](docs/production_run.md).

Running on the UW Hyak cluster (`klone`): environment, build, tests,
submitting and restarting jobs, OpenMP scaling on one node and the
measured runs are in [docs/cluster.md](docs/cluster.md); the scripts are
in `tools/hyak/`.

How the Fortran maps onto the C++, the ambiguous constructs and the bugs
found in the original are described in
[docs/port_notes.md](docs/port_notes.md).
