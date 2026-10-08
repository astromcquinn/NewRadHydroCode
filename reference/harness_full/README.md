# Fortran reference harness: the full radhydro.f90

This directory builds Hy Trac's complete RadHydro program (hydro + PM
gravity + halos + ray-tracing RT) with gfortran on one thread, with FFTW3f
instead of MKL, and runs it on the problems of `tests/problems_full`. Its
outputs are the reference for `apps/radhydro` (see `docs/full_run.md`).

```
reference/harness_full/
  src/radhydro_full.f90       radhydro.f90 + HARNESS_FULL patches (PATCHES.md 2)
  src/harness_full_tools.f90  parameters, on-the-fly sources (C interface),
                              compact output (PATCHES.md 3)
  src/cosmo_parameters.f90    box by the preprocessor (PATCHES.md 1)
  Makefile                    everything else comes from ../harness/src,
                              ../harness_hydro/src, ../harness_pm/src,
                              ../harness_common and ../Radhydro
  run_problem.sh              build for the problem's sizes and run one problem
  run_all.sh                  run (and --store) every validation case
  store_outputs.py            compact copy of a run -> outputs/<case>
  outputs/<case>/             stored reference outputs (CTest full_ref_*)
```

## Build and run

Needs gfortran (13.3 tested), FFTW3 single precision and, for the
on-the-fly source mode, the C++ libraries of the CMake build
(`cmake --build build --target radray_halos`).

```
make N=32 NDM=16 NRF=4 BOX=4          # build/N32_D16_R4_B4_F1/radhydro_full.x
make debug N=32 NDM=16 NRF=4 BOX=4    # -O0 -g -fcheck=all
./run_problem.sh ../../tests/problems_full/legacy runs/legacy [--debug] \
                 [--set legacy_rt=0 --set legacy_physics=0]
./run_all.sh [runs] [--store]         # all cases, 4 in parallel
```

`run_problem.sh` reads `N`, `Ndm`, `nray_max_factor` and `box` from the
problem (and from `--set`), builds the matching executable if needed, raises
the stack limit, runs, and deletes the `D2D/` slices.

## Outputs

`<outdir>` is also the original's working directory: `Log/hydro.log`
(`summarize_gas` records), `Log/rt.log`, `Log/ray.log`, `Pk/pk.<ztag>`,
`DM/xv_dm.<ztag>`, `checkpoint.<ztag>`, `timing.txt` as the original writes
them, plus the compact validation output `step_%04d/` (`stats.txt` every
step; `fHI fHeI fHeII temp rhog ngam nu fgam tau gamma lambda dx u_gas w_gas
v_gas xv_dm` `.bin` at the zdata1 output steps and the last step).

## Cases

| case | problem | flags | steps | RT substeps | Fortran / C++ time |
|---|---|---|---|---|---|
| smoke_rt1, smoke_rt0 | smoke (N=16) | rt1: legacy_rt = legacy_physics = 1; rt0: both 0 | 15 | 1479 | see docs/full_run.md |
| legacy_rt1, legacy_rt0 | legacy (N=32, legacy catalogs) | as above | 24 | 1669 | ~4 min |
| onthefly_rt1, onthefly_rt0 | onthefly (N=32, C++ FoF) | as above | 24 | 1669 | ~4 min |

All cases are bit-identical between this harness and `apps/radhydro`
(`compare.py --rtol 0 --atol-frac 0`, plus `hydro.log`, `Pk` and `DM` files
byte for byte).
