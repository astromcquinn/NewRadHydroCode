# reference/: the Fortran specification and its stored outputs

The C++ port is validated bit for bit against Hy Trac's RadHydro and the
gfortran harnesses built from it (`CLAUDE.md`, rule 2). Since 2026-10-26 the
Fortran itself lives in the **private** repository
`astromcquinn/RadHydro-fortran`, carried here as the git submodule
`reference/fortran` (Hy Trac's condition for this repository being public:
his code stays private). What is here:

| Path | What |
|---|---|
| `fortran/` | the submodule: `Radhydro/` (the untouched original), `harness/`, `harness_hydro/`, `harness_pm/`, `harness_full/` (the gfortran harnesses: patched sources, Makefiles, run scripts, `PATCHES.md`), `harness_common/`, `tests/` (the Fortran test helpers). Empty without access. |
| `harness/outputs*/` | stored RT reference outputs (stats for all steps, fields at the last step) for the legacy, `legacy_rt = 0`, `legacy_rt = legacy_physics = 0` and `c_reduced = 0.1` sets |
| `harness_hydro/outputs/`, `harness_pm/outputs/`, `harness_full/outputs/` | stored hydro, PM and full-loop reference outputs |

The stored outputs are data the harnesses produced; the CTests compare the
C++ against them on every build (`ref_*`, `pm_ref_*`, `hydro_ref_*`,
`full_ref_*`). The live comparisons that build and run the harnesses
(labels `fortran`) need the submodule:

```
git submodule update --init reference/fortran        # with access to the private repository
cmake -S . -B build && cmake --build build && ctest --test-dir build -L fortran
```

`-DRADRAY_FORTRAN_DIR=<path>` points CMake at a checkout elsewhere. CI
populates the submodule from the repository secret `FORTRAN_TOKEN` (a
fine-grained token with read access to `RadHydro-fortran`); without it the
live comparisons are skipped and the stored-reference comparisons still
guard the port. The regeneration workflow (`regen_references.yml`) runs the
harnesses from the submodule and stores the outputs here.
