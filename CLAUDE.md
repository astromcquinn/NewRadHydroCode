# NewRadHydroCode — guide for Claude Code sessions

This repo holds a C++ port of Hy Trac's cosmological RadHydro, the adaptive
ray-tracing RT of Trac & Cen 2007 and Trac, Cen & Loeb 2008. It comes with
gfortran reference builds of the original Fortran. The goals are a faithful
port first, then MPI+GPU parallelism, then compact ray encodings and adaptive
RT.

At the start of a session, read:
- `docs/STATUS.md`: what's done, what's in progress, open questions.
- `docs/PLAN.md`: every user decision so far.
- `docs/ROADMAP.md`: the long-term plan (AMR, tree/multigrid gravity,
  OpenMP/MPI/GPU, ray compression), the proposed order, and suggested
  threads.

## Ground rules (non-negotiable)

1. **Never edit `reference/Radhydro/`.** It is the untouched original. Changes
   to the Fortran go only into the harness copies
   (`reference/harness*/src`) as documented patches in that harness's
   `PATCHES.md`.
2. **Bit-for-bit validation.**
   - The gfortran harnesses hold the original specification, for
     comparison. With the default (original) settings the C++ must stay
     byte-identical to them on every validation problem.
   - The existing mirrored flags (`legacy_rt`, `legacy_physics`) are also
     in the Fortran harnesses; keep the C++ byte-identical in every
     combination of them.
   - Check with `tools/compare/compare.py REF NEW --rtol 0 --atol-frac 0`.
     `--rtol 0` alone still allows absolute differences.
   - Changing results on purpose requires a C++ flag whose default keeps
     the original behaviour. New flags do **not** need a mirrored
     implementation in the Fortran harness (user decision, 2026-10-08):
     the Fortran only has to provide the original for comparison. Example:
     the research-only `ray_quant_*` options (`docs/ray_compression.md`).
3. **Legacy flags always reproduce the original.**
   - `legacy_rt = 1` and `legacy_physics = 1` (the defaults) give the
     original behaviour.
   - `= 0` turns on the RT bookkeeping fixes B1–B4, B8, B9 and the physics
     fixes P1–P10 (`docs/port_notes.md` §5a/§5b).
   - Every new correction to original behaviour goes behind a flag. Use
     `legacy_physics` for physics and `legacy_rt` for RT bookkeeping, unless
     the user says otherwise.
4. **Fortran semantics in C++.**
   - Match the precision exactly: real(4) is float, real(8) is double.
   - Fortran single-precision literals such as `0.99` and `1E-6` stay float.
   - Use the helpers in `include/radray/fortran.hpp`, including gfortran
     MAX/MIN tie and NaN behaviour (`f::gmax`/`f::gmin`).
   - Build with `-ffp-contract=off`, no fast-math.
5. **Determinism.**
   - Results must not depend on thread count.
   - The original is nondeterministic when multithreaded; it is validated
     with Ncpu = 1 only.
6. **Flag, don't silently "fix", original quirks.** Each one is listed with
   an ID in the port notes: B*, H*, PM, F* and halo quirks.

## Layout

| Path | What |
|---|---|
| `include/radray/`, `src/` | C++17 library, namespace `radray` |
| `apps/rt_driver.cpp` | RT only (`docs/validation_spec.md`) |
| `apps/hydro_driver.cpp` | hydro only |
| `apps/pm_driver.cpp` | PM gravity only |
| `apps/radhydro.cpp` | full `radhydro.f90` loop (`docs/full_run.md`) |
| `apps/find_halos.cpp` | FoF on a snapshot (`docs/halos_and_sources.md`) |
| `reference/harness` | gfortran RT reference |
| `reference/harness_hydro` | gfortran hydro reference |
| `reference/harness_pm` | gfortran PM reference |
| `reference/harness_full` | gfortran full-loop reference |
| `reference/harness_common` | FFTW-based MKL DFTI shim, shared timestep |
| `tests/problems*` | validation problems |
| `reference/*/outputs*` | stored reference outputs: stats for all steps, fields at the last step |
| `tools/compare` | field/stat comparison (with self-tests) |
| `tools/ics` | Zel'dovich ICs |
| `tools/analysis` | figures |
| `docs/` | plan, status, port notes per component, design doc, figures |

## Build and test

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build -j4 -LE "long|nightly"   # fast suite
ctest --test-dir build -L long                   # ~5 min each, full-loop references
```

- Needs `g++`, `gfortran` (for the live Fortran comparison tests), `libfftw3-dev`
  (FFTW3 single precision; MKL is replaced by FFTW in both codes) and python3
  with numpy (and matplotlib for figures).
- Each harness has a `Makefile`/`run_all.sh`; see its README.
- The RT harness needs `ulimit -s unlimited`, which its run scripts set.
- Runtime parameters: `params.txt` plus `--set key=value` on every driver.

## Working conventions

- **Branches:** a feature branch and a PR per workstream.
- **Commits:** end each commit message with the attribution lines the
  harness asks for.
- **Parallel agents:** use git worktrees, one per agent. Keep their edits to
  shared files (`CMakeLists.txt`, `tests/CMakeLists.txt`, `simulation.hpp`,
  `config.hpp`) minimal and additive. Merge, rebuild and run the fast
  suite before pushing.
- **API rate limits** can interrupt long agent runs. Commit often, and resume
  an interrupted agent instead of restarting it.
- **Network:** arXiv, ADS and journal sites are not reachable. Don't write
  published coefficients from memory without flagging them; ask the user for
  the paper (they supplied Glover & Jappsen 2007 for P10).
- **End of a work session:** update `docs/STATUS.md`.
