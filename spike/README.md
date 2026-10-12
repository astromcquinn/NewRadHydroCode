# Framework spike (Stage C): AMReX vs Parthenon

Result and recommendation: `docs/framework_decision.md`.

| dir | what |
|---|---|
| `common/` | the host driver logic of the RT prototypes (`rt_host.hpp`) and the frozen hydro kernel (`hydro_mh.hpp`; its maintained copy is `include/radray/kernels/hydro_mh.hpp`, PR #15). The ray and RT-cell kernels (`ray_kernels.hpp`, `rt_cell.hpp`, `healpix_hd.hpp`, `consts_hd.hpp`, `rr_hd.hpp`) moved to `include/radray/kernels/` in Stage D (`docs/stage_d.md`) and are maintained there |
| `serial/` | `rt_serial`, `hydro_serial`: the reference (OpenMP, no framework); `rt_serial_fx36` round-off sensitivity build. `rt_serial` builds against the library kernels and is the bitwise reference of the AMReX backend (`ctest -L amrex`, also built as `tests/rt_serial` with `RADRAY_AMREX=ON`) |
| `amrex/` | AMReX prototype (`rt_amrex`, `hydro_amrex`, `validate.sh`, README with notes). **Frozen**: `rt_amrex.cpp` predates the moved kernels and no longer compiles; its successor is the library component `src/amrex/rt_level.cpp` (`rt_backend = amrex`) |
| `parthenon/` | Parthenon prototype (`rt_parthenon`, `hydro_parthenon`, input decks, README with notes). **Frozen** likewise |
| `bench/` | benchmark problem generator, `bench.sh` matrix runner, Sod check, `results/` |
| `hyak/` | Hyak environment, build, single-node and 1–4 node benchmark jobs, GPU build |
| `gpu/` | CUDA compile test of `common/` |
| `cmake/` | `radray_core` (the host parts of the radray library each project links) |

Only `legacy_rt = 0` is implemented: sorted ray arrays cannot reproduce the
original's zombie-ray thresholds (design doc §0.1).

Quick start (this machine):

```
cmake -S spike/serial -B build-spike-serial && cmake --build build-spike-serial -j
build-spike-serial/rt_serial tests/problems/lognormal_16src_smallbuf out_ser --set legacy_rt=0
build/rt_driver tests/problems/lognormal_16src_smallbuf out_ref --set legacy_rt=0
python3 tools/compare/compare.py out_ref out_ser   # stats identical, fields to round-off
```

Framework builds and validation: see `amrex/README.md`, `parthenon/README.md`.
