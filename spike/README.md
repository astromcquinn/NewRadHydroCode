# Framework spike (Stage C): AMReX vs Parthenon

Result and recommendation: `docs/framework_decision.md`.

| dir | what |
|---|---|
| `common/` | framework-independent, header-only host/device kernels: rays (`ray_kernels.hpp`), per-cell RT physics (`rt_cell.hpp`), MUSCL-Hancock hydro (`hydro_mh.hpp`), host driver logic (`rt_host.hpp`) |
| `serial/` | `rt_serial`, `hydro_serial`: the reference (OpenMP, no framework); `rt_serial_fx36` round-off sensitivity build |
| `amrex/` | AMReX prototype (`rt_amrex`, `hydro_amrex`, `validate.sh`, README with notes) |
| `parthenon/` | Parthenon prototype (`rt_parthenon`, `hydro_parthenon`, input decks, README with notes) |
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
