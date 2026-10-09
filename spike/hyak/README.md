# Running the framework spike on Hyak (klone)

Everything here is for the `astro` account on `cpu-g2` (AMD Genoa,
2 × 96 cores, 1.5 TB per node). It can be run by you or by the Hyak
Claude thread (`claude remote-control` in tmux on a login node).

```
cd NewRadHydroCode && git checkout framework-spike
# 1. environment: adjust module names if they differ (module avail gcc ompi cmake)
export SPIKE_MODULES="gcc ompi cmake"          # e.g. "gcc/13.2.0 ompi/4.1.6 cmake/3.28.3"
export SPIKE_WORK=/gscratch/astro/$USER/spike   # deps, builds, problems, results
source spike/hyak/env.sh
# 2. build (clones AMReX 26.10 and Parthenon v25.12 into $SPIKE_WORK/deps);
#    do it on a compute node: ~3 min with 32 cores
salloc -A astro -p cpu-g2 -c 32 --mem=64G --time=1:00:00 bash spike/hyak/build.sh
# 3. quick correctness check on the compute node (bitwise vs the serial reference)
bash spike/amrex/validate.sh $SPIKE_BUILD/serial $SPIKE_BUILD/amrex $SPIKE_WORK/val quick
# 4. benchmarks
sbatch spike/hyak/bench_node.sbatch               # one node, MPI x OpenMP layouts
sbatch --nodes=4 spike/hyak/bench_scale.sbatch    # 1, 2, 4 nodes; LAYOUT=24:8 default
```

The CSVs land in `$SPIKE_WORK/results/` (one line per run:
`kernel,code,ranks,threads,block,wall_s,work,ns_per_unit_core`).
Please copy them into `spike/bench/results/` on the branch, so the
decision document can be updated with them.

What the jobs measure:

| job | problems | layouts (ranks:threads) |
|---|---|---|
| `bench_node.sbatch` | RT log-normal 128³ (512 sources) and 256³; hydro blast 256³ | 1:12, 1:48, 8:6, 24:2, then 192 cores as 1:192, 2:96, 8:24 (1 rank per NUMA domain with NPS4), 24:8 (1 per CCD), 48:4, 96:2, 192:1; block sizes 32 and 16 |
| `bench_scale.sbatch` | strong: RT 256³; weak: RT + hydro 192³ / 240³ / 304³ on 1 / 2 / 4 nodes | `LAYOUT` per node (default 24:8; set it to the node job's winner) |

Notes:
- AMReX 26.x needs C++20, so a gcc ≥ 11 module; its GPU build needs CUDA
  ≥ 12.2. Check `module avail gcc cuda` and set `SPIKE_MODULES`.
- `env.sh` sets `OMP_PROC_BIND=close OMP_PLACES=cores`; the scripts bind each
  MPI rank to a contiguous range of cores (`--map-by slot:PE=<threads>` /
  `ppr:<n>:node:PE=<threads>`). If your MPI is not OpenMPI, set `MPIFLAGS`.
- `-march=znver4` is used by `build.sh`; it changes vector width only. The
  projects keep `-ffp-contract=off`, so results stay bitwise identical to
  the serial reference (the quick validation in step 3 checks this on the
  node).
- GPU (`ckpt` partitions, preemptible): `build_gpu.sh` builds AMReX with CUDA
  and the Parthenon prototype with Kokkos CUDA. Nothing has been run on a
  GPU yet; see the header of that script.
