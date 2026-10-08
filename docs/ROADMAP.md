# Roadmap after the faithful port

*Written 2026-10-08 at the end of the first thread, from the user's list
of goals. `docs/PLAN.md` holds the decisions made so far, and
`docs/STATUS.md` holds the current state.*

The user's goals, in their words:

1. Adaptive mesh refinement on overdense regions for the gas and RT. For
   now, keep turning off cooling in overdense gas or quickly forming that
   gas into stars. Ramses and Enzo are the models.
2. A tree or other adaptive structure for the N-body force: TreePM as in
   GADGET/AREPO, or gravity on the AMR grid.
3. OpenMP, MPI and GPU.
4. (Wish list, possibly early.) Cut the memory and time cost of the rays,
   for example by learning the distribution of ray values and
   entropy-coding them, sacrificing flux precision where that is
   harmless, or other tricks.

The question is the order, and in particular whether MPI should come first.

---

## 1. The key decision: build the mesh ourselves or adopt a framework

Items 1–3 all depend on the same thing: the data structure that holds cells,
particles and rays across levels, MPI ranks and GPUs. Writing our own MPI for
the current uniform grid, and then rewriting it for AMR and again for GPUs,
would triple the work. So the first step is a short **framework spike**:
prototype on two or three options and decide.

| Option | What it gives us | Concerns |
|---|---|---|
| **AMReX** (C++, MPI + CUDA/HIP/SYCL + OpenMP) | Block-structured AMR, particles on AMR, multigrid Poisson, load balancing, I/O. Nyx (cosmological hydro + N-body PM on AMR) is built on it, so it is the closest existing analogue. | Large dependency, with its own idioms. Our ray-tracing data structure must be built on its particle containers or beside them. |
| **Parthenon** (C++, Kokkos + MPI) | Block-based octree AMR (MeshBlocks), Kokkos portability, swarms (particles), used by AthenaPK and KHARMA. Lighter than AMReX. | Less cosmology heritage, and no built-in multigrid gravity (it is being added). |
| **In-house** octree of fixed-size blocks (e.g. 16³), Kokkos + MPI | Full control, fits the ray-tracing algorithm exactly. | Several person-months of infrastructure that the frameworks already have. |

Cell-based octrees like Ramses suit CPUs but are awkward on GPUs.
Block-structured AMR (AMReX patches or Parthenon/Athena++ meshblocks) is the
GPU- and MPI-friendly choice. **Recommendation: pick AMReX or Parthenon after
the spike.** The test is porting the two kernels we understand best: the
hydro sweep and RT ray casting with split/merge.

Whatever we choose, `reference/` and the current serial uniform-grid code stay
as the **bit-for-bit legacy reference**. New code is validated against them on
uniform grids, with refinement off, before AMR is switched on.

## 2. Does MPI come sooner?

**For testing on the cluster: not necessarily.** A single cluster node with
OpenMP has 32–128 cores and 256 GB or more. Once the production-readiness
work (OpenMP, reduced c, ray-pool limits) lands, such a node can run 256³ RT.
That is enough for the next set of science tests and convergence studies, and
a single-node job only needs a job script.

**For production scale: yes, MPI matters.** Do it once, in the chosen
framework, right after the spike and the ray data-structure rewrite (§3). A
quick standalone MPI for the uniform grid is possible: slab decomposition fits
the original's transposes and slab FFT, plus ray migration between neighbours.
It is worth doing only if the cluster runs are needed before the framework
port is ready. It would be thrown away later.

## 3. The ray data structure, a prerequisite for everything

The original keeps rays in per-cell linked lists, as a serial or OpenMP-domain
algorithm. GPU, MPI, AMR and compression all need instead:

- rays in flat arrays (structure of arrays), sorted by (cell, cone) every substep;
- merge as a segmented reduction, split and emission as prefix-sum expansion;
- deterministic ordering, so results don't depend on the thread count, rank
  layout or backend.

`docs/design_parallel_and_encoding.md` §2–3 has the full design. This is the
first engineering step after the production run. It can be validated against
the legacy code: identical physics, with differences only from summation order,
plus a deterministic mode that is bitwise reproducible across thread counts.

## 4. Ray memory and speed (wish-list item 4): can start now

This is largely independent research, so it can run in parallel with
everything else, using ray snapshots dumped from the current code. Proposed
stages:

1. **Measure.**
   - Dump all rays at several times from production-like runs: 1 source,
     16 sources, a full box before and after overlap.
   - Report the empirical entropy of each field: direction offset from the
     HEALPix pixel centre, position within the cell, r and dr, the photon
     count in each of the 5 groups, and the spectral shape.
   - Find which fields are nearly deterministic. For example, a ray from a
     single source that has never been merged has its direction and position
     fully set by (origin, nside, ipix, r).
2. **Cheap structural compression**, about 3× at no accuracy cost. This is
   the design doc §4 encoding:
   - drop `dr`, which is derivable;
   - store the direction as a HEALPix index plus a quantised offset;
   - store position as an offset within the cell;
   - store the photon groups as a float total plus a log-quantised spectral
     shape;
   - flag never-merged rays so they need no offset at all.
3. **Learned or entropy coding.**
   - Use conditional probability models of the remaining fields, for example
     the photon count given the parent cell's optical depth and the HEALPix
     level. Code them with rANS, a fast arithmetic coder.
   - A small neural network is the general form of the model. Simple
     parametric or histogram models conditioned on a few features may get
     most of the gain at far lower cost.
   - **Key constraint:** rays are read and written every substep. Decoding
     must cost much less than the RT work per ray, about 1–3 µs on CPU,
     which favours block-wise coding of ray batches.
4. **Lossy precision with error control.**
   - Quantise photon counts with dithering (unbiased), and check that
     photon conservation and the ionisation history stay within a set
     tolerance.
   - Bit budget: fluxes to about 1%, which the user is happy to trade.
5. **Algorithmic speed-ups**, compared in `docs/design_parallel_and_encoding.md` §5:
   - merge and split criteria based on estimated error, not memory pressure;
   - an adaptive number of cones per cell (12, 48 or 192) driven by an
     I-front indicator;
   - fewer rays in fully ionised, optically thin regions, for example by
     handing radiation over to a cheap moment (M1) or background field
     after overlap;
   - casting multiple cells per substep where the gas is transparent;
   - the reduced speed of light (now in progress);
   - frequency-group experiments.

**Metrics:** bytes per ray, ns per ray per substep, and errors in x_HII(z),
T, I-front positions and the photon budget. Measure them on the Iliev et al.
2006 tests and on the production box.

## 5. AMR for gas and RT (user item 1)

- **Hydro:**
  - Port the hydro update to the framework's patches.
  - Refine on overdensity, e.g. a Lagrangian criterion with a mass per cell
    of about 8× the base level.
  - Keep the original's overcooling suppression (no cooling above 200× mean
    density), or switch to an explicit star-formation sink, as the user
    suggested. This avoids needing a full ISM model.
  - The original's directionally split TVD scheme with a moving frame does
    not carry over naturally to AMR. Expect to move to an unsplit Godunov
    scheme (MUSCL-Hancock + HLLC, as in Ramses/Athena) and validate it with
    the hydro test suite now being built.
- **RT on AMR:**
  - Rays move one local cell per substep, so fine levels need substeps
    sub-cycled within the coarse step.
  - Angular resolution (`ray_per_cell`) is set relative to the local cell
    size, and splitting continues into fine levels.
  - Merging happens within cells at the ray's current level.
  - This part is new research. Validate on the Iliev tests with refinement,
    and against the uniform-grid code at matched finest resolution.

## 6. Gravity (user item 2)

Two routes:

- **Multigrid Poisson on the AMR grid** (PM on AMR), as in Ramses, Enzo and
  Nyx. It shares the mesh with gas and RT, so it is the simplest once AMR
  exists, and it is GPU-friendly. Its force resolution follows the
  refinement.
- **TreePM** (GADGET, AREPO): better force accuracy on small scales
  independent of the gas mesh, but a separate, harder-to-GPU tree code.

**Recommendation:**
- Do multigrid on AMR first, because it comes with the framework.
- Keep TreePM as a later option if the DM halo structure on small scales
  needs it.
- Before either, test the force accuracy of the current PM (the halo mass
  function deficit seen at low resolution is a resolution effect).

## 7. Proposed order

| Stage | Work | Unlocks |
|---|---|---|
| **A (now)** | Production readiness (OpenMP, c_reduced, ray pool), hydro test suite, abundance-matched sources, 10 Mpc/h production run and movie. Cluster single-node job scripts. | First science run; cluster use without MPI. |
| **B (parallel with A/C)** | Ray measurement and compression study (§4 steps 1–3) on dumped rays. | Evidence for the encoding choice. |
| **C** | Framework spike (§1): hydro and RT kernels on AMReX vs Parthenon. Decision. | Everything below. |
| **D** | Ray data-structure rewrite (sorted arrays, deterministic), with the chosen encoding, in the framework on a uniform grid, with Kokkos/GPU. Validated against the legacy code. | GPU, MPI, AMR. |
| **E** | MPI on uniform grid in the framework: hydro, PM/FFT, RT ray migration, load balancing. Scaling tests on the cluster. | Large production runs. |
| **F** | AMR for gas (new unsplit solver), then RT on AMR. Overcooling switch or star-formation sink. | Item 1. |
| **G** | Multigrid gravity on AMR; TreePM optional. | Item 2. |
| **H** | Adaptive RT: error-driven split/merge, cones per cell, hybrid optically thin treatment. | Item 4, speed-ups. |

## 8. How to run the work: threads

The repo carries the context: `CLAUDE.md` is loaded automatically, and
`docs/PLAN.md`, `docs/STATUS.md` and this file hold the rest.

- One thread per workstream, each on its own branch and PR.
- Start each thread with: *"Read CLAUDE.md, docs/PLAN.md, docs/STATUS.md and
  docs/ROADMAP.md, then …"*.
- At the end of a thread, update `docs/STATUS.md` and record any decisions in
  `docs/PLAN.md`.

Suggested threads and opening prompts:

1. **Production run and movie (Stage A).** "… finish the follow-up PR from the
   first thread (OpenMP, c_reduced, hydro test suite) if it is not merged, then
   do abundance matching to the high-z UV luminosity functions and the
   10 Mpc/h, 0.01c production run with a movie."
2. **Cluster setup (Stage A).** Run it on the cluster with
   `claude remote-control` inside tmux. "… set up the build on this cluster:
   modules, a job-script template with restarts, a 128³ test job. Write
   docs/cluster.md."
3. **Ray compression study (Stage B).** "… implement ray snapshot dumps,
   measure per-field entropy, prototype the structural encoding and a
   learned/rANS coder, and report bytes per ray, speed and accuracy, per
   docs/ROADMAP.md §4."
4. **Framework spike (Stage C).** "… prototype the hydro sweep and RT ray
   casting on AMReX and on Parthenon, measure them, and write a decision
   document."
5. Later threads for each of Stages D–H.

Threads 1, 2 and 3 can run at the same time. Thread 4 can also start soon,
since it only produces prototypes and a decision.

## Open questions for the user

- Is there a preference or existing experience in the group with AMReX, Parthenon or Kokkos?
- The cluster's node size and GPU type (NVIDIA or AMD), and the typical queue limits.
- For AMR: overcooling suppression as in the original, or a simple
  star-formation sink? Can a sink consume dense gas with a constant
  efficiency per free-fall time?
- What precision is acceptable for ray fluxes (e.g. 1%), and which outputs
  must be preserved (x_HII(z), τ_es, I-front positions)?
