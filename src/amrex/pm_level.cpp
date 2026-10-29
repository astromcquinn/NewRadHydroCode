// The AMReX backend of the PM gravity (pm_backend = amrex; docs/stage_e_pm.md):
// one PM step (link, density, Poisson solve, kick) on one uniform AMReX
// level of the PM mesh, MPI + OpenMP, with the particle and cell arithmetic
// of the serial port (include/radray/kernels/pm_kernels.hpp). Stage E2 of
// docs/framework_decision.md section 6.
//
//   * particles  ParticleContainerPureSoA<9, 0>: AMReX position = the centre
//                of the particle's home mesh cell (link_dm's cell), id = the
//                host index ip, SoA reals = the position in DM-grid units
//                and the velocity (the original's float xv_dm);
//   * cells      FabArray<BaseFab<float>> pot (density, then the potential;
//                3 ghost layers for the 4-point gradient), frc (the mesh
//                force, 3 components, 1 ghost layer: the CIC reach), gas
//                (u(5) + v(3) or e on the gas grid, the mesh BoxArray
//                coarsened by pm_refine, no ghosts);
//   * deposits   FabArray<BaseFab<Long>> acc: the CIC mass in int64 fixed
//                point (2^-kDensBits), per-thread buffers of the grown box,
//                then SumBoundary (exact, associative);
//   * sums       the kick statistics arms and vrms in 128-bit fixed point,
//                amax and vmax as maxima; all reductions are integers or
//                maxima, so the result is bitwise identical for any rank
//                count, box size and thread count;
//   * Poisson    amrex_pm_fft = pencil (default): 1-D FFTW transforms per
//                pencil on z- and y-slab layouts moved with ParallelCopy
//                (every line gets the same plan whatever the layout: rank-
//                and thread-independent); amrex: amrex::FFT::R2C.
//
// Stage E2's MPI model: the hydro and the RT hold the full state on every
// rank, so the level copies the host state in at the start of the step
// (every rank creates the particles of its own boxes from s.xv and fills
// its gas boxes from u_gas/v_gas) and gathers the results back onto every
// rank at the end (the kicked velocities into s.xv, the kicked gas into
// u_gas/v_gas: each rank writes its own cells into a zeroed array, then an
// MPI bitwise-OR all-reduce; exact). The level holds no state between
// steps, so checkpoints and restarts are the host's.
//
// Stage G (pm_level_solve = 1, docs/stage_g.md): with a refined gas
// hierarchy the particles live on the finest level covering them (one
// ParticleContainer with the hierarchy's levels, Redistribute by level),
// every level deposits its own particles and gas on its own mesh (the
// level's gas boxes refined by pm_refine, dx = 2^-l level-0 mesh cells),
// the deposits are restricted onto the parents (so level 0 sees all the
// mass), level 0 keeps the FFT and each refined level solves
// lap phi = (2a/3)(rho - mean) with AMReX's MLMG (MLPoisson, Dirichlet
// values from the parent's potential, a fixed number of V-cycles, the
// smoother as the bottom solver: every operation is per cell, so the solve
// is bitwise rank-independent), then the 4-point gradient and the kicks on
// each level. See FineMesh below.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <AMReX.H>
#include <AMReX_FFT.H>
#include <AMReX_FabArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuComplex.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLPoisson.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParticleLocator.H>
#include <AMReX_Particles.H>
#include <AMReX_iMultiFab.H>
#include <fftw3.h>
#ifdef AMREX_USE_OMP
#include <omp.h>
#endif

#include "exact_sum.hpp"
#include "gas_state.hpp"
#include "radray/constants.hpp"
#include "radray/fortran.hpp"
#include "radray/kernels/hydro_mh.hpp"
#include "radray/kernels/pm_kernels.hpp"
#include "radray/pm.hpp"
#include "radray/pm_backend.hpp"
#include "radray/rt_backend.hpp"

namespace radray {

namespace {

namespace ax = amrex;
using ax::Long;
using FFab = ax::FabArray<ax::BaseFab<float>>;
using LFab = ax::FabArray<ax::BaseFab<Long>>;
using CFab = ax::FabArray<ax::BaseFab<ax::GpuComplex<float>>>;
// 3 positions (the home cell) + x(3) + v(3) + the star's mass and formation
// scale factor (0 for dark matter); one int: the type (0 dark matter, 1 a
// star of the sink, docs/stage_f.md section 2)
using PC = ax::ParticleContainerPureSoA<11, 1>;
using PTD = PC::ParticleTileType::ParticleTileDataType;
constexpr int kX = 3, kV = 6, kM = 9, kA = 10;  // SoA real components
constexpr int kType = 0;                       // SoA int component

static_assert(std::is_same_v<ax::ParticleReal, float>, "build AMReX with AMReX_PARTICLES_PRECISION=SINGLE");
static_assert(sizeof(Long) == 8, "int64 accumulators");
static_assert(sizeof(ax::GpuComplex<float>) == sizeof(fftwf_complex), "complex layout");

// The CIC mass in int64 fixed point: resolution 2^-38 = 3.6e-12 (below the
// float resolution of any density above 2^-14), range 2^25 = 3.4e7 per
// cell (the densest cell of a 1024^3 run holds ~1e5 particles of mass
// ~0.8 in mesh units).
constexpr int kDensBits = 38;
// The kick statistics (diagnostics: arms = sum |f|^2 over the particles,
// vrms = sum |v|^2) in 128-bit fixed point: resolution 2^-56 = 1.4e-17,
// range 2^71.
constexpr int kStatBits = 56;
// (exact_sum.hpp)
using xs::allreduce_i128;
using xs::from_fx128;
using xs::i128;
using xs::to_fx128;

inline Long to_fx(double x, int bits) { return static_cast<Long>(std::floor(std::ldexp(x, bits) + 0.5)); }

// In-place bitwise OR of an array over all ranks (exact gather of disjoint
// pieces: every rank holds its own values and zeros elsewhere).
void allreduce_bor(float* data, std::size_t n) {
  if (ax::ParallelDescriptor::NProcs() == 1) return;
  auto* u = reinterpret_cast<std::uint32_t*>(data);
  const std::size_t chunk = std::size_t(1) << 28;
  for (std::size_t off = 0; off < n; off += chunk) {
    const int cnt = static_cast<int>(std::min(chunk, n - off));
    MPI_Allreduce(MPI_IN_PLACE, u + off, cnt, MPI_UINT32_T, MPI_BOR, ax::ParallelDescriptor::Communicator());
  }
}

// Wrapped global 0-based index g -> its periodic image nearest to the home
// index h
inline int near_image(int g, int h, int N) {
  int d = g - h;
  if (2 * d > N) d -= N;
  else if (2 * d < -N) d += N;
  return h + d;
}

double now() {
  return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

#ifdef AMREX_USE_OMP
inline int nthreads_max() { return omp_get_max_threads(); }
inline int thread_id() { return omp_get_thread_num(); }
#else
inline int nthreads_max() { return 1; }
inline int thread_id() { return 0; }
#endif

struct Timers {
  double prep = 0, drift = 0, deposit = 0, sumbnd = 0, fft = 0, force = 0, kick_dm = 0, kick_gas = 0, gather = 0;
  double d_fine = 0, k_fine = 0;  // pm_gas_levels: the refined levels' deposit and kick (+ average-down)
  double sync = 0, d_lev = 0, solve = 0, f_lev = 0, k_lev = 0;  // pm_level_solve: the refined levels
  double g_vwrite = 0, g_vred = 0, g_vcopy = 0, g_gas = 0, g_pot = 0, p_locate = 0, p_fill = 0;
  Long steps = 0;
};

// ---------------------------------------------------------------------------
// The pencil FFT: real (N^3) <-> half-complex ((N/2+1) x N x N) with 1-D
// FFTW transforms, one plan per direction, executed line by line
// (FFTW_ESTIMATE | FFTW_UNALIGNED, howmany = 1, single-threaded plans), so
// every line gets exactly the same arithmetic whatever the rank and thread
// layout. The x and y transforms run on z-slabs (full x and y per box), the
// z transform on y-slabs; the two layouts are connected by ParallelCopy.
// ---------------------------------------------------------------------------
struct PencilFFT {
  int N = 0, H = 0;
  ax::BoxArray baz_r, baz_c, bay_c;  // z-slabs (real, complex), y-slabs (complex)
  ax::DistributionMapping dmz, dmy;
  FFab rz;
  CFab cz, cy;
  fftwf_plan p_r2c = nullptr, p_c2r = nullptr, p_fwd = nullptr, p_bwd = nullptr;

  explicit PencilFFT(int n) : N(n), H(n / 2) {
    const int nranks = ax::ParallelDescriptor::NProcs();
    const ax::Box rdom(ax::IntVect(0), ax::IntVect(N - 1));
    const ax::Box cdom(ax::IntVect(0), ax::IntVect(H, N - 1, N - 1));
    const int slab = std::max(1, (N + nranks - 1) / nranks);
    baz_r = ax::BoxArray(rdom);
    baz_r.maxSize(ax::IntVect(N, N, slab));
    baz_c = ax::BoxArray(cdom);
    baz_c.maxSize(ax::IntVect(H + 1, N, slab));
    bay_c = ax::BoxArray(cdom);
    bay_c.maxSize(ax::IntVect(H + 1, slab, N));
    dmz = ax::DistributionMapping(baz_r);
    dmy = ax::DistributionMapping(bay_c);
    rz.define(baz_r, dmz, 1, 0);
    cz.define(baz_c, dmz, 1, 0);
    cy.define(bay_c, dmy, 1, 0);
    // the plans (the planner is not thread-safe: made once, here)
    std::vector<float> rb(static_cast<std::size_t>(N) + 2);
    std::vector<fftwf_complex> cb(static_cast<std::size_t>(N));
    const unsigned flags = FFTW_ESTIMATE | FFTW_UNALIGNED;
    p_r2c = fftwf_plan_dft_r2c_1d(N, rb.data(), cb.data(), flags);
    p_c2r = fftwf_plan_dft_c2r_1d(N, cb.data(), rb.data(), flags);
    p_fwd = fftwf_plan_dft_1d(N, cb.data(), cb.data(), FFTW_FORWARD, flags);
    p_bwd = fftwf_plan_dft_1d(N, cb.data(), cb.data(), FFTW_BACKWARD, flags);
    if (!p_r2c || !p_c2r || !p_fwd || !p_bwd) throw std::runtime_error("pm_backend = amrex: FFTW planning failed");
  }
  ~PencilFFT() {
    fftwf_destroy_plan(p_r2c);
    fftwf_destroy_plan(p_c2r);
    fftwf_destroy_plan(p_fwd);
    fftwf_destroy_plan(p_bwd);
  }

  // c2c along y (dir = 1) or z (dir = 2) of every line of a FabArray whose
  // boxes span the full extent in that direction
  void c2c_lines(CFab& c, int dir, fftwf_plan plan) {
    for (ax::MFIter mfi(c); mfi.isValid(); ++mfi) {
      const ax::Box bx = mfi.validbox();
      auto a = c.array(mfi);
      const auto lo = ax::lbound(bx), hi = ax::ubound(bx);
      // the lines: (i, k) for y, (i, j) for z
      const int n1 = hi.x - lo.x + 1;
      const int n2 = dir == 1 ? hi.z - lo.z + 1 : hi.y - lo.y + 1;
      const Long nl = static_cast<Long>(n1) * n2;
#pragma omp parallel
      {
        std::vector<ax::GpuComplex<float>> buf(static_cast<std::size_t>(N));
        auto* fb = reinterpret_cast<fftwf_complex*>(buf.data());
#pragma omp for schedule(static)
        for (Long l = 0; l < nl; ++l) {
          const int i = lo.x + static_cast<int>(l % n1);
          const int o = static_cast<int>(l / n1);
          if (dir == 1) {
            const int k = lo.z + o;
            for (int j = 0; j < N; ++j) buf[j] = a(i, lo.y + j, k);
            fftwf_execute_dft(plan, fb, fb);
            for (int j = 0; j < N; ++j) a(i, lo.y + j, k) = buf[j];
          } else {
            const int j = lo.y + o;
            for (int k = 0; k < N; ++k) buf[k] = a(i, j, lo.z + k);
            fftwf_execute_dft(plan, fb, fb);
            for (int k = 0; k < N; ++k) a(i, j, lo.z + k) = buf[k];
          }
        }
      }
    }
  }

  // in: the valid region of a real FabArray on any layout; post(i, j, k,
  // c) on every spectral coefficient (global indices); the backward
  // transform (unnormalised, as MKL/FFTW) into out's valid region
  template <class F>
  void forward_then_backward(const FFab& in, FFab& out, F&& post) {
    rz.ParallelCopy(in, 0, 0, 1);
    // x: r2c line by line (contiguous in i), into cz
    for (ax::MFIter mfi(rz); mfi.isValid(); ++mfi) {
      const ax::Box bx = mfi.validbox();
      auto r = rz.array(mfi);
      auto c = cz.array(mfi);
      const auto lo = ax::lbound(bx), hi = ax::ubound(bx);
      const int ny = hi.y - lo.y + 1, nz = hi.z - lo.z + 1;
      const Long nl = static_cast<Long>(ny) * nz;
#pragma omp parallel
      {
        std::vector<float> rb(static_cast<std::size_t>(N));
        std::vector<ax::GpuComplex<float>> cb(static_cast<std::size_t>(H) + 1);
#pragma omp for schedule(static)
        for (Long l = 0; l < nl; ++l) {
          const int j = lo.y + static_cast<int>(l % ny), k = lo.z + static_cast<int>(l / ny);
          for (int i = 0; i < N; ++i) rb[i] = r(i, j, k);
          fftwf_execute_dft_r2c(p_r2c, rb.data(), reinterpret_cast<fftwf_complex*>(cb.data()));
          for (int i = 0; i <= H; ++i) c(i, j, k) = cb[i];
        }
      }
    }
    c2c_lines(cz, 1, p_fwd);     // y
    cy.ParallelCopy(cz, 0, 0, 1);
    c2c_lines(cy, 2, p_fwd);     // z
    for (ax::MFIter mfi(cy); mfi.isValid(); ++mfi) {
      const ax::Box bx = mfi.validbox();
      auto a = cy.array(mfi);
      ax::ParallelForOMP(bx, [=](int i, int j, int k) { post(i, j, k, a(i, j, k)); });
    }
    c2c_lines(cy, 2, p_bwd);
    cz.ParallelCopy(cy, 0, 0, 1);
    c2c_lines(cz, 1, p_bwd);
    for (ax::MFIter mfi(rz); mfi.isValid(); ++mfi) {
      const ax::Box bx = mfi.validbox();
      auto r = rz.array(mfi);
      auto c = cz.array(mfi);
      const auto lo = ax::lbound(bx), hi = ax::ubound(bx);
      const int ny = hi.y - lo.y + 1, nz = hi.z - lo.z + 1;
      const Long nl = static_cast<Long>(ny) * nz;
#pragma omp parallel
      {
        std::vector<float> rb(static_cast<std::size_t>(N));
        std::vector<ax::GpuComplex<float>> cb(static_cast<std::size_t>(H) + 1);
#pragma omp for schedule(static)
        for (Long l = 0; l < nl; ++l) {
          const int j = lo.y + static_cast<int>(l % ny), k = lo.z + static_cast<int>(l / ny);
          for (int i = 0; i <= H; ++i) cb[i] = c(i, j, k);
          fftwf_execute_dft_c2r(p_c2r, reinterpret_cast<fftwf_complex*>(cb.data()), rb.data());
          for (int i = 0; i < N; ++i) r(i, j, k) = rb[i];
        }
      }
    }
    out.ParallelCopy(rz, 0, 0, 1);
  }
};

// ---------------------------------------------------------------------------
// pm_level_solve (Stage G): a refined level's mesh
// ---------------------------------------------------------------------------
//
// Level l (1..max_level) of the gas hierarchy (gs::GasState::fine[l - 1])
// gets a PM mesh of its gas boxes refined by pm_refine: Nm 2^l cells per
// dimension in the periodic domain, dx = 2^-l level-0 mesh cells. Its
// particles are those of the hierarchy's ParticleContainer at level l (the
// finest level whose boxes hold the particle's home cell). Deposits are in
// int64 fixed point as on level 0, in density units (mass per level-0
// mesh-cell volume: a particle of mass m in a level-l cell is a density
// m 8^l), with `bits` fractional bits: 38 to level 2 and 3 fewer per level
// below (26 at most), so that the range (2^(63-bits) per cell) holds the
// deep levels' densities; the parent's value is the sum of the 8 children
// (exact when the bits differ by 3). own = the CIC deposit of the level's
// own particles (after SumBoundary), tot = own + the restricted deposits of
// the finer levels + the finer levels' ghost deposits that fell outside
// them (their mass lumped into the adjacent parent cell). The density the
// level solves is gas + tot + own(parent) prolonged (the parent-level
// particles whose CIC cloud reaches into the level).
struct FineMesh {
  int lev = 0;
  int bits = kDensBits;
  ax::Box domain;
  ax::Geometry geom;
  ax::Periodicity period;
  ax::BoxArray ba, bag;  // the mesh (the level's gas boxes refined by pm_refine), the gas boxes
  ax::DistributionMapping dm;
  FFab pot, frc, gas;
  LFab acc, own, tot;
  ax::iMultiFab cov;  // 1 where a mesh cell (valid or one of the 3 ghost layers) lies in the level
  ax::ParticleLocator<ax::DenseBins<ax::Box>> locator;
  bool aggable = false;  // the boxes fill their bounding box (MLMG may agglomerate: deeper multigrid)
  int mg_depth = 0;      // the multigrid coarsenings (layout-independent)
  int bottom = 0;        // the coarsest multigrid grid's widest side (aggable levels; 0 = not known)
  int sweeps = 0;        // the bottom solve's red-black sweeps (solve_level)
  int iters = 0;         // the last solve: V-cycles, final / initial residual
  double resid = 0;
  Long capped = 0;       // solves that ended at the V-cycle cap above kCapResid
};

// A level solve that ends at pm_mlmg_iters with a relative residual above
// this is reported (the float potential's rounding is ~1e-7)
constexpr double kCapResid = 1e-6;

// The fractional bits of level l's fixed-point deposits (density units)
inline int bits_of(int l) { return std::max(26, kDensBits - 3 * std::max(0, l - 2)); }

// Interpolation of a parent-level array at the centre of fine cell (i, j,
// k), which sits at (i + 1/2)/2 - 1/2 in parent index space: 1/4 past the
// parent cell's centre for an odd i, 1/4 before it for an even i.
// Trilinear (the initial guess of the level solve; 2 parent ghosts):
AMREX_FORCE_INLINE double interp_parent(const ax::Array4<const float>& q, int i, int j, int k) {
  const int i0 = (i >> 1) - (1 - (i & 1)), j0 = (j >> 1) - (1 - (j & 1)), k0 = (k >> 1) - (1 - (k & 1));
  const double wx = (i & 1) ? 0.25 : 0.75, wy = (j & 1) ? 0.25 : 0.75, wz = (k & 1) ? 0.25 : 0.75;
  double v = 0;
  for (int dk = 0; dk < 2; ++dk)
    for (int dj = 0; dj < 2; ++dj)
      for (int di = 0; di < 2; ++di) {
        const double w = (di ? wx : 1 - wx) * (dj ? wy : 1 - wy) * (dk ? wz : 1 - wz);
        v += w * static_cast<double>(q(i0 + di, j0 + dj, k0 + dk));
      }
  return v;
}
// Tri-quadratic (the potential's ghost cells outside a level, which the
// 4-point gradient of the boundary cells and the force of the ghost layer
// read: the trilinear error, h^2 phi''/8, is a few per cent of the force
// 8 cells from a mass; the quadratic one is O(h^3 phi''')); the 3 parent
// cells nearest in each direction: 3 parent ghosts around the coarsened box
AMREX_FORCE_INLINE double interp_parent_quad(const ax::Array4<const float>& q, int i, int j, int k) {
  const int c[3] = {i >> 1, j >> 1, k >> 1};  // the parent cell holding the fine cell
  double w[3][3];
  const int idx[3] = {i, j, k};
  for (int d = 0; d < 3; ++d) {
    const double t = (idx[d] & 1) ? 0.25 : -0.25;  // the fine centre relative to the parent centre
    w[d][0] = t * (t - 1) / 2;
    w[d][1] = 1 - t * t;
    w[d][2] = t * (t + 1) / 2;
  }
  double v = 0;
  for (int dk = 0; dk < 3; ++dk)
    for (int dj = 0; dj < 3; ++dj)
      for (int di = 0; di < 3; ++di)
        v += w[0][di] * w[1][dj] * w[2][dk] * static_cast<double>(q(c[0] + di - 1, c[1] + dj - 1, c[2] + dk - 1));
  return v;
}

// ---------------------------------------------------------------------------
// The level
// ---------------------------------------------------------------------------
struct PMLevel {
  int Nm = 0, Ng = 0, r = 1, mgs = 16;
  std::int64_t Npart = 0;
  bool mf = true;
  int fft_kind = 0;
  rpk::MeshGeom mg;
  double x_dm_m = 1;
  double x_dm_gas = 1;  // Ngas/Ndm (DM to gas grid velocity units)
  ax::Box domain;
  ax::Geometry geom;
  ax::Periodicity period;
  ax::BoxArray ba, bag;
  ax::DistributionMapping dm;
  FFab pot, frc, gas;
  LFab acc;
  std::unique_ptr<PC> pc;
  ax::ParticleLocator<ax::DenseBins<ax::Box>> locator;
  std::unique_ptr<PencilFFT> pfft;
  std::unique_ptr<ax::FFT::R2C<float>> afft;
  gs::GasState* state = nullptr;  // amrex_state = 1: the gas on the level (gas_state.hpp)
  bool loaded = false;            // amrex_state = 1: the particles live on the level
  bool gas_levels = false;        // pm_gas_levels: the state's refined levels enter the PM
  int lr = 0;                     // log2(pm_refine) (a power of 2 with gas_levels)
  // pm_level_solve (Stage G): the hierarchy's depth (the particle container
  // has max_level + 1 levels), the refined meshes in use (fm[l - 1] = level
  // l), level 0's own and total deposits, the mean density (the k = 0 mode
  // the FFT drops; subtracted on the refined levels), the solver settings
  bool level_solve = false;
  int max_level = 0;
  int bf = 8;  // the hierarchy's blocking factor (the multigrid depth)
  int mlmg_iters = 20, mlmg_bottom_sweeps = 8, mlmg_verbose = 0;
  double mlmg_tol = 1e-10, mlmg_bottom_scale = 1.0;
  std::vector<FineMesh> fm;
  LFab own0, tot0;
  double rho_mean = 0;
  Timers T;

  explicit PMLevel(const pm::PMConfig& pc_) {
    Nm = pc_.Nm;
    Ng = pc_.Ngas;
    r = pc_.refine;
    mgs = pc_.amrex_max_grid_size;
    Npart = pc_.Npart;
    mf = pc_.mf;
    fft_kind = pc_.amrex_fft;
    mg = rpk::MeshGeom{pc_.Nm, pc_.nmr, pc_.x_dm_m};
    x_dm_m = pc_.x_dm_m;
    x_dm_gas = pc_.x_dm_gas;
    if (mgs < 4 || mgs % r != 0 || Nm % mgs != 0)
      throw std::runtime_error("pm_backend = amrex: the box size must be a multiple of pm_refine and divide the mesh");
    domain = ax::Box(ax::IntVect(0), ax::IntVect(Nm - 1));
    const ax::RealBox rb({0., 0., 0.}, {double(Nm), double(Nm), double(Nm)});  // dx = 1
    const ax::Array<int, 3> per{1, 1, 1};
    geom = ax::Geometry(domain, rb, ax::CoordSys::cartesian, per);
    period = geom.periodicity();
    ba = ax::BoxArray(domain);
    ba.maxSize(mgs);
    dm = ax::DistributionMapping(ba);
    bag = ba;
    bag.coarsen(r);  // the gas grid: each gas cell's r^3 mesh cells lie in its box
    if (pc_.amrex_state) {
      // the gas state's decomposition, refined by pm_refine for the mesh
      // (gas_state.hpp): the gas is read and written locally
      state = gs::active();
      if (!state) throw std::runtime_error("amrex_state = 1: amrex_state_init was not called");
      if (state->N != Ng) throw std::runtime_error("pm_backend = amrex: the gas state has another N");
      bag = state->ba;
      ba = state->ba;
      ba.refine(r);
      dm = state->dm;
      gas_levels = pc_.gas_levels;
      if (gas_levels && (r & (r - 1)) != 0)
        throw std::runtime_error("pm_gas_levels = 1 needs pm_refine a power of 2");
      while ((1 << lr) < r) ++lr;
      level_solve = pc_.level_solve;
      max_level = level_solve ? pc_.max_level : 0;
      bf = pc_.blocking_factor;
      mlmg_iters = pc_.mlmg_iters;
      mlmg_tol = pc_.mlmg_tol;
      mlmg_bottom_sweeps = pc_.mlmg_bottom_sweeps;
      mlmg_bottom_scale = pc_.mlmg_bottom_scale;
      mlmg_verbose = pc_.mlmg_verbose;
      if (level_solve && static_cast<std::int64_t>(Nm) * (std::int64_t(1) << max_level) > (std::int64_t(1) << 23))
        throw std::runtime_error("pm_level_solve = 1: the finest mesh exceeds 2^23 cells per dimension (the "
                                 "particle positions are floats)");
    }
    pot.define(ba, dm, 1, 3);
    frc.define(ba, dm, 3, 1);
    acc.define(ba, dm, 1, 1);
    gas.define(bag, dm, 8, 0);
    if (level_solve) {
      // the hierarchy's levels in one container: level l's domain is the
      // mesh refined 2^l times in the same box; the refined levels get their
      // boxes from the gas state at every step (sync_levels)
      own0.define(ba, dm, 1, 0);
      tot0.define(ba, dm, 1, 0);
      ax::Vector<ax::Geometry> geoms(static_cast<std::size_t>(max_level) + 1);
      ax::Vector<ax::DistributionMapping> dms(static_cast<std::size_t>(max_level) + 1);
      ax::Vector<ax::BoxArray> bas(static_cast<std::size_t>(max_level) + 1);
      ax::Vector<int> rr(static_cast<std::size_t>(max_level), 2);
      for (int l = 0; l <= max_level; ++l)
        geoms[static_cast<std::size_t>(l)] =
            ax::Geometry(ax::Box(ax::IntVect(0), ax::IntVect((Nm << l) - 1)), rb, ax::CoordSys::cartesian, per);
      bas[0] = ba;
      dms[0] = dm;
      pc = std::make_unique<PC>(geoms, dms, bas, rr);
    } else {
      pc = std::make_unique<PC>(geom, dm, ba);
    }
    for (ax::MFIter mfi(pot); mfi.isValid(); ++mfi)
      pc->DefineAndReturnParticleTile(0, mfi.index(), mfi.LocalTileIndex());
    locator.build(ba, geom);
  }

  // ---- pm_level_solve: the levels ---------------------------------------------

  int finest() const { return static_cast<int>(fm.size()); }  // the refined levels in use (0 without level_solve)
  FFab& POT(int l) { return l == 0 ? pot : fm[static_cast<std::size_t>(l) - 1].pot; }
  FFab& FRC(int l) { return l == 0 ? frc : fm[static_cast<std::size_t>(l) - 1].frc; }
  FFab& GAS(int l) { return l == 0 ? gas : fm[static_cast<std::size_t>(l) - 1].gas; }
  LFab& ACC(int l) { return l == 0 ? acc : fm[static_cast<std::size_t>(l) - 1].acc; }
  LFab& OWN(int l) { return l == 0 ? own0 : fm[static_cast<std::size_t>(l) - 1].own; }
  LFab& TOT(int l) { return l == 0 ? tot0 : fm[static_cast<std::size_t>(l) - 1].tot; }
  const ax::Periodicity& PER(int l) const { return l == 0 ? period : fm[static_cast<std::size_t>(l) - 1].period; }
  const ax::BoxArray& BA(int l) const { return l == 0 ? ba : fm[static_cast<std::size_t>(l) - 1].ba; }
  const ax::DistributionMapping& DM(int l) const { return l == 0 ? dm : fm[static_cast<std::size_t>(l) - 1].dm; }
  int BITS(int l) const { return l == 0 ? kDensBits : fm[static_cast<std::size_t>(l) - 1].bits; }
  const ax::ParticleLocator<ax::DenseBins<ax::Box>>& LOC(int l) const {
    return l == 0 ? locator : fm[static_cast<std::size_t>(l) - 1].locator;
  }
  // The CIC geometry of level l's mesh (Nm 2^l cells, x_dm_m 2^l mesh cells
  // per DM-grid cell: exact scalings of level 0's)
  rpk::MeshGeom mesh_geom(int l) const {
    return rpk::MeshGeom{Nm << l, static_cast<double>(Nm << l), x_dm_m * static_cast<double>(1 << l)};
  }
  // The AMReX position of a particle at x (DM-grid units): the centre of
  // its home cell at the finest level's resolution (max_level = 0: level
  // 0's, as before), so that floor(pos 2^l) is the home cell at every
  // level l; exact in float while Nm 2^max_level <= 2^23
  AMREX_FORCE_INLINE float fpos(float x) const {
    if (max_level == 0) return static_cast<float>(home(x)) + 0.5f;
    const int L = max_level;
    const int c = f::mod(f::int4(static_cast<double>(x) * x_dm_m * static_cast<double>(1 << L)), Nm << L);
    return (static_cast<float>(c) + 0.5f) / static_cast<float>(1 << L);
  }
  // The home cell at level l of an AMReX position
  static AMREX_FORCE_INLINE int home_l(float pos, int l) {
    return static_cast<int>(std::floor(pos * static_cast<float>(1 << l)));
  }
  // Every particle tile of every level: f(lev, mfi) with mfi over the
  // level's mesh FabArray
  template <class F>
  void for_tiles(F&& f) {
    for (int l = 0; l <= finest(); ++l)
      for (ax::MFIter mfi(POT(l)); mfi.isValid(); ++mfi) f(l, mfi);
  }
  PC::ParticleTileType& tile(int l, const ax::MFIter& mfi) {
    return pc->DefineAndReturnParticleTile(l, mfi.index(), mfi.LocalTileIndex());
  }

  // The refined meshes from the gas state's levels (after the hydro's
  // regrid): redefined when a level's boxes or mapping changed, and the
  // particle container's levels set accordingly (the caller's Redistribute
  // then moves the particles onto them). Collective.
  bool sync_levels() {
    if (!level_solve) return false;
    const double t0 = now();
    const int nf = state->finest;
    if (nf > max_level) throw std::runtime_error("pm_level_solve: the gas state has more levels than amrex_max_level");
    bool changed = finest() != nf;
    for (int l = 1; l <= nf && !changed; ++l) {
      const gs::GasState::FineLevel& F = state->fine[static_cast<std::size_t>(l) - 1];
      const FineMesh& M = fm[static_cast<std::size_t>(l) - 1];
      if (M.lev != l || !(M.bag == F.ba) || !(M.dm == F.dm)) changed = true;
    }
    if (!changed) {
      T.sync += now() - t0;
      return false;
    }
    const ax::RealBox rb({0., 0., 0.}, {double(Nm), double(Nm), double(Nm)});
    const ax::Array<int, 3> per{1, 1, 1};
    fm.resize(static_cast<std::size_t>(nf));
    for (int l = 1; l <= nf; ++l) {
      FineMesh& M = fm[static_cast<std::size_t>(l) - 1];
      const gs::GasState::FineLevel& F = state->fine[static_cast<std::size_t>(l) - 1];
      if (M.lev == l && M.bag == F.ba && M.dm == F.dm) continue;
      M.lev = l;
      M.bits = bits_of(l);
      M.domain = ax::Box(ax::IntVect(0), ax::IntVect((Nm << l) - 1));
      M.geom = ax::Geometry(M.domain, rb, ax::CoordSys::cartesian, per);
      M.period = M.geom.periodicity();
      M.bag = F.ba;
      M.ba = F.ba;
      M.ba.refine(r);
      M.dm = F.dm;
      M.pot.define(M.ba, M.dm, 1, 3);
      M.frc.define(M.ba, M.dm, 3, 1);
      M.gas.define(M.bag, M.dm, 8, 0);
      M.acc.define(M.ba, M.dm, 1, 1);
      M.own.define(M.ba, M.dm, 1, 0);
      M.tot.define(M.ba, M.dm, 1, 0);
      M.cov.define(M.ba, M.dm, 1, 3);
      M.cov.setVal(0);
      M.cov.setVal(1, 0, 1, 0);
      M.cov.FillBoundary(M.period);
      M.locator.build(M.ba, M.geom);
      M.aggable = M.ba.minimalBox().numPts() == M.ba.numPts();
      // the boxes are multiples of the blocking factor (bf gas cells = bf r
      // mesh cells) however the rank count chops them, so log2(bf r) - 1
      // coarsenings (to a width of 2) are possible on any layout; a level
      // whose boxes fill their bounding box can be agglomerated and coarsened
      // further by MLMG, again independently of the layout
      int depth = -1;
      for (int w = bf * r; w > 1; w >>= 1) ++depth;
      M.mg_depth = M.aggable ? 30 : std::max(0, depth);
      // An aggable level's multigrid halves its bounding box while every side
      // stays even (MLLinOp's agglomeration: down to a width of 2), so the
      // bottom grid's side is the odd part of each width: 9 for 144 = 16 x 9,
      // 65 for a padded 1040. The smoother bottom solve needs of order
      // width^2 sweeps there (solve_level). Layout-independent.
      M.bottom = 0;
      if (M.aggable) {
        const ax::Box mb = M.ba.minimalBox();
        for (int d = 0; d < 3; ++d) {
          int w = mb.length(d);
          while (w % 2 == 0 && w / 2 >= 2) w /= 2;
          M.bottom = std::max(M.bottom, w);
        }
      }
      // proper nesting: the parent's cells must cover 3 parent cells around
      // the level (the Dirichlet data and the ghost potential come from them;
      // AmrMesh nests by blocking_factor / 2 = 4 cells)
      {
        ax::iMultiFab m2(ax::coarsen(M.ba, 2), M.dm, 1, 3), ones(BA(l - 1), DM(l - 1), 1, 0);
        m2.setVal(0);
        ones.setVal(1);
        m2.ParallelCopy(ones, 0, 0, 1, 0, 3, PER(l - 1));
        if (m2.min(0, 3) < 1)
          throw std::runtime_error("pm_level_solve: level " + std::to_string(l) +
                                   " is not nested 3 parent cells inside level " + std::to_string(l - 1));
      }
      pc->SetParticleBoxArray(l, M.ba);
      pc->SetParticleDistributionMap(l, M.dm);
    }
    for (int l = nf + 1; l <= max_level; ++l) {
      pc->SetParticleBoxArray(l, ax::BoxArray());
      pc->SetParticleDistributionMap(l, ax::DistributionMapping());
    }
    T.sync += now() - t0;
    return true;
  }

  // ---- host state -> level ------------------------------------------------

  // The particles of s.xv whose home cell (link_dm's cell) lies in a local
  // box, and the gas of the local gas boxes (e_gas: the dual-energy field in
  // place of v_gas when the kick is the unsplit one)
  // ---- amrex_balance = 1: the particles onto a new mapping ------------------

  // The persistent particles move by a Redistribute under the new mapping;
  // the mesh and gas FabArrays are work arrays, rebuilt
  void remap(const ax::DistributionMapping& dm2) {
    const double t0 = now();
    dm = dm2;
    pot.define(ba, dm, 1, 3);
    frc.define(ba, dm, 3, 1);
    acc.define(ba, dm, 1, 1);
    gas.define(bag, dm, 8, 0);
    if (level_solve) {
      own0.define(ba, dm, 1, 0);
      tot0.define(ba, dm, 1, 0);
    }
    pc->SetParticleDistributionMap(0, dm);
    pc->Redistribute();
    for (ax::MFIter mfi(pot); mfi.isValid(); ++mfi)
      pc->DefineAndReturnParticleTile(0, mfi.index(), mfi.LocalTileIndex());
    pfft.reset();
    afft.reset();
    T.prep += now() - t0;
  }

  // ---- amrex_state = 1: the particles on the level between steps -----------

  // link_dm's home cell of a position (mesh units): int(x * x_dm_m) mod Nm
  AMREX_FORCE_INLINE int home(float x) const { return f::mod(f::int4(static_cast<double>(x) * x_dm_m), Nm); }

  // Rank 0's xv onto the level: broadcast in chunks, every rank keeps the
  // particles of its own boxes (id = the host index ip). Collective.
  void load_particles(const pm::PMState& s) {
    const double t0 = now();
    const int myproc = ax::ParallelDescriptor::MyProc();
    const int root = ax::ParallelDescriptor::IOProcessorNumber();
    const bool io = ax::ParallelDescriptor::IOProcessor();
    if (io && static_cast<std::int64_t>(s.xv.size()) != 6 * Npart)
      throw std::runtime_error("amrex_state = 1: rank 0's xv is not allocated");
    for (ax::MFIter mfi(pot); mfi.isValid(); ++mfi)
      pc->DefineAndReturnParticleTile(0, mfi.index(), mfi.LocalTileIndex()).resize(0);
    const auto assign = locator.getGridAssignor();
    constexpr std::int64_t kChunk = std::int64_t(1) << 22;  // 4 M particles, 100 MB
    std::vector<float> buf;
    Long nmine = 0;
    for (std::int64_t ip0 = 1; ip0 <= Npart; ip0 += kChunk) {
      const std::int64_t n = std::min(kChunk, Npart - ip0 + 1);
      buf.resize(static_cast<std::size_t>(6 * n));
      if (io) std::memcpy(buf.data(), s.xv_p(ip0), static_cast<std::size_t>(6 * n) * sizeof(float));
      if (ax::ParallelDescriptor::NProcs() > 1) ax::ParallelDescriptor::Bcast(buf.data(), 6 * n, root);
      // the local ones, per box, in index order (deterministic tiles)
      const int nlocal = pot.local_size();
      std::vector<std::vector<std::int64_t>> lists(nlocal);
      for (std::int64_t q = 0; q < n; ++q) {
        const float* p = &buf[static_cast<std::size_t>(6 * q)];
        const auto gt = assign(ax::IntVect(home(p[0]), home(p[1]), home(p[2])), ax::IntVect(0));
        if (gt.first < 0) throw std::runtime_error("pm_backend = amrex: a particle has no home box");
        if (dm[gt.first] == myproc) lists[pot.localindex(gt.first)].push_back(q);
      }
      for (ax::MFIter mfi(pot); mfi.isValid(); ++mfi) {
        const auto& L = lists[mfi.LocalIndex()];
        if (L.empty()) continue;
        auto& pt = pc->DefineAndReturnParticleTile(0, mfi.index(), mfi.LocalTileIndex());
        const std::size_t n0 = pt.numParticles();
        pt.resize(n0 + L.size());
        const PTD pd = pt.getParticleTileData();
        for (std::size_t q = 0; q < L.size(); ++q) {
          const float* p = &buf[static_cast<std::size_t>(6 * L[q])];
          const int n_ = static_cast<int>(n0 + q);
          for (int c = 0; c < 3; ++c) {
            pd.m_rdata[c][n_] = fpos(p[c]);
            pd.m_rdata[kX + c][n_] = p[c];
            pd.m_rdata[kV + c][n_] = p[3 + c];
          }
          pd.m_rdata[kM][n_] = 0.0f;
          pd.m_rdata[kA][n_] = 0.0f;
          pd.m_idata[kType][n_] = 0;
          pd.m_idcpu[n_] = ax::SetParticleIDandCPU(ip0 + L[q], myproc);
        }
        nmine += static_cast<Long>(L.size());
      }
    }
    ax::ParallelDescriptor::ReduceLongSum(nmine);
    if (nmine != Npart) throw std::runtime_error("pm_backend = amrex: the particle count does not add up");
    loaded = true;
    T.prep += now() - t0;
  }

  // calc_drift on the level: x = mod(x + v dt/2 + ndmr, ndmr) per particle
  // (the host arithmetic), then the home cells and a Redistribute
  void drift(double dt_sim, double ndmr) {
    sync_levels();  // the hydro may have regridded: the Redistribute below moves the particles onto the levels
    const double t0 = now();
    for_tiles([&](int l, const ax::MFIter& mfi) {
      auto& pt = tile(l, mfi);
      const int np = static_cast<int>(pt.numParticles());
      const PTD pd = pt.getParticleTileData();
#pragma omp parallel for schedule(static)
      for (int n = 0; n < np; ++n) {
        for (int c = 0; c < 3; ++c) {
          const double t = static_cast<double>(pd.m_rdata[kX + c][n]) +
                           static_cast<double>(pd.m_rdata[kV + c][n]) * dt_sim / 2 + ndmr;
          const float x = static_cast<float>(f::mod(t, ndmr));
          pd.m_rdata[kX + c][n] = x;
          pd.m_rdata[c][n] = fpos(x);
        }
      }
    });
    pc->Redistribute();
    T.drift += now() - t0;
  }

  // The level's particles into rank 0's xv by id (the mirror): every rank
  // sends (id, x, v) records to rank 0 in chunks. Collective.
  void gather_particles(pm::PMState& s) {
    const double t0 = now();
    struct Rec {
      Long id;
      float xv[6];
    };
    std::vector<Rec> mine;
    for_tiles([&](int l, const ax::MFIter& mfi) {
      auto& pt = tile(l, mfi);
      const int np = static_cast<int>(pt.numParticles());
      const PTD pd = pt.getParticleTileData();
      const std::size_t n0 = mine.size();
      mine.resize(n0 + np);
      std::size_t nn = n0;
      for (int n = 0; n < np; ++n) {
        if (pd.m_idata[kType][n] != 0) continue;  // the stars have no host mirror
        Rec& r = mine[nn++];
        r.id = static_cast<Long>(ax::ParticleIDWrapper(pd.m_idcpu[n]));
        for (int c = 0; c < 3; ++c) {
          r.xv[c] = pd.m_rdata[kX + c][n];
          r.xv[3 + c] = pd.m_rdata[kV + c][n];
        }
      }
      mine.resize(nn);
    });
    const bool io = ax::ParallelDescriptor::IOProcessor();
    if (io && static_cast<std::int64_t>(s.xv.size()) != 6 * Npart)
      throw std::runtime_error("amrex_state = 1: rank 0's xv is not allocated");
    auto unpack = [&](const Rec* r, std::size_t n) {
      for (std::size_t q = 0; q < n; ++q) {
        if (r[q].id < 1 || r[q].id > Npart) throw std::runtime_error("pm_backend = amrex: a particle id out of range");
        std::memcpy(s.xv_p(r[q].id), r[q].xv, sizeof r[q].xv);
      }
    };
    const int nranks = ax::ParallelDescriptor::NProcs();
    const int myproc = ax::ParallelDescriptor::MyProc();
    const int root = ax::ParallelDescriptor::IOProcessorNumber();
    if (nranks == 1) {
      unpack(mine.data(), mine.size());
    } else {
#ifdef AMREX_USE_MPI
      const MPI_Comm comm = ax::ParallelDescriptor::Communicator();
      constexpr std::size_t kChunk = std::size_t(1) << 24;  // records per message (512 MB)
      Long count = static_cast<Long>(mine.size());
      std::vector<Long> counts(nranks);
      MPI_Gather(&count, 1, MPI_LONG_LONG, counts.data(), 1, MPI_LONG_LONG, root, comm);
      if (io) unpack(mine.data(), mine.size());
      std::vector<Rec> rbuf;
      for (int r = 0; r < nranks; ++r) {
        if (r == root) continue;
        if (myproc == r) {
          for (std::size_t off = 0; off < mine.size(); off += kChunk) {
            const int n = static_cast<int>(std::min(kChunk, mine.size() - off));
            MPI_Send(mine.data() + off, n * static_cast<int>(sizeof(Rec)), MPI_BYTE, root, 7, comm);
          }
        } else if (io) {
          const std::size_t total = static_cast<std::size_t>(counts[r]);
          for (std::size_t off = 0; off < total; off += kChunk) {
            const int n = static_cast<int>(std::min(kChunk, total - off));
            rbuf.resize(static_cast<std::size_t>(n));
            MPI_Recv(rbuf.data(), n * static_cast<int>(sizeof(Rec)), MPI_BYTE, r, 7, comm, MPI_STATUS_IGNORE);
            unpack(rbuf.data(), static_cast<std::size_t>(n));
          }
        }
      }
#endif
    }
    T.gather += now() - t0;
  }

  // ---- the star-formation sink's particles (docs/stage_f.md section 2) ------

  // One star per sink event: x and v from level-0 gas-cell units to the
  // particle's (DM grid) units (x ndmr/ngasr), the mass as given, the home
  // cell, an id Npart + 1 + (the global count before) and a type 1; placed
  // on this rank's first tile and moved by a Redistribute. Collective
  // (the prefix count). Returns the global number added.
  Long add_stars(pm::PMState& s, const std::vector<SinkEvent>& ev, double a_form) {
    if (!state) throw std::runtime_error("the star-formation sink needs amrex_state = 1");
    const double t0 = now();
    const Long nmine = static_cast<Long>(ev.size());
    // the ids are independent of the rank count: every rank gathers this
    // call's event positions (one event per cell, so they are distinct),
    // sorts them, and numbers its own events by their rank in that order
    std::vector<float> pos(static_cast<std::size_t>(3 * nmine));
    for (Long q = 0; q < nmine; ++q)
      for (int c = 0; c < 3; ++c) pos[static_cast<std::size_t>(3 * q + c)] = ev[static_cast<std::size_t>(q)].x[c];
    std::vector<float> all = pos;
#ifdef AMREX_USE_MPI
    if (ax::ParallelDescriptor::NProcs() > 1) {
      const MPI_Comm comm = ax::ParallelDescriptor::Communicator();
      const int nranks = ax::ParallelDescriptor::NProcs();
      const int cnt = static_cast<int>(3 * nmine);
      std::vector<int> counts(nranks), displs(nranks);
      MPI_Allgather(&cnt, 1, MPI_INT, counts.data(), 1, MPI_INT, comm);
      int off = 0;
      for (int r = 0; r < nranks; ++r) {
        displs[r] = off;
        off += counts[r];
      }
      all.resize(static_cast<std::size_t>(off));
      MPI_Allgatherv(pos.data(), cnt, MPI_FLOAT, all.data(), counts.data(), displs.data(), MPI_FLOAT, comm);
    }
#endif
    const Long total = static_cast<Long>(all.size() / 3);
    if (total == 0) return 0;
    std::vector<Long> order(static_cast<std::size_t>(total));
    for (Long q = 0; q < total; ++q) order[static_cast<std::size_t>(q)] = q;
    auto less = [&](Long a, Long b) {
      for (int c = 0; c < 3; ++c) {
        const float xa = all[static_cast<std::size_t>(3 * a + c)], xb = all[static_cast<std::size_t>(3 * b + c)];
        if (xa != xb) return xa < xb;
      }
      return false;
    };
    std::sort(order.begin(), order.end(), less);
    std::vector<Long> rank_of(static_cast<std::size_t>(total));  // position in the sorted order
    for (Long q = 0; q < total; ++q) rank_of[static_cast<std::size_t>(order[static_cast<std::size_t>(q)])] = q;
    Long before = 0;  // this rank's events sit at [before, before + nmine) of `all`
#ifdef AMREX_USE_MPI
    if (ax::ParallelDescriptor::NProcs() > 1) {
      MPI_Exscan(&nmine, &before, 1, MPI_LONG_LONG, MPI_SUM, ax::ParallelDescriptor::Communicator());
      if (ax::ParallelDescriptor::MyProc() == 0) before = 0;
    }
#endif
    const double scale = s.pc.ndmr / s.pc.ngasr;
    const int myproc = ax::ParallelDescriptor::MyProc();
    if (nmine > 0) {
      ax::MFIter mfi(pot);
      if (!mfi.isValid()) throw std::runtime_error("add_stars: a rank with events but no box");
      auto& pt = pc->DefineAndReturnParticleTile(0, mfi.index(), mfi.LocalTileIndex());
      const std::size_t n0 = pt.numParticles();
      pt.resize(n0 + static_cast<std::size_t>(nmine));
      const PTD pd = pt.getParticleTileData();
      for (Long q = 0; q < nmine; ++q) {
        const SinkEvent& e = ev[static_cast<std::size_t>(q)];
        const int n_ = static_cast<int>(n0 + static_cast<std::size_t>(q));
        for (int c = 0; c < 3; ++c) {
          const float x = static_cast<float>(f::mod(static_cast<double>(e.x[c]) * scale, s.pc.ndmr));
          pd.m_rdata[kX + c][n_] = x;
          pd.m_rdata[kV + c][n_] = static_cast<float>(static_cast<double>(e.v[c]) * scale);
          pd.m_rdata[c][n_] = fpos(x);
        }
        pd.m_rdata[kM][n_] = e.mass;
        pd.m_rdata[kA][n_] = static_cast<float>(a_form);
        pd.m_idata[kType][n_] = 1;
        pd.m_idcpu[n_] = ax::SetParticleIDandCPU(Npart + 1 + s.nstar + rank_of[static_cast<std::size_t>(before + q)], myproc);
      }
    }
    s.nstar += total;
    pc->Redistribute();
    T.prep += now() - t0;
    return total;
  }

  // Every star into rank 0's s.stars (sorted by id). Collective.
  Long gather_stars(pm::PMState& s) {
    std::vector<pm::Star> mine;
    for_tiles([&](int l, const ax::MFIter& mfi) {
      auto& pt = tile(l, mfi);
      const int np = static_cast<int>(pt.numParticles());
      const PTD pd = pt.getParticleTileData();
      for (int n = 0; n < np; ++n) {
        if (pd.m_idata[kType][n] == 0) continue;
        pm::Star st;
        st.id = static_cast<Long>(ax::ParticleIDWrapper(pd.m_idcpu[n]));
        for (int c = 0; c < 3; ++c) {
          st.x[c] = pd.m_rdata[kX + c][n];
          st.v[c] = pd.m_rdata[kV + c][n];
        }
        st.mass = pd.m_rdata[kM][n];
        st.a_form = pd.m_rdata[kA][n];
        mine.push_back(st);
      }
    });
    const bool io = ax::ParallelDescriptor::IOProcessor();
    const int nranks = ax::ParallelDescriptor::NProcs();
    Long total = static_cast<Long>(mine.size());
    if (nranks == 1) {
      s.stars = std::move(mine);
    } else {
#ifdef AMREX_USE_MPI
      const MPI_Comm comm = ax::ParallelDescriptor::Communicator();
      const int root = ax::ParallelDescriptor::IOProcessorNumber();
      const int cnt = static_cast<int>(mine.size() * sizeof(pm::Star));
      std::vector<int> counts(nranks), displs(nranks);
      MPI_Gather(&cnt, 1, MPI_INT, counts.data(), 1, MPI_INT, root, comm);
      std::vector<pm::Star> all;
      if (io) {
        int off = 0;
        for (int r = 0; r < nranks; ++r) {
          displs[r] = off;
          off += counts[r];
        }
        all.resize(static_cast<std::size_t>(off) / sizeof(pm::Star));
      }
      MPI_Gatherv(mine.data(), cnt, MPI_BYTE, all.data(), counts.data(), displs.data(), MPI_BYTE, root, comm);
      if (io) s.stars = std::move(all);
      MPI_Allreduce(MPI_IN_PLACE, &total, 1, MPI_LONG_LONG, MPI_SUM, comm);
#endif
    }
    if (io) std::sort(s.stars.begin(), s.stars.end(), [](const pm::Star& a, const pm::Star& b) { return a.id < b.id; });
    return total;
  }

  // s.stars (the same list on every rank) onto the level: each rank keeps
  // the stars whose home cell lies in its boxes
  void load_stars(const pm::PMState& s) {
    const int myproc = ax::ParallelDescriptor::MyProc();
    const auto assign = locator.getGridAssignor();
    const double xdg = x_dm_m;
    const int NN = Nm;
    for (const pm::Star& st : s.stars) {
      const int h[3] = {f::mod(f::int4(static_cast<double>(st.x[0]) * xdg), NN),
                        f::mod(f::int4(static_cast<double>(st.x[1]) * xdg), NN),
                        f::mod(f::int4(static_cast<double>(st.x[2]) * xdg), NN)};
      const auto gt = assign(ax::IntVect(h[0], h[1], h[2]), ax::IntVect(0));
      if (gt.first < 0) throw std::runtime_error("pm_backend = amrex: a star has no home box");
      if (dm[gt.first] != myproc) continue;
      auto& pt = pc->DefineAndReturnParticleTile(0, gt.first, 0);
      const std::size_t n0 = pt.numParticles();
      pt.resize(n0 + 1);
      const PTD pd = pt.getParticleTileData();
      const int n_ = static_cast<int>(n0);
      for (int c = 0; c < 3; ++c) {
        pd.m_rdata[c][n_] = fpos(st.x[c]);
        pd.m_rdata[kX + c][n_] = st.x[c];
        pd.m_rdata[kV + c][n_] = st.v[c];
      }
      pd.m_rdata[kM][n_] = st.mass;
      pd.m_rdata[kA][n_] = st.a_form;
      pd.m_idata[kType][n_] = 1;
      pd.m_idcpu[n_] = ax::SetParticleIDandCPU(st.id, myproc);
    }
  }

  void prep(const pm::PMState& s, const float* u_gas, const float* v_gas, const float* e_gas, int kind) {
    const double t0 = now();
    const int myproc = ax::ParallelDescriptor::MyProc();
    const int nlocal = pot.local_size();
    if (state) {
      // amrex_state = 1: the particles are on the level already (load_particles,
      // drift); only the gas below
      if (!loaded) throw std::runtime_error("amrex_state = 1: amrex_pm_load_particles was not called");
      if (kind != 2) throw std::runtime_error("amrex_state = 1 needs hydro_scheme = unsplit");
      for (int l = 0; l <= finest(); ++l) {
        FFab& G = GAS(l);
        const FFab& S = l == 0 ? state->G : state->fine[static_cast<std::size_t>(l) - 1].G;
        for (ax::MFIter mfi(G); mfi.isValid(); ++mfi) {
          auto g = G.array(mfi);
          auto ga = S.const_array(mfi);
          ax::ParallelForOMP(mfi.validbox(), [=](int i, int j, int k) {
            for (int m = 0; m < 5; ++m) g(i, j, k, m) = ga(i, j, k, gs::kU + m);
            g(i, j, k, 5) = ga(i, j, k, gs::kE);
            g(i, j, k, 6) = g(i, j, k, 7) = 0.0f;
          });
        }
      }
      T.prep += now() - t0;
      return;
    }
    for (ax::MFIter mfi(pot); mfi.isValid(); ++mfi)
      pc->DefineAndReturnParticleTile(0, mfi.index(), mfi.LocalTileIndex()).resize(0);
    const auto assign = locator.getGridAssignor();
    const int nt = nthreads_max();
    std::vector<std::vector<std::vector<std::int64_t>>> lists(nt, std::vector<std::vector<std::int64_t>>(nlocal));
    const double xdg = x_dm_m;
    const int NN = Nm;
    Long bad = 0;
#pragma omp parallel reduction(+ : bad)
    {
      auto& mine = lists[thread_id()];
#pragma omp for schedule(static)
      for (std::int64_t ip = 1; ip <= Npart; ++ip) {
        const float* p = s.xv_p(ip);
        const int i = f::mod(f::int4(static_cast<double>(p[0]) * xdg), NN);
        const int j = f::mod(f::int4(static_cast<double>(p[1]) * xdg), NN);
        const int k = f::mod(f::int4(static_cast<double>(p[2]) * xdg), NN);
        const auto gt = assign(ax::IntVect(i, j, k), ax::IntVect(0));
        if (gt.first < 0) {
          bad += 1;
          continue;
        }
        if (dm[gt.first] == myproc) mine[pot.localindex(gt.first)].push_back(ip);
      }
    }
    if (bad) throw std::runtime_error("pm_backend = amrex: a particle has no home box");
    const double t1 = now();
    T.p_locate += t1 - t0;
    Long nmine = 0;
    for (ax::MFIter mfi(pot); mfi.isValid(); ++mfi) {
      const int b = mfi.LocalIndex();
      std::size_t n = 0;
      for (int t = 0; t < nt; ++t) n += lists[t][b].size();
      auto& pt = pc->DefineAndReturnParticleTile(0, mfi.index(), mfi.LocalTileIndex());
      pt.resize(n);
      const PTD pd = pt.getParticleTileData();
      std::size_t n0 = 0;
      for (int t = 0; t < nt; ++t) {
        const auto& L = lists[t][b];
#pragma omp parallel for schedule(static)
        for (std::size_t q = 0; q < L.size(); ++q) {
          const std::int64_t ip = L[q];
          const float* p = s.xv_p(ip);
          const int n_ = static_cast<int>(n0 + q);
          const int i = f::mod(f::int4(static_cast<double>(p[0]) * xdg), NN);
          const int j = f::mod(f::int4(static_cast<double>(p[1]) * xdg), NN);
          const int k = f::mod(f::int4(static_cast<double>(p[2]) * xdg), NN);
          pd.m_rdata[0][n_] = static_cast<float>(i) + 0.5f;
          pd.m_rdata[1][n_] = static_cast<float>(j) + 0.5f;
          pd.m_rdata[2][n_] = static_cast<float>(k) + 0.5f;
          for (int c = 0; c < 3; ++c) {
            pd.m_rdata[kX + c][n_] = p[c];
            pd.m_rdata[kV + c][n_] = p[3 + c];
          }
          pd.m_rdata[kM][n_] = 0.0f;
          pd.m_rdata[kA][n_] = 0.0f;
          pd.m_idata[kType][n_] = 0;
          pd.m_idcpu[n_] = ax::SetParticleIDandCPU(ip, myproc);
        }
        n0 += L.size();
      }
      nmine += static_cast<Long>(n);
    }
    ax::ParallelDescriptor::ReduceLongSum(nmine);
    if (nmine != Npart) throw std::runtime_error("pm_backend = amrex: the particle count does not add up");
    T.p_fill += now() - t1;
    // the gas (every rank holds the full host grid)
    const int NG = Ng;
    const bool mf_ = mf;
    for (ax::MFIter mfi(gas); mfi.isValid(); ++mfi) {
      auto g = gas.array(mfi);
      ax::ParallelForOMP(mfi.validbox(), [=](int i, int j, int k) {
        const std::int64_t gc = i + static_cast<std::int64_t>(NG) * (j + static_cast<std::int64_t>(NG) * k);
        for (int m = 0; m < 5; ++m) g(i, j, k, m) = u_gas[gc * 5 + m];
        if (kind == 2) {
          g(i, j, k, 5) = e_gas ? e_gas[gc] : 0.0f;
          g(i, j, k, 6) = g(i, j, k, 7) = 0.0f;
        } else if (mf_ && v_gas) {
          for (int m = 0; m < 3; ++m) g(i, j, k, 5 + m) = v_gas[gc * 3 + m];
        } else {
          g(i, j, k, 5) = g(i, j, k, 6) = g(i, j, k, 7) = 0.0f;
        }
      });
    }
    T.prep += now() - t0;
  }

  // ---- calc_density_field: CIC in int64 fixed point + the gas -------------

  // The CIC deposit of level l's particles into ACC(l) (its own mesh, in
  // density units: mass 8^l, with the level's fixed point), then
  // SumBoundary. bad counts the corners outside the ghost layer.
  void deposit_level(int l, double mass, Long& bad) {
    LFab& A = ACC(l);
    A.setVal(0);
    const int NN = Nm << l;
    const int bits = BITS(l);
    const double scale = std::ldexp(1.0, 3 * l);  // mass per level-l cell -> density per level-0 mesh cell
    const double r3 = static_cast<double>(r) * r * r;
    const rpk::MeshGeom g = l == 0 ? mg : mesh_geom(l);
    for (ax::MFIter mfi(A); mfi.isValid(); ++mfi) {
      const ax::Box bx = mfi.validbox();
      const ax::Box gbx = ax::grow(bx, 1);
      const auto glo = ax::lbound(gbx), ghi = ax::ubound(gbx), glen = ax::length(gbx);
      const Long ngc = gbx.numPts();
      auto& pt = tile(l, mfi);
      const int np = static_cast<int>(pt.numParticles());
      const PTD pd = pt.getParticleTileData();
      const int nt = nthreads_max();
      std::vector<Long> buf(static_cast<std::size_t>(nt) * ngc, 0);
      auto a = A.array(mfi);
#pragma omp parallel reduction(+ : bad)
      {
        Long* mine = buf.data() + static_cast<std::size_t>(thread_id()) * ngc;
#pragma omp for schedule(static)
        for (int n = 0; n < np; ++n) {
          const float x[3] = {pd.m_rdata[kX][n], pd.m_rdata[kX + 1][n], pd.m_rdata[kX + 2][n]};
          const int h[3] = {home_l(pd.m_rdata[0][n], l), home_l(pd.m_rdata[1][n], l), home_l(pd.m_rdata[2][n], l)};
          // dark matter: the particle mass; a star: its gas mass (mean-matter
          // x level-0 gas-cell volumes) spread over the r^3 mesh cells of a
          // gas cell, i.e. mass_gas r^3 per mesh cell of full weight
          const double m = (pd.m_idata[kType][n] == 0 ? mass : static_cast<double>(pd.m_rdata[kM][n]) * r3) * scale;
          const rpk::Cic c = rpk::cic_weights(g, x);
          rpk::cic_corners(c, [&](int i, int j, int k, double w) {
            const int ii = near_image(i - 1, h[0], NN), jj = near_image(j - 1, h[1], NN),
                      kk = near_image(k - 1, h[2], NN);
            if (ii < glo.x || ii > ghi.x || jj < glo.y || jj > ghi.y || kk < glo.z || kk > ghi.z) {
              bad += 1;
              return;
            }
            mine[(ii - glo.x) + glen.x * ((jj - glo.y) + static_cast<Long>(glen.y) * (kk - glo.z))] +=
                to_fx(w * m, bits);
          });
        }
      }
      ax::ParallelForOMP(gbx, [=, &buf](int i, int j, int k) {
        const Long idx = (i - glo.x) + glen.x * ((j - glo.y) + static_cast<Long>(glen.y) * (k - glo.z));
        Long sum = 0;
        for (int t = 0; t < nt; ++t) sum += buf[static_cast<std::size_t>(t) * ngc + idx];
        a(i, j, k) += sum;
      });
    }
    const double t1 = now();
    A.SumBoundary(PER(l));  // fold the ghost deposits into their owners (int64: exact)
    T.sumbnd += now() - t1;
  }

  // pm_level_solve: level l's deposits restricted onto its parent and added
  // into TOT(l - 1): the total of its valid cells and, in its uncovered
  // ghost cells, the CIC mass of its own particles that fell outside the
  // level. The restriction is the [1, 3, 3, 1]/8 stencil per dimension
  // (the parent cell's two children and their outer neighbours), not the
  // plain mean of the 8 children: a fine CIC cloud restricted by the mean
  // is a sharper source than the parent's own CIC (NGP-like at the parent's
  // resolution), whose 7-point Laplacian then gives 5-8 % force errors 8
  // parent cells away (the lattice anisotropy that the CIC smoothing
  // cancels); the parent CIC's spectrum is the fine one's times
  // cos^2(k h_f / 2), i.e. the fine field filtered with [1, 2, 1]/4, and
  // sampling that at the parent centres (half a fine cell off the fine
  // centres) with the two-point mean gives [1, 3, 3, 1]/8: the parent sees
  // every fine level's mass as its own CIC would, times one more
  // half-fine-cell top-hat. Mass is conserved (each fine cell's weights
  // sum to 1/8 in 3-D: the mean density). Integer arithmetic: exact and
  // layout-independent. Collective.
  void restrict_and_spill(int l) {
    FineMesh& M = fm[static_cast<std::size_t>(l) - 1];
    const int db = BITS(l - 1) - BITS(l);  // 0 or 3
    const Long mul = Long(1) << db;
    LFab cr(ax::coarsen(M.ba, 2), M.dm, 1, 1);
    cr.setVal(0);
    for (ax::MFIter mfi(cr); mfi.isValid(); ++mfi) {
      const ax::Box cbx = mfi.validbox();
      const ax::Box fbx = M.ba[mfi.index()];
      const ax::Box fgb = ax::grow(fbx, 1);
      auto c = cr.array(mfi);
      auto t = M.tot.const_array(mfi);
      auto a = M.acc.const_array(mfi);
      auto cv = M.cov.const_array(mfi);
      ax::ParallelForOMP(ax::grow(cbx, 1), [=](int i, int j, int k) {
        static constexpr int kw[4] = {1, 3, 3, 1};
        i128 sum = 0;
        for (int dk = -1; dk <= 2; ++dk)
          for (int dj = -1; dj <= 2; ++dj)
            for (int di = -1; di <= 2; ++di) {
              const ax::IntVect f(2 * i + di, 2 * j + dj, 2 * k + dk);
              Long v;
              if (fbx.contains(f))
                v = t(f);  // this box's own cells: the level's total there
              else if (fgb.contains(f) && cv(f) == 0)
                v = a(f);  // outside the level: this box's particles' spill
              else
                continue;  // another box's cell (it adds it) or beyond the reach of any CIC cloud
              sum += static_cast<i128>(v) * (kw[di + 1] * kw[dj + 1] * kw[dk + 1]);
            }
        c(i, j, k) = static_cast<Long>((sum * mul) / 512);
      });
    }
    TOT(l - 1).ParallelAdd(cr, 0, 0, 1, 1, 0, PER(l - 1));
  }

  // pm_level_solve: the density of level l into POT(l): the level's gas,
  // its total deposit and, below level 0, the parent level's own deposit
  // prolonged (piecewise constant), rounded to float once
  void rho_level(int l) {
    FFab& P = POT(l);
    const LFab& Tt = TOT(l);
    const FFab& G = GAS(l);
    const int rr = r;
    const int bits = BITS(l);
    if (l == 0) {
      for (ax::MFIter mfi(P); mfi.isValid(); ++mfi) {
        auto p = P.array(mfi);
        auto t = Tt.const_array(mfi);
        auto gsa = G.const_array(mfi);
        ax::ParallelForOMP(mfi.validbox(), [=](int i, int j, int k) {
          const double rho = gsa(i / rr, j / rr, k / rr, 0);
          p(i, j, k) = static_cast<float>(rho + std::ldexp(static_cast<double>(t(i, j, k)), -bits));
        });
      }
      return;
    }
    FineMesh& M = fm[static_cast<std::size_t>(l) - 1];
    LFab po(ax::coarsen(M.ba, 2), M.dm, 1, 0);
    po.setVal(0);
    po.ParallelCopy(OWN(l - 1), 0, 0, 1, 0, 0, PER(l - 1));
    const int bp = BITS(l - 1);
    for (ax::MFIter mfi(P); mfi.isValid(); ++mfi) {
      auto p = P.array(mfi);
      auto t = Tt.const_array(mfi);
      auto o = po.const_array(mfi);
      auto gsa = G.const_array(mfi);
      ax::ParallelForOMP(mfi.validbox(), [=](int i, int j, int k) {
        const double rho = gsa(i / rr, j / rr, k / rr, 0);
        p(i, j, k) = static_cast<float>(rho + std::ldexp(static_cast<double>(t(i, j, k)), -bits) +
                                        std::ldexp(static_cast<double>(o(i >> 1, j >> 1, k >> 1)), -bp));
      });
    }
  }

  void density(double mass) {
    const double t0 = now();
    Long bad = 0;
    const int nf = finest();
    for (int l = 0; l <= nf; ++l) {
      if (l == 1) T.deposit += now() - t0;
      // the fixed point's range: a particle of full weight must fit the level
      if (l > 0 && mass * std::ldexp(1.0, 3 * l) > std::ldexp(1.0, 62 - BITS(l)))
        throw std::runtime_error("pm_level_solve: the particle mass overflows level " + std::to_string(l) +
                                 "'s fixed point (lighter particles or fewer levels)");
      deposit_level(l, mass, bad);
    }
    if (bad) throw std::runtime_error("pm_backend = amrex: a CIC corner fell outside the ghost layer");
    if (level_solve) {
      // own and total deposits per level, the finer levels restricted onto
      // their parents, the densities, the mean of level 0's (fixed point:
      // exact and layout-independent)
      const double t2 = now();
      for (int l = 0; l <= nf; ++l) {
        ax::Copy(OWN(l), ACC(l), 0, 0, 1, 0);
        ax::Copy(TOT(l), ACC(l), 0, 0, 1, 0);
      }
      for (int l = nf; l >= 1; --l) restrict_and_spill(l);
      for (int l = 0; l <= nf; ++l) rho_level(l);
      i128 S = 0;
      for (ax::MFIter mfi(pot); mfi.isValid(); ++mfi) {
        auto p = pot.const_array(mfi);
        const ax::Box bx = mfi.validbox();
        const int nt = nthreads_max();
        std::vector<i128> part(static_cast<std::size_t>(nt), 0);
#pragma omp parallel
        {
          i128 mine = 0;
          const auto lo = ax::lbound(bx), hi = ax::ubound(bx);
#pragma omp for schedule(static) collapse(2)
          for (int k = lo.z; k <= hi.z; ++k)
            for (int j = lo.y; j <= hi.y; ++j)
              for (int i = lo.x; i <= hi.x; ++i) mine += to_fx128(static_cast<double>(p(i, j, k)), kDensBits);
          part[static_cast<std::size_t>(thread_id())] = mine;
        }
        for (int t = 0; t < nt; ++t) S += part[static_cast<std::size_t>(t)];
      }
      allreduce_i128(&S, 1);
      rho_mean = from_fx128(S, kDensBits) / (static_cast<double>(Nm) * Nm * Nm);
      if (std::getenv("RADRAY_PM_LEVEL_DEBUG")) {
        // the mass of every level's deposits (density x cell volume): the
        // total on level 0 must be the particles' (DM + stars)
        for (int l = 0; l <= nf; ++l) {
          i128 so = 0, st = 0;
          for (ax::MFIter mfi(TOT(l)); mfi.isValid(); ++mfi) {
            auto t = TOT(l).const_array(mfi);
            auto o = OWN(l).const_array(mfi);
            ax::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
              st += t(i, j, k);
              so += o(i, j, k);
            });
          }
          i128 v[2] = {so, st};
          allreduce_i128(v, 2);
          if (ax::ParallelDescriptor::IOProcessor())
            std::fprintf(stderr, "PM level debug: level %d: own mass %.9e, total mass %.9e (mesh cells %lld, mean density %.9e)\n",
                         l, from_fx128(v[0], BITS(l)) / std::ldexp(1.0, 3 * l),
                         from_fx128(v[1], BITS(l)) / std::ldexp(1.0, 3 * l), static_cast<long long>(BA(l).numPts()),
                         rho_mean);
        }
      }
      T.d_lev += now() - t2;
      if (nf == 0) T.deposit += now() - t0;
      return;
    }
    // density = the gas density of the gas cell the mesh cell lies in + the
    // DM mass (the original: float, added particle by particle; here one
    // rounding of the exact sum)
    const int rr = r;
    if (fine_levels() > 0) {
      // pm_gas_levels: the gas cell is the finest level's cell holding the
      // mesh cell (composite_density), the same rounding of the sum
      FFab gm(ba, dm, 1, 0);
      composite_density(gm);
      for (ax::MFIter mfi(pot); mfi.isValid(); ++mfi) {
        auto p = pot.array(mfi);
        auto a = acc.const_array(mfi);
        auto gc = gm.const_array(mfi);
        ax::ParallelForOMP(mfi.validbox(), [=](int i, int j, int k) {
          const double rho = gc(i, j, k);
          p(i, j, k) = static_cast<float>(rho + std::ldexp(static_cast<double>(a(i, j, k)), -kDensBits));
        });
      }
      T.deposit += now() - t0;
      return;
    }
    for (ax::MFIter mfi(pot); mfi.isValid(); ++mfi) {
      auto p = pot.array(mfi);
      auto a = acc.const_array(mfi);
      auto gs = gas.const_array(mfi);
      ax::ParallelForOMP(mfi.validbox(), [=](int i, int j, int k) {
        const double rho = gs(i / rr, j / rr, k / rr, 0);
        p(i, j, k) = static_cast<float>(rho + std::ldexp(static_cast<double>(a(i, j, k)), -kDensBits));
      });
    }
    T.deposit += now() - t0;
  }

  // ---- pm_gas_levels: the refined levels of the gas state -------------------
  //
  // A level-l cell spans q = r / 2^l mesh cells per dimension when l <= lr
  // (log2 r), or lies inside one mesh cell (2^(l - lr) level cells per mesh
  // cell) when l > lr. The deposit gives every mesh cell the density of
  // the finest level with l <= lr covering it: a level finer than the mesh
  // is already in its average-down onto level lr (the hydro and the sink
  // average every step down, and the RT and the PM do the same after
  // changing the fine cells). The force on a level-l cell is the mean mesh
  // force over its q^3 mesh cells, or the force of the mesh cell holding
  // it, divided by r (mesh to level-0 cell units, the units of every
  // level's velocities). Spread and average are again transposes of each
  // other, so the gas-DM forces stay pairwise antisymmetric (total momentum
  // conserved to round-off, docs/pm_refine.md section 1). Afterwards the
  // kicked levels are averaged down onto their parents (momenta and E),
  // so level 0 is the composite and the hydro's push-down sees no change.

  // The refined levels in use (0 without pm_gas_levels)
  int fine_levels() const { return gas_levels && state ? state->finest : 0; }

  // gm (mesh layout, 1 component): the gas density of the finest level
  // (<= lr) covering each mesh cell
  void composite_density(FFab& gm) {
    const double t0 = now();
    const int rr = r;
    for (ax::MFIter mfi(gm); mfi.isValid(); ++mfi) {
      auto o = gm.array(mfi);
      auto gs = gas.const_array(mfi);
      ax::ParallelForOMP(mfi.validbox(), [=](int i, int j, int k) { o(i, j, k) = gs(i / rr, j / rr, k / rr, 0); });
    }
    const int nl = std::min(fine_levels(), lr);
    for (int l = 1; l <= nl; ++l) {
      const gs::GasState::FineLevel& F = state->fine[static_cast<std::size_t>(l) - 1];
      const int q = r >> l;
      ax::BoxArray bq = F.ba;
      bq.refine(q);
      FFab t(bq, F.dm, 1, 0);
      for (ax::MFIter mfi(t); mfi.isValid(); ++mfi) {
        auto o = t.array(mfi);
        auto g = F.G.const_array(mfi);
        ax::ParallelForOMP(mfi.validbox(), [=](int i, int j, int k) { o(i, j, k) = g(i / q, j / q, k / q, gs::kU); });
      }
      gm.ParallelCopy(t, 0, 0, 1);  // the finer level last: it overwrites its parent
    }
    T.d_fine += now() - t0;
  }

  // The unsplit kick of the refined levels by the mesh force over their own
  // cells (kind 2 of kick_gas), then the average-down of u(2:5) onto the
  // parents, finest first. After gather (level 0's own kick is in the
  // state). Collective.
  void kick_fine(double dt_sim, const Units& u, double eta1) {
    const double t0 = now();
    const int nf = fine_levels();
    const double Dmin = u.Dgas_min;
    for (int l = 1; l <= nf; ++l) {
      gs::GasState::FineLevel& F = state->fine[static_cast<std::size_t>(l) - 1];
      const bool spans = l <= lr;
      const int q = spans ? r >> l : 1;           // mesh cells per level cell (l <= lr)
      const int sh = spans ? 0 : l - lr;          // level cells per mesh cell 2^sh (l > lr)
      ax::BoxArray bq = F.ba;
      if (spans)
        bq.refine(q);
      else
        bq.coarsen(1 << sh);
      FFab ft(bq, F.dm, 3, 0);
      ft.ParallelCopy(frc, 0, 0, 3);
      const double norm = static_cast<double>(q) * q * q * r;
      for (ax::MFIter mfi(F.G); mfi.isValid(); ++mfi) {
        auto g = F.G.array(mfi);
        auto fa = ft.const_array(mfi);
        ax::ParallelForOMP(mfi.validbox(), [=](int i, int j, int k) {
          double Fc[3] = {0, 0, 0};
          if (spans) {
            for (int kk = k * q; kk < (k + 1) * q; ++kk)
              for (int jj = j * q; jj < (j + 1) * q; ++jj)
                for (int ii = i * q; ii < (i + 1) * q; ++ii)
                  for (int m = 0; m < 3; ++m) Fc[m] += static_cast<double>(fa(ii, jj, kk, m));
          } else {
            for (int m = 0; m < 3; ++m) Fc[m] = fa(i >> sh, j >> sh, k >> sh, m);
          }
          for (int m = 0; m < 3; ++m) Fc[m] /= norm;
          double U[5];
          for (int m = 0; m < 5; ++m) U[m] = g(i, j, k, gs::kU + m);
          const double eD = g(i, j, k, gs::kE);
          const double eth = rrh::thermal_energy(U, eD, eta1);
          const double rho = f::max(U[0], Dmin);
          for (int m = 0; m < 3; ++m) U[1 + m] = U[1 + m] + rho * Fc[m] * dt_sim;
          U[4] = rrh::kinetic_energy(U) + eth;
          for (int m = 1; m < 5; ++m) g(i, j, k, gs::kU + m) = static_cast<float>(U[m]);
        });
      }
    }
    average_down_momenta();
    T.k_fine += now() - t0;
  }

  // pm_level_solve: the unsplit kick of every refined level's gas by its
  // own mesh force (gas_force over its r^3 mesh cells), then the
  // average-down onto the parents. After gather (level 0's kick is in the
  // state). Collective.
  void kick_levels(double dt_sim, const Units& u, double eta1) {
    const double t0 = now();
    const double Dmin = u.Dgas_min;
    for (int l = 1; l <= finest(); ++l) {
      gs::GasState::FineLevel& F = state->fine[static_cast<std::size_t>(l) - 1];
      const FFab& Fr = FRC(l);
      for (ax::MFIter mfi(F.G); mfi.isValid(); ++mfi) {
        auto g = F.G.array(mfi);
        auto fa = Fr.const_array(mfi);
        ax::ParallelForOMP(mfi.validbox(), [=, this](int i, int j, int k) {
          double Fc[3];
          gas_force(fa, i, j, k, Fc);
          double U[5];
          for (int m = 0; m < 5; ++m) U[m] = g(i, j, k, gs::kU + m);
          const double eD = g(i, j, k, gs::kE);
          const double eth = rrh::thermal_energy(U, eD, eta1);
          const double rho = f::max(U[0], Dmin);
          for (int m = 0; m < 3; ++m) U[1 + m] = U[1 + m] + rho * Fc[m] * dt_sim;
          U[4] = rrh::kinetic_energy(U) + eth;
          for (int m = 1; m < 5; ++m) g(i, j, k, gs::kU + m) = static_cast<float>(U[m]);
        });
      }
    }
    average_down_momenta();
    T.k_lev += now() - t0;
  }

  // u(2:5) of every refined level averaged onto its parent's covered cells
  // (the plain mean of the 8 children in double, rounded once, as the RT
  // hierarchy's average_down), finest first, so level 0 is the composite;
  // then the hydro is told not to push the level-0 kick down again
  void average_down_momenta() {
    const int nf = state->finest;
    constexpr int nc = 4;  // u(2:5)
    for (int l = nf; l >= 1; --l) {
      const gs::GasState::FineLevel& F = state->fine[static_cast<std::size_t>(l) - 1];
      FFab& PG = l == 1 ? state->G : state->fine[static_cast<std::size_t>(l) - 2].G;
      const ax::BoxArray& pba = l == 1 ? state->ba : state->fine[static_cast<std::size_t>(l) - 2].ba;
      const ax::DistributionMapping& pdm = l == 1 ? state->dm : state->fine[static_cast<std::size_t>(l) - 2].dm;
      FFab Fc(ax::coarsen(F.ba, 2), F.dm, nc + 1, 0);  // the averages and a "covered" flag
      for (ax::MFIter mfi(Fc); mfi.isValid(); ++mfi) {
        auto a = Fc.array(mfi);
        auto g = F.G.const_array(mfi);
        ax::ParallelForOMP(mfi.validbox(), [=](int i, int j, int k) {
          for (int m = 0; m < nc; ++m) {
            double sum = 0;
            for (int kk = 0; kk < 2; ++kk)
              for (int jj = 0; jj < 2; ++jj)
                for (int ii = 0; ii < 2; ++ii) sum += g(2 * i + ii, 2 * j + jj, 2 * k + kk, gs::kU + 1 + m);
            a(i, j, k, m) = static_cast<float>(sum / 8);
          }
          a(i, j, k, nc) = 1.0f;
        });
      }
      FFab Pc(pba, pdm, nc + 1, 0);
      Pc.setVal(0.0f);
      Pc.ParallelCopy(Fc, 0, 0, nc + 1);
      for (ax::MFIter mfi(PG); mfi.isValid(); ++mfi) {
        auto g = PG.array(mfi);
        auto a = Pc.const_array(mfi);
        ax::ParallelForOMP(mfi.validbox(), [=](int i, int j, int k) {
          if (a(i, j, k, nc) == 0.0f) return;
          for (int m = 0; m < nc; ++m) g(i, j, k, gs::kU + 1 + m) = a(i, j, k, m);
        });
      }
    }
    amrex_hydro_level0_is_current();  // the hydro must not push the level-0 kick down again
  }

  // RADRAY_PM_MOMENTUM_CHECK=1 (amrex_state = 1; a diagnostic for the tests):
  // the total momentum of the particles and of the gas before and after the
  // kicks, in mesh-cell masses x gas-grid velocities. The gas total is level
  // 0's, the composite (every refined level is averaged down onto it before
  // and after the PM). The DM-gas and self forces are pairwise
  // antisymmetric, so the change must vanish to the float round-off of the
  // stored velocities; a deposit and a kick that do not use the same cells
  // break it. Printed by the I/O rank: the change relative to the summed
  // |change| of the gas cells. Order-dependent sums: not bitwise.
  void momentum(double out[6], double mass) {
    for (int m = 0; m < 6; ++m) out[m] = 0;
    const double xg = x_dm_gas, r3 = static_cast<double>(r) * r * r;
    for_tiles([&](int l, const ax::MFIter& mfi) {
      auto& pt = tile(l, mfi);
      const int np = static_cast<int>(pt.numParticles());
      const PTD pd = pt.getParticleTileData();
      for (int n = 0; n < np; ++n) {
        const double mp = pd.m_idata[kType][n] == 0 ? mass : static_cast<double>(pd.m_rdata[kM][n]) * r3;
        for (int m = 0; m < 3; ++m) out[m] += mp * xg * pd.m_rdata[kV + m][n];
      }
    });
    for (ax::MFIter mfi(state->G); mfi.isValid(); ++mfi) {
      auto g = state->G.const_array(mfi);
      const ax::Box bx = mfi.validbox();
      ax::LoopOnCpu(bx, [&](int i, int j, int k) {
        for (int m = 0; m < 3; ++m) out[3 + m] += r3 * g(i, j, k, gs::kU + 1 + m);
      });
    }
    ax::ParallelDescriptor::ReduceRealSum(out, 6);
  }
  double gas_change_abs(const FFab& before) const {
    double s = 0;
    const double r3 = static_cast<double>(r) * r * r;
    for (ax::MFIter mfi(state->G); mfi.isValid(); ++mfi) {
      auto g = state->G.const_array(mfi);
      auto b = before.const_array(mfi);
      ax::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
        for (int m = 0; m < 3; ++m) s += r3 * std::fabs(static_cast<double>(g(i, j, k, gs::kU + 1 + m)) - b(i, j, k, m));
      });
    }
    ax::ParallelDescriptor::ReduceRealSum(s);
    return s;
  }

  // ---- calc_gravity_field: the Poisson solve --------------------------------

  void gravity(double scalefactor) {
    const double t0 = now();
    const int NN = Nm;
    const double Ak = 2 * pc::pi / static_cast<double>(Nm);
    // Aphi = -4*pi*scalefactor/(6*pi)  (unary minus applies to the quotient)
    const double Aphi = -(4 * pc::pi * scalefactor / (6 * pc::pi));
    auto post = [=](int i, int j, int k, ax::GpuComplex<float>& c) {
      const double wphi = rpk::green_wphi(NN, i, j, k, Ak, Aphi);
      c.m_real = static_cast<float>(static_cast<double>(c.m_real) * wphi);
      c.m_imag = static_cast<float>(static_cast<double>(c.m_imag) * wphi);
    };
    if (fft_kind == 1) {
      if (!afft) afft = std::make_unique<ax::FFT::R2C<float>>(domain);
      afft->forwardThenBackward(pot, pot, post);
      // amrex::FFT sets FFTW's planner to its thread count when AMReX uses
      // the threaded FFTW: restore the serial plans of the host code (hydro,
      // the serial PM). Only then: AMReX built with AMReX_FFTW_IGNORE_OMP
      // (klone) never initialises FFTW's threads, and plan_with_nthreads
      // before init_threads runs fftwf_cleanup, which leaves the solver's
      // plans undefined (their destruction crashed at 1024^3)
      if (fftwf_planner_nthreads() > 1) fftwf_plan_with_nthreads(1);
    } else {
      if (!pfft) pfft = std::make_unique<PencilFFT>(Nm);
      pfft->forward_then_backward(pot, pot, post);
    }
    // bscale = 1E0/(int(n1,8)*int(n2,8)*int(n3,8)) on the real outputs
    const float bscale = 1.0f / static_cast<float>(static_cast<std::int64_t>(Nm) * Nm * Nm);
    for (ax::MFIter mfi(pot); mfi.isValid(); ++mfi) {
      auto p = pot.array(mfi);
      ax::ParallelForOMP(mfi.validbox(), [=](int i, int j, int k) { p(i, j, k) = p(i, j, k) * bscale; });
    }
    pot.FillBoundary(period);
    T.fft += now() - t0;
    for (int l = 1; l <= finest(); ++l) solve_level(l, scalefactor);
  }

  // ---- pm_level_solve: the Poisson solve of a refined level ------------------
  //
  // lap phi = (2a/3)(rho - mean) on the level's mesh (dx = 2^-l), the
  // 7-point Laplacian that the FFT's Green's function inverts on level 0,
  // with MLMG: Dirichlet values on the level's boundary from the parent's
  // potential (setCoarseFineBC), periodic domain, the parent's potential
  // interpolated as the initial guess, mlmg_iters V-cycles (or until the
  // residual is below mlmg_tol of the initial one), red-black Gauss-Seidel
  // and the smoother as the bottom solver: every operation is per cell and
  // the norms are maxima, so the result is bitwise rank-, thread- and
  // layout-independent (the multigrid depth is fixed by sync_levels). The
  // potential's ghost cells outside the level take the parent's potential
  // interpolated (the 4-point gradient of the boundary cells and the force
  // of the ghost layer the kicks read).
  void solve_level(int l, double scalefactor) {
    const double t0 = now();
    FineMesh& M = fm[static_cast<std::size_t>(l) - 1];
    const int lp = l - 1;
    FFab pc2(ax::coarsen(M.ba, 2), M.dm, 1, 3);
    pc2.setVal(0.0f);
    pc2.ParallelCopy(POT(lp), 0, 0, 1, 0, 3, PER(lp));
    ax::MultiFab sol(M.ba, M.dm, 1, 1), rhs(M.ba, M.dm, 1, 0);
    sol.setVal(0.0);
    const double c = 2 * scalefactor / 3;  // -Aphi
    const double mean = rho_mean;
    for (ax::MFIter mfi(rhs); mfi.isValid(); ++mfi) {
      auto b = rhs.array(mfi);
      auto sa = sol.array(mfi);
      auto p = M.pot.const_array(mfi);
      auto q = pc2.const_array(mfi);
      ax::ParallelForOMP(mfi.validbox(), [=](int i, int j, int k) {
        b(i, j, k) = c * (static_cast<double>(p(i, j, k)) - mean);
        sa(i, j, k) = interp_parent(q, i, j, k);
      });
    }
    ax::MultiFab phic(BA(lp), DM(lp), 1, 0);
    {
      const FFab& PP = POT(lp);
      for (ax::MFIter mfi(phic); mfi.isValid(); ++mfi) {
        auto d = phic.array(mfi);
        auto p = PP.const_array(mfi);
        ax::ParallelForOMP(mfi.validbox(), [=](int i, int j, int k) { d(i, j, k) = p(i, j, k); });
      }
    }
    ax::LPInfo info;
    info.setMaxCoarseningLevel(M.mg_depth);
    ax::MLPoisson mlp({M.geom}, {M.ba}, {M.dm}, info);
    const auto P = ax::LinOpBCType::Periodic;
    mlp.setDomainBC({P, P, P}, {P, P, P});
    mlp.setCoarseFineBC(&phic, 2);
    mlp.setLevelBC(0, nullptr);
    mlp.setEnforceSingularSolvable(false);
    ax::MLMG mlmg(mlp);
    mlmg.setVerbose(mlmg_verbose);
    mlmg.setBottomVerbose(0);
    mlmg.setFixedIter(mlmg_iters);
    mlmg.setBottomSolver(ax::BottomSolver::smoother);
    // the smoother bottom solve: this many red-black sweeps on the coarsest
    // grid, at least the square of its width (a 9^3 bottom grid stalled the
    // V-cycles at a factor 0.7 each with AMReX's 8 sweeps; 64 gave 20x per
    // cycle; the sweeps one cycle needs grow as width^2), and otherwise
    // pm_mlmg_bottom_sweeps (8, AMReX's default: a bottom grid of width 2
    // needs no more, and every extra sweep is a boundary exchange; the 64 of
    // before cost more than it saved on 24 ranks), docs/stage_g.md 5.1
    M.sweeps = std::max(mlmg_bottom_sweeps,
                        static_cast<int>(std::ceil(mlmg_bottom_scale * M.bottom * M.bottom)));
    mlmg.setFinalSmooth(M.sweeps);
    mlmg.solve({&sol}, {&rhs}, mlmg_tol, 0.0);
    M.iters = mlmg.getNumIters();
    const double r0 = mlmg.getInitResidual();
    M.resid = r0 > 0 ? mlmg.getFinalResidual() / r0 : 0.0;
    // a level that ran out of V-cycles far from the target: say so (the
    // first time per level, then every 100th), the force is then only as
    // good as that residual
    if (M.iters >= mlmg_iters && M.resid > kCapResid) {
      M.capped += 1;
      if (ax::ParallelDescriptor::IOProcessor() && (M.capped == 1 || M.capped % 100 == 0))
        std::fprintf(stderr,
                     "pm_level_solve WARNING: level %d (%lld cells, multigrid bottom width %d) reached pm_mlmg_iters = %d "
                     "at relative residual %.2e > %.0e (%lld such solves); raise pm_mlmg_iters or "
                     "pm_mlmg_bottom_sweeps, or check the level's shape (docs/stage_g.md section 5.1)\n",
                     l, static_cast<long long>(M.ba.numPts()), M.bottom, mlmg_iters, M.resid, kCapResid,
                     static_cast<long long>(M.capped));
    }
    for (ax::MFIter mfi(M.pot); mfi.isValid(); ++mfi) {
      auto p = M.pot.array(mfi);
      auto sa = sol.const_array(mfi);
      ax::ParallelForOMP(mfi.validbox(), [=](int i, int j, int k) { p(i, j, k) = static_cast<float>(sa(i, j, k)); });
    }
    M.pot.FillBoundary(M.period);
    for (ax::MFIter mfi(M.pot); mfi.isValid(); ++mfi) {
      auto p = M.pot.array(mfi);
      auto q = pc2.const_array(mfi);
      auto cv = M.cov.const_array(mfi);
      ax::ParallelForOMP(ax::grow(mfi.validbox(), 3), [=](int i, int j, int k) {
        if (cv(i, j, k) == 0) p(i, j, k) = static_cast<float>(interp_parent_quad(q, i, j, k));
      });
    }
    T.solve += now() - t0;
  }

  // ---- the mesh force on the box and its ghost layer ------------------------

  // Level l: the stencil's differences are per level-l cell; divided by
  // dx = 2^-l (an exact float scaling) the force is in level-0 mesh units
  // on every level
  void force_level(int l) {
    const int NN = Nm << l;
    FFab& Fr = FRC(l);
    const FFab& P = POT(l);
    const float scale = static_cast<float>(1 << l);
    for (ax::MFIter mfi(Fr); mfi.isValid(); ++mfi) {
      auto f = Fr.array(mfi);
      auto p = P.const_array(mfi);
      ax::ParallelForOMP(ax::grow(mfi.validbox(), 1), [=](int i, int j, int k) {
        double F[3];
        rpk::force_stencil(
            NN, i + 1, j + 1, k + 1,
            [&](int ii, int jj, int kk) -> float {
              return p(near_image(ii - 1, i, NN), near_image(jj - 1, j, NN), near_image(kk - 1, k, NN));
            },
            F);
        if (l == 0)
          for (int m = 0; m < 3; ++m) f(i, j, k, m) = static_cast<float>(F[m]);  // exact: F is a float
        else
          for (int m = 0; m < 3; ++m) f(i, j, k, m) = static_cast<float>(F[m]) * scale;
      });
    }
  }

  void force() {
    const double t0 = now();
    force_level(0);
    T.force += now() - t0;
    const double t1 = now();
    for (int l = 1; l <= finest(); ++l) force_level(l);
    T.f_lev += now() - t1;
  }

  // ---- calc_kick: accelerate_dm ---------------------------------------------

  // The particles of level l kicked by its mesh force (the CIC at the
  // level's resolution); the statistics accumulated into the sums and maxima
  void kick_dm_level(int l, double dt_sim, i128& s_a, i128& s_v, double& mx_a, double& mx_v, Long& bad) {
    const int NN = Nm << l;
    const double xdg = x_dm_m;
    const rpk::MeshGeom g = l == 0 ? mg : mesh_geom(l);
    const FFab& Fr = FRC(l);
    for (ax::MFIter mfi(Fr); mfi.isValid(); ++mfi) {
      const ax::Box gbx = ax::grow(mfi.validbox(), 1);
      const auto glo = ax::lbound(gbx), ghi = ax::ubound(gbx);
      auto fa = Fr.const_array(mfi);
      auto& pt = tile(l, mfi);
      const int np = static_cast<int>(pt.numParticles());
      const PTD pd = pt.getParticleTileData();
      const int nt = nthreads_max();
      std::vector<i128> ta(nt, 0), tv(nt, 0);
      std::vector<double> tma(nt, 0.0), tmv(nt, 0.0);
#pragma omp parallel reduction(+ : bad)
      {
        const int t = thread_id();
        i128 la = 0, lv = 0;
        double lma = 0, lmv = 0;
#pragma omp for schedule(static)
        for (int n = 0; n < np; ++n) {
          const float x[3] = {pd.m_rdata[kX][n], pd.m_rdata[kX + 1][n], pd.m_rdata[kX + 2][n]};
          const int h[3] = {home_l(pd.m_rdata[0][n], l), home_l(pd.m_rdata[1][n], l), home_l(pd.m_rdata[2][n], l)};
          const rpk::Cic c = rpk::cic_weights(g, x);
          double fv[3];
          rpk::particle_force(
              c,
              [&](int ii, int jj, int kk, double F[3]) {
                const int i = near_image(ii - 1, h[0], NN), j = near_image(jj - 1, h[1], NN),
                          k = near_image(kk - 1, h[2], NN);
                if (i < glo.x || i > ghi.x || j < glo.y || j > ghi.y || k < glo.z || k > ghi.z) {
                  bad += 1;
                  F[0] = F[1] = F[2] = 0;
                  return;
                }
                for (int m = 0; m < 3; ++m) F[m] = fa(i, j, k, m);
              },
              fv);
          // xv_dm(4:6,ip) = xv_dm(4:6,ip) + (f/x_dm_gas)*dt_sim  (x_dm_m)
          float v[3];
          for (int m = 0; m < 3; ++m) {
            v[m] = static_cast<float>(static_cast<double>(pd.m_rdata[kV + m][n]) + (fv[m] / xdg) * dt_sim);
            pd.m_rdata[kV + m][n] = v[m];
          }
          if (pd.m_idata[kType][n] != 0) continue;  // the statistics are the dark matter's
          // arms = arms + sum(f**2); amax = max(amax,maxval(abs(f)))
          const double a2 = (0.0 + fv[0] * fv[0]) + fv[1] * fv[1] + fv[2] * fv[2];
          la += to_fx128(a2, kStatBits);
          lma = f::max(lma, f::max(f::max(std::fabs(fv[0]), std::fabs(fv[1])), std::fabs(fv[2])));
          // vrms = vrms + sum(xv_dm(4:6,ip)**2)   (real(4) sum)
          const float vs = ((0.0f + v[0] * v[0]) + v[1] * v[1]) + v[2] * v[2];
          lv += to_fx128(static_cast<double>(vs), kStatBits);
          const float va = f::max(f::max(std::fabs(v[0]), std::fabs(v[1])), std::fabs(v[2]));
          lmv = f::max(lmv, static_cast<double>(va));
        }
        ta[t] = la;
        tv[t] = lv;
        tma[t] = lma;
        tmv[t] = lmv;
      }
      for (int t = 0; t < nt; ++t) {
        s_a += ta[t];
        s_v += tv[t];
        mx_a = f::max(mx_a, tma[t]);
        mx_v = f::max(mx_v, tmv[t]);
      }
    }
  }

  void kick_dm(double dt_sim, double& amax, double& arms, double& vmax, double& vrms) {
    const double t0 = now();
    i128 s_a = 0, s_v = 0;
    double mx_a = 0, mx_v = 0;
    Long bad = 0;
    for (int l = 0; l <= finest(); ++l) kick_dm_level(l, dt_sim, s_a, s_v, mx_a, mx_v, bad);
    if (bad) throw std::runtime_error("pm_backend = amrex: a CIC corner of the kick fell outside the ghost layer");
    i128 sums[2] = {s_a, s_v};
    allreduce_i128(sums, 2);
    ax::ParallelDescriptor::ReduceRealMax(mx_a);
    ax::ParallelDescriptor::ReduceRealMax(mx_v);
    amax = mx_a;
    arms = from_fx128(sums[0], kStatBits);
    vmax = mx_v;
    vrms = from_fx128(sums[1], kStatBits);
    T.kick_dm += now() - t0;
  }

  // ---- calc_kick: accelerate_gas (the three variants of src/pm.cpp) -------

  // The force on a gas cell in gas grid units: pm_refine = 1 the mesh
  // force; otherwise the average of the mesh forces over the r^3 mesh cells
  // (k, j, i order, double sums), divided by r (gas_force in src/pm.cpp)
  template <class FA>
  void gas_force(const FA& fa, int i, int j, int k, double F[3]) const {
    if (r == 1) {
      for (int m = 0; m < 3; ++m) F[m] = fa(i, j, k, m);
      return;
    }
    double sum[3] = {0, 0, 0};
    for (int kk = k * r; kk < (k + 1) * r; ++kk)
      for (int jj = j * r; jj < (j + 1) * r; ++jj)
        for (int ii = i * r; ii < (i + 1) * r; ++ii)
          for (int m = 0; m < 3; ++m) sum[m] += static_cast<double>(fa(ii, jj, kk, m));
    const double norm = static_cast<double>(r) * r * r * r;
    for (int m = 0; m < 3; ++m) F[m] = sum[m] / norm;
  }

  // kind: 0 MF (v_gas += F dt), 1 non-MF (u_gas(2:5)), 2 unsplit
  void kick_gas(int kind, double dt_sim, const Units& u, double eta1) {
    const double t0 = now();
    for (ax::MFIter mfi(gas); mfi.isValid(); ++mfi) {
      auto g = gas.array(mfi);
      auto fa = frc.const_array(mfi);
      ax::ParallelForOMP(mfi.validbox(), [=, this](int i, int j, int k) {
        double F[3];
        gas_force(fa, i, j, k, F);
        if (kind == 0) {
          // v_gas(:,i,j,k) = v_gas(:,i,j,k) + force(i,j,k)*dt_sim
          for (int m = 0; m < 3; ++m)
            g(i, j, k, 5 + m) = static_cast<float>(static_cast<double>(g(i, j, k, 5 + m)) + F[m] * dt_sim);
        } else if (kind == 1) {
          float ug[5];
          for (int m = 0; m < 5; ++m) ug[m] = g(i, j, k, m);
          const double rho = f::max(static_cast<double>(ug[0]), u.Dgas_min);
          double dv[3];
          for (int m = 0; m < 3; ++m) dv[m] = static_cast<double>(ug[1 + m]) / rho;
          double eth = static_cast<double>(ug[4]) -
                       rho * ((0.0 + dv[0] * dv[0]) + dv[1] * dv[1] + dv[2] * dv[2]) / 2;
          eth = f::min(f::max(eth, rho * u.Tgas_min), rho * u.Tgas_max);
          ug[0] = static_cast<float>(rho);
          for (int m = 0; m < 3; ++m) ug[1 + m] = static_cast<float>(rho * (dv[m] + F[m] * dt_sim));
          const float s2 = ((0.0f + ug[1] * ug[1]) + ug[2] * ug[2]) + ug[3] * ug[3];
          ug[4] = static_cast<float>(eth + static_cast<double>(s2) / rho / 2);
          for (int m = 0; m < 5; ++m) g(i, j, k, m) = ug[m];
        } else {
          double U[5];
          for (int m = 0; m < 5; ++m) U[m] = g(i, j, k, m);
          const double eD = g(i, j, k, 5);
          const double eth = rrh::thermal_energy(U, eD, eta1);
          const double rho = f::max(U[0], u.Dgas_min);
          for (int m = 0; m < 3; ++m) U[1 + m] = U[1 + m] + rho * F[m] * dt_sim;
          U[4] = rrh::kinetic_energy(U) + eth;
          for (int m = 1; m < 5; ++m) g(i, j, k, m) = static_cast<float>(U[m]);
        }
      });
    }
    T.kick_gas += now() - t0;
  }

  // ---- level -> host state on every rank ------------------------------------

  // The kicked velocities into s.xv, the kicked gas components into
  // u_gas / v_gas, the potential into s.density when asked: each rank
  // writes its own values into zeroed arrays, then a bitwise OR all-reduce
  void gather(pm::PMState& s, float* u_gas, float* v_gas, int kind) {
    const double t0 = now();
    if (state) {
      // amrex_state = 1: the particles stay on the level (their velocities
      // are kicked in place); the kicked u(2:5) into the level's gas state
      for (ax::MFIter mfi(gas); mfi.isValid(); ++mfi) {
        auto g = gas.const_array(mfi);
        auto ga = state->G.array(mfi);
        ax::ParallelForOMP(mfi.validbox(), [=](int i, int j, int k) {
          for (int m = 1; m < 5; ++m) ga(i, j, k, gs::kU + m) = g(i, j, k, m);
        });
      }
      T.g_gas += now() - t0;
      T.gather += now() - t0;
      return;
    }
    // DM velocities
    std::vector<float> vb(static_cast<std::size_t>(3) * Npart, 0.0f);
    for (ax::MFIter mfi(pot); mfi.isValid(); ++mfi) {
      auto& pt = pc->DefineAndReturnParticleTile(0, mfi.index(), mfi.LocalTileIndex());
      const int np = static_cast<int>(pt.numParticles());
      const PTD pd = pt.getParticleTileData();
#pragma omp parallel for schedule(static)
      for (int n = 0; n < np; ++n) {
        const std::int64_t ip = static_cast<std::int64_t>(ax::ParticleIDWrapper(pd.m_idcpu[n]));
        float* v = &vb[static_cast<std::size_t>(ip - 1) * 3];
        for (int m = 0; m < 3; ++m) v[m] = pd.m_rdata[kV + m][n];
      }
    }
    const double t1 = now();
    T.g_vwrite += t1 - t0;
    allreduce_bor(vb.data(), vb.size());
    const double t2 = now();
    T.g_vred += t2 - t1;
#pragma omp parallel for schedule(static)
    for (std::int64_t ip = 1; ip <= Npart; ++ip) {
      float* p = s.xv_p(ip);
      const float* v = &vb[static_cast<std::size_t>(ip - 1) * 3];
      for (int m = 0; m < 3; ++m) p[3 + m] = v[m];
    }
    vb.clear();
    vb.shrink_to_fit();
    const double t3 = now();
    T.g_vcopy += t3 - t2;
    // the gas: the changed components (MF: v; non-MF: u; unsplit: u(2:5))
    const std::int64_t ncell = static_cast<std::int64_t>(Ng) * Ng * Ng;
    const int NG = Ng;
    float* arr = kind == 0 ? v_gas : u_gas;
    const int nc = kind == 0 ? 3 : 5;
    const int c0 = kind == 2 ? 1 : 0;
#pragma omp parallel for schedule(static)
    for (std::int64_t c = 0; c < ncell; ++c)
      for (int m = c0; m < nc; ++m) arr[c * nc + m] = 0.0f;
    for (ax::MFIter mfi(gas); mfi.isValid(); ++mfi) {
      auto g = gas.const_array(mfi);
      ax::ParallelForOMP(mfi.validbox(), [=](int i, int j, int k) {
        const std::int64_t gc = i + static_cast<std::int64_t>(NG) * (j + static_cast<std::int64_t>(NG) * k);
        for (int m = c0; m < nc; ++m) arr[gc * nc + m] = g(i, j, k, kind == 0 ? 5 + m : m);
      });
    }
    allreduce_bor(arr, static_cast<std::size_t>(ncell) * nc);
    const double t4 = now();
    T.g_gas += t4 - t3;
    // the potential (pm_driver's phi output; the host's density array keeps
    // the in-place FFT padding)
    if (s.gather_potential) {
      std::fill(s.density.begin(), s.density.end(), 0.0f);
      float* d = s.density.data();
      const int NM = Nm;
      for (ax::MFIter mfi(pot); mfi.isValid(); ++mfi) {
        auto p = pot.const_array(mfi);
        ax::ParallelForOMP(mfi.validbox(), [=](int i, int j, int k) {
          d[i + static_cast<std::int64_t>(NM + 2) * (j + static_cast<std::int64_t>(NM) * k)] = p(i, j, k);
        });
      }
      allreduce_bor(d, s.density.size());
    }
    T.g_pot += now() - t4;
    T.gather += now() - t0;
  }

  // ---- the fields for the tests: potential and mesh force onto the host -------

  void fields(pm::PMState& s, const float* u_gas, double scalefactor) {
    const double mass = (1 - cp::omegab / cp::omegam) * s.pc.dm_mass_ratio;
    if (state) {
      sync_levels();
      pc->Redistribute();
    }
    prep(s, u_gas, nullptr, nullptr, state ? 2 : (mf ? 0 : 1));
    density(mass);
    gravity(scalefactor);
    force();
    const bool keep = s.gather_potential;
    s.gather_potential = true;
    // the potential (gather's last part) and the force mesh
    {
      std::fill(s.density.begin(), s.density.end(), 0.0f);
      float* d = s.density.data();
      const int NM = Nm;
      for (ax::MFIter mfi(pot); mfi.isValid(); ++mfi) {
        auto p = pot.const_array(mfi);
        ax::ParallelForOMP(mfi.validbox(), [=](int i, int j, int k) {
          d[i + static_cast<std::int64_t>(NM + 2) * (j + static_cast<std::int64_t>(NM) * k)] = p(i, j, k);
        });
      }
      allreduce_bor(d, s.density.size());
      s.fmesh.assign(3 * static_cast<std::size_t>(s.pc.Nmcell), 0.0f);
      float* fm = s.fmesh.data();
      for (ax::MFIter mfi(frc); mfi.isValid(); ++mfi) {
        auto f = frc.const_array(mfi);
        ax::ParallelForOMP(mfi.validbox(), [=](int i, int j, int k) {
          float* o = fm + 3 * (i + static_cast<std::int64_t>(NM) * (j + static_cast<std::int64_t>(NM) * k));
          for (int m = 0; m < 3; ++m) o[m] = f(i, j, k, m);
        });
      }
      allreduce_bor(fm, s.fmesh.size());
      s.fmesh_valid = true;
    }
    s.gather_potential = keep;
  }

  // ---- the force at given positions (tests) -----------------------------------

  // The force a particle at each position of xs (DM-grid units) would get
  // from the finest level covering its home cell: the CIC interpolation of
  // that level's mesh force, in level-0 mesh units. Collective.
  void force_at(const std::vector<float>& xs, std::vector<double>& out, std::vector<int>* levels) {
    const int n = static_cast<int>(xs.size() / 3);
    out.assign(static_cast<std::size_t>(3) * n, 0.0);
    std::vector<int> lv(static_cast<std::size_t>(n), 0);
    const int myproc = ax::ParallelDescriptor::MyProc();
    for (int q = 0; q < n; ++q) {
      const float* x = &xs[static_cast<std::size_t>(3) * q];
      for (int l = finest(); l >= 0; --l) {
        const int NN = Nm << l;
        const double xdg = x_dm_m * static_cast<double>(1 << l);
        const int h[3] = {f::mod(f::int4(static_cast<double>(x[0]) * xdg), NN),
                          f::mod(f::int4(static_cast<double>(x[1]) * xdg), NN),
                          f::mod(f::int4(static_cast<double>(x[2]) * xdg), NN)};
        const auto gt = LOC(l).getGridAssignor()(ax::IntVect(h[0], h[1], h[2]), ax::IntVect(0));
        if (gt.first < 0) continue;
        if (DM(l)[gt.first] == myproc) {
          const FFab& Fr = FRC(l);
          const ax::Box gbx = ax::grow(BA(l)[gt.first], 1);
          const auto glo = ax::lbound(gbx), ghi = ax::ubound(gbx);
          auto fa = Fr.const_array(gt.first);
          const rpk::Cic c = rpk::cic_weights(l == 0 ? mg : mesh_geom(l), x);
          double fv[3];
          bool bad = false;
          rpk::particle_force(
              c,
              [&](int ii, int jj, int kk, double F[3]) {
                const int i = near_image(ii - 1, h[0], NN), j = near_image(jj - 1, h[1], NN),
                          k = near_image(kk - 1, h[2], NN);
                if (i < glo.x || i > ghi.x || j < glo.y || j > ghi.y || k < glo.z || k > ghi.z) {
                  bad = true;
                  F[0] = F[1] = F[2] = 0;
                  return;
                }
                for (int m = 0; m < 3; ++m) F[m] = fa(i, j, k, m);
              },
              fv);
          if (bad) throw std::runtime_error("amrex_pm_force_at: a CIC corner fell outside the ghost layer");
          for (int m = 0; m < 3; ++m) out[static_cast<std::size_t>(3) * q + m] = fv[m];
          lv[static_cast<std::size_t>(q)] = l;
        }
        break;
      }
    }
    ax::ParallelDescriptor::ReduceRealSum(out.data(), 3 * n);
    ax::ParallelDescriptor::ReduceIntMax(lv.data(), n);
    if (levels) *levels = lv;
  }

  // ---- one PM step ------------------------------------------------------------

  // kind: 0 MF, 1 non-MF, 2 unsplit (e_gas may be null: dual_energy = 0)
  void step(pm::PMState& s, float* u_gas, float* v_gas, const float* e_gas, int kind, double eta1,
            SimClock& clk, const Units& u) {
    const double mass = (1 - cp::omegab / cp::omegam) * s.pc.dm_mass_ratio;  // dm_mass of src/pm.cpp
    const bool pcheck = state && std::getenv("RADRAY_PM_MOMENTUM_CHECK") != nullptr;
    if (sync_levels()) pc->Redistribute();  // (a no-op after the drift's; the tests call the step without a drift)
    double P0[6];
    FFab before;
    if (pcheck) {
      momentum(P0, mass);
      before.define(state->ba, state->dm, 3, 0);
      before.ParallelCopy(state->G, gs::kU + 1, 0, 3);
    }
    prep(s, u_gas, v_gas, e_gas, kind);
    density(mass);
    gravity(clk.scalefactor);
    force();
    double amax = 0, arms = 0, vmax = 0, vrms = 0;
    kick_dm(clk.dt_sim, amax, arms, vmax, vrms);
    kick_gas(kind, clk.dt_sim, u, eta1);
    gather(s, u_gas, v_gas, kind);
    if (fine_levels() > 0) kick_fine(clk.dt_sim, u, eta1);
    if (finest() > 0) kick_levels(clk.dt_sim, u, eta1);
    if (pcheck) {
      double P1[6];
      momentum(P1, mass);
      const double sabs = gas_change_abs(before);
      double rel = 0;
      for (int m = 0; m < 3; ++m) {
        const double d = (P1[m] - P0[m]) + (P1[3 + m] - P0[3 + m]);
        rel = std::max(rel, sabs > 0 ? std::fabs(d) / sabs : 0.0);
      }
      if (ax::ParallelDescriptor::IOProcessor())
        std::fprintf(stderr, "PM momentum check: levels %d, |dP_dm + dP_gas| / sum|dp_gas| = %.3e (sum|dp_gas| %.6e)\n",
                     std::max(fine_levels(), finest()), rel, sabs);
    }
    // set_kick_globals (PM grid units)
    const double npart = static_cast<double>(s.pc.Npart);
    clk.amax_dm = amax;
    clk.arms_dm = std::sqrt(arms / npart);
    clk.vmax_dm = vmax;
    clk.vrms_dm = std::sqrt(vrms / npart);
    T.steps += 1;
  }

  std::string timers() const {
    char buf[768];
    std::snprintf(buf, sizeof buf,
                  "AMReX PM timers [s]: prep %.3f (locate %.3f fill %.3f) drift %.3f deposit %.3f (SumBoundary %.3f) "
                  "fft %.3f force %.3f kick_dm %.3f kick_gas %.3f gather %.3f (v write %.3f reduce %.3f "
                  "copy %.3f, gas %.3f, potential %.3f); steps %lld",
                  T.prep, T.p_locate, T.p_fill, T.drift, T.deposit, T.sumbnd, T.fft, T.force, T.kick_dm, T.kick_gas,
                  T.gather, T.g_vwrite, T.g_vred, T.g_vcopy, T.g_gas, T.g_pot, static_cast<long long>(T.steps));
    if (gas_levels) {
      char b2[128];
      std::snprintf(b2, sizeof b2, "; gas levels: deposit %.3f kick %.3f", T.d_fine, T.k_fine);
      return std::string(buf) + b2;
    }
    if (level_solve) {
      char b2[256];
      std::snprintf(b2, sizeof b2, "; level solve: sync %.3f deposit %.3f solve %.3f force %.3f kick %.3f; levels %d",
                    T.sync, T.d_lev, T.solve, T.f_lev, T.k_lev, finest());
      std::string out = std::string(buf) + b2;
      for (const FineMesh& M : fm) {
        char b3[128];
        std::snprintf(b3, sizeof b3,  " [%d: %lld cells, mg depth %d, bottom %d, sweeps %d, %d iters, resid %.2e%s]", M.lev,
                      static_cast<long long>(M.ba.numPts()), M.mg_depth, M.bottom, M.sweeps, M.iters, M.resid,
                      M.capped ? (", capped " + std::to_string(M.capped) + "x").c_str() : "");
        out += b3;
      }
      return out;
    }
    return buf;
  }
};

}  // namespace

struct AmrexPM {
  std::unique_ptr<PMLevel> level;
};

namespace {
PMLevel& level_of(pm::PMState& s) {
  if (!s.amrex) {
    if (!amrex_available() || backend_nranks() < 1)
      throw std::runtime_error("pm_backend = amrex: backend_init was not called");
    s.amrex = std::make_shared<AmrexPM>();
    s.amrex->level = std::make_unique<PMLevel>(s.pc);
  }
  PMLevel& L = *s.amrex->level;
  if (L.Nm != s.pc.Nm || L.Npart != s.pc.Npart)
    throw std::runtime_error("pm_backend = amrex: the PM configuration changed after the level was built");
  return L;
}
}  // namespace

void amrex_pm_step(pm::PMState& s, float* u_gas, float* v_gas, SimClock& clk, const Units& u) {
  level_of(s).step(s, u_gas, v_gas, nullptr, s.pc.mf ? 0 : 1, 0.0, clk, u);
}

void amrex_pm_step_unsplit(pm::PMState& s, float* u_gas, const float* e_gas, double eta1, SimClock& clk,
                           const Units& u) {
  // the dual-energy field may be absent (dual_energy = 0): the kernel then
  // takes E - KE
  level_of(s).step(s, u_gas, nullptr, e_gas, 2, eta1, clk, u);
}

void amrex_pm_fields(pm::PMState& s, const float* u_gas, double scalefactor) {
  level_of(s).fields(s, u_gas, scalefactor);
}
void amrex_pm_force_at(pm::PMState& s, const std::vector<float>& xs, std::vector<double>& F,
                       std::vector<int>* levels) {
  level_of(s).force_at(xs, F, levels);
}

void amrex_pm_load_particles(pm::PMState& s) { level_of(s).load_particles(s); }
std::int64_t amrex_pm_add_stars(pm::PMState& s, const std::vector<SinkEvent>& ev, double a_form) {
  return level_of(s).add_stars(s, ev, a_form);
}
std::int64_t amrex_pm_gather_stars(pm::PMState& s) { return level_of(s).gather_stars(s); }
void amrex_pm_load_stars(pm::PMState& s) { level_of(s).load_stars(s); }
void pm_remap(pm::PMState& s, const amrex::DistributionMapping& dm) { level_of(s).remap(dm); }
void amrex_pm_gather_particles(pm::PMState& s) { level_of(s).gather_particles(s); }
void amrex_pm_drift(pm::PMState& s, double dt_sim) { level_of(s).drift(dt_sim, s.pc.ndmr); }

std::string amrex_pm_timers(const pm::PMState& s) { return s.amrex ? s.amrex->level->timers() : ""; }

}  // namespace radray
