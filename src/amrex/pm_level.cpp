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
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
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
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParticleLocator.H>
#include <AMReX_Particles.H>
#include <fftw3.h>
#ifdef AMREX_USE_OMP
#include <omp.h>
#endif

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
using PC = ax::ParticleContainerPureSoA<9, 0>;  // 3 positions + x(3) + v(3)
using PTD = PC::ParticleTileType::ParticleTileDataType;
constexpr int kX = 3, kV = 6;  // SoA real components

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
using i128 = __int128;

inline Long to_fx(double x, int bits) { return static_cast<Long>(std::floor(std::ldexp(x, bits) + 0.5)); }
inline i128 to_fx128(double x, int bits) { return static_cast<i128>(std::floor(std::ldexp(x, bits) + 0.5)); }
inline double from_fx128(i128 v, int bits) { return std::ldexp(static_cast<double>(v), -bits); }

// MPI sum of 128-bit integers as four 32-bit limbs in int64 (rt_level.cpp)
void allreduce_i128(i128* v, int n) {
  if (ax::ParallelDescriptor::NProcs() == 1) return;
  std::vector<Long> limbs(static_cast<std::size_t>(4 * n));
  for (int i = 0; i < n; ++i) {
    const unsigned __int128 u = static_cast<unsigned __int128>(v[i]);
    for (int k = 0; k < 4; ++k) limbs[4 * i + k] = static_cast<Long>((u >> (32 * k)) & 0xffffffffULL);
  }
  ax::ParallelDescriptor::ReduceLongSum(limbs.data(), 4 * n);
  for (int i = 0; i < n; ++i) {
    unsigned __int128 r = 0;
    for (int k = 0; k < 4; ++k)
      r += static_cast<unsigned __int128>(static_cast<unsigned long long>(limbs[4 * i + k])) << (32 * k);
    v[i] = static_cast<i128>(r);
  }
}

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
// The level
// ---------------------------------------------------------------------------
struct PMLevel {
  int Nm = 0, Ng = 0, r = 1, mgs = 16;
  std::int64_t Npart = 0;
  bool mf = true;
  int fft_kind = 0;
  rpk::MeshGeom mg;
  double x_dm_m = 1;
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
    }
    pot.define(ba, dm, 1, 3);
    frc.define(ba, dm, 3, 1);
    acc.define(ba, dm, 1, 1);
    gas.define(bag, dm, 8, 0);
    pc = std::make_unique<PC>(geom, dm, ba);
    for (ax::MFIter mfi(pot); mfi.isValid(); ++mfi)
      pc->DefineAndReturnParticleTile(0, mfi.index(), mfi.LocalTileIndex());
    locator.build(ba, geom);
  }

  // ---- host state -> level ------------------------------------------------

  // The particles of s.xv whose home cell (link_dm's cell) lies in a local
  // box, and the gas of the local gas boxes (e_gas: the dual-energy field in
  // place of v_gas when the kick is the unsplit one)
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
          pd.m_rdata[0][n_] = static_cast<float>(home(p[0])) + 0.5f;
          pd.m_rdata[1][n_] = static_cast<float>(home(p[1])) + 0.5f;
          pd.m_rdata[2][n_] = static_cast<float>(home(p[2])) + 0.5f;
          for (int c = 0; c < 3; ++c) {
            pd.m_rdata[kX + c][n_] = p[c];
            pd.m_rdata[kV + c][n_] = p[3 + c];
          }
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
    const double t0 = now();
    const double xdg = x_dm_m;
    const int NN = Nm;
    for (ax::MFIter mfi(pot); mfi.isValid(); ++mfi) {
      auto& pt = pc->DefineAndReturnParticleTile(0, mfi.index(), mfi.LocalTileIndex());
      const int np = static_cast<int>(pt.numParticles());
      const PTD pd = pt.getParticleTileData();
#pragma omp parallel for schedule(static)
      for (int n = 0; n < np; ++n) {
        for (int c = 0; c < 3; ++c) {
          const double t = static_cast<double>(pd.m_rdata[kX + c][n]) +
                           static_cast<double>(pd.m_rdata[kV + c][n]) * dt_sim / 2 + ndmr;
          const float x = static_cast<float>(f::mod(t, ndmr));
          pd.m_rdata[kX + c][n] = x;
          pd.m_rdata[c][n] = static_cast<float>(f::mod(f::int4(static_cast<double>(x) * xdg), NN)) + 0.5f;
        }
      }
    }
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
    for (ax::MFIter mfi(pot); mfi.isValid(); ++mfi) {
      auto& pt = pc->DefineAndReturnParticleTile(0, mfi.index(), mfi.LocalTileIndex());
      const int np = static_cast<int>(pt.numParticles());
      const PTD pd = pt.getParticleTileData();
      const std::size_t n0 = mine.size();
      mine.resize(n0 + np);
      for (int n = 0; n < np; ++n) {
        Rec& r = mine[n0 + n];
        r.id = static_cast<Long>(ax::ParticleIDWrapper(pd.m_idcpu[n]));
        for (int c = 0; c < 3; ++c) {
          r.xv[c] = pd.m_rdata[kX + c][n];
          r.xv[3 + c] = pd.m_rdata[kV + c][n];
        }
      }
    }
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

  void prep(const pm::PMState& s, const float* u_gas, const float* v_gas, const float* e_gas, int kind) {
    const double t0 = now();
    const int myproc = ax::ParallelDescriptor::MyProc();
    const int nlocal = pot.local_size();
    if (state) {
      // amrex_state = 1: the particles are on the level already (load_particles,
      // drift); only the gas below
      if (!loaded) throw std::runtime_error("amrex_state = 1: amrex_pm_load_particles was not called");
      if (kind != 2) throw std::runtime_error("amrex_state = 1 needs hydro_scheme = unsplit");
      for (ax::MFIter mfi(gas); mfi.isValid(); ++mfi) {
        auto g = gas.array(mfi);
        auto ga = state->G.const_array(mfi);
        ax::ParallelForOMP(mfi.validbox(), [=](int i, int j, int k) {
          for (int m = 0; m < 5; ++m) g(i, j, k, m) = ga(i, j, k, gs::kU + m);
          g(i, j, k, 5) = ga(i, j, k, gs::kE);
          g(i, j, k, 6) = g(i, j, k, 7) = 0.0f;
        });
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

  void density(double mass) {
    const double t0 = now();
    acc.setVal(0);
    const int NN = Nm;
    const rpk::MeshGeom g = mg;
    Long bad = 0;
    for (ax::MFIter mfi(acc); mfi.isValid(); ++mfi) {
      const ax::Box bx = mfi.validbox();
      const ax::Box gbx = ax::grow(bx, 1);
      const auto glo = ax::lbound(gbx), ghi = ax::ubound(gbx), glen = ax::length(gbx);
      const Long ngc = gbx.numPts();
      auto& pt = pc->DefineAndReturnParticleTile(0, mfi.index(), mfi.LocalTileIndex());
      const int np = static_cast<int>(pt.numParticles());
      const PTD pd = pt.getParticleTileData();
      const int nt = nthreads_max();
      std::vector<Long> buf(static_cast<std::size_t>(nt) * ngc, 0);
      auto a = acc.array(mfi);
#pragma omp parallel reduction(+ : bad)
      {
        Long* mine = buf.data() + static_cast<std::size_t>(thread_id()) * ngc;
#pragma omp for schedule(static)
        for (int n = 0; n < np; ++n) {
          const float x[3] = {pd.m_rdata[kX][n], pd.m_rdata[kX + 1][n], pd.m_rdata[kX + 2][n]};
          const int h[3] = {static_cast<int>(pd.m_rdata[0][n]), static_cast<int>(pd.m_rdata[1][n]),
                            static_cast<int>(pd.m_rdata[2][n])};
          const rpk::Cic c = rpk::cic_weights(g, x);
          rpk::cic_corners(c, [&](int i, int j, int k, double w) {
            const int ii = near_image(i - 1, h[0], NN), jj = near_image(j - 1, h[1], NN),
                      kk = near_image(k - 1, h[2], NN);
            if (ii < glo.x || ii > ghi.x || jj < glo.y || jj > ghi.y || kk < glo.z || kk > ghi.z) {
              bad += 1;
              return;
            }
            mine[(ii - glo.x) + glen.x * ((jj - glo.y) + static_cast<Long>(glen.y) * (kk - glo.z))] +=
                to_fx(w * mass, kDensBits);
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
    if (bad) throw std::runtime_error("pm_backend = amrex: a CIC corner fell outside the ghost layer");
    const double t1 = now();
    acc.SumBoundary(period);  // fold the ghost deposits into their owners (int64: exact)
    T.sumbnd += now() - t1;
    // density = the gas density of the gas cell the mesh cell lies in + the
    // DM mass (the original: float, added particle by particle; here one
    // rounding of the exact sum)
    const int rr = r;
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
      // amrex::FFT sets FFTW's planner to its thread count: restore the
      // serial plans of the host code (hydro, the serial PM)
      fftwf_plan_with_nthreads(1);
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
  }

  // ---- the mesh force on the box and its ghost layer ------------------------

  void force() {
    const double t0 = now();
    const int NN = Nm;
    for (ax::MFIter mfi(frc); mfi.isValid(); ++mfi) {
      auto f = frc.array(mfi);
      auto p = pot.const_array(mfi);
      ax::ParallelForOMP(ax::grow(mfi.validbox(), 1), [=](int i, int j, int k) {
        double F[3];
        rpk::force_stencil(
            NN, i + 1, j + 1, k + 1,
            [&](int ii, int jj, int kk) -> float {
              return p(near_image(ii - 1, i, NN), near_image(jj - 1, j, NN), near_image(kk - 1, k, NN));
            },
            F);
        for (int m = 0; m < 3; ++m) f(i, j, k, m) = static_cast<float>(F[m]);  // exact: F is a float
      });
    }
    T.force += now() - t0;
  }

  // ---- calc_kick: accelerate_dm ---------------------------------------------

  void kick_dm(double dt_sim, double& amax, double& arms, double& vmax, double& vrms) {
    const double t0 = now();
    const int NN = Nm;
    const double xdg = x_dm_m;
    const rpk::MeshGeom g = mg;
    i128 s_a = 0, s_v = 0;
    double mx_a = 0, mx_v = 0;
    Long bad = 0;
    for (ax::MFIter mfi(frc); mfi.isValid(); ++mfi) {
      const ax::Box gbx = ax::grow(mfi.validbox(), 1);
      const auto glo = ax::lbound(gbx), ghi = ax::ubound(gbx);
      auto fa = frc.const_array(mfi);
      auto& pt = pc->DefineAndReturnParticleTile(0, mfi.index(), mfi.LocalTileIndex());
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
          const int h[3] = {static_cast<int>(pd.m_rdata[0][n]), static_cast<int>(pd.m_rdata[1][n]),
                            static_cast<int>(pd.m_rdata[2][n])};
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
    const double mass = (1 - cp::omegab / cp::omegam) * s.pc.N_dm_m_c;
    prep(s, u_gas, nullptr, nullptr, mf ? 0 : 1);
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

  // ---- one PM step ------------------------------------------------------------

  // kind: 0 MF, 1 non-MF, 2 unsplit (e_gas may be null: dual_energy = 0)
  void step(pm::PMState& s, float* u_gas, float* v_gas, const float* e_gas, int kind, double eta1,
            SimClock& clk, const Units& u) {
    const double mass = (1 - cp::omegab / cp::omegam) * s.pc.N_dm_m_c;  // dm_mass of src/pm.cpp
    prep(s, u_gas, v_gas, e_gas, kind);
    density(mass);
    gravity(clk.scalefactor);
    force();
    double amax = 0, arms = 0, vmax = 0, vrms = 0;
    kick_dm(clk.dt_sim, amax, arms, vmax, vrms);
    kick_gas(kind, clk.dt_sim, u, eta1);
    gather(s, u_gas, v_gas, kind);
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

void amrex_pm_load_particles(pm::PMState& s) { level_of(s).load_particles(s); }
void amrex_pm_gather_particles(pm::PMState& s) { level_of(s).gather_particles(s); }
void amrex_pm_drift(pm::PMState& s, double dt_sim) { level_of(s).drift(dt_sim, s.pc.ndmr); }

std::string amrex_pm_timers(const pm::PMState& s) { return s.amrex ? s.amrex->level->timers() : ""; }

}  // namespace radray
