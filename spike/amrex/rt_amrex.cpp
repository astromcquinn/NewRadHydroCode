// AMReX prototype of the sorted-ray-array RT substep: the same phases, in
// the same order, as spike/serial/rt_serial.cpp, with the same kernels
// (spike/common) and the same host logic (rrhost::RTRun), but distributed:
//
//   * cell fields   FabArray<BaseFab<float>> C (all 45 floats of
//                   rrk::CellData, no ghosts) and TF (tau(5), fgam(5),
//                   4 ghosts, FillBoundary before cast_1 / cast_3);
//   * deposits      FabArray<BaseFab<Long>> A (13 int64 fixed-point
//                   accumulators, 4 ghosts), atomic adds into the box
//                   including its ghosts, then SumBoundary (periodic);
//   * rays          ParticleContainerPureSoA (position = centre of the home
//                   cell; the RayP payload in 18 float + 2 int components);
//                   gathered per box into a cell-sorted rrk::RayP scratch
//                   array (DenseBins), migrated by Redistribute() after
//                   cast_3;
//   * reductions    integer histogram / stats / ngam_tot via ReduceLongSum.
//
// Bitwise identical to rt_serial for any rank count, max_grid_size and
// thread count (all cross-cell sums are integer).
//
// usage: rt_amrex <problem_dir> <outdir> [--nsteps n] [--set key=value]...
//                 [--no-output] [--max_grid_size M] [--redist_local L] [--ghosts G]
#include <chrono>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <type_traits>
#include <vector>

#include <AMReX.H>
#include <AMReX_DenseBins.H>
#include <AMReX_FabArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuAtomic.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Particles.H>
#include <AMReX_Reduce.H>
#include <AMReX_Scan.H>

#include "rt_host.hpp"

using namespace amrex;
using RTCell = rrk::CellData;  // not CellData: amrex::CellData exists too
using rrk::RayP;

namespace {

// ---- cell layout: RTCell as 45 float components ------------------------
constexpr int NC = static_cast<int>(sizeof(RTCell) / sizeof(float));
static_assert(sizeof(RTCell) == 45 * sizeof(float), "RTCell must be 45 packed floats");
constexpr int kLuv = offsetof(RTCell, Luv) / sizeof(float);
constexpr int kSfr = offsetof(RTCell, sfr) / sizeof(float);
constexpr int kTau = offsetof(RTCell, tau) / sizeof(float);
constexpr int kNGdefault = 4;  // ghost cells of TF and A (--ghosts)

using FFab = FabArray<BaseFab<float>>;
using LFab = FabArray<BaseFab<Long>>;
static_assert(sizeof(Long) == sizeof(long long), "int64 accumulators");

// Cell storage. Default: AoS, one rrk::CellData per cell in a
// FabArray<BaseFab<rrk::CellData>> (1 component); the per-cell kernels take
// the whole struct, and on the CPU this is 3-5x faster than the SoA layout
// (45 float components, 45 memory streams per cell). -DRT_CELL_SOA selects
// the SoA layout (FabArray<BaseFab<float>>, 45 components), the natural one
// for coalesced GPU access.
#ifdef RT_CELL_SOA
using CFab = FFab;
constexpr int kCComp = NC;
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE RTCell load_cell(Array4<float const> const& a, int i, int j, int k) {
  const auto cd = a.cellData(i, j, k);  // amrex::CellData: pointer + component stride
  float b[NC];
  for (int m = 0; m < NC; ++m) b[m] = cd[m];
  RTCell d;
  std::memcpy(&d, b, sizeof d);
  return d;
}
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void store_cell(Array4<float> const& a, int i, int j, int k,
                                                         const RTCell& d) {
  const auto cd = a.cellData(i, j, k);
  float b[NC];
  std::memcpy(b, &d, sizeof d);
  for (int m = 0; m < NC; ++m) cd[m] = b[m];
}
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE float cf(Array4<float const> const& a, int i, int j, int k, int m) {
  return a(i, j, k, m);
}
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void cset(Array4<float> const& a, int i, int j, int k, int m, float v) {
  a(i, j, k, m) = v;
}
#else
using CFab = FabArray<BaseFab<RTCell>>;
constexpr int kCComp = 1;
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE RTCell load_cell(Array4<RTCell const> const& a, int i, int j, int k) {
  return a(i, j, k);
}
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void store_cell(Array4<RTCell> const& a, int i, int j, int k,
                                                         const RTCell& d) {
  a(i, j, k) = d;
}
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE float cf(Array4<RTCell const> const& a, int i, int j, int k, int m) {
  float v;
  std::memcpy(&v, reinterpret_cast<const char*>(&a(i, j, k)) + m * sizeof(float), sizeof v);
  return v;
}
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void cset(Array4<RTCell> const& a, int i, int j, int k, int m, float v) {
  std::memcpy(reinterpret_cast<char*>(&a(i, j, k)) + m * sizeof(float), &v, sizeof v);
}
#endif

// ---- rays: PureSoA particles ----------------------------------------------
// real comps: 0-2 position (home-cell centre), 3 r, 4 dr, 5-7 x, 8-10 dx,
// 11-15 nu, 16-20 ngam; int comps: 0 nside, 1 ipix.
constexpr int kNR = 3 + 18, kNI = 2;
constexpr int kR = 3, kDR = 4, kX = 5, kDX = 8, kNU = 11, kNGAM = 16;
static_assert(std::is_same_v<ParticleReal, float>, "build AMReX with AMReX_PARTICLES_PRECISION=SINGLE");
using RayPC = ParticleContainerPureSoA<kNR, kNI>;
using PTD = RayPC::ParticleTileType::ParticleTileDataType;

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE RayP load_ray(PTD const& p, int n) {
  RayP r;
  r.nside = p.m_idata[0][n];
  r.ipix = p.m_idata[1][n];
  r.r = p.m_rdata[kR][n];
  r.dr = p.m_rdata[kDR][n];
  for (int d = 0; d < 3; ++d) {
    r.x[d] = p.m_rdata[kX + d][n];
    r.dx[d] = p.m_rdata[kDX + d][n];
  }
  for (int q = 0; q < rrk::NF; ++q) {
    r.nu[q] = p.m_rdata[kNU + q][n];
    r.ngam[q] = p.m_rdata[kNGAM + q][n];
  }
  return r;
}
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void store_payload(PTD const& p, int n, const RayP& r) {
  p.m_idata[0][n] = r.nside;
  p.m_idata[1][n] = r.ipix;
  p.m_rdata[kR][n] = r.r;
  p.m_rdata[kDR][n] = r.dr;
  for (int d = 0; d < 3; ++d) {
    p.m_rdata[kX + d][n] = r.x[d];
    p.m_rdata[kDX + d][n] = r.dx[d];
  }
  for (int q = 0; q < rrk::NF; ++q) {
    p.m_rdata[kNU + q][n] = r.nu[q];
    p.m_rdata[kNGAM + q][n] = r.ngam[q];
  }
}
// Wrapped global 0-based index g -> index of its periodic image nearest to
// the home index h (h is a valid cell of the box).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE int near_image(int g, int h, int N) {
  int d = g - h;
  if (2 * d > N) d -= N;
  else if (2 * d < -N) d += N;
  return h + d;
}

// Full store: payload, valid id, and AMReX position = centre of the new home
// cell cell_of(x), taken as the periodic image nearest to the emitting cell
// h (so it may lie just outside the domain). Redistribute() applies the
// periodic shift itself; a pre-wrapped position would look like a jump of
// ~N cells and break the local Redistribute algorithm. All values are exact
// in float, so the owner of cell_of(x) is found exactly.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void store_ray(PTD const& p, int n, const RayP& r, int N, int cpu,
                                                        const int h[3]) {
  store_payload(p, n, r);
  for (int d = 0; d < 3; ++d)
    p.m_rdata[d][n] = static_cast<float>(near_image(rrk::cell_of(r.x[d], N), h[d], N)) + 0.5f;
  p.m_idcpu[n] = SetParticleIDandCPU(1, cpu);
}

struct Timers {
  double cell = 0, hist = 0, step1 = 0, cast1 = 0, phys = 0, cast3 = 0, sort = 0, migrate = 0, post = 0;
  double fillbnd = 0, sumbnd = 0, zero = 0;  // parts of cell+phys, cast1, step1
};
double now() {
  return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

void zero(Gpu::DeviceVector<Long>& v) {
  Long* p = v.data();
  ParallelFor(static_cast<int>(v.size()), [=] AMREX_GPU_DEVICE(int n) { p[n] = 0; });
}

// Per-box scratch, kept across substeps
struct BoxScratch {
  DenseBins<int> bins;
  Gpu::DeviceVector<int> hcell;           // local home cell of each particle
  Gpu::DeviceVector<RayP> rays, out;      // cell-sorted input / output rays
  Gpu::DeviceVector<int> idx;             // step1_write scratch
  Gpu::DeviceVector<rrk::MergeKey> key;   // step1_write scratch
  Gpu::DeviceVector<Long> cnt, ostart;    // per-cell output counts / offsets
  Gpu::DeviceVector<int> ocell;           // emitting (home) cell of each output ray
  Gpu::DeviceVector<int> keep, wpos;      // cast_3 survivors / compaction offsets
  Long nout = 0;
};

}  // namespace

int main(int argc, char** argv) {
  amrex::Initialize(argc, argv, false, MPI_COMM_WORLD, []() {
    ParmParse pp("amrex");
    pp.add("verbose", 0);
  });
  int rc = 0;
  if (argc < 3) {
    std::fprintf(stderr,
                 "usage: %s <problem_dir> <outdir> [--nsteps n] [--set k=v] [--no-output] "
                 "[--max_grid_size M] [--redist_local L] [--ghosts G]\n",
                 argv[0]);
    amrex::Finalize();
    return 2;
  }
  try {
    const std::string pdir = argv[1], odir = argv[2];
    std::vector<std::string> sets;
    int nsteps = -1, mgs = 16, redist_local = 3, ng = kNGdefault;
    bool output = true;
    for (int a = 3; a < argc; ++a) {
      const std::string s = argv[a];
      if (s == "--set" && a + 1 < argc) sets.emplace_back(argv[++a]);
      else if (s == "--nsteps" && a + 1 < argc) nsteps = std::atoi(argv[++a]);
      else if (s == "--no-output") output = false;
      else if (s == "--max_grid_size" && a + 1 < argc) mgs = std::atoi(argv[++a]);
      else if (s == "--redist_local" && a + 1 < argc) redist_local = std::atoi(argv[++a]);
      else if (s == "--ghosts" && a + 1 < argc) ng = std::atoi(argv[++a]);
      else throw std::runtime_error("unknown argument " + s);
    }
    const bool io = ParallelDescriptor::IOProcessor();
    const int myproc = ParallelDescriptor::MyProc();
    const int nranks = ParallelDescriptor::NProcs();

    rrhost::RTRun run;
    run.init(pdir, sets);
    if (nsteps >= 0) run.P.nsteps = nsteps;
    const int N = run.P.N;

    // ---- grid, decomposition, fields --------------------------------------
    const Box domain(IntVect(0), IntVect(N - 1));
    const RealBox rb({0., 0., 0.}, {double(N), double(N), double(N)});  // dx = 1
    const Array<int, 3> per{1, 1, 1};
    Geometry geom(domain, rb, CoordSys::cartesian, per);
    const Periodicity period = geom.periodicity();
    BoxArray ba(domain);
    ba.maxSize(mgs);
    DistributionMapping dm(ba);
    CFab C(ba, dm, kCComp, 0);  // RTCell
    FFab TF(ba, dm, 10, ng);  // tau(5), fgam(5) with ghosts
    LFab A(ba, dm, rrk::kNdep, ng);
    RayPC pc(geom, dm, ba);
    std::vector<BoxScratch> scratch(C.local_size());

    // ---- initial cells (host: the gas input is a host array) --------------
    {
      radray::GasGrid g;
      run.read_gas(pdir, g);
      for (MFIter mfi(C); mfi.isValid(); ++mfi) {
        auto c = C.array(mfi);
        LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
          const long long gc = i + static_cast<long long>(N) * (j + static_cast<long long>(N) * k);
          store_cell(c, i, j, k, rrhost::initial_cell(g, gc));
        });
        for (const auto& s : run.sources) {
          if (mfi.validbox().contains(IntVect(s.i, s.j, s.k))) {
            cset(c, s.i, s.j, s.k, kLuv, s.Luv);
            cset(c, s.i, s.j, s.k, kSfr, s.sfr);
          }
        }
      }
    }

    Gpu::DeviceVector<Long> hist(rrk::kNpdfMax + 1);  // bins + npos
    Gpu::DeviceVector<Long> cstat(4);                 // ncst, nspl, nmrg, ndel
    Gpu::DeviceVector<int> err(1, 0);
    Long* const hist_p = hist.data();
    Long* const cstat_p = cstat.data();
    int* const err_p = err.data();
    auto check_err = [&](const char* where) {
      int e = 0;
      Gpu::copy(Gpu::deviceToHost, err.begin(), err.end(), &e);
      ParallelDescriptor::ReduceIntMax(e);
      if (e) throw std::runtime_error(std::string("rt_amrex: index outside the ghost region / count mismatch in ") + where);
    };

    Timers T;
    Long ray_substeps = 0;
    ParallelDescriptor::Barrier();
    const double t_start = now();
    int nthreads = 1;
#ifdef AMREX_USE_OMP
    nthreads = omp_get_max_threads();
#endif

    for (int step = 1; step <= run.P.nsteps; ++step) {
      const int nsub = run.begin_step();
      const double sfr_min_step = run.units.sfr_min;  // is_src of rt_serial: fixed per hydro step
      {
        const radray::Units u = run.units;
        for (MFIter mfi(C); mfi.isValid(); ++mfi) {
          auto c = C.array(mfi);
          ParallelForOMP(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            RTCell d = load_cell(c, i, j, k);
            rrk::prep_cell(d, u);
            store_cell(c, i, j, k, d);
          });
        }
      }
      rrhost::StepStats st;
      for (int sub = 1; sub <= nsub; ++sub) {
        double t0 = now();
        // ---- 1a. reset and tau; tau -> TF ghosts ----
        {
          const radray::Units u = run.units;
          double sig[3][5];
          std::memcpy(sig, run.sig_init, sizeof sig);
          for (MFIter mfi(C); mfi.isValid(); ++mfi) {
            auto c = C.array(mfi);
            auto tf = TF.array(mfi);
            ParallelForOMP(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
              RTCell d = load_cell(c, i, j, k);
              rrk::calc_rt_1_cell(d, u, sig);
              store_cell(c, i, j, k, d);
              for (int q = 0; q < 5; ++q) tf(i, j, k, q) = d.tau[q];
            });
          }
          const double tf0 = now();
          TF.FillBoundary(0, 5, period);
          T.fillbnd += now() - tf0;
        }
        double t1 = now();
        T.cell += t1 - t0;
        // ---- 1b. global histogram -> thresholds ----
        {
          const int Np = run.N_pdf;
          zero(hist);
          for (MFIter mfi(C); mfi.isValid(); ++mfi) {
            auto& pt = pc.DefineAndReturnParticleTile(0, mfi.index(), mfi.LocalTileIndex());
            const int np = static_cast<int>(pt.numParticles());
            const PTD p = pt.getParticleTileData();
            ParallelForOMP(np, [=] AMREX_GPU_DEVICE(int n) {
              float ngv[rrk::NF];
              for (int q = 0; q < rrk::NF; ++q) ngv[q] = p.m_rdata[kNGAM + q][n];
              const double ntot = static_cast<double>(rrf::sum4<rrk::NF>(ngv));
              HostDevice::Atomic::Add(hist_p + rrk::pdf_bin(ntot, Np), Long(1));
              if (ntot > 0) HostDevice::Atomic::Add(hist_p + rrk::kNpdfMax, Long(1));
            });
          }
          std::vector<Long> h(rrk::kNpdfMax + 1);
          Gpu::copy(Gpu::deviceToHost, hist.begin(), hist.end(), h.begin());
          ParallelDescriptor::ReduceLongSum(h.data(), rrk::kNpdfMax + 1);
          std::vector<long long> hh(h.begin(), h.end() - 1);
          run.thresholds(hh.data(), h[rrk::kNpdfMax]);
        }
        const rrk::RTP P = run.rtp();
        double t2 = now();
        T.hist += t2 - t1;
        // ---- gather rays per box into cell-sorted scratch (DenseBins) ----
        for (MFIter mfi(C); mfi.isValid(); ++mfi) {
          BoxScratch& S = scratch[mfi.LocalIndex()];
          auto& pt = pc.DefineAndReturnParticleTile(0, mfi.index(), mfi.LocalTileIndex());
          const int np = static_cast<int>(pt.numParticles());
          const PTD p = pt.getParticleTileData();
          const Box bx = mfi.validbox();
          const auto lo = lbound(bx), len = length(bx);
          S.hcell.resize(np);
          int* hc = S.hcell.data();
          ParallelForOMP(np, [=] AMREX_GPU_DEVICE(int n) {
            const int i = static_cast<int>(p.m_rdata[0][n]) - lo.x;
            const int j = static_cast<int>(p.m_rdata[1][n]) - lo.y;
            const int k = static_cast<int>(p.m_rdata[2][n]) - lo.z;
            if (i < 0 || i >= len.x || j < 0 || j >= len.y || k < 0 || k >= len.z) {
              HostDevice::Atomic::Add(err_p, 1);
              hc[n] = 0;
            } else {
              hc[n] = i + len.x * (j + len.y * k);
            }
          });
#if defined(AMREX_USE_GPU) || !defined(AMREX_USE_OMP)
          S.bins.build(np, S.hcell.data(), static_cast<int>(bx.numPts()), [=] AMREX_GPU_DEVICE(int h) { return h; });
#else
          S.bins.build(BinPolicy::OpenMP, np, S.hcell.data(), static_cast<int>(bx.numPts()),
                       [=](int h) { return h; });
#endif
          S.rays.resize(np);
          RayP* rays = S.rays.data();
          const int* perm = S.bins.permutationPtr();
          ParallelForOMP(np, [=] AMREX_GPU_DEVICE(int n) { rays[n] = load_ray(p, perm[n]); });
        }
        check_err("gather");
        double t2b = now();
        T.sort += t2b - t2;
        // ---- 1d. per-cell classify/emit/split/merge: count, scan, write ----
        zero(cstat);
        const double tz = now();
        A.setVal(0);
        T.zero += now() - tz;
        for (MFIter mfi(C); mfi.isValid(); ++mfi) {
          BoxScratch& S = scratch[mfi.LocalIndex()];
          const Box bx = mfi.validbox();
          const auto lo = lbound(bx), len = length(bx);
          const int ncb = static_cast<int>(bx.numPts());
          const int np = static_cast<int>(S.rays.size());
          auto c = C.const_array(mfi);
          auto acc = A.array(mfi);
          const int* off = S.bins.offsetsPtr();
          const RayP* rays = S.rays.data();
          S.cnt.resize(ncb + 1);
          S.ostart.resize(ncb + 1);
          Long* cnt = S.cnt.data();
          ParallelForOMP(ncb + 1, [=] AMREX_GPU_DEVICE(int q) {
            if (q == ncb) {
              cnt[q] = 0;
              return;
            }
            const int i = lo.x + q % len.x, j = lo.y + (q / len.x) % len.y, k = lo.z + q / (len.x * len.y);
            float tau[5];
            for (int m = 0; m < 5; ++m) tau[m] = cf(c, i, j, k, kTau + m);
            const bool src = static_cast<double>(cf(c, i, j, k, kSfr)) > sfr_min_step;
            const bool thick = src && rrk::source_is_thick(tau);
            rrk::CellStats cs;
            cnt[q] = rrk::step1_count(rays + off[q], off[q + 1] - off[q], src, thick, i, j, k, cf(c, i, j, k, kLuv), P, cs);
            if (cs.ncst) HostDevice::Atomic::Add(cstat_p + 0, Long(cs.ncst));
            if (cs.nspl) HostDevice::Atomic::Add(cstat_p + 1, Long(cs.nspl));
            if (cs.nmrg) HostDevice::Atomic::Add(cstat_p + 2, Long(cs.nmrg));
            if (cs.ndel) HostDevice::Atomic::Add(cstat_p + 3, Long(cs.ndel));
          });
          S.nout = Scan::ExclusiveSum(ncb + 1, S.cnt.data(), S.ostart.data(), Scan::retSum);
          S.out.resize(S.nout);
          S.ocell.resize(S.nout);
          S.idx.resize(np);
          S.key.resize(np);
          const Long* ost = S.ostart.data();
          RayP* out = S.out.data();
          int* ocell = S.ocell.data();
          int* idx = S.idx.data();
          rrk::MergeKey* key = S.key.data();
          ParallelForOMP(ncb, [=] AMREX_GPU_DEVICE(int q) {
            const int i = lo.x + q % len.x, j = lo.y + (q / len.x) % len.y, k = lo.z + q / (len.x * len.y);
            float tau[5];
            for (int m = 0; m < 5; ++m) tau[m] = cf(c, i, j, k, kTau + m);
            const float Luv = cf(c, i, j, k, kLuv);
            const bool src = static_cast<double>(cf(c, i, j, k, kSfr)) > sfr_min_step;
            const bool thick = src && rrk::source_is_thick(tau);
            if (thick) {
              double v[rrk::kNdep];
              rrk::source_thick_deposit(Luv, P, v);
              for (int m = 0; m < rrk::kNdep; ++m)
                HostDevice::Atomic::Add(&acc(i, j, k, m), Long(rrk::to_fx(v[m], rrk::dep_scale(m))));
            }
            const Long n0 = ost[q], nc = ost[q + 1] - ost[q];
            if (nc == 0) return;
            const int w = rrk::step1_write(rays + off[q], off[q + 1] - off[q], src, thick, i, j, k, Luv, P,
                                           idx + off[q], key + off[q], out + n0);
            if (w != nc) HostDevice::Atomic::Add(err_p, 1);
            for (Long r = n0; r < n0 + nc; ++r) ocell[r] = q;
          });
        }
        check_err("step1");
        double t3 = now();
        T.step1 += t3 - t2b;
        // ---- 1d'. cast every output ray once (deposits incl. ghosts) ----
        for (MFIter mfi(C); mfi.isValid(); ++mfi) {
          BoxScratch& S = scratch[mfi.LocalIndex()];
          const Box bx = mfi.validbox();
          const Box gbx = amrex::grow(bx, ng);
          const auto lo = lbound(bx), len = length(bx), glo = lbound(gbx), ghi = ubound(gbx);
          auto tfa = TF.const_array(mfi);
          auto acc = A.array(mfi);
          const RayP* out = S.out.data();
          const int* ocell = S.ocell.data();
          ParallelForOMP(S.nout, [=] AMREX_GPU_DEVICE(Long r) {
            const int q = ocell[r];
            const int hi = lo.x + q % len.x, hj = lo.y + (q / len.x) % len.y, hk = lo.z + q / (len.x * len.y);
            bool bad = false;
            auto map = [&](int i, int j, int k, int& ii, int& jj, int& kk) {
              ii = near_image(i, hi, N);
              jj = near_image(j, hj, N);
              kk = near_image(k, hk, N);
              if (ii < glo.x || ii > ghi.x || jj < glo.y || jj > ghi.y || kk < glo.z || kk > ghi.z) {
                bad = true;
                ii = hi;
                jj = hj;
                kk = hk;
              }
            };
            float tb[5];
            auto tau = [&](int i, int j, int k) -> const float* {
              int ii, jj, kk;
              map(i, j, k, ii, jj, kk);
              for (int m = 0; m < 5; ++m) tb[m] = tfa(ii, jj, kk, m);
              return tb;
            };
            auto dep = [&](int i, int j, int k, const double* v) {
              int ii, jj, kk;
              map(i, j, k, ii, jj, kk);
              for (int m = 0; m < rrk::kNdep; ++m)
                HostDevice::Atomic::Add(&acc(ii, jj, kk, m), Long(rrk::to_fx(v[m], rrk::dep_scale(m))));
            };
            rrk::cast_1(out[r], P, tau, dep);
            if (bad) HostDevice::Atomic::Add(err_p, 1);
          });
          ray_substeps += S.nout;
        }
        const double tsb = now();
        A.SumBoundary(period);  // fold ghost deposits into their owners (int64: exact)
        T.sumbnd += now() - tsb;
        check_err("cast_1");
        double t4 = now();
        T.cast1 += t4 - t3;
        // ---- 2. physics; fgam -> TF ghosts ----
        {
          const rrk::PhysCtx X = run.phys();
          for (MFIter mfi(C); mfi.isValid(); ++mfi) {
            auto c = C.array(mfi);
            auto tf = TF.array(mfi);
            auto acc = A.const_array(mfi);
            ParallelForOMP(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
              RTCell d = load_cell(c, i, j, k);
              for (int m = 0; m < 5; ++m) d.ngam[m] = static_cast<float>(rrk::from_fx(acc(i, j, k, m), rrk::kScaleNgam));
              for (int m = 0; m < 5; ++m) d.nu[m] = static_cast<float>(rrk::from_fx(acc(i, j, k, 5 + m), rrk::kScaleNu));
              for (int m = 0; m < 3; ++m) d.dx[m] = static_cast<float>(rrk::from_fx(acc(i, j, k, 10 + m), rrk::kScaleNgam));
              rrk::physics_cell(d, X);
              store_cell(c, i, j, k, d);
              for (int q = 0; q < 5; ++q) tf(i, j, k, 5 + q) = d.fgam[q];
            });
          }
          const double tf0 = now();
          TF.FillBoundary(5, 5, period);
          T.fillbnd += now() - tf0;
        }
        double t5 = now();
        T.phys += t5 - t4;
        // ---- 3. advance ----
        Long nout_loc = 0, nkeep_loc = 0;
        for (MFIter mfi(C); mfi.isValid(); ++mfi) {
          BoxScratch& S = scratch[mfi.LocalIndex()];
          const Box bx = mfi.validbox();
          const Box gbx = amrex::grow(bx, ng);
          const auto lo = lbound(bx), len = length(bx), glo = lbound(gbx), ghi = ubound(gbx);
          auto tfa = TF.const_array(mfi);
          RayP* out = S.out.data();
          const int* ocell = S.ocell.data();
          S.keep.resize(S.nout + 1);
          S.wpos.resize(S.nout + 1);
          int* keep = S.keep.data();
          const Long nout = S.nout;
          ParallelForOMP(nout + 1, [=] AMREX_GPU_DEVICE(Long r) {
            if (r == nout) {
              keep[r] = 0;
              return;
            }
            const int q = ocell[r];
            const int hi = lo.x + q % len.x, hj = lo.y + (q / len.x) % len.y, hk = lo.z + q / (len.x * len.y);
            bool bad = false;
            float tb[5], fb[5];
            auto tf = [&](int i, int j, int k, const float*& tau, const float*& fgam) {
              int ii = near_image(i, hi, N), jj = near_image(j, hj, N), kk = near_image(k, hk, N);
              if (ii < glo.x || ii > ghi.x || jj < glo.y || jj > ghi.y || kk < glo.z || kk > ghi.z) {
                bad = true;
                ii = hi;
                jj = hj;
                kk = hk;
              }
              for (int m = 0; m < 5; ++m) {
                tb[m] = tfa(ii, jj, kk, m);
                fb[m] = tfa(ii, jj, kk, 5 + m);
              }
              tau = tb;
              fgam = fb;
            };
            keep[r] = rrk::cast_3(out[r], P, tf) ? 1 : 0;
            if (bad) HostDevice::Atomic::Add(err_p, 1);
          });
          nout_loc += nout;
        }
        check_err("cast_3");
        double t6 = now();
        T.cast3 += t6 - t5;
        // ---- write survivors back into the particle tiles ----
        for (MFIter mfi(C); mfi.isValid(); ++mfi) {
          BoxScratch& S = scratch[mfi.LocalIndex()];
          const Long nout = S.nout;
          const Long nkeep = Scan::ExclusiveSum(nout + 1, S.keep.data(), S.wpos.data(), Scan::retSum);
          nkeep_loc += nkeep;
          auto& pt = pc.DefineAndReturnParticleTile(0, mfi.index(), mfi.LocalTileIndex());
          pt.resize(nkeep);
          const PTD p = pt.getParticleTileData();
          const int* keep = S.keep.data();
          const int* wpos = S.wpos.data();
          const RayP* out = S.out.data();
          const int* ocell = S.ocell.data();
          const auto lo = lbound(mfi.validbox()), len = length(mfi.validbox());
          ParallelForOMP(nout, [=] AMREX_GPU_DEVICE(Long r) {
            if (!keep[r]) return;
            const int q = ocell[r];
            const int h[3] = {lo.x + q % len.x, lo.y + (q / len.x) % len.y, lo.z + q / (len.x * len.y)};
            store_ray(p, wpos[r], out[r], N, myproc, h);
          });
        }
        double t7 = now();
        T.sort += t7 - t6;
        // ---- migrate rays to the owners of their new home cells ----
        if (redist_local > 0)
          pc.Redistribute(0, 0, 0, redist_local);
        else
          pc.Redistribute();
        T.migrate += now() - t7;
        {
          Long s[6];
          Gpu::copy(Gpu::deviceToHost, cstat.begin(), cstat.end(), s);
          s[4] = nout_loc - nkeep_loc;  // deleted in cast_3
          s[5] = 0;
          ParallelDescriptor::ReduceLongSum(s, 5);
          st.ncst1 = s[0];
          st.nspl1 = s[1];
          st.nmrg1 = s[2];
          st.ndel1 = s[3] + s[4];
        }
        run.end_substep();
      }
      // ---- post_rt: redshift rays, final gas/tau, ngam_tot ----
      const double t8 = now();
      {
        const rrk::RedshiftTable RT = run.redshift_table();
        const radray::Units u = run.units;
        double sig[3][5];
        std::memcpy(sig, run.sig_init, sizeof sig);
        Long s[2] = {0, 0};
        for (MFIter mfi(C); mfi.isValid(); ++mfi) {
          auto& pt = pc.DefineAndReturnParticleTile(0, mfi.index(), mfi.LocalTileIndex());
          const int np = static_cast<int>(pt.numParticles());
          const PTD p = pt.getParticleTileData();
          ParallelForOMP(np, [=] AMREX_GPU_DEVICE(int n) {
            RayP r = load_ray(p, n);
            rrk::redshift_ray(r, RT);
            store_payload(p, n, r);
          });
          s[0] += Reduce::Sum<Long>(np, [=] AMREX_GPU_DEVICE(int n) -> Long {
            Long t = 0;
            for (int m = 0; m < 5; ++m)
              t += rrk::to_fx(static_cast<double>(p.m_rdata[kNGAM + m][n]), rrk::kScaleTot);
            return t;
          });
          s[1] += np;
          auto c = C.array(mfi);
          ParallelForOMP(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            RTCell d = load_cell(c, i, j, k);
            rrk::post_cell(d, u, sig);
            store_cell(c, i, j, k, d);
          });
        }
        ParallelDescriptor::ReduceLongSum(s, 2);
        st.ngam_tot_fx = s[0];
        st.nray_active = s[1];
      }
      T.post += now() - t8;
      if (io)
        std::printf(" step %3d  z = %.5f  nstep_rt = %d  rays = %lld  cast = %lld  merge = %lld  split = %lld\n",
                    run.clk.istep, run.clk.redshift, nsub, st.nray_active, st.ncst1, st.nmrg1, st.nspl1);
      if (output) {
        // gather all cells on rank 0 in global order
        BoxArray ba0(domain);
        DistributionMapping dm0(Vector<int>{0});
        // ParallelCopy needs an arithmetic FAB (it instantiates an ADD path
        // that requires operator+= on the element type), so the AoS cells
        // go through a 45-component float FabArray.
        FFab Cf(ba, dm, NC, 0);
        for (MFIter mfi(C); mfi.isValid(); ++mfi) {
          auto c = C.const_array(mfi);
          auto f = Cf.array(mfi);
          ParallelForOMP(mfi.validbox(), NC, [=] AMREX_GPU_DEVICE(int i, int j, int k, int m) {
            f(i, j, k, m) = cf(c, i, j, k, m);
          });
        }
        FFab C0(ba0, dm0, NC, 0);
        C0.ParallelCopy(Cf, 0, 0, NC);
        if (io) {
          const auto a = C0.const_array(0);
          std::vector<RTCell> cells(static_cast<std::size_t>(N) * N * N);
          for (int k = 0; k < N; ++k)
            for (int j = 0; j < N; ++j)
              for (int i = 0; i < N; ++i) {
                float b[NC];
                for (int m = 0; m < NC; ++m) b[m] = a(i, j, k, m);
                std::memcpy(&cells[i + static_cast<std::size_t>(N) * (j + static_cast<std::size_t>(N) * k)], b,
                            sizeof b);
              }
          run.write_step(odir, cells, st);
        }
      }
    }
    ParallelDescriptor::Barrier();
    const double tt = now() - t_start;
    ParallelDescriptor::ReduceLongSum(ray_substeps);
    if (io) {
      std::printf(" ranks %d  threads %d  boxes %d  wall %.3f s  ray-substeps %lld  ns/ray-substep*threads*ranks %.1f\n",
                  nranks, nthreads, static_cast<int>(ba.size()), tt, static_cast<long long>(ray_substeps),
                  1e9 * tt * nthreads * nranks / (ray_substeps > 0 ? ray_substeps : 1));
      std::printf(" timers: cell %.3f hist %.3f step1 %.3f cast1 %.3f phys %.3f cast3 %.3f sort %.3f migrate %.3f post %.3f\n",
                  T.cell, T.hist, T.step1, T.cast1, T.phys, T.cast3, T.sort, T.migrate, T.post);
      std::printf(" of which: FillBoundary(tau,fgam) %.3f SumBoundary(dep) %.3f setVal(dep) %.3f\n", T.fillbnd,
                  T.sumbnd, T.zero);
    }
  } catch (const std::exception& e) {
    std::fprintf(stderr, "rt_amrex: error: %s\n", e.what());
    if (ParallelDescriptor::NProcs() > 1) amrex::Abort(e.what());
    rc = 1;
  }
  amrex::Finalize();
  return rc;
}
