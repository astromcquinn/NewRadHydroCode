// Parthenon prototype of the sorted-ray-array RT substep: exactly the phase
// sequence of spike/serial/rt_serial.cpp, distributed over meshblocks, with
// the same kernels (spike/common) and host logic (rrhost::RTRun), and
// bit-identical results for any rank count, meshblock size and thread count.
//
//   cells   : Parthenon variables "rtc" (35 comps: CellData without tau,
//             fgam), "tau" (5) and "fgam" (5); tau/fgam have
//             Metadata::FillGhost (nghost = 4) and are exchanged with
//             Parthenon's boundary tasks before cast_1 / cast_3.
//   deposits: per block int64 fixed-point accumulators incl. a 4-cell ghost
//             zone (Kokkos::View, Kokkos::atomic_add), then a hand-written
//             reverse ("ghost -> owner, summed") exchange over the block's
//             neighbour list (same-rank: direct add, else MPI).
//   rays    : Parthenon Swarm "rays" (int "ri"(2) = nside, ipix; Real
//             "rf"(18) = the 18 floats of rrk::RayP); swarm position = centre
//             of the ray's home cell (image nearest the block); migration by
//             swarm boundary communication after cast_3.
//   global  : MPI_Allreduce on integers only (histogram, counts, fixed point).
//
// usage: rt_parthenon -i rt.in <problem_dir> <outdir> [--nsteps n] [--set k=v]...
//                     [--no-output] [--mb 16] [block/par=value ...]
#include <cstddef>

#include "parth_common.hpp"

#include "rt_host.hpp"

using namespace parthenon;
using rrk::CellData;
using rrk::RayP;

namespace {

constexpr int G = 4;  // ghost cells (tau/fgam reads and deposits reach <= 3)
constexpr int NDEP = rrk::kNdep;

// CellData as 45 floats: [0,26) and [36,45) go to "rtc", fgam = [26,31),
// tau = [31,36).
constexpr int kNcd = 45, kFg = 26, kTau = 31, kNrtc = 35;
static_assert(sizeof(CellData) == kNcd * sizeof(float), "CellData layout");
static_assert(offsetof(CellData, fgam) == kFg * sizeof(float), "CellData layout");
static_assert(offsetof(CellData, tau) == kTau * sizeof(float), "CellData layout");
static_assert(offsetof(CellData, u) == (kTau + 5) * sizeof(float), "CellData layout");

using PVar = ParArrayND<Real, VariableState>;

// Bit masks over the 45 floats of CellData: the kernels touch only part of
// a cell, and with Parthenon's component-slowest layout every component is a
// separate memory stream, so each phase loads/stores only what it uses.
constexpr unsigned long long rng(int a, int b) { return ((1ULL << b) - 1) & ~((1ULL << a) - 1); }
constexpr unsigned long long kAll = rng(0, kNcd);
constexpr unsigned long long kGas = rng(0, 5);           // rhog fHI fHeI fHeII temp
constexpr unsigned long long kSrc = rng(5, 7);           // Luv sfr
constexpr unsigned long long kDx = rng(7, 10), kGamLam = rng(10, 16), kNu = rng(16, 21), kNgam = rng(21, 26);
constexpr unsigned long long kFgam = rng(kFg, kFg + 5), kTauM = rng(kTau, kTau + 5);
constexpr unsigned long long kUW = rng(kTau + 5, kNcd);  // u_gas, w_gas

template <class A>
struct CellVarsT {
  A rtc, tau, fgam;
  int is, js, ks;  // first interior index
  KOKKOS_FORCEINLINE_FUNCTION Real& at(int q, int K, int J, int I) const {
    if (q < kFg) return rtc(q, K, J, I);
    if (q < kTau) return fgam(q - kFg, K, J, I);
    if (q < kTau + 5) return tau(q - kTau, K, J, I);
    return rtc(q - 10, K, J, I);
  }
  template <unsigned long long M = kAll>
  KOKKOS_FORCEINLINE_FUNCTION CellData load(int i, int j, int k) const {
    CellData c{};
    float* f = reinterpret_cast<float*>(&c);
    const int I = i + is, J = j + js, K = k + ks;
    for (int q = 0; q < kNcd; ++q)
      if ((M >> q) & 1ULL) f[q] = static_cast<float>(at(q, K, J, I));
    return c;
  }
  template <unsigned long long M = kAll>
  KOKKOS_FORCEINLINE_FUNCTION void store(int i, int j, int k, const CellData& c) const {
    const float* f = reinterpret_cast<const float*>(&c);
    const int I = i + is, J = j + js, K = k + ks;
    for (int q = 0; q < kNcd; ++q)
      if ((M >> q) & 1ULL) at(q, K, J, I) = f[q];
  }
};
using CellVars = CellVarsT<PVar>;
template <class A>
CellVarsT<A> make_cv(A rtc, A tau, A fgam, int is, int js, int ks) {
  return CellVarsT<A>{rtc, tau, fgam, is, js, ks};
}

// Ray payload <-> swarm variables
using IVar = ParArrayND<int>;
using RVar = ParArrayND<Real>;
KOKKOS_INLINE_FUNCTION RayP load_ray(const IVar& ri, const RVar& rf, int n) {
  RayP r;
  r.nside = ri(0, n);
  r.ipix = ri(1, n);
  float* f = reinterpret_cast<float*>(reinterpret_cast<char*>(&r) + 2 * sizeof(std::int32_t));
  for (int q = 0; q < 18; ++q) f[q] = static_cast<float>(rf(q, n));
  return r;
}
KOKKOS_INLINE_FUNCTION void store_ray(const IVar& ri, const RVar& rf, int n, const RayP& r) {
  ri(0, n) = r.nside;
  ri(1, n) = r.ipix;
  const float* f = reinterpret_cast<const float*>(reinterpret_cast<const char*>(&r) + 2 * sizeof(std::int32_t));
  for (int q = 0; q < 18; ++q) rf(q, n) = f[q];
}

// Wrapped global cell index -> index in the block's ghost-grown range
// [0, n + 2G) (nearest periodic image), or -1 if no image is in range.
KOKKOS_INLINE_FUNCTION int local_img(int gi, int off, int n, int N) {
  const int d = gi - off;
  if (d >= -G && d < n + G) return d + G;
  if (d + N >= -G && d + N < n + G) return d + N + G;
  if (d - N >= -G && d - N < n + G) return d - N + G;
  return -1;
}

struct Sig {
  double s[3][5];
};

template <class T>
using V1 = Kokkos::View<T*>;
using Acc4 = Kokkos::View<long long****, Kokkos::LayoutRight>;  // (k,j,i,13) ghost-grown
using TF4 = Kokkos::View<float****, Kokkos::LayoutRight>;       // (k,j,i,10) tau|fgam

template <class T>
void ensure(V1<T>& v, std::size_t n, const char* name) {
  if (v.extent(0) < n) v = V1<T>(Kokkos::view_alloc(Kokkos::WithoutInitializing, std::string(name)), n + n / 4 + 16);
}

struct Blk {
  MeshBlock* pmb;
  sp::BlockGeom g;
  int nx, ny, nz, ncell;
  CellVars cv;
  Acc4 acc;
  TF4 tf;
  V1<RayP> rays, out;
  V1<long long> start, ostart, cnt, fill;
  V1<int> pcell, idx, keep, kpos;
  V1<rrk::MergeKey> key;
  V1<long long> hist, stats;  // hist[kNpdfMax] + npos; stats: ncst nspl nmrg ndel ndel3
  long long nr = 0, nout = 0;
};

Packages_t process_packages(std::unique_ptr<ParameterInput>& pin) {
  Packages_t packages;
  auto pkg = std::make_shared<StateDescriptor>("spike_rt");
  pkg->AddField("rtc", Metadata({Metadata::Cell, Metadata::Derived, Metadata::OneCopy}, std::vector<int>{kNrtc}));
  Metadata mg({Metadata::Cell, Metadata::Derived, Metadata::OneCopy, Metadata::FillGhost}, std::vector<int>{5});
  pkg->AddField("tau", mg);
  pkg->AddField("fgam", mg);
  pkg->AddSwarm("rays", Metadata({Metadata::Provides, Metadata::None, Metadata::NoPersistentParticleIds}));
  pkg->AddSwarmValue("ri", "rays", Metadata({Metadata::Integer}, std::vector<int>{2}));
  pkg->AddSwarmValue("rf", "rays", Metadata({Metadata::Real}, std::vector<int>{18}));
  packages.Add(pkg);
  return packages;
}

std::shared_ptr<Swarm> get_swarm(MeshBlock* pmb) {
  return pmb->meshblock_data.Get()->GetSwarmData()->Get("rays");
}

// tau | fgam float mirror (incl. ghosts) of the Parthenon variables
void mirror_tf(Blk& B, bool do_tau, bool do_fgam) {
  auto tf = B.tf;
  auto tau = B.cv.tau, fgam = B.cv.fgam;
  const int mz = B.nz + 2 * G, my = B.ny + 2 * G, mx = B.nx + 2 * G;
  B.pmb->par_for(
      "rt_mirror_tf", 0, mz - 1, 0, my - 1, 0, mx - 1, KOKKOS_LAMBDA(const int k, const int j, const int i) {
        if (do_tau)
          for (int q = 0; q < 5; ++q) tf(k, j, i, q) = static_cast<float>(tau(q, k, j, i));
        if (do_fgam)
          for (int q = 0; q < 5; ++q) tf(k, j, i, 5 + q) = static_cast<float>(fgam(q, k, j, i));
      });
}

// ---------------------------------------------------------------------------
// Reverse ghost exchange: add every block's deposit ghost zone into the
// interior of the neighbour that owns it. Integer adds, so the result does not
// depend on the order of messages or blocks.
// ---------------------------------------------------------------------------
struct Box {
  int lo[3], hi[3];
  int vol() const { return (hi[0] - lo[0]) * (hi[1] - lo[1]) * (hi[2] - lo[2]); }
};
// sender's ghost slab in direction o (sender's ghost-grown coordinates)
Box send_box(const Blk& B, const int o[3]) {
  const int n[3] = {B.nx, B.ny, B.nz};
  Box b;
  for (int d = 0; d < 3; ++d) {
    b.lo[d] = o[d] < 0 ? 0 : (o[d] == 0 ? G : G + n[d]);
    b.hi[d] = o[d] < 0 ? G : (o[d] == 0 ? G + n[d] : 2 * G + n[d]);
  }
  return b;
}
// receiver's interior slab that the sender's ghost slab o maps onto
Box recv_box(const Blk& B, const int o[3]) {
  const int n[3] = {B.nx, B.ny, B.nz};
  Box b;
  for (int d = 0; d < 3; ++d) {
    b.lo[d] = o[d] < 0 ? n[d] : G;
    b.hi[d] = o[d] < 0 ? n[d] + G : (o[d] == 0 ? G + n[d] : 2 * G);
  }
  return b;
}
int dir_index(const int o[3]) { return (o[0] + 1) + 3 * (o[1] + 1) + 9 * (o[2] + 1); }

void pack_box(const Acc4& acc, const Box& b, const V1<long long>& buf) {
  const int lx = b.hi[0] - b.lo[0], ly = b.hi[1] - b.lo[1];
  const int x0 = b.lo[0], y0 = b.lo[1], z0 = b.lo[2];
  parthenon::par_for(
      "rt_fold_pack", b.lo[2], b.hi[2] - 1, b.lo[1], b.hi[1] - 1, b.lo[0], b.hi[0] - 1,
      KOKKOS_LAMBDA(const int k, const int j, const int i) {
        const long long c = (i - x0) + lx * ((j - y0) + static_cast<long long>(ly) * (k - z0));
        for (int m = 0; m < NDEP; ++m) buf(c * NDEP + m) = acc(k, j, i, m);
      });
}
void add_box(const Acc4& acc, const Box& b, const V1<long long>& buf) {
  const int lx = b.hi[0] - b.lo[0], ly = b.hi[1] - b.lo[1];
  const int x0 = b.lo[0], y0 = b.lo[1], z0 = b.lo[2];
  parthenon::par_for(
      "rt_fold_add", b.lo[2], b.hi[2] - 1, b.lo[1], b.hi[1] - 1, b.lo[0], b.hi[0] - 1,
      KOKKOS_LAMBDA(const int k, const int j, const int i) {
        const long long c = (i - x0) + lx * ((j - y0) + static_cast<long long>(ly) * (k - z0));
        for (int m = 0; m < NDEP; ++m) acc(k, j, i, m) += buf(c * NDEP + m);
      });
}

struct Fold {
  std::vector<V1<long long>> sbuf, rbuf;  // one per (block, neighbour)
  std::vector<int> rblk, sblk, snb, rnb;
};

void fold_deposits(std::vector<Blk>& blks, const std::vector<int>& lid_of_gid, int rank, MPI_Comm comm,
                   Fold& F) {
  std::vector<MPI_Request> req;
  // post receives
  std::size_t ir = 0;
  for (std::size_t b = 0; b < blks.size(); ++b) {
    for (auto& nb : blks[b].pmb->neighbors) {
      if (nb.rank == rank) continue;
      int o[3];
      for (int d = 0; d < 3; ++d) o[d] = -static_cast<int>(nb.offsets[d]);  // sender's offset to me
      const Box rb = recv_box(blks[b], o);
      if (F.rbuf.size() <= ir) F.rbuf.emplace_back("rt_fold_rbuf", 1);
      ensure(F.rbuf[ir], static_cast<std::size_t>(rb.vol()) * NDEP, "rt_fold_rbuf");
      req.emplace_back();
      MPI_Irecv(F.rbuf[ir].data(), rb.vol() * NDEP, MPI_LONG_LONG, nb.rank, blks[b].pmb->gid * 27 + dir_index(o),
                comm, &req.back());
      ++ir;
    }
  }
  const std::size_t nrecv = req.size();
  // pack and send / add locally
  std::size_t is = 0;
  for (std::size_t b = 0; b < blks.size(); ++b) {
    for (auto& nb : blks[b].pmb->neighbors) {
      int o[3];
      for (int d = 0; d < 3; ++d) o[d] = static_cast<int>(nb.offsets[d]);
      const Box sb = send_box(blks[b], o);
      if (F.sbuf.size() <= is) F.sbuf.emplace_back("rt_fold_sbuf", 1);
      ensure(F.sbuf[is], static_cast<std::size_t>(sb.vol()) * NDEP, "rt_fold_sbuf");
      pack_box(blks[b].acc, sb, F.sbuf[is]);
      if (nb.rank == rank) {
        const int t = lid_of_gid[nb.gid];
        if (t < 0) sp::die("fold: local neighbour not found");
        add_box(blks[t].acc, recv_box(blks[t], o), F.sbuf[is]);
      } else {
        Kokkos::fence();
        req.emplace_back();
        MPI_Isend(F.sbuf[is].data(), sb.vol() * NDEP, MPI_LONG_LONG, nb.rank, nb.gid * 27 + dir_index(o), comm,
                  &req.back());
      }
      ++is;
    }
  }
  // wait for the receives and add them (in neighbour order; integer adds)
  MPI_Waitall(static_cast<int>(nrecv), req.data(), MPI_STATUSES_IGNORE);
  ir = 0;
  for (std::size_t b = 0; b < blks.size(); ++b) {
    for (auto& nb : blks[b].pmb->neighbors) {
      if (nb.rank == rank) continue;
      int o[3];
      for (int d = 0; d < 3; ++d) o[d] = -static_cast<int>(nb.offsets[d]);
      add_box(blks[b].acc, recv_box(blks[b], o), F.rbuf[ir]);
      ++ir;
    }
  }
  MPI_Waitall(static_cast<int>(req.size() - nrecv), req.data() + nrecv, MPI_STATUSES_IGNORE);
  Kokkos::fence();
}

// Swarm migration: Parthenon swarm boundary communication. Rays move at most
// ~3 cells, so one hop to a neighbour is always enough.
void migrate_rays(Mesh* pmesh) {
  TaskCollection tc;
  TaskID none(0);
  auto& blocks = pmesh->block_list;
  const int nb = static_cast<int>(blocks.size());
  TaskRegion& r0 = tc.AddRegion(nb);  // reset all flags before any block sends
  for (int b = 0; b < nb; ++b) {
    auto sc = blocks[b]->meshblock_data.Get()->GetSwarmData();
    r0[b].AddTask(none, &SwarmContainer::ResetCommunication, sc.get());
  }
  TaskRegion& r1 = tc.AddRegion(nb);
  for (int b = 0; b < nb; ++b) {
    auto sc = blocks[b]->meshblock_data.Get()->GetSwarmData();
    auto send = r1[b].AddTask(none, &SwarmContainer::Send, sc.get(), BoundaryCommSubset::all);
    r1[b].AddTask(send, &SwarmContainer::Receive, sc.get(), BoundaryCommSubset::all);
  }
  if (tc.Execute() != TaskListStatus::complete) sp::die("ray migration failed");
}

struct Timers {
  double cell = 0, halo = 0, hist = 0, step1 = 0, cast1 = 0, fold = 0, phys = 0, cast3 = 0, sort = 0,
         migrate = 0, post = 0;
};

}  // namespace

int main(int argc, char** argv) {
  std::vector<std::string> pos, sets;
  int nsteps = -1, mb = 16;
  bool output = true;
  for (int a = 1; a < argc; ++a) {
    const std::string s = argv[a];
    if (s == "-i") ++a;
    else if (s == "--set" && a + 1 < argc) sets.emplace_back(argv[++a]);
    else if (s == "--nsteps" && a + 1 < argc) nsteps = std::atoi(argv[++a]);
    else if (s == "--mb" && a + 1 < argc) mb = std::atoi(argv[++a]);
    else if (s == "--no-output") output = false;
    else if (s.find('=') != std::string::npos && s.find('/') < s.find('=')) continue;  // block/par=value
    else if (!s.empty() && s[0] != '-') pos.push_back(s);
    else {
      std::fprintf(stderr, "unknown argument %s\n", s.c_str());
      return 2;
    }
  }
  if (pos.size() != 2) {
    std::fprintf(stderr,
                 "usage: %s -i rt.in <problem_dir> <outdir> [--nsteps n] [--set k=v] [--no-output] [--mb 16]\n",
                 argv[0]);
    return 2;
  }
  const std::string pdir = pos[0], odir = pos[1];

  sp::ParthArgs pargs = sp::parthenon_args(argc, argv);
  ParthenonManager pman;
  auto pst = pman.ParthenonInitEnv(pargs.argc(), pargs.argv());
  if (pst != ParthenonStatus::ok) {
    pman.ParthenonFinalize();
    return pst == ParthenonStatus::complete ? 0 : 1;
  }
  int status = 0;
  try {
    rrhost::RTRun run;
    run.init(pdir, sets);
    if (nsteps >= 0) run.P.nsteps = nsteps;
    const int N = run.P.N;
    {
      auto& pin = pman.pinput;
      const char* dn[3] = {"1", "2", "3"};
      for (int d = 0; d < 3; ++d) {
        pin->SetInteger("parthenon/mesh", std::string("nx") + dn[d], N);
        pin->SetReal("parthenon/mesh", std::string("x") + dn[d] + "min", 0.0);
        pin->SetReal("parthenon/mesh", std::string("x") + dn[d] + "max", static_cast<double>(N));
        pin->SetString("parthenon/mesh", std::string("ix") + dn[d] + "_bc", "periodic");
        pin->SetString("parthenon/mesh", std::string("ox") + dn[d] + "_bc", "periodic");
        pin->SetInteger("parthenon/meshblock", std::string("nx") + dn[d], mb);
      }
      pin->SetInteger("parthenon/mesh", "nghost", G);
      pin->SetString("parthenon/mesh", "refinement", "none");
    }
    if (mb < 2 * G) throw std::runtime_error("meshblocks must have >= 8 cells per side");
    pman.app_input->ProcessPackages = process_packages;
    pman.ParthenonInitPackagesAndMesh();
    Mesh* pmesh = pman.pmesh.get();

    MPI_Comm comm;
    MPI_Comm_dup(MPI_COMM_WORLD, &comm);
    int rank, nranks;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nranks);
    const long long ncell_glob = static_cast<long long>(N) * N * N;

    // rate tables in Kokkos views (device-visible pointers for PhysCtx)
    V1<double> rt_cool("cool", run.rates.coef_cool_.size()), rt_alpha("alpha", run.rates.alpha_rec_.size()),
        rt_gcol("gcol", run.rates.gamma_col_.size());
    {
      auto cp = [](V1<double>& v, const std::vector<double>& h) {
        auto m = Kokkos::create_mirror_view(v);
        for (std::size_t q = 0; q < h.size(); ++q) m(q) = h[q];
        Kokkos::deep_copy(v, m);
      };
      cp(rt_cool, run.rates.coef_cool_);
      cp(rt_alpha, run.rates.alpha_rec_);
      cp(rt_gcol, run.rates.gamma_col_);
    }
    Sig sig;
    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 5; ++b) sig.s[a][b] = run.sig_init[a][b];

    // blocks
    auto& blocks = pmesh->block_list;
    const int nb = static_cast<int>(blocks.size());
    std::vector<Blk> blks(nb);
    std::vector<int> lid_of_gid(pmesh->nbtotal, -1);
    {
      radray::GasGrid gas;
      run.read_gas(pdir, gas);
      std::vector<float> srcL(ncell_glob, 0.f), srcS(ncell_glob, 0.f);
      for (const auto& s : run.sources) {
        const long long c = s.i + static_cast<long long>(N) * (s.j + static_cast<long long>(N) * s.k);
        srcL[c] = s.Luv;
        srcS[c] = s.sfr;
      }
      for (int b = 0; b < nb; ++b) {
        Blk& B = blks[b];
        B.pmb = blocks[b].get();
        lid_of_gid[B.pmb->gid] = b;
        B.g = sp::block_geom(B.pmb, N, 0.0, static_cast<double>(N));
        B.nx = B.g.n[0];
        B.ny = B.g.n[1];
        B.nz = B.g.n[2];
        B.ncell = B.nx * B.ny * B.nz;
        auto rc = B.pmb->meshblock_data.Get();
        B.cv = CellVars{rc->Get("rtc").data, rc->Get("tau").data, rc->Get("fgam").data, B.g.is, B.g.js, B.g.ks};
        B.acc = Acc4("rt_acc", B.nz + 2 * G, B.ny + 2 * G, B.nx + 2 * G, NDEP);
        B.tf = TF4("rt_tf", B.nz + 2 * G, B.ny + 2 * G, B.nx + 2 * G, 10);
        B.start = V1<long long>("start", B.ncell + 1);
        B.ostart = V1<long long>("ostart", B.ncell + 1);
        B.cnt = V1<long long>("cnt", B.ncell);
        B.fill = V1<long long>("fill", B.ncell);
        B.hist = V1<long long>("hist", rrk::kNpdfMax + 1);
        B.stats = V1<long long>("stats", 5);
        B.rays = V1<RayP>("rays", 16);
        B.out = V1<RayP>("out", 16);
        B.pcell = V1<int>("pcell", 16);
        B.idx = V1<int>("idx", 16);
        B.keep = V1<int>("keep", 16);
        B.kpos = V1<int>("kpos", 16);
        B.key = V1<rrk::MergeKey>("key", 16);
        // initial cells (host fill, then copy into the Parthenon variables)
        auto hrtc = B.cv.rtc.GetHostMirror();
        auto htau = B.cv.tau.GetHostMirror();
        auto hfg = B.cv.fgam.GetHostMirror();
        auto hv = make_cv(hrtc, htau, hfg, B.g.is, B.g.js, B.g.ks);
        for (int k = 0; k < B.nz; ++k)
          for (int j = 0; j < B.ny; ++j)
            for (int i = 0; i < B.nx; ++i) {
              const long long c = (B.g.off[0] + i) + static_cast<long long>(N) * ((B.g.off[1] + j) + static_cast<long long>(N) * (B.g.off[2] + k));
              CellData d = rrhost::initial_cell(gas, c);
              d.Luv = srcL[c];  // 0 except in source cells, as in rt_serial
              d.sfr = srcS[c];
              hv.store(i, j, k, d);
            }
        B.cv.rtc.DeepCopy(hrtc);
        B.cv.tau.DeepCopy(htau);
        B.cv.fgam.DeepCopy(hfg);
      }
    }

    Timers T;
    Fold F;
    long long ray_substeps = 0;
    const int nthreads = sp::num_threads();
    MPI_Barrier(comm);
    const double t_start = sp::now();

    for (int step = 1; step <= run.P.nsteps; ++step) {
      const int nsub = run.begin_step();
      {
        const radray::Units u = run.units;
        for (auto& B : blks) {
          const CellVars cv = B.cv;
          const int nx = B.nx, ny = B.ny;
          B.pmb->par_for(
              "rt_prep", 0, B.ncell - 1, KOKKOS_LAMBDA(const int c) {
                const int i = c % nx, j = (c / nx) % ny, k = c / (nx * ny);
                CellData d = cv.load<kUW>(i, j, k);
                rrk::prep_cell(d, u);
                cv.store<kGas>(i, j, k, d);
              });
        }
      }
      rrhost::StepStats st;
      for (int sub = 1; sub <= nsub; ++sub) {
        // 0. rays: swarm -> cell-sorted array (counting sort), swarm emptied
        double t0 = sp::now();
        for (auto& B : blks) {
          auto swarm = get_swarm(B.pmb);
          const int nmax = swarm->GetMaxActiveIndex();
          const long long nact = swarm->GetNumActive();
          ensure(B.rays, nact, "rays");
          ensure(B.pcell, nmax + 1, "pcell");
          Kokkos::deep_copy(B.cnt, 0);
          Kokkos::deep_copy(B.fill, 0);
          auto ri = swarm->Get<int>("ri").Get();
          auto rf = swarm->Get<Real>("rf").Get();
          auto sd = swarm->GetDeviceContext();
          auto cnt = B.cnt, fill = B.fill, start = B.start;
          auto pcell = B.pcell;
          auto rays = B.rays;
          const int nx = B.nx, ny = B.ny, nz = B.nz, ox = B.g.off[0], oy = B.g.off[1], oz = B.g.off[2];
          if (nmax >= 0) {
            B.pmb->par_for(
                "rt_home", 0, nmax, KOKKOS_LAMBDA(const int n) {
                  if (!sd.IsActive(n)) return;
                  const RayP r = load_ray(ri, rf, n);
                  const int li = local_img(rrk::cell_of(r.x[0], N), ox, nx, N) - G;
                  const int lj = local_img(rrk::cell_of(r.x[1], N), oy, ny, N) - G;
                  const int lk = local_img(rrk::cell_of(r.x[2], N), oz, nz, N) - G;
                  if (li < 0 || li >= nx || lj < 0 || lj >= ny || lk < 0 || lk >= nz)
                    Kokkos::abort("ray not in its home block after migration");
                  const int c = li + nx * (lj + ny * lk);
                  pcell(n) = c;
                  Kokkos::atomic_add(&cnt(c), 1LL);
                });
          }
          const int ncell = B.ncell;
          Kokkos::parallel_scan(
              "rt_scan_home", ncell, KOKKOS_LAMBDA(const int c, long long& upd, const bool fin) {
                if (fin) start(c) = upd;
                upd += cnt(c);
                if (fin && c == ncell - 1) start(ncell) = upd;
              });
          if (nmax >= 0) {
            B.pmb->par_for(
                "rt_scatter", 0, nmax, KOKKOS_LAMBDA(const int n) {
                  if (!sd.IsActive(n)) return;
                  const int c = pcell(n);
                  const long long p = start(c) + Kokkos::atomic_fetch_add(&fill(c), 1LL);
                  rays(p) = load_ray(ri, rf, n);
                  sd.MarkParticleForRemoval(n);
                });
          }
          swarm->RemoveMarkedParticles();
          auto hs = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), Kokkos::subview(start, ncell));
          B.nr = hs();
          if (B.nr != nact) sp::die("ray count mismatch in gather");
        }
        double t1 = sp::now();
        T.sort += t1 - t0;
        t0 = t1;
        // 1a. reset and tau, tau ghosts
        {
          const radray::Units u = run.units;
          for (auto& B : blks) {
            const CellVars cv = B.cv;
            const int nx = B.nx, ny = B.ny;
            B.pmb->par_for(
                "rt_calc_rt_1", 0, B.ncell - 1, KOKKOS_LAMBDA(const int c) {
                  const int i = c % nx, j = (c / nx) % ny, k = c / (nx * ny);
                  CellData d = cv.load<rng(0, 4)>(i, j, k);  // rhog fHI fHeI fHeII
                  rrk::calc_rt_1_cell(d, u, sig.s);
                  cv.store<rng(7, kTau + 5)>(i, j, k, d);    // dx .. tau
                });
          }
          Kokkos::fence();
        }
        t1 = sp::now();
        T.cell += t1 - t0;
        t0 = t1;
        sp::exchange_ghosts(pmesh, "tau_only", {"tau"});
        for (auto& B : blks) mirror_tf(B, true, false);
        Kokkos::fence();
        t1 = sp::now();
        T.halo += t1 - t0;
        t0 = t1;
        // 1b. histogram -> thresholds (exact integer allreduce)
        std::vector<long long> hist(rrk::kNpdfMax + 1, 0);
        {
          const int Np = run.N_pdf;
          for (auto& B : blks) {
            Kokkos::deep_copy(B.hist, 0);
            auto h = B.hist;
            auto rays = B.rays;
            if (B.nr > 0)
              Kokkos::parallel_for(
                  "rt_hist", B.nr, KOKKOS_LAMBDA(const long long r) {
                    const double ntot = rrk::sum_ngam(rays(r));
                    Kokkos::atomic_add(&h(rrk::pdf_bin(ntot, Np)), 1LL);
                    if (ntot > 0) Kokkos::atomic_add(&h(rrk::kNpdfMax), 1LL);
                  });
            auto hh = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), h);
            for (int q = 0; q <= rrk::kNpdfMax; ++q) hist[q] += hh(q);
          }
          MPI_Allreduce(MPI_IN_PLACE, hist.data(), rrk::kNpdfMax + 1, MPI_LONG_LONG, MPI_SUM, comm);
        }
        run.thresholds(hist.data(), hist[rrk::kNpdfMax]);
        const rrk::RTP P = run.rtp();
        t1 = sp::now();
        T.hist += t1 - t0;
        t0 = t1;
        // 1d. per-cell classify/emit/split/merge: count, scan, write
        const double sfr_min = run.units.sfr_min;
        for (auto& B : blks) {
          const CellVars cv = B.cv;
          const int nx = B.nx, ny = B.ny, ox = B.g.off[0], oy = B.g.off[1], oz = B.g.off[2];
          auto rays = B.rays;
          auto start = B.start, ostart = B.ostart, cnt = B.cnt, stats = B.stats;
          Kokkos::deep_copy(stats, 0);
          B.pmb->par_for(
              "rt_step1_count", 0, B.ncell - 1, KOKKOS_LAMBDA(const int c) {
                const int i = c % nx, j = (c / nx) % ny, k = c / (nx * ny);
                const CellData d = cv.load<kSrc | kTauM>(i, j, k);
                const bool src = static_cast<double>(d.sfr) > sfr_min;
                const bool thick = src && rrk::source_is_thick(d.tau);
                rrk::CellStats cs;
                cnt(c) = rrk::step1_count(&rays(0) + start(c), static_cast<int>(start(c + 1) - start(c)), src, thick,
                                          ox + i, oy + j, oz + k, d.Luv, P, cs);
                if (cs.ncst) Kokkos::atomic_add(&stats(0), cs.ncst);
                if (cs.nspl) Kokkos::atomic_add(&stats(1), cs.nspl);
                if (cs.nmrg) Kokkos::atomic_add(&stats(2), cs.nmrg);
                if (cs.ndel) Kokkos::atomic_add(&stats(3), cs.ndel);
              });
          const int ncell = B.ncell;
          Kokkos::parallel_scan(
              "rt_scan_out", ncell, KOKKOS_LAMBDA(const int c, long long& upd, const bool fin) {
                if (fin) ostart(c) = upd;
                upd += cnt(c);
                if (fin && c == ncell - 1) ostart(ncell) = upd;
              });
          auto hs = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), Kokkos::subview(ostart, ncell));
          B.nout = hs();
          ensure(B.out, B.nout, "out");
          ensure(B.idx, B.nr, "idx");
          ensure(B.key, B.nr, "key");
          Kokkos::deep_copy(B.acc, 0);
          auto out = B.out;
          auto idx = B.idx;
          auto key = B.key;
          auto acc = B.acc;
          B.pmb->par_for(
              "rt_step1_write", 0, B.ncell - 1, KOKKOS_LAMBDA(const int c) {
                const int i = c % nx, j = (c / nx) % ny, k = c / (nx * ny);
                const CellData d = cv.load<kSrc | kTauM>(i, j, k);
                const bool src = static_cast<double>(d.sfr) > sfr_min;
                const bool thick = src && rrk::source_is_thick(d.tau);
                if (thick) {
                  double v[NDEP];
                  rrk::source_thick_deposit(d.Luv, P, v);
                  for (int m = 0; m < NDEP; ++m)
                    Kokkos::atomic_add(&acc(k + G, j + G, i + G, m), rrk::to_fx(v[m], rrk::dep_scale(m)));
                }
                if (cnt(c) == 0) return;
                const int w = rrk::step1_write(&rays(0) + start(c), static_cast<int>(start(c + 1) - start(c)), src,
                                               thick, ox + i, oy + j, oz + k, d.Luv, P, &idx(0) + start(c),
                                               &key(0) + start(c), &out(0) + ostart(c));
                if (w != cnt(c)) Kokkos::abort("step1_write count mismatch");
              });
        }
        Kokkos::fence();
        t1 = sp::now();
        T.step1 += t1 - t0;
        t0 = t1;
        // 1d'. cast every output ray once (deposits incl. ghost zone)
        for (auto& B : blks) {
          if (B.nout == 0) continue;
          auto out = B.out;
          auto acc = B.acc;
          auto tf = B.tf;
          const int nx = B.nx, ny = B.ny, nz = B.nz, ox = B.g.off[0], oy = B.g.off[1], oz = B.g.off[2];
          Kokkos::parallel_for(
              "rt_cast_1", B.nout, KOKKOS_LAMBDA(const long long r) {
                auto tau = [&](int i, int j, int k) -> const float* {
                  const int li = local_img(i, ox, nx, N), lj = local_img(j, oy, ny, N), lk = local_img(k, oz, nz, N);
                  if (li < 0 || lj < 0 || lk < 0) Kokkos::abort("cast_1 read outside the ghost zone");
                  return &tf(lk, lj, li, 0);
                };
                auto dep = [&](int i, int j, int k, const double* v) {
                  const int li = local_img(i, ox, nx, N), lj = local_img(j, oy, ny, N), lk = local_img(k, oz, nz, N);
                  if (li < 0 || lj < 0 || lk < 0) Kokkos::abort("cast_1 deposit outside the ghost zone");
                  for (int m = 0; m < NDEP; ++m)
                    Kokkos::atomic_add(&acc(lk, lj, li, m), rrk::to_fx(v[m], rrk::dep_scale(m)));
                };
                rrk::cast_1(out(r), P, tau, dep);
              });
          ray_substeps += B.nout;
        }
        Kokkos::fence();
        t1 = sp::now();
        T.cast1 += t1 - t0;
        t0 = t1;
        // deposits: ghost zones -> owners (reverse exchange)
        fold_deposits(blks, lid_of_gid, rank, comm, F);
        t1 = sp::now();
        T.fold += t1 - t0;
        t0 = t1;
        // 2. physics
        {
          rrk::PhysCtx X = run.phys();
          X.rates.cool = rt_cool.data();
          X.rates.alpha = rt_alpha.data();
          X.rates.gcol = rt_gcol.data();
          for (auto& B : blks) {
            const CellVars cv = B.cv;
            const int nx = B.nx, ny = B.ny;
            auto acc = B.acc;
            B.pmb->par_for(
                "rt_physics", 0, B.ncell - 1, KOKKOS_LAMBDA(const int c) {
                  const int i = c % nx, j = (c / nx) % ny, k = c / (nx * ny);
                  CellData d = cv.load<kGas | kGamLam | kFgam | kUW>(i, j, k);
                  const long long* a = &acc(k + G, j + G, i + G, 0);
                  for (int m = 0; m < 5; ++m) d.ngam[m] = static_cast<float>(rrk::from_fx(a[m], rrk::kScaleNgam));
                  for (int m = 0; m < 5; ++m) d.nu[m] = static_cast<float>(rrk::from_fx(a[5 + m], rrk::kScaleNu));
                  for (int m = 0; m < 3; ++m) d.dx[m] = static_cast<float>(rrk::from_fx(a[10 + m], rrk::kScaleNgam));
                  rrk::physics_cell(d, X);
                  cv.store<kGas | kDx | kGamLam | kNu | kNgam | kFgam | kUW>(i, j, k, d);
                });
          }
          Kokkos::fence();
        }
        t1 = sp::now();
        T.phys += t1 - t0;
        t0 = t1;
        sp::exchange_ghosts(pmesh, "fgam_only", {"fgam"});
        for (auto& B : blks) mirror_tf(B, false, true);
        Kokkos::fence();
        t1 = sp::now();
        T.halo += t1 - t0;
        t0 = t1;
        // 3. advance
        for (auto& B : blks) {
          ensure(B.keep, B.nout, "keep");
          ensure(B.kpos, B.nout, "kpos");
          if (B.nout == 0) continue;
          auto out = B.out;
          auto tf = B.tf;
          auto keep = B.keep;
          auto stats = B.stats;
          const int nx = B.nx, ny = B.ny, nz = B.nz, ox = B.g.off[0], oy = B.g.off[1], oz = B.g.off[2];
          Kokkos::parallel_for(
              "rt_cast_3", B.nout, KOKKOS_LAMBDA(const long long r) {
                auto tfn = [&](int i, int j, int k, const float*& tau, const float*& fgam) {
                  const int li = local_img(i, ox, nx, N), lj = local_img(j, oy, ny, N), lk = local_img(k, oz, nz, N);
                  if (li < 0 || lj < 0 || lk < 0) Kokkos::abort("cast_3 read outside the ghost zone");
                  tau = &tf(lk, lj, li, 0);
                  fgam = &tf(lk, lj, li, 5);
                };
                keep(r) = rrk::cast_3(out(r), P, tfn) ? 1 : 0;
                if (!keep(r)) Kokkos::atomic_add(&stats(4), 1LL);
              });
        }
        Kokkos::fence();
        t1 = sp::now();
        T.cast3 += t1 - t0;
        t0 = t1;
        // surviving rays -> swarm (position = centre of the home cell image
        // nearest the block)
        for (auto& B : blks) {
          if (B.nout == 0) continue;
          auto keep = B.keep, kpos = B.kpos;
          int nkeep = 0;
          Kokkos::parallel_scan(
              "rt_scan_keep", static_cast<int>(B.nout),
              KOKKOS_LAMBDA(const int r, int& upd, const bool fin) {
                if (fin) kpos(r) = upd;
                upd += keep(r);
              },
              nkeep);
          auto swarm = get_swarm(B.pmb);
          auto ctx = swarm->AddEmptyParticles(nkeep);
          auto ri = swarm->Get<int>("ri").Get();
          auto rf = swarm->Get<Real>("rf").Get();
          auto& xs = swarm->Get<Real>(swarm_position::x::name()).Get();
          auto& ys = swarm->Get<Real>(swarm_position::y::name()).Get();
          auto& zs = swarm->Get<Real>(swarm_position::z::name()).Get();
          auto out = B.out;
          const int nx = B.nx, ny = B.ny, nz = B.nz, ox = B.g.off[0], oy = B.g.off[1], oz = B.g.off[2];
          Kokkos::parallel_for(
              "rt_to_swarm", static_cast<int>(B.nout), KOKKOS_LAMBDA(const int r) {
                if (!keep(r)) return;
                const int n = ctx.GetNewParticleIndex(kpos(r));
                const RayP& ry = out(r);
                store_ray(ri, rf, n, ry);
                const int li = local_img(rrk::cell_of(ry.x[0], N), ox, nx, N);
                const int lj = local_img(rrk::cell_of(ry.x[1], N), oy, ny, N);
                const int lk = local_img(rrk::cell_of(ry.x[2], N), oz, nz, N);
                if (li < 0 || lj < 0 || lk < 0) Kokkos::abort("ray moved beyond the ghost zone");
                xs(n) = ox + (li - G) + 0.5;
                ys(n) = oy + (lj - G) + 0.5;
                zs(n) = oz + (lk - G) + 0.5;
              });
        }
        Kokkos::fence();
        t1 = sp::now();
        T.sort += t1 - t0;
        t0 = t1;
        migrate_rays(pmesh);
        Kokkos::fence();
        t1 = sp::now();
        T.migrate += t1 - t0;
        // substep statistics (integers)
        long long sv[5] = {0, 0, 0, 0, 0};
        for (auto& B : blks) {
          auto h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), B.stats);
          for (int q = 0; q < 5; ++q) sv[q] += h(q);
        }
        MPI_Allreduce(MPI_IN_PLACE, sv, 5, MPI_LONG_LONG, MPI_SUM, comm);
        st.ncst1 = sv[0];
        st.nspl1 = sv[1];
        st.nmrg1 = sv[2];
        st.ndel1 = sv[3] + sv[4];
        run.end_substep();
      }
      // post_rt: redshift rays, final gas/tau
      const double t7 = sp::now();
      {
        const rrk::RedshiftTable RT = run.redshift_table();
        const radray::Units u = run.units;
        long long tot[2] = {0, 0};
        for (auto& B : blks) {
          auto swarm = get_swarm(B.pmb);
          const int nmax = swarm->GetMaxActiveIndex();
          if (nmax >= 0) {
            auto ri = swarm->Get<int>("ri").Get();
            auto rf = swarm->Get<Real>("rf").Get();
            auto sd = swarm->GetDeviceContext();
            long long bt = 0;
            Kokkos::parallel_reduce(
                "rt_redshift", nmax + 1,
                KOKKOS_LAMBDA(const int n, long long& acc_) {
                  if (!sd.IsActive(n)) return;
                  RayP r = load_ray(ri, rf, n);
                  rrk::redshift_ray(r, RT);
                  store_ray(ri, rf, n, r);
                  for (int m = 0; m < 5; ++m) acc_ += rrk::to_fx(static_cast<double>(r.ngam[m]), rrk::kScaleTot);
                },
                bt);
            tot[0] += bt;
          }
          tot[1] += swarm->GetNumActive();
          const CellVars cv = B.cv;
          const int nx = B.nx, ny = B.ny;
          B.pmb->par_for(
              "rt_post_cell", 0, B.ncell - 1, KOKKOS_LAMBDA(const int c) {
                const int i = c % nx, j = (c / nx) % ny, k = c / (nx * ny);
                CellData d = cv.load<kUW>(i, j, k);
                rrk::post_cell(d, u, sig.s);
                cv.store<kGas | kTauM>(i, j, k, d);
              });
        }
        Kokkos::fence();
        MPI_Allreduce(MPI_IN_PLACE, tot, 2, MPI_LONG_LONG, MPI_SUM, comm);
        st.ngam_tot_fx = tot[0];
        st.nray_active = tot[1];
      }
      T.post += sp::now() - t7;
      if (rank == 0)
        std::printf(" step %3d  z = %.5f  nstep_rt = %d  rays = %lld  cast = %lld  merge = %lld  split = %lld\n",
                    run.clk.istep, run.clk.redshift, nsub, st.nray_active, st.ncst1, st.nmrg1, st.nspl1);
      if (output) {
        std::vector<long long> gidx;
        std::vector<char> rec;
        for (auto& B : blks) {
          auto hv = make_cv(B.cv.rtc.GetHostMirrorAndCopy(), B.cv.tau.GetHostMirrorAndCopy(),
                            B.cv.fgam.GetHostMirrorAndCopy(), B.g.is, B.g.js, B.g.ks);
          for (int k = 0; k < B.nz; ++k)
            for (int j = 0; j < B.ny; ++j)
              for (int i = 0; i < B.nx; ++i) {
                gidx.push_back((B.g.off[0] + i) + static_cast<long long>(N) * ((B.g.off[1] + j) + static_cast<long long>(N) * (B.g.off[2] + k)));
                const CellData d = hv.load(i, j, k);
                const char* p = reinterpret_cast<const char*>(&d);
                rec.insert(rec.end(), p, p + sizeof d);
              }
        }
        std::vector<char> all = sp::gather_cells(gidx, rec, sizeof(CellData), ncell_glob, comm);
        if (rank == 0) {
          std::vector<CellData> cells(ncell_glob);
          std::memcpy(cells.data(), all.data(), all.size());
          run.write_step(odir, cells, st);
        }
      }
    }
    MPI_Barrier(comm);
    const double tt = sp::now() - t_start;
    MPI_Allreduce(MPI_IN_PLACE, &ray_substeps, 1, MPI_LONG_LONG, MPI_SUM, comm);
    // per-phase timers: max over ranks
    double tv[11] = {T.cell, T.halo, T.hist, T.step1, T.cast1, T.fold, T.phys, T.cast3, T.sort, T.migrate, T.post};
    MPI_Allreduce(MPI_IN_PLACE, tv, 11, MPI_DOUBLE, MPI_MAX, comm);
    if (rank == 0) {
      std::printf(" ranks %d  threads %d  blocks %d (mb %d)  wall %.3f s  ray-substeps %lld  "
                  "ns/ray-substep*threads*ranks %.1f\n",
                  nranks, nthreads, static_cast<int>(pmesh->nbtotal), mb, tt, ray_substeps,
                  1e9 * tt * nthreads * nranks / (ray_substeps > 0 ? ray_substeps : 1));
      std::printf(" timers (max over ranks): cell %.3f halo %.3f hist %.3f step1 %.3f cast1 %.3f fold %.3f "
                  "phys %.3f cast3 %.3f sort %.3f migrate %.3f post %.3f\n",
                  tv[0], tv[1], tv[2], tv[3], tv[4], tv[5], tv[6], tv[7], tv[8], tv[9], tv[10]);
    }
    MPI_Comm_free(&comm);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "rt_parthenon: error: %s\n", e.what());
    status = 1;
  }
  pman.ParthenonFinalize();
  return status;
}
