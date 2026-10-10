// Complete binary checkpoint (see checkpoint.hpp).
#include "radray/checkpoint.hpp"

#include "radray/recomb_rad.hpp"
#include "radray/rt.hpp"
#include "radray/rt_backend.hpp"
#include "radray/rt_moments.hpp"

#include <cstdio>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <type_traits>

namespace radray {

namespace {

constexpr char kMagic[8] = {'R', 'A', 'D', 'R', 'A', 'Y', 'C', 'K'};
// Version 2 (hydro_scheme = unsplit, docs/hydro_unsplit.md): the internal
// energy density gas.e (empty with the legacy scheme) and the unsplit
// counters after the hydro scalars. Version 1 files are still read (e empty).
// Version 3 (rt_method = m1 | otvet, docs/rt_moments.md): the moment RT
// state before the trailer (empty with the rays). Versions 1-2 are read as
// ray runs.
// Version 4 (recomb_rad = 1, docs/recomb_radiation.md): the recombination
// photon buffers, on-the-spot counts, cell modes and budget before the
// trailer (empty when off).
// Version 5 (Stage D, docs/stage_d.md): the ray storage layout is
// rrk::RayP + links (the links moved from between ipix and r to the end;
// versions <= 4 are converted on reading), and the AMReX backend's section
// (rt_backend, then its rays gathered on rank 0) before the trailer.
constexpr std::int32_t kVersion = 5;

// The ray layout of versions <= 4
struct RayV4 {
  std::int32_t nside, ipix;
  std::int32_t next, prev;
  float r, dr;
  float x[3], dx[3];
  float nu[kNfreq], ngam[kNfreq];
};
static_assert(sizeof(RayV4) == sizeof(Ray), "RayV4 layout");

struct Out {
  std::ofstream f;
  template <class T>
  void pod(const T& x) {
    static_assert(std::is_trivially_copyable<T>::value, "POD only");
    f.write(reinterpret_cast<const char*>(&x), sizeof x);
  }
  template <class T>
  void vec(const std::vector<T>& v) {
    const std::int64_t n = static_cast<std::int64_t>(v.size());
    pod(n);
    if (n) f.write(reinterpret_cast<const char*>(v.data()), static_cast<std::streamsize>(n * sizeof(T)));
  }
};

struct In {
  std::ifstream f;
  std::string path;
  template <class T>
  void pod(T& x) {
    static_assert(std::is_trivially_copyable<T>::value, "POD only");
    f.read(reinterpret_cast<char*>(&x), sizeof x);
    if (!f) throw std::runtime_error("checkpoint truncated: " + path);
  }
  template <class T>
  void vec(std::vector<T>& v, std::int64_t expect = -1) {
    std::int64_t n = 0;
    pod(n);
    if (expect >= 0 && n != expect)
      throw std::runtime_error("checkpoint does not match this run (array size " +
                               std::to_string(n) + " != " + std::to_string(expect) + "): " + path);
    v.resize(static_cast<std::size_t>(n));
    if (n) f.read(reinterpret_cast<char*>(v.data()), static_cast<std::streamsize>(n * sizeof(T)));
    if (!f) throw std::runtime_error("checkpoint truncated: " + path);
  }
};

template <class IO, class RS>
void ray_stat_io(IO& io, RS& r) {
  io.pod(r.ncst1);
  io.pod(r.ncst2);
  io.pod(r.nmrg1);
  io.pod(r.nmrg2);
  io.pod(r.nspl1);
  io.pod(r.nspl2);
  io.pod(r.ndel1);
  io.pod(r.nfree);
  io.vec(r.ifree);
}

}  // namespace

void write_checkpoint(const std::string& path, const Simulation& s, const HydroState& hs,
                      const pm::PMState& pm, const DriverState& d) {
  // rt_backend = amrex: the rays live on the level; gather them on rank 0
  // (collective: every rank calls this function, only io_rank() writes)
  std::vector<rrk::RayP> brays;
  if (use_amrex(s.cfg)) amrex_rt_gather_rays(s, brays);
  if (!io_rank()) return;
  const std::string tmp = path + ".tmp";
  {
    Out o;
    o.f.open(tmp, std::ios::binary);
    if (!o.f) throw std::runtime_error("cannot write " + tmp);
    o.f.write(kMagic, 8);
    o.pod(kVersion);
    o.pod(s.cfg.Ngas);
    o.pod(s.cfg.Ndm);
    o.pod(s.cfg.Nray_max);
    o.pod(s.clk);
    // RT grid
    const RTGrid& G = s.rt;
    for (const auto* v : {&G.rhog, &G.fHI, &G.fHeI, &G.fHeII, &G.temp, &G.Luv, &G.sfr, &G.dx,
                          &G.gamma, &G.lambda, &G.nu, &G.ngam, &G.fgam, &G.tau})
      o.vec(*v);
    o.vec(G.htoc);
    // gas, grid velocity, hydro scalars
    o.vec(s.gas.u);
    o.vec(s.gas.w);
    o.vec(hs.v);
    o.pod(hs.c_limit);
    o.pod(hs.dt_hy);
    o.pod(hs.vmax_last);
    o.pod(hs.counters);
    o.vec(s.gas.e);
    o.pod(hs.ucounters);
    // rays: storage up to the high-water mark, pool and thresholds
    const RayState& R = s.rays;
    o.pod(R.Nray);
    o.pod(R.Nrbuf);
    const std::int64_t hwm = R.hwm();
    o.pod(hwm);
    if (hwm > 0)
      o.f.write(reinterpret_cast<const char*>(&R.ray[1]),
                static_cast<std::streamsize>(hwm * sizeof(Ray)));
    o.vec(R.ray_index_);
    ray_stat_io(o, const_cast<RayStat&>(R.stat));
    ray_stat_io(o, const_cast<RayStat&>(R.save));
    o.pod(R.ngam_ray);
    o.pod(R.npool_limited);
    o.pod(R.npool_forced);
    // deterministic_parallel pool: all free slots in one list
    std::vector<std::int32_t> freel = R.pool_free;
    for (const ThreadRays& t : R.tl) freel.insert(freel.end(), t.free.begin(), t.free.end());
    o.vec(freel);
    o.pod(R.next_fresh);
    // RT statistics, sources
    o.pod(s.rt_stat);
    o.pod(s.rt_save);
    o.pod(s.Nsource);
    // particles
    o.vec(pm.xv);
    o.pod(pm.dm_cpu);
    // driver
    // Copy d into a zeroed struct first: DriverState has padding bytes, and a
    // raw write of d would copy whatever the stack held there, so identical
    // runs would give checkpoints that differ in those bytes.
    DriverState dz;
    std::memset(&dz, 0, sizeof dz);
    dz.idata1 = d.idata1;
    dz.idata2 = d.idata2;
    dz.data1_out = d.data1_out;
    dz.data2_out = d.data2_out;
    dz.Nhalo = d.Nhalo;
    dz.movie_frame = d.movie_frame;
    dz.movie_next = d.movie_next;
    o.pod(dz);
    // moment RT (version 3)
    const std::int32_t method = static_cast<std::int32_t>(s.cfg.rt_method);
    o.pod(method);
    const std::vector<float> none;
    const MomentState* M = s.moments.get();
    for (int q = 0; q < 5; ++q) {
      const std::vector<float>* v[5] = {&none, &none, &none, &none, &none};
      if (M) {
        v[0] = &M->n;
        v[1] = &M->fx;
        v[2] = &M->fy;
        v[3] = &M->fz;
        v[4] = &M->e;
      }
      o.vec(*v[q]);
    }
    const double budget[3] = {M ? M->injected : 0, M ? M->absorbed : 0, M ? M->lost : 0};
    o.pod(budget);
    // recombination radiation (version 4)
    const RecombState* Rc = s.recomb.get();
    const std::vector<std::uint8_t> none8;
    o.vec(Rc ? Rc->ngam : none);
    o.vec(Rc ? Rc->egam : none);
    o.vec(Rc ? Rc->ots : none);
    o.vec(Rc ? Rc->mode : none8);
    const double rbud[3] = {Rc ? Rc->launched : 0, Rc ? Rc->dropped : 0, Rc ? Rc->onspot : 0};
    o.pod(rbud);
    // AMReX backend (version 5): its rays in the float layout
    const std::int32_t backend = static_cast<std::int32_t>(s.cfg.rt_backend);
    o.pod(backend);
    o.vec(brays);
    o.pod(s.ngam_total_backend);
    o.f.write(kMagic, 8);
    if (!o.f) throw std::runtime_error("error writing " + tmp);
  }
  if (std::rename(tmp.c_str(), path.c_str()) != 0)
    throw std::runtime_error("cannot rename " + tmp + " to " + path);
}

void read_checkpoint(const std::string& path, Simulation& s, HydroState& hs, pm::PMState& pm,
                     DriverState& d) {
  In in;
  in.path = path;
  in.f.open(path, std::ios::binary);
  if (!in.f) throw std::runtime_error("cannot open checkpoint " + path);
  char m[8];
  in.f.read(m, 8);
  if (!in.f || std::memcmp(m, kMagic, 8) != 0) throw std::runtime_error("not a checkpoint: " + path);
  std::int32_t ver = 0, Ngas = 0, Ndm = 0;
  std::int64_t nrmax = 0;
  in.pod(ver);
  in.pod(Ngas);
  in.pod(Ndm);
  in.pod(nrmax);
  if (ver < 1 || ver > kVersion)
    throw std::runtime_error("checkpoint version mismatch: " + path);
  if (Ngas != s.cfg.Ngas || Ndm != s.cfg.Ndm || nrmax != s.cfg.Nray_max)
    throw std::runtime_error("checkpoint is for N = " + std::to_string(Ngas) + ", Ndm = " +
                             std::to_string(Ndm) + ", Nray_max = " + std::to_string(nrmax) +
                             "; this run has different sizes");
  in.pod(s.clk);
  RTGrid& G = s.rt;
  for (auto* v : {&G.rhog, &G.fHI, &G.fHeI, &G.fHeII, &G.temp, &G.Luv, &G.sfr, &G.dx, &G.gamma,
                  &G.lambda, &G.nu, &G.ngam, &G.fgam, &G.tau})
    in.vec(*v, static_cast<std::int64_t>(v->size()));
  in.vec(G.htoc, static_cast<std::int64_t>(G.htoc.size()));
  in.vec(s.gas.u, static_cast<std::int64_t>(s.gas.u.size()));
  in.vec(s.gas.w, static_cast<std::int64_t>(s.gas.w.size()));
  in.vec(hs.v, static_cast<std::int64_t>(hs.v.size()));
  in.pod(hs.c_limit);
  in.pod(hs.dt_hy);
  in.pod(hs.vmax_last);
  in.pod(hs.counters);
  if (ver >= 2) {
    in.vec(s.gas.e);
    if (!s.gas.e.empty() && static_cast<std::int64_t>(s.gas.e.size()) != s.gas.ncell)
      throw std::runtime_error("checkpoint: internal-energy field of the wrong size: " + path);
    in.pod(hs.ucounters);
  } else {
    s.gas.e.clear();
    hs.ucounters = UnsplitCounters{};
  }
  RayState& R = s.rays;
  in.pod(R.Nray);
  in.pod(R.Nrbuf);
  std::int64_t hwm = 0;
  in.pod(hwm);
  if (hwm > R.Nray_max) throw std::runtime_error("checkpoint: ray high-water mark > Nray_max");
  if (hwm > 0) {
    in.f.read(reinterpret_cast<char*>(&R.ray[1]), static_cast<std::streamsize>(hwm * sizeof(Ray)));
    if (!in.f) throw std::runtime_error("checkpoint truncated: " + path);
    if (ver <= 4) {
      // the links sat between ipix and r
      for (std::int64_t ir = 1; ir <= hwm; ++ir) {
        RayV4 o;
        std::memcpy(&o, &R.ray[ir], sizeof o);
        Ray& r = R.ray[ir];
        r.nside = o.nside;
        r.ipix = o.ipix;
        r.next = o.next;
        r.prev = o.prev;
        r.r = o.r;
        r.dr = o.dr;
        for (int d = 0; d < 3; ++d) {
          r.x[d] = o.x[d];
          r.dx[d] = o.dx[d];
        }
        for (int n = 0; n < kNfreq; ++n) {
          r.nu[n] = o.nu[n];
          r.ngam[n] = o.ngam[n];
        }
      }
    }
  }
  R.det_set_hwm(hwm);
  in.vec(R.ray_index_);
  ray_stat_io(in, R.stat);
  ray_stat_io(in, R.save);
  in.pod(R.ngam_ray);
  in.pod(R.npool_limited);
  in.pod(R.npool_forced);
  in.vec(R.pool_free);
  for (ThreadRays& t : R.tl) t.free.clear();
  in.pod(R.next_fresh);
  in.pod(s.rt_stat);
  in.pod(s.rt_save);
  in.pod(s.Nsource);
  in.vec(pm.xv, static_cast<std::int64_t>(pm.xv.size()));
  in.pod(pm.dm_cpu);
  in.pod(d);
  std::int32_t method = 0;
  if (ver >= 3) in.pod(method);
  if (method != static_cast<std::int32_t>(s.cfg.rt_method))
    throw std::runtime_error(std::string("checkpoint was written with rt_method = ") +
                             rt_method_name(static_cast<RTMethod>(method)) +
                             ", this run has rt_method = " + rt_method_name(s.cfg.rt_method));
  if (ver >= 3) {
    std::vector<float> v[5];
    for (int q = 0; q < 5; ++q) in.vec(v[q]);
    double budget[3];
    in.pod(budget);
    if (use_moments(s.cfg)) {
      moment_init(s);
      MomentState& M = *s.moments;
      const std::int64_t n = kNfreq * M.ncell;
      std::vector<float>* dst[5] = {&M.n, &M.fx, &M.fy, &M.fz, &M.e};
      for (int q = 0; q < 5; ++q) {
        if (static_cast<std::int64_t>(v[q].size()) != n)
          throw std::runtime_error("checkpoint: moment RT state of the wrong size: " + path);
        *dst[q] = std::move(v[q]);
      }
      M.injected = budget[0];
      M.absorbed = budget[1];
      M.lost = budget[2];
    }
  }
  if (ver >= 4) {
    std::vector<float> ng, eg, ots;
    std::vector<std::uint8_t> mode;
    in.vec(ng);
    in.vec(eg);
    in.vec(ots);
    in.vec(mode);
    double rbud[3];
    in.pod(rbud);
    if (s.cfg.recomb_rad != !ng.empty())
      throw std::runtime_error("checkpoint and run differ in recomb_rad: " + path);
    if (s.cfg.recomb_rad) {
      recomb_init(s);
      RecombState& Rc = *s.recomb;
      if (ng.size() != Rc.ngam.size() || eg.size() != Rc.egam.size() ||
          ots.size() != Rc.ots.size() || mode.size() != Rc.mode.size())
        throw std::runtime_error("checkpoint: recombination buffers of the wrong size: " + path);
      Rc.ngam = std::move(ng);
      Rc.egam = std::move(eg);
      Rc.ots = std::move(ots);
      Rc.mode = std::move(mode);
      Rc.launched = rbud[0];
      Rc.dropped = rbud[1];
      Rc.onspot = rbud[2];
    }
  } else if (s.cfg.recomb_rad) {
    throw std::runtime_error("checkpoint has no recombination radiation state: " + path);
  }
  std::int32_t backend = 0;
  std::vector<rrk::RayP> brays;
  if (ver >= 5) {
    in.pod(backend);
    in.vec(brays);
    in.pod(s.ngam_total_backend);
  }
  if (backend == 0 && use_amrex(s.cfg)) {
    // A serial (legacy_rt = 0) checkpoint restarted on the backend: the
    // live rays of the serial chains go to the level (the timing and
    // scaling runs restart production checkpoints this way; the result is
    // the backend's, not the serial continuation). The serial pool is then
    // released.
    if (s.cfg.legacy_rt) throw std::runtime_error("rt_backend = amrex needs a legacy_rt = 0 checkpoint");
    brays.clear();
    for (std::int64_t c = 0; c < G.ncell; ++c)
      for (std::int64_t ir = G.htoc_ray(3, c); ir > 0; ir = R.at(ir).next) brays.push_back(R.at(ir));
    s.ngam_total_backend = ngam_rays_total_serial(s);
    R.init(R.Nray_max, R.Nfree_max);
    std::fill(G.htoc.begin(), G.htoc.end(), 0);
    R.Nray = static_cast<std::int64_t>(brays.size());
  } else if (backend != static_cast<std::int32_t>(s.cfg.rt_backend)) {
    throw std::runtime_error("checkpoint was written with rt_backend = " +
                             std::string(backend ? "amrex" : "serial") + ", this run has " +
                             rt_backend_summary(s.cfg));
  }
  if (use_amrex(s.cfg)) amrex_rt_scatter_rays(s, brays);
  in.f.read(m, 8);
  if (!in.f || std::memcmp(m, kMagic, 8) != 0)
    throw std::runtime_error("checkpoint trailer missing: " + path);
}

}  // namespace radray
