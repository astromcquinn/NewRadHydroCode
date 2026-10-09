// Complete binary checkpoint (see checkpoint.hpp).
#include "radray/checkpoint.hpp"

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
constexpr std::int32_t kVersion = 2;

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
    o.pod(d);
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
  if (ver != 1 && ver != kVersion)
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
  in.f.read(m, 8);
  if (!in.f || std::memcmp(m, kMagic, 8) != 0)
    throw std::runtime_error("checkpoint trailer missing: " + path);
}

}  // namespace radray
