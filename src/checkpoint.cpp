// Complete binary checkpoint (see checkpoint.hpp).
#include "radray/checkpoint.hpp"

#include "radray/amrex_state.hpp"
#include "radray/recomb_rad.hpp"
#include "radray/rt.hpp"
#include "radray/rt_backend.hpp"
#include "radray/rt_moments.hpp"

#include <algorithm>
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
// Version 6 (docs/rt_multicell.md): rt_cells_per_substep before the
// trailer (versions <= 5: 1). The substep carry rt_carry of the clock is in
// substeps of that many cells; a restart with another value converts it.
// Version 7 (fof_dt_myr, docs/pm_refine.md): the on-the-fly FoF cache
// before the trailer (empty when off). Versions <= 6 are read as "no cache".
constexpr std::int32_t kVersion = 7;

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
  // skips a vector written by Out::vec
  template <class T>
  void skip_vec() {
    std::int64_t n = 0;
    pod(n);
    f.seekg(static_cast<std::streamoff>(n * static_cast<std::int64_t>(sizeof(T))), std::ios::cur);
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
  if (use_amrex(s.cfg) && use_moments(s.cfg))
    throw std::runtime_error("checkpoint: rt_method = m1 | otvet with rt_backend = amrex is not implemented yet");
  // rt_backend = amrex: the rays live on the level; rank 0 receives them in
  // chunks while it writes (collective: every rank takes part in the count
  // and the gather, only io_rank() writes; rank 0's transient is one chunk)
  const std::int64_t nbrays = use_amrex(s.cfg) ? amrex_rt_nrays(s) : 0;
  if (!io_rank()) {
    if (use_amrex(s.cfg)) amrex_rt_gather_rays_chunks(s, [](const std::vector<rrk::RayP>&) {});
    return;
  }
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
    o.pod(nbrays);  // the vec layout: the count, then the rays, streamed
    if (use_amrex(s.cfg))
      amrex_rt_gather_rays_chunks(s, [&](const std::vector<rrk::RayP>& chunk) {
        o.f.write(reinterpret_cast<const char*>(chunk.data()),
                  static_cast<std::streamsize>(chunk.size() * sizeof(rrk::RayP)));
      });
    o.pod(s.ngam_total_backend);
    // rt_cells_per_substep (version 6)
    const std::int32_t cells = s.cfg.rt_cells_per_substep;
    o.pod(cells);
    // FoF cache (version 7)
    const FofCache& fc = s.fof;
    const std::int32_t fvalid = fc.valid ? 1 : 0;
    o.pod(fvalid);
    const double fd[3] = {fc.year, fc.vel_unit, fc.dtau};
    o.pod(fd);
    const std::int64_t fn[2] = {fc.nfof, fc.nreuse};
    o.pod(fn);
    o.vec(fc.halos);
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
  if (use_amrex(s.cfg) && use_moments(s.cfg))
    throw std::runtime_error("checkpoint: rt_method = m1 | otvet with rt_backend = amrex is not implemented yet");
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
  // amrex_state = 1: the grids are rank 0's mirror (scattered onto the
  // levels by the driver afterwards); the other ranks skip their blocks,
  // Luv and sfr too (their sources come from the next calc_halo)
  const bool skip_grids = use_amrex_state(s.cfg) && !io_rank();
  RTGrid& G = s.rt;
  for (auto* v : {&G.rhog, &G.fHI, &G.fHeI, &G.fHeII, &G.temp, &G.Luv, &G.sfr, &G.dx, &G.gamma,
                  &G.lambda, &G.nu, &G.ngam, &G.fgam, &G.tau}) {
    if (skip_grids)
      in.skip_vec<float>();
    else
      in.vec(*v, static_cast<std::int64_t>(v->size()));
  }
  if (skip_grids) {
    in.skip_vec<std::int64_t>();
    in.skip_vec<float>();
    in.skip_vec<float>();
    in.skip_vec<float>();
  } else {
    in.vec(G.htoc, static_cast<std::int64_t>(G.htoc.size()));
    in.vec(s.gas.u, static_cast<std::int64_t>(s.gas.u.size()));
    in.vec(s.gas.w, static_cast<std::int64_t>(s.gas.w.size()));
    in.vec(hs.v, static_cast<std::int64_t>(hs.v.size()));
  }
  in.pod(hs.c_limit);
  in.pod(hs.dt_hy);
  in.pod(hs.vmax_last);
  in.pod(hs.counters);
  if (ver >= 2) {
    if (skip_grids) {
      in.skip_vec<float>();
      s.gas.e_on = true;
    } else {
      in.vec(s.gas.e);
      if (!s.gas.e.empty() && static_cast<std::int64_t>(s.gas.e.size()) != s.gas.ncell)
        throw std::runtime_error("checkpoint: internal-energy field of the wrong size: " + path);
    }
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
  // rt_backend = amrex: only rank 0 needs the serial rays (they go to the
  // level through rank 0); the other ranks skip them
  const bool skip_rays = use_amrex(s.cfg) && !io_rank();
  if (hwm > 0 && skip_rays) {
    in.f.seekg(static_cast<std::streamoff>(hwm * static_cast<std::int64_t>(sizeof(Ray))), std::ios::cur);
    if (!in.f) throw std::runtime_error("checkpoint truncated: " + path);
  } else if (hwm > 0) {
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
  if (skip_grids)
    in.skip_vec<float>();  // the particles are rank 0's mirror too
  else
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
    if (!s.cfg.recomb_rad && !ng.empty())
      throw std::runtime_error("checkpoint and run differ in recomb_rad: " + path);
    // recomb_rad = 1 with empty buffers: written before the RT started
    // (prep_rt allocates them), so the state is left unallocated as in a
    // fresh run
    if (s.cfg.recomb_rad && !ng.empty()) {
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
  // rt_backend = amrex: the rays go onto the level in chunks as rank 0
  // reads them (collective: the same number of scatters on every rank)
  constexpr std::int64_t kRayChunk = std::int64_t(1) << 25;  // 32 M rays, 2.7 GB
  bool scattered = false;
  if (ver >= 5) {
    in.pod(backend);
    if (use_amrex(s.cfg) && backend != 0) {
      std::int64_t n = 0;
      if (!skip_rays)
        in.pod(n);
      else
        in.skip_vec<rrk::RayP>();  // past the block; rank 0 streams it
      double nd = static_cast<double>(n);
      backend_bcast(&nd, 1);
      n = static_cast<std::int64_t>(nd);
      const std::int64_t nchunk = std::max<std::int64_t>(1, (n + kRayChunk - 1) / kRayChunk);
      for (std::int64_t c = 0; c < nchunk; ++c) {
        brays.clear();
        if (!skip_rays) {
          const std::int64_t m = std::min(kRayChunk, n - c * kRayChunk);
          brays.resize(static_cast<std::size_t>(std::max<std::int64_t>(m, 0)));
          if (m > 0) {
            in.f.read(reinterpret_cast<char*>(brays.data()), static_cast<std::streamsize>(m * sizeof(rrk::RayP)));
            if (!in.f) throw std::runtime_error("checkpoint truncated: " + path);
          }
        }
        amrex_rt_scatter_rays(s, brays, c == 0);
      }
      brays.clear();
      brays.shrink_to_fit();
      scattered = true;
    } else if (skip_rays) {
      in.skip_vec<rrk::RayP>();
    } else {
      in.vec(brays);
    }
    in.pod(s.ngam_total_backend);
  }
  if (backend == 0 && use_amrex(s.cfg)) {
    // A serial (legacy_rt = 0) checkpoint restarted on the backend: the
    // live rays of the serial chains go to the level (the timing and
    // scaling runs restart production checkpoints this way; the result is
    // the backend's, not the serial continuation). The serial pool is then
    // released.
    if (s.cfg.legacy_rt) throw std::runtime_error("rt_backend = amrex needs a legacy_rt = 0 checkpoint");
    // the live rays of the chains, streamed onto the level in chunks (rank
    // 0 holds the serial pool and one chunk, not a second copy of the pool)
    double tot = 0;
    std::int64_t nlive = 0;
    if (!skip_rays) {
      for (std::int64_t c = 0; c < G.ncell; ++c)
        for (std::int64_t ir = G.htoc_ray(3, c); ir > 0; ir = R.at(ir).next) ++nlive;
      tot = ngam_rays_total_serial(s);
    }
    double nd = static_cast<double>(nlive);
    backend_bcast(&nd, 1);
    nlive = static_cast<std::int64_t>(nd);
    const std::int64_t nchunk = std::max<std::int64_t>(1, (nlive + kRayChunk - 1) / kRayChunk);
    std::int64_t c = 0, ir = 0;  // the traversal state on rank 0
    if (!skip_rays && G.ncell > 0) ir = G.htoc_ray(3, 0);
    for (std::int64_t k = 0; k < nchunk; ++k) {
      brays.clear();
      if (!skip_rays) {
        while (static_cast<std::int64_t>(brays.size()) < kRayChunk && c < G.ncell) {
          if (ir > 0) {
            brays.push_back(R.at(ir));
            ir = R.at(ir).next;
          } else if (++c < G.ncell) {
            ir = G.htoc_ray(3, c);
          }
        }
      }
      amrex_rt_scatter_rays(s, brays, k == 0);
    }
    brays.clear();
    brays.shrink_to_fit();
    scattered = true;
    // (rank 0 writes the outputs and checkpoints; the other ranks keep 0
    // until the next post_rt sets the reduced total on every rank)
    s.ngam_total_backend = tot;
    R.init(R.Nray_max, R.Nfree_max);
    std::fill(G.htoc.begin(), G.htoc.end(), 0);
    R.Nray = nlive;
  } else if (backend != static_cast<std::int32_t>(s.cfg.rt_backend)) {
    throw std::runtime_error("checkpoint was written with rt_backend = " +
                             std::string(backend ? "amrex" : "serial") + ", this run has " +
                             rt_backend_summary(s.cfg));
  }
  if (use_amrex(s.cfg) && !scattered) amrex_rt_scatter_rays(s, brays, true);
  std::int32_t cells = 1;
  if (ver >= 6) in.pod(cells);
  if (cells < 1) throw std::runtime_error("checkpoint: bad rt_cells_per_substep: " + path);
  if (cells != s.cfg.rt_cells_per_substep) {
    // rt_carry counts light substeps of `cells` RT cells; this run's
    // substeps cross rt_cells_per_substep cells
    s.clk.rt_carry = s.clk.rt_carry * cells / s.cfg.rt_cells_per_substep;
    if (io_rank())
      std::fprintf(stderr,
                   " read_checkpoint: rt_cells_per_substep changes from %d to %d at this restart "
                   "(substep carry converted)\n",
                   static_cast<int>(cells), s.cfg.rt_cells_per_substep);
  }
  s.fof = FofCache{};
  if (ver >= 7) {
    std::int32_t fvalid = 0;
    double fd[3];
    std::int64_t fn[2];
    in.pod(fvalid);
    in.pod(fd);
    in.pod(fn);
    in.vec(s.fof.halos);
    s.fof.valid = fvalid != 0;
    s.fof.year = fd[0];
    s.fof.vel_unit = fd[1];
    s.fof.dtau = fd[2];
    s.fof.nfof = fn[0];
    s.fof.nreuse = fn[1];
  }
  in.f.read(m, 8);
  if (!in.f || std::memcmp(m, kMagic, 8) != 0)
    throw std::runtime_error("checkpoint trailer missing: " + path);
}

}  // namespace radray
