// Optional ray snapshots for the ray-compression study
// (include/radray/ray_dump.hpp, docs/ray_compression.md). Off by default;
// params.txt ray_dump_dir with ray_dump_every and/or ray_dump_at turn it on
// in rt_driver and radhydro (ray_dump_enable, before the first RT step).
// radiative_transfer (src/rt_tools.cpp) then advances the global substep
// counter before rt_step_1 and calls ray_dump_after_substep after rt_step_3,
// on the serial ray backend only (not the moment RT, not rt_backend = amrex).
// A snapshot is every live ray in the cell chains, cell order then chain
// order, with its provenance (RayProv, filled by ray_source, ray_cast_split
// and ray_cast_merge when tracking is on) and its cell's state, as a numpy
// structured array <stem>.npy plus a JSON sidecar <stem>.json with the run's
// parameters and the step's bookkeeping; tools/ray_compression reads them.
// Nothing here is read by the RT, so results are bit-identical with the
// dump on or off (the ray_dump_* CTests).
#include "radray/ray_dump.hpp"

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "radray/simulation.hpp"

namespace radray {

namespace {

// numpy dtype of RayDumpRecord (little endian, no padding): field by field
// the struct of ray_dump.hpp, whose static_assert pins the 180-byte size
const char* kDescr =
    "[('cell', '<i4'), ('chainpos', '<i4'), ('nside', '<i4'), ('ipix', '<i4'), "
    "('r', '<f4'), ('dr', '<f4'), ('x', '<f4', (3,)), ('dx', '<f4', (3,)), "
    "('nu', '<f4', (5,)), ('ngam', '<f4', (5,)), ('src', '<i4'), ('birth', '<i4'), "
    "('nmerge', '<u2'), ('nclamp', 'u1'), ('flags', 'u1'), "
    "('c_tau', '<f4', (5,)), ('c_fgam', '<f4', (5,)), ('c_ngam', '<f4', (5,)), "
    "('c_fHI', '<f4'), ('c_fHeI', '<f4'), ('c_fHeII', '<f4'), ('c_temp', '<f4'), "
    "('c_rhog', '<f4')]";

// .npy version 1.0 header for a 1-D structured array of n records: magic,
// version, little-endian header length, then the dict padded with spaces so
// that the data start on a 64-byte boundary (the format's alignment rule)
std::string npy_header(std::int64_t n) {
  std::ostringstream d;
  d << "{'descr': " << kDescr << ", 'fortran_order': False, 'shape': (" << n << ",), }";
  std::string h = d.str();
  const std::size_t pre = 10;  // magic(6) + version(2) + len(2)
  std::size_t total = pre + h.size() + 1;
  const std::size_t pad = (64 - total % 64) % 64;
  h.append(pad, ' ');
  h.push_back('\n');
  if (h.size() > 65535) throw std::runtime_error("npy header too long");
  std::string out("\x93NUMPY\x01\x00", 8);
  out.push_back(static_cast<char>(h.size() & 0xff));
  out.push_back(static_cast<char>((h.size() >> 8) & 0xff));
  return out + h;
}

// a JSON array of doubles; %.17g round-trips every value exactly
void json_doubles(std::ostream& o, const double* v, int n) {
  o << "[";
  char buf[40];
  for (int i = 0; i < n; ++i) {
    std::snprintf(buf, sizeof buf, "%.17g", v[i]);
    o << (i ? ", " : "") << buf;
  }
  o << "]";
}

}  // namespace

// ray_dump_at: substep numbers separated by blanks and/or commas
std::vector<std::int64_t> parse_substep_list(const std::string& s0) {
  std::string s = s0;
  for (auto& ch : s)
    if (ch == ',') ch = ' ';
  std::istringstream in(s);
  std::vector<std::int64_t> v;
  std::int64_t x;
  while (in >> x) v.push_back(x);
  return v;
}

// Turns the dump on: the provenance side array gets one slot per ray storage
// slot (1-based like the pool, hence Nray_max + 1) and is filled from now
// on, so rays created earlier would have no origin.
void ray_dump_enable(Simulation& s, const RayDumpConfig& cfg) {
  if (!cfg.enabled()) return;
  if (s.rays.Nray > 0)
    throw std::runtime_error("ray_dump_enable: must be called before rays exist");
  std::filesystem::create_directories(cfg.dir);
  s.raydump_state = std::make_shared<RayDump>();
  s.raydump_state->cfg = cfg;
  s.rays.track_prov = true;
  s.rays.prov.reset(s.rays.Nray_max + 1);
}

// After rt_step_3 of every substep: write a snapshot if this global substep
// is selected (every n-th, or listed in ray_dump_at), named by its number.
void ray_dump_after_substep(Simulation& s) {
  RayDump* d = s.raydump_state.get();
  if (!d) return;
  // The counter was advanced at the start of the substep (see
  // radiative_transfer), so prov_substep is the current substep.
  const std::int64_t k = d->substep;
  const bool hit = (d->cfg.every > 0 && k % d->cfg.every == 0) ||
                   std::find(d->cfg.at.begin(), d->cfg.at.end(), k) != d->cfg.at.end();
  if (!hit) return;
  char stem[64];
  std::snprintf(stem, sizeof stem, "rays_s%07lld", static_cast<long long>(k));
  ray_dump_write(s, d->cfg.dir + "/" + stem);
  ++d->ndumped;
}

// One snapshot: walks every cell's live chain (htoc_ray(3, c) the head,
// htoc_ray(4, c) the tail, grid.hpp) and copies each ray's fields, its
// provenance and the cell's state after the substep; on the way it collects
// the source cells (the RT's criterion sfr > sfr_min) and the volume
// fraction of HII for the sidecar.
void ray_dump_write(const Simulation& s, const std::string& stem) {
  const RTGrid& G = s.rt;
  const RayState& R = s.rays;
  const bool have_prov = R.track_prov && !R.prov.empty();

  std::vector<RayDumpRecord> rec;
  std::int64_t nsrc = 0;
  double xhii_vol = 0;
  std::vector<std::int64_t> src_cells;
  for (std::int64_t c = 0; c < G.ncell; ++c) {
    xhii_vol += 1.0 - static_cast<double>(G.fHI[c]);
    if (static_cast<double>(G.sfr[c]) > s.units.sfr_min) {
      ++nsrc;
      src_cells.push_back(c);
    }
    std::int64_t ir = G.htoc_ray(3, c);
    const std::int64_t tail = G.htoc_ray(4, c);
    std::int32_t pos = 0;
    while (ir > 0) {
      const Ray& ry = R.at(ir);
      RayDumpRecord q{};
      q.cell = static_cast<std::int32_t>(c);
      q.chainpos = pos++;
      q.nside = ry.nside;
      q.ipix = ry.ipix;
      q.r = ry.r;
      q.dr = ry.dr;
      std::memcpy(q.x, ry.x, sizeof q.x);
      std::memcpy(q.dx, ry.dx, sizeof q.dx);
      std::memcpy(q.nu, ry.nu, sizeof q.nu);
      std::memcpy(q.ngam, ry.ngam, sizeof q.ngam);
      if (have_prov && ir < R.prov.cap) {
        const RayProv& p = R.prov[ir];
        // birth >= 1 for every ray created while tracking; 0 means the ray
        // predates it (e.g. read from a restart checkpoint): unknown origin
        q.src = p.birth > 0 ? p.src : -2;
        q.birth = p.birth;
        q.nmerge = p.nmerge;
        q.nclamp = p.nclamp;
        q.flags = p.flags;
      } else {
        q.src = -2;  // unknown
      }
      std::memcpy(q.c_tau, &G.tau[c * kNfreq], sizeof q.c_tau);
      std::memcpy(q.c_fgam, &G.fgam[c * kNfreq], sizeof q.c_fgam);
      std::memcpy(q.c_ngam, &G.ngam[c * kNfreq], sizeof q.c_ngam);
      q.c_fHI = G.fHI[c];
      q.c_fHeI = G.fHeI[c];
      q.c_fHeII = G.fHeII[c];
      q.c_temp = G.temp[c];
      q.c_rhog = G.rhog[c];
      rec.push_back(q);
      if (ir == tail) break;
      ir = ry.next;
    }
  }
  xhii_vol /= static_cast<double>(G.ncell);

  // the records as one unpadded block after the numpy header
  {
    std::ofstream out(stem + ".npy", std::ios::binary);
    if (!out) throw std::runtime_error("cannot write " + stem + ".npy");
    const std::string h = npy_header(static_cast<std::int64_t>(rec.size()));
    out.write(h.data(), static_cast<std::streamsize>(h.size()));
    out.write(reinterpret_cast<const char*>(rec.data()),
              static_cast<std::streamsize>(rec.size() * sizeof(RayDumpRecord)));
  }
  // the sidecar: the grid and the flags that shape the rays, the clock and
  // the substep of the snapshot, the thresholds of the ray bookkeeping, the
  // sources and the ngam_ray threshold table of merging and splitting (64
  // entries, rays.hpp)
  {
    std::ofstream o(stem + ".json");
    char buf[64];
    auto num = [&](double v) {
      std::snprintf(buf, sizeof buf, "%.17g", v);
      return std::string(buf);
    };
    o << "{\n";
    o << "  \"N\": " << s.cfg.Nrt << ",\n";
    o << "  \"box\": " << num(s.cfg.box) << ",\n";
    o << "  \"nray_max\": " << s.cfg.Nray_max << ",\n";
    o << "  \"legacy_rt\": " << (s.cfg.legacy_rt ? 1 : 0) << ",\n";
    o << "  \"ray_rotate\": " << s.cfg.ray_rotate << ",\n";
    o << "  \"redshift_per_distance\": " << (s.cfg.redshift_per_distance ? 1 : 0) << ",\n";
    o << "  \"legacy_physics\": " << (s.cfg.legacy_physics ? 1 : 0) << ",\n";
    o << "  \"istep\": " << s.clk.istep << ",\n";
    o << "  \"istep_rt\": " << s.clk.istep_rt << ",\n";
    o << "  \"nstep_rt\": " << s.clk.nstep_rt << ",\n";
    o << "  \"substep\": " << (s.raydump_state ? s.raydump_state->substep : 0) << ",\n";
    o << "  \"redshift\": " << num(s.clk.redshift) << ",\n";
    o << "  \"dtsec_rt\": " << num(s.clk.dtsec_rt) << ",\n";
    o << "  \"ngam_min\": " << num(s.units.ngam_min) << ",\n";
    o << "  \"fgam_min\": " << num(s.units.fgam_min) << ",\n";
    o << "  \"nrays\": " << rec.size() << ",\n";
    o << "  \"xHII_vol\": " << num(xhii_vol) << ",\n";
    o << "  \"nsources\": " << nsrc << ",\n";
    o << "  \"source_cells\": [";
    for (std::size_t i = 0; i < src_cells.size(); ++i) o << (i ? ", " : "") << src_cells[i];
    o << "],\n";
    o << "  \"ngam_ray\": ";
    json_doubles(o, R.ngam_ray.data(), 64);
    o << "\n}\n";
  }
}

}  // namespace radray
