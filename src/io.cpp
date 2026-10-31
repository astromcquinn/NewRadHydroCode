#include "radray/io.hpp"

#include "radray/cosmology.hpp"
#include "radray/plane_parallel.hpp"

#include "radray/rt_backend.hpp"
#include "radray/rt_moments.hpp"

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace radray {

namespace {
std::string trim(const std::string& s) {
  const auto a = s.find_first_not_of(" \t\r\n");
  if (a == std::string::npos) return "";
  const auto b = s.find_last_not_of(" \t\r\n");
  return s.substr(a, b - a + 1);
}

bool little_endian() {
  const std::uint32_t one = 1;
  unsigned char b;
  std::memcpy(&b, &one, 1);
  return b == 1;
}
}  // namespace

void set_param(ProblemParams& p, const std::string& key_in, std::string val) {
  const std::string key = trim(key_in);
  val = trim(val);
  if (key == "ray_dump_dir") {
    p.ray_dump_dir = val;
    return;
  }
  if (set_rt_method_param(p.rtm, key, val)) return;  // rt_moments.hpp
  if (set_rt_backend_param(p.rtm, key, val)) return;  // rt_backend.hpp
  if (set_plane_param(p.rtm, key, val)) return;       // plane_parallel.hpp
  for (auto& ch : val)
    if (ch == 'd' || ch == 'D') ch = 'e';  // accept Fortran exponents
  // omegam, omegab, hubble0, sigma8, ns, T_min, z_floor (cosmology.hpp)
  if (set_cosmology_param(p.rtm, key, val)) return;
  if (key == "N")
    p.N = std::stoi(val);
  else if (key == "z_init")
    p.z_init = std::stod(val);
  else if (key == "da_sim")
    p.da_sim = std::stod(val);
  else if (key == "nsteps")
    p.nsteps = std::stoi(val);
  else if (key == "nray_max_factor")
    p.nray_max_factor = std::stoi(val);
  else if (key == "box") {
    p.box = std::stod(val);
    if (!(p.box > 0)) throw std::runtime_error("box must be > 0");
  }
  else if (key == "ray_dump_every")
    p.ray_dump_every = std::stoll(val);
  else if (key == "ray_dump_at")
    p.ray_dump_at = val;
  else if (key == "ray_quant_bits")
    p.ray_quant_bits = std::stoi(val);
  else if (key == "ray_quant_dither")
    p.ray_quant_dither = std::stoi(val);
  else if (key == "ray_quant_struct")
    p.ray_quant_struct = std::stoi(val);
  else if (key == "ray_quant_offbits")
    p.ray_quant_offbits = std::stoi(val);
  else if (key == "legacy_rt") {
    p.legacy_rt = std::stoi(val);
    if (p.legacy_rt != 0 && p.legacy_rt != 1)
      throw std::runtime_error("legacy_rt must be 0 or 1");
  } else if (key == "deterministic_parallel") {
    p.deterministic_parallel = std::stoi(val);
    if (p.deterministic_parallel != 0 && p.deterministic_parallel != 1)
      throw std::runtime_error("deterministic_parallel must be 0 or 1");
  } else if (key == "threads") {
    p.threads = std::stoi(val);
  } else if (key == "c_reduced") {
    p.c_reduced = std::stod(val);
    if (!(p.c_reduced > 0)) throw std::runtime_error("c_reduced must be > 0");
  } else if (key == "ray_rotate") {
    p.ray_rotate = std::stoi(val);
    if (p.ray_rotate < 0 || p.ray_rotate > 3) throw std::runtime_error("ray_rotate must be 0..3");
  } else if (key == "ray_rotate_seed") {
    p.ray_rotate_seed = std::stoll(val);
  } else if (key == "redshift_per_distance") {
    p.redshift_per_distance = std::stoi(val);
    if (p.redshift_per_distance != 0 && p.redshift_per_distance != 1)
      throw std::runtime_error("redshift_per_distance must be 0 or 1");
  } else if (key == "rt_cells_per_substep") {
    p.rt_cells_per_substep = std::stoi(val);
    if (p.rt_cells_per_substep < 1) throw std::runtime_error("rt_cells_per_substep must be >= 1");
  } else if (key == "ray_deposit") {
    if (val == "share" || val == "0")
      p.ray_deposit = 0;
    else if (val == "exposure" || val == "1")
      p.ray_deposit = 1;
    else
      throw std::runtime_error("ray_deposit must be share or exposure");
  } else if (key == "cool_suppress_overdensity") {
    p.cool_suppress_overdensity = std::stod(val);
    if (!(p.cool_suppress_overdensity >= 0)) throw std::runtime_error("cool_suppress_overdensity must be >= 0 (0 = off)");
  } else if (key == "sn_heating") {
    p.sn_heating = std::stoi(val);
    if (p.sn_heating != 0 && p.sn_heating != 1) throw std::runtime_error("sn_heating must be 0 or 1");
  } else if (key == "photon_conserve") {
    p.photon_conserve = std::stoi(val);  // -1: follow legacy_physics (the default)
    if (p.photon_conserve < -1 || p.photon_conserve > 1)
      throw std::runtime_error("photon_conserve must be 0 or 1");
  } else if (key == "photon_consume_dense") {
    p.photon_consume_dense = std::stoi(val);  // -1: follow legacy_physics (the default)
    if (p.photon_consume_dense < -1 || p.photon_consume_dense > 1)
      throw std::runtime_error("photon_consume_dense must be 0 or 1");
  } else if (key == "rt_substep_carry") {
    p.rt_substep_carry = std::stoi(val);
    if (p.rt_substep_carry != 0 && p.rt_substep_carry != 1)
      throw std::runtime_error("rt_substep_carry must be 0 or 1");
  } else if (key == "legacy_physics") {
    p.legacy_physics = std::stoi(val);
    if (p.legacy_physics != 0 && p.legacy_physics != 1)
      throw std::runtime_error("legacy_physics must be 0 or 1");
  }
}

ProblemParams read_params(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open " + path);
  ProblemParams p;
  std::string line;
  while (std::getline(in, line)) {
    const auto h = line.find('#');
    if (h != std::string::npos) line = line.substr(0, h);
    const auto eq = line.find('=');
    if (eq == std::string::npos) continue;
    set_param(p, line.substr(0, eq), line.substr(eq + 1));
  }
  return p;
}

void read_gas(const std::string& path, GasGrid& gas) {
  if (!little_endian()) throw std::runtime_error("big-endian hosts not supported");
  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::runtime_error("cannot open " + path);
  std::vector<float> buf(static_cast<std::size_t>(9 * gas.ncell));
  in.read(reinterpret_cast<char*>(buf.data()), buf.size() * sizeof(float));
  if (!in || in.gcount() != static_cast<std::streamsize>(buf.size() * sizeof(float)))
    throw std::runtime_error("gas.bin has wrong size: " + path);
  in.peek();
  if (!in.eof()) throw std::runtime_error("gas.bin larger than (9,N,N,N): " + path);
  for (std::int64_t c = 0; c < gas.ncell; ++c) {
    for (int m = 0; m < 5; ++m) gas.u[c * 5 + m] = buf[c * 9 + m];
    for (int m = 0; m < 4; ++m) gas.w[c * 4 + m] = buf[c * 9 + 5 + m];
  }
}

std::vector<Source> read_sources(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open " + path);
  std::vector<Source> v;
  std::string line;
  while (std::getline(in, line)) {
    const auto h = line.find('#');
    if (h != std::string::npos) line = line.substr(0, h);
    if (trim(line).empty()) continue;
    for (auto& ch : line)
      if (ch == 'd' || ch == 'D') ch = 'e';
    std::istringstream ss(line);
    Source s;
    double L, sfr;
    if (!(ss >> s.i >> s.j >> s.k >> L >> sfr))
      throw std::runtime_error("bad line in " + path + ": " + line);
    // source_type%L and %sfr are real(4)
    s.Luv = static_cast<float>(L);
    s.sfr = static_cast<float>(sfr);
    v.push_back(s);
  }
  return v;
}

void write_raw_f32(const std::string& path, const std::vector<float>& v) {
  if (!little_endian()) throw std::runtime_error("big-endian hosts not supported");
  std::ofstream out(path, std::ios::binary);
  if (!out) throw std::runtime_error("cannot write " + path);
  out.write(reinterpret_cast<const char*>(v.data()), v.size() * sizeof(float));
}

void write_step_output(const std::string& dir, const Simulation& s) {
  const RTGrid& G = s.rt;
  const auto P = [&](const char* name) { return dir + "/" + name + ".bin"; };
  write_raw_f32(P("fHI"), G.fHI);
  write_raw_f32(P("fHeI"), G.fHeI);
  write_raw_f32(P("fHeII"), G.fHeII);
  write_raw_f32(P("temp"), G.temp);
  write_raw_f32(P("rhog"), G.rhog);
  write_raw_f32(P("ngam"), G.ngam);
  write_raw_f32(P("nu"), G.nu);
  write_raw_f32(P("fgam"), G.fgam);
  write_raw_f32(P("tau"), G.tau);
  write_raw_f32(P("gamma"), G.gamma);
  write_raw_f32(P("lambda"), G.lambda);
  write_raw_f32(P("dx"), G.dx);
  write_raw_f32(P("u_gas"), s.gas.u);
  write_raw_f32(P("w_gas"), s.gas.w);

  std::FILE* fp = std::fopen((dir + "/stats.txt").c_str(), "w");
  if (!fp) throw std::runtime_error("cannot write stats.txt in " + dir);
  const RayStat& rv = s.rays.save;
  std::fprintf(fp, "istep = %d\n", s.clk.istep);
  std::fprintf(fp, "redshift = %.17g\n", s.clk.redshift);
  std::fprintf(fp, "scalefactor = %.17g\n", s.clk.scalefactor);
  std::fprintf(fp, "nstep_rt = %d\n", s.clk.nstep_rt);
  std::fprintf(fp, "nray_active = %lld\n",
               static_cast<long long>(s.rays.Nray - s.rays.save.nfree));
  std::fprintf(fp, "ncst1 = %lld\n", static_cast<long long>(rv.ncst1));
  std::fprintf(fp, "nmrg1 = %lld\n", static_cast<long long>(rv.nmrg1));
  std::fprintf(fp, "nspl1 = %lld\n", static_cast<long long>(rv.nspl1));
  std::fprintf(fp, "ndel1 = %lld\n", static_cast<long long>(rv.ndel1));
  std::fprintf(fp, "ngam_rays_total = %.17g\n", ngam_rays_total(s));
  if (s.plane.on) {  // plane_parallel = 1: the step's photon bookkeeping (plane_parallel.hpp)
    std::fprintf(fp, "plane_emitted = %.17g\n", s.plane.emitted);
    std::fprintf(fp, "plane_absorbed = %.17g\n", s.plane.absorbed());
    std::fprintf(fp, "plane_escaped = %.17g\n", s.plane.escaped());
  }
  std::fclose(fp);

  // rt_stats.txt: the extra diagnostics written by the Fortran harness
  // (reference/harness/src/rt_harness.f90: write_step_output), same keys.
  fp = std::fopen((dir + "/rt_stats.txt").c_str(), "w");
  if (!fp) throw std::runtime_error("cannot write rt_stats.txt in " + dir);
  const RTStat& v = s.rt_save;
  const Units& u = s.units;
  auto d1 = [&](const char* k, double x) { std::fprintf(fp, "%s = %.17g\n", k, x); };
  auto d3 = [&](const char* k, const std::array<double, 3>& a) {
    std::fprintf(fp, "%s = %.17g %.17g %.17g\n", k, a[0], a[1], a[2]);
  };
  d1("year", s.clk.year);
  d1("dt_sim", s.clk.dt_sim);
  d1("dt_yr", s.clk.dt_yr);
  d1("dtsec_rt", s.clk.dtsec_rt);
  d1("temp_min", u.temp_min);
  d1("nH_avg", u.nH_avg);
  d1("nHe_avg", u.nHe_avg);
  d1("vol_rt", u.vol_rt);
  std::fprintf(fp, "Nray = %lld\n", static_cast<long long>(s.rays.Nray));
  std::fprintf(fp, "Nsfill = %lld\n", static_cast<long long>(v.nsfill));
  std::fprintf(fp, "Ngfill = %lld\n", static_cast<long long>(v.ngfill));
  d1("Luv", v.Luv);
  d1("sfr", v.sfr);
  d3("ngam", v.ngam);
  d3("nion", v.nion);
  d3("ncol", v.ncol);
  d3("nrec", v.nrec);
  d3("nrab", v.nrab);
  d3("gamma", v.gamma);
  d3("lambda", v.lambda);
  d1("tau", v.tau);
  d1("dtau", v.dtau);
  std::fprintf(fp, "ngam_ray =");
  for (double x : s.rays.ngam_ray) std::fprintf(fp, " %.17g", x);
  std::fprintf(fp, "\n");
  std::fclose(fp);
}

}  // namespace radray
