#include "radray/io.hpp"

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
  for (auto& ch : val)
    if (ch == 'd' || ch == 'D') ch = 'e';  // accept Fortran exponents
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
  else if (key == "legacy_rt") {
    p.legacy_rt = std::stoi(val);
    if (p.legacy_rt != 0 && p.legacy_rt != 1)
      throw std::runtime_error("legacy_rt must be 0 or 1");
  } else if (key == "c_reduced") {
    p.c_reduced = std::stod(val);
    if (!(p.c_reduced > 0)) throw std::runtime_error("c_reduced must be > 0");
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
