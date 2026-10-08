// Compatibility with the original external halo catalogs:
//  * read_halo port (sfr.txt scan, redshift tag, rescaling) against the
//    verbatim Fortran read_halo in tests/halo_catalog_ref.f90 (gfortran),
//    bit for bit;
//  * reading catalogs written by gfortran (stream = Intel form='binary',
//    and sequential with record markers);
//  * writing catalogs read back by gfortran;
//  * the redshift tag (Fortran F editing) against gfortran;
//  * link_halo + calc_halo + add_halo_to_rt metal injection on a grid.
//
// usage: test_halo_catalog <halo_catalog_ref exe or ""> <work dir>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "check.hpp"
#include "radray/constants.hpp"
#include "radray/halo_sources.hpp"
#include "radray/rt.hpp"

using namespace radray;
namespace fs = std::filesystem;

namespace {

std::string g_ref;
std::string g_work;

std::vector<std::string> run_ref(const std::string& args) {
  const std::string out = g_work + "/ref_out.txt";
  const std::string cmd = "\"" + g_ref + "\" " + args + " > \"" + out + "\"";
  const int rc = std::system(cmd.c_str());
  CHECK(rc == 0);
  std::vector<std::string> lines;
  std::ifstream in(out);
  std::string l;
  while (std::getline(in, l)) lines.push_back(l);
  return lines;
}

std::uint32_t bits(float f) {
  std::uint32_t u;
  std::memcpy(&u, &f, 4);
  return u;
}
std::int64_t bits(double d) {
  std::int64_t u;
  std::memcpy(&u, &d, 8);
  return u;
}

void write_sfr_txt(const std::string& path) {
  std::ofstream o(path);
  // z, 8 dummies, sfr. Exercises list-directed features: a row continued on
  // the next line, commas, a 'd' exponent, extra values (discarded), and an
  // exact decimal tie in the tag (8.03125 -> f6.4).
  o << "14.5  1 2 3 4 5 6 7 8  1.0e-4\n";
  o << "12.25 1 2 3 4\n 5 6 7 8  3.5e-4\n";
  o << "10.0,1,2,3,4,5,6,7,8,1.2d-3\n";
  o << "9.0   1 2 3 4 5 6 7 8  2.0e-3   99 99\n";
  o << "8.03125 0 0 0 0 0 0 0 0 3.1e-3\n";
  o << "7.5   0 0 0 0 0 0 0 0 4.4e-3\n";
  o << "6.0   0 0 0 0 0 0 0 0 9.0e-3\n";
  o << "5.2   0 0 0 0 0 0 0 0 1.3e-2\n";
}

void test_read_halo_vs_fortran() {
  const std::string dir = g_work + "/Halo/";
  fs::create_directories(dir);
  write_sfr_txt(dir + "sfr.txt");
  const auto hist = read_sfr_history(dir + "sfr.txt");
  CHECK(hist.size() == 8);
  CHECK(hist[1].z == 12.25 && hist[1].sfr == 3.5e-4 && hist[1].dummy[7] == 8);
  CHECK(hist[2].sfr == 1.2e-3);
  CHECK(hist[3].sfr == 2.0e-3);
  if (g_ref.empty()) {
    std::printf("  (no gfortran: Fortran comparisons skipped)\n");
    return;
  }
  run_ref("write " + dir);
  int ncompared = 0;
  for (const char* zs : {"15.0", "14.0", "12.25", "11.0", "9.5", "8.6", "8.0", "7.0", "5.5",
                         "4.0", "13.0001"}) {
    const double z = std::strtod(zs, nullptr);
    const auto ref = run_ref(std::string("read ") + dir + " " + zs);
    const SourceCatalog c = read_halo(dir, z);
    CHECK(ref.size() >= 2);
    if (ref.size() < 2) continue;
    std::istringstream h(ref[1]);
    long long nh = -1, ns = -1, sc = 0;
    h >> nh >> ns >> sc;
    if (ref[0] == "none") {
      CHECK(c.halo.empty() && c.Nsource == 0);
      std::printf("  z=%-8s none\n", zs);
      continue;
    }
    CHECK(c.tag == ref[0]);
    CHECK(static_cast<long long>(c.halo.size()) == nh);
    CHECK(c.Nsource == ns);
    CHECK(bits(c.sfr_scale) == sc);
    CHECK(static_cast<long long>(ref.size()) == 2 + ns);
    for (long long i = 0; i < ns && 2 + i < static_cast<long long>(ref.size()); ++i) {
      std::istringstream r(ref[2 + i]);
      long long L = 0, s = 0;
      r >> L >> s;
      CHECK(static_cast<std::int32_t>(bits(c.halo[i].L)) == L);
      CHECK(static_cast<std::int32_t>(bits(c.halo[i].sfr)) == s);
    }
    std::printf("  z=%-8s %s Nhalo=%lld Nsource=%lld scale=%.6f\n", zs, c.tag.c_str(), nh, ns,
                c.sfr_scale);
    ++ncompared;
    // The sequential copy written by gfortran reads back identically
    const auto seq = read_halo_catalog(dir + "seq.halo." + c.tag);
    const auto str = read_halo_catalog(dir + "halo." + c.tag, CatalogFormat::Stream);
    CHECK(seq.size() == str.size());
    CHECK(!seq.empty() && std::memcmp(seq.data(), str.data(), seq.size() * 56) == 0);
  }
  CHECK(ncompared >= 8);
}

void test_write_read_back() {
  std::mt19937 rng(5);
  std::uniform_real_distribution<float> U(-1e3f, 1e3f);
  std::vector<SourceRecord> h(37);
  for (std::size_t i = 0; i < h.size(); ++i) {
    auto& r = h[i];
    for (int d = 0; d < 3; ++d) {
      r.x[d] = std::fabs(U(rng)) / 200;
      r.v[d] = U(rng);
    }
    r.N = static_cast<std::int32_t>(1000 - i);
    r.M = 1e9f * (1 + std::fabs(U(rng)));
    r.R = 0.01f;
    r.T = 2e4f;
    r.L = i < 20 ? 1e27f * (i + 1) : 0.0f;
    r.sfr = i < 20 ? 0.1f : 0.0f;
    r.mstar = 1e6f;
    r.fgas = 0.15f;
  }
  const std::string fs_ = g_work + "/cpp.stream", fq = g_work + "/cpp.seq";
  write_halo_catalog(fs_, h, CatalogFormat::Stream);
  write_halo_catalog(fq, h, CatalogFormat::Sequential);
  CHECK(fs::file_size(fs_) == 8 + 56 * h.size());
  CHECK(fs::file_size(fq) == 16 + 64 * h.size());
  for (const auto& f : {fs_, fq}) {
    const auto b = read_halo_catalog(f);
    CHECK(b.size() == h.size() && std::memcmp(b.data(), h.data(), 56 * h.size()) == 0);
  }
  if (g_ref.empty()) return;
  const auto ref = run_ref("dump " + fs_);
  CHECK(ref.size() == 1 + h.size());
  if (ref.size() != 1 + h.size()) return;
  CHECK(std::stoll(ref[0]) == static_cast<long long>(h.size()));
  for (std::size_t i = 0; i < h.size(); ++i) {
    std::istringstream r(ref[1 + i]);
    long long v[14];
    for (auto& x : v) r >> x;
    const auto& s = h[i];
    const float fl[13] = {s.x[0], s.x[1], s.x[2], s.v[0], s.v[1], s.v[2], 0,
                          s.M,    s.R,    s.T,    s.L,    s.sfr,  s.mstar};
    for (int k = 0; k < 13; ++k) {
      if (k == 6) {
        CHECK(v[6] == s.N);
      } else {
        CHECK(v[k] == static_cast<std::int32_t>(bits(fl[k])));
      }
    }
    CHECK(v[13] == static_cast<std::int32_t>(bits(s.fgas)));
  }
}

void test_tags() {
  CHECK(redshift_tag(8) == "z=08.0000");
  CHECK(redshift_tag(12.5) == "z=12.5000");
  CHECK(redshift_tag(0.5) == "z=00.5000");
  bool threw = false;
  try {
    redshift_tag(150);
  } catch (const std::exception&) {
    threw = true;
  }
  CHECK(threw);
  if (g_ref.empty()) return;
  int n = 0;
  for (const char* zs : {"8", "6.03125", "6.03135", "0.00005", "9.99996", "9.99994", "0.5",
                         "12.000049", "99.99996", "0.15625", "0.28125", "0.46875", "1.00005",
                         "10.00005", "10.00015", "7.123456789", "19.99995", "5.999949999"}) {
    const auto ref = run_ref(std::string("tag ") + zs);
    const std::string t = redshift_tag(std::strtod(zs, nullptr));
    if (ref.empty() || t != ref[0])
      std::fprintf(stderr, "tag mismatch for z=%s: C++ '%s' gfortran '%s'\n", zs, t.c_str(),
                   ref.empty() ? "" : ref[0].c_str());
    CHECK(!ref.empty() && t == ref[0]);
    ++n;
  }
  std::printf("  %d redshift tags identical to gfortran\n", n);
}

void test_link_and_metals() {
  Simulation s;
  init_simulation(s, 8, 4, 9.0);
  s.gas.allocate(8);
  s.clk.dt_yr = 2.5e6;
  calc_sim_var(s.cfg, s.clk, s.units);
  const double dx = s.cfg.x_rt_box;  // box/Nrt Mpc/h
  std::vector<SourceRecord> src(5);
  // positions (Mpc/h) -> expected 1-based cells
  const float pos[5][3] = {{0.1f * static_cast<float>(dx), 0.0f, 0.0f},
                           {static_cast<float>(3.5 * dx), static_cast<float>(7.99 * dx), 0.0f},
                           {static_cast<float>(3.2 * dx), static_cast<float>(7.5 * dx), 0.0f},
                           {static_cast<float>(5 * dx), static_cast<float>(5 * dx), static_cast<float>(5 * dx)},
                           {1, 1, 1}};
  for (int n = 0; n < 5; ++n) {
    for (int d = 0; d < 3; ++d) src[n].x[d] = pos[n][d];
    src[n].L = 1e27f * (n + 1);
    src[n].sfr = 0.1f * (n + 1);
  }
  src[4].L = 0;  // not a source
  const auto cells = link_halo(src, 4, s.cfg);
  CHECK(cells.size() == 4);
  CHECK(cells[0].i == 1 && cells[0].j == 1 && cells[0].k == 1);
  CHECK(cells[1].i == 4 && cells[1].j == 8 && cells[1].k == 1);
  CHECK(cells[2].i == 4 && cells[2].j == 8 && cells[2].k == 1);
  CHECK(cells[3].k >= 5 && cells[3].k <= 6);
  calc_halo_records(s, src, 4, true);
  const auto c48 = s.rt.cell(4, 8, 1);
  CHECK(s.rt.Luv[c48] == static_cast<float>(2e27 + 3e27));
  CHECK_CLOSE(s.rt_stat.sfr, 0.1 + 0.2 + 0.3 + 0.4, 1e-6);
  CHECK(s.Nsource == 4);
  // Metals: rhoZ/rho_unit per source (single cell, single source)
  const double fZ = static_cast<double>(0.10f) * static_cast<double>(0.02f);
  const double expect = fZ * (static_cast<double>(0.1f) * 2.5e6 * pc::Msun_cgs) /
                        s.units.vol_gas / s.units.rho_unit;
  CHECK(s.gas.w[s.gas.cell(1, 1, 1) * 4 + 3] == static_cast<float>(expect));
  double tot = 0;
  for (std::int64_t c = 0; c < s.gas.ncell; ++c) tot += s.gas.w[c * 4 + 3];
  CHECK_CLOSE(tot, expect * (1 + 2 + 3 + 4), 1e-5);
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 3) {
    std::fprintf(stderr, "usage: %s <halo_catalog_ref|\"\"> <workdir>\n", argv[0]);
    return 2;
  }
  g_ref = argv[1];
  g_work = argv[2];
  fs::remove_all(g_work);
  fs::create_directories(g_work);
  std::printf("read_halo vs Fortran\n");
  test_read_halo_vs_fortran();
  std::printf("write/read back\n");
  test_write_read_back();
  std::printf("tags\n");
  test_tags();
  std::printf("link_halo / calc_halo / metals\n");
  test_link_and_metals();
  TEST_MAIN_END();
}
