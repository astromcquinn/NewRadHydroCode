// Fortran-semantics helpers and folded constants vs a gfortran dump.
#include <cstdint>
#include <cstring>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "check.hpp"
#include "radray/constants.hpp"
#include "radray/fortran.hpp"
#include "radray/ray_tracing.hpp"

using namespace radray;

static std::uint64_t bits(double v) {
  std::uint64_t u;
  std::memcpy(&u, &v, 8);
  return u;
}

int main(int argc, char** argv) {
  // INT truncates toward zero
  CHECK(f::int4(2.9) == 2);
  CHECK(f::int4(-2.9) == -2);
  CHECK(f::int4(-0.5) == 0);
  CHECK(f::int4(std::nan("")) == INT32_MIN);
  CHECK(f::int4(1e300) == INT32_MIN);
  CHECK(f::int4(-std::numeric_limits<double>::infinity()) == INT32_MIN);
  // NINT rounds half away from zero
  CHECK(f::nint(2.5) == 3);
  CHECK(f::nint(-2.5) == -3);
  CHECK(f::nint(2.4999999) == 2);
  CHECK(f::nint(0.5) == 1);
  // FLOOR / CEILING
  CHECK(f::floor4(-0.5) == -1);
  CHECK(f::floor4(3.0) == 3);
  CHECK(f::ceiling4(3.0) == 3);
  CHECK(f::ceiling4(-0.5) == 0);
  CHECK(f::ceiling4(3.0000001) == 4);
  // MOD has the sign of the first argument
  CHECK(f::mod(7, 3) == 1);
  CHECK(f::mod(-7, 3) == -1);
  CHECK(f::mod(7, -3) == 1);
  CHECK(f::mod(-7.5, 2.0) == -1.5);
  CHECK(f::mod(33.0, 32.0) == 1.0);
  // Single-precision literals promoted
  CHECK(static_cast<double>(0.99f) != 0.99);
  CHECK(cp::omegam == static_cast<double>(0.27f));
  // real(4) sums accumulate in float
  {
    const float a[5] = {1.0f, 1e-8f, 1e-8f, 1e-8f, 1e-8f};
    CHECK(f::sum4<5>(a) == 1.0f);  // float accumulation loses the small terms
  }
  // integer powers
  CHECK(f::pow4(3.0) == 81.0);

  // Compare with the gfortran dump (tests/fortran_consts/expected_consts.txt)
  if (argc > 1) {
    std::ifstream in(argv[1]);
    CHECK(static_cast<bool>(in));
    std::map<std::string, std::vector<std::string>> m;
    std::string line;
    while (std::getline(in, line)) {
      std::istringstream ss(line);
      std::string k, v;
      ss >> k >> v;
      if (!k.empty()) m[k].push_back(v);
    }
    auto hex = [](double v) {
      char b[32];
      std::snprintf(b, sizeof b, "%016llX", static_cast<unsigned long long>(bits(v)));
      return std::string(b);
    };
    auto expect = [&](const char* key, int idx, double v) {
      auto it = m.find(key);
      if (it == m.end() || static_cast<int>(it->second.size()) <= idx) {
        std::fprintf(stderr, "missing key %s\n", key);
        ++radray_test::failures();
        return;
      }
      if (it->second[idx] != hex(v)) {
        std::fprintf(stderr, "constant %s[%d]: gfortran %s, C++ %s (%.17g)\n", key, idx,
                     it->second[idx].c_str(), hex(v).c_str(), v);
        ++radray_test::failures();
      }
    };
    expect("pi", 0, pc::pi);
    expect("omegam", 0, cp::omegam);
    expect("omegal", 0, cp::omegal);
    expect("omegab", 0, cp::omegab);
    expect("omegar", 0, cp::omegar);
    expect("hubble0", 0, cp::hubble0);
    expect("XHy", 0, cp::XHy);
    expect("YHe", 0, cp::YHe);
    expect("gamma_ideal", 0, cp::gamma_ideal);
    expect("Tcmb0", 0, cp::Tcmb0);
    for (int i = 0; i < 5; ++i) {
      expect("freq_init", i, cp::freq_init[i]);
      expect("pop2_init", i, cp::pop2_init[i]);
    }
    expect("fResc", 0, cp::fResc);
    expect("eSF", 0, cp::eSF);
    expect("eSN", 0, cp::eSN);
    expect("lgT_del", 0, cp::lgT_del);
    expect("lgn_del", 0, folded::lgn_del);
    expect("xT", 0, folded::temp_x);
    expect("log4", 0, folded::log4_r4);
    expect("Dmin", 0, folded::ten_lgD_min);
    expect("Dmax", 0, folded::ten_lgD_max);
    expect("Tmax", 0, folded::ten_lgT_max);
    // Expressions the C++ evaluates at run/compile time (no hardcoding):
    const double om = cp::omegam, ol = cp::omegal, h0 = cp::hubble0;
    expect("yrpre", 0, 2 / (3 * pc::H0_cgs * h0 * std::sqrt(ol)));
    expect("dlpre", 0, 2 * pc::c_cgs / (pc::H0_cgs * h0 * std::sqrt(om)));
    expect("s3pi", 0, std::sqrt(3 / (pc::pi * rayp::ray_per_cell)));
    expect("cc14", 0,
           4 * pc::sS_cgs / pc::c_cgs * static_cast<double>(8.0f / 3.0f) *
               static_cast<double>(3.36e-24f) / (cp::gamma_ideal - 1));
    expect("freq_edge2", 0, std::sqrt(pc::HI_eV * pc::HeI_eV));
    expect("freq_edge4", 0, std::sqrt(pc::HeI_eV * pc::HeII_eV));
    expect("HI_eV", 0, pc::HI_eV);
    expect("HI_K", 0, pc::HI_K);
    expect("HeI_K", 0, pc::HeI_K);
    expect("HeII_K", 0, pc::HeII_K);
    expect("sec2yr", 0, pc::sec2yr);
    expect("tdpre", 0, 3 * pc::pi);
    expect("gpre", 0, 32 * pc::G_cgs);
    expect("c2", 0, pc::c_cgs * pc::c_cgs);
    expect("sqom3", 0, 3 * std::sqrt(om));
    expect("bdel", 0, static_cast<double>(0.1f));
    // integer parameters
    const int N_pdf = f::nint((3.0 - (-6.0)) / folded::lgn_del);
    CHECK(m["N_pdf"].size() == 1 && std::stoi(m["N_pdf"][0]) == N_pdf);
    const int Nbin = 1 + f::nint((6.0 - (-6.0)) / static_cast<double>(0.1f));
    CHECK(m["Nbin"].size() == 1 && std::stoi(m["Nbin"][0]) == Nbin);
  }
  TEST_MAIN_END();
}
