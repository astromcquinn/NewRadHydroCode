// Stand-alone halo finder + source model on a DM snapshot.
//
// usage: find_halos <xv_dm file | --synthetic Ndm> <redshift> <out catalog>
//                   [--b 0.2] [--nmin 20] [--so] [--potmin]
//                   [--model global|dynamical] [--sfr Halo/sfr.txt]
//                   [--tvir 1e4] [--seq]
//
// xv_dm: raw float32 xv_dm(6,Npart) as written by io_tools.f90
// (form='binary', no header), positions in DM-grid units [0,Ndm). Ndm is
// taken from the file size, Ngas = 2*Ndm (radhydro_global.f90). The output
// is a legacy halo catalog (stream format unless --seq), readable by the
// original read_halo. --synthetic generates a clustered test distribution
// (for timing).
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <random>
#include <string>
#include <vector>

#include "radray/constants.hpp"
#include "radray/cosmology.hpp"
#include "radray/halo_sources.hpp"
#include "radray/halos.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace radray;

static std::vector<float> synthetic(int Ndm, unsigned seed) {
  // Uniform background plus Gaussian blobs with a power-law mass spectrum
  const double L = Ndm;
  const std::int64_t N = static_cast<std::int64_t>(Ndm) * Ndm * Ndm;
  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<double> U(0, 1);
  std::normal_distribution<double> G(0, 1);
  std::vector<float> xv;
  xv.reserve(6 * N);
  auto add = [&](double x, double y, double z) {
    for (double a : {x, y, z}) {
      a = std::fmod(a, L);
      if (a < 0) a += L;
      float f = static_cast<float>(a);
      xv.push_back(f >= static_cast<float>(L) ? 0.0f : f);
    }
    for (int d = 0; d < 3; ++d) xv.push_back(static_cast<float>(G(rng)));
  };
  std::int64_t used = 0;
  while (used < N / 3) {
    const double n = 20 * std::pow(1 - U(rng), -1.0);  // dN/dn ~ n^-2
    const std::int64_t np = static_cast<std::int64_t>(std::min(n, 2e5));
    const double c[3] = {U(rng) * L, U(rng) * L, U(rng) * L};
    const double s = 0.08 * std::cbrt(static_cast<double>(np));
    for (std::int64_t q = 0; q < np; ++q) add(c[0] + s * G(rng), c[1] + s * G(rng), c[2] + s * G(rng));
    used += np;
  }
  while (static_cast<std::int64_t>(xv.size() / 6) < N) add(U(rng) * L, U(rng) * L, U(rng) * L);
  return xv;
}

int main(int argc, char** argv) {
  if (argc < 4) {
    std::fprintf(stderr,
                 "usage: %s <xv_dm|--synthetic Ndm> <redshift> <out> [--b x] [--nmin n] [--so]\n"
                 "       [--potmin] [--model global|dynamical] [--sfr file] [--tvir T] [--seq]\n"
                 "       [--box Mpc/h (default 5)] [--ngas Ngas (default 2*Ndm)]\n"
                 "       [--cosmo omegam=..|omegab=..|hubble0=..|sigma8=..|ns=.. (repeatable)]\n",
                 argv[0]);
    return 2;
  }
  try {
    int a = 1;
    std::vector<float> xv;
    int Ndm = 0;
    if (std::string(argv[a]) == "--synthetic") {
      Ndm = std::atoi(argv[a + 1]);
      xv = synthetic(Ndm, 1);
      a += 2;
    } else {
      const std::string fn = argv[a++];
      const auto size = std::filesystem::file_size(fn);
      const std::int64_t np = static_cast<std::int64_t>(size / 24);
      Ndm = static_cast<int>(std::llround(std::cbrt(static_cast<double>(np))));
      if (static_cast<std::int64_t>(Ndm) * Ndm * Ndm != np || size % 24)
        throw std::runtime_error(fn + ": size is not 24*Ndm^3 bytes");
      xv.resize(6 * np);
      std::ifstream in(fn, std::ios::binary);
      in.read(reinterpret_cast<char*>(xv.data()), static_cast<std::streamsize>(size));
      if (!in) throw std::runtime_error("cannot read " + fn);
    }
    const double z = std::atof(argv[a++]);
    const std::string out = argv[a++];
    HaloFinderParams fp;
    SourceModelParams mp;
    CatalogFormat fmt = CatalogFormat::Stream;
    std::vector<SfrHistoryRow> table;
    double box = cp::box_default;
    int ngas = 0;
    Cosmology cosmo;  // --cosmo key=value (the original's by default)
    for (; a < argc; ++a) {
      const std::string s = argv[a];
      if (s == "--b") fp.b = std::atof(argv[++a]);
      else if (s == "--nmin") fp.min_npart = std::atoll(argv[++a]);
      else if (s == "--so") fp.so_refine = true;
      else if (s == "--potmin") fp.centre = HaloCentre::PotentialMinimum;
      else if (s == "--model") mp.model = std::string(argv[++a]) == "dynamical" ? SfrModel::Dynamical : SfrModel::GlobalMassProportional;
      else if (s == "--sfr") table = read_sfr_history(argv[++a]);
      else if (s == "--tvir") mp.Tvir_min = std::atof(argv[++a]);
      else if (s == "--seq") fmt = CatalogFormat::Sequential;
      else if (s == "--box") box = std::atof(argv[++a]);
      else if (s == "--ngas") ngas = std::atoi(argv[++a]);
      else if (s == "--cosmo") {
        const std::string kv = argv[++a];
        const auto eq = kv.find('=');
        if (eq == std::string::npos || !set_cosmology_key(cosmo, kv.substr(0, eq), kv.substr(eq + 1)))
          throw std::runtime_error("--cosmo omegam|omegab|hubble0|sigma8|ns=value");
      } else throw std::runtime_error("unknown option " + s);
    }
    Config cfg = Config::make(ngas > 0 ? ngas : 2 * Ndm, 64, Ndm, box);
    cosmo.derive();
    cfg.cosmo = cosmo;
    SimClock clk;
    init_variables(cfg.cosmo, z, clk);
    Units u;
    calc_sim_var(cfg, clk, u);
    const HaloUnits hu = make_halo_units(cfg, clk, u);
    const std::int64_t np = static_cast<std::int64_t>(xv.size() / 6);
    const auto t0 = std::chrono::steady_clock::now();
    const auto halos = find_halos(xv.data(), np, hu, fp);
    const auto t1 = std::chrono::steady_clock::now();
    const SourceCatalog cat = assign_sources(halos, hu, mp, table.empty() ? nullptr : &table);
    write_halo_catalog(out, cat.halo, fmt);
    int nt = 1;
#ifdef _OPENMP
    nt = omp_get_max_threads();
#endif
    double sfr = 0;
    for (std::int64_t i = 0; i < cat.Nsource; ++i) sfr += cat.halo[i].sfr;
    std::printf(" Npart %lld (Ndm %d), z %g, threads %d\n", static_cast<long long>(np), Ndm, z, nt);
    std::printf(" Halos %zu, sources %lld, total SFR %.4g Msun/yr\n", cat.halo.size(),
                static_cast<long long>(cat.Nsource), sfr);
    if (!halos.empty())
      std::printf(" Largest: N %lld M %.4g Msun T %.4g K\n", static_cast<long long>(halos[0].N),
                  halos[0].M, halos[0].T);
    std::printf(" find_halos %.3f s\n", std::chrono::duration<double>(t1 - t0).count());
    std::printf(" Wrote %s (tag for this redshift: halo.%s)\n", out.c_str(), cat.tag.c_str());
  } catch (const std::exception& e) {
    std::fprintf(stderr, "find_halos: error: %s\n", e.what());
    return 1;
  }
  return 0;
}
