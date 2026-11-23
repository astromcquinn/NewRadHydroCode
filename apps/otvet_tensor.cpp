// The alternative OTVET tensors (otvet_closure.hpp) on a snapshot with a
// tensor dump (radhydro tensor_dump = 1, rt_method = rays): the a priori
// check of docs/otvet_closure.md, against the rays' true tensor and the
// numpy prototypes of tools/analysis/tensor_compare.py.
//
// usage: otvet_tensor STEPDIR OUTDIR --box BOX [--set key=value ...] [--repeat n]
//
// --repeat n computes the tensor n times and prints each wall time (the
// FFT kernel spectra are cached after the first call, as in a run).
//
// Reads STEPDIR/{sources.txt, rhog.bin, fHI.bin, tensor.bin}; the
// solution's photons N are the rays' S0 (Gamma-weighted, the units of the
// dump), the HI column follows calc_rt_1. Writes OUTDIR/tensor.bin in the
// dump's layout (10 float32 per cell: S0, S1 = 0, S2 = S0 h), with
// otvet_split = 1 applied to the Gamma-weighted S0
// (P = S0 I/3 + min(S0, N_dir)(h - I/3)), for tensor_compare.py --moment.
//
// --set otvet_background=ga01 (docs/otvet_background.md section 4c): the
// pressure-only form of the ga01 background, the a priori proxy for it,
// P = S0 h + min(S0, Lambda <S0>)(I/3 - h) with the thin h, <S0> the box
// mean and Lambda the source-spectrum Gamma-weighted mean of Lambda_g =
// exp(-0.611 N <tau_cell,g>) (tau_cell,g from rhog, fHI, fHeI and fHeII as
// calc_rt_1), or otvet_background_fix. Also reads fHeI.bin, fHeII.bin.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "radray/config.hpp"
#include "radray/constants.hpp"
#include "radray/cosmology.hpp"
#include "radray/otvet_closure.hpp"
#include "radray/parallel.hpp"
#include "radray/rt_moments.hpp"
#include "radray/rt_rates.hpp"

using namespace radray;

namespace {

std::vector<float> read_floats(const std::string& path, std::size_t n) {
  std::ifstream f(path, std::ios::binary);
  if (!f) throw std::runtime_error("cannot read " + path);
  std::vector<float> v(n);
  f.read(reinterpret_cast<char*>(v.data()), static_cast<std::streamsize>(n * sizeof(float)));
  if (!f) throw std::runtime_error("short file " + path);
  return v;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc < 3) {
      std::cerr << "usage: otvet_tensor STEPDIR OUTDIR --box BOX [--set key=value ...]\n";
      return 2;
    }
    const std::string dir = argv[1], outdir = argv[2];
    double box = 0;
    Config cfg;
    cfg.rt_method = RTMethod::OTVET;
    std::vector<std::pair<std::string, std::string>> sets;
    int method = 0;
    int repeat = 1;
    for (int a = 3; a < argc; ++a) {
      const std::string s = argv[a];
      if (s == "--box" && a + 1 < argc) {
        box = std::stod(argv[++a]);
      } else if (s == "--repeat" && a + 1 < argc) {
        repeat = std::max(1, std::stoi(argv[++a]));
      } else if (s == "--set" && a + 1 < argc) {
        const std::string kv = argv[++a];
        const auto eq = kv.find('=');
        if (eq == std::string::npos) throw std::runtime_error("--set key=value");
        sets.push_back({kv.substr(0, eq), kv.substr(eq + 1)});
      } else {
        throw std::runtime_error("unknown argument " + s);
      }
    }
    if (!(box > 0)) throw std::runtime_error("--box is required");

    // sources.txt: header qfac, sigma18, redshift z N n; rows i j k Luv
    OtvetClosureIn in;
    double z = 0;
    {
      std::ifstream f(dir + "/sources.txt");
      if (!f) throw std::runtime_error("cannot read " + dir + "/sources.txt");
      std::string line;
      while (std::getline(f, line)) {
        std::istringstream is(line);
        if (line.rfind("#", 0) == 0) {
          std::string hash, key;
          is >> hash >> key;
          if (key == "qfac")
            for (int g = 0; g < kNfreq; ++g) is >> in.q[g];
          else if (key == "sigma18")
            for (int g = 0; g < kNfreq; ++g) is >> in.sigma18[g];
          else if (key == "redshift") {
            std::string nkey;
            is >> z >> nkey >> in.N;
          }
          continue;
        }
        int i, j, k;
        double L;
        if (is >> i >> j >> k >> L) in.src.push_back({i - 1, j - 1, k - 1, L});
      }
    }
    if (in.N <= 0) throw std::runtime_error("sources.txt: no '# redshift z N n' line");
    const int N = in.N;
    cfg = Config::make(N, 64, -1, box);
    cfg.rt_method = RTMethod::OTVET;
    for (const auto& kv : sets)
      if (!set_otvet_closure_param(cfg, kv.first, kv.second)) {
        if (kv.first == "otvet_tensor")
          method = std::stoi(kv.second);
        else
          throw std::runtime_error("unknown key " + kv.first);
      }
    check_otvet_closure_param(cfg);
    const OtvetClosureParams p = otvet_closure_params(cfg);

    const std::size_t nc = static_cast<std::size_t>(N) * N * N;
    const std::vector<float> rhog = read_floats(dir + "/rhog.bin", nc);
    const std::vector<float> fHI = read_floats(dir + "/fHI.bin", nc);
    const std::vector<float> S = read_floats(dir + "/tensor.bin", 10 * nc);
    std::vector<double> nw(nc);
    for (std::size_t c = 0; c < nc; ++c) nw[c] = S[10 * c];
    in.nw = nw.data();

    SimClock clk;
    clk.redshift = z;
    clk.scalefactor = 1 / (1 + z);
    Units u;
    calc_sim_var(cfg, clk, u);
    std::vector<float> col(nc);
    for (std::size_t c = 0; c < nc; ++c) {
      const double nH = static_cast<double>(rhog[c]) * u.rhob_avg * cp::XHy / pc::mH_cgs;
      col[c] = static_cast<float>(nH * static_cast<double>(fHI[c]) * u.len_rt * 1e-18);
    }
    in.col = col.data();
    std::vector<float> xion(nc);
    for (std::size_t c = 0; c < nc; ++c) xion[c] = 1.0f - fHI[c];
    in.xion = xion.data();

    // ga01's pressure-only form: the background share Lambda <S0> per cell
    double bg_share = -1;
    if (cfg.otvet_background) {
      const std::vector<float> fHeI = read_floats(dir + "/fHeI.bin", nc);
      const std::vector<float> fHeII = read_floats(dir + "/fHeII.bin", nc);
      const SigmaTable sig = photoion_sigma(cp::freq_init);
      double tau[kNfreq] = {}, lam[kNfreq], s0 = 0;
      for (std::size_t c = 0; c < nc; ++c) {
        const double rg = static_cast<double>(rhog[c]) * u.rhob_avg;
        const double nH = rg * cp::XHy / pc::mH_cgs, nHe = rg * cp::YHe / pc::mHe_cgs;
        const double a = nH * static_cast<double>(fHI[c]), b = nHe * static_cast<double>(fHeI[c]),
                     d = nHe * static_cast<double>(fHeII[c]);
        for (int g = 0; g < kNfreq; ++g)
          tau[g] += static_cast<double>(
              static_cast<float>((((0.0 + a * sig[0][g]) + b * sig[1][g]) + d * sig[2][g]) * u.len_rt));
        s0 += nw[c];
      }
      double wl = 0, ws = 0;
      std::printf("ga01: Lambda_g");
      for (int g = 0; g < kNfreq; ++g) {
        tau[g] /= static_cast<double>(nc);
        lam[g] = cfg.otvet_background_fix >= 0 ? cfg.otvet_background_fix : std::exp(-kGa01Path * N * tau[g]);
        wl += in.q[g] * in.sigma18[g] * lam[g];
        ws += in.q[g] * in.sigma18[g];
        std::printf(" %.4g", lam[g]);
      }
      const double L = ws > 0 ? wl / ws : 0;
      bg_share = L * s0 / static_cast<double>(nc);
      std::printf("; Gamma-weighted Lambda %.4g, <S0> %.4g\n", L, s0 / static_cast<double>(nc));
    }

    OtvetClosureOut out;
    for (int it = 0; it < repeat; ++it) {
      otvet_closure_tensor(p, method, in, out);
      if (repeat > 1) std::printf("call %d: %.2f s\n", it + 1, out.seconds);
    }
    double qw = 0;
    for (int g = 0; g < kNfreq; ++g) qw += in.q[g] * in.sigma18[g];
    std::printf("z = %.3f N = %d sources = %zu: %s; lambda %.1f lambda_x %.1f N_dir/N %.3f, %.1f s\n",
                z, N, in.src.size(), otvet_closure_summary(cfg).c_str(), out.lambda, out.lam_x,
                out.ratio, out.seconds);

    std::vector<float> o(10 * nc, 0.0f);
    for (std::size_t c = 0; c < nc; ++c) {
      const double S0 = nw[c];
      const float* h = &out.h[6 * c];
      double P[6];
      if (bg_share >= 0 && S0 > 0) {
        const double m = std::min(S0, bg_share);
        for (int k = 0; k < 3; ++k) P[k] = S0 * static_cast<double>(h[k]) + m * (1.0 / 3 - static_cast<double>(h[k]));
        for (int k = 3; k < 6; ++k) P[k] = (S0 - m) * static_cast<double>(h[k]);
      } else if (!out.ndir.empty() && S0 > 0) {
        const double m = std::min(S0, qw * static_cast<double>(out.ndir[c]));
        for (int k = 0; k < 3; ++k) P[k] = (S0 - m) / 3 + m * static_cast<double>(h[k]);
        for (int k = 3; k < 6; ++k) P[k] = m * static_cast<double>(h[k]);
      } else {
        const double w = S0 > 0 ? S0 : 1.0;
        for (int k = 0; k < 6; ++k) P[k] = w * static_cast<double>(h[k]);
      }
      o[10 * c] = static_cast<float>(S0);
      for (int k = 0; k < 6; ++k) o[10 * c + 4 + k] = static_cast<float>(P[k]);
    }
    std::string mk = "mkdir -p '" + outdir + "'";
    if (std::system(mk.c_str()) != 0) throw std::runtime_error("cannot create " + outdir);
    std::ofstream f(outdir + "/tensor.bin", std::ios::binary);
    f.write(reinterpret_cast<const char*>(o.data()), static_cast<std::streamsize>(o.size() * sizeof(float)));
    if (!f) throw std::runtime_error("cannot write " + outdir + "/tensor.bin");
    {
      // a copy of the source list with the header (tensor_compare reads N, z)
      std::ifstream a(dir + "/sources.txt");
      std::ofstream b(outdir + "/sources.txt");
      b << a.rdbuf();
    }
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "otvet_tensor: " << e.what() << "\n";
    return 1;
  }
}
