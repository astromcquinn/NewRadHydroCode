// On-the-fly IGM diagnostics (igm_diag.hpp).
#include "radray/igm_diag.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>

#include "radray/fft.hpp"

namespace radray {

namespace {

// |k|^2 index of one FFT axis (numpy fftfreq: 0..N/2-1, then -N/2..-1)
inline int freq(int i, int N) { return i <= N / 2 ? i : i - N; }

}  // namespace

void velocity_spectrum(const float* u, int N, double box, double vel_kms, int nbins, int nthreads,
                       VelocitySpectrum& out) {
  const double box_kpc = box * 1e3, kf = 2 * M_PI / box_kpc;
  const std::int64_t n3 = static_cast<std::int64_t>(N) * N * N;
  const int nx = N / 2 + 1, n1p = N + 2;
  out.edges.resize(nbins + 1);
  const double l0 = std::log10(kf * 0.999), l1 = std::log10(M_PI * N / box_kpc * 1.001);
  for (int b = 0; b <= nbins; ++b) out.edges[b] = std::pow(10.0, l0 + (l1 - l0) * b / nbins);
  // the shell of every (ix, iy, iz) is the same for each transform: one table
  std::vector<int> bin(static_cast<std::size_t>(nx) * N * N);
  std::vector<double> kmag(bin.size());
#pragma omp parallel for schedule(static)
  for (int iz = 0; iz < N; ++iz)
    for (int iy = 0; iy < N; ++iy)
      for (int ix = 0; ix < nx; ++ix) {
        const int fy = freq(iy, N), fz = freq(iz, N);
        const double k = kf * std::sqrt(static_cast<double>(ix * ix + fy * fy + fz * fz));
        const std::size_t m = ix + static_cast<std::size_t>(nx) * (iy + static_cast<std::size_t>(N) * iz);
        kmag[m] = k;
        bin[m] = static_cast<int>(std::upper_bound(out.edges.begin(), out.edges.end(), k) - out.edges.begin()) - 1;
      }
  auto weight = [&](int ix) { return (ix == 0 || (N % 2 == 0 && ix == N / 2)) ? 1.0 : 2.0; };
  // sums per k-plane, then in plane order (independent of the thread count)
  std::vector<double> plane(static_cast<std::size_t>(N) * nbins);
  auto shell_sum = [&](const float* a, std::vector<double>& acc) {
    std::fill(plane.begin(), plane.end(), 0.0);
#pragma omp parallel for schedule(static)
    for (int iz = 0; iz < N; ++iz) {
      double* p = &plane[static_cast<std::size_t>(iz) * nbins];
      for (int iy = 0; iy < N; ++iy)
        for (int ix = 0; ix < nx; ++ix) {
          const std::size_t m = ix + static_cast<std::size_t>(nx) * (iy + static_cast<std::size_t>(N) * iz);
          const int b = bin[m];
          if (b < 0 || b >= nbins) continue;
          const std::size_t c = 2 * (ix + static_cast<std::size_t>(nx) * (iy + static_cast<std::size_t>(N) * iz));
          const double re = a[c], im = a[c + 1];
          p[b] += weight(ix) * (re * re + im * im);
        }
    }
    for (int iz = 0; iz < N; ++iz)
      for (int b = 0; b < nbins; ++b) acc[b] += plane[static_cast<std::size_t>(iz) * nbins + b];
  };
  out.cnt.assign(nbins, 0.0);
  out.k.assign(nbins, 0.0);
  {
    std::vector<double> ksum(nbins, 0.0);
    for (int iz = 0; iz < N; ++iz)
      for (int iy = 0; iy < N; ++iy)
        for (int ix = 0; ix < nx; ++ix) {
          const std::size_t m = ix + static_cast<std::size_t>(nx) * (iy + static_cast<std::size_t>(N) * iz);
          const int b = bin[m];
          if (b < 0 || b >= nbins) continue;
          out.cnt[b] += weight(ix);
          ksum[b] += weight(ix) * kmag[m];
        }
    for (int b = 0; b < nbins; ++b) out.k[b] = ksum[b] / std::max(out.cnt[b], 1.0);
  }
  // <rho> for the mass weighting, per plane then in order
  std::vector<double> rp(N, 0.0);
#pragma omp parallel for schedule(static)
  for (int k = 0; k < N; ++k) {
    double s = 0;
    for (std::int64_t c = static_cast<std::int64_t>(k) * N * N; c < static_cast<std::int64_t>(k + 1) * N * N; ++c)
      s += u[c * 5];
    rp[k] = s;
  }
  double rho_avg = 0;
  for (int k = 0; k < N; ++k) rho_avg += rp[k];
  rho_avg /= static_cast<double>(n3);

  std::vector<float> a(static_cast<std::size_t>(n1p) * N * N);
  std::vector<double> sv(nbins, 0.0), sm(nbins, 0.0);
  for (int pass = 0; pass < 2; ++pass)
    for (int comp = 0; comp < 3; ++comp) {
#pragma omp parallel for schedule(static)
      for (int k = 0; k < N; ++k)
        for (int j = 0; j < N; ++j)
          for (int i = 0; i < N; ++i) {
            const std::int64_t c = i + static_cast<std::int64_t>(N) * (j + static_cast<std::int64_t>(N) * k);
            const float* p = u + c * 5;
            double v = static_cast<double>(p[1 + comp]) / p[0] * vel_kms;
            if (pass == 1) v *= std::sqrt(p[0] / rho_avg);
            a[i + static_cast<std::size_t>(n1p) * (j + static_cast<std::size_t>(N) * k)] = static_cast<float>(v);
          }
      fft::fft_3d('f', N, N, N, a.data(), nthreads);
      shell_sum(a.data(), pass == 0 ? sv : sm);
    }
  const double norm = 2.0 * std::pow(static_cast<double>(N), 6);
  out.E.resize(nbins);
  out.Em.resize(nbins);
  out.v2 = out.v2m = 0;
  for (int b = 0; b < nbins; ++b) {
    const double dk = out.edges[b + 1] - out.edges[b];
    out.E[b] = sv[b] / norm / dk;
    out.Em[b] = sm[b] / norm / dk;
    out.v2 += 2 * out.E[b] * dk;
    out.v2m += 2 * out.Em[b] * dk;
  }
}

void write_velocity_spectrum(const std::string& path, double z, int N, double box,
                             const VelocitySpectrum& s) {
  std::FILE* fp = std::fopen(path.c_str(), "w");
  if (!fp) throw std::runtime_error("cannot write " + path);
  std::fprintf(fp,
               "# E(k) of the gas velocity (igm_diag.hpp; tools/analysis/velocity_power.py conventions)\n"
               "# z = %.6f N = %d box = %g Mpc/h\n# <v^2> = %.6e (km/s)^2, mass-weighted <rho v^2>/<rho> = %.6e\n"
               "# k_lo k_hi [h/kpc] k_mean count E [(km/s)^2 kpc/h] E_mass\n",
               z, N, box, s.v2, s.v2m);
  for (std::size_t b = 0; b < s.E.size(); ++b)
    std::fprintf(fp, "%.6e %.6e %.6e %.0f %.6e %.6e\n", s.edges[b], s.edges[b + 1], s.k[b], s.cnt[b], s.E[b],
                 s.Em[b]);
  std::fclose(fp);
}

void tdelta_hist(const float* u, const float* temp, std::int64_t ncell, int N, TDeltaHist& out) {
  if (ncell != static_cast<std::int64_t>(N) * N * N) throw std::runtime_error("tdelta_hist: ncell != N^3");
  const std::int64_t nb = static_cast<std::int64_t>(kTDeltaND) * kTDeltaNT;
  const std::int64_t nplane = static_cast<std::int64_t>(N) * N;
  std::vector<double> rp(N, 0.0);
#pragma omp parallel for schedule(static)
  for (int k = 0; k < N; ++k) {
    double s = 0;
    for (std::int64_t c = k * nplane; c < (k + 1) * nplane; ++c) s += u[c * 5];
    rp[k] = s;
  }
  double rho_avg = 0;
  for (int k = 0; k < N; ++k) rho_avg += rp[k];
  rho_avg /= static_cast<double>(ncell);
  out.count.assign(nb, 0);
  out.mass.assign(nb, 0.0);
  auto idx = [](double x, double x0, int n) {
    const int i = static_cast<int>(std::floor((x - x0) / kTDeltaBin));
    return std::min(std::max(i, 0), n - 1);
  };
  // fixed chunks of planes (not per thread), each with its own histogram,
  // summed in chunk order: the masses do not depend on the thread count
  const int chunk = 16, nchunk = (N + chunk - 1) / chunk;
  std::vector<std::int64_t> cc(static_cast<std::size_t>(nchunk) * nb, 0);
  std::vector<double> cm(cc.size(), 0.0);
#pragma omp parallel for schedule(dynamic)
  for (int q = 0; q < nchunk; ++q) {
    std::int64_t* hc = &cc[static_cast<std::size_t>(q) * nb];
    double* hm = &cm[static_cast<std::size_t>(q) * nb];
    for (int k = q * chunk; k < std::min(N, (q + 1) * chunk); ++k)
      for (std::int64_t c = k * nplane; c < (k + 1) * nplane; ++c) {
        const double d = u[c * 5] / rho_avg;
        const double t = std::max(static_cast<double>(temp[c]), 1e-30);
        const std::int64_t m =
            static_cast<std::int64_t>(idx(std::log10(std::max(d, 1e-30)), kTDeltaD0, kTDeltaND)) * kTDeltaNT +
            idx(std::log10(t), kTDeltaT0, kTDeltaNT);
        hc[m] += 1;
        hm[m] += d;
      }
  }
  for (int q = 0; q < nchunk; ++q)
    for (std::int64_t m = 0; m < nb; ++m) {
      out.count[m] += cc[static_cast<std::size_t>(q) * nb + m];
      out.mass[m] += cm[static_cast<std::size_t>(q) * nb + m];
    }
}

void write_tdelta(const std::string& path, double z, const TDeltaHist& h) {
  std::FILE* fp = std::fopen(path.c_str(), "w");
  if (!fp) throw std::runtime_error("cannot write " + path);
  std::fprintf(fp,
               "# T-Delta histogram of the gas (igm_diag.hpp)\n# z = %.6f\n"
               "# bins: log10 Delta from %g, log10 T [K] from %g, width %g; %d x %d (edge bins hold the overflow)\n"
               "# iDelta iT cells mass (sum of Delta), non-empty bins only\n",
               z, kTDeltaD0, kTDeltaT0, kTDeltaBin, kTDeltaND, kTDeltaNT);
  for (int i = 0; i < kTDeltaND; ++i)
    for (int j = 0; j < kTDeltaNT; ++j) {
      const std::int64_t m = static_cast<std::int64_t>(i) * kTDeltaNT + j;
      if (h.count[m] > 0)
        std::fprintf(fp, "%d %d %lld %.9e\n", i, j, static_cast<long long>(h.count[m]), h.mass[m]);
    }
  std::fclose(fp);
}

bool diag_due(double z0, double dz, double z_prev, double z) {
  if (dz <= 0) return false;
  auto n = [&](double zz) {
    if (zz > z0 + 1e-9) return -1.0;
    return std::floor((z0 - zz) / dz + 1e-9);
  };
  return n(z) > n(z_prev);
}

}  // namespace radray
