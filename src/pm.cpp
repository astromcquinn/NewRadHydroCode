// Port of pm_tools.f90 (particle-mesh gravity for dark matter and gas) and
// the DM/gravity parts of radhydro.f90, hydro_tools.f90, cosmo_tools.f90 and
// io_tools.f90. Every expression keeps the Fortran operation order and
// precision; see docs/port_notes_pm.md.
#include "radray/pm.hpp"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>
#include <cstdio>
#include <fstream>
#include <ostream>
#include <stdexcept>

#include "radray/constants.hpp"
#include "radray/fft.hpp"
#include "radray/fortran.hpp"
#include "radray/parallel.hpp"

namespace radray {
namespace pm {

using namespace pc;
using namespace cp;

namespace {

// pm_tools.f90: cic_mass_assignment, local parameter
//   mass = (1-omegab/omegam)*N_dm_gas_c
inline double dm_mass(const PMConfig& pc) { return (1 - omegab / omegam) * pc.N_dm_gas_c; }

// CIC weights shared by cic_mass_assignment and accelerate_dm:
//   x = mod(xv_dm(1,ip)*x_dm_gas - 0.5 + ngasr, ngasr), i1 = 1+int(x), ...
struct Cic {
  int i1, i2, j1, j2, k1, k2;
  double dx1, dx2, dy1, dy2, dz1, dz2;
};
inline Cic cic_weights(const PMConfig& pc, const float* xv) {
  const int N = pc.Ngas;
  Cic c;
  const double x = f::mod(static_cast<double>(xv[0]) * pc.x_dm_gas - 0.5 + pc.ngasr, pc.ngasr);
  const double y = f::mod(static_cast<double>(xv[1]) * pc.x_dm_gas - 0.5 + pc.ngasr, pc.ngasr);
  const double z = f::mod(static_cast<double>(xv[2]) * pc.x_dm_gas - 0.5 + pc.ngasr, pc.ngasr);
  c.i1 = 1 + f::int4(x);
  c.i2 = 1 + f::mod(c.i1, N);
  c.dx1 = c.i1 - x;
  c.dx2 = 1 - c.dx1;
  c.j1 = 1 + f::int4(y);
  c.j2 = 1 + f::mod(c.j1, N);
  c.dy1 = c.j1 - y;
  c.dy2 = 1 - c.dy1;
  c.k1 = 1 + f::int4(z);
  c.k2 = 1 + f::mod(c.k1, N);
  c.dz1 = c.k1 - z;
  c.dz2 = 1 - c.dz1;
  return c;
}

// calc_kick: function force(i,j,k). The whole finite difference is a real(4)
// expression (2./3 and 1./12 are real(4) literals, density is real(4));
// only the result is stored into the real(8) function value.
inline void force(const PMState& s, int i, int j, int k, double F[3]) {
  const int N = s.pc.Ngas;
  const int im2 = 1 + f::mod(i - 3 + N, N);
  const int im1 = 1 + f::mod(i - 2 + N, N);
  const int ip1 = 1 + f::mod(i, N);
  const int ip2 = 1 + f::mod(i + 1 + N, N);
  const int jm2 = 1 + f::mod(j - 3 + N, N);
  const int jm1 = 1 + f::mod(j - 2 + N, N);
  const int jp1 = 1 + f::mod(j, N);
  const int jp2 = 1 + f::mod(j + 1 + N, N);
  const int km2 = 1 + f::mod(k - 3 + N, N);
  const int km1 = 1 + f::mod(k - 2 + N, N);
  const int kp1 = 1 + f::mod(k, N);
  const int kp2 = 1 + f::mod(k + 1 + N, N);
  const float c23 = 2.0f / 3.0f;
  const float c112 = 1.0f / 12.0f;
  const float* d = s.density.data();
  const float fx = -(c23 * (d[s.dcell(ip1, j, k)] - d[s.dcell(im1, j, k)])) +
                   c112 * (d[s.dcell(ip2, j, k)] - d[s.dcell(im2, j, k)]);
  const float fy = -(c23 * (d[s.dcell(i, jp1, k)] - d[s.dcell(i, jm1, k)])) +
                   c112 * (d[s.dcell(i, jp2, k)] - d[s.dcell(i, jm2, k)]);
  const float fz = -(c23 * (d[s.dcell(i, j, kp1)] - d[s.dcell(i, j, km1)])) +
                   c112 * (d[s.dcell(i, j, kp2)] - d[s.dcell(i, j, km2)]);
  F[0] = fx;
  F[1] = fy;
  F[2] = fz;
}

// general_tools.f90: sinc
inline double sinc(double x) {
  if (std::fabs(x) > 1e-6) return std::sin(x) / x;
  return 1 - f::pow2(x) / 6;
}

}  // namespace

PMConfig PMConfig::make(int Ngas, int Ndm, bool mf) {
  if (Ngas <= 0 || Ndm <= 0 || Ngas % Ndm != 0)
    throw std::runtime_error("PMConfig: Ngas must be a positive multiple of Ndm");
  if (Ngas % 2 != 0) throw std::runtime_error("PMConfig: Ngas must be even");
  PMConfig c;
  c.Ngas = Ngas;
  c.Hgas = Ngas / 2;
  c.Ndm = Ndm;
  c.Npart = static_cast<std::int64_t>(Ndm) * Ndm * Ndm;
  c.Ncell = static_cast<std::int64_t>(Ngas) * Ngas * Ngas;
  c.ngasr = Ngas;
  c.ndmr = Ndm;
  c.N_dm_gas = Ngas / Ndm;
  c.N_dm_gas_c = c.N_dm_gas * c.N_dm_gas * c.N_dm_gas;
  c.x_dm_gas = c.N_dm_gas;
  c.mf = mf;
  return c;
}

// pm_tools.f90: init_pm
void init_pm(PMState& s, const PMConfig& pc, const Domains& dom) {
  s.pc = pc;
  s.dom = dom;
  const std::int64_t N = pc.Ngas;
  // npt = Npart/Ncpu + min(1,mod(Npart,Ncpu)) with Ncpu = 1
  s.dm_cpu[0] = 1;
  s.dm_cpu[1] = pc.Npart;
  s.xv.assign(6 * pc.Npart, 0.0f);
  s.ll.assign(pc.Npart + 1, 0);
  s.hoc.assign(N * N * N, 0);
  s.hoc_cpu.assign(N * N, 0);
  s.density.assign((N + 2) * N * N, 0.0f);
}

// pm_tools.f90: calc_drift
//   xv_dm(1:3,ip) = mod(xv_dm(1:3,ip) + xv_dm(4:6,ip)*dt_sim/2 + ndmr,ndmr)
void calc_drift(PMState& s, double dt_sim) {
  const double ndmr = s.pc.ndmr;
#pragma omp parallel for schedule(static) if (par::worth(s.pc.Npart))
  for (std::int64_t ip = s.dm_cpu[0]; ip <= s.dm_cpu[1]; ++ip) {
    float* p = s.xv_p(ip);
    for (int c = 0; c < 3; ++c) {
      const double t = static_cast<double>(p[c]) + static_cast<double>(p[3 + c]) * dt_sim / 2 + ndmr;
      p[c] = static_cast<float>(f::mod(t, ndmr));
    }
  }
}

// pm_tools.f90: link_dm (Ncpu = 1)
void link_dm(PMState& s) {
  const int N = s.pc.Ngas;
  const double xdg = s.pc.x_dm_gas;
  std::fill(s.hoc_cpu.begin(), s.hoc_cpu.end(), 0);
  auto hc = [&](int j, int k) -> std::int64_t& {
    return s.hoc_cpu[(j - 1) + static_cast<std::int64_t>(N) * (k - 1)];
  };
  for (std::int64_t ip = s.dm_cpu[0]; ip <= s.dm_cpu[1]; ++ip) {
    const float* p = s.xv_p(ip);
    const int j = 1 + f::mod(f::int4(static_cast<double>(p[1]) * xdg), N);
    const int k = 1 + f::mod(f::int4(static_cast<double>(p[2]) * xdg), N);
    s.ll[ip] = hc(j, k);
    hc(j, k) = ip;
  }
  // the rows (j,k) are independent (each relinks only its own particles)
#pragma omp parallel for schedule(static) if (par::worth(s.pc.Ncell))
  for (int k = 1; k <= N; ++k) {
    for (int j = 1; j <= N; ++j) {
      for (int i = 1; i <= N; ++i) s.hoc[s.gcell(i, j, k)] = 0;
      // loop over cpus (Ncpu..1): one
      std::int64_t ip = hc(j, k);
      while (ip > 0) {
        const int i = 1 + f::mod(f::int4(static_cast<double>(s.xv_p(ip)[0]) * xdg), N);
        const std::int64_t next = s.ll[ip];
        s.ll[ip] = s.hoc[s.gcell(i, j, k)];
        s.hoc[s.gcell(i, j, k)] = ip;
        ip = next;
      }
    }
  }
}

// pm_tools.f90: cic_mass_assignment(indx)
void cic_mass_assignment(PMState& s, const std::array<int, 6>& b) {
  const double mass = dm_mass(s.pc);
  float* d = s.density.data();
  auto add = [&](int i, int j, int k, double w) {
    float& r = d[s.dcell(i, j, k)];
    r = static_cast<float>(static_cast<double>(r) + w * mass);
  };
  for (int k = b[4]; k <= b[5]; ++k)
    for (int j = b[2]; j <= b[3]; ++j)
      for (int i = b[0]; i <= b[1]; ++i) {
        std::int64_t ip = s.hoc[s.gcell(i, j, k)];
        while (ip > 0) {
          const Cic c = cic_weights(s.pc, s.xv_p(ip));
          // density(..) = density(..) + dx*dy*dz*mass
          add(c.i1, c.j1, c.k1, c.dx1 * c.dy1 * c.dz1);
          add(c.i2, c.j1, c.k1, c.dx2 * c.dy1 * c.dz1);
          add(c.i1, c.j2, c.k1, c.dx1 * c.dy2 * c.dz1);
          add(c.i2, c.j2, c.k1, c.dx2 * c.dy2 * c.dz1);
          add(c.i1, c.j1, c.k2, c.dx1 * c.dy1 * c.dz2);
          add(c.i2, c.j1, c.k2, c.dx2 * c.dy1 * c.dz2);
          add(c.i1, c.j2, c.k2, c.dx1 * c.dy2 * c.dz2);
          add(c.i2, c.j2, c.k2, c.dx2 * c.dy2 * c.dz2);
          ip = s.ll[ip];
        }
      }
}

// cic_mass_assignment over all gas domains in the order 1..Ndomain (what
// set_domain_gas hands out with one thread), in parallel: thread t owns the
// target planes [ka,kb], walks the source cells in exactly the serial order
// (restricted to the planes that can deposit into its slab) and adds only
// into its own planes. Every density cell therefore receives the same
// addends in the same order as in the serial loop: bit-identical for any
// thread count.
void cic_mass_assignment_all(PMState& s) {
  const int N = s.pc.Ngas;
  const int nt = par::worth(s.pc.Ncell) ? std::min(par::max_threads(), N) : 1;
  if (nt <= 1) {
    for (const auto& b : s.dom.box) cic_mass_assignment(s, b);
    return;
  }
  const double mass = dm_mass(s.pc);
  float* d = s.density.data();
#pragma omp parallel num_threads(nt)
  {
    const int t = par::thread_id(), nth = par::num_threads();
    const int ka = 1 + static_cast<int>(static_cast<std::int64_t>(t) * N / nth);
    const int kb = static_cast<int>(static_cast<std::int64_t>(t + 1) * N / nth);
    std::vector<char> src(N + 1, 0), own(N + 1, 0);
    for (int k = ka; k <= kb; ++k) own[k] = 1;
    // CIC from source plane k reaches planes k-1..k+1; margin of 2
    for (int k = ka - 2; k <= kb + 2; ++k) src[1 + ((k - 1) % N + N) % N] = 1;
    auto add = [&](int i, int j, int k, double w) {
      if (!own[k]) return;
      float& r = d[s.dcell(i, j, k)];
      r = static_cast<float>(static_cast<double>(r) + w * mass);
    };
    for (const auto& b : s.dom.box)
      for (int k = b[4]; k <= b[5]; ++k) {
        if (!src[k]) continue;
        for (int j = b[2]; j <= b[3]; ++j)
          for (int i = b[0]; i <= b[1]; ++i) {
            std::int64_t ip = s.hoc[s.gcell(i, j, k)];
            while (ip > 0) {
              const Cic c = cic_weights(s.pc, s.xv_p(ip));
              add(c.i1, c.j1, c.k1, c.dx1 * c.dy1 * c.dz1);
              add(c.i2, c.j1, c.k1, c.dx2 * c.dy1 * c.dz1);
              add(c.i1, c.j2, c.k1, c.dx1 * c.dy2 * c.dz1);
              add(c.i2, c.j2, c.k1, c.dx2 * c.dy2 * c.dz1);
              add(c.i1, c.j1, c.k2, c.dx1 * c.dy1 * c.dz2);
              add(c.i2, c.j1, c.k2, c.dx2 * c.dy1 * c.dz2);
              add(c.i1, c.j2, c.k2, c.dx1 * c.dy2 * c.dz2);
              add(c.i2, c.j2, c.k2, c.dx2 * c.dy2 * c.dz2);
              ip = s.ll[ip];
            }
          }
      }
  }
}

// pm_tools.f90: calc_density_field
void calc_density_field(PMState& s, const float* u_gas) {
  const int N = s.pc.Ngas;
#pragma omp parallel for schedule(static) if (par::worth(s.pc.Ncell))
  for (int k = 1; k <= N; ++k)
    for (int j = 1; j <= N; ++j) {
      for (int i = 1; i <= N; ++i) s.density[s.dcell(i, j, k)] = u_gas[s.gcell(i, j, k) * 5];
      for (int i = N + 1; i <= N + 2; ++i) s.density[s.dcell(i, j, k)] = 0;
    }
  // set_domain_gas/end_domain_gas hand out domains 1..Ndomain in order
  cic_mass_assignment_all(s);
}

// pm_tools.f90: calc_gravity_field
void calc_gravity_field(PMState& s, double scalefactor) {
  const int N = s.pc.Ngas;
  const int H = s.pc.Hgas;
  const double Ak = 2 * pi / s.pc.ngasr;
  // Aphi = -4*pi*scalefactor/(6*pi)  (unary minus applies to the quotient)
  const double Aphi = -(4 * pi * scalefactor / (6 * pi));

  fft::fft_3d('f', N, N, N, s.density.data());

#pragma omp parallel for schedule(static) if (par::worth(s.pc.Ncell))
  for (int k = 1; k <= N; ++k) {
    const double kz = (k <= H + 1) ? Ak * (k - 1) : Ak * (k - 1 - N);
    const double kc = 2 * std::sin(kz / 2);
    for (int j = 1; j <= N; ++j) {
      const double ky = (j <= H + 1) ? Ak * (j - 1) : Ak * (j - 1 - N);
      const double kb = 2 * std::sin(ky / 2);
      for (int i = 1; i <= N + 2; i += 2) {
        const double kx = Ak * ((i - 1) / 2);
        const double ka = 2 * std::sin(kx / 2);
        const double ksq = f::pow2(kx) + f::pow2(ky) + f::pow2(kz);
        double wphi;
        if (ksq > 0)
          wphi = Aphi / (f::pow2(ka) + f::pow2(kb) + f::pow2(kc));
        else
          wphi = 0;
        float* d = &s.density[s.dcell(i, j, k)];
        d[0] = static_cast<float>(static_cast<double>(d[0]) * wphi);
        d[1] = static_cast<float>(static_cast<double>(d[1]) * wphi);
      }
    }
  }

  fft::fft_3d('b', N, N, N, s.density.data());
}

// pm_tools.f90: calc_kick (accelerate_dm for all slabs, then accelerate_gas)
void calc_kick(PMState& s, float* u_gas, float* v_gas, SimClock& clk, const Units& u) {
  const double dt_sim = clk.dt_sim;
  const int N = s.pc.Ngas;
  const double xdg = s.pc.x_dm_gas;
  double amax = 0, arms = 0, vmax = 0, vrms = 0;

  // accelerate_dm(k). Particles are independent; the double sums arms and
  // vrms are recorded per plane and added in plane order afterwards (same
  // addends, same order as the serial loop), amax/vmax are maxima.
  const int nth = par::worth(s.pc.Ncell) ? par::max_threads() : 1;
  const int batch = nth > 1 ? std::min(N, 4 * nth) : 1;
  std::vector<std::vector<double>> add_a(batch), add_v(batch);
  std::vector<double> amax_t(nth > 1 ? batch : 1, 0.0), vmax_t(nth > 1 ? batch : 1, 0.0);
  for (int k0 = 1; k0 <= N; k0 += batch) {
  const int nb = std::min(batch, N - k0 + 1);
#pragma omp parallel for schedule(dynamic, 1) if (nth > 1)
  for (int b = 0; b < nb; ++b) {
    const int k = k0 + b;
    std::vector<double>& la = add_a[b];
    std::vector<double>& lv = add_v[b];
    double amax = 0, vmax = 0;
    for (int j = 1; j <= N; ++j)
      for (int i = 1; i <= N; ++i) {
        std::int64_t ip = s.hoc[s.gcell(i, j, k)];
        while (ip > 0) {
          float* p = s.xv_p(ip);
          const Cic c = cic_weights(s.pc, p);
          double fv[3] = {0, 0, 0};
          double F[3];
          auto acc = [&](int ii, int jj, int kk, double wx, double wy, double wz) {
            force(s, ii, jj, kk, F);
            for (int m = 0; m < 3; ++m) fv[m] = fv[m] + F[m] * wx * wy * wz;
          };
          acc(c.i1, c.j1, c.k1, c.dx1, c.dy1, c.dz1);
          acc(c.i2, c.j1, c.k1, c.dx2, c.dy1, c.dz1);
          acc(c.i1, c.j2, c.k1, c.dx1, c.dy2, c.dz1);
          acc(c.i2, c.j2, c.k1, c.dx2, c.dy2, c.dz1);
          acc(c.i1, c.j1, c.k2, c.dx1, c.dy1, c.dz2);
          acc(c.i2, c.j1, c.k2, c.dx2, c.dy1, c.dz2);
          acc(c.i1, c.j2, c.k2, c.dx1, c.dy2, c.dz2);
          acc(c.i2, c.j2, c.k2, c.dx2, c.dy2, c.dz2);

          // xv_dm(4:6,ip) = xv_dm(4:6,ip) + (f/x_dm_gas)*dt_sim
          for (int m = 0; m < 3; ++m)
            p[3 + m] = static_cast<float>(static_cast<double>(p[3 + m]) + (fv[m] / xdg) * dt_sim);

          // arms = arms + sum(f**2); amax = max(amax,maxval(abs(f)))
          la.push_back((0.0 + fv[0] * fv[0]) + fv[1] * fv[1] + fv[2] * fv[2]);
          amax = f::max(amax, f::max(f::max(std::fabs(fv[0]), std::fabs(fv[1])), std::fabs(fv[2])));
          // vrms = vrms + sum(xv_dm(4:6,ip)**2)   (real(4) sum)
          const float vs = ((0.0f + p[3] * p[3]) + p[4] * p[4]) + p[5] * p[5];
          lv.push_back(static_cast<double>(vs));
          const float va = f::max(f::max(std::fabs(p[3]), std::fabs(p[4])), std::fabs(p[5]));
          vmax = f::max(vmax, static_cast<double>(va));

          ip = s.ll[ip];
        }
      }
    amax_t[nth > 1 ? b : 0] = f::max(amax_t[nth > 1 ? b : 0], amax);
    vmax_t[nth > 1 ? b : 0] = f::max(vmax_t[nth > 1 ? b : 0], vmax);
  }
  for (int b = 0; b < nb; ++b) {
    for (double a : add_a[b]) arms = arms + a;
    for (double v : add_v[b]) vrms = vrms + v;
    add_a[b].clear();
    add_v[b].clear();
  }
  }
  for (double a : amax_t) amax = f::max(amax, a);
  for (double v : vmax_t) vmax = f::max(vmax, v);

  // accelerate_gas(k): cells are independent
#pragma omp parallel for schedule(static) if (nth > 1)
  for (int k = 1; k <= N; ++k)
    for (int j = 1; j <= N; ++j)
      for (int i = 1; i <= N; ++i) {
        double F[3];
        force(s, i, j, k, F);
        const std::int64_t c = s.gcell(i, j, k);
        if (s.pc.mf) {
          // v_gas(:,i,j,k) = v_gas(:,i,j,k) + force(i,j,k)*dt_sim
          float* v = &v_gas[c * 3];
          for (int m = 0; m < 3; ++m)
            v[m] = static_cast<float>(static_cast<double>(v[m]) + F[m] * dt_sim);
        } else {
          float* ug = &u_gas[c * 5];
          const double rho = f::max(static_cast<double>(ug[0]), u.Dgas_min);
          double dv[3];
          for (int m = 0; m < 3; ++m) dv[m] = static_cast<double>(ug[1 + m]) / rho;
          double eth = static_cast<double>(ug[4]) -
                       rho * ((0.0 + dv[0] * dv[0]) + dv[1] * dv[1] + dv[2] * dv[2]) / 2;
          eth = f::min(f::max(eth, rho * u.Tgas_min), rho * u.Tgas_max);
          ug[0] = static_cast<float>(rho);
          for (int m = 0; m < 3; ++m) ug[1 + m] = static_cast<float>(rho * (dv[m] + F[m] * dt_sim));
          const float s2 = ((0.0f + ug[1] * ug[1]) + ug[2] * ug[2]) + ug[3] * ug[3];
          ug[4] = static_cast<float>(eth + static_cast<double>(s2) / rho / 2);
        }
      }

  // Update global variables in PM grid units
  const double npart = static_cast<double>(s.pc.Npart);
  clk.amax_dm = amax;
  clk.arms_dm = std::sqrt(arms / npart);
  clk.vmax_dm = vmax;
  clk.vrms_dm = std::sqrt(vrms / npart);
}

// pm_tools.f90: particle_mesh
void particle_mesh(PMState& s, float* u_gas, float* v_gas, SimClock& clk, const Units& u) {
  calc_drift(s, clk.dt_sim);
  link_dm(s);
  calc_density_field(s, u_gas);
  calc_gravity_field(s, clk.scalefactor);
  calc_kick(s, u_gas, v_gas, clk, u);
  calc_drift(s, clk.dt_sim);
}

// cosmo_tools.f90: calc_pk_3d
void calc_pk_3d(int nc, const float* field, std::vector<double>& Dk, bool deconvolve,
                double box) {
  std::vector<double> pk(5 * static_cast<std::size_t>(nc), 0.0);
  auto P = [&](int c, int k) -> double& { return pk[(k - 1) * 5 + (c - 1)]; };
  const int hc = nc / 2;
  const double ncr = nc;
  const double ncube = f::pow3(ncr);
  const std::int64_t ld = nc + 2;
  for (int k = 1; k <= nc; ++k) {
    const double kz = (k <= hc + 1) ? 2 * pi / ncr * (k - 1) : 2 * pi / ncr * (k - 1 - nc);
    for (int j = 1; j <= nc; ++j) {
      const double ky = (j <= hc + 1) ? 2 * pi / ncr * (j - 1) : 2 * pi / ncr * (j - 1 - nc);
      for (int i = 1; i <= nc + 2; i += 2) {
        const double kx = 2 * pi / ncr * ((i - 1) / 2);
        double kr = std::sqrt(f::pow2(kx) + f::pow2(ky) + f::pow2(kz));
        if (kr > 0) {
          double wcic;
          if (deconvolve)
            wcic = f::pow2(sinc(kx / 2) * sinc(ky / 2) * sinc(kz / 2));
          else
            wcic = 1;
          const float* fp = field + (i - 1) + ld * ((j - 1) + static_cast<std::int64_t>(nc) * (k - 1));
          const double a0 = static_cast<double>(fp[0]) / ncube;
          const double a1 = static_cast<double>(fp[1]) / ncube;
          const double p = (0.0 + a0 * a0) + a1 * a1;
          const double pcor = p / f::pow2(wcic);
          kr = kr / (2 * pi / ncr);
          const int kk = 1 + f::nint(kr);
          P(1, kk) = P(1, kk) + 1;
          P(2, kk) = P(2, kk) + p;
          P(3, kk) = P(3, kk) + f::pow2(p);
          P(4, kk) = P(4, kk) + pcor;
          P(5, kk) = P(5, kk) + f::pow2(pcor);
        }
      }
    }
  }
  Dk.assign(5 * static_cast<std::size_t>(nc), 0.0);
  for (int k = 1; k <= nc; ++k) {
    if (P(1, k) > 0) {
      for (int c = 2; c <= 5; ++c) P(c, k) = P(c, k) / P(1, k);
      P(3, k) = f::max((P(3, k) - f::pow2(P(2, k))) / f::max(P(1, k) - 1, 1.0), 0.0);
      P(3, k) = std::sqrt(P(3, k));
      P(5, k) = f::max((P(5, k) - f::pow2(P(4, k))) / f::max(P(1, k) - 1, 1.0), 0.0);
      P(5, k) = std::sqrt(P(5, k));
    } else {
      for (int c = 2; c <= 5; ++c) P(c, k) = 0;
    }
    const double km1 = k - 1.0;
    Dk[(k - 1) * 5 + 0] = 2 * pi * km1 / box;
    for (int c = 2; c <= 5; ++c) Dk[(k - 1) * 5 + (c - 1)] = 4 * pi * f::pow3(km1) * P(c, k);
  }
}

// pm_tools.f90: calc_power_spectrum (computation; file name/IO in the caller)
void calc_power_spectrum(PMState& s, const float* u_gas, double box,
                         std::vector<double>& pk_dm, std::vector<double>& pk_gas) {
  const int N = s.pc.Ngas;
  const double fb = omegab / omegam;
  const double fd = 1 - fb;

  link_dm(s);

  std::fill(s.density.begin(), s.density.end(), 0.0f);
  cic_mass_assignment_all(s);
#pragma omp parallel for schedule(static) if (par::worth(s.pc.Ncell))
  for (int k = 1; k <= N; ++k)
    for (int j = 1; j <= N; ++j)
      for (int i = 1; i <= N; ++i) {
        float& d = s.density[s.dcell(i, j, k)];
        d = static_cast<float>(static_cast<double>(d) / fd - 1);
      }
  fft::fft_3d('f', N, N, N, s.density.data());
  calc_pk_3d(N, s.density.data(), pk_dm, true, box);

  for (int k = 1; k <= N; ++k)
    for (int j = 1; j <= N; ++j) {
      for (int i = 1; i <= N; ++i)
        s.density[s.dcell(i, j, k)] =
            static_cast<float>(static_cast<double>(u_gas[s.gcell(i, j, k) * 5]) / fb - 1);
      for (int i = N + 1; i <= N + 2; ++i) s.density[s.dcell(i, j, k)] = 0;
    }
  fft::fft_3d('f', N, N, N, s.density.data());
  calc_pk_3d(N, s.density.data(), pk_gas, false, box);
}

void write_power_spectrum(const std::string& path, const PMState& s,
                          const std::vector<double>& pk_dm,
                          const std::vector<double>& pk_gas) {
  std::FILE* fp = std::fopen(path.c_str(), "w");
  if (!fp) throw std::runtime_error("cannot write " + path);
  // write(11,"(7es14.6)") pk_dm(:,k),pk_gas(2:3,k), k = 2..Hgas+1
  for (int k = 2; k <= s.pc.Hgas + 1; ++k) {
    const double* a = &pk_dm[(k - 1) * 5];
    const double* g = &pk_gas[(k - 1) * 5];
    std::fprintf(fp, "%14.6E%14.6E%14.6E%14.6E%14.6E%14.6E%14.6E\n", a[0], a[1], a[2], a[3],
                 a[4], g[1], g[2]);
  }
  std::fclose(fp);
}

// pm_tools.f90: summarize_dm
DMSummary summarize_dm(const PMState& s, const Units& u, std::ostream* out) {
  DMSummary r;
  for (int c = 0; c < 3; ++c) {
    r.xavg[c] = 0;
    r.xmax[c] = 0;
    r.xmin[c] = std::numeric_limits<double>::max();
    r.vavg[c] = 0;
    r.vmax[c] = 0;
    r.vmin[c] = std::numeric_limits<double>::max();
  }
  for (std::int64_t ip = s.dm_cpu[0]; ip <= s.dm_cpu[1]; ++ip) {
    const float* p = s.xv_p(ip);
    for (int c = 0; c < 3; ++c) {
      const double x = p[c], v = p[3 + c];
      r.xavg[c] = r.xavg[c] + x;
      r.xmax[c] = f::max(r.xmax[c], x);
      r.xmin[c] = f::min(r.xmin[c], x);
      r.vavg[c] = r.vavg[c] + v;
      r.vmax[c] = f::max(r.vmax[c], v);
      r.vmin[c] = f::min(r.vmin[c], v);
    }
  }
  const double npart = static_cast<double>(s.pc.Npart);
  const double c5 = static_cast<double>(1e5f);
  for (int c = 0; c < 3; ++c) {
    r.xavg[c] = r.xavg[c] / npart;
    r.vavg[c] = r.vavg[c] / npart;
    r.vavg[c] = r.vavg[c] * s.pc.x_dm_gas * u.vel_unit / c5;
    r.vmin[c] = r.vmin[c] * s.pc.x_dm_gas * u.vel_unit / c5;
    r.vmax[c] = r.vmax[c] * s.pc.x_dm_gas * u.vel_unit / c5;
  }
  if (out) {
    char buf[256];
    *out << " DM\n";
    const char* lab[6] = {" x  : ", " y  : ", " z  : ", " vx : ", " vy : ", " vz : "};
    for (int c = 0; c < 3; ++c) {
      std::snprintf(buf, sizeof buf, "%s %15.8g %15.8g %15.8g\n", lab[c],
                    static_cast<float>(r.xavg[c]), static_cast<float>(r.xmin[c]),
                    static_cast<float>(r.xmax[c]));
      *out << buf;
    }
    for (int c = 0; c < 3; ++c) {
      std::snprintf(buf, sizeof buf, "%s %15.8g %15.8g %15.8g\n", lab[3 + c],
                    static_cast<float>(r.vavg[c]), static_cast<float>(r.vmin[c]),
                    static_cast<float>(r.vmax[c]));
      *out << buf;
    }
  }
  return r;
}

namespace {
void read_f32(const std::string& path, float* dst, std::size_t n) {
  std::ifstream f(path, std::ios::binary);
  if (!f) throw std::runtime_error("cannot open " + path);
  f.read(reinterpret_cast<char*>(dst), static_cast<std::streamsize>(n * sizeof(float)));
  if (static_cast<std::size_t>(f.gcount()) != n * sizeof(float))
    throw std::runtime_error("short read on " + path);
}
}  // namespace

// io_tools.f90: read_init_conditions
void read_init_conditions(const std::string& dir, PMState& s, std::vector<float>& u_gas,
                          std::vector<float>& w_gas, const Units& u) {
  read_f32(dir + "/xv_dm.init", s.xv.data(), s.xv.size());
  read_f32(dir + "/u_gas.init", u_gas.data(), u_gas.size());
  const std::int64_t ncell = s.pc.Ncell;
  for (std::int64_t c = 0; c < ncell; ++c) {
    float* ug = &u_gas[c * 5];
    ug[4] = static_cast<float>(static_cast<double>(ug[4]) + static_cast<double>(ug[0]) * u.Tgas_min);
  }
  for (std::int64_t c = 0; c < ncell; ++c) {
    const double rho = u_gas[c * 5];
    w_gas[c * 4 + 0] = static_cast<float>(XHy * rho);
    w_gas[c * 4 + 1] = static_cast<float>(YHe * rho);
    w_gas[c * 4 + 2] = 0;
    w_gas[c * 4 + 3] = 0;
  }
}

// cosmo_tools.f90: calc_sim_var, redshift file tag
std::string ztag(double z) {
  char buf[32];
  if (z < 10)
    std::snprintf(buf, sizeof buf, "z=0%6.4f", z);
  else if (z < 100)
    std::snprintf(buf, sizeof buf, "z=%7.4f", z);
  else if (z < 1000)
    std::snprintf(buf, sizeof buf, "z=%7.3f", z);
  else
    buf[0] = '\0';  // the original keeps the previous tag
  return buf;
}

}  // namespace pm
}  // namespace radray
