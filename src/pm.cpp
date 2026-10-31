// Port of pm_tools.f90 (particle-mesh gravity for dark matter and gas) and
// the DM/gravity parts of radhydro.f90, hydro_tools.f90, cosmo_tools.f90 and
// io_tools.f90. Every expression keeps the Fortran operation order and
// precision; see docs/port_notes_pm.md.
#include "radray/pm.hpp"
#include "radray/kernels/hydro_mh.hpp"
#include "radray/kernels/pm_kernels.hpp"

#include <algorithm>
#include <limits>
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
#include "radray/pm_backend.hpp"

namespace radray {
namespace pm {

using namespace pc;
using namespace cp;

namespace {

// pm_tools.f90: cic_mass_assignment, local parameter
//   mass = (1-omegab/omegam)*N_dm_gas_c
// (N_dm_gas_c = (Nm/Ndm)^3: the PM mesh is the gas grid in the original)
inline double dm_mass(const PMConfig& pc) { return (1 - pc.cosmo.omegab / pc.cosmo.omegam) * pc.dm_mass_ratio; }

// CIC weights shared by cic_mass_assignment and accelerate_dm
// (rpk::cic_weights, kernels/pm_kernels.hpp: the same arithmetic as the
// Fortran, shared with the AMReX backend)
using rpk::Cic;
inline rpk::MeshGeom mesh_geom(const PMConfig& pc) { return rpk::MeshGeom{pc.Nm, pc.nmr, pc.x_dm_m}; }
inline Cic cic_weights(const PMConfig& pc, const float* xv) { return rpk::cic_weights(mesh_geom(pc), xv); }

// calc_kick: function force(i,j,k) (rpk::force_stencil: the whole finite
// difference is a real(4) expression; only the result is stored into the
// real(8) function value). (i,j,k) is a mesh cell.
inline void force(const PMState& s, int i, int j, int k, double F[3]) {
  const float* d = s.density.data();
  rpk::force_stencil(s.pc.Nm, i, j, k,
                     [&](int ii, int jj, int kk) -> float { return d[s.dcell(ii, jj, kk)]; }, F);
}

// force(i,j,k), from the stored mesh force when it is valid (pm_refine > 1,
// during the kick): the same float values
inline void mesh_force(const PMState& s, int i, int j, int k, double F[3]) {
  if (s.fmesh_valid) {
    const float* f = &s.fmesh[3 * s.mcell(i, j, k)];
    F[0] = f[0];
    F[1] = f[1];
    F[2] = f[2];
  } else {
    force(s, i, j, k, F);
  }
}

// The force on gas cell (i,j,k) in gas grid units. pm_refine = 1: the
// original force(i,j,k). Otherwise the average of the mesh forces over the
// r^3 mesh cells of the gas cell (k, j, i order, double sums), divided by r
// to convert mesh cells to gas cells: the adjoint of the spread in
// calc_density_field.
inline void gas_force(const PMState& s, int i, int j, int k, double F[3]) {
  const int r = s.pc.refine;
  if (r == 1) {
    force(s, i, j, k, F);
    return;
  }
  double sum[3] = {0, 0, 0};
  for (int kk = (k - 1) * r + 1; kk <= k * r; ++kk)
    for (int jj = (j - 1) * r + 1; jj <= j * r; ++jj)
      for (int ii = (i - 1) * r + 1; ii <= i * r; ++ii) {
        double Fm[3];
        mesh_force(s, ii, jj, kk, Fm);
        for (int m = 0; m < 3; ++m) sum[m] += Fm[m];
      }
  const double norm = static_cast<double>(r) * r * r * r;
  for (int m = 0; m < 3; ++m) F[m] = sum[m] / norm;
}

// CIC interpolation of the mesh force to a particle (accelerate_dm;
// rpk::particle_force)
inline void particle_force(const PMState& s, const Cic& c, double fv[3]) {
  rpk::particle_force(
      c, [&](int ii, int jj, int kk, double F[3]) { mesh_force(s, ii, jj, kk, F); }, fv);
}

// FFT threads for the pm_refine > 1 mesh. FFTW's OpenMP plans of 1024^3
// on a klone cpu-g2 node (2x96 cores): 1 thread 7.2 s per transform, 8
// threads 1.0 s, 32 threads 0.28 s, but 48-192 threads 3-10 s. The result
// is byte-identical for every count (test_fft_threads, 1..192), so the cap
// costs nothing. pm_refine = 1: the original serial plans.
inline int pm_fft_threads(const PMConfig& pc) {
  return pc.refine > 1 ? std::min(par::max_threads(), 32) : 1;
}

}  // namespace

PMConfig PMConfig::make(int Ngas, int Ndm, bool mf, int refine) {
  if (refine < 1) throw std::runtime_error("PMConfig: pm_refine must be >= 1");
  // the mesh a multiple of Ndm (the original: N_dm_gas particles' spacings
  // per cell), or Ndm a multiple of the mesh (more particles than mesh
  // cells: the zoom tests of pm_level_solve, docs/stage_g.md, whose refined
  // levels resolve them)
  const std::int64_t nm = static_cast<std::int64_t>(refine) * Ngas;
  if (Ngas <= 0 || Ndm <= 0 || (nm % Ndm != 0 && Ndm % nm != 0))
    throw std::runtime_error(refine == 1 ? "PMConfig: Ngas must be a positive multiple of Ndm, or Ndm of Ngas"
                                         : "PMConfig: the PM mesh pm_refine*Ngas must be a positive "
                                           "multiple of Ndm, or Ndm of the mesh");
  if (Ngas % 2 != 0) throw std::runtime_error("PMConfig: Ngas must be even");
  PMConfig c;
  c.Ngas = Ngas;
  c.Hgas = Ngas / 2;
  c.Ndm = Ndm;
  c.refine = refine;
  c.Nm = refine * Ngas;
  c.Hm = c.Nm / 2;
  c.Npart = static_cast<std::int64_t>(Ndm) * Ndm * Ndm;
  c.Ncell = static_cast<std::int64_t>(Ngas) * Ngas * Ngas;
  c.Nmcell = static_cast<std::int64_t>(c.Nm) * c.Nm * c.Nm;
  c.ngasr = Ngas;
  c.ndmr = Ndm;
  c.nmr = c.Nm;
  c.N_dm_m = c.Nm / Ndm;  // 0 with more particles than mesh cells
  c.N_dm_m_c = c.N_dm_m * c.N_dm_m * c.N_dm_m;
  c.x_dm_m = static_cast<double>(c.Nm) / Ndm;
  c.dm_mass_ratio = c.x_dm_m * c.x_dm_m * c.x_dm_m;
  c.x_dm_gas = static_cast<double>(Ngas) / Ndm;
  c.mf = mf;
  return c;
}

// pm_tools.f90: init_pm
void init_pm(PMState& s, const PMConfig& pc, const Domains& dom, bool host_arrays) {
  s.pc = pc;
  s.dom = dom;
  // pm_refine > 1: the same domains, in the same order, on the mesh
  if (pc.refine > 1)
    for (auto& b : s.dom.box)
      for (int e = 0; e < 6; e += 2) {
        b[e] = (b[e] - 1) * pc.refine + 1;
        b[e + 1] = b[e + 1] * pc.refine;
      }
  const std::int64_t N = pc.Nm;
  // npt = Npart/Ncpu + min(1,mod(Npart,Ncpu)) with Ncpu = 1
  s.dm_cpu[0] = 1;
  s.dm_cpu[1] = pc.Npart;
  if (!host_arrays) return;
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
  const int N = s.pc.Nm;
  const double xdg = s.pc.x_dm_m;
  std::fill(s.hoc_cpu.begin(), s.hoc_cpu.end(), 0);
  auto hc = [&](int j, int k) -> std::int64_t& {
    return s.hoc_cpu[(j - 1) + static_cast<std::int64_t>(N) * (k - 1)];
  };
  const int nchunk = std::min(par::max_threads(), 32);
  if (s.pc.refine > 1 && nchunk > 1 && par::worth(s.pc.Npart)) {
    // pm_refine > 1: the same row chains built in parallel. Chunk c of the
    // particles (consecutive ip) chains its own particles per row (head: its
    // largest ip, tail: its smallest); the chunks are then joined per row in
    // ip order. Every ll(ip) and hoc_cpu(j,k) is the serial one (the next
    // smaller ip of the same row), for any thread count.
    const std::int64_t rows = static_cast<std::int64_t>(N) * N;
    const std::int64_t p0 = s.dm_cpu[0], np = s.dm_cpu[1] - s.dm_cpu[0] + 1;
    std::vector<std::int64_t> head(nchunk * rows, 0), tail(nchunk * rows, 0);
#pragma omp parallel for schedule(static, 1) num_threads(nchunk)
    for (int c = 0; c < nchunk; ++c) {
      std::int64_t* H = &head[c * rows];
      std::int64_t* T = &tail[c * rows];
      const std::int64_t a = p0 + np * c / nchunk, b = p0 + np * (c + 1) / nchunk - 1;
      for (std::int64_t ip = a; ip <= b; ++ip) {
        const float* p = s.xv_p(ip);
        const int j = 1 + f::mod(f::int4(static_cast<double>(p[1]) * xdg), N);
        const int k = 1 + f::mod(f::int4(static_cast<double>(p[2]) * xdg), N);
        const std::int64_t row = (j - 1) + static_cast<std::int64_t>(N) * (k - 1);
        s.ll[ip] = H[row];
        if (H[row] == 0) T[row] = ip;
        H[row] = ip;
      }
    }
#pragma omp parallel for schedule(static)
    for (std::int64_t row = 0; row < rows; ++row) {
      std::int64_t h = 0;
      for (int c = 0; c < nchunk; ++c)
        if (head[c * rows + row] > 0) {
          s.ll[tail[c * rows + row]] = h;
          h = head[c * rows + row];
        }
      s.hoc_cpu[row] = h;
    }
  } else {
    for (std::int64_t ip = s.dm_cpu[0]; ip <= s.dm_cpu[1]; ++ip) {
      const float* p = s.xv_p(ip);
      const int j = 1 + f::mod(f::int4(static_cast<double>(p[1]) * xdg), N);
      const int k = 1 + f::mod(f::int4(static_cast<double>(p[2]) * xdg), N);
      s.ll[ip] = hc(j, k);
      hc(j, k) = ip;
    }
  }
  // the rows (j,k) are independent (each relinks only its own particles)
#pragma omp parallel for schedule(static) if (par::worth(s.pc.Nmcell))
  for (int k = 1; k <= N; ++k) {
    for (int j = 1; j <= N; ++j) {
      for (int i = 1; i <= N; ++i) s.hoc[s.mcell(i, j, k)] = 0;
      // loop over cpus (Ncpu..1): one
      std::int64_t ip = hc(j, k);
      while (ip > 0) {
        const int i = 1 + f::mod(f::int4(static_cast<double>(s.xv_p(ip)[0]) * xdg), N);
        const std::int64_t next = s.ll[ip];
        s.ll[ip] = s.hoc[s.mcell(i, j, k)];
        s.hoc[s.mcell(i, j, k)] = ip;
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
        std::int64_t ip = s.hoc[s.mcell(i, j, k)];
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
  const int N = s.pc.Nm;
  const int nt = par::worth(s.pc.Nmcell) ? std::min(par::max_threads(), N) : 1;
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
            std::int64_t ip = s.hoc[s.mcell(i, j, k)];
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
//   (pm_refine = r > 1: every mesh cell takes the density of the gas cell
//   it lies in, so each gas cell's mass is spread evenly over its r^3 mesh
//   cells; gas_force is the adjoint)
void calc_density_field(PMState& s, const float* u_gas) {
  const int N = s.pc.Nm;
  const int r = s.pc.refine;
#pragma omp parallel for schedule(static) if (par::worth(s.pc.Nmcell))
  for (int k = 1; k <= N; ++k)
    for (int j = 1; j <= N; ++j) {
      if (r == 1) {
        for (int i = 1; i <= N; ++i) s.density[s.dcell(i, j, k)] = u_gas[s.gcell(i, j, k) * 5];
      } else {
        const int jg = (j - 1) / r + 1, kg = (k - 1) / r + 1;
        for (int i = 1; i <= N; ++i)
          s.density[s.dcell(i, j, k)] = u_gas[s.gcell((i - 1) / r + 1, jg, kg) * 5];
      }
      for (int i = N + 1; i <= N + 2; ++i) s.density[s.dcell(i, j, k)] = 0;
    }
  // set_domain_gas/end_domain_gas hand out domains 1..Ndomain in order
  cic_mass_assignment_all(s);
}

// pm_tools.f90: calc_gravity_field
void calc_gravity_field(PMState& s, double scalefactor) {
  const int N = s.pc.Nm;
  const double Ak = 2 * pi / s.pc.nmr;
  // pm_refine > 1: threaded FFTW (fft.hpp); 1: the original serial plans
  const int fft_threads = pm_fft_threads(s.pc);
  s.fmesh_valid = false;
  // Aphi = -4*pi*scalefactor/(6*pi)  (unary minus applies to the quotient)
  const double Aphi = -(4 * pi * scalefactor / (6 * pi));

  fft::fft_3d('f', N, N, N, s.density.data(), fft_threads);

#pragma omp parallel for schedule(static) if (par::worth(s.pc.Nmcell))
  for (int k = 1; k <= N; ++k) {
    for (int j = 1; j <= N; ++j) {
      for (int i = 1; i <= N + 2; i += 2) {
        // rpk::green_wphi: kx = Ak*ix, ky/kz folded above the Nyquist mode
        const double wphi = rpk::green_wphi(N, (i - 1) / 2, j - 1, k - 1, Ak, Aphi);
        float* d = &s.density[s.dcell(i, j, k)];
        d[0] = static_cast<float>(static_cast<double>(d[0]) * wphi);
        d[1] = static_cast<float>(static_cast<double>(d[1]) * wphi);
      }
    }
  }

  fft::fft_3d('b', N, N, N, s.density.data(), fft_threads);
}

namespace {

// pm_refine > 1: force(i,j,k) of every mesh cell into s.fmesh (valid until
// the next calc_gravity_field)
void fill_force_mesh(PMState& s) {
  const int N = s.pc.Nm;
  s.fmesh_valid = false;
  s.fmesh.resize(3 * static_cast<std::size_t>(s.pc.Nmcell));
#pragma omp parallel for schedule(static) if (par::worth(s.pc.Nmcell))
  for (int k = 1; k <= N; ++k)
    for (int j = 1; j <= N; ++j)
      for (int i = 1; i <= N; ++i) {
        double F[3];
        force(s, i, j, k, F);
        float* f = &s.fmesh[3 * s.mcell(i, j, k)];
        for (int m = 0; m < 3; ++m) f[m] = static_cast<float>(F[m]);  // exact: F is a float
      }
  s.fmesh_valid = true;
}

// calc_kick, DM part: accelerate_dm(k) for all slabs. Particles are
// independent; the double sums arms and vrms are recorded per plane and
// added in plane order afterwards (same addends, same order as the serial
// loop), amax/vmax are maxima.
void accelerate_dm_all(PMState& s, double dt_sim, double& amax, double& arms, double& vmax,
                       double& vrms) {
  const int N = s.pc.Nm;
  const double xdg = s.pc.x_dm_m;
  const int nth = par::worth(s.pc.Nmcell) ? par::max_threads() : 1;
  const int batch = nth > 1 ? std::min(N, 4 * nth) : 1;
  std::vector<std::vector<double>> add_a(batch), add_v(batch);
  std::vector<double> amax_t(nth > 1 ? batch : 1, 0.0), vmax_t(nth > 1 ? batch : 1, 0.0);
  // pm_refine > 1: arms and vrms (diagnostics only) summed per plane, the
  // planes then in order: thread-count independent, not the serial order
  const bool plane_sums = s.pc.refine > 1;
  std::vector<double> sum_a(batch, 0.0), sum_v(batch, 0.0);
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
        std::int64_t ip = s.hoc[s.mcell(i, j, k)];
        while (ip > 0) {
          float* p = s.xv_p(ip);
          const Cic c = cic_weights(s.pc, p);
          double fv[3];
          particle_force(s, c, fv);

          // xv_dm(4:6,ip) = xv_dm(4:6,ip) + (f/x_dm_gas)*dt_sim  (x_dm_m)
          for (int m = 0; m < 3; ++m)
            p[3 + m] = static_cast<float>(static_cast<double>(p[3 + m]) + (fv[m] / xdg) * dt_sim);

          // arms = arms + sum(f**2); amax = max(amax,maxval(abs(f)))
          const double a2 = (0.0 + fv[0] * fv[0]) + fv[1] * fv[1] + fv[2] * fv[2];
          if (plane_sums)
            sum_a[b] += a2;
          else
            la.push_back(a2);
          amax = f::max(amax, f::max(f::max(std::fabs(fv[0]), std::fabs(fv[1])), std::fabs(fv[2])));
          // vrms = vrms + sum(xv_dm(4:6,ip)**2)   (real(4) sum)
          const float vs = ((0.0f + p[3] * p[3]) + p[4] * p[4]) + p[5] * p[5];
          if (plane_sums)
            sum_v[b] += static_cast<double>(vs);
          else
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
    if (plane_sums) {
      arms = arms + sum_a[b];
      vrms = vrms + sum_v[b];
      sum_a[b] = sum_v[b] = 0;
      continue;
    }
    for (double a : add_a[b]) arms = arms + a;
    for (double v : add_v[b]) vrms = vrms + v;
    add_a[b].clear();
    add_v[b].clear();
  }
  }
  for (double a : amax_t) amax = f::max(amax, a);
  for (double v : vmax_t) vmax = f::max(vmax, v);
}

}  // namespace

void force_at(const PMState& s, const float* xv, double F[3]) {
  particle_force(s, cic_weights(s.pc, xv), F);
}

namespace {

// Update global variables in PM grid units (end of calc_kick)
void set_kick_globals(const PMState& s, SimClock& clk, double amax, double arms, double vmax,
                      double vrms) {
  const double npart = static_cast<double>(s.pc.Npart);
  clk.amax_dm = amax;
  clk.arms_dm = std::sqrt(arms / npart);
  clk.vmax_dm = vmax;
  clk.vrms_dm = std::sqrt(vrms / npart);
}

}  // namespace

// pm_tools.f90: calc_kick (accelerate_dm for all slabs, then accelerate_gas)
void calc_kick(PMState& s, float* u_gas, float* v_gas, SimClock& clk, const Units& u) {
  const double dt_sim = clk.dt_sim;
  const int N = s.pc.Ngas;
  double amax = 0, arms = 0, vmax = 0, vrms = 0;
  if (s.pc.refine > 1) fill_force_mesh(s);
  accelerate_dm_all(s, dt_sim, amax, arms, vmax, vrms);
  const int nth = par::worth(s.pc.Ncell) ? par::max_threads() : 1;

  // accelerate_gas(k): cells are independent
#pragma omp parallel for schedule(static) if (nth > 1)
  for (int k = 1; k <= N; ++k)
    for (int j = 1; j <= N; ++j)
      for (int i = 1; i <= N; ++i) {
        double F[3];
        gas_force(s, i, j, k, F);
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

  set_kick_globals(s, clk, amax, arms, vmax, vrms);
}

// pm_tools.f90: particle_mesh (pm_backend = amrex: link, density, gravity
// and kick on the AMReX level, docs/stage_e_pm.md)
void particle_mesh(PMState& s, float* u_gas, float* v_gas, SimClock& clk, const Units& u) {
  calc_drift(s, clk.dt_sim);
  if (s.pc.backend == PMBackend::AMReX) {
    amrex_pm_step(s, u_gas, v_gas, clk, u);
  } else {
    link_dm(s);
    calc_density_field(s, u_gas);
    calc_gravity_field(s, clk.scalefactor);
    calc_kick(s, u_gas, v_gas, clk, u);
  }
  calc_drift(s, clk.dt_sim);
}

// hydro_scheme = unsplit (docs/hydro_unsplit.md section 3): the DM part of
// calc_kick unchanged; the gas kick acts on the conserved variables of the
// lab frame (v_gas = 0): rho v += rho g dt and E = KE_new + e, where e is
// the thermal energy selected by the dual-energy rule (E - KE when it is a
// reliable fraction of E, else the internal-energy field, which gravity
// does not change). Computed in double from the float state.
void calc_kick_unsplit(PMState& s, float* u_gas, const float* e_gas, double eta1, SimClock& clk,
                       const Units& u) {
  const double dt_sim = clk.dt_sim;
  const int N = s.pc.Ngas;
  double amax = 0, arms = 0, vmax = 0, vrms = 0;
  if (s.pc.refine > 1) fill_force_mesh(s);
  accelerate_dm_all(s, dt_sim, amax, arms, vmax, vrms);
  const int nth = par::worth(s.pc.Ncell) ? par::max_threads() : 1;
#pragma omp parallel for schedule(static) if (nth > 1)
  for (int k = 1; k <= N; ++k)
    for (int j = 1; j <= N; ++j)
      for (int i = 1; i <= N; ++i) {
        double F[3];
        gas_force(s, i, j, k, F);
        const std::int64_t c = s.gcell(i, j, k);
        float* ug = &u_gas[c * 5];
        double U[5];
        for (int m = 0; m < 5; ++m) U[m] = ug[m];
        const double eD = e_gas ? static_cast<double>(e_gas[c]) : 0.0;
        const double eth = rrh::thermal_energy(U, eD, eta1);
        const double rho = f::max(U[0], u.Dgas_min);
        for (int m = 0; m < 3; ++m) U[1 + m] = U[1 + m] + rho * F[m] * dt_sim;
        U[4] = rrh::kinetic_energy(U) + eth;
        for (int m = 1; m < 5; ++m) ug[m] = static_cast<float>(U[m]);
      }
  set_kick_globals(s, clk, amax, arms, vmax, vrms);
}

void particle_mesh_unsplit(PMState& s, float* u_gas, const float* e_gas, double eta1,
                           SimClock& clk, const Units& u) {
  if (s.pc.amrex_state) {
    // amrex_state = 1: the particles and the gas live on the level
    // (pm_backend.hpp); the drifts there too
    amrex_pm_drift(s, clk.dt_sim);
    amrex_pm_step_unsplit(s, u_gas, e_gas, eta1, clk, u);
    amrex_pm_drift(s, clk.dt_sim);
    return;
  }
  calc_drift(s, clk.dt_sim);
  if (s.pc.backend == PMBackend::AMReX) {
    amrex_pm_step_unsplit(s, u_gas, e_gas, eta1, clk, u);
  } else {
    link_dm(s);
    calc_density_field(s, u_gas);
    calc_gravity_field(s, clk.scalefactor);
    calc_kick_unsplit(s, u_gas, e_gas, eta1, clk, u);
  }
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
            wcic = f::pow2(pk_sinc(kx / 2) * pk_sinc(ky / 2) * pk_sinc(kz / 2));
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
  calc_pk_finish(nc, pk, Dk, box);
}


// pm_tools.f90: calc_power_spectrum (computation; file name/IO in the caller)
// pm_refine > 1: the DM power on the mesh (rows k = 1..Nm of pk_dm), the gas
// power on the gas grid in a separate buffer; density then holds the DM
// transform.
void calc_power_spectrum(PMState& s, const float* u_gas, double box,
                         std::vector<double>& pk_dm, std::vector<double>& pk_gas, int stride) {
  const int N = s.pc.Ngas;
  const int Nm = s.pc.Nm;
  const double fb = s.pc.cosmo.omegab / s.pc.cosmo.omegam;
  const double fd = 1 - fb;
  const int fft_threads = pm_fft_threads(s.pc);
  s.fmesh_valid = false;

  link_dm(s);

  std::fill(s.density.begin(), s.density.end(), 0.0f);
  cic_mass_assignment_all(s);
#pragma omp parallel for schedule(static) if (par::worth(s.pc.Nmcell))
  for (int k = 1; k <= Nm; ++k)
    for (int j = 1; j <= Nm; ++j)
      for (int i = 1; i <= Nm; ++i) {
        float& d = s.density[s.dcell(i, j, k)];
        d = static_cast<float>(static_cast<double>(d) / fd - 1);
      }
  fft::fft_3d('f', Nm, Nm, Nm, s.density.data(), fft_threads);
  calc_pk_3d(Nm, s.density.data(), pk_dm, true, box);

  std::vector<float> gbuf;
  float* g = s.density.data();  // density(Ngas+2,Ngas,Ngas) when the mesh is the gas grid
  if (s.pc.refine > 1) {
    gbuf.assign(static_cast<std::size_t>(N + 2) * N * N, 0.0f);
    g = gbuf.data();
  }
  const std::int64_t ld = N + 2;
  for (int k = 1; k <= N; ++k)
    for (int j = 1; j <= N; ++j) {
      float* row = g + ld * ((j - 1) + static_cast<std::int64_t>(N) * (k - 1));
      for (int i = 1; i <= N; ++i)
        row[i - 1] = static_cast<float>(static_cast<double>(u_gas[s.gcell(i, j, k) * stride]) / fb - 1);
      for (int i = N + 1; i <= N + 2; ++i) row[i - 1] = 0;
    }
  fft::fft_3d('f', N, N, N, g, fft_threads);
  calc_pk_3d(N, g, pk_gas, false, box);
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
namespace {
void dm_summary_init(DMSummary& r) {
  for (int c = 0; c < 3; ++c) {
    r.xavg[c] = 0;
    r.xmax[c] = 0;
    r.xmin[c] = std::numeric_limits<double>::max();
    r.vavg[c] = 0;
    r.vmax[c] = 0;
    r.vmin[c] = std::numeric_limits<double>::max();
  }
}
void dm_summary_add(const PMState& s, std::int64_t a, std::int64_t b, DMSummary& r) {
  for (std::int64_t ip = a; ip <= b; ++ip) {
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
}
}  // namespace

DMSummary summarize_dm(const PMState& s, const Units& u, std::ostream* out) {
  DMSummary r;
  dm_summary_init(r);
  if (s.pc.refine > 1 && par::worth(s.pc.Npart)) {
    // pm_refine > 1: sums over fixed chunks of 2^20 particles, added in
    // chunk order (independent of the thread count; not the serial order)
    const std::int64_t ch = std::int64_t(1) << 20;
    const std::int64_t np = s.dm_cpu[1] - s.dm_cpu[0] + 1, nch = (np + ch - 1) / ch;
    std::vector<DMSummary> part(nch);
#pragma omp parallel for schedule(dynamic, 16)
    for (std::int64_t q = 0; q < nch; ++q) {
      dm_summary_init(part[q]);
      dm_summary_add(s, s.dm_cpu[0] + q * ch, std::min(s.dm_cpu[1], s.dm_cpu[0] + (q + 1) * ch - 1),
                     part[q]);
    }
    for (const DMSummary& p : part)
      for (int c = 0; c < 3; ++c) {
        r.xavg[c] = r.xavg[c] + p.xavg[c];
        r.xmax[c] = f::max(r.xmax[c], p.xmax[c]);
        r.xmin[c] = f::min(r.xmin[c], p.xmin[c]);
        r.vavg[c] = r.vavg[c] + p.vavg[c];
        r.vmax[c] = f::max(r.vmax[c], p.vmax[c]);
        r.vmin[c] = f::min(r.vmin[c], p.vmin[c]);
      }
  } else {
    dm_summary_add(s, s.dm_cpu[0], s.dm_cpu[1], r);
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
                          std::vector<float>& w_gas, const Units& u, double T_init_K) {
  read_f32(dir + "/xv_dm.init", s.xv.data(), s.xv.size());
  if (u_gas.empty()) return;  // amrex_state = 1 on a rank other than 0: no host gas grid
  read_f32(dir + "/u_gas.init", u_gas.data(), u_gas.size());
  const std::int64_t ncell = s.pc.Ncell;
  const double Tinit = ics_tinit(u, T_init_K);
  for (std::int64_t c = 0; c < ncell; ++c) ics_gas_cell(&u_gas[c * 5], &w_gas[c * 4], Tinit);
}

// thermal_init = recfast: the RECFAST temperature (neutral mu, temp1_unit)
double ics_tinit(const Units& u, double T_init_K) { return T_init_K > 0 ? T_init_K / u.temp1_unit : u.Tgas_min; }

// u(5) += rho Tinit; w neutral (each cell its own: the same values as the
// original's two loops over the grid)
void ics_gas_cell(float* ug, float* wg, double Tinit) {
  ug[4] = static_cast<float>(static_cast<double>(ug[4]) + static_cast<double>(ug[0]) * Tinit);
  const double rho = ug[0];
  wg[0] = static_cast<float>(XHy * rho);
  wg[1] = static_cast<float>(YHe * rho);
  wg[2] = 0;
  wg[3] = 0;
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
