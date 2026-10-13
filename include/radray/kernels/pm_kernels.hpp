// The per-particle and per-cell arithmetic of the PM gravity
// (pm_tools.f90: cic_mass_assignment, calc_kick's force function and
// accelerate_dm), in one place for the serial port (src/pm.cpp) and the
// AMReX backend (src/amrex/pm_level.cpp, pm_backend = amrex,
// docs/stage_e_pm.md). Every expression keeps the Fortran operation order
// and precision (docs/port_notes_pm.md); the serial port calling these is
// byte-identical to the gfortran reference.
//
// Units: the PM mesh of Nm cells per dimension (the gas grid when
// pm_refine = 1); particle positions xv[0..2] in DM-grid units
// (x_dm_m = Nm/Ndm mesh cells per DM cell), velocities xv[3..5].
// Mesh cells are 1-based (i, j, k), as in the Fortran.
#pragma once

#include <cstdint>

#include "radray/kernels/rr_hd.hpp"

namespace radray {
namespace rpk {

// Geometry of the CIC on the mesh: the fields of pm::PMConfig that the
// kernels read.
struct MeshGeom {
  int Nm = 0;        // mesh cells per dimension
  double nmr = 0;    // Nm as a real
  double x_dm_m = 1; // mesh cells per DM-grid cell (Nm/Ndm)
};

// CIC weights shared by cic_mass_assignment and accelerate_dm:
//   x = mod(xv_dm(1,ip)*x_dm_gas - 0.5 + ngasr, ngasr), i1 = 1+int(x),
//   i2 = 1+mod(i1,N), dx1 = i1-x, dx2 = 1-dx1   (and y, z)
// on the PM mesh (x_dm_m, nmr; x_dm_gas, ngasr with pm_refine = 1).
struct Cic {
  int i1, i2, j1, j2, k1, k2;
  double dx1, dx2, dy1, dy2, dz1, dz2;
};
RR_HD Cic cic_weights(const MeshGeom& g, const float* xv) {
  const int N = g.Nm;
  Cic c;
  const double x = rrf::mod(static_cast<double>(xv[0]) * g.x_dm_m - 0.5 + g.nmr, g.nmr);
  const double y = rrf::mod(static_cast<double>(xv[1]) * g.x_dm_m - 0.5 + g.nmr, g.nmr);
  const double z = rrf::mod(static_cast<double>(xv[2]) * g.x_dm_m - 0.5 + g.nmr, g.nmr);
  c.i1 = 1 + rrf::int4(x);
  c.i2 = 1 + rrf::mod(c.i1, N);
  c.dx1 = c.i1 - x;
  c.dx2 = 1 - c.dx1;
  c.j1 = 1 + rrf::int4(y);
  c.j2 = 1 + rrf::mod(c.j1, N);
  c.dy1 = c.j1 - y;
  c.dy2 = 1 - c.dy1;
  c.k1 = 1 + rrf::int4(z);
  c.k2 = 1 + rrf::mod(c.k1, N);
  c.dz1 = c.k1 - z;
  c.dz2 = 1 - c.dz1;
  return c;
}

// The eight CIC corners (i, j, k, weight) of the particle, in the order of
// cic_mass_assignment / accelerate_dm. f(i, j, k, w) is called for each.
template <class F>
RR_HD void cic_corners(const Cic& c, F&& f) {
  f(c.i1, c.j1, c.k1, c.dx1 * c.dy1 * c.dz1);
  f(c.i2, c.j1, c.k1, c.dx2 * c.dy1 * c.dz1);
  f(c.i1, c.j2, c.k1, c.dx1 * c.dy2 * c.dz1);
  f(c.i2, c.j2, c.k1, c.dx2 * c.dy2 * c.dz1);
  f(c.i1, c.j1, c.k2, c.dx1 * c.dy1 * c.dz2);
  f(c.i2, c.j1, c.k2, c.dx2 * c.dy1 * c.dz2);
  f(c.i1, c.j2, c.k2, c.dx1 * c.dy2 * c.dz2);
  f(c.i2, c.j2, c.k2, c.dx2 * c.dy2 * c.dz2);
}

// calc_kick: function force(i,j,k), the 4-point gradient of the potential
// on the periodic mesh. The whole finite difference is a real(4)
// expression (2./3 and 1./12 are real(4) literals, the potential is
// real(4)); only the result is stored into the real(8) function value.
// phi(i, j, k) returns the potential of a mesh cell (1-based, wrapped
// indices are passed; the caller may map them to its own storage).
template <class Phi>
RR_HD void force_stencil(int N, int i, int j, int k, Phi&& phi, double F[3]) {
  const int im2 = 1 + rrf::mod(i - 3 + N, N);
  const int im1 = 1 + rrf::mod(i - 2 + N, N);
  const int ip1 = 1 + rrf::mod(i, N);
  const int ip2 = 1 + rrf::mod(i + 1 + N, N);
  const int jm2 = 1 + rrf::mod(j - 3 + N, N);
  const int jm1 = 1 + rrf::mod(j - 2 + N, N);
  const int jp1 = 1 + rrf::mod(j, N);
  const int jp2 = 1 + rrf::mod(j + 1 + N, N);
  const int km2 = 1 + rrf::mod(k - 3 + N, N);
  const int km1 = 1 + rrf::mod(k - 2 + N, N);
  const int kp1 = 1 + rrf::mod(k, N);
  const int kp2 = 1 + rrf::mod(k + 1 + N, N);
  const float c23 = 2.0f / 3.0f;
  const float c112 = 1.0f / 12.0f;
  const float fx = -(c23 * (phi(ip1, j, k) - phi(im1, j, k))) + c112 * (phi(ip2, j, k) - phi(im2, j, k));
  const float fy = -(c23 * (phi(i, jp1, k) - phi(i, jm1, k))) + c112 * (phi(i, jp2, k) - phi(i, jm2, k));
  const float fz = -(c23 * (phi(i, j, kp1) - phi(i, j, km1))) + c112 * (phi(i, j, kp2) - phi(i, j, km2));
  F[0] = fx;
  F[1] = fy;
  F[2] = fz;
}

// CIC interpolation of the mesh force to a particle (accelerate_dm):
// fv = sum over the corners of force(corner) * wx*wy*wz, in the order of
// the Fortran. force(i, j, k, F) fills the three components of a cell.
template <class Force>
RR_HD void particle_force(const Cic& c, Force&& force, double fv[3]) {
  double F[3];
  fv[0] = fv[1] = fv[2] = 0;
  auto acc = [&](int ii, int jj, int kk, double wx, double wy, double wz) {
    force(ii, jj, kk, F);
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
}

// calc_gravity_field: the Green's function factor of the spectral mode
// with 0-based indices (ix = 0..Nm/2 of the half-complex x axis, jy and kz
// = 0..Nm-1), Ak = 2 pi / nmr, Aphi = -(4 pi a / (6 pi)). Zero for the
// k = 0 mode. Multiplies real and imaginary parts as the Fortran does
// (double product, stored as real(4)).
RR_HD double green_wphi(int Nm, int ix, int jy, int kz, double Ak, double Aphi) {
  const int H = Nm / 2;
  const double kx = Ak * ix;
  const double ky = (jy <= H) ? Ak * jy : Ak * (jy - Nm);
  const double kzz = (kz <= H) ? Ak * kz : Ak * (kz - Nm);
  const double ka = 2 * std::sin(kx / 2);
  const double kb = 2 * std::sin(ky / 2);
  const double kc = 2 * std::sin(kzz / 2);
  const double ksq = rrf::pow2(kx) + rrf::pow2(ky) + rrf::pow2(kzz);
  if (ksq > 0) return Aphi / (rrf::pow2(ka) + rrf::pow2(kb) + rrf::pow2(kc));
  return 0;
}

}  // namespace rpk
}  // namespace radray
