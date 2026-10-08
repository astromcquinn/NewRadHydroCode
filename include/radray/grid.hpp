// Grid fields, structure-of-arrays.
//
// rt_grid (radhydro_global.f90: rt_grid_type, array rt_grid(Nrt,Nrt,Nrt)),
// u_gas(5,Ngas^3), w_gas(4,Ngas^3) and htoc_ray(4,Nrt^3).
//
// Indexing: cells are addressed with Fortran 1-based (i,j,k), i fastest:
//   cell(i,j,k) = (i-1) + (j-1)*N + (k-1)*N*N    (0-based linear index)
// Multi-component fields keep the component index fastest, exactly as in
// the Fortran derived type / arrays and the output format:
//   ngam(f, cell) -> ngam[cell*kNfreq + (f-1)]
#pragma once

#include <cstdint>
#include <vector>

#include "radray/constants.hpp"

namespace radray {

struct RTGrid {
  int N = 0;
  std::int64_t ncell = 0;

  // real(4) :: rhog,fHI,fHeI,fHeII,temp
  std::vector<float> rhog, fHI, fHeI, fHeII, temp;
  // real(4) :: Luv,sfr
  std::vector<float> Luv, sfr;
  // real(4), dimension(3) :: dx,gamma,lambda   (3 per cell)
  std::vector<float> dx, gamma, lambda;
  // real(4), dimension(Nfreq) :: nu,ngam,fgam,tau (kNfreq per cell)
  std::vector<float> nu, ngam, fgam, tau;

  // htoc_ray(4,Nrt,Nrt,Nrt): (1:2) = head/tail snapshot, (3:4) = live
  // head/tail of each cell's ray chain (ray storage indices, 0 = none).
  std::vector<std::int64_t> htoc;

  void allocate(int n);

  // 1-based Fortran cell indices -> linear index
  std::int64_t cell(int i, int j, int k) const {
    return (i - 1) + static_cast<std::int64_t>(N) * ((j - 1) + static_cast<std::int64_t>(N) * (k - 1));
  }
  // htoc_ray(m, cell), m = 1..4
  std::int64_t& htoc_ray(int m, std::int64_t c) { return htoc[c * 4 + (m - 1)]; }
  std::int64_t htoc_ray(int m, std::int64_t c) const { return htoc[c * 4 + (m - 1)]; }
  float* ngam_p(std::int64_t c) { return &ngam[c * kNfreq]; }
  float* nu_p(std::int64_t c) { return &nu[c * kNfreq]; }
  float* fgam_p(std::int64_t c) { return &fgam[c * kNfreq]; }
  float* tau_p(std::int64_t c) { return &tau[c * kNfreq]; }
  float* dx_p(std::int64_t c) { return &dx[c * 3]; }
  float* gamma_p(std::int64_t c) { return &gamma[c * 3]; }
  float* lambda_p(std::int64_t c) { return &lambda[c * 3]; }
};

struct GasGrid {
  int N = 0;
  std::int64_t ncell = 0;
  // u_gas(1:5,i,j,k): rho, rho*v (3), total energy density (code units)
  std::vector<float> u;
  // w_gas(1:4,i,j,k): rho_HI, rho_HeI, rho_HeII, rho_Z
  std::vector<float> w;

  void allocate(int n) {
    N = n;
    ncell = static_cast<std::int64_t>(n) * n * n;
    u.assign(5 * ncell, 0.0f);
    w.assign(4 * ncell, 0.0f);
  }
  std::int64_t cell(int i, int j, int k) const {
    return (i - 1) + static_cast<std::int64_t>(N) * ((j - 1) + static_cast<std::int64_t>(N) * (k - 1));
  }
  float* u_p(std::int64_t c) { return &u[c * 5]; }
  float* w_p(std::int64_t c) { return &w[c * 4]; }
};

}  // namespace radray
