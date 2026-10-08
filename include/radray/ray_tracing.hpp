// Port of ray_tracing_tools.f90 (RAYFAST undefined, FREQ1 undefined).
//
// Conventions: ray storage indices `ir` are 1-based (0 = none); RT cell
// indices (irt,jrt,krt) are 1-based Fortran indices; `c` denotes the
// 0-based linear cell index RTGrid::cell(irt,jrt,krt).
#pragma once

#include <array>
#include <cstdint>

#include "radray/simulation.hpp"

namespace radray {

// ray_tracing_tools.f90: init_ray_tracing
void init_ray_tracing(Simulation& s);

// ray_tracing_tools.f90: ray_tracing_1 / ray_tracing_3 over one domain box
// {i1,i2,j1,j2,k1,k2} (1-based inclusive).
void ray_tracing_1(Simulation& s, const std::array<int, 6>& box);
void ray_tracing_3(Simulation& s, const std::array<int, 6>& box);

// ray_tracing_tools.f90: ray_source
void ray_source(Simulation& s, int irt, int jrt, int krt);
// ray_tracing_tools.f90: ray_cast_split
void ray_cast_split(Simulation& s, std::int64_t ir, int irt, int jrt, int krt);
// ray_tracing_tools.f90: ray_cast_merge (irc: nrc ray indices of the cone)
void ray_cast_merge(Simulation& s, int nrc, const std::int64_t* irc, int irt,
                    int jrt, int krt);
// ray_tracing_tools.f90: ray_cast_1 / ray_cast_3 (non-RAYFAST versions).
// ray_cast_3 sets ir = 0 if the ray was deleted.
void ray_cast_1(Simulation& s, std::int64_t ir, int irt, int jrt, int krt);
void ray_cast_3(Simulation& s, std::int64_t& ir, int irt, int jrt, int krt);

// ray_tracing_tools.f90: ray_segment -- lengths of the cell crossings of a
// unit-length step from x0 along dx. The original computes nseg = 3
// segments, which can truncate the step (B8); nseg = 4 (legacy_rt = 0)
// always covers it. Elements beyond nseg are 0.
std::array<double, 4> ray_segment(const double x0[3], const double dx[3],
                                  double seg_min, int nseg = 3);

// ray_tracing_tools.f90: ray_create / ray_remove / ray_insert_head/tail
void ray_create(Simulation& s, int nr, std::int64_t* indx);
void ray_remove(Simulation& s, std::int64_t ir, int irt, int jrt, int krt,
                bool flag);
void ray_insert_head(Simulation& s, std::int64_t ir, int irt, int jrt, int krt);
void ray_insert_tail(Simulation& s, std::int64_t ir, int irt, int jrt, int krt);

// ray_tracing_tools.f90: lside/nside/nray/rmin/rmax helpers
int lside_from_radius(double r);
int nside_from_radius(double r);
int nray_from_nside(int nside);
double rmin_from_nside(int nside);
double rmax_from_nside(int nside);

// ray_tracing_tools.f90: ray_redshifting (non-FREQ1 redshift_ray)
void ray_redshifting(Simulation& s);

// ray_tracing_tools.f90: ray_buffer_allocate / _supplement / _deallocate
void ray_buffer_allocate(Simulation& s);
void ray_buffer_supplement(Simulation& s);
void ray_buffer_deallocate(Simulation& s);

// general_tools.f90: interpolate(x,f,a,b,'lin') for f(2,n) stored as
// pairs f[2*i+0] = f(1,i+1), f[2*i+1] = f(2,i+1); a,b in {1,2}.
double interpolate_lin(double x, const double* f, int n, int a, int b);

}  // namespace radray
