// plane_transport = sweep (plane_parallel.hpp, docs/plane_parallel.md): the
// serial RT's photon packets of the plane-parallel RT-domain mode as a
// field per lit face instead of rays. A packet "about to enter cell e
// along the face's axis" has travelled m(e) cells from its source face
// (m = e mod L for a + face, L - 1 - e mod L for a - face); packets exist
// only at m = 0 mod k (k = rt_cells_per_substep, which divides L), the
// packet at m = 0 is the source packet (never stored). Every substep the
// cast (rt_step_1) deposits each packet along its k-cell path with
// rrk::cast_1, as the ray it stands for would, and the advance
// (rt_step_3) attenuates it with rrk::cast_3 and moves it k cells on, or
// counts it as escaped at m + k = L, or drops it when faint. The same
// kernels and arithmetic as the rays, so the deposits agree with the ray
// transport to round-off (one packet per cell and direction: the float
// sums are in direction order either way).
#include <cmath>
#include <vector>

#include "radray/kernels/ray_kernels.hpp"
#include "radray/parallel.hpp"
#include "radray/plane_parallel.hpp"
#include "radray/ray_tracing.hpp"
#include "radray/simulation.hpp"

namespace radray {

namespace {

using rrk::plane_face_sign;
using rrk::plane_packet_m;
using rrk::plane_packet_ray;

inline bool all_zero(const float* v) {
  for (int q = 0; q < rrk::NF; ++q)
    if (v[q] != 0) return false;
  return true;
}

}  // namespace

void plane_sweep_cast(Simulation& s) {
  PlaneState& ps = s.plane;
  if (!ps.on || ps.transport != 1) return;
  RTGrid& G = s.rt;
  const rrk::RTP P = rt_params(s);
  const int N = s.cfg.Nrt, L = ps.L, k = P.cells;
  const std::int64_t ncol = static_cast<std::int64_t>(N) * N;
  std::int64_t ncst = 0;
  for (std::size_t fi = 0; fi < ps.face_list.size(); ++fi) {
    const int f = ps.face_list[fi], d = f / 2, sgn = plane_face_sign(f);
    const int e1 = (d + 1) % 3, e2 = (d + 2) % 3;
    const float* F = ps.field[fi].data();
    // one column of cells along the axis per thread: a packet deposits
    // only into its own column, so the cell sums are in a fixed order
#pragma omp parallel for schedule(static) reduction(+ : ncst) if (par::worth(ncol))
    for (std::int64_t q = 0; q < ncol; ++q) {
      int idx[3];
      idx[e1] = static_cast<int>(q % N);
      idx[e2] = static_cast<int>(q / N);
      auto tau = [&](int i, int j, int kk) -> const float* { return G.tau_p(G.cell(i + 1, j + 1, kk + 1)); };
      auto dep = [&](int i, int j, int kk, const double* v) {
        const std::int64_t c = G.cell(i + 1, j + 1, kk + 1);
        float* gng = G.ngam_p(c);
        float* gnu = G.nu_p(c);
        float* gdx = G.dx_p(c);
        for (int m = 0; m < kNfreq; ++m) gng[m] = static_cast<float>(static_cast<double>(gng[m]) + v[m]);
        for (int m = 0; m < kNfreq; ++m) gnu[m] = static_cast<float>(static_cast<double>(gnu[m]) + v[5 + m]);
        for (int m = 0; m < 3; ++m) gdx[m] = static_cast<float>(static_cast<double>(gdx[m]) + v[10 + m]);
      };
      for (int e = 0; e < N; ++e) {
        const int m = plane_packet_m(e, L, sgn);
        if (m % k != 0) continue;
        idx[d] = e;
        const std::int64_t c = G.cell(idx[0] + 1, idx[1] + 1, idx[2] + 1);
        const float* pk = m == 0 ? P.plane.ngam : F + static_cast<std::size_t>(kNfreq) * c;
        if (m > 0 && all_zero(pk)) continue;
        const rrk::RayP ry = plane_packet_ray(idx, f, m, P, pk);
        if (rrk::classify(ry, P) == rrk::kRemove) continue;  // fainter than ngam_min
        rrk::cast_1(ry, P, tau, dep);
        ++ncst;
      }
    }
  }
  s.rays.stat.ncst1 += ncst;
  s.rays.stat.ncst2 = ncst > s.rays.stat.ncst2 ? ncst : s.rays.stat.ncst2;
}

void plane_sweep_advance(Simulation& s) {
  PlaneState& ps = s.plane;
  if (!ps.on || ps.transport != 1) return;
  RTGrid& G = s.rt;
  const rrk::RTP P = rt_params(s);
  const int N = s.cfg.Nrt, L = ps.L, k = P.cells;
  const std::int64_t ncol = static_cast<std::int64_t>(N) * N;
  const int nt = par::max_threads();
  std::vector<__int128> esc_t(nt, 0), inf_t(nt, 0);
  std::int64_t ndel = 0;
  for (std::size_t fi = 0; fi < ps.face_list.size(); ++fi) {
    const int f = ps.face_list[fi], d = f / 2, sgn = plane_face_sign(f);
    const int e1 = (d + 1) % 3, e2 = (d + 2) % 3;
    float* F = ps.field[fi].data();
#pragma omp parallel for schedule(static) reduction(+ : ndel) if (par::worth(ncol))
    for (std::int64_t q = 0; q < ncol; ++q) {
      const int t = par::thread_id();
      int idx[3];
      idx[e1] = static_cast<int>(q % N);
      idx[e2] = static_cast<int>(q / N);
      auto tf = [&](int i, int j, int kk, const float*& tau, const float*& fgam) {
        const std::int64_t c = G.cell(i + 1, j + 1, kk + 1);
        tau = G.tau_p(c);
        fgam = G.fgam_p(c);
      };
      // in place along the column, against the direction of motion: the
      // destination of a packet (k cells on) was consumed before it is
      // written; no packet crosses a domain boundary (it escapes there)
      for (int n = 0; n < N; ++n) {
        const int e = sgn > 0 ? N - 1 - n : n;
        const int m = plane_packet_m(e, L, sgn);
        if (m % k != 0) continue;
        idx[d] = e;
        const std::int64_t c = G.cell(idx[0] + 1, idx[1] + 1, idx[2] + 1);
        float* slot = F + static_cast<std::size_t>(kNfreq) * c;
        float src[rrk::NF];
        const float* pk;
        if (m == 0) {
          pk = P.plane.ngam;
        } else {
          for (int qq = 0; qq < rrk::NF; ++qq) src[qq] = slot[qq];
          pk = src;
          for (int qq = 0; qq < rrk::NF; ++qq) slot[qq] = 0;  // consumed
          if (all_zero(src)) continue;
        }
        rrk::RayP ry = plane_packet_ray(idx, f, m, P, pk);
        if (rrk::classify(ry, P) == rrk::kRemove || !rrk::cast_3(ry, P, tf)) {
          ++ndel;  // faint
          continue;
        }
        if (m + k == L) {
          esc_t[t] += plane_to_fx(rrk::sum_ngam(ry));  // leaves at the far face
          ++ndel;
          continue;
        }
        idx[d] = e + sgn * k;  // within the domain (m + k < L)
        float* dst = F + static_cast<std::size_t>(kNfreq) * G.cell(idx[0] + 1, idx[1] + 1, idx[2] + 1);
        for (int qq = 0; qq < rrk::NF; ++qq) dst[qq] = ry.ngam[qq];
        inf_t[t] += plane_to_fx(rrk::sum_ngam(ry));
      }
    }
  }
  __int128 esc = 0, inf = 0;
  for (int t = 0; t < nt; ++t) {
    esc += esc_t[t];
    inf += inf_t[t];
  }
  ps.escaped_fx += esc;
  ps.inflight_fx = inf;
  s.rays.stat.ndel1 += ndel;
}

}  // namespace radray
