// Port of hydro_tools.f90 (Hy Trac's RadHydro hydro solver) with -DMF,
// plus the hydro parts of radhydro.f90 (calc_timestep, radhydro).
//
// Faithfulness rules (docs/port_notes.md section 3, docs/port_notes_hydro.md):
//   * real(4) -> float, real(8) -> double; every expression keeps the
//     Fortran operation order and the Fortran promotion points. In
//     particular:
//       - sum(u(2:4,i)**2) of a real(4) array is a float sum of float squares;
//       - (v(:,ip)-v(:,im))/2 and (u(:,im)-u(:,i))/4 of real(4) arrays are
//         float operations, widened afterwards;
//       - u(:,i)*v(1,i) (advect, real(4) * real(4)) is a float product;
//       - assignments to real(4) arrays round to float at that point.
//   * MAX/MIN follow gfortran exactly: NaN propagates, and on ties MAX
//     returns its first and MIN its second argument (gmax/gmin below;
//     checked against gfortran 13 -O2, see port_notes_hydro.md).
//   * Single-precision literals (0.1, 2E8, 1E5, 2.740, ...) are float.
//   * The 1-D FFTs (MKL DFTI in the original) use FFTW3 single precision
//     with FFTW_ESTIMATE plans on an aligned n+2 buffer, exactly like the
//     MKL_DFTI shim of reference/harness_hydro, so both agree bit for bit.
#include "radray/hydro.hpp"
#include "radray/hydro_unsplit.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <memory>
#include <string>
#include <stdexcept>

#include "radray/constants.hpp"
#include "radray/fft.hpp"
#include "radray/fortran.hpp"
#include "radray/parallel.hpp"

namespace radray {

using namespace pc;
using namespace cp;

namespace {

constexpr double g = gamma_ideal;

// gfortran MAX/MIN semantics (two arguments, same kind after promotion)
inline double gmax(double a, double b) { return f::gmax(a, b); }
inline double gmin(double a, double b) { return f::gmin(a, b); }
inline float gmaxf(float a, float b) { return f::gmax(a, b); }

// sum(x(1:3)**2): real(4) array -> float accumulation; real(8) -> double
inline float sumsq3(const float* x) {
  float s = 0.0f;
  s = s + x[0] * x[0];
  s = s + x[1] * x[1];
  s = s + x[2] * x[2];
  return s;
}
inline double sumsq3(const double* x) {
  double s = 0.0;
  s = s + x[0] * x[0];
  s = s + x[1] * x[1];
  s = s + x[2] * x[2];
  return s;
}

// MAXVAL(ABS(x(1:3)))
inline double maxabs3(const double* x) {
  double m = -std::numeric_limits<double>::max();
  for (int k = 0; k < 3; ++k) {
    const double a = std::fabs(x[k]);
    if (a > m) m = a;
  }
  return m;
}

// Fortran index helpers (1-based, periodic)
inline int im_of(int i, int N) { return 1 + (i - 2 + N) % N; }
inline int ip_of(int i, int N) { return 1 + i % N; }

using fft::Fft1d;  // MKL DFTI 1-D descriptor equivalent (fft.hpp)

// Parameters of one hydro_step and scratch space for one pencil.
struct Sweep {
  int N = 0;
  double Dgas_min = 0, Dgas_frz = 0;
  double c_limit = 0, Tmf_min = 0, Tmf_max = 0;
  double temp1_unit = 0, dt_hy = 0, ngasr = 0;
  HydroCounters* cnt = nullptr;
  Fft1d* fft = nullptr;

  // scratch (sized in init)
  std::vector<double> a, c, p;            // (N)
  std::vector<double> uh, fr, fl, uv, uvh, fu;  // (5,N)
  std::vector<double> vh;                 // (3,N)
  std::vector<double> wh, fw;             // (4,N)
  std::vector<float> u0;                  // (5,N)
  std::vector<float> wt, vx, vy, vz;      // (N+2)

  void init(int n) {
    N = n;
    a.assign(n, 0);
    c.assign(n, 0);
    p.assign(n, 0);
    uh.assign(5 * n, 0);
    fr.assign(5 * n, 0);
    fl.assign(5 * n, 0);
    uv.assign(5 * n, 0);
    uvh.assign(5 * n, 0);
    fu.assign(5 * n, 0);
    vh.assign(3 * n, 0);
    wh.assign(4 * n, 0);
    fw.assign(4 * n, 0);
    u0.assign(5 * n, 0);
    wt.assign(n + 2, 0);
    vx.assign(n + 2, 0);
    vy.assign(n + 2, 0);
    vz.assign(n + 2, 0);
  }
};

// --------------------------------------------------------------------------
// smooth_cells_1 (real(4) state) / smooth_cells_2 (real(8) state)
// Mixes cell i with its neighbours; the v argument of the original is unused.
// --------------------------------------------------------------------------
template <class T>
void smooth_cells(int i, T* u, T* w, const Sweep& S) {
  const int N = S.N;
  const int im = im_of(i, N), ip = ip_of(i, N);
  T* um = u + 5 * (im - 1);
  T* ui = u + 5 * (i - 1);
  T* up = u + 5 * (ip - 1);
  T* wm = w + 4 * (im - 1);
  T* wi = w + 4 * (i - 1);
  T* wp = w + 4 * (ip - 1);
  double du1[5], du2[5];
  for (int m = 0; m < 5; ++m) {
    du1[m] = static_cast<double>((um[m] - ui[m]) / static_cast<T>(4));
    du2[m] = static_cast<double>((up[m] - ui[m]) / static_cast<T>(4));
  }
  double r = 1 - du1[0] / gmax(static_cast<double>(um[0]), S.Dgas_min);
  for (int m = 0; m < 5; ++m) um[m] = static_cast<T>(um[m] - du1[m]);
  um[0] = static_cast<T>(gmax(static_cast<double>(um[0]), S.Dgas_min));
  for (int m = 0; m < 4; ++m) wm[m] = static_cast<T>(r * wm[m]);
  r = 1 + (du1[0] + du2[0]) / gmax(static_cast<double>(ui[0]), S.Dgas_min);
  for (int m = 0; m < 5; ++m) ui[m] = static_cast<T>(ui[m] + du1[m] + du2[m]);
  ui[0] = static_cast<T>(gmax(static_cast<double>(ui[0]), S.Dgas_min));
  for (int m = 0; m < 4; ++m) wi[m] = static_cast<T>(r * wi[m]);
  r = 1 - du2[0] / gmax(static_cast<double>(up[0]), S.Dgas_min);
  for (int m = 0; m < 5; ++m) up[m] = static_cast<T>(up[m] - du2[m]);
  up[0] = static_cast<T>(gmax(static_cast<double>(up[0]), S.Dgas_min));
  for (int m = 0; m < 4; ++m) wp[m] = static_cast<T>(r * wp[m]);
}

inline void smooth_cells_1(int i, float* u, float* w, Sweep& S) {
  ++S.cnt->smooth_cells_1;
  smooth_cells<float>(i, u, w, S);
}
inline void smooth_cells_2(int i, double* u, double* w, Sweep& S) {
  ++S.cnt->smooth_cells_2;
  smooth_cells<double>(i, u, w, S);
}

// --------------------------------------------------------------------------
// check_temp: temperature floor Tmin = 1/temp1_unit (1 K), energy
// redistributed with the two neighbours in proportion to the old T.
// --------------------------------------------------------------------------
void check_temp(float* u, const float* u0, Sweep& S) {
  const int N = S.N;
  const double Tmin = 1 / S.temp1_unit;
  for (int i = 1; i <= N; ++i) {
    float* ui = u + 5 * (i - 1);
    const double d = ui[0];
    const double k = static_cast<double>(sumsq3(ui + 1)) / d / 2;
    const double e = ui[4] - k;
    double T = e / d;
    if (T < Tmin) {
      ++S.cnt->check_temp_fix;
      const int im = im_of(i, N), ip = ip_of(i, N);
      float* um = u + 5 * (im - 1);
      float* up = u + 5 * (ip - 1);
      const float* u0m = u0 + 5 * (im - 1);
      const float* u0i = u0 + 5 * (i - 1);
      const float* u0p = u0 + 5 * (ip - 1);
      const double dm = um[0];
      const double dp = up[0];
      const double km = static_cast<double>(sumsq3(um + 1)) / dm / 2;
      const double kp = static_cast<double>(sumsq3(up + 1)) / dp / 2;
      const double em = um[4] - km;
      const double ep = up[4] - kp;
      const double d0m = u0m[0];
      const double d0 = u0i[0];
      const double d0p = u0p[0];
      const double k0m = static_cast<double>(sumsq3(u0m + 1)) / d0m / 2;
      const double k0 = static_cast<double>(sumsq3(u0i + 1)) / d0 / 2;
      const double k0p = static_cast<double>(sumsq3(u0p + 1)) / d0p / 2;
      const double e0m = u0m[4] - k0m;
      const double e0 = u0i[4] - k0;
      const double e0p = u0p[4] - k0p;
      const double T0m = gmax(e0m / d0m, Tmin);
      const double T0 = gmax(e0 / d0, Tmin);
      const double T0p = gmax(e0p / d0p, Tmin);
      const double etot = em + e + ep;
      double Tm, Tp;
      if (etot > 0) {
        const double x = etot / (dm * T0m + d * T0 + dp * T0p);
        Tm = gmax(T0m * x, Tmin);
        T = gmax(T0 * x, Tmin);
        Tp = gmax(T0p * x, Tmin);
      } else {
        Tm = T0m;
        T = T0;
        Tp = T0p;
      }
      um[4] = static_cast<float>(km + dm * Tm);
      ui[4] = static_cast<float>(k + d * T);
      up[4] = static_cast<float>(kp + dp * Tp);
    }
  }
}

// --------------------------------------------------------------------------
// euler: 2nd-order relaxing-TVD Euler sweep in the moving frame
// (the `if (.false.)` 1st-order half update of the original is dead code
// and not ported).
// --------------------------------------------------------------------------

// Freezing speed and left/right fluxes of cell j (1-based) from state U
// (float u or double uh). Returns rho.
template <class T>
inline double cell_flux(int j, const T* U, Sweep& S) {
  const T* uj = U + 5 * (j - 1);
  const double rho = uj[0];
  const double vel0 = uj[1] / rho;
  const double eth = uj[4] - static_cast<double>(sumsq3(uj + 1)) / rho / 2;
  const int jj = j - 1;
  S.p[jj] = (g - 1) * eth;
  S.a[jj] = vel0;
  S.c[jj] = std::fabs(S.a[jj]) + std::sqrt(g * std::fabs(S.p[jj] / rho));
  return rho;
}
template <class T>
inline void cell_lr(int j, const T* U, Sweep& S) {
  const T* uj = U + 5 * (j - 1);
  const int jj = j - 1;
  double* fr = &S.fr[5 * jj];
  double* fl = &S.fl[5 * jj];
  const double cpa = S.c[jj] + S.a[jj];
  const double cma = S.c[jj] - S.a[jj];
  for (int m = 0; m < 5; ++m) fr[m] = uj[m] * cpa;
  fr[1] = fr[1] + S.p[jj];
  fr[4] = fr[4] + S.p[jj] * S.a[jj];
  for (int m = 0; m < 5; ++m) fl[m] = uj[m] * cma;
  fl[1] = fl[1] - S.p[jj];
  fl[4] = fl[4] - S.p[jj] * S.a[jj];
}

// interface flux fri-fli with the van-Leer-type limiter of the original;
// dfr1/dfl1 are the "previous" differences, dfr2/dfl2 the new ones.
inline void tvd_flux(int i, int ip, const double* dfr1, const double* dfr2,
                     const double* dfl1, const double* dfl2, const Sweep& S,
                     double* fu) {
  const double* fri0 = &S.fr[5 * (i - 1)];
  const double* flp0 = &S.fl[5 * (ip - 1)];
  for (int m = 0; m < 5; ++m) {
    double dd = dfr1[m] * dfr2[m];
    const double fri = (dd > 0) ? fri0[m] + dd / (dfr1[m] + dfr2[m]) : fri0[m];
    dd = dfl1[m] * dfl2[m];
    const double fli = (dd > 0) ? flp0[m] + dd / (dfl1[m] + dfl2[m]) : flp0[m];
    fu[m] = fri - fli;
  }
}

void euler(float* u, float* v, float* w, double vmax[2], Sweep& S) {
  const int N = S.N;
  const double m01 = static_cast<double>(0.1f);
  std::copy(u, u + 5 * N, S.u0.begin());
  double dt = S.dt_hy;
  double dth = dt / 2;
  double dtq = dt / 4;
  double cmax = 0;
  int it = 0, nt = 0;
  double rho = 0;
  double fu[2][5], fw[2][4], dfr[2][5], dfl[2][5], dvg[3];
  ++S.cnt->euler_pencils;

  for (;;) {
    it = it + 1;

    // 1st order fluxes
    for (int i = 1; i <= N; ++i)
      if (u[5 * (i - 1)] < S.Dgas_min) smooth_cells_1(i, u, w, S);
    for (int i = 1; i <= N; ++i) {
      rho = cell_flux(i, u, S);
      if (S.c[i - 1] < S.c_limit) {
        cell_lr(i, u, S);
      } else {
        ++S.cnt->climit_branch;
        smooth_cells_1(i, u, w, S);
        const int n3[3] = {im_of(i, N), i, ip_of(i, N)};
        for (int k = 0; k < 3; ++k) {
          rho = cell_flux(n3[k], u, S);
          cell_lr(n3[k], u, S);
        }
      }
      cmax = gmax(cmax, S.c[i - 1]);
      if (rho > S.Dgas_frz) vmax[0] = gmax(vmax[0], S.c[i - 1]);
    }

    // number of iterations
    if (it == 1) {
      nt = f::ceiling4(cmax * dt / hydrop::cfl);
      if (nt > 1) {
        dt = dt / nt;
        dth = dt / 2;
        dtq = dt / 4;
        ++S.cnt->euler_subcycled;
      }
      S.cnt->euler_max_nt = std::max<std::int64_t>(S.cnt->euler_max_nt, nt);
      if (nt < 1)
        throw std::runtime_error(
            "euler: nt = ceiling(cmax*dt/cfl) < 1; the original loops forever here "
            "(port_notes_hydro.md H1)");
    }

    // Half update (2nd order branch)
    for (int m = 0; m < 5; ++m) fu[0][m] = 0;
    for (int m = 0; m < 4; ++m) fw[0][m] = 0;
    for (int m = 0; m < 5; ++m) {
      dfr[0][m] = S.fr[m] - S.fr[5 * (N - 1) + m];
      dfl[0][m] = S.fl[m] - S.fl[5 * 1 + m];
    }
    {
      int i, im, ip, ipp;
      for (int n = 1; n <= N; n += 2) {
        // odd cell
        i = n;
        ip = ip_of(i, N);
        ipp = ip_of(ip, N);
        im = im_of(i, N);
        for (int m = 0; m < 5; ++m) {
          dfr[1][m] = S.fr[5 * (ip - 1) + m] - S.fr[5 * (i - 1) + m];
          dfl[1][m] = S.fl[5 * (ip - 1) + m] - S.fl[5 * (ipp - 1) + m];
        }
        tvd_flux(i, ip, dfr[0], dfr[1], dfl[0], dfl[1], S, fu[1]);
        {
          const int ic = (fu[1][0] > 0) ? i : ip;
          for (int m = 0; m < 4; ++m)
            fw[1][m] = fu[1][0] * w[4 * (ic - 1) + m] / u[5 * (ic - 1)];
        }
        for (int m = 0; m < 3; ++m)
          dvg[m] = static_cast<double>((v[3 * (ip - 1) + m] - v[3 * (im - 1) + m]) / 2.0f);
        vmax[1] = gmax(vmax[1], maxabs3(dvg));
        for (int m = 0; m < 3; ++m) dvg[m] = gmin(gmax(dvg[m], -(m01 / dt)), m01 / dt) * dth;
        for (int m = 0; m < 3; ++m)
          S.vh[3 * (i - 1) + m] = v[3 * (i - 1) + m] - S.a[i - 1] * dvg[m];
        for (int m = 0; m < 5; ++m)
          S.uh[5 * (i - 1) + m] = u[5 * (i - 1) + m] - (fu[1][m] - fu[0][m]) * dtq;
        S.uh[5 * (i - 1) + 4] = S.uh[5 * (i - 1) + 4] -
                                S.p[i - 1] / (g - 1) * (1 - std::exp(-((g - 1) * dvg[0])));
        for (int m = 0; m < 4; ++m)
          S.wh[4 * (i - 1) + m] = w[4 * (i - 1) + m] - (fw[1][m] - fw[0][m]) * dtq;

        // even cell
        i = i + 1;
        im = 1 + im % N;
        ip = 1 + ip % N;
        ipp = 1 + ipp % N;
        for (int m = 0; m < 5; ++m) {
          dfr[0][m] = S.fr[5 * (ip - 1) + m] - S.fr[5 * (i - 1) + m];
          dfl[0][m] = S.fl[5 * (ip - 1) + m] - S.fl[5 * (ipp - 1) + m];
        }
        tvd_flux(i, ip, dfr[0], dfr[1], dfl[0], dfl[1], S, fu[0]);
        {
          const int ic = (fu[0][0] > 0) ? i : ip;
          for (int m = 0; m < 4; ++m)
            fw[0][m] = fu[0][0] * w[4 * (ic - 1) + m] / u[5 * (ic - 1)];
        }
        for (int m = 0; m < 3; ++m)
          dvg[m] = static_cast<double>((v[3 * (ip - 1) + m] - v[3 * (im - 1) + m]) / 2.0f);
        vmax[1] = gmax(vmax[1], maxabs3(dvg));
        for (int m = 0; m < 3; ++m) dvg[m] = gmin(gmax(dvg[m], -(m01 / dt)), m01 / dt) * dth;
        for (int m = 0; m < 3; ++m)
          S.vh[3 * (i - 1) + m] = v[3 * (i - 1) + m] - S.a[i - 1] * dvg[m];
        for (int m = 0; m < 5; ++m)
          S.uh[5 * (i - 1) + m] = u[5 * (i - 1) + m] - (fu[0][m] - fu[1][m]) * dtq;
        S.uh[5 * (i - 1) + 4] = S.uh[5 * (i - 1) + 4] -
                                S.p[i - 1] / (g - 1) * (1 - std::exp(-((g - 1) * dvg[0])));
        for (int m = 0; m < 4; ++m)
          S.wh[4 * (i - 1) + m] = w[4 * (i - 1) + m] - (fw[0][m] - fw[1][m]) * dtq;
      }
    }
    // boundary update
    for (int m = 0; m < 5; ++m) S.uh[m] = S.uh[m] + fu[0][m] * dtq;
    for (int m = 0; m < 4; ++m) S.wh[m] = S.wh[m] + fw[0][m] * dtq;

    // 2nd order fluxes
    double* uh = S.uh.data();
    double* wh = S.wh.data();
    for (int i = 1; i <= N; ++i)
      if (uh[5 * (i - 1)] < S.Dgas_min) smooth_cells_2(i, uh, wh, S);
    for (int i = 1; i <= N; ++i) {
      rho = cell_flux(i, uh, S);
      if (S.c[i - 1] < S.c_limit) {
        cell_lr(i, uh, S);
      } else {
        ++S.cnt->climit_branch;
        smooth_cells_2(i, uh, wh, S);
        const int n3[3] = {im_of(i, N), i, ip_of(i, N)};
        for (int k = 0; k < 3; ++k) {
          rho = cell_flux(n3[k], uh, S);
          cell_lr(n3[k], uh, S);
        }
      }
      if (rho > S.Dgas_frz) vmax[0] = gmax(vmax[0], S.c[i - 1]);
    }

    // Full update
    for (int m = 0; m < 5; ++m) fu[0][m] = 0;
    for (int m = 0; m < 4; ++m) fw[0][m] = 0;
    for (int m = 0; m < 5; ++m) {
      dfr[0][m] = S.fr[m] - S.fr[5 * (N - 1) + m];
      dfl[0][m] = S.fl[m] - S.fl[5 * 1 + m];
    }
    {
      int i, im, ip, ipp;
      for (int n = 1; n <= N; n += 2) {
        // odd cell
        i = n;
        ip = ip_of(i, N);
        ipp = ip_of(ip, N);
        im = im_of(i, N);
        for (int m = 0; m < 5; ++m) {
          dfr[1][m] = S.fr[5 * (ip - 1) + m] - S.fr[5 * (i - 1) + m];
          dfl[1][m] = S.fl[5 * (ip - 1) + m] - S.fl[5 * (ipp - 1) + m];
        }
        tvd_flux(i, ip, dfr[0], dfr[1], dfl[0], dfl[1], S, fu[1]);
        {
          const int ic = (fu[1][0] > 0) ? i : ip;
          for (int m = 0; m < 4; ++m) fw[1][m] = fu[1][0] * wh[4 * (ic - 1) + m] / uh[5 * (ic - 1)];
        }
        for (int m = 0; m < 3; ++m) dvg[m] = (S.vh[3 * (ip - 1) + m] - S.vh[3 * (im - 1) + m]) / 2;
        vmax[1] = gmax(vmax[1], maxabs3(dvg));
        for (int m = 0; m < 3; ++m) dvg[m] = gmin(gmax(dvg[m], -(m01 / dt)), m01 / dt) * dt;
        for (int m = 0; m < 3; ++m)
          v[3 * (i - 1) + m] = static_cast<float>(v[3 * (i - 1) + m] - S.a[i - 1] * dvg[m]);
        for (int m = 0; m < 5; ++m)
          u[5 * (i - 1) + m] =
              static_cast<float>(u[5 * (i - 1) + m] - (fu[1][m] - fu[0][m]) * dth);
        u[5 * (i - 1) + 4] = static_cast<float>(
            u[5 * (i - 1) + 4] - S.p[i - 1] / (g - 1) * (1 - std::exp(-((g - 1) * dvg[0]))));
        for (int m = 0; m < 4; ++m)
          w[4 * (i - 1) + m] =
              static_cast<float>(w[4 * (i - 1) + m] - (fw[1][m] - fw[0][m]) * dth);

        // even cell
        i = i + 1;
        im = 1 + im % N;
        ip = 1 + ip % N;
        ipp = 1 + ipp % N;
        for (int m = 0; m < 5; ++m) {
          dfr[0][m] = S.fr[5 * (ip - 1) + m] - S.fr[5 * (i - 1) + m];
          dfl[0][m] = S.fl[5 * (ip - 1) + m] - S.fl[5 * (ipp - 1) + m];
        }
        tvd_flux(i, ip, dfr[0], dfr[1], dfl[0], dfl[1], S, fu[0]);
        {
          const int ic = (fu[0][0] > 0) ? i : ip;
          for (int m = 0; m < 4; ++m) fw[0][m] = fu[0][0] * wh[4 * (ic - 1) + m] / uh[5 * (ic - 1)];
        }
        for (int m = 0; m < 3; ++m) dvg[m] = (S.vh[3 * (ip - 1) + m] - S.vh[3 * (im - 1) + m]) / 2;
        vmax[1] = gmax(vmax[1], maxabs3(dvg));
        for (int m = 0; m < 3; ++m) dvg[m] = gmin(gmax(dvg[m], -(m01 / dt)), m01 / dt) * dt;
        for (int m = 0; m < 3; ++m)
          v[3 * (i - 1) + m] = static_cast<float>(v[3 * (i - 1) + m] - S.a[i - 1] * dvg[m]);
        for (int m = 0; m < 5; ++m)
          u[5 * (i - 1) + m] =
              static_cast<float>(u[5 * (i - 1) + m] - (fu[0][m] - fu[1][m]) * dth);
        u[5 * (i - 1) + 4] = static_cast<float>(
            u[5 * (i - 1) + 4] - S.p[i - 1] / (g - 1) * (1 - std::exp(-((g - 1) * dvg[0]))));
        for (int m = 0; m < 4; ++m)
          w[4 * (i - 1) + m] =
              static_cast<float>(w[4 * (i - 1) + m] - (fw[0][m] - fw[1][m]) * dth);
      }
    }
    // boundary update
    for (int m = 0; m < 5; ++m) u[m] = static_cast<float>(u[m] + fu[0][m] * dth);
    for (int m = 0; m < 4; ++m) w[m] = static_cast<float>(w[m] + fw[0][m] * dth);

    if (it == nt) break;
  }

  check_temp(u, S.u0.data(), S);
}

// --------------------------------------------------------------------------
// advect (MF): advection of u and w with the grid velocity v(1)
// --------------------------------------------------------------------------
void advect(float* u, float* v, float* w, double vmax[2], Sweep& S) {
  const int N = S.N;
  const double m01 = static_cast<double>(0.1f);
  std::copy(u, u + 5 * N, S.u0.begin());
  double dt = S.dt_hy;
  double dth = dt / 2;
  double amax = 0;
  int it = 0, nt = 0;
  double f[5], dfl[5], dfr[5], dvg[3];
  double* uh = S.uh.data();
  double* vh = S.vh.data();
  double* wh = S.wh.data();
  double* uv = S.uv.data();
  double* uvh = S.uvh.data();
  double* fu = S.fu.data();
  double* fw = S.fw.data();
  ++S.cnt->advect_pencils;

  // interface flux at ip (between i and ip) from UV, velocity VX (stride 3)
  auto iface = [&](int i, int im, int ip, int ipp, const double* UV, auto vi, auto vip) {
    double* fup = fu + 5 * (ip - 1);
    if (vi * vip > 0) {
      if (vi > 0) {
        for (int m = 0; m < 5; ++m) {
          f[m] = UV[5 * (i - 1) + m];
          dfl[m] = UV[5 * (i - 1) + m] - UV[5 * (im - 1) + m];
          dfr[m] = UV[5 * (ip - 1) + m] - UV[5 * (i - 1) + m];
        }
      } else {
        for (int m = 0; m < 5; ++m) {
          f[m] = UV[5 * (ip - 1) + m];
          dfl[m] = UV[5 * (i - 1) + m] - UV[5 * (ip - 1) + m];
          dfr[m] = UV[5 * (ip - 1) + m] - UV[5 * (ipp - 1) + m];
        }
      }
      for (int m = 0; m < 5; ++m) {
        const double dd = dfl[m] * dfr[m];
        fup[m] = (dd > 0) ? f[m] + dd / (dfl[m] + dfr[m]) : f[m];
      }
    } else {
      for (int m = 0; m < 5; ++m) fup[m] = (UV[5 * (i - 1) + m] + UV[5 * (ip - 1) + m]) / 2;
    }
  };

  for (;;) {
    it = it + 1;

    for (int i = 1; i <= N; ++i)
      if (u[5 * (i - 1)] < S.Dgas_min) smooth_cells_1(i, u, w, S);
    for (int i = 1; i <= N; ++i) {
      const float v1 = v[3 * (i - 1)];
      for (int m = 0; m < 5; ++m) uv[5 * (i - 1) + m] = static_cast<double>(u[5 * (i - 1) + m] * v1);
      const double a = std::fabs(v1);
      amax = gmax(amax, a);
      if (u[5 * (i - 1)] > S.Dgas_frz) vmax[0] = gmax(vmax[0], a);
    }

    if (it == 1) {
      nt = f::ceiling4(amax * dt / hydrop::cfl);
      if (nt > 1) {
        dt = dt / nt;
        dth = dt / 2;
        ++S.cnt->advect_subcycled;
      }
      S.cnt->advect_max_nt = std::max<std::int64_t>(S.cnt->advect_max_nt, nt);
      if (nt < 1)
        throw std::runtime_error(
            "advect: grid velocity of a whole pencil is exactly zero (nt = 0); the "
            "original loops forever here (port_notes_hydro.md H1)");
    }

    // Half update
    for (int m = 0; m < 5; ++m) fu[m] = 0;
    for (int m = 0; m < 4; ++m) fw[m] = 0;
    for (int i = 1; i <= N; ++i) {
      const int im = im_of(i, N), ip = ip_of(i, N), ipp = ip_of(ip, N);
      iface(i, im, ip, ipp, uv, v[3 * (i - 1)], v[3 * (ip - 1)]);
      {
        const double fu1 = fu[5 * (ip - 1)];
        const int ic = (fu1 > 0) ? i : ip;
        for (int m = 0; m < 4; ++m) fw[4 * (ip - 1) + m] = fu1 * w[4 * (ic - 1) + m] / u[5 * (ic - 1)];
      }
      for (int m = 0; m < 3; ++m)
        dvg[m] = static_cast<double>((v[3 * (ip - 1) + m] - v[3 * (im - 1) + m]) / 2.0f);
      vmax[1] = gmax(vmax[1], maxabs3(dvg));
      for (int m = 0; m < 3; ++m) dvg[m] = gmin(gmax(dvg[m], -(m01 / dt)), m01 / dt) * dth;
      const float v1 = v[3 * (i - 1)];
      for (int m = 0; m < 3; ++m) vh[3 * (i - 1) + m] = v[3 * (i - 1) + m] - v1 * dvg[m];
      for (int m = 0; m < 5; ++m)
        uh[5 * (i - 1) + m] = u[5 * (i - 1) + m] - (fu[5 * (ip - 1) + m] - fu[5 * (i - 1) + m]) * dth;
      for (int m = 0; m < 4; ++m)
        wh[4 * (i - 1) + m] = w[4 * (i - 1) + m] - (fw[4 * (ip - 1) + m] - fw[4 * (i - 1) + m]) * dth;
    }
    for (int m = 0; m < 5; ++m) uh[m] = uh[m] + fu[m] * dth;
    for (int m = 0; m < 4; ++m) wh[m] = wh[m] + fw[m] * dth;

    // 2nd order fluxes
    for (int i = 1; i <= N; ++i)
      if (uh[5 * (i - 1)] < S.Dgas_min) smooth_cells_2(i, uh, wh, S);
    for (int i = 1; i <= N; ++i) {
      const double v1 = vh[3 * (i - 1)];
      for (int m = 0; m < 5; ++m) uvh[5 * (i - 1) + m] = uh[5 * (i - 1) + m] * v1;
      const double a = std::fabs(v1);
      if (uh[5 * (i - 1)] > S.Dgas_frz) vmax[0] = gmax(vmax[0], a);
    }

    // Full update
    for (int m = 0; m < 5; ++m) fu[m] = 0;
    for (int m = 0; m < 4; ++m) fw[m] = 0;
    for (int i = 1; i <= N; ++i) {
      const int im = im_of(i, N), ip = ip_of(i, N), ipp = ip_of(ip, N);
      iface(i, im, ip, ipp, uvh, vh[3 * (i - 1)], vh[3 * (ip - 1)]);
      {
        const double fu1 = fu[5 * (ip - 1)];
        const int ic = (fu1 > 0) ? i : ip;
        for (int m = 0; m < 4; ++m) fw[4 * (ip - 1) + m] = fu1 * wh[4 * (ic - 1) + m] / uh[5 * (ic - 1)];
      }
      for (int m = 0; m < 3; ++m) dvg[m] = (vh[3 * (ip - 1) + m] - vh[3 * (im - 1) + m]) / 2;
      vmax[1] = gmax(vmax[1], maxabs3(dvg));
      for (int m = 0; m < 3; ++m) dvg[m] = gmin(gmax(dvg[m], -(m01 / dt)), m01 / dt) * dt;
      const double v1 = vh[3 * (i - 1)];
      for (int m = 0; m < 3; ++m)
        v[3 * (i - 1) + m] = static_cast<float>(v[3 * (i - 1) + m] - v1 * dvg[m]);
      for (int m = 0; m < 5; ++m)
        u[5 * (i - 1) + m] = static_cast<float>(
            u[5 * (i - 1) + m] - (fu[5 * (ip - 1) + m] - fu[5 * (i - 1) + m]) * dt);
      for (int m = 0; m < 4; ++m)
        w[4 * (i - 1) + m] = static_cast<float>(
            w[4 * (i - 1) + m] - (fw[4 * (ip - 1) + m] - fw[4 * (i - 1) + m]) * dt);
    }
    for (int m = 0; m < 5; ++m) u[m] = static_cast<float>(u[m] + fu[m] * dt);
    for (int m = 0; m < 4; ++m) w[m] = static_cast<float>(w[m] + fw[m] * dt);

    if (it == nt) break;
  }

  check_temp(u, S.u0.data(), S);
}

// --------------------------------------------------------------------------
// smooth (MF): grid velocity = Gaussian-smoothed (rmf cells) total velocity,
// weighted by 1/sqrt(T) clamped to [Tmf_min,Tmf_max]; frame change of u.
// --------------------------------------------------------------------------
void smooth(float* u, float* v, float* w, Sweep& S) {
  const int N = S.N;
  float* wt = S.wt.data();
  float* vx = S.vx.data();
  float* vy = S.vy.data();
  float* vz = S.vz.data();

  auto weights = [&](int j, double tmp, const double* vel) {
    wt[j - 1] = static_cast<float>(1 / std::sqrt(tmp));
    vx[j - 1] = static_cast<float>(vel[0] * wt[j - 1]);
    vy[j - 1] = static_cast<float>(vel[1] * wt[j - 1]);
    vz[j - 1] = static_cast<float>(vel[2] * wt[j - 1]);
  };

  for (int i = 1; i <= N; ++i)
    if (u[5 * (i - 1)] < S.Dgas_min) smooth_cells_1(i, u, w, S);
  for (int i = 1; i <= N; ++i) {
    const float* ui = u + 5 * (i - 1);
    const double rho = ui[0];
    double vel[3];
    for (int m = 0; m < 3; ++m) vel[m] = ui[1 + m] / rho + v[3 * (i - 1) + m];
    const double eth = ui[4] - static_cast<double>(sumsq3(ui + 1)) / rho / 2;
    double tmp = eth / rho;
    const double a = std::fabs(vel[0]);
    if (a < S.c_limit) {
      tmp = gmin(gmax(tmp, S.Tmf_min), S.Tmf_max);
      weights(i, tmp, vel);
    } else {
      ++S.cnt->climit_branch;
      smooth_cells_1(i, u, w, S);
      const int n3[3] = {im_of(i, N), i, ip_of(i, N)};
      for (int k = 0; k < 3; ++k) {
        const int j = n3[k];
        const float* uj = u + 5 * (j - 1);
        const double rhoj = uj[0];
        double velj[3];
        for (int m = 0; m < 3; ++m) velj[m] = uj[1 + m] / rhoj + v[3 * (j - 1) + m];
        const double ethj = uj[4] - static_cast<double>(sumsq3(uj + 1)) / rhoj / 2;
        const double tmpj = gmin(gmax(ethj / rhoj, S.Tmf_min), S.Tmf_max);
        weights(j, tmpj, velj);
      }
    }
  }
  for (int i = N; i < N + 2; ++i) wt[i] = vx[i] = vy[i] = vz[i] = 0;

  S.fft->forward(wt);
  S.fft->forward(vx);
  S.fft->forward(vy);
  S.fft->forward(vz);

  // Gaussian kernel; pi*rmf is a compile-time constant in Fortran (exact)
  const double pirmf = pi * hydrop::rmf;
  for (int i = 1; i <= N + 1; i += 2) {
    const double t = pirmf * (i - 1) / S.ngasr;
    const double wg = std::exp(-(t * t) / 2);
    for (float* arr : {wt, vx, vy, vz}) {
      arr[i - 1] = static_cast<float>(arr[i - 1] * wg);
      arr[i] = static_cast<float>(arr[i] * wg);
    }
  }

  S.fft->backward(wt);
  S.fft->backward(vx);
  S.fft->backward(vy);
  S.fft->backward(vz);

  // Frame change
  for (int i = 1; i <= N; ++i) {
    float* ui = u + 5 * (i - 1);
    float* vi = v + 3 * (i - 1);
    vx[i - 1] = vx[i - 1] / wt[i - 1];
    vy[i - 1] = vy[i - 1] / wt[i - 1];
    vz[i - 1] = vz[i - 1] / wt[i - 1];
    double dvg[3];
    dvg[0] = static_cast<double>(vx[i - 1] - vi[0]);
    dvg[1] = static_cast<double>(vy[i - 1] - vi[1]);
    dvg[2] = static_cast<double>(vz[i - 1] - vi[2]);
    vi[0] = vx[i - 1];
    vi[1] = vy[i - 1];
    vi[2] = vz[i - 1];

    const double rho = ui[0];
    double dv[3];
    for (int m = 0; m < 3; ++m) dv[m] = ui[1 + m] / rho;
    const double eth = ui[4] - rho * sumsq3(dv) / 2;
    ui[0] = static_cast<float>(rho);
    for (int m = 0; m < 3; ++m) ui[1 + m] = static_cast<float>(rho * (dv[m] - dvg[m]));
    ui[4] = static_cast<float>(static_cast<double>(sumsq3(ui + 1)) / rho / 2 + eth);
  }
}

// --------------------------------------------------------------------------
// freeze: smooth cells that are below Dgas_min or faster than c_limit
// --------------------------------------------------------------------------
void freeze(float* u, float* w, Sweep& S) {
  const int N = S.N;
  for (int i = 1; i <= N; ++i)
    if (u[5 * (i - 1)] < S.Dgas_min) smooth_cells_1(i, u, w, S);
  for (int i = 1; i <= N; ++i) {
    const float* ui = u + 5 * (i - 1);
    const double rho = ui[0];
    const double a = ui[1] / rho;
    const double eth = ui[4] - static_cast<double>(sumsq3(ui + 1)) / rho / 2;
    const double p = (g - 1) * eth;
    const double c = std::fabs(a) + std::sqrt(g * std::fabs(p / rho));
    if (c > S.c_limit) {
      ++S.cnt->freeze_smooth;
      smooth_cells_1(i, u, w, S);
    }
  }
}

// --------------------------------------------------------------------------
// Fortran ESw.d formatting
// --------------------------------------------------------------------------
std::string right(const std::string& s, int w) {
  if (static_cast<int>(s.size()) > w) return std::string(w, '*');
  return std::string(w - s.size(), ' ') + s;
}

}  // namespace

std::string format_es(double x, int w, int d) {
  if (std::isnan(x)) return right("NaN", w);
  if (std::isinf(x)) return right(x > 0 ? "Infinity" : "-Infinity", w);
  char buf[64];
  std::snprintf(buf, sizeof buf, "%.*E", d, x);
  std::string s(buf);
  const auto e = s.find('E');
  std::string mant = s.substr(0, e);
  int ex = std::atoi(s.c_str() + e + 1);
  char eb[16];
  if (std::abs(ex) <= 99)
    std::snprintf(eb, sizeof eb, "E%c%02d", ex < 0 ? '-' : '+', std::abs(ex));
  else
    std::snprintf(eb, sizeof eb, "%c%03d", ex < 0 ? '-' : '+', std::abs(ex));
  return right(mant + eb, w);
}

std::string format_gas_record(const std::vector<double>& rec) {
  std::string s;
  for (double x : rec) s += format_es(x, 13, 5);
  return s;
}

// --------------------------------------------------------------------------
// init_hydro
// --------------------------------------------------------------------------
void init_hydro(const Config& cfg, GasGrid& gas, HydroState& hs) {
  gas.allocate(cfg.Ngas);
  hs.allocate(cfg.Ngas);
}

// --------------------------------------------------------------------------
// hydro_step
// --------------------------------------------------------------------------
void hydro_step(const Config& cfg, const Units& un, const SimClock& clk,
                SweepOrder order, GasGrid& gas, HydroState& hs) {
  const int N = cfg.Ngas;
  if (N % 2 != 0)
    throw std::runtime_error("hydro_step: Ngas must be even (odd/even cell pairs in euler)");

  // There are two hydro steps per PM step
  hs.dt_hy = clk.dt_sim / 2;

  double vmax[2] = {0, 0};
  const float bscale = 1.0f / static_cast<float>(N);

  // Calc freezing speed limit
  hs.c_limit = gmax(2 * clk.cmax_gas, static_cast<double>(2e8f) / un.vel_unit);

  // The pencils (j,k) are independent: each touches only its own cells.
  // One sweep workspace, FFT descriptor (planned serially: the FFTW planner
  // is not thread-safe) and set of counters per thread; the counters are
  // integer sums/maxima and vmax a maximum, so the result is the same for
  // any thread count.
  const int nth = par::max_threads();
  std::vector<std::unique_ptr<Fft1d>> ffts;
  std::vector<Sweep> sweeps(nth);
  std::vector<HydroCounters> cnts(nth);
  std::vector<std::array<double, 2>> vmaxs(nth, {0.0, 0.0});
  for (int t = 0; t < nth; ++t) {
    ffts.emplace_back(new Fft1d(N, bscale));
    Sweep& S = sweeps[t];
    S.init(N);
    S.N = N;
    S.Dgas_min = un.Dgas_min;
    S.Dgas_frz = un.Dgas_frz;
    S.c_limit = hs.c_limit;
    S.Tmf_min = un.Tmf_min;
    S.Tmf_max = un.Tmf_max;
    S.temp1_unit = un.temp1_unit;
    S.dt_hy = hs.dt_hy;
    S.ngasr = cfg.ngasr;
    S.cnt = &cnts[t];
    S.fft = ffts[t].get();
  }

  const std::int64_t npencil = static_cast<std::int64_t>(N) * N;
  std::string err;
#pragma omp parallel for schedule(dynamic, 16) if (par::worth(npencil * N))
  for (std::int64_t p = 0; p < npencil; ++p) {
    const int t = par::thread_id();
    Sweep& S = sweeps[t];
    double* vm = vmaxs[t].data();
    const int k = static_cast<int>(p / N) + 1, j = static_cast<int>(p % N) + 1;
    try {
      const std::int64_t c0 = gas.cell(1, j, k);
      float* u = gas.u_p(c0);
      float* w = gas.w_p(c0);
      float* v = hs.v_p(c0);
      smooth(u, v, w, S);
      if (order == SweepOrder::Forward) {
        euler(u, v, w, vm, S);
        advect(u, v, w, vm, S);
      } else {
        advect(u, v, w, vm, S);
        euler(u, v, w, vm, S);
      }
      freeze(u, w, S);
    } catch (const std::exception& e) {
#pragma omp critical(radray_hydro_err)
      if (err.empty()) err = e.what();
    }
  }
  if (!err.empty()) throw std::runtime_error(err);
  HydroCounters& C = hs.counters;
  for (int t = 0; t < nth; ++t) {
    vmax[0] = gmax(vmax[0], vmaxs[t][0]);
    vmax[1] = gmax(vmax[1], vmaxs[t][1]);
    const HydroCounters& c = cnts[t];
    C.smooth_cells_1 += c.smooth_cells_1;
    C.smooth_cells_2 += c.smooth_cells_2;
    C.climit_branch += c.climit_branch;
    C.freeze_smooth += c.freeze_smooth;
    C.check_temp_fix += c.check_temp_fix;
    C.euler_pencils += c.euler_pencils;
    C.euler_subcycled += c.euler_subcycled;
    C.euler_max_nt = std::max(C.euler_max_nt, c.euler_max_nt);
    C.advect_pencils += c.advect_pencils;
    C.advect_subcycled += c.advect_subcycled;
    C.advect_max_nt = std::max(C.advect_max_nt, c.advect_max_nt);
  }
  hs.vmax_last = {vmax[0], vmax[1]};
}

// --------------------------------------------------------------------------
// transpose_ij / transpose_ik (blocked loops of the original, nb = 8)
// --------------------------------------------------------------------------
namespace {
constexpr int nb = 8;

// swap a(:,i,j,k) and a(:,j,i,k) with component permutation perm (1-based)
template <int NC>
void swap_perm(float* x, float* y, const int (&perm)[NC]) {
  float xi[NC], yj[NC];
  for (int m = 0; m < NC; ++m) {
    xi[m] = x[m];
    yj[m] = y[m];
  }
  for (int m = 0; m < NC; ++m) x[m] = yj[perm[m] - 1];
  for (int m = 0; m < NC; ++m) y[m] = xi[perm[m] - 1];
}

// generic blocked transpose of a(NC,N,N,N) in the (i,j) or (i,k) plane
template <int NC, bool IK>
void transpose_field(std::vector<float>& a, int N, const int (&perm)[NC]) {
  auto idx = [&](int i, int j, int k) {
    return static_cast<std::size_t>(NC) *
           ((i - 1) + static_cast<std::size_t>(N) * ((j - 1) + static_cast<std::size_t>(N) * (k - 1)));
  };
  // (i,j,k) <-> (j,i,k)  or  (i,j,k) <-> (k,j,i);  outer loop over the
  // untouched index l, block loops over the swapped pair (b = j or k)
  auto pair = [&](int i, int b, int l) {
    if (IK)
      swap_perm<NC>(&a[idx(i, l, b)], &a[idx(b, l, i)], perm);
    else
      swap_perm<NC>(&a[idx(i, b, l)], &a[idx(b, i, l)], perm);
  };
  // each plane l is independent
#pragma omp parallel for schedule(static) if (par::worth(static_cast<long long>(N) * N * N))
  for (int l = 1; l <= N; ++l) {
    for (int b1 = 1; b1 <= N; b1 += nb) {
      for (int b = b1; b <= b1 + nb - 1; ++b)
        for (int i = b; i <= b1 + nb - 1; ++i) pair(i, b, l);
      for (int i1 = b1 + nb; i1 <= N; i1 += nb)
        for (int b = b1; b <= b1 + nb - 1; ++b)
          for (int i = i1; i <= i1 + nb - 1; ++i) pair(i, b, l);
    }
  }
}

void check_nb(int N) {
  if (N % nb != 0)
    throw std::runtime_error(
        "transpose: Ngas must be a multiple of 8 (blocked loops of the original "
        "run out of bounds otherwise)");
}
}  // namespace

void transpose_ij(GasGrid& gas, HydroState& hs) {
  check_nb(gas.N);
  static const int pu[5] = {1, 3, 2, 4, 5};
  static const int pv[3] = {2, 1, 3};
  static const int pw[4] = {1, 2, 3, 4};
  transpose_field<5, false>(gas.u, gas.N, pu);
  transpose_field<3, false>(hs.v, gas.N, pv);
  transpose_field<4, false>(gas.w, gas.N, pw);
}

void transpose_ik(GasGrid& gas, HydroState& hs) {
  check_nb(gas.N);
  static const int pu[5] = {1, 4, 3, 2, 5};
  static const int pv[3] = {3, 2, 1};
  static const int pw[4] = {1, 2, 3, 4};
  transpose_field<5, true>(gas.u, gas.N, pu);
  transpose_field<3, true>(hs.v, gas.N, pv);
  transpose_field<4, true>(gas.w, gas.N, pw);
}

// --------------------------------------------------------------------------
// summarize_gas
// --------------------------------------------------------------------------
namespace {
// HII case-B recombination coefficient as written in summarize_gas
inline double alpha_B(double x) {
  return 2.753e-14 * std::pow(x, static_cast<double>(1.5f)) /
         std::pow(1.0 + std::pow(x / static_cast<double>(2.740f), static_cast<double>(0.407f)),
                  static_cast<double>(2.242f));
}
}  // namespace

GasSummary summarize_gas(const Config& cfg, const Units& un, SimClock& clk,
                         GasGrid& gas, HydroState& hs, std::ostream* out) {
  const double ob = cfg.cosmo.omegab, om = cfg.cosmo.omegam;
  const double fb = ob / om;
  const double Zsolar = static_cast<double>(0.02f);
  const int N = cfg.Ngas;
  const double Ncell = static_cast<double>(cfg.Ncell);

  double x1[11], x2[11], x3[11], x4[11];
  for (int m = 0; m < 11; ++m) {
    x1[m] = 0;
    x2[m] = 0;
    x3[m] = std::numeric_limits<double>::max();
    x4[m] = 0;
  }
  x3[1] = x3[2] = x3[3] = 0;
  double cmax = 0, vrms = 0;
  double clump[21];  // 1-based
  for (double& c : clump) c = 0;

  auto acc = [&](int m, double val, double rho) {
    x1[m] = x1[m] + val;
    x2[m] = x2[m] + val * rho;
    x3[m] = gmin(x3[m], val);
    x4[m] = gmax(x4[m], val);
  };

  // Per-cell values are computed in parallel (one k-plane batch at a time,
  // cells independent; the clamps of u_gas/w_gas are per cell) and summed
  // serially in the original k,j,i order with the original expressions, so
  // every sum is bit-identical for any thread count. cmax is a maximum.
  struct Rec {
    double rho, vel[3], cs, tmp, fHI, fHeI, fHeII, fe, Z, nH, nHII, ne, alpha, x, cmax;
  };
  const std::int64_t plane = static_cast<std::int64_t>(N) * N;
  const int nth = par::worth(cfg.Ncell) ? par::max_threads() : 1;
  const int batch = nth > 1 ? std::min(N, 4 * nth) : 1;
  std::vector<Rec> recs(static_cast<std::size_t>(batch * plane));
  for (int k0 = 1; k0 <= N; k0 += batch) {
    const int nb = std::min(batch, N - k0 + 1);
#pragma omp parallel for schedule(static) if (nth > 1)
    for (int b = 0; b < nb; ++b) {
      const int k = k0 + b;
      for (int j = 1; j <= N; ++j)
        for (int i = 1; i <= N; ++i) {
          const std::int64_t c = gas.cell(i, j, k);
          float* u = gas.u_p(c);
          float* w = gas.w_p(c);
          const float* v = hs.v_p(c);
          Rec& R = recs[b * plane + (i - 1) + static_cast<std::int64_t>(N) * (j - 1)];

          // Hydro u_gas variables
          u[0] = static_cast<float>(gmax(static_cast<double>(u[0]), un.Dgas_min));
          const double rho = u[0];
          double vel[3];
          for (int m = 0; m < 3; ++m) vel[m] = u[1 + m] / rho + v[m];
          double eth = u[4] - static_cast<double>(sumsq3(u + 1)) / rho / 2;
          // hydro_scheme = unsplit: the dual-energy field (empty otherwise)
          if (!gas.e.empty()) eth = gas.e[c];
          eth = gmax(eth, 0.0);
          const double cs = std::sqrt(g * (g - 1) * eth / rho);

          // Chemistry w_gas variables
          for (int m = 0; m < 4; ++m) w[m] = gmaxf(w[m], 0.0f);
          const double fHI = w[0] / (XHy * rho);
          const double fHeI = w[1] / (YHe * rho);
          const double fHeII = w[2] / (YHe * rho);
          const double Z = w[3] / rho;

          // Additional variables
          const double fHII = gmin(gmax(1 - fHI, 0.0), 1.0);
          const double fHeIII = gmin(gmax(1 - fHeI - fHeII, 0.0), 1.0);
          const double fe = (fHII * XHy + (fHeII + 2 * fHeIII) * YHe / 4) / (XHy + 2 * YHe / 4.0);
          const double feff = ((1 + fHII) * XHy + (1 + fHeII + 2 * fHeIII) * YHe / 4) /
                              (2 * XHy + 3 * YHe / 4.0);
          const double tmp = eth / (rho * feff);

          // Max freezing speed (of this cell)
          double cm = 0;
          if (rho > un.Dgas_frz) {
            const double c1 = maxabs3(vel);
            double dvv[3];
            for (int m = 0; m < 3; ++m) dvv[m] = vel[m] - v[m];
            const double c2 = maxabs3(dvv) + cs;
            cm = gmax(gmax(cm, c1), c2);
          }

          // HII recombination case B
          const double x = 2 * HI_K / (tmp * un.temp3_unit);
          const double alpha = alpha_B(x);

          // Clumping factor
          const double nH = rho * un.H_unit;
          const double nHe = rho * un.He_unit;
          const double nHII = fHII * nH;
          const double ne = fe * (nH + 2 * nHe);
          R = Rec{rho, {vel[0], vel[1], vel[2]}, cs, tmp, fHI, fHeI, fHeII, fe, Z, nH, nHII, ne,
                  alpha, gmin(std::exp(-(rho / fb - rho_source) / rho_source), 1.0), cm};
        }
    }
    for (std::int64_t n = 0; n < nb * plane; ++n) {
      const Rec& R = recs[n];
      const double rho = R.rho;
      const double* vel = R.vel;
      acc(0, rho, rho);
      for (int m = 0; m < 3; ++m) acc(1 + m, vel[m], rho);
      acc(4, R.cs, rho);
      acc(5, R.tmp, rho);
      acc(6, R.fHI, rho);
      acc(7, R.fHeI, rho);
      acc(8, R.fHeII, rho);
      acc(9, R.fe, rho);
      acc(10, R.Z, rho);

      // Velocity dispersion
      vrms = vrms + sumsq3(vel);

      // Max freezing speed
      if (rho > un.Dgas_frz) cmax = gmax(cmax, R.cmax);

      const double nH = R.nH, nHII = R.nHII, ne = R.ne, alpha = R.alpha;
      clump[1] += nHII;
      clump[2] += nH * nH;
      clump[3] += nHII * nHII;
      clump[4] += nHII * ne * alpha;
      if (rho / fb < 50) {
        clump[5] += nHII;
        clump[6] += nH * nH;
        clump[7] += nHII * nHII;
        clump[8] += nHII * ne * alpha;
      }
      if (rho / fb < 100) {
        clump[9] += nHII;
        clump[10] += nH * nH;
        clump[11] += nHII * nHII;
        clump[12] += nHII * ne * alpha;
      }
      if (rho / fb < 200) {
        clump[13] += nHII;
        clump[14] += nH * nH;
        clump[15] += nHII * nHII;
        clump[16] += nHII * ne * alpha;
      }
      const double x = R.x;
      clump[17] += x * nHII;
      clump[18] += x * (nH * nH);
      clump[19] += x * (nHII * nHII);
      clump[20] += x * nHII * ne * alpha;
    }
  }

  // Update global variables
  clk.cmax_gas = cmax;
  clk.vmax_gas = gmax(gmax(x4[1], x4[2]) /* maxval(x4(2:3)): vz omitted */,
                     gmax(gmax(std::fabs(x3[1]), std::fabs(x3[2])), std::fabs(x3[3])));
  clk.vrms_gas = std::sqrt(vrms / Ncell);

  const double v5 = static_cast<double>(1e5f);
  GasSummary S;
  auto& t = S.table;
  // Density in rho/<rho>
  t[0] = {x1[0] / fb / Ncell, x2[0] / fb / x1[0], x3[0] / fb, x4[0] / fb};
  // Velocity in km/s (normalisations as in the original: swapped)
  for (int m = 0; m < 3; ++m)
    t[1 + m] = {x1[1 + m] * un.vel_unit / v5 / x1[0], x2[1 + m] * un.vel_unit / v5 / Ncell,
                x3[1 + m] * un.vel_unit / v5, x4[1 + m] * un.vel_unit / v5};
  // Sound speed in km/s
  t[4] = {x1[4] * un.cs_unit / v5 / Ncell, x2[4] * un.cs_unit / v5 / x1[0],
          x3[4] * un.cs_unit / v5, x4[4] * un.cs_unit / v5};
  // Temperature in K
  t[5] = {x1[5] * un.temp3_unit / Ncell, x2[5] * un.temp3_unit / x1[0], x3[5] * un.temp3_unit,
          x4[5] * un.temp3_unit};
  // fractions
  for (int m = 6; m <= 9; ++m) t[m] = {x1[m] / Ncell, x2[m] / x1[0], x3[m], x4[m]};
  // Metallicity in solar units
  t[10] = {x1[10] / Zsolar / Ncell, x2[10] / Zsolar / x1[0], x3[10] / Zsolar, x4[10] / Zsolar};

  // Clumping factors
  {
    const double nH = un.nH_avg;
    const double ne = un.nH_avg + un.nHe_avg;
    const double x = 2 * HI_K / static_cast<double>(1e4f);
    const double alpha = alpha_B(x);
    for (int m = 1; m <= 20; ++m) clump[m] = clump[m] / Ncell;
    for (int b = 0; b < 5; ++b) {
      const int o = 4 * b;  // clump(1+o) = nHII sum of this density cut
      clump[2 + o] = clump[2 + o] / (nH * nH);
      clump[3 + o] = clump[3 + o] / (nH * clump[1 + o]);
      clump[4 + o] = clump[4 + o] / (alpha * ne * clump[1 + o]);
    }
  }

  S.record = {clk.redshift, t[6][0], t[6][1], t[7][0], t[7][1], t[8][0], t[8][1],
              t[9][0], t[9][1], t[5][0], t[5][1], t[10][0], t[10][1]};
  for (int b = 0; b < 5; ++b)
    for (int m = 2; m <= 4; ++m) S.record.push_back(clump[m + 4 * b]);

  if (out) {
    static const char* names[11] = {"rho   : ", "vx    : ", "vy    : ", "vz    : ",
                                    "cs    : ", "T     : ", "HI    : ", "HeI   : ",
                                    "HeII  : ", "e     : ", "Z     : "};
    char buf[256];
    *out << " Gas\n";
    for (int m = 0; m < 11; ++m) {
      std::snprintf(buf, sizeof buf, " %s %14.7g %14.7g %14.7g %14.7g\n", names[m],
                    static_cast<float>(t[m][0]), static_cast<float>(t[m][1]),
                    static_cast<float>(t[m][2]), static_cast<float>(t[m][3]));
      *out << buf;
    }
    static const char* cn[5] = {"C   0 : ", "C  50 : ", "C 100 : ", "C 200 : ", "C eff : "};
    for (int b = 0; b < 5; ++b) {
      std::snprintf(buf, sizeof buf, " %s %14.7g %14.7g %14.7g\n", cn[b],
                    static_cast<float>(clump[2 + 4 * b]), static_cast<float>(clump[3 + 4 * b]),
                    static_cast<float>(clump[4 + 4 * b]));
      *out << buf;
    }
  }
  return S;
}

// --------------------------------------------------------------------------
// radhydro.f90: radhydro (hydro part)
// --------------------------------------------------------------------------
void radhydro_hydro_step(const Config& cfg, SimClock& clk, Units& un, GasGrid& gas,
                         HydroState& hs) {
  // Forward hydro step with transposes in between (legacy), or one unsplit
  // step over dt_sim/2
  hydro_half_step(cfg, un, clk, HydroHalf::First, gas, hs);

  clk.year = clk.year + clk.dt_yr / 2;
  clk.scalefactor = clk.scalefactor + clk.da_sim / 2;
  clk.redshift = 1 / clk.scalefactor - 1;
  calc_sim_var(cfg, clk, un);

  // particle_mesh: not part of this port (no gravity)

  // Reverse hydro step with transposes in between (legacy), or the second
  // unsplit step over dt_sim/2
  hydro_half_step(cfg, un, clk, HydroHalf::Second, gas, hs);

  clk.year = clk.year + clk.dt_yr / 2;
  clk.scalefactor = clk.scalefactor + clk.da_sim / 2;
  clk.redshift = 1 / clk.scalefactor - 1;
  calc_sim_var(cfg, clk, un);
}

}  // namespace radray
