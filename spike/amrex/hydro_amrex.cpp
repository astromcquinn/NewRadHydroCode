// AMReX prototype of the spike hydro step: the unsplit MUSCL-Hancock + HLLC
// update of spike/serial/hydro_serial.cpp (kernels from
// spike/common/hydro_mh.hpp) on a periodic amrex::Geometry, with the state
// in a MultiFab (5 components, 2 ghosts). Must be bitwise identical to
// hydro_serial for any rank count, max_grid_size and thread count.
//
// usage: hydro_amrex [--n 64] [--steps 20] [--ic blast|sod] [--cfl 0.4]
//                    [--out U.bin] [--max_grid_size 32]
// Output: raw float64, (5, n, n, n), component fastest, i fastest (rank 0).
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include <AMReX.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Reduce.H>

#include "hydro_mh.hpp"

using namespace amrex;

// Predictor states: by default AoS, one QCell (20 doubles) per cell in a
// FabArray<BaseFab<QCell>>, which is what the kernels (contiguous Q[NQ]) want;
// -DHYDRO_Q_SOA stores them as a 20-component MultiFab (AMReX's native SoA).
struct QCell {
  double v[rrh::NQ];
};

int main(int argc, char** argv) {
  // Our own command line; keep AMReX quiet and do not let it parse argv.
  amrex::Initialize(argc, argv, false, MPI_COMM_WORLD, []() {
    ParmParse pp("amrex");
    pp.add("verbose", 0);
  });
  int rc = 0;
  {
    int n = 64, steps = 20, ic = 0, mgs = 32;
    double cfl = 0.4;
    std::string out;
    for (int a = 1; a < argc; ++a) {
      const std::string s = argv[a];
      if (s == "--n") n = std::atoi(argv[++a]);
      else if (s == "--steps") steps = std::atoi(argv[++a]);
      else if (s == "--ic") ic = std::string(argv[++a]) == "sod" ? 1 : 0;
      else if (s == "--cfl") cfl = std::atof(argv[++a]);
      else if (s == "--out") out = argv[++a];
      else if (s == "--max_grid_size") mgs = std::atoi(argv[++a]);
      else amrex::Abort("unknown argument " + s);
    }
    constexpr int G = rrh::NG, NV = rrh::NV;
    [[maybe_unused]] constexpr int NQ = rrh::NQ;
    const Box domain(IntVect(0), IntVect(n - 1));
    const RealBox rb({0., 0., 0.}, {1., 1., 1.});
    const Array<int, 3> per{1, 1, 1};
    Geometry geom(domain, rb, CoordSys::cartesian, per);
    BoxArray ba(domain);
    ba.maxSize(mgs);
    DistributionMapping dm(ba);
    MultiFab U(ba, dm, NV, G);
#ifdef HYDRO_Q_SOA
    MultiFab Q(ba, dm, NQ, 1);
#else
    FabArray<BaseFab<QCell>> Q(ba, dm, 1, 1);
#endif
    Array<MultiFab, 3> F;
    for (int d = 0; d < 3; ++d) F[d].define(amrex::convert(ba, IntVect::TheDimensionVector(d)), dm, NV, 0);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(U, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      auto u = U.array(mfi);
      ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        double v[NV];
        rrh::initial_state(ic, i, j, k, n, v);
        for (int m = 0; m < NV; ++m) u(i, j, k, m) = v[m];
      });
    }

    const double dx = 1.0 / n;
    double t = 0;
    ParallelDescriptor::Barrier();
    const auto t0 = std::chrono::steady_clock::now();
    auto clk = []() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); };
    double Tfill = 0, Tred = 0, Tpred = 0, Tflux = 0, Tupd = 0;
    for (int s = 0; s < steps; ++s) {
      double c0 = clk();
      U.FillBoundary(geom.periodicity());
      double c1 = clk();
      Tfill += c1 - c0;
      // global max signal speed (max is exact: order-independent)
      ReduceOps<ReduceOpMax> rop;
      ReduceData<double> rdat(rop);
      for (MFIter mfi(U, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto u = U.const_array(mfi);
        rop.eval(mfi.tilebox(), rdat, [=] AMREX_GPU_DEVICE(int i, int j, int k) -> GpuTuple<double> {
          double v[NV];
          for (int m = 0; m < NV; ++m) v[m] = u(i, j, k, m);
          return {rrh::max_speed(v)};
        });
      }
      double smax = amrex::get<0>(rdat.value(rop));
      ParallelDescriptor::ReduceRealMax(smax);
      const double dt = cfl * dx / smax;
      const double dtdx = dt / dx;
      double c2 = clk();
      Tred += c2 - c1;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(U, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto u = U.const_array(mfi);
        auto q = Q.array(mfi);
        // predictor on the tile plus one ghost layer
        ParallelFor(mfi.growntilebox(1), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
          auto UA = [&](int ii, int jj, int kk, int m) { return u(ii, jj, kk, m); };
#ifdef HYDRO_Q_SOA
          double qv[NQ];
          rrh::predict_cell(UA, i, j, k, dtdx, qv);
          for (int m = 0; m < NQ; ++m) q(i, j, k, m) = qv[m];
#else
          rrh::predict_cell(UA, i, j, k, dtdx, q(i, j, k).v);
#endif
        });
      }
      double c3 = clk();
      Tpred += c3 - c2;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(U, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto q = Q.const_array(mfi);
        for (int d = 0; d < 3; ++d) {
          auto f = F[d].array(mfi);
          const int ei = d == 0, ej = d == 1, ek = d == 2;
          ParallelFor(mfi.nodaltilebox(d), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            double fv[NV];
#ifdef HYDRO_Q_SOA
            double QL[NQ], QR[NQ];
            for (int m = 0; m < NQ; ++m) {
              QL[m] = q(i - ei, j - ej, k - ek, m);
              QR[m] = q(i, j, k, m);
            }
            rrh::face_flux(QL, QR, d, fv);
#else
            rrh::face_flux(q(i - ei, j - ej, k - ek).v, q(i, j, k).v, d, fv);
#endif
            for (int m = 0; m < NV; ++m) f(i, j, k, m) = fv[m];
          });
        }
      }
      double c4 = clk();
      Tflux += c4 - c3;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(U, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto u = U.array(mfi);
        auto fx = F[0].const_array(mfi), fy = F[1].const_array(mfi), fz = F[2].const_array(mfi);
        ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
          double uv[NV], f[6][NV];
          for (int m = 0; m < NV; ++m) {
            uv[m] = u(i, j, k, m);
            f[0][m] = fx(i, j, k, m);
            f[1][m] = fx(i + 1, j, k, m);
            f[2][m] = fy(i, j, k, m);
            f[3][m] = fy(i, j + 1, k, m);
            f[4][m] = fz(i, j, k, m);
            f[5][m] = fz(i, j, k + 1, m);
          }
          rrh::update_cell(uv, f[0], f[1], f[2], f[3], f[4], f[5], dtdx);
          for (int m = 0; m < NV; ++m) u(i, j, k, m) = uv[m];
        });
      }
      Tupd += clk() - c4;
      t += dt;
    }
    ParallelDescriptor::Barrier();
    const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

    // gather to rank 0 (single box), then diagnostics/output in global order
    BoxArray ba0(domain);
    DistributionMapping dm0(Vector<int>{0});
    MultiFab U0(ba0, dm0, NV, 0);
    U0.ParallelCopy(U, 0, 0, NV);
    if (ParallelDescriptor::IOProcessor()) {
      const auto u = U0.const_array(0);
      double mass = 0, etot = 0;
      for (int k = 0; k < n; ++k)
        for (int j = 0; j < n; ++j)
          for (int i = 0; i < n; ++i) {
            mass += u(i, j, k, 0);
            etot += u(i, j, k, 4);
          }
      int nth = 1;
#ifdef AMREX_USE_OMP
      nth = omp_get_max_threads();
#endif
      const int nr = ParallelDescriptor::NProcs();
      std::printf(" n %d steps %d t %.6e mass %.17g energy %.17g\n", n, steps, t, mass * dx * dx * dx,
                  etot * dx * dx * dx);
      std::printf(" ranks %d threads %d boxes %d wall %.3f s  ns/cell-update %.2f (core-ns %.2f)\n", nr, nth,
                  static_cast<int>(ba.size()), wall, 1e9 * wall / (double(n) * n * n * steps),
                  1e9 * wall * nth * nr / (double(n) * n * n * steps));
      std::printf(" timers: FillBoundary %.3f maxspeed %.3f predict %.3f flux %.3f update %.3f\n", Tfill, Tred, Tpred,
                  Tflux, Tupd);
      if (!out.empty()) {
        std::vector<double> o(static_cast<long long>(n) * n * n * NV);
        for (int k = 0; k < n; ++k)
          for (int j = 0; j < n; ++j)
            for (int i = 0; i < n; ++i)
              for (int c = 0; c < NV; ++c)
                o[((i + static_cast<long long>(n) * (j + static_cast<long long>(n) * k)) * NV) + c] =
                    u(i, j, k, c);
        std::FILE* fp = std::fopen(out.c_str(), "wb");
        if (!fp) {
          std::fprintf(stderr, "cannot write %s\n", out.c_str());
          rc = 1;
        } else {
          std::fwrite(o.data(), sizeof(double), o.size(), fp);
          std::fclose(fp);
        }
      }
    }
  }
  amrex::Finalize();
  return rc;
}
