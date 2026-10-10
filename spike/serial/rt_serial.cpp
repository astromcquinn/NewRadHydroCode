// Serial (optionally OpenMP) reference of the sorted-ray-array RT substep
// (docs/design_parallel_and_encoding.md section 2.2), built only from the
// framework-independent kernels of spike/common. It is
//   * the reference the AMReX and Parthenon prototypes must match bit for
//     bit (any rank count, any thread count), and
//   * itself checked against the serial port (apps/rt_driver,
//     legacy_rt = 0) to round-off: same ray statistics, fields within
//     float round-off (differences only from deposit/merge summation order).
// It is also the "in-house, no framework" baseline of the benchmarks.
//
// usage: rt_serial <problem_dir> <outdir> [--nsteps n] [--set key=value]...
//                  [--no-output]
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "rt_host.hpp"

using rrk::CellData;
using rrk::RayP;

namespace {

struct Timers {
  double cell = 0, hist = 0, step1 = 0, cast1 = 0, phys = 0, cast3 = 0, sort = 0, post = 0;
};
double now() {
  return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

inline void atomic_add(long long* p, long long v) {
#ifdef _OPENMP
#pragma omp atomic
#endif
  *p += v;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 3) {
    std::fprintf(stderr, "usage: %s <problem_dir> <outdir> [--nsteps n] [--set k=v] [--no-output]\n",
                 argv[0]);
    return 2;
  }
  const std::string pdir = argv[1], odir = argv[2];
  std::vector<std::string> sets;
  int nsteps = -1;
  bool output = true;
  for (int a = 3; a < argc; ++a) {
    const std::string s = argv[a];
    if (s == "--set" && a + 1 < argc) sets.emplace_back(argv[++a]);
    else if (s == "--nsteps" && a + 1 < argc) nsteps = std::atoi(argv[++a]);
    else if (s == "--no-output") output = false;
    else {
      std::fprintf(stderr, "unknown argument %s\n", s.c_str());
      return 2;
    }
  }
  try {
    rrhost::RTRun run;
    run.init(pdir, sets);
    if (nsteps >= 0) run.P.nsteps = nsteps;
    const int N = run.P.N;
    const long long ncell = static_cast<long long>(N) * N * N;
    auto cidx = [N](int i, int j, int k) {
      return i + static_cast<long long>(N) * (j + static_cast<long long>(N) * k);
    };

    std::vector<CellData> cells(ncell);
    {
      radray::GasGrid g;
      run.read_gas(pdir, g);
      for (long long c = 0; c < ncell; ++c) cells[c] = rrhost::initial_cell(g, c);
    }
    std::vector<char> is_src(ncell, 0);
    for (const auto& s : run.sources) {
      const long long c = cidx(s.i, s.j, s.k);
      cells[c].Luv = s.Luv;
      cells[c].sfr = s.sfr;
    }

    std::vector<RayP> rays, out;                 // rays sorted by home cell
    std::vector<long long> start(ncell + 1, 0);  // CSR offsets
    std::vector<long long> acc(ncell * rrk::kNdep);
    std::vector<long long> ostart(ncell + 1);
    std::vector<int> idx;
    std::vector<rrk::MergeKey> key;
    std::vector<int> home;
    std::vector<char> keep;
    Timers T;
    long long ray_substeps = 0;
    const double t_start = now();
    int nthreads = 1;
#ifdef _OPENMP
    nthreads = omp_get_max_threads();
#endif

    for (int step = 1; step <= run.P.nsteps; ++step) {
      const int nsub = run.begin_step();
      for (long long c = 0; c < ncell; ++c) {
        is_src[c] = static_cast<double>(cells[c].sfr) > run.units.sfr_min;
      }
      {
        const radray::Units u = run.units;
#pragma omp parallel for schedule(static)
        for (long long c = 0; c < ncell; ++c) rrk::prep_cell(cells[c], u);
      }
      rrhost::StepStats st;
      for (int sub = 1; sub <= nsub; ++sub) {
        double t0 = now();
        // 1a. reset and tau
        {
          const radray::Units u = run.units;
#pragma omp parallel for schedule(static)
          for (long long c = 0; c < ncell; ++c) rrk::calc_rt_1_cell(rrk::cell_ref(cells[c]), u, run.sig_init, false);
        }
        double t1 = now();
        T.cell += t1 - t0;
        // 1b. histogram -> thresholds, pool guard
        std::vector<long long> hist(rrk::kNpdfMax, 0), spl(rrk::kNpdfMax, 0);
        long long npos = 0;
        {
          const int Np = run.N_pdf;
          const long long nr = static_cast<long long>(rays.size());
#pragma omp parallel
          {
            std::vector<long long> h(rrk::kNpdfMax, 0), hs(rrk::kNpdfMax, 0);
            long long np = 0;
#pragma omp for schedule(static) nowait
            for (long long r = 0; r < nr; ++r) {
              const double ntot = rrk::sum_ngam(rays[r]);
              const int b = rrk::pdf_bin(ntot, Np);
              h[b] += 1;
              if (rays[r].dr < 0 && rays[r].nside < radray::rayp::nside_max) hs[b] += 1;
              if (ntot > 0) np += 1;
            }
#pragma omp critical
            {
              for (int b = 0; b < rrk::kNpdfMax; ++b) hist[b] += h[b];
              for (int b = 0; b < rrk::kNpdfMax; ++b) spl[b] += hs[b];
              npos += np;
            }
          }
        }
        long long nsrc = 0;
        for (long long c = 0; c < ncell; ++c) nsrc += is_src[c] ? 1 : 0;
        run.thresholds(hist.data(), npos, spl.data(), nsrc);
        const rrk::RTP P = run.rtp();
        double t2 = now();
        T.hist += t2 - t1;
        // 1d. per-cell classify/emit/split/merge: count, scan, write
        std::vector<long long> cnt(ncell);
        long long ncst = 0, nspl = 0, nmrg = 0, ndel = 0;
#pragma omp parallel for schedule(dynamic, 256) reduction(+ : ncst, nspl, nmrg, ndel)
        for (long long c = 0; c < ncell; ++c) {
          const int i = static_cast<int>(c % N), j = static_cast<int>((c / N) % N),
                    k = static_cast<int>(c / (static_cast<long long>(N) * N));
          const bool src = is_src[c];
          const bool thick = src && rrk::source_is_thick(cells[c].tau);
          rrk::SourceSpec spec;
          if (src && !thick) spec = rrk::source_spec(i, j, k, cells[c].Luv, P);
          rrk::CellStats cs;
          cnt[c] = rrk::step1_count(rays.data() + start[c], static_cast<int>(start[c + 1] - start[c]),
                                    (src && !thick) ? &spec : nullptr, P, cs);
          ncst += cs.ncst;
          nspl += cs.nspl;
          nmrg += cs.nmrg;
          ndel += cs.ndel;
        }
        ostart[0] = 0;
        for (long long c = 0; c < ncell; ++c) ostart[c + 1] = ostart[c] + cnt[c];
        out.resize(ostart[ncell]);
        idx.resize(rays.size());
        key.resize(rays.size());
        std::fill(acc.begin(), acc.end(), 0);
#pragma omp parallel for schedule(dynamic, 256)
        for (long long c = 0; c < ncell; ++c) {
          const int i = static_cast<int>(c % N), j = static_cast<int>((c / N) % N),
                    k = static_cast<int>(c / (static_cast<long long>(N) * N));
          const bool src = is_src[c];
          const bool thick = src && rrk::source_is_thick(cells[c].tau);
          if (thick) {
            double v[rrk::kNdep];
            rrk::source_thick_deposit(cells[c].Luv, P, v);
            for (int m = 0; m < rrk::kNdep; ++m)
              atomic_add(&acc[c * rrk::kNdep + m], rrk::to_fx(v[m], rrk::dep_scale(m)));
          }
          if (cnt[c] == 0) continue;
          rrk::SourceSpec spec;
          if (src && !thick) spec = rrk::source_spec(i, j, k, cells[c].Luv, P);
          const int w = rrk::step1_write(rays.data() + start[c], static_cast<int>(start[c + 1] - start[c]),
                                         (src && !thick) ? &spec : nullptr, i, j, k, P, idx.data() + start[c],
                                         key.data() + start[c], out.data() + ostart[c]);
          if (w != cnt[c]) std::abort();
        }
        double t3 = now();
        T.step1 += t3 - t2;
        // 1d'. cast every ray once (deposits)
        {
          const long long nr = static_cast<long long>(out.size());
          auto tau = [&](int i, int j, int k) { return cells[cidx(i, j, k)].tau; };
          auto dep = [&](int i, int j, int k, const double* v) {
            long long* a = &acc[cidx(i, j, k) * rrk::kNdep];
            for (int m = 0; m < rrk::kNdep; ++m) atomic_add(a + m, rrk::to_fx(v[m], rrk::dep_scale(m)));
          };
#pragma omp parallel for schedule(static)
          for (long long r = 0; r < nr; ++r) rrk::cast_1(out[r], P, tau, dep);
          ray_substeps += nr;
        }
        double t4 = now();
        T.cast1 += t4 - t3;
        // 2. physics
        {
          const rrk::PhysCtx X = run.phys();
#pragma omp parallel for schedule(dynamic, 1024)
          for (long long c = 0; c < ncell; ++c) {
            CellData& d = cells[c];
            const long long* a = &acc[c * rrk::kNdep];
            for (int m = 0; m < 5; ++m) d.ngam[m] = static_cast<float>(rrk::from_fx(a[m], rrk::kScaleNgam));
            for (int m = 0; m < 5; ++m) d.nu[m] = static_cast<float>(rrk::from_fx(a[5 + m], rrk::kScaleNu));
            for (int m = 0; m < 3; ++m) d.dx[m] = static_cast<float>(rrk::from_fx(a[10 + m], rrk::kScaleNgam));
            rrk::physics_cell(d, X);
          }
        }
        double t5 = now();
        T.phys += t5 - t4;
        // 3. advance
        long long ndel3 = 0;
        {
          const long long nr = static_cast<long long>(out.size());
          home.resize(nr);
          keep.resize(nr);
          auto tf = [&](int i, int j, int k, const float*& tau, const float*& fgam) {
            const CellData& d = cells[cidx(i, j, k)];
            tau = d.tau;
            fgam = d.fgam;
          };
#pragma omp parallel for schedule(static) reduction(+ : ndel3)
          for (long long r = 0; r < nr; ++r) {
            keep[r] = rrk::cast_3(out[r], P, tf);
            if (keep[r]) {
              home[r] = static_cast<int>(cidx(rrk::cell_of(out[r].x[0], N), rrk::cell_of(out[r].x[1], N),
                                              rrk::cell_of(out[r].x[2], N)));
            } else {
              ndel3 += 1;
            }
          }
        }
        double t6 = now();
        T.cast3 += t6 - t5;
        // re-sort by home cell (stable counting sort, compacting deletions)
        {
          std::fill(start.begin(), start.end(), 0);
          const long long nr = static_cast<long long>(out.size());
          for (long long r = 0; r < nr; ++r)
            if (keep[r]) start[home[r] + 1] += 1;
          for (long long c = 0; c < ncell; ++c) start[c + 1] += start[c];
          rays.resize(start[ncell]);
          std::vector<long long> pos(start.begin(), start.end() - 1);
          for (long long r = 0; r < nr; ++r)
            if (keep[r]) rays[pos[home[r]]++] = out[r];
        }
        T.sort += now() - t6;
        st.ncst1 = ncst;
        st.nspl1 = nspl;
        st.nmrg1 = nmrg;
        st.ndel1 = ndel + ndel3;
        run.end_substep();
      }
      // post_rt: redshift rays, final gas/tau
      const double t7 = now();
      {
        const rrk::RedshiftTable RT = run.redshift_table();
        const double fscale = run.redshift_fscale(), ngam_min = run.units.ngam_min;
        const long long nr = static_cast<long long>(rays.size());
#pragma omp parallel for schedule(static)
        for (long long r = 0; r < nr; ++r) rrk::redshift_ray(rays[r], RT, fscale, ngam_min);
        const radray::Units u = run.units;
#pragma omp parallel for schedule(static)
        for (long long c = 0; c < ncell; ++c) rrk::post_cell(cells[c], u, run.sig_init);
        long long tot = 0;
#pragma omp parallel for schedule(static) reduction(+ : tot)
        for (long long r = 0; r < nr; ++r)
          for (int m = 0; m < 5; ++m) tot += rrk::to_fx(static_cast<double>(rays[r].ngam[m]), rrk::kScaleTot);
        st.ngam_tot_fx = tot;
        st.nray_active = nr;
      }
      T.post += now() - t7;
      std::printf(" step %3d  z = %.5f  nstep_rt = %d  rays = %lld  cast = %lld  merge = %lld  split = %lld\n",
                  run.clk.istep, run.clk.redshift, nsub, st.nray_active, st.ncst1, st.nmrg1, st.nspl1);
      if (output) run.write_step(odir, cells, st);
    }
    const double tt = now() - t_start;
    std::printf(" threads %d  wall %.3f s  ray-substeps %lld  ns/ray-substep %.1f\n", nthreads, tt,
                ray_substeps, 1e9 * tt * nthreads / (ray_substeps > 0 ? ray_substeps : 1));
    std::printf(" timers: cell %.3f hist %.3f step1 %.3f cast1 %.3f phys %.3f cast3 %.3f sort %.3f post %.3f\n",
                T.cell, T.hist, T.step1, T.cast1, T.phys, T.cast3, T.sort, T.post);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "rt_serial: error: %s\n", e.what());
    return 1;
  }
  return 0;
}
