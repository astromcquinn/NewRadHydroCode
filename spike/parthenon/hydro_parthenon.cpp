// Parthenon prototype of the spike hydro step: the unsplit MUSCL-Hancock +
// HLLC update of spike/common/hydro_mh.hpp on a periodic n^3 mesh of
// meshblocks (2 ghosts), bit-identical to spike/serial/hydro_serial.cpp.
//
//   * conserved state: Parthenon variable "cons" (5 components,
//     Metadata::FillGhost), ghosts filled by Parthenon's boundary
//     communication tasks (task collection, one list per block partition);
//   * per block: Kokkos par_for for predict / face flux / update, scratch
//     Kokkos views for the predictor states Q and the face fluxes F;
//   * dt from a per-block Kokkos max reduction + MPI_Allreduce(MAX) (exact);
//   * our own fixed number of steps and dt (no Parthenon integrator/driver).
//
// usage: hydro_parthenon -i hydro.in [--n 64] [--steps 20] [--ic blast|sod]
//                        [--cfl 0.4] [--mb 16] [--out U.bin] [block/par=value ...]
// Output: raw float64 (5, n, n, n), component fastest, i fastest (rank 0).
#include "parth_common.hpp"

#include "hydro_mh.hpp"

using namespace parthenon;

namespace {

Packages_t process_packages(std::unique_ptr<ParameterInput>& pin) {
  Packages_t packages;
  auto pkg = std::make_shared<StateDescriptor>("spike_hydro");
  Metadata m({Metadata::Cell, Metadata::Independent, Metadata::FillGhost}, std::vector<int>{rrh::NV});
  pkg->AddField("cons", m);
  packages.Add(pkg);
  return packages;
}

// U(i,j,k,m) accessor over the Parthenon variable, local interior indices
struct UAcc {
  ParArrayND<Real, VariableState> v;
  int is, js, ks;
  KOKKOS_INLINE_FUNCTION double operator()(int i, int j, int k, int m) const {
    return v(m, k + ks, j + js, i + is);
  }
};

struct BlockScratch {
  Kokkos::View<double****, Kokkos::LayoutRight> Q;     // (n+2)^3 x NQ, offset 1
  Kokkos::View<double****, Kokkos::LayoutRight> F[3];  // (n+1)^3 x NV
};

}  // namespace

int main(int argc, char** argv) {
  int n = 64, steps = 20, ic = 0, mb = 16;
  double cfl = 0.4;
  std::string out;
  for (int a = 1; a < argc; ++a) {
    const std::string s = argv[a];
    auto next = [&]() -> const char* {
      if (a + 1 >= argc) { std::fprintf(stderr, "missing value for %s\n", s.c_str()); std::exit(2); }
      return argv[++a];
    };
    if (s == "--n") n = std::atoi(next());
    else if (s == "--steps") steps = std::atoi(next());
    else if (s == "--ic") ic = std::string(next()) == "sod" ? 1 : 0;
    else if (s == "--cfl") cfl = std::atof(next());
    else if (s == "--mb") mb = std::atoi(next());
    else if (s == "--out") out = next();
    else if (s == "-i") ++a;  // Parthenon input deck
    else if (s.find('=') != std::string::npos && s.find('/') < s.find('=')) continue;  // block/par=value
    else { std::fprintf(stderr, "unknown argument %s\n", s.c_str()); return 2; }
  }

  sp::ParthArgs pargs = sp::parthenon_args(argc, argv);
  ParthenonManager pman;
  auto st = pman.ParthenonInitEnv(pargs.argc(), pargs.argv());
  if (st != ParthenonStatus::ok) {
    pman.ParthenonFinalize();
    return st == ParthenonStatus::complete ? 0 : 1;
  }
  {
    auto& pin = pman.pinput;
    const char* dn[3] = {"1", "2", "3"};
    for (int d = 0; d < 3; ++d) {
      pin->SetInteger("parthenon/mesh", std::string("nx") + dn[d], n);
      pin->SetReal("parthenon/mesh", std::string("x") + dn[d] + "min", 0.0);
      pin->SetReal("parthenon/mesh", std::string("x") + dn[d] + "max", 1.0);
      pin->SetString("parthenon/mesh", std::string("ix") + dn[d] + "_bc", "periodic");
      pin->SetString("parthenon/mesh", std::string("ox") + dn[d] + "_bc", "periodic");
      pin->SetInteger("parthenon/meshblock", std::string("nx") + dn[d], mb);
    }
    pin->SetInteger("parthenon/mesh", "nghost", rrh::NG);
    pin->SetString("parthenon/mesh", "refinement", "none");
  }
  pman.app_input->ProcessPackages = process_packages;
  pman.ParthenonInitPackagesAndMesh();
  {
    Mesh* pmesh = pman.pmesh.get();
    MPI_Comm comm;
    MPI_Comm_dup(MPI_COMM_WORLD, &comm);
    int rank, nranks;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nranks);
    constexpr int NV = rrh::NV, NQ = rrh::NQ;
    auto& blocks = pmesh->block_list;
    const int nb = static_cast<int>(blocks.size());
    std::vector<sp::BlockGeom> geo(nb);
    std::vector<BlockScratch> scr(nb);
    for (int b = 0; b < nb; ++b) {
      MeshBlock* pmb = blocks[b].get();
      geo[b] = sp::block_geom(pmb, n, 0.0, 1.0);
      const auto& g = geo[b];
      scr[b].Q = Kokkos::View<double****, Kokkos::LayoutRight>("Q", g.n[2] + 2, g.n[1] + 2, g.n[0] + 2, NQ);
      for (int d = 0; d < 3; ++d)
        scr[b].F[d] = Kokkos::View<double****, Kokkos::LayoutRight>("F", g.n[2] + 1, g.n[1] + 1, g.n[0] + 1, NV);
      // initial condition on the interior
      auto u = pmb->meshblock_data.Get()->Get("cons").data;
      const int is = g.is, js = g.js, ks = g.ks, ox = g.off[0], oy = g.off[1], oz = g.off[2];
      const int icl = ic, nn = n;
      pmb->par_for(
          "hydro_ic", 0, g.n[2] - 1, 0, g.n[1] - 1, 0, g.n[0] - 1,
          KOKKOS_LAMBDA(const int k, const int j, const int i) {
            double U[NV];
            rrh::initial_state(icl, ox + i, oy + j, oz + k, nn, U);
            for (int m = 0; m < NV; ++m) u(m, k + ks, j + js, i + is) = U[m];
          });
    }
    const double dx = 1.0 / n;
    double t = 0, t_halo = 0;
    MPI_Barrier(comm);
    const double t0 = sp::now();
    for (int s = 0; s < steps; ++s) {
      const double th = sp::now();
      sp::exchange_ghosts(pmesh, "cons", {});
      t_halo += sp::now() - th;
      double smax = 0;
      for (int b = 0; b < nb; ++b) {
        MeshBlock* pmb = blocks[b].get();
        const auto& g = geo[b];
        auto u = pmb->meshblock_data.Get()->Get("cons").data;
        const int is = g.is, js = g.js, ks = g.ks;
        double bmax = 0;
        pmb->par_reduce(
            "hydro_smax", 0, g.n[2] - 1, 0, g.n[1] - 1, 0, g.n[0] - 1,
            KOKKOS_LAMBDA(const int k, const int j, const int i, double& lmax) {
              double U[NV];
              for (int m = 0; m < NV; ++m) U[m] = u(m, k + ks, j + js, i + is);
              lmax = rrf::max(lmax, rrh::max_speed(U));
            },
            Kokkos::Max<double>(bmax));
        smax = rrf::max(smax, bmax);
      }
      MPI_Allreduce(MPI_IN_PLACE, &smax, 1, MPI_DOUBLE, MPI_MAX, comm);
      const double dt = cfl * dx / smax;
      const double dtdx = dt / dx;
      for (int b = 0; b < nb; ++b) {
        MeshBlock* pmb = blocks[b].get();
        const auto& g = geo[b];
        const UAcc U{pmb->meshblock_data.Get()->Get("cons").data, g.is, g.js, g.ks};
        auto Q = scr[b].Q;
        pmb->par_for(
            "hydro_predict", -1, g.n[2], -1, g.n[1], -1, g.n[0],
            KOKKOS_LAMBDA(const int k, const int j, const int i) {
              rrh::predict_cell(U, i, j, k, dtdx, &Q(k + 1, j + 1, i + 1, 0));
            });
        for (int d = 0; d < 3; ++d) {
          const int ei = d == 0, ej = d == 1, ek = d == 2;
          auto F = scr[b].F[d];
          pmb->par_for(
              "hydro_flux", 0, g.n[2] - 1 + ek, 0, g.n[1] - 1 + ej, 0, g.n[0] - 1 + ei,
              KOKKOS_LAMBDA(const int k, const int j, const int i) {
                rrh::face_flux(&Q(k + 1 - ek, j + 1 - ej, i + 1 - ei, 0), &Q(k + 1, j + 1, i + 1, 0), d,
                               &F(k, j, i, 0));
              });
        }
        auto u = U.v;
        const int is = g.is, js = g.js, ks = g.ks;
        auto Fx = scr[b].F[0], Fy = scr[b].F[1], Fz = scr[b].F[2];
        pmb->par_for(
            "hydro_update", 0, g.n[2] - 1, 0, g.n[1] - 1, 0, g.n[0] - 1,
            KOKKOS_LAMBDA(const int k, const int j, const int i) {
              double Uc[NV];
              for (int m = 0; m < NV; ++m) Uc[m] = u(m, k + ks, j + js, i + is);
              rrh::update_cell(Uc, &Fx(k, j, i, 0), &Fx(k, j, i + 1, 0), &Fy(k, j, i, 0), &Fy(k, j + 1, i, 0),
                               &Fz(k, j, i, 0), &Fz(k + 1, j, i, 0), dtdx);
              for (int m = 0; m < NV; ++m) u(m, k + ks, j + js, i + is) = Uc[m];
            });
      }
      Kokkos::fence();
      t += dt;
    }
    MPI_Barrier(comm);
    const double wall = sp::now() - t0;

    // gather the interior to rank 0, (component fastest, i fastest)
    std::vector<long long> gidx;
    std::vector<char> rec;
    for (int b = 0; b < nb; ++b) {
      MeshBlock* pmb = blocks[b].get();
      const auto& g = geo[b];
      auto u = pmb->meshblock_data.Get()->Get("cons").data.GetHostMirrorAndCopy();
      for (int k = 0; k < g.n[2]; ++k)
        for (int j = 0; j < g.n[1]; ++j)
          for (int i = 0; i < g.n[0]; ++i) {
            gidx.push_back((g.off[0] + i) + static_cast<long long>(n) * ((g.off[1] + j) + static_cast<long long>(n) * (g.off[2] + k)));
            double U[NV];
            for (int m = 0; m < NV; ++m) U[m] = u(m, k + g.ks, j + g.js, i + g.is);
            const char* p = reinterpret_cast<const char*>(U);
            rec.insert(rec.end(), p, p + sizeof U);
          }
    }
    const long long ncell = static_cast<long long>(n) * n * n;
    std::vector<char> all = sp::gather_cells(gidx, rec, NV * sizeof(double), ncell, comm);
    if (rank == 0) {
      const double* o = reinterpret_cast<const double*>(all.data());
      double mass = 0, etot = 0;  // same summation order as hydro_serial
      for (long long c = 0; c < ncell; ++c) {
        mass += o[c * NV];
        etot += o[c * NV + 4];
      }
      const int nth = sp::num_threads();
      std::printf(" n %d steps %d t %.6e mass %.17g energy %.17g\n", n, steps, t, mass * dx * dx * dx,
                  etot * dx * dx * dx);
      std::printf(" ranks %d threads %d blocks %d (mb %d) wall %.3f s  halo %.3f s  ns/cell-update %.2f (core-ns %.2f)\n",
                  nranks, nth, static_cast<int>(pmesh->nbtotal), mb, wall, t_halo,
                  1e9 * wall / (double(n) * n * n * steps),
                  1e9 * wall * nth * nranks / (double(n) * n * n * steps));
      if (!out.empty()) {
        std::FILE* fp = std::fopen(out.c_str(), "wb");
        if (!fp) sp::die("cannot open output");
        std::fwrite(o, sizeof(double), ncell * NV, fp);
        std::fclose(fp);
      }
    }
    MPI_Comm_free(&comm);
  }
  pman.ParthenonFinalize();
  return 0;
}
