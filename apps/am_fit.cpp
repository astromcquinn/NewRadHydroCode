// am_fit: abundance-match halos to observed UV luminosity functions and
// fit the parametric L_UV(M_halo, z) used by sfr_model = abundance
// (docs/abundance_matching.md).
//
//   am_fit --lf data/uv_lf/lf_params.txt --out am_out [options]
//
// Options:
//   --hmf st|tinker        analytic halo mass function (default tinker)
//   --tinker-zmax Z        freeze the Tinker z-evolution above Z
//   --st-params A,a,p      Sheth-Tormen coefficients (default 0.3222,0.707,0.3;
//                          production: the fit to the FoF halos,
//                          data/hmf/st_fof_*.txt)
//   --sim-hmf FILE         match to the simulated FoF halos: the analytic HMF
//                          (--hmf; use st with fit_fof_hmf.py) times the
//                          fitted correction in FILE (data/hmf/*.txt); the
//                          matched masses start at nmin particles
//   --sim-hmf FILE         match to the simulated FoF halos: the analytic HMF
//                          (--hmf, use st with fit_fof_hmf.py) times the
//                          fitted correction of FILE (data/hmf/*.txt); the
//                          matched masses start at nmin particles
//   --dust 0|1             Meurer IRX-beta dust correction (default 0)
//   --beta0 B --dbeta D    beta_UV(M) = B + D (M + 19.5)
//   --extend-faint M       the shape fit also uses matched points down to
//                          M_UV = M (default -12; 0 = observed range only)
//   --logM-min X --logM-max X   halo mass range of the matched points
//   --free l0,l1,a0,...    free parameters (default: all eight)
//   --set am_b=0.3         starting/fixed value of a parameter (repeatable)
//   --shape-zmax Z         entries with z <= Z fix the shape of L_UV(M, z)
//                          (all parameters fitted); entries above Z only set
//                          the normalisation (one offset per redshift, the
//                          shape held at its z_hi form). Default 9.5.
//   --zlo Z --zhi Z        clamp of the shape's redshift dependence
//                          (default: the range of the shape entries)
//   --only-confirmed       ignore entries marked unconfirmed
//
// Writes <out>/am_points.txt, <out>/am_params.txt (params.txt lines) and
// <out>/lf_check.txt (input vs implied LF), and prints a summary.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "radray/abundance_matching.hpp"

using namespace radray;

namespace {

void usage() {
  std::cerr << "usage: am_fit --lf FILE --out DIR [--hmf st|tinker] [--dust 0|1] "
               "[--beta0 B] [--dbeta D] [--extend-faint M] [--logM-min X] [--logM-max X] "
               "[--free l0,l1,...] [--set am_x=v] [--zlo Z] [--zhi Z] [--tinker-zmax Z] [--st-params A,a,p] [--sim-hmf FILE] "
               "[--only-confirmed] [--cosmo key=value]\n"
               "  --cosmo omegam=..|omegab=..|hubble0=..|sigma8=..|ns=..: the runtime cosmology of\n"
               "  the simulation (params.txt keys, docs/production_run.md section 10; repeatable;\n"
               "  default the original's)\n";
}

int par_index(const std::string& name) {
  for (int i = 0; i < LuvModel::kNpar; ++i) {
    const std::string n = LuvModel::name(i);
    if (name == n || "am_" + name == n) return i;
  }
  throw std::runtime_error("unknown parameter " + name);
}

}  // namespace

int main(int argc, char** argv) {
  try {
    std::string lf_path, out_dir = "am_out";
    AMOptions opt;
    opt.extend_faint_to = -12;
    LuvModel start;
    std::array<bool, LuvModel::kNpar> free{};
    free.fill(true);
    bool only_confirmed = false;
    double zlo = -1, zhi = -1, shape_zmax = 9.5;
    Cosmology cosmo;  // --cosmo key=value (the original's by default)
    for (int i = 1; i < argc; ++i) {
      const std::string a = argv[i];
      auto next = [&]() -> std::string {
        if (i + 1 >= argc) throw std::runtime_error("missing value for " + a);
        return argv[++i];
      };
      if (a == "--lf") lf_path = next();
      else if (a == "--out") out_dir = next();
      else if (a == "--hmf") {
        const std::string v = next();
        if (v == "st") opt.hmf = MassFunction::Kind::ShethTormen;
        else if (v == "tinker") opt.hmf = MassFunction::Kind::Tinker08;
        else if (v == "ps") opt.hmf = MassFunction::Kind::PressSchechter;
        else throw std::runtime_error("--hmf st|tinker|ps");
      } else if (a == "--tinker-zmax") opt.tinker_zmax_evol = std::stod(next());
      else if (a == "--st-params") {  // A,a,p (data/hmf/st_fof_*.txt)
        const std::string v = next();
        if (std::sscanf(v.c_str(), "%lf,%lf,%lf", &opt.st.A, &opt.st.a, &opt.st.p) != 3)
          throw std::runtime_error("--st-params A,a,p");
      }
      else if (a == "--sim-hmf") opt.sim_hmf = SimHMF::read(next());
      else if (a == "--dust") opt.dust.enabled = std::stoi(next()) != 0;
      else if (a == "--beta0") opt.dust.beta0 = std::stod(next());
      else if (a == "--dbeta") opt.dust.dbeta_dM = std::stod(next());
      else if (a == "--extend-faint") opt.extend_faint_to = std::stod(next());
      else if (a == "--logM-min") opt.logM_min = std::stod(next());
      else if (a == "--logM-max") opt.logM_max = std::stod(next());
      else if (a == "--zlo") zlo = std::stod(next());
      else if (a == "--zhi") zhi = std::stod(next());
      else if (a == "--shape-zmax") shape_zmax = std::stod(next());
      else if (a == "--only-confirmed") only_confirmed = true;
      else if (a == "--free") {
        free.fill(false);
        std::stringstream ss(next());
        std::string t;
        while (std::getline(ss, t, ',')) free[par_index(t)] = true;
      } else if (a == "--set") {
        const std::string kv = next();
        const auto eq = kv.find('=');
        if (eq == std::string::npos) throw std::runtime_error("--set name=value");
        start.p[par_index(kv.substr(0, eq))] = std::stod(kv.substr(eq + 1));
      } else if (a == "--cosmo") {
        const std::string kv = next();
        const auto eq = kv.find('=');
        if (eq == std::string::npos || !set_cosmology_key(cosmo, kv.substr(0, eq), kv.substr(eq + 1)))
          throw std::runtime_error("--cosmo omegam|omegab|hubble0|sigma8|ns=value");
      } else if (a == "-h" || a == "--help") {
        usage();
        return 0;
      } else {
        usage();
        throw std::runtime_error("unknown option " + a);
      }
    }
    cosmo.derive();
    opt.cosmo = MassFunctionCosmology::from(cosmo);
    if (opt.sim_hmf.enabled)
      opt.logM_min = std::max(opt.logM_min,
                              std::log10(opt.sim_hmf.nmin * opt.sim_hmf.mpart_msun) + 1e-9);
    if (lf_path.empty()) {
      usage();
      return 2;
    }

    std::vector<LFEntry> lf;
    int n_unconf = 0;
    for (const auto& e : read_lf_table(lf_path)) {
      if (!e.confirmed && e.weight > 0) ++n_unconf;
      if (only_confirmed && !e.confirmed) continue;
      if (!(e.weight > 0)) continue;  // kept in the table for reference only
      lf.push_back(e);
    }
    if (lf.empty()) throw std::runtime_error("no LF entries to fit");
    if (n_unconf > 0)
      std::cout << "WARNING: " << n_unconf
                << " LF entries are marked unconfirmed (parameters not yet checked "
                   "against the papers)\n";

    double zmin = 1e30, zmax = -1e30;
    std::vector<AMPoint> pts;
    for (const auto& e : lf) {
      const auto p = abundance_match(e, opt);
      pts.insert(pts.end(), p.begin(), p.end());
      if (e.z <= shape_zmax && e.weight > 0) {
        zmin = std::min(zmin, e.z);
        zmax = std::max(zmax, e.z);
      }
    }
    if (zmax < zmin) throw std::runtime_error("no LF entries at z <= --shape-zmax");
    start.z_lo = zlo > 0 ? zlo : zmin;
    start.z_hi = zhi > 0 ? zhi : zmax;

    FitResult fr = fit_luv_model(pts, start, free);
    fr.model = fit_highz_norm(pts, fr.model);

    std::filesystem::create_directories(out_dir);
    {
      std::ofstream o(out_dir + "/am_points.txt");
      o << "# z logM[Msun] n_halo[Mpc^-3] M_obs M_int logL_AM logL_fit in_range in_shape\n";
      for (const auto& p : pts) {
        char buf[256];
        std::snprintf(buf, sizeof buf, "%.3f %.3f %.6e %.4f %.4f %.5f %.5f %d %d\n", p.z,
                      p.logM, p.n_halo, p.M_obs, p.M_int, p.logL, fr.model.log_luv(p.logM, p.z),
                      p.in_range ? 1 : 0, p.in_shape ? 1 : 0);
        o << buf;
      }
    }
    {
      std::ofstream o(out_dir + "/am_params.txt");
      o << "# sfr_model = abundance parameters from am_fit (" << lf_path << ")\n"
        << "# hmf = "
        << (opt.hmf == MassFunction::Kind::Tinker08 ? "tinker" : opt.hmf == MassFunction::Kind::ShethTormen ? "st" : "ps")
        << (opt.hmf == MassFunction::Kind::ShethTormen
                ? " (A " + std::to_string(opt.st.A) + ", a " + std::to_string(opt.st.a) + ", p " +
                      std::to_string(opt.st.p) + ")"
                : std::string())
        << (opt.sim_hmf.enabled ? " x simulated-FoF correction" : "")
        << ", dust = " << (opt.dust.enabled ? 1 : 0) << ", rms = " << fr.rms_dex
        << " dex over " << fr.npoints << " points"
        << (n_unconf > 0 && !only_confirmed ? ", INCLUDES UNCONFIRMED LF ENTRIES" : "")
        << "\nsfr_model = abundance\n"
        << fr.model.to_params();
    }
    {
      std::ofstream o(out_dir + "/lf_check.txt");
      o << "# z ref M_UV phi_input phi_model [Mpc^-3 mag^-1] in_range\n";
      for (const auto& e : lf) {
        std::vector<double> M;
        for (double m = -23.5; m <= -12.0 + 1e-9; m += 0.25) M.push_back(m);
        const auto phm = model_lf(fr.model, e.z, M, opt);
        for (std::size_t k = 0; k < M.size(); ++k) {
          char buf[256];
          std::snprintf(buf, sizeof buf, "%.3f %s %.2f %.6e %.6e %d\n", e.z, e.ref.c_str(), M[k],
                        e.phi(M[k]), phm[k], (M[k] >= e.M_bright && M[k] <= e.M_faint) ? 1 : 0);
          o << buf;
        }
      }
    }

    std::cout << "entries: " << lf.size() << ", fitted points: " << fr.npoints
              << ", rms = " << fr.rms_dex << " dex, iterations " << fr.iterations << "\n";
    std::cout << fr.model.to_params();
    if (fr.model.norm_z.empty())
      std::cout << "(no LF entries above z_hi = " << fr.model.z_hi
                << ": L(M) is held at its z_hi form at all higher z)\n";
    // 20 particles at 128^3 in 10 Mpc/h: 7.1e8 Msun/h = 1.0e9 Msun (logM = 9.0)
    std::cout << "M_UV of the fit (logM in Msun; 128^3 in 10 Mpc/h resolves logM >= 9.0):\n";
    for (double z : {6.0, 8.0, 10.0, 12.0, 14.0, 16.0}) {
      std::cout << "  z = " << z << ":";
      for (double lm : {8.5, 9.0, 10.0, 11.0}) {
        char buf[64];
        std::snprintf(buf, sizeof buf, "  logM=%.1f -> M_UV=%.2f", lm,
                      muv_from_luv(fr.model.luv(std::pow(10.0, lm), z)));
        std::cout << buf;
      }
      std::cout << "\n";
    }
    return 0;
  } catch (const std::exception& ex) {
    std::cerr << "am_fit: " << ex.what() << "\n";
    return 1;
  }
}
