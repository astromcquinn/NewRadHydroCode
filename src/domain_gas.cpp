// Port of domain_tools.f90: init_domain_gas, set_domain_gas, end_domain_gas.
//
// The gas grid is split into n1^3 cubic domains (n1 = 2 for Ncpu = 1) with a
// 26-neighbour table; set_domain_gas hands out a domain none of whose
// neighbours is in progress (status 1), end_domain_gas marks it done (2) and
// resets all statuses once every domain is done. In the original these
// schedule the OpenMP threads of data_tools.f90 and pm_tools.f90 (gas
// density assignment / force interpolation); the hydro sweeps themselves
// use plain k-loops.
#include <cmath>
#include <ostream>
#include <stdexcept>

#include "radray/fortran.hpp"
#include "radray/hydro.hpp"

namespace radray {

void init_domain_gas(const Config& cfg, Domains& dom, std::ostream* out) {
  const int Ngas = cfg.Ngas;
  // n = nint(Ncpu**(1./3)) : integer ** real(4)
  const int n0 = f::nint(static_cast<double>(
      std::pow(static_cast<float>(cfg.Ncpu), 1.0f / 3.0f)));
  int n1 = 0;
  for (int i = n0; i <= Ngas; ++i) {
    if (Ngas % (2 * i) == 0) {
      n1 = 2 * i;
      break;
    }
  }
  // The original leaves n1 undefined in this case.
  if (n1 == 0) throw std::runtime_error("init_domain_gas: no valid subdivision (Ngas odd?)");
  const int n2 = n1, n3 = n1;
  dom.Nx_dom = Ngas / n1;
  dom.Ndomain = n1 * n2 * n3;
  if (out)
    *out << " " << cfg.Ncpu << " " << dom.Ndomain << " " << n1 << " " << dom.Nx_dom << " "
         << Ngas << "\n";
  dom.box.assign(dom.Ndomain, {});
  dom.domain.assign(dom.Ndomain, {});
  for (int n = 1; n <= dom.Ndomain; ++n) {
    const int i = 1 + (n - 1) % n1;
    const int j = 1 + ((n - 1) / n1) % n2;
    const int k = 1 + ((n - 1) / n1 / n2);
    auto& b = dom.box[n - 1];
    b[0] = 1 + (i - 1) * dom.Nx_dom;
    b[1] = i * dom.Nx_dom;
    b[2] = 1 + (j - 1) * dom.Nx_dom;
    b[3] = j * dom.Nx_dom;
    b[4] = 1 + (k - 1) * dom.Nx_dom;
    b[5] = k * dom.Nx_dom;
    auto& d = dom.domain[n - 1];
    int m = 1;
    d[0] = 0;
    for (int cc = -1; cc <= 1; ++cc) {
      const int k1 = 1 + (k + cc - 1 + n3) % n3;
      for (int bb = -1; bb <= 1; ++bb) {
        const int j1 = 1 + (j + bb - 1 + n2) % n2;
        for (int aa = -1; aa <= 1; ++aa) {
          const int i1 = 1 + (i + aa - 1 + n1) % n1;
          const int l = i1 + (j1 - 1) * n1 + (k1 - 1) * n1 * n2;
          if (l != n) {
            m = m + 1;
            d[m - 1] = l;
          }
        }
      }
    }
  }
}

int set_domain_gas(Domains& dom) {
  for (int n = 1; n <= dom.Ndomain; ++n) {
    auto& d = dom.domain[n - 1];
    if (d[0] == 0) {
      bool eval = true;
      for (int m = 2; m <= 27; ++m) {
        if (dom.domain[d[m - 1] - 1][0] == 1) {
          eval = false;
          break;
        }
      }
      if (eval) {
        d[0] = 1;
        return n - 1;
      }
    }
  }
  // The original spins in `do while (loop)` until another thread releases a
  // domain; with one thread that would never terminate.
  throw std::runtime_error("set_domain_gas: no eligible domain");
}

void end_domain_gas(Domains& dom, int idom) {
  dom.domain[idom][0] = 2;
  long sum = 0;
  for (auto& d : dom.domain) sum += d[0];
  if (sum == 2L * dom.Ndomain)
    for (auto& d : dom.domain) d[0] = 0;
}

}  // namespace radray
