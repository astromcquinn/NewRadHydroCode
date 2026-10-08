// Helper for healpix_vs_healpy.py.
//   healpix_tool vec2pix <in.f64> <out.i32>   in: records (nside, x, y, z) as f64
//   healpix_tool pix2vec <in.i32> <out.f64>   in: records (nside, ipix) as i32
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <vector>

#include "radray/healpix.hpp"

template <class T>
static std::vector<T> slurp(const char* path) {
  std::ifstream in(path, std::ios::binary | std::ios::ate);
  const auto n = static_cast<std::size_t>(in.tellg());
  in.seekg(0);
  std::vector<T> v(n / sizeof(T));
  in.read(reinterpret_cast<char*>(v.data()), v.size() * sizeof(T));
  return v;
}

int main(int argc, char** argv) {
  if (argc != 4) return 2;
  if (std::strcmp(argv[1], "vec2pix") == 0) {
    const auto in = slurp<double>(argv[2]);
    std::vector<std::int32_t> out(in.size() / 4);
    for (std::size_t n = 0; n < out.size(); ++n) {
      const double v[3] = {in[4 * n + 1], in[4 * n + 2], in[4 * n + 3]};
      out[n] = radray::healpix::vec2pix_nest(static_cast<int>(in[4 * n]), v);
    }
    std::ofstream(argv[3], std::ios::binary)
        .write(reinterpret_cast<const char*>(out.data()), out.size() * 4);
  } else if (std::strcmp(argv[1], "pix2vec") == 0) {
    const auto in = slurp<std::int32_t>(argv[2]);
    std::vector<double> out(3 * (in.size() / 2));
    for (std::size_t n = 0; n < in.size() / 2; ++n)
      radray::healpix::pix2vec_nest(in[2 * n], in[2 * n + 1], &out[3 * n]);
    std::ofstream(argv[3], std::ios::binary)
        .write(reinterpret_cast<const char*>(out.data()), out.size() * 8);
  } else {
    return 2;
  }
  return 0;
}
