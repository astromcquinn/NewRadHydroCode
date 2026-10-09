// Compile-only check that the spike's framework-independent kernels
// (spike/common) build as CUDA device code. Each kernel instantiates one
// family of kernels; nothing is launched.
//   nvcc -std=c++17 --expt-relaxed-constexpr --fmad=false -arch=sm_80 \
//        -I include -I spike/common -c spike/gpu/kernels_compile_test.cu
#include "hydro_mh.hpp"
#include "ray_kernels.hpp"
#include "rt_cell.hpp"

using namespace rrk;

struct TauView {
  const float* tau;
  int n;
  __device__ const float* operator()(int i, int j, int k) const {
    return tau + 5 * (i + n * (j + n * k));
  }
};
struct DepView {
  unsigned long long* acc;
  int n;
  __device__ void operator()(int i, int j, int k, const double* v) const {
    unsigned long long* a = acc + kNdep * (i + n * (j + n * k));
    for (int m = 0; m < kNdep; ++m)
      atomicAdd(a + m, static_cast<unsigned long long>(to_fx(v[m], dep_scale(m))));
  }
};
struct TFView {
  const float* tau;
  const float* fgam;
  int n;
  __device__ void operator()(int i, int j, int k, const float*& t, const float*& f) const {
    const int c = 5 * (i + n * (j + n * k));
    t = tau + c;
    f = fgam + c;
  }
};

__global__ void k_cast(RayP* rays, int nr, RTP P, TauView tv, DepView dv, TFView tf, int* keep) {
  const int r = blockIdx.x * blockDim.x + threadIdx.x;
  if (r >= nr) return;
  cast_1(rays[r], P, tv, dv);
  keep[r] = cast_3(rays[r], P, tf);
}

__global__ void k_step1(const RayP* in, const long long* start, int ncell, RTP P, int* idx,
                        MergeKey* key, RayP* out, const long long* ostart, long long* cnt) {
  const int c = blockIdx.x * blockDim.x + threadIdx.x;
  if (c >= ncell) return;
  const int n = static_cast<int>(start[c + 1] - start[c]);
  CellStats st;
  cnt[c] = step1_count(in + start[c], n, false, false, 0, 0, 0, 0.f, P, st);
  step1_write(in + start[c], n, false, false, 0, 0, 0, 0.f, P, idx + start[c], key + start[c],
              out + ostart[c]);
}

__global__ void k_cell(CellData* cells, int ncell, PhysCtx X) {
  const int c = blockIdx.x * blockDim.x + threadIdx.x;
  if (c >= ncell) return;
  physics_cell(cells[c], X);
}

__global__ void k_redshift(RayP* rays, int nr, const RedshiftTable* T) {
  const int r = blockIdx.x * blockDim.x + threadIdx.x;
  if (r < nr) redshift_ray(rays[r], *T);
}

struct UView {
  const double* u;
  int m;
  __device__ double operator()(int i, int j, int k, int c) const {
    return u[5 * ((i + 2) + m * ((j + 2) + m * (k + 2))) + c];
  }
};
__global__ void k_hydro(UView U, double* Q, int n, double dtdx) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) return;
  rrh::predict_cell(U, i, 0, 0, dtdx, Q + rrh::NQ * i);
  double F[rrh::NV];
  rrh::face_flux(Q, Q + rrh::NQ, 0, F);
}
