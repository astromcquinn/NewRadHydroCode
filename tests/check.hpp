// Minimal test helpers (no external framework).
#pragma once

#include <cmath>
#include <cstdio>
#include <cstdlib>

namespace radray_test {
inline int& failures() {
  static int n = 0;
  return n;
}
}  // namespace radray_test

#define CHECK(cond)                                                          \
  do {                                                                       \
    if (!(cond)) {                                                           \
      std::fprintf(stderr, "%s:%d: CHECK failed: %s\n", __FILE__, __LINE__, \
                   #cond);                                                   \
      ++radray_test::failures();                                             \
    }                                                                        \
  } while (0)

#define CHECK_CLOSE(a, b, rtol)                                                \
  do {                                                                         \
    const double _a = (a), _b = (b);                                           \
    const double _d = std::fabs(_a - _b);                                      \
    if (!(_d <= (rtol) * std::fmax(std::fabs(_a), std::fabs(_b)) || _d == 0)) { \
      std::fprintf(stderr, "%s:%d: CHECK_CLOSE failed: %s=%.17g %s=%.17g\n",  \
                   __FILE__, __LINE__, #a, _a, #b, _b);                        \
      ++radray_test::failures();                                               \
    }                                                                          \
  } while (0)

#define TEST_MAIN_END()                                         \
  do {                                                          \
    if (radray_test::failures()) {                              \
      std::fprintf(stderr, "%d failure(s)\n", radray_test::failures()); \
      return 1;                                                 \
    }                                                           \
    std::printf("OK\n");                                        \
    return 0;                                                   \
  } while (0)
