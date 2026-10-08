cmake_minimum_required(VERSION 3.16)
project(radray CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

option(RADRAY_SANITIZE "Build with -fsanitize=address,undefined" OFF)

# Bit-for-bit agreement with the gfortran reference requires: no FMA
# contraction, no fast-math, IEEE SSE arithmetic (default on x86-64).
add_compile_options(-Wall -Wextra -ffp-contract=off -fno-fast-math)
if(RADRAY_SANITIZE)
  add_compile_options(-fsanitize=address,undefined -fno-omit-frame-pointer -g)
  add_link_options(-fsanitize=address,undefined)
endif()

add_library(radray
  src/healpix.cpp
  src/cosmology.cpp
  src/rt_rates.cpp
  src/ray_tracing.cpp
  src/rt_tools.cpp
  src/io.cpp)
target_include_directories(radray PUBLIC include)

# On-the-fly halo finder, source model, legacy halo catalogs
# (docs/halos_and_sources.md). OpenMP is optional; results are identical
# for any thread count.
option(RADRAY_HALOS_OPENMP "Build the halo finder with OpenMP" OFF)
add_library(radray_halos
  src/halos.cpp
  src/halo_sources.cpp
  src/mass_function.cpp)
target_link_libraries(radray_halos PUBLIC radray)
if(RADRAY_HALOS_OPENMP)
  find_package(OpenMP REQUIRED)
  target_link_libraries(radray_halos PUBLIC OpenMP::OpenMP_CXX)
else()
  target_compile_options(radray_halos PRIVATE -Wno-unknown-pragmas)
endif()

add_executable(rt_driver apps/rt_driver.cpp)
target_link_libraries(rt_driver PRIVATE radray)

enable_testing()
add_subdirectory(tests)
