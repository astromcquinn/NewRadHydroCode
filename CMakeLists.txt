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

# OpenMP (docs/parallel.md): hot loops of RT, hydro, PM, FoF and the
# summaries. Results are bit-identical for any thread count; the ray tracing
# itself runs in parallel only with deterministic_parallel = 1.
option(RADRAY_OPENMP "Build with OpenMP" ON)
set(RADRAY_HAVE_OPENMP OFF)
if(RADRAY_OPENMP)
  find_package(OpenMP)
  if(OpenMP_CXX_FOUND)
    set(RADRAY_HAVE_OPENMP ON)
    target_link_libraries(radray PUBLIC OpenMP::OpenMP_CXX)
    message(STATUS "OpenMP: enabled")
  else()
    message(STATUS "OpenMP: not found, serial build")
  endif()
endif()
if(NOT RADRAY_HAVE_OPENMP)
  target_compile_options(radray PUBLIC -Wno-unknown-pragmas)
endif()

# On-the-fly halo finder, source model, legacy halo catalogs
# (docs/halos_and_sources.md). OpenMP is optional; results are identical
# for any thread count.
add_library(radray_halos
  src/halos.cpp
  src/halo_sources.cpp
  src/halo_capi.cpp
  src/mass_function.cpp)
target_link_libraries(radray_halos PUBLIC radray)
add_executable(find_halos apps/find_halos.cpp)
target_link_libraries(find_halos PRIVATE radray_halos)

add_executable(rt_driver apps/rt_driver.cpp)
target_link_libraries(rt_driver PRIVATE radray)

# PM gravity (pm_tools.f90 / fft_tools.f90) needs FFTW3 single precision.
# Optional: without it only the RT part is built.
option(RADRAY_WITH_PM "Build the PM gravity and hydro ports (need FFTW3f)" ON)
set(RADRAY_HAVE_PM OFF)
if(RADRAY_WITH_PM)
  find_path(FFTW3F_INCLUDE_DIR fftw3.h)
  find_library(FFTW3F_LIBRARY fftw3f)
  if(FFTW3F_INCLUDE_DIR AND FFTW3F_LIBRARY)
    set(RADRAY_HAVE_PM ON)
    # One FFT layer (MKL DFTI semantics on FFTW3f) for the hydro (1-D) and
    # the PM (3-D), see include/radray/fft.hpp.
    add_library(radray_fft src/fft.cpp)
    target_include_directories(radray_fft PRIVATE ${FFTW3F_INCLUDE_DIR})
    target_link_libraries(radray_fft PUBLIC radray ${FFTW3F_LIBRARY})
    # Hydro solver (hydro_tools.f90 -DMF), summarize_gas and the gas
    # domains, see docs/port_notes_hydro.md.
    add_library(radray_hydro src/hydro.cpp src/domain_gas.cpp)
    target_link_libraries(radray_hydro PUBLIC radray_fft)
    add_executable(hydro_driver apps/hydro_driver.cpp)
    target_link_libraries(hydro_driver PRIVATE radray_hydro)
    # PM gravity (pm_tools.f90); uses the gas domains and summarize_gas.
    add_library(radray_pm src/pm.cpp src/checkpoint.cpp)
    target_link_libraries(radray_pm PUBLIC radray_hydro)
    add_executable(pm_driver apps/pm_driver.cpp)
    target_link_libraries(pm_driver PRIVATE radray_pm)
    # Full radhydro.f90 time loop (docs/full_run.md)
    add_executable(radhydro apps/radhydro.cpp)
    target_link_libraries(radhydro PRIVATE radray_pm radray_halos)
    message(STATUS "PM gravity and hydro: enabled (FFTW3f: ${FFTW3F_LIBRARY})")
  else()
    message(STATUS "PM gravity: disabled (FFTW3f not found)")
  endif()
endif()

# The ~5-minute tests (label nightly) are registered but disabled unless
# this is ON: cmake -DRADRAY_NIGHTLY_TESTS=ON build && ctest -L nightly
option(RADRAY_NIGHTLY_TESTS "Enable the long (nightly) CTests" OFF)

enable_testing()
add_subdirectory(tests)
