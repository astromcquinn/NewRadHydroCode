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

add_executable(rt_driver apps/rt_driver.cpp)
target_link_libraries(rt_driver PRIVATE radray)

# Hydro solver (hydro_tools.f90 -DMF) and gas domains: needs FFTW3 single
# precision (replaces MKL DFTI, see docs/port_notes_hydro.md).
find_path(FFTW3_INCLUDE_DIR fftw3.h REQUIRED)
find_library(FFTW3F_LIBRARY fftw3f REQUIRED)
add_library(radray_hydro src/hydro.cpp src/domain_gas.cpp)
target_include_directories(radray_hydro PUBLIC ${FFTW3_INCLUDE_DIR})
target_link_libraries(radray_hydro PUBLIC radray ${FFTW3F_LIBRARY})
add_executable(hydro_driver apps/hydro_driver.cpp)
target_link_libraries(hydro_driver PRIVATE radray_hydro)

enable_testing()
add_subdirectory(tests)
