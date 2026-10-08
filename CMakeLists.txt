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

# PM gravity (pm_tools.f90 / fft_tools.f90) needs FFTW3 single precision.
# Optional: without it only the RT part is built.
option(RADRAY_WITH_PM "Build the PM gravity port (needs FFTW3f)" ON)
set(RADRAY_HAVE_PM OFF)
if(RADRAY_WITH_PM)
  find_path(FFTW3F_INCLUDE_DIR fftw3.h)
  find_library(FFTW3F_LIBRARY fftw3f)
  if(FFTW3F_INCLUDE_DIR AND FFTW3F_LIBRARY)
    set(RADRAY_HAVE_PM ON)
    add_library(radray_pm src/fft.cpp src/pm.cpp)
    target_include_directories(radray_pm PRIVATE ${FFTW3F_INCLUDE_DIR})
    target_link_libraries(radray_pm PUBLIC radray ${FFTW3F_LIBRARY})
    add_executable(pm_driver apps/pm_driver.cpp)
    target_link_libraries(pm_driver PRIVATE radray_pm)
    message(STATUS "PM gravity: enabled (FFTW3f: ${FFTW3F_LIBRARY})")
  else()
    message(STATUS "PM gravity: disabled (FFTW3f not found)")
  endif()
endif()

enable_testing()
add_subdirectory(tests)
