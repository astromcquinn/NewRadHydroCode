# AMReX prototype of the spike kernels (hydro_amrex, rt_amrex).
#
#   cmake -S spike/amrex -B build-spike-amrex \
#         -DAMReX_DIR=/path/to/amrex-install/lib/cmake/AMReX
#   cmake --build build-spike-amrex -j2
#
# AMReX must be built with MPI, OpenMP, particles in SINGLE precision
# (AMReX_PARTICLES_PRECISION=SINGLE) and 3 dimensions; see README.md.
# Without an installed AMReX, -DSPIKE_FETCH_AMREX=ON downloads and builds
# AMReX 26.10 as a subproject (needs network access).
cmake_minimum_required(VERSION 3.20)
project(spike_amrex CXX C)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

option(SPIKE_FETCH_AMREX "Download and build AMReX 26.10 if find_package fails" OFF)
find_package(AMReX 26.10 CONFIG QUIET COMPONENTS 3D MPI OMP PARTICLES PSINGLE)
if(NOT AMReX_FOUND)
  if(NOT SPIKE_FETCH_AMREX)
    message(FATAL_ERROR "AMReX (3D, MPI, OMP, PARTICLES, PSINGLE) not found: set "
                        "-DAMReX_DIR=<install>/lib/cmake/AMReX or -DSPIKE_FETCH_AMREX=ON")
  endif()
  include(FetchContent)
  set(AMReX_SPACEDIM 3 CACHE STRING "" FORCE)
  set(AMReX_MPI ON CACHE BOOL "" FORCE)
  set(AMReX_OMP ON CACHE BOOL "" FORCE)
  set(AMReX_PARTICLES ON CACHE BOOL "" FORCE)
  set(AMReX_PARTICLES_PRECISION SINGLE CACHE STRING "" FORCE)
  set(AMReX_FORTRAN OFF CACHE BOOL "" FORCE)
  set(AMReX_FASTMATH OFF CACHE BOOL "" FORCE)
  FetchContent_Declare(amrex GIT_REPOSITORY https://github.com/AMReX-Codes/amrex.git GIT_TAG 26.10)
  FetchContent_MakeAvailable(amrex)
else()
  message(STATUS "AMReX ${AMReX_VERSION} from ${AMReX_DIR}")
endif()

include(${CMAKE_CURRENT_LIST_DIR}/../cmake/radray_core.cmake)

foreach(app hydro_amrex rt_amrex)
  add_executable(${app} ${app}.cpp)
  target_link_libraries(${app} PRIVATE spike_common AMReX::amrex)
  set_target_properties(${app} PROPERTIES CXX_STANDARD 20 CXX_STANDARD_REQUIRED ON)
  target_compile_options(${app} PRIVATE -Wall -Wextra -ffp-contract=off -fno-fast-math)
  # AMReX builds itself with -finline-limit=43210 (GCC) but does not pass it
  # to consumers; its header templates (FillBoundary/SumBoundary copy loops)
  # are instantiated in our translation units and run 2x slower without it.
  # Inlining does not change floating-point results (-ffp-contract=off).
  if(TARGET AMReX::Flags_INLINE)
    target_link_libraries(${app} PRIVATE AMReX::Flags_INLINE)
  endif()
endforeach()
