# Parthenon prototype of the spike (RT substep + MUSCL-Hancock hydro).
#
#   cmake -S spike/parthenon -B build-spike-parthenon \
#         -DPARTHENON_SOURCE_DIR=/path/to/parthenon      # a v25.12 checkout
#   cmake --build build-spike-parthenon -j2
#
# Parthenon (and its bundled Kokkos) is built as part of this project via
# add_subdirectory, the way AthenaPK consumes it. Without PARTHENON_SOURCE_DIR
# it is fetched with FetchContent at PARTHENON_TAG (default v25.12); with
# FETCHCONTENT_SOURCE_DIR_PARTHENON the fetch uses a local checkout instead.
cmake_minimum_required(VERSION 3.16)
project(spike_parthenon C CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS_RELEASE "-O2")  # same as spike/serial
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

# --- Parthenon configuration (forced: the spike needs exactly this) ---------
set(Kokkos_ENABLE_OPENMP ON CACHE BOOL "" FORCE)
set(PARTHENON_DISABLE_HDF5 ON CACHE BOOL "" FORCE)
set(PARTHENON_DISABLE_EXAMPLES ON CACHE BOOL "" FORCE)
set(PARTHENON_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
set(PARTHENON_ENABLE_PYTHON_MODULE_CHECK OFF CACHE BOOL "" FORCE)
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(PARTHENON_SINGLE_PRECISION OFF CACHE BOOL "" FORCE)  # Real = double

set(PARTHENON_SOURCE_DIR "" CACHE PATH "Local Parthenon checkout (else FetchContent)")
set(PARTHENON_TAG "v25.12" CACHE STRING "Parthenon tag for FetchContent")
if(PARTHENON_SOURCE_DIR)
  add_subdirectory(${PARTHENON_SOURCE_DIR} parthenon)
else()
  include(FetchContent)
  FetchContent_Declare(parthenon
    GIT_REPOSITORY https://github.com/parthenon-hpc-lab/parthenon.git
    GIT_TAG ${PARTHENON_TAG}
    GIT_SUBMODULES_RECURSE ON)
  FetchContent_MakeAvailable(parthenon)
endif()

# --- our code -----------------------------------------------------------------
include(${CMAKE_CURRENT_LIST_DIR}/../cmake/radray_core.cmake)
find_package(MPI REQUIRED COMPONENTS CXX)

foreach(app rt_parthenon hydro_parthenon)
  add_executable(${app} ${app}.cpp)
  target_link_libraries(${app} PRIVATE parthenon spike_common MPI::MPI_CXX)
  target_compile_options(${app} PRIVATE -ffp-contract=off -fno-fast-math)
endforeach()
