# Builds the parts of the radray library the spike needs for host-side setup
# (clock, units, rate tables, problem I/O) as `radray_core`, and the
# interface target `spike_common` (spike/common headers). Included by every
# spike project so each one is a standalone CMake build.
get_filename_component(RADRAY_ROOT "${CMAKE_CURRENT_LIST_DIR}/../.." ABSOLUTE)
if(NOT TARGET radray_core)
  # The source list is read from the top-level CMakeLists.txt's
  # add_library(radray ...) block, so it cannot drift from the main library.
  file(READ "${RADRAY_ROOT}/CMakeLists.txt" _radray_top)
  string(REGEX MATCH "add_library\\(radray\n([^)]*)\\)" _radray_m "${_radray_top}")
  string(REGEX MATCHALL "src/[A-Za-z0-9_]+\\.cpp" _radray_srcs "${CMAKE_MATCH_1}")
  if(NOT _radray_srcs)
    message(FATAL_ERROR "radray_core: could not read add_library(radray ...) from ${RADRAY_ROOT}/CMakeLists.txt")
  endif()
  list(TRANSFORM _radray_srcs PREPEND "${RADRAY_ROOT}/")
  # the AMReX-free stub of the rt_backend = amrex functions (rt_backend.hpp)
  list(APPEND _radray_srcs "${RADRAY_ROOT}/src/amrex/rt_level_stub.cpp")
  add_library(radray_core STATIC ${_radray_srcs})
  target_include_directories(radray_core PUBLIC ${RADRAY_ROOT}/include)
  target_compile_options(radray_core PRIVATE -ffp-contract=off -fno-fast-math)
  set_target_properties(radray_core PROPERTIES CXX_STANDARD 17 POSITION_INDEPENDENT_CODE ON)
  # the moment RT (src/rt_moments_fft.cpp) uses FFTW3f when available
  find_path(FFTW3F_INCLUDE_DIR fftw3.h)
  find_library(FFTW3F_LIBRARY fftw3f)
  if(FFTW3F_INCLUDE_DIR AND FFTW3F_LIBRARY)
    target_include_directories(radray_core PRIVATE ${FFTW3F_INCLUDE_DIR})
    target_link_libraries(radray_core PUBLIC ${FFTW3F_LIBRARY})
  endif()
  add_library(spike_common INTERFACE)
  target_include_directories(spike_common INTERFACE ${RADRAY_ROOT}/spike/common)
  target_link_libraries(spike_common INTERFACE radray_core)
endif()
