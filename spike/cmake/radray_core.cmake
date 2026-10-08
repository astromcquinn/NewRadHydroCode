# Builds the parts of the radray library the spike needs for host-side setup
# (clock, units, rate tables, problem I/O) as `radray_core`, and the
# interface target `spike_common` (spike/common headers). Included by every
# spike project so each one is a standalone CMake build.
get_filename_component(RADRAY_ROOT "${CMAKE_CURRENT_LIST_DIR}/../.." ABSOLUTE)
if(NOT TARGET radray_core)
  add_library(radray_core STATIC
    ${RADRAY_ROOT}/src/healpix.cpp
    ${RADRAY_ROOT}/src/cosmology.cpp
    ${RADRAY_ROOT}/src/rt_rates.cpp
    ${RADRAY_ROOT}/src/ray_tracing.cpp
    ${RADRAY_ROOT}/src/rt_tools.cpp
    ${RADRAY_ROOT}/src/io.cpp)
  target_include_directories(radray_core PUBLIC ${RADRAY_ROOT}/include)
  target_compile_options(radray_core PRIVATE -ffp-contract=off -fno-fast-math)
  set_target_properties(radray_core PROPERTIES CXX_STANDARD 17 POSITION_INDEPENDENT_CODE ON)
  add_library(spike_common INTERFACE)
  target_include_directories(spike_common INTERFACE ${RADRAY_ROOT}/spike/common)
  target_link_libraries(spike_common INTERFACE radray_core)
endif()
