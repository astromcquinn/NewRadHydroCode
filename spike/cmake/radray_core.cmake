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
  add_library(radray_core STATIC ${_radray_srcs})
  target_include_directories(radray_core PUBLIC ${RADRAY_ROOT}/include)
  target_compile_options(radray_core PRIVATE -ffp-contract=off -fno-fast-math)
  set_target_properties(radray_core PROPERTIES CXX_STANDARD 17 POSITION_INDEPENDENT_CODE ON)
  add_library(spike_common INTERFACE)
  target_include_directories(spike_common INTERFACE ${RADRAY_ROOT}/spike/common)
  target_link_libraries(spike_common INTERFACE radray_core)
endif()
