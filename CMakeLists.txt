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
  src/ray_dump.cpp
  src/ray_quant.cpp
  src/rt_moments.cpp
  src/rt_tensor_diag.cpp
  src/rt_moments_fft.cpp
  src/otvet_closure.cpp
  src/recomb_rad.cpp
  src/rt_backend.cpp
  src/pm_backend.cpp
  src/hydro_backend.cpp
  src/amrex_state.cpp
  src/thermal_history.cpp
  src/io.cpp)
target_include_directories(radray PUBLIC include)

# AMReX backend of the ray tracing (rt_backend = amrex, docs/stage_d.md):
# the sorted-ray-array RT on an AMReX level, MPI + OpenMP. Optional: with
# RADRAY_AMREX off the backend is a stub that throws. AMReX (3D, MPI, OMP,
# particles in SINGLE precision) is found with find_package (set
# -DAMReX_DIR=<install>/lib/cmake/AMReX) or, with RADRAY_FETCH_AMREX=ON,
# downloaded and built as a subproject (AMReX 26.10, needs network access).
# AMReX 26.x needs C++20: the backend is its own target; the library stays
# C++17.
option(RADRAY_AMREX "Build the AMReX RT backend (rt_backend = amrex)" OFF)
option(RADRAY_FETCH_AMREX "Download and build AMReX 26.10 if find_package fails" OFF)
set(RADRAY_HAVE_AMREX OFF)
if(RADRAY_AMREX)
  # (the PM backend, pm_backend = amrex, docs/stage_e_pm.md, adds the FFT
  # component: AMReX built with -DAMReX_FFT=ON, which needs FFTW3 on the CPU)
  find_package(AMReX 26.10 CONFIG QUIET COMPONENTS 3D MPI OMP PARTICLES PSINGLE FFT)
  if(NOT AMReX_FOUND)
    if(NOT RADRAY_FETCH_AMREX)
      message(FATAL_ERROR "RADRAY_AMREX=ON but AMReX (3D, MPI, OMP, PARTICLES, PSINGLE, FFT) was not "
                          "found: set -DAMReX_DIR=<install>/lib/cmake/AMReX or -DRADRAY_FETCH_AMREX=ON")
    endif()
    include(FetchContent)
    set(AMReX_SPACEDIM 3 CACHE STRING "" FORCE)
    set(AMReX_MPI ON CACHE BOOL "" FORCE)
    set(AMReX_OMP ON CACHE BOOL "" FORCE)
    set(AMReX_PARTICLES ON CACHE BOOL "" FORCE)
    set(AMReX_PARTICLES_PRECISION SINGLE CACHE STRING "" FORCE)
    set(AMReX_FORTRAN OFF CACHE BOOL "" FORCE)
    set(AMReX_FASTMATH OFF CACHE BOOL "" FORCE)
    set(AMReX_LINEAR_SOLVERS OFF CACHE BOOL "" FORCE)
    set(AMReX_FFT ON CACHE BOOL "" FORCE)
    set(AMReX_PROBINIT OFF CACHE BOOL "" FORCE)
    set(AMReX_AMRLEVEL OFF CACHE BOOL "" FORCE)
    if(RADRAY_AMREX_GPU_BACKEND)
      set(AMReX_GPU_BACKEND ${RADRAY_AMREX_GPU_BACKEND} CACHE STRING "" FORCE)
    endif()
    FetchContent_Declare(amrex
      GIT_REPOSITORY https://github.com/AMReX-Codes/amrex.git
      GIT_TAG 26.10
      GIT_SHALLOW TRUE)
    FetchContent_MakeAvailable(amrex)
  endif()
  set(RADRAY_HAVE_AMREX ON)
  add_library(radray_amrex src/amrex/rt_level.cpp src/amrex/pm_level.cpp src/amrex/hydro_level.cpp
                           src/amrex/gas_state.cpp src/amrex/moment_level.cpp)
  target_link_libraries(radray_amrex PUBLIC radray AMReX::amrex)
  # the PM level's pencil FFT uses FFTW3f directly (the library of the
  # serial port and of AMReX's own FFT module)
  find_library(RADRAY_AMREX_FFTW3F fftw3f REQUIRED)
  target_link_libraries(radray_amrex PRIVATE ${RADRAY_AMREX_FFTW3F})
  # rt_tools.cpp (in radray) calls the backend: every consumer of radray
  # links radray_amrex too (a cycle between static libraries, which CMake
  # resolves by repeating them on the link line)
  target_link_libraries(radray INTERFACE radray_amrex)
  if(TARGET AMReX::Flags_INLINE)
    # AMReX does not pass its inlining flag to consumers; the communication
    # loops are header templates instantiated here (spike/amrex/README.md)
    target_link_libraries(radray_amrex PRIVATE AMReX::Flags_INLINE)
  endif()
  set_target_properties(radray_amrex PROPERTIES CXX_STANDARD 20 CXX_STANDARD_REQUIRED ON)
  target_compile_definitions(radray_amrex PUBLIC RADRAY_HAVE_AMREX=1)
  if(AMReX_GPU_BACKEND STREQUAL "CUDA")
    # the ray kernels are compiled as CUDA for the device (docs/stage_d.md 5)
    set_source_files_properties(src/amrex/rt_level.cpp PROPERTIES LANGUAGE CUDA)
    set_target_properties(radray_amrex PROPERTIES CUDA_STANDARD 20 CUDA_SEPARABLE_COMPILATION ON)
    target_compile_options(radray_amrex PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--fmad=false>)
  endif()
  set(RADRAY_RT_BACKEND_LIB radray_amrex)
  # The stubs as a library of their own: the gfortran harness of the live
  # full-loop tests links libradray.a (whose rt_tools.cpp and halo code
  # reference the backend functions) without AMReX, so it links these
  # after it (reference/harness_full/Makefile, when the file exists).
  add_library(radray_backend_stubs STATIC src/amrex/rt_level_stub.cpp src/amrex/pm_level_stub.cpp
                                          src/amrex/hydro_level_stub.cpp src/amrex/gas_state_stub.cpp)
  target_link_libraries(radray_backend_stubs PUBLIC radray)
  message(STATUS "AMReX RT backend: enabled (AMReX ${AMReX_VERSION})")
else()
  target_sources(radray PRIVATE src/amrex/rt_level_stub.cpp src/amrex/pm_level_stub.cpp src/amrex/hydro_level_stub.cpp
                                src/amrex/gas_state_stub.cpp)
  set(RADRAY_RT_BACKEND_LIB radray)
  message(STATUS "AMReX RT backend: disabled (RADRAY_AMREX=OFF)")
endif()

# Moment RT (docs/rt_moments.md): the OTVET tensor uses an FFT convolution
# when FFTW3f is available (otherwise only the direct sum).
find_path(FFTW3F_INCLUDE_DIR fftw3.h)
find_library(FFTW3F_LIBRARY fftw3f)
if(FFTW3F_INCLUDE_DIR AND FFTW3F_LIBRARY)
  target_include_directories(radray PRIVATE ${FFTW3F_INCLUDE_DIR})
  target_link_libraries(radray PUBLIC ${FFTW3F_LIBRARY})
endif()

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
  src/mass_function.cpp
  src/abundance_matching.cpp
  src/xray_background.cpp)
target_link_libraries(radray_halos PUBLIC radray)
add_executable(find_halos apps/find_halos.cpp)
target_link_libraries(find_halos PRIVATE radray_halos)

add_executable(rt_driver apps/rt_driver.cpp)
target_link_libraries(rt_driver PRIVATE ${RADRAY_RT_BACKEND_LIB})

# The alternative OTVET tensors on a tensor-dump snapshot (docs/otvet_closure.md)
add_executable(otvet_tensor apps/otvet_tensor.cpp)
target_link_libraries(otvet_tensor PRIVATE radray)

# Ray-compression study prototype (docs/ray_compression.md); reads ray dumps
add_executable(ray_codec tools/ray_compression/ray_codec.cpp)
target_link_libraries(ray_codec PRIVATE radray)
target_compile_options(ray_codec PRIVATE -O3 -fno-math-errno)

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
    # Threaded FFTW for the refined PM mesh (pm_refine > 1, docs/pm_refine.md):
    # libfftw3f_omp from the same installation, if there is one. The default
    # (pm_refine = 1) never uses it.
    get_filename_component(_fftw_libdir ${FFTW3F_LIBRARY} DIRECTORY)
    find_library(FFTW3F_OMP_LIBRARY fftw3f_omp HINTS ${_fftw_libdir} NO_DEFAULT_PATH)
    if(FFTW3F_OMP_LIBRARY AND OpenMP_CXX_FOUND)
      target_link_libraries(radray_fft PUBLIC ${FFTW3F_OMP_LIBRARY} ${FFTW3F_LIBRARY})
      target_compile_definitions(radray_fft PRIVATE RADRAY_FFTW_OMP)
      # the PM level's pencil FFT (pm_level.cpp) calls the threaded API; the
      # level is reached from the AMReX library's own objects (the state's
      # load balancing), so every consumer of radray_amrex needs it too
      if(TARGET radray_amrex)
        target_link_libraries(radray_amrex PRIVATE ${FFTW3F_OMP_LIBRARY} ${FFTW3F_LIBRARY})
        # the moment level's threaded OTVET FFT plans (moment_level.cpp)
        target_compile_definitions(radray_amrex PRIVATE RADRAY_FFTW_OMP)
      endif()
      message(STATUS "Threaded FFTW (pm_refine > 1): ${FFTW3F_OMP_LIBRARY}")
    else()
      message(STATUS "Threaded FFTW: not found (pm_refine > 1 uses serial FFTs)")
    endif()
    # Hydro solver (hydro_tools.f90 -DMF), summarize_gas and the gas
    # domains, see docs/port_notes_hydro.md.
    add_library(radray_hydro src/hydro.cpp src/hydro_unsplit.cpp src/domain_gas.cpp)
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
    target_link_libraries(radhydro PRIVATE radray_pm radray_halos ${RADRAY_RT_BACKEND_LIB})
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

# Abundance matching of halo UV luminosities to observed UV LFs
# (docs/abundance_matching.md)
add_executable(am_fit apps/am_fit.cpp)
target_link_libraries(am_fit PRIVATE radray_halos)
add_executable(hmf_table apps/hmf_table.cpp)
target_link_libraries(hmf_table PRIVATE radray_halos)
# The X-ray background of xray_heating = 1, and a one-zone thermal history
# (docs/xray_heating.md)
add_executable(xray_bkg apps/xray_bkg.cpp)
target_link_libraries(xray_bkg PRIVATE radray_halos)
