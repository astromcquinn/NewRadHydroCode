# Fortran reference build of the FULL original radhydro.f90 (gfortran,
# 1 thread, FFTW3f instead of MKL): hydro + PM + halos + RT.
#
#   make                       # N=32 NDM=16 NRF=64 BOX=5 -> build/N32_D16_R64_B5/radhydro_full.x
#   make N=32 NDM=16 NRF=1 BOX=2
#   make debug ...             # -O0 -g -fcheck=all
#   make ONTHEFLY=0            # without the on-the-fly halo finder (no C++ needed)
#   make clean
#
# Compile-time parameters of the original (preprocessor, see PATCHES.md):
#   N         -> -DNGRID            Ngas = Nrt
#   NDM       -> -DNDM              Ndm (default N/2, as the original)
#   NRF       -> -DNRAY_MAX_FACTOR  Nray_max = NRF*N^3 (default 64)
#   BOX       -> -DBOX              box in Mpc/h (default 5; a literal exact in real(4))
#   NHALO_MAX -> -DNHALO_MAX        source(:) capacity (default 16384; the original
#                                   Npart/1024 is 4 for Ndm = 16)
#
# The on-the-fly source mode (ONTHEFLY=1, default) links the C++ halo finder
# and source model (src/halo_capi.cpp) from the CMake build directory
# RADRAY_BUILD (default ../../build; libradray_halos.a, libradray.a).
#
# Sources: unchanged originals from ../Radhydro (pm_tools.f90); the RT
# harness copies ../harness/src (radhydro_global, rt_tools, rt_rates,
# ray_tracing_tools, ... with their PATCHES.md); hydro_tools.f90 from
# ../harness_hydro/src; io_tools.f90 from ../harness_pm/src (Ncpu = 1
# read_init_conditions fix); the MKL_DFTI shim and fft_tools.f90 from
# ../harness_common; cosmo_parameters.f90 (BOX) and the driver
# radhydro_full.f90 (= radhydro.f90 + HARNESS_FULL patches) in src/.

ifeq ($(origin FC),default)
  FC = gfortran
endif
N         ?= 32
NDM       ?= $(shell echo $$(( $(N) / 2 )))
NRF       ?= 64
BOX       ?= 5
NHALO_MAX ?= 16384
ONTHEFLY  ?= 1
MODE      ?= opt
RADRAY_BUILD ?= ../../build

FFTW_INC ?= /usr/include
FFTW_LIB ?= -lfftw3f

ORIG = ../Radhydro
RTH  = ../harness/src
COM  = ../harness_common

DEFS = -DMF -DNGRID=$(N) -DNDM=$(NDM) -DNCPU=1 -DNRAY_MAX_FACTOR=$(NRF) -DBOX=$(BOX) \
       -DNHALO_MAX=$(NHALO_MAX) -DHARNESS_FULL_F1
ifeq ($(ONTHEFLY),1)
  DEFS    += -DONTHEFLY
  # libradray_backend_stubs.a exists in a build with the AMReX backends,
  # whose libradray.a references them (resolved here by the stubs: the
  # harness never runs a backend)
  CXXDEPS  = $(RADRAY_BUILD)/libradray_halos.a $(RADRAY_BUILD)/libradray.a \
             $(wildcard $(RADRAY_BUILD)/libradray_backend_stubs.a)
  CXXLIBS  = $(CXXDEPS) -lstdc++ -lm
endif

COMMON = -cpp -fopenmp -ffree-line-length-none -mcmodel=medium -I$(FFTW_INC) -I$(COM)/include
OPT    = -O2 -g
DBG    = -O0 -g -fcheck=all -fbacktrace -Wall -Wno-conversion -Wno-unused-variable \
         -Wno-unused-dummy-argument -Wno-unused-function -Wno-compare-reals \
         -Wno-unused-label -Wno-maybe-uninitialized -Wno-character-truncation
HPX_FFLAGS = -cpp -DGFORTRAN -DNOCFITSIO -O2 -g -ffree-line-length-none -fno-range-check

TAG = N$(N)_D$(NDM)_R$(NRF)_B$(BOX)_F$(ONTHEFLY)
ifeq ($(MODE),debug)
  FFLAGS = $(COMMON) $(DBG)
  BUILD  = build/debug_$(TAG)
else
  FFLAGS = $(COMMON) $(OPT)
  BUILD  = build/$(TAG)
endif

# HEALPix 3.11 subset (as the RT harness)
HPX_SRC  = ../harness/healpix
HPX_MODS = healpix_types extension long_intrinsic misc_utils bit_manipulation \
           indmed num_rec statistics pix_tools
HPX_OBJ  = $(BUILD)/hpx/cgetEnvironment.o $(addprefix $(BUILD)/hpx/,$(addsuffix .o,$(HPX_MODS)))

# module dependency order: source path of each object
SRCS = $(RTH)/ifport_shim.f90 $(COM)/mkl_dfti_shim.f90 \
       $(RTH)/physical_constants.f90 src/cosmo_parameters.f90 \
       $(RTH)/radhydro_global.f90 $(RTH)/general_tools.f90 $(RTH)/cosmo_tools.f90 \
       $(RTH)/domain_tools.f90 $(RTH)/data_tools.f90 ../harness_pm/src/io_tools.f90 \
       $(RTH)/halo_tools.f90 $(COM)/fft_tools.f90 ../harness_hydro/src/hydro_tools.f90 \
       $(ORIG)/pm_tools.f90 $(RTH)/ray_tracing_tools.f90 $(RTH)/rt_rates.f90 \
       $(RTH)/rt_tools.f90 src/harness_full_tools.f90 src/radhydro_full.f90
OBJS = $(addprefix $(BUILD)/,$(notdir $(SRCS:.f90=.o)))

EXE = $(BUILD)/radhydro_full.x

.PHONY: all debug clean print-exe
.NOTPARALLEL:

all: $(EXE)

debug:
	$(MAKE) MODE=debug N=$(N) NDM=$(NDM) NRF=$(NRF) BOX=$(BOX) NHALO_MAX=$(NHALO_MAX) \
	        ONTHEFLY=$(ONTHEFLY) RADRAY_BUILD=$(RADRAY_BUILD)

print-exe:
	@echo $(EXE)

$(EXE): $(HPX_OBJ) $(OBJS) $(CXXDEPS)
	$(FC) $(FFLAGS) -o $@ $(HPX_OBJ) $(OBJS) $(FFTW_LIB) $(CXXLIBS)

$(BUILD)/hpx/%.o: $(HPX_SRC)/%.F90 | $(BUILD)/hpx
	$(FC) $(HPX_FFLAGS) -I$(HPX_SRC) -J$(BUILD)/hpx -c $< -o $@

$(BUILD)/hpx/%.o: $(HPX_SRC)/%.f90 | $(BUILD)/hpx
	$(FC) $(HPX_FFLAGS) -I$(HPX_SRC) -J$(BUILD)/hpx -c $< -o $@

$(BUILD)/hpx/cgetEnvironment.o: $(HPX_SRC)/cgetEnvironment.c | $(BUILD)/hpx
	gcc -O2 -c $< -o $@

# one explicit rule per source (sources live in several directories)
define compile_rule
$(BUILD)/$(notdir $(1:.f90=.o)): $(1) $(HPX_OBJ) | $(BUILD)
	$(FC) $(FFLAGS) $(DEFS) -I$(BUILD)/hpx -J$(BUILD) -c $(1) -o $$@
endef
$(foreach s,$(SRCS),$(eval $(call compile_rule,$(s))))

# every object depends on the modules compiled before it
$(foreach i,$(shell seq 2 $(words $(OBJS))),$(eval $(word $(i),$(OBJS)): $(word $(shell expr $(i) - 1),$(OBJS))))

$(BUILD) $(BUILD)/hpx:
	mkdir -p $@

clean:
	rm -rf build
