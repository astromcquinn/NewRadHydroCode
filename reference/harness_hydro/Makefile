# Fortran reference harness for the RadHydro HYDRO solver (gfortran, 1 thread).
#
#   make                 # optimised build, N=32  -> build/N32/hydro_harness.x
#   make N=64            # other grid size        -> build/N64/hydro_harness.x
#   make debug           # -O0 -g -fcheck=all     -> build/debug_N32/hydro_harness.x
#   make clean
#
# Sources (see PATCHES.md): hydro_tools.f90, cosmo_tools.f90, domain_tools.f90,
# general_tools.f90, physical_constants.f90 and cosmo_parameters.f90 are
# compiled UNMODIFIED straight from ../Radhydro.  radhydro_global.f90 (sizes
# settable by -DNGRID/-DNCPU) and ifport_shim.f90 are reused from the RT
# harness (../harness/src).  ../harness_common/mkl_dfti_shim.f90 replaces Intel
# MKL DFTI by FFTW3 (single precision) and ../harness_common/harness_timestep.f90
# holds calc_timestep (both shared with the PM harness).  src/hydro_harness.f90
# is the driver.
#
# -DMF as in the original Makefile (RADHY = -DMF): moving-frame code paths.

ifeq ($(origin FC),default)
  FC = gfortran
endif
N       ?= 32
NCPU    ?= 1
MODE    ?= opt

ORIG    = ../Radhydro
RTH     = ../harness/src
COM     = ../harness_common

DEFS    = -DMF -DNGRID=$(N) -DNCPU=$(NCPU)
COMMON  = -cpp -fopenmp -ffree-line-length-none -mcmodel=medium
OPT     = -O2 -g
DBG     = -O0 -g -fcheck=all -fbacktrace -Wall -Wno-conversion -Wno-unused-variable \
          -Wno-unused-dummy-argument -Wno-unused-function -Wno-compare-reals \
          -Wno-unused-label -Wno-maybe-uninitialized
ifdef FPE
  DBG += -ffpe-trap=$(FPE)
endif

ifeq ($(MODE),debug)
  FFLAGS = $(COMMON) $(DBG)
  BUILD  = build/debug_N$(N)
else
  FFLAGS = $(COMMON) $(OPT)
  BUILD  = build/N$(N)
endif

FFTW_INC ?= /usr/include
FFTW_LIB ?= -lfftw3f

# module dependency order: source path of each object
SRCS = $(RTH)/ifport_shim.f90 $(COM)/mkl_dfti_shim.f90 \
       $(ORIG)/physical_constants.f90 $(ORIG)/cosmo_parameters.f90 \
       $(RTH)/radhydro_global.f90 $(ORIG)/general_tools.f90 \
       $(ORIG)/cosmo_tools.f90 $(COM)/harness_timestep.f90 \
       $(ORIG)/domain_tools.f90 src/hydro_tools.f90 \
       src/hydro_harness.f90
OBJS = $(addprefix $(BUILD)/,$(notdir $(SRCS:.f90=.o)))

EXE = $(BUILD)/hydro_harness.x

.PHONY: all debug clean print-exe
.NOTPARALLEL:

all: $(EXE)

debug:
	$(MAKE) MODE=debug N=$(N) NCPU=$(NCPU) FPE=$(FPE)

print-exe:
	@echo $(EXE)

$(EXE): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^ $(FFTW_LIB)

# one explicit rule per source (sources live in three directories)
define compile_rule
$(BUILD)/$(notdir $(1:.f90=.o)): $(1) | $(BUILD)
	$(FC) $(FFLAGS) $(DEFS) -I$(FFTW_INC) -J$(BUILD) -c $(1) -o $$@
endef
$(foreach s,$(SRCS),$(eval $(call compile_rule,$(s))))

# every object depends on the modules compiled before it
$(foreach i,$(shell seq 2 $(words $(OBJS))),$(eval $(word $(i),$(OBJS)): $(word $(shell expr $(i) - 1),$(OBJS))))

$(BUILD):
	mkdir -p $@

clean:
	rm -rf build
