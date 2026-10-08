# Fortran reference harness for the RadHydro PM gravity (gfortran, 1 thread,
# FFTW3 single precision instead of MKL DFTI).
#
#   make                      # NGAS=32 NDM=16 MF=1  -> build/G32_D16_MF1/pm_harness.x
#   make NGAS=64 NDM=32       # other sizes
#   make MF=0                 # build without -DMF (gas kick on u_gas instead of v_gas)
#   make debug                # -O0 -g -fcheck=all build
#   make run-all              # build what is needed and run all tests/problems_pm
#   make clean
#
# Compile-time parameters (preprocessor; see PATCHES.md):
#   NGAS -> -DNGRID  Ngas (= Nrt)                 (default 32)
#   NDM  -> -DNDM    Ndm  (Npart = Ndm^3)          (default NGAS/2, as the original)
#   MF   -> -DMF     moving-frame branch, as the original Makefile (default 1)
#   Ncpu is fixed to 1 (-DNCPU=1).
#
# Sources: unchanged modules come from reference/Radhydro (pm_tools.f90,
# fft_tools.f90) and from reference/harness/src (the RT harness copies, whose
# patches are listed in reference/harness/PATCHES.md); hydro_tools.f90 (for
# init_hydro and summarize_gas) from reference/harness_hydro/src; the MKL_DFTI
# shim and calc_timestep from reference/harness_common; io_tools.f90 is a
# patched copy in src/ (PATCHES.md here).

ifeq ($(origin FC),default)
  FC = gfortran
endif
NGAS ?= 32
NDM  ?= $(shell echo $$(( $(NGAS) / 2 )))
MF   ?= 1
MODE ?= opt

FFTW_INC ?= /usr/include
FFTW_LIB ?= -lfftw3f

DEFS = -DNGRID=$(NGAS) -DNDM=$(NDM) -DNCPU=1 -DNRAY_MAX_FACTOR=1
ifeq ($(MF),1)
  DEFS += -DMF
endif

COMMON  = -cpp -fopenmp -ffree-line-length-none -mcmodel=medium -I$(FFTW_INC) -I$(COM)/include
OPT     = -O2 -g
DBG     = -O0 -g -fcheck=all -fbacktrace -Wall -Wno-conversion -Wno-unused-variable \
          -Wno-unused-dummy-argument -Wno-unused-function -Wno-maybe-uninitialized

ifeq ($(MODE),debug)
  FFLAGS = $(COMMON) $(DBG)
  BUILD  = build/debug_G$(NGAS)_D$(NDM)_MF$(MF)
else
  FFLAGS = $(COMMON) $(OPT)
  BUILD  = build/G$(NGAS)_D$(NDM)_MF$(MF)
endif

RTH  = ../harness/src
ORIG = ../Radhydro
COM  = ../harness_common

# module dependency order; directory of each source
SRC_ifport_shim       = $(RTH)/ifport_shim.f90
SRC_physical_constants = $(RTH)/physical_constants.f90
SRC_cosmo_parameters  = $(RTH)/cosmo_parameters.f90
SRC_radhydro_global   = $(RTH)/radhydro_global.f90
SRC_general_tools     = $(RTH)/general_tools.f90
SRC_cosmo_tools       = $(RTH)/cosmo_tools.f90
SRC_domain_tools      = $(RTH)/domain_tools.f90
SRC_mkl_dfti_shim     = $(COM)/mkl_dfti_shim.f90
SRC_harness_timestep  = $(COM)/harness_timestep.f90
SRC_fft_tools         = $(COM)/fft_tools.f90
SRC_hydro_tools       = ../harness_hydro/src/hydro_tools.f90
SRC_io_tools          = src/io_tools.f90
SRC_pm_tools          = $(ORIG)/pm_tools.f90
SRC_pm_harness        = src/pm_harness.f90

MODS = ifport_shim mkl_dfti_shim physical_constants cosmo_parameters radhydro_global \
       general_tools cosmo_tools harness_timestep domain_tools fft_tools io_tools \
       hydro_tools pm_tools pm_harness
OBJ  = $(addprefix $(BUILD)/,$(addsuffix .o,$(MODS)))

EXE = $(BUILD)/pm_harness.x

.PHONY: all debug clean run-all print-exe
.NOTPARALLEL:

all: $(EXE)

debug:
	$(MAKE) MODE=debug NGAS=$(NGAS) NDM=$(NDM) MF=$(MF)

print-exe:
	@echo $(EXE)

$(EXE): $(OBJ)
	$(FC) $(FFLAGS) -o $@ $^ $(FFTW_LIB)

define compile_rule
$(BUILD)/$(1).o: $(SRC_$(1)) | $(BUILD)
	$(FC) $(FFLAGS) $(DEFS) -J$(BUILD) -c $$< -o $$@
endef
$(foreach m,$(MODS),$(eval $(call compile_rule,$(m))))

$(BUILD):
	mkdir -p $@

run-all:
	./run_all.sh

clean:
	rm -rf build
