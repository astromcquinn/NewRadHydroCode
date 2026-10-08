# Fortran reference harness for the RadHydro RT module (gfortran, 1 thread).
#
#   make                 # optimised build, N=32           -> build/N32_R64/rt_harness.x
#   make N=64            # other grid size                 -> build/N64_R64/rt_harness.x
#   make NRF=4           # Nray_max = 4*N^3                -> build/N32_R4/rt_harness.x
#   make debug           # -O0 -g -fcheck=all build        -> build/debug_N32_R64/rt_harness.x
#   make run-all         # build what is needed and run all test problems
#   make clean
#
# Compile-time parameters (preprocessor, see PATCHES.md):
#   N         -> -DNGRID   Ngas = Nrt = N           (default 32)
#   NRF       -> -DNRAY_MAX_FACTOR  Nray_max = NRF*N^3 (default 64)
#   NCPU      -> -DNCPU    Ncpu (OpenMP threads)    (default 1; validation uses 1)
#   NDM       -> -DNDM     Ndm (only sizes unused particle units; default N/2)
#   NHALO_MAX -> -DNHALO_MAX (unused halo buffer;    default (N/2)^3/1024)

ifeq ($(origin FC),default)
  FC = gfortran
endif
N       ?= 32
NCPU    ?= 1
NRF     ?= 64
MODE    ?= opt

DEFS = -DNGRID=$(N) -DNCPU=$(NCPU) -DNRAY_MAX_FACTOR=$(NRF)
ifdef NDM
  DEFS += -DNDM=$(NDM)
endif
ifdef NHALO_MAX
  DEFS += -DNHALO_MAX=$(NHALO_MAX)
endif

# Common flags: C preprocessor (was ifort -fpp), OpenMP (was -openmp),
# no line-length limit, -mcmodel=medium as in the original Makefile (large
# static arrays: ray_stat%ifree has Nray_max/(4*Ncpu) entries per thread).
COMMON  = -cpp -fopenmp -ffree-line-length-none -mcmodel=medium
OPT     = -O2 -g
DBG     = -O0 -g -fcheck=all -fbacktrace -Wall -Wno-conversion -Wno-unused-variable -Wno-unused-dummy-argument -Wno-unused-function
# FPE traps are off by default: photo_cell's ne iteration legitimately
# evaluates ne/ne_prev with ne_prev = 0 (Inf/NaN, IEEE semantics) in fully
# neutral cells.  `make debug FPE=overflow` (or invalid,zero,...) to trap.
ifdef FPE
  DBG += -ffpe-trap=$(FPE)
endif

# build directory: build/[debug_]N<N>_R<NRF>  (default build/N32_R64)
ifeq ($(MODE),debug)
  FFLAGS = $(COMMON) $(DBG)
  BUILD  = build/debug_N$(N)_R$(NRF)
else
  FFLAGS = $(COMMON) $(OPT)
  BUILD  = build/N$(N)_R$(NRF)
endif

# HEALPix 3.11 subset (pix2vec_nest / vec2pix_nest)
HPX_SRC  = healpix
HPX_MODS = healpix_types extension long_intrinsic misc_utils bit_manipulation \
           indmed num_rec statistics pix_tools
HPX_FFLAGS = -cpp -DGFORTRAN -DNOCFITSIO -O2 -g -ffree-line-length-none -fno-range-check
HPX_OBJ  = $(BUILD)/hpx/cgetEnvironment.o $(addprefix $(BUILD)/hpx/,$(addsuffix .o,$(HPX_MODS)))

# RT sources (order = module dependency order)
RT_MODS = ifport_shim physical_constants cosmo_parameters radhydro_global \
          general_tools cosmo_tools domain_tools data_tools io_tools \
          halo_tools ray_tracing_tools rt_rates rt_tools
RT_OBJ  = $(addprefix $(BUILD)/,$(addsuffix .o,$(RT_MODS)))

EXE = $(BUILD)/rt_harness.x

.PHONY: all debug clean run-all test-healpix print-exe
.NOTPARALLEL:

all: $(EXE)

debug:
	$(MAKE) MODE=debug N=$(N) NCPU=$(NCPU) NRF=$(NRF) FPE=$(FPE)

print-exe:
	@echo $(EXE)

$(EXE): $(HPX_OBJ) $(RT_OBJ) $(BUILD)/rt_harness.o
	$(FC) $(FFLAGS) -o $@ $^

$(BUILD)/hpx/%.o: $(HPX_SRC)/%.F90 | $(BUILD)/hpx
	$(FC) $(HPX_FFLAGS) -I$(HPX_SRC) -J$(BUILD)/hpx -c $< -o $@

$(BUILD)/hpx/%.o: $(HPX_SRC)/%.f90 | $(BUILD)/hpx
	$(FC) $(HPX_FFLAGS) -I$(HPX_SRC) -J$(BUILD)/hpx -c $< -o $@

$(BUILD)/hpx/cgetEnvironment.o: $(HPX_SRC)/cgetEnvironment.c | $(BUILD)/hpx
	gcc -O2 -c $< -o $@

$(BUILD)/%.o: src/%.f90 $(HPX_OBJ) | $(BUILD)
	$(FC) $(FFLAGS) $(DEFS) -I$(BUILD)/hpx -J$(BUILD) -c $< -o $@

$(BUILD) $(BUILD)/hpx:
	mkdir -p $@

# HEALPix check program (compares against healpy, see tests/check_healpix.py)
$(BUILD)/hpx/hpx_check.x: $(HPX_OBJ) tests/hpx_check.f90
	$(FC) $(HPX_FFLAGS) -I$(BUILD)/hpx -J$(BUILD)/hpx -o $@ $(HPX_OBJ) tests/hpx_check.f90

test-healpix: $(BUILD)/hpx/hpx_check.x
	python3 tests/check_healpix.py $(BUILD)/hpx/hpx_check.x

run-all:
	./run_all.sh

clean:
	rm -rf build
