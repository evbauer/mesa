CODE := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
WORKSPACE := $(abspath $(CODE)/../../../../..)
MESA := $(WORKSPACE)/mesa
SDK := $(WORKSPACE)/mesasdk-x86_64-linux-26.6.1
PRIVATE := $(WORKSPACE)/freeeos_biquintic_progress/experiments/052_stellar_api/build
BUILD := $(WORKSPACE)/freeeos_biquintic_progress/experiments/052_stellar_api/runs/driver_build
FC := $(SDK)/bin/gfortran
empty :=
space := $(empty) $(empty)
PRIVATE_PC := $(subst $(space),:,$(wildcard $(PRIVATE)/*/lib/pkgconfig)):$(SDK)/lib/pkgconfig
MODULE_FLAGS := $(shell PKG_CONFIG_PATH='$(PRIVATE_PC)' pkg-config --define-variable=math_slot=crmath --cflags mesa-star)
LINK_FLAGS := $(shell PKG_CONFIG_PATH='$(PRIVATE_PC)' pkg-config --define-variable=math_slot=crmath --libs --static mesa-star)
FLAGS := -O2 -g -fopenmp -ftrampoline-impl=heap -ffp-contract=off -ffree-line-length-none -fcheck=all
COMPILE = $(FC) $(FLAGS) -I$(BUILD) -J$(BUILD) $(MODULE_FLAGS) -I$(PRIVATE)/auto_diff/modules
LIBS := -L$(SDK)/lib -L$(SDK)/lib64 -Wl,-rpath,$(SDK)/lib -Wl,-rpath,$(SDK)/lib64
CORES := $(MESA)/eos/eosFreeEOS_builder/biquintic
OBJECTS := $(addprefix $(BUILD)/,helm_polynomials.o freeeos_biquintic.o freeeos_table.o runtime_sha256.o freeeos_runtime_adapter.o run_star_extras.o run_star.o run.o)

.PHONY: all
all: $(BUILD)/star

$(BUILD):
	mkdir -p $@

$(BUILD)/helm_polynomials.o: $(MESA)/eos/private/helm_polynomials.f90 | $(BUILD)
	$(COMPILE) -c $< -o $@

$(BUILD)/freeeos_biquintic.o: $(CORES)/test_fortran/freeeos_biquintic.f90 $(BUILD)/helm_polynomials.o
	$(COMPILE) -c $< -o $@

$(BUILD)/freeeos_table.o: $(CORES)/table_fortran/freeeos_table.f90 $(BUILD)/freeeos_biquintic.o
	$(COMPILE) -c $< -o $@

$(BUILD)/runtime_sha256.o: $(CORES)/runtime_adapter/runtime_sha256.f90 | $(BUILD)
	$(COMPILE) -c $< -o $@

$(BUILD)/freeeos_runtime_adapter.o: $(CORES)/runtime_adapter/freeeos_runtime_adapter.f90 $(BUILD)/runtime_sha256.o $(BUILD)/freeeos_table.o
	$(COMPILE) -c $< -o $@

$(BUILD)/run_star_extras.o: $(CODE)/run_star_extras.f90 $(BUILD)/freeeos_runtime_adapter.o
	$(COMPILE) -c $< -o $@

$(BUILD)/run_star.o: $(MESA)/star/job/run_star.f90 $(BUILD)/run_star_extras.o
	$(COMPILE) -c $< -o $@

$(BUILD)/run.o: $(CORES)/stellar_control/run.f90 $(BUILD)/run_star.o
	$(COMPILE) -c $< -o $@

$(BUILD)/star: $(OBJECTS)
	$(FC) $(FLAGS) $^ -o $@ $(LIBS) $(LINK_FLAGS) -Wl,-Map,$(BUILD)/link.map
