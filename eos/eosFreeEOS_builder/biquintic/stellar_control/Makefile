CODE := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
WORKSPACE := $(abspath $(CODE)/../../../../..)
MESA := $(WORKSPACE)/mesa
SDK := $(WORKSPACE)/mesasdk-x86_64-linux-26.6.1
BUILD ?= $(WORKSPACE)/freeeos_biquintic_progress/experiments/049_stellar_preparation/build
FC := $(SDK)/bin/gfortran
empty :=
space := $(empty) $(empty)
BASE_PC := $(subst $(space),:,$(wildcard $(MESA)/build/*/lib/pkgconfig)):$(SDK)/lib/pkgconfig
PACKAGES := mesa-star
BASE_FLAGS := $(shell PKG_CONFIG_PATH='$(BASE_PC)' pkg-config --define-variable=math_slot=crmath --cflags $(PACKAGES))
BASE_LIBS := $(shell PKG_CONFIG_PATH='$(BASE_PC)' pkg-config --define-variable=math_slot=crmath --libs --static $(PACKAGES))
FLAGS := -O2 -g -fopenmp -ftrampoline-impl=heap -ffp-contract=off -ffree-line-length-none -fcheck=all
LIBS := -L$(SDK)/lib -L$(SDK)/lib64 -Wl,-rpath,$(SDK)/lib -Wl,-rpath,$(SDK)/lib64

.PHONY: all
all: $(BUILD)/star $(BUILD)/guard_driver

$(BUILD):
	mkdir -p $@

$(BUILD)/run_star_extras.o: $(CODE)/run_star_extras.f90 | $(BUILD)
	$(FC) $(FLAGS) $(BASE_FLAGS) -J$(BUILD) -I$(BUILD) -c $< -o $@

$(BUILD)/run_star.o: $(MESA)/star/job/run_star.f90 $(BUILD)/run_star_extras.o
	$(FC) $(FLAGS) $(BASE_FLAGS) -J$(BUILD) -I$(BUILD) -c $< -o $@

$(BUILD)/run.o: $(CODE)/run.f90 $(BUILD)/run_star.o
	$(FC) $(FLAGS) $(BASE_FLAGS) -J$(BUILD) -I$(BUILD) -c $< -o $@

$(BUILD)/star: $(BUILD)/run_star_extras.o $(BUILD)/run_star.o $(BUILD)/run.o
	$(FC) $(FLAGS) $^ -o $@ $(LIBS) $(BASE_LIBS)

$(BUILD)/guard_driver: $(CODE)/guard_driver.f90 $(BUILD)/run_star_extras.o
	$(FC) $(FLAGS) $(BASE_FLAGS) -I$(BUILD) -J$(BUILD) $^ -o $@ $(LIBS) $(BASE_LIBS)
