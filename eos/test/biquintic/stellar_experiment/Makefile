CODE := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
WORKSPACE := $(abspath $(CODE)/../../../../..)
MESA := $(WORKSPACE)/mesa
SDK := $(WORKSPACE)/mesasdk-x86_64-linux-26.6.1
PRIVATE := $(WORKSPACE)/freeeos_biquintic_progress/experiments/065_native_integration/build
BUILD := $(WORKSPACE)/freeeos_biquintic_progress/experiments/064_stellar_native/driver_build
FC := $(SDK)/bin/gfortran
empty :=
space := $(empty) $(empty)
PRIVATE_PC := $(subst $(space),:,$(wildcard $(PRIVATE)/*/lib/pkgconfig)):$(SDK)/lib/pkgconfig
MODULE_FLAGS := $(shell PKG_CONFIG_PATH='$(PRIVATE_PC)' pkg-config --define-variable=math_slot=crmath --cflags mesa-star)
LINK_FLAGS := $(shell PKG_CONFIG_PATH='$(PRIVATE_PC)' pkg-config --define-variable=math_slot=crmath --libs --static mesa-star)
FLAGS := -O2 -g -fopenmp -ftrampoline-impl=heap -ffp-contract=off -ffree-line-length-none -fcheck=all
COMPILE = $(FC) $(FLAGS) -I$(BUILD) -J$(BUILD) -I$(PRIVATE)/eos/modules $(MODULE_FLAGS)
LIBS := -L$(SDK)/lib -L$(SDK)/lib64 -Wl,-rpath,$(SDK)/lib -Wl,-rpath,$(SDK)/lib64
OBJECTS := $(addprefix $(BUILD)/,run_star_extras.o run_star.o run.o)

.PHONY: all
all: $(BUILD)/star

$(BUILD):
	mkdir -p $@

$(BUILD)/run_star_extras.o: $(CODE)/run_star_extras.f90 | $(BUILD)
	$(COMPILE) -c $< -o $@

$(BUILD)/run_star.o: $(MESA)/star/job/run_star.f90 $(BUILD)/run_star_extras.o
	$(COMPILE) -c $< -o $@

$(BUILD)/run.o: $(CODE)/run.f90 $(BUILD)/run_star.o
	$(COMPILE) -c $< -o $@

$(BUILD)/star: $(OBJECTS)
	$(FC) $(FLAGS) $^ -o $@ $(LIBS) $(LINK_FLAGS) -Wl,-Map,$(BUILD)/link.map
