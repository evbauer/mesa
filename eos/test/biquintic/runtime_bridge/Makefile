BRIDGE := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
WORKSPACE := $(abspath $(BRIDGE)/../../../../..)
MESA_DIR := $(WORKSPACE)/mesa
SDK_DIR := $(WORKSPACE)/mesasdk-x86_64-linux-26.6.1
PRIVATE_EOS := $(WORKSPACE)/freeeos_biquintic_progress/experiments/050_runtime_api/build/eos
BUILD := $(PRIVATE_EOS)/../driver
FC := $(SDK_DIR)/bin/gfortran
empty :=
space := $(empty) $(empty)
BASE_PC := $(subst $(space),:,$(wildcard $(MESA_DIR)/build/*/lib/pkgconfig)):$(SDK_DIR)/lib/pkgconfig
PACKAGES := mesa-eos mesa-const mesa-math mesa-chem mesa-utils
BASE_FLAGS := $(shell PKG_CONFIG_PATH='$(BASE_PC)' pkg-config --define-variable=math_slot=crmath --cflags $(PACKAGES))
BASE_LIBS := $(shell PKG_CONFIG_PATH='$(BASE_PC)' pkg-config --define-variable=math_slot=crmath --libs --static $(PACKAGES))
FLAGS := -O2 -g -fopenmp -ftrampoline-impl=heap -ffp-contract=off -ffree-line-length-none -fcheck=all
LIBS := -L$(SDK_DIR)/lib -L$(SDK_DIR)/lib64 -Wl,-rpath,$(SDK_DIR)/lib -Wl,-rpath,$(SDK_DIR)/lib64

.PHONY: all
all: $(BUILD)/hook_driver $(BUILD)/baseline_driver $(BUILD)/candidate_driver

ADAPTER_SOURCES := $(MESA_DIR)/eos/private/helm_polynomials.f90 \
 $(BRIDGE)/../test_fortran/freeeos_biquintic.f90 \
 $(MESA_DIR)/eos/private/freeeos_table.f90 \
 $(MESA_DIR)/eos/private/runtime_sha256.f90 \
 $(BRIDGE)/../runtime_adapter/freeeos_runtime_adapter.f90

$(BUILD)/candidate_driver: $(ADAPTER_SOURCES) bridge_support.f90 candidate_driver.f90 $(PRIVATE_EOS)/lib/libeos.a
	mkdir -p $(BUILD)/candidate_modules
	$(FC) $(FLAGS) -I$(PRIVATE_EOS)/include $(BASE_FLAGS) -I$(MESA_DIR)/build/auto_diff/modules -I$(BUILD)/candidate_modules -J$(BUILD)/candidate_modules $(ADAPTER_SOURCES) bridge_support.f90 candidate_driver.f90 -o $@ -L$(PRIVATE_EOS)/lib $(LIBS) $(BASE_LIBS)

$(BUILD)/hook_driver: bridge_support.f90 hook_driver.f90 $(PRIVATE_EOS)/lib/libeos.a
	mkdir -p $(BUILD)/private_modules
	$(FC) $(FLAGS) -I$(PRIVATE_EOS)/include $(BASE_FLAGS) -J$(BUILD)/private_modules $< hook_driver.f90 -o $@ -L$(PRIVATE_EOS)/lib $(LIBS) $(BASE_LIBS)

$(BUILD)/baseline_driver: bridge_support.f90 baseline_driver.f90
	mkdir -p $(BUILD)/baseline_modules
	$(FC) $(FLAGS) $(BASE_FLAGS) -J$(BUILD)/baseline_modules $< baseline_driver.f90 -o $@ $(LIBS) $(BASE_LIBS)
