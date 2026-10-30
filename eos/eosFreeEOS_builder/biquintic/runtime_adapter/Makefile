ADAPTER_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
MESA_DIR ?= $(abspath $(ADAPTER_DIR)/../../../..)
SDK_DIR ?= $(abspath $(MESA_DIR)/../mesasdk-x86_64-linux-26.6.1)
BUILD_DIR ?= $(abspath $(MESA_DIR)/../freeeos_biquintic_progress/experiments/047_runtime_adapter/build)
FC := $(SDK_DIR)/bin/gfortran
empty :=
space := $(empty) $(empty)
PKG_CONFIG_PATH := $(subst $(space),:,$(wildcard $(MESA_DIR)/build/*/lib/pkgconfig)):$(SDK_DIR)/lib/pkgconfig
PACKAGES := mesa-auto_diff mesa-chem mesa-const mesa-math mesa-utils
MODULE_FLAGS := $(shell PKG_CONFIG_PATH='$(PKG_CONFIG_PATH)' pkg-config --define-variable=math_slot=crmath --cflags $(PACKAGES))
LINK_FLAGS := $(shell PKG_CONFIG_PATH='$(PKG_CONFIG_PATH)' pkg-config --define-variable=math_slot=crmath --libs --static $(PACKAGES))
FFLAGS ?= -O2 -g -fopenmp -ftrampoline-impl=heap -ffp-contract=off -ffree-line-length-none -fcheck=all -Wall -Wextra -Wno-compare-reals
COMPILE = MESASDK_ROOT='$(SDK_DIR)' LD_LIBRARY_PATH='$(SDK_DIR)/lib:$(SDK_DIR)/lib64' $(FC) $(FFLAGS) $(MODULE_FLAGS) -I$(MESA_DIR)/build/eos/modules -I$(MESA_DIR)/build/auto_diff/modules -I$(BUILD_DIR) -J$(BUILD_DIR)
OBJECTS := $(addprefix $(BUILD_DIR)/,helm_polynomials.o freeeos_biquintic.o freeeos_table.o runtime_sha256.o freeeos_runtime_adapter.o)

.PHONY: all check-deps
all: $(BUILD_DIR)/runtime_driver

check-deps:
	@PKG_CONFIG_PATH='$(PKG_CONFIG_PATH)' pkg-config --exists $(PACKAGES)

$(BUILD_DIR):
	mkdir -p "$@"

$(BUILD_DIR)/helm_polynomials.o: $(MESA_DIR)/eos/private/helm_polynomials.f90 | $(BUILD_DIR) check-deps
	$(COMPILE) -c $< -o $@

$(BUILD_DIR)/freeeos_biquintic.o: $(ADAPTER_DIR)/../test_fortran/freeeos_biquintic.f90 $(BUILD_DIR)/helm_polynomials.o
	$(COMPILE) -c $< -o $@

$(BUILD_DIR)/freeeos_table.o: $(ADAPTER_DIR)/../table_fortran/freeeos_table.f90 $(BUILD_DIR)/freeeos_biquintic.o
	$(COMPILE) -c $< -o $@

$(BUILD_DIR)/runtime_sha256.o: $(ADAPTER_DIR)/runtime_sha256.f90 | $(BUILD_DIR)
	$(COMPILE) -c $< -o $@

$(BUILD_DIR)/freeeos_runtime_adapter.o: $(ADAPTER_DIR)/freeeos_runtime_adapter.f90 $(BUILD_DIR)/freeeos_table.o $(BUILD_DIR)/runtime_sha256.o
	$(COMPILE) -c $< -o $@

$(BUILD_DIR)/runtime_driver: $(ADAPTER_DIR)/runtime_driver.f90 $(OBJECTS)
	$(COMPILE) $< $(OBJECTS) -o $@ -L$(SDK_DIR)/lib -L$(SDK_DIR)/lib64 -Wl,-rpath,$(SDK_DIR)/lib -Wl,-rpath,$(SDK_DIR)/lib64 $(LINK_FLAGS)
