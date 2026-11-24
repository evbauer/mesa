MESA_DIR ?= $(abspath ../../../..)
FC := gfortran
PKG_CONFIG ?= pkg-config
INTERNAL_DEPENDS_ON := const utils math auto_diff
empty :=
space := $(empty) $(empty)
PKG_CONFIG_PATH := $(subst $(space),:,$(foreach package,$(INTERNAL_DEPENDS_ON),$(MESA_DIR)/build/$(package)/lib/pkgconfig)):$(PKG_CONFIG_PATH)
export PKG_CONFIG_PATH

FFLAGS ?= -O2 -g -std=f2008 -Wall -Wextra -Wno-compare-reals -Wno-unused-dummy-argument -fcheck=all -ffp-contract=off -fopenmp
BUILD_DIR ?= build/optimized
BIN_DIR ?= .
INCLUDE_DIRS = $(shell PKG_CONFIG_PATH='$(PKG_CONFIG_PATH)' $(PKG_CONFIG) --cflags mesa-auto_diff mesa-const) \
  -I$(MESA_DIR)/build/auto_diff/modules
LIB_DEP_ARGS = $(shell PKG_CONFIG_PATH='$(PKG_CONFIG_PATH)' $(PKG_CONFIG) --libs --static mesa-auto_diff mesa-const)
SDK_LIB_DIR = $(if $(MESASDK_ROOT),-L$(MESASDK_ROOT)/lib)
FCOMPILE = $(FC) $(FFLAGS) $(INCLUDE_DIRS) -I$(BUILD_DIR) -J$(BUILD_DIR) -c
EXECUTABLE = $(FC) $(FFLAGS) $(SDK_LIB_DIR)
KERNEL_OBJECTS := $(BUILD_DIR)/helm_polynomials.o $(BUILD_DIR)/freeeos_biquintic.o
TEST_EXECUTABLE := $(BIN_DIR)/test_biquintic
DRIVER_EXECUTABLE := $(BIN_DIR)/cell_driver

.DEFAULT_GOAL := all
.DELETE_ON_ERROR:
.PHONY: all test test-kernel test-driver check-deps
all: $(TEST_EXECUTABLE) $(DRIVER_EXECUTABLE)

check-deps:
	@$(PKG_CONFIG) --exists mesa-auto_diff mesa-const || \
	  (echo 'Existing MESA const/utils/math/auto_diff libraries are required; no MESA rebuild is attempted.' && exit 1)

$(sort $(BUILD_DIR) $(BIN_DIR)):
	mkdir -p $@

$(BUILD_DIR)/helm_polynomials.o: $(MESA_DIR)/eos/private/helm_polynomials.f90 | $(BUILD_DIR) check-deps
	$(FCOMPILE) $< -o $@

$(BUILD_DIR)/freeeos_biquintic.o: freeeos_biquintic.f90 $(BUILD_DIR)/helm_polynomials.o
	$(FCOMPILE) $< -o $@

$(BUILD_DIR)/test_biquintic.o: test_biquintic.f90 $(KERNEL_OBJECTS)
	$(FCOMPILE) $< -o $@

$(BUILD_DIR)/cell_driver.o: cell_driver.f90 $(KERNEL_OBJECTS)
	$(FCOMPILE) $< -o $@

$(TEST_EXECUTABLE): $(BUILD_DIR)/test_biquintic.o $(KERNEL_OBJECTS) | $(BIN_DIR)
	$(EXECUTABLE) -o $@ $^ $(LIB_DEP_ARGS)

$(DRIVER_EXECUTABLE): $(BUILD_DIR)/cell_driver.o $(KERNEL_OBJECTS) | $(BIN_DIR)
	$(EXECUTABLE) -o $@ $^ $(LIB_DEP_ARGS)

test: test-kernel test-driver

test-kernel: $(TEST_EXECUTABLE)
	$(TEST_EXECUTABLE)

test-driver: $(DRIVER_EXECUTABLE)
	bash test_driver.sh $(DRIVER_EXECUTABLE) $(BUILD_DIR)/driver_smoke
