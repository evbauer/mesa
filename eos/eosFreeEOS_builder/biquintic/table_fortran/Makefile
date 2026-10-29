MESA_DIR ?= $(abspath ../../../..)
FC := gfortran
PKG_CONFIG ?= pkg-config
INTERNAL_DEPENDS_ON := const utils math auto_diff
empty :=
space := $(empty) $(empty)
PKG_CONFIG_PATH := $(subst $(space),:,$(foreach package,$(INTERNAL_DEPENDS_ON),$(MESA_DIR)/build/$(package)/lib/pkgconfig)):$(PKG_CONFIG_PATH)
export PKG_CONFIG_PATH

FFLAGS ?= -O2 -g -std=f2008 -Wall -Wextra -Wno-compare-reals -Wno-unused-dummy-argument -fcheck=all -ffp-contract=off -fopenmp
BUILD_DIR ?= build
INCLUDE_DIRS = $(shell PKG_CONFIG_PATH='$(PKG_CONFIG_PATH)' $(PKG_CONFIG) --cflags mesa-auto_diff mesa-const) -I$(MESA_DIR)/build/auto_diff/modules
LIB_DEP_ARGS = $(shell PKG_CONFIG_PATH='$(PKG_CONFIG_PATH)' $(PKG_CONFIG) --libs --static mesa-auto_diff mesa-const)
SDK_LIB_DIR = $(if $(MESASDK_ROOT),-L$(MESASDK_ROOT)/lib)
FCOMPILE = $(FC) $(FFLAGS) $(INCLUDE_DIRS) -I$(BUILD_DIR) -J$(BUILD_DIR) -c
OBJECTS := $(BUILD_DIR)/helm_polynomials.o $(BUILD_DIR)/freeeos_biquintic.o $(BUILD_DIR)/freeeos_table.o

.DEFAULT_GOAL := all
.DELETE_ON_ERROR:
.PHONY: all test check-deps
all: $(BUILD_DIR)/table_query $(BUILD_DIR)/test_table

check-deps:
	@PKG_CONFIG_PATH='$(PKG_CONFIG_PATH)' $(PKG_CONFIG) --exists mesa-auto_diff mesa-const || \
	  (echo 'Installed MESA dependencies required; no dependency rebuild is attempted.' && exit 1)

$(BUILD_DIR):
	mkdir -p $@

$(BUILD_DIR)/helm_polynomials.o: $(MESA_DIR)/eos/private/helm_polynomials.f90 | $(BUILD_DIR) check-deps
	$(FCOMPILE) $< -o $@

$(BUILD_DIR)/freeeos_biquintic.o: ../test_fortran/freeeos_biquintic.f90 $(BUILD_DIR)/helm_polynomials.o
	$(FCOMPILE) $< -o $@

$(BUILD_DIR)/freeeos_table.o: freeeos_table.f90 $(BUILD_DIR)/freeeos_biquintic.o
	$(FCOMPILE) $< -o $@

$(BUILD_DIR)/table_query.o: table_query.f90 $(OBJECTS)
	$(FCOMPILE) $< -o $@

$(BUILD_DIR)/test_table.o: test_table.f90 $(OBJECTS)
	$(FCOMPILE) $< -o $@

$(BUILD_DIR)/table_query: $(BUILD_DIR)/table_query.o $(OBJECTS)
	$(FC) $(FFLAGS) $(SDK_LIB_DIR) -o $@ $^ $(LIB_DEP_ARGS)

$(BUILD_DIR)/test_table: $(BUILD_DIR)/test_table.o $(OBJECTS)
	$(FC) $(FFLAGS) $(SDK_LIB_DIR) -o $@ $^ $(LIB_DEP_ARGS)

test: all
	$(BUILD_DIR)/test_table $(BUILD_DIR)
	python3 test_driver.py $(BUILD_DIR)/table_query $(BUILD_DIR)
