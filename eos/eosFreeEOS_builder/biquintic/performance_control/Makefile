CONTROL_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
MESA_DIR ?= $(abspath $(CONTROL_DIR)/../../../..)
SDK_DIR ?= $(abspath $(MESA_DIR)/../mesasdk-x86_64-linux-26.6.1)
TABLE_BUILD ?= $(MESA_DIR)/eos/eosFreeEOS_builder/biquintic/table_fortran/build
OUTPUT_DIR := $(abspath $(MESA_DIR)/../freeeos_biquintic_progress/experiments/017_performance_control)
BUILD_DIR ?= $(OUTPUT_DIR)/build
FC := $(SDK_DIR)/bin/gfortran
empty :=
space := $(empty) $(empty)
PKG_CONFIG_PATH := $(subst $(space),:,$(wildcard $(MESA_DIR)/build/*/lib/pkgconfig)):$(SDK_DIR)/lib/pkgconfig:$(shell printenv PKG_CONFIG_PATH)
export PKG_CONFIG_PATH
FLAGS := -O2 -g -std=f2008 -Wall -Wextra -Wno-compare-reals -Wno-unused-dummy-argument -ffree-line-length-none -ffp-contract=off -fopenmp -ftrampoline-impl=heap
ENV := MESASDK_ROOT="$(SDK_DIR)" LD_LIBRARY_PATH="$(SDK_DIR)/lib:$(SDK_DIR)/lib64:$$LD_LIBRARY_PATH"
TABLE_PACKAGES := mesa-auto_diff mesa-const
MESA_PACKAGES := mesa-eos mesa-chem mesa-math mesa-utils mesa-const
TABLE_MODULES := $(shell PKG_CONFIG_PATH='$(PKG_CONFIG_PATH)' pkg-config --cflags $(TABLE_PACKAGES)) -I$(MESA_DIR)/build/auto_diff/modules -I$(TABLE_BUILD)
MESA_MODULES := $(shell PKG_CONFIG_PATH='$(PKG_CONFIG_PATH)' pkg-config --define-variable=math_slot=crmath --cflags $(MESA_PACKAGES))
TABLE_LIBS := $(shell PKG_CONFIG_PATH='$(PKG_CONFIG_PATH)' pkg-config --libs --static $(TABLE_PACKAGES))
MESA_LIBS := $(shell PKG_CONFIG_PATH='$(PKG_CONFIG_PATH)' pkg-config --define-variable=math_slot=crmath --libs --static $(MESA_PACKAGES))
SDK_LIBS := -L$(SDK_DIR)/lib -L$(SDK_DIR)/lib64 -Wl,-rpath,$(SDK_DIR)/lib -Wl,-rpath,$(SDK_DIR)/lib64
TABLE_OBJECTS := $(addprefix $(TABLE_BUILD)/,freeeos_table.o freeeos_biquintic.o helm_polynomials.o)

.DEFAULT_GOAL := all
.DELETE_ON_ERROR:
.PHONY: all table mesa check-deps
all: table mesa
table: $(BUILD_DIR)/table_benchmark
mesa: $(BUILD_DIR)/mesa_benchmark

check-deps:
	@PKG_CONFIG_PATH='$(PKG_CONFIG_PATH)' pkg-config --exists $(TABLE_PACKAGES) || (echo 'Installed libraries required; no dependency rebuild.' && exit 1)
	@test -f "$(TABLE_BUILD)/freeeos_table.mod" || (echo 'Existing table_fortran build required; do not regenerate it here.' && exit 1)

$(BUILD_DIR):
	mkdir -p "$@"

$(BUILD_DIR)/benchmark_support.o: $(CONTROL_DIR)benchmark_support.f90 | $(BUILD_DIR) check-deps
	$(ENV) $(FC) $(FLAGS) $(TABLE_MODULES) -J$(BUILD_DIR) -c "$<" -o "$@"

$(BUILD_DIR)/table_benchmark: $(CONTROL_DIR)table_benchmark.f90 $(BUILD_DIR)/benchmark_support.o $(TABLE_OBJECTS)
	$(ENV) $(FC) $(FLAGS) $(TABLE_MODULES) -I$(BUILD_DIR) -J$(BUILD_DIR) "$<" $(BUILD_DIR)/benchmark_support.o $(TABLE_OBJECTS) -o "$@" $(SDK_LIBS) $(TABLE_LIBS)

$(BUILD_DIR)/mesa_benchmark: $(CONTROL_DIR)mesa_benchmark.f90 $(BUILD_DIR)/benchmark_support.o
	$(ENV) $(FC) $(FLAGS) $(MESA_MODULES) -I$(BUILD_DIR) -J$(BUILD_DIR) "$<" $(BUILD_DIR)/benchmark_support.o -o "$@" $(SDK_LIBS) $(MESA_LIBS)
