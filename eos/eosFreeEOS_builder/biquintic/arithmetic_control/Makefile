CONTROL := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
WORKSPACE := $(abspath $(CONTROL)/../../../../..)
MESA_DIR := $(WORKSPACE)/mesa
SDK_DIR := $(WORKSPACE)/mesasdk-x86_64-linux-26.6.1
PRIVATE_EOS := $(WORKSPACE)/freeeos_biquintic_progress/experiments/050_runtime_api/build/eos
BUILD := $(WORKSPACE)/freeeos_biquintic_progress/experiments/053_arithmetic_control/build/driver
BRIDGE := $(CONTROL)/../runtime_bridge
FC := $(SDK_DIR)/bin/gfortran
empty :=
space := $(empty) $(empty)
BASE_PC := $(subst $(space),:,$(wildcard $(MESA_DIR)/build/*/lib/pkgconfig)):$(SDK_DIR)/lib/pkgconfig
PACKAGES := mesa-eos mesa-const mesa-math mesa-chem mesa-utils
BASE_FLAGS := $(shell PKG_CONFIG_PATH='$(BASE_PC)' pkg-config --define-variable=math_slot=crmath --cflags $(PACKAGES))
BASE_LIBS := $(shell PKG_CONFIG_PATH='$(BASE_PC)' pkg-config --define-variable=math_slot=crmath --libs --static $(PACKAGES))
FLAGS := -O2 -g -fopenmp -ftrampoline-impl=heap -ffp-contract=off -ffree-line-length-none -fcheck=all
LIBS := -L$(SDK_DIR)/lib -L$(SDK_DIR)/lib64 -Wl,-rpath,$(SDK_DIR)/lib -Wl,-rpath,$(SDK_DIR)/lib64
SOURCES := $(MESA_DIR)/eos/private/helm_polynomials.f90 $(CONTROL)/arithmetic_weights.f90 \
 $(CONTROL)/freeeos_biquintic.f90 $(CONTROL)/../table_fortran/freeeos_table.f90 \
 $(CONTROL)/../runtime_adapter/runtime_sha256.f90 $(CONTROL)/../runtime_adapter/freeeos_runtime_adapter.f90

.PHONY: all
all: $(BUILD)/candidate_driver $(BUILD)/test_biquintic $(BUILD)/cell_driver $(BUILD)/baseline_cell_driver $(BUILD)/runtime_driver

$(BUILD)/candidate_driver: $(SOURCES) $(BRIDGE)/bridge_support.f90 $(BRIDGE)/candidate_driver.f90
	mkdir -p $(BUILD)/modules
	$(FC) $(FLAGS) -I$(PRIVATE_EOS)/include $(BASE_FLAGS) -I$(MESA_DIR)/build/auto_diff/modules -I$(BUILD)/modules -J$(BUILD)/modules $(SOURCES) $(BRIDGE)/bridge_support.f90 $(BRIDGE)/candidate_driver.f90 -o $@ -L$(PRIVATE_EOS)/lib $(LIBS) $(BASE_LIBS)

$(BUILD)/test_biquintic: $(SOURCES) $(CONTROL)/../test_fortran/test_biquintic.f90
	mkdir -p $(BUILD)/test_modules
	$(FC) $(FLAGS) $(BASE_FLAGS) -I$(MESA_DIR)/build/auto_diff/modules -J$(BUILD)/test_modules $(wordlist 1,3,$(SOURCES)) $(CONTROL)/../test_fortran/test_biquintic.f90 -o $@ $(LIBS) $(BASE_LIBS)

$(BUILD)/cell_driver: $(SOURCES) $(CONTROL)/../test_fortran/cell_driver.f90
	mkdir -p $(BUILD)/cell_modules
	$(FC) $(FLAGS) $(BASE_FLAGS) -I$(MESA_DIR)/build/auto_diff/modules -J$(BUILD)/cell_modules $(wordlist 1,3,$(SOURCES)) $(CONTROL)/../test_fortran/cell_driver.f90 -o $@ $(LIBS) $(BASE_LIBS)

$(BUILD)/baseline_cell_driver: $(CONTROL)/../test_fortran/freeeos_biquintic.f90 $(CONTROL)/../test_fortran/cell_driver.f90
	mkdir -p $(BUILD)/baseline_cell_modules
	$(FC) $(FLAGS) $(BASE_FLAGS) -I$(MESA_DIR)/build/auto_diff/modules -J$(BUILD)/baseline_cell_modules $(MESA_DIR)/eos/private/helm_polynomials.f90 $^ -o $@ $(LIBS) $(BASE_LIBS)

$(BUILD)/runtime_driver: $(SOURCES) $(CONTROL)/../runtime_adapter/runtime_driver.f90
	mkdir -p $(BUILD)/runtime_modules
	$(FC) $(FLAGS) $(BASE_FLAGS) -I$(MESA_DIR)/build/auto_diff/modules -J$(BUILD)/runtime_modules $(SOURCES) $(CONTROL)/../runtime_adapter/runtime_driver.f90 -o $@ $(LIBS) $(BASE_LIBS)
