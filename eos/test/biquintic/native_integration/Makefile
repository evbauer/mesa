SOURCE := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
WORKSPACE := $(abspath $(SOURCE)/../../../../..)
MESA_DIR := $(WORKSPACE)/mesa
SDK_DIR := $(WORKSPACE)/mesasdk-x86_64-linux-26.6.1
export MESASDK_ROOT := $(SDK_DIR)
BUILD_ROOT := $(WORKSPACE)/freeeos_biquintic_progress/experiments/065_native_integration/build
DRIVER_DIR := $(BUILD_ROOT)/driver
FC := $(SDK_DIR)/bin/gfortran
empty :=
space := $(empty) $(empty)
PC := $(subst $(space),:,$(wildcard $(BUILD_ROOT)/*/lib/pkgconfig)):$(SDK_DIR)/lib/pkgconfig
PACKAGES := mesa-eos mesa-const mesa-math mesa-chem mesa-utils
CFLAGS := $(shell PKG_CONFIG_PATH='$(PC)' pkg-config --define-variable=math_slot=crmath --cflags $(PACKAGES))
LIBS := $(shell PKG_CONFIG_PATH='$(PC)' pkg-config --define-variable=math_slot=crmath --libs --static $(PACKAGES))
FLAGS := -O2 -g -fopenmp -ftrampoline-impl=heap -ffp-contract=off -ffree-line-length-none -fcheck=all

.PHONY: all
all: $(DRIVER_DIR)/native_driver

$(DRIVER_DIR)/native_driver: native_driver.f90 ../runtime_bridge/bridge_support.f90 Makefile $(BUILD_ROOT)/eos/lib/libeos.a
	mkdir -p $(DRIVER_DIR)/modules
	$(FC) $(FLAGS) $(CFLAGS) -I$(BUILD_ROOT)/eos/modules -J$(DRIVER_DIR)/modules ../runtime_bridge/bridge_support.f90 native_driver.f90 -o $@ -L$(SDK_DIR)/lib -L$(SDK_DIR)/lib64 -Wl,-rpath,$(SDK_DIR)/lib -Wl,-rpath,$(SDK_DIR)/lib64 $(LIBS)
