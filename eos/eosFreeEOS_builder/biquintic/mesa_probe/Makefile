PROBE_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
MESA_DIR ?= $(abspath $(PROBE_DIR)/../../../..)
SDK_DIR ?= $(abspath $(MESA_DIR)/../mesasdk-x86_64-linux-26.6.1)
FC := $(SDK_DIR)/bin/gfortran
empty :=
space := $(empty) $(empty)
PKG_CONFIG_PATH := $(subst $(space),:,$(wildcard $(MESA_DIR)/build/*/lib/pkgconfig)):$(SDK_DIR)/lib/pkgconfig:$(shell printenv PKG_CONFIG_PATH)
export PKG_CONFIG_PATH
PACKAGES := mesa-eos mesa-const mesa-math mesa-chem mesa-utils
MODULE_FLAGS := $(shell PKG_CONFIG_PATH='$(PKG_CONFIG_PATH)' pkg-config --define-variable=math_slot=crmath --cflags $(PACKAGES))
LINK_FLAGS := $(shell PKG_CONFIG_PATH='$(PKG_CONFIG_PATH)' pkg-config --define-variable=math_slot=crmath --libs --static $(PACKAGES))
FLAGS := -O2 -g -fopenmp -ftrampoline-impl=heap -ffp-contract=off -ffree-line-length-none

.PHONY: all
all: $(PROBE_DIR)build/mesa_probe

$(PROBE_DIR)build/mesa_probe: $(PROBE_DIR)probe.f90
	mkdir -p "$(PROBE_DIR)build"
	MESASDK_ROOT="$(SDK_DIR)" LD_LIBRARY_PATH="$(SDK_DIR)/lib:$(SDK_DIR)/lib64:$$LD_LIBRARY_PATH" $(FC) $(FLAGS) $(MODULE_FLAGS) "$<" -o "$@" -L$(SDK_DIR)/lib -L$(SDK_DIR)/lib64 -Wl,-rpath,$(SDK_DIR)/lib -Wl,-rpath,$(SDK_DIR)/lib64 $(LINK_FLAGS)
