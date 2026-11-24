WORKSPACE := $(abspath ../../../../../)
MESA := $(WORKSPACE)/mesa
SDK := $(WORKSPACE)/mesasdk-x86_64-linux-26.6.1
SOURCE := $(MESA)/eos/test/biquintic
EXPERIMENT := $(WORKSPACE)/freeeos_biquintic_progress/experiments/075_co_table_support
BUILD := $(EXPERIMENT)/build
FC := $(SDK)/bin/gfortran
empty :=
space := $(empty) $(empty)
PKG_CONFIG_PATH := $(subst $(space),:,$(wildcard $(MESA)/build/*/lib/pkgconfig)):$(SDK)/lib/pkgconfig
PACKAGES := mesa-eos mesa-const mesa-math mesa-chem mesa-utils
INCLUDES := -I$(BUILD) $(shell PKG_CONFIG_PATH='$(PKG_CONFIG_PATH)' pkg-config --define-variable=math_slot=crmath --cflags $(PACKAGES))
LIBRARIES := $(shell PKG_CONFIG_PATH='$(PKG_CONFIG_PATH)' pkg-config --define-variable=math_slot=crmath --libs --static $(PACKAGES))
NATIVE_BUILD := $(EXPERIMENT)/native_build
NATIVE_PC := $(subst $(space),:,$(wildcard $(NATIVE_BUILD)/*/lib/pkgconfig)):$(SDK)/lib/pkgconfig
NATIVE_INCLUDES := -I$(NATIVE_BUILD)/eos/modules $(shell PKG_CONFIG_PATH='$(NATIVE_PC)' pkg-config --define-variable=math_slot=crmath --cflags $(PACKAGES))
NATIVE_LIBRARIES := $(shell PKG_CONFIG_PATH='$(NATIVE_PC)' pkg-config --define-variable=math_slot=crmath --libs --static $(PACKAGES))
FLAGS := -O0 -g -fopenmp -ffree-line-length-none -fcheck=all -ffp-contract=off
OBJECTS := $(BUILD)/helm_polynomials.o $(BUILD)/arithmetic_weights.o $(BUILD)/freeeos_biquintic.o \
   $(BUILD)/freeeos_table.o $(BUILD)/freeeos_family.o $(BUILD)/test_co.o

all: $(BUILD)/test_co $(BUILD)/test_family

native: $(NATIVE_BUILD)/driver/native_driver

dual: $(NATIVE_BUILD)/driver/dual_handle_driver

$(BUILD):
	mkdir -p $@

$(BUILD)/helm_polynomials.o: $(MESA)/eos/private/helm_polynomials.f90 | $(BUILD)
	$(FC) $(FLAGS) $(INCLUDES) -I$(BUILD) -J$(BUILD) -c $< -o $@

$(BUILD)/arithmetic_weights.o: $(MESA)/eos/private/arithmetic_weights.f90 | $(BUILD)
	$(FC) $(FLAGS) $(INCLUDES) -I$(BUILD) -J$(BUILD) -c $< -o $@

$(BUILD)/freeeos_biquintic.o: $(MESA)/eos/private/freeeos_biquintic.f90 $(BUILD)/helm_polynomials.o $(BUILD)/arithmetic_weights.o
	$(FC) $(FLAGS) $(INCLUDES) -I$(MESA)/build/auto_diff/modules -I$(BUILD) -J$(BUILD) -c $< -o $@

$(BUILD)/freeeos_table.o: $(MESA)/eos/private/freeeos_table.f90 $(BUILD)/freeeos_biquintic.o
	$(FC) $(FLAGS) $(INCLUDES) -I$(MESA)/build/auto_diff/modules -I$(BUILD) -J$(BUILD) -c $< -o $@

$(BUILD)/freeeos_family.o: $(MESA)/eos/private/freeeos_family.f90 $(BUILD)/freeeos_table.o
	$(FC) $(FLAGS) $(INCLUDES) -I$(MESA)/build/auto_diff/modules -I$(BUILD) -J$(BUILD) -c $< -o $@

$(BUILD)/test_co.o: test_co.f90 $(BUILD)/freeeos_family.o
	$(FC) $(FLAGS) $(INCLUDES) -I$(BUILD) -J$(BUILD) -c $< -o $@

$(EXPERIMENT)/test_family.f90: $(WORKSPACE)/freeeos_biquintic_progress/experiments/063_family_runtime/test_family.f90
	cp $< $@

$(BUILD)/test_family.o: $(EXPERIMENT)/test_family.f90 $(BUILD)/freeeos_family.o
	$(FC) $(FLAGS) $(INCLUDES) -I$(BUILD) -J$(BUILD) -c $< -o $@

$(BUILD)/test_co: $(OBJECTS)
	$(FC) $(FLAGS) -Wl,--enable-new-dtags -o $@ $(OBJECTS) -L$(SDK)/lib -L$(SDK)/lib64 -Wl,-rpath,$(SDK)/lib -Wl,-rpath,$(SDK)/lib64 $(LIBRARIES)

$(BUILD)/test_family: $(filter-out $(BUILD)/test_co.o,$(OBJECTS)) $(BUILD)/test_family.o
	$(FC) $(FLAGS) -Wl,--enable-new-dtags -o $@ $^ -L$(SDK)/lib -L$(SDK)/lib64 -Wl,-rpath,$(SDK)/lib -Wl,-rpath,$(SDK)/lib64 $(LIBRARIES)

$(NATIVE_BUILD)/driver/native_driver: $(SOURCE)/native_integration/native_driver.f90 $(SOURCE)/runtime_bridge/bridge_support.f90 $(NATIVE_BUILD)/eos/lib/libeos.a
	mkdir -p $(NATIVE_BUILD)/driver/modules
	$(FC) -O2 -g -fopenmp -ffree-line-length-none -fcheck=all -ffp-contract=off \
	   $(NATIVE_INCLUDES) -J$(NATIVE_BUILD)/driver/modules \
	   $(SOURCE)/runtime_bridge/bridge_support.f90 $(SOURCE)/native_integration/native_driver.f90 \
	   -o $@ -Wl,--enable-new-dtags -L$(SDK)/lib -L$(SDK)/lib64 -Wl,-rpath,$(SDK)/lib \
	   -Wl,-rpath,$(SDK)/lib64 $(NATIVE_LIBRARIES)

$(NATIVE_BUILD)/driver/dual_handle_driver: dual_handle_driver.f90 $(SOURCE)/runtime_bridge/bridge_support.f90 $(NATIVE_BUILD)/eos/lib/libeos.a
	mkdir -p $(NATIVE_BUILD)/driver/modules
	$(FC) -O2 -g -fopenmp -ffree-line-length-none -fcheck=all -ffp-contract=off \
	   $(NATIVE_INCLUDES) -J$(NATIVE_BUILD)/driver/modules \
	   $(SOURCE)/runtime_bridge/bridge_support.f90 dual_handle_driver.f90 \
	   -o $@ -Wl,--enable-new-dtags -L$(SDK)/lib -L$(SDK)/lib64 -Wl,-rpath,$(SDK)/lib \
	   -Wl,-rpath,$(SDK)/lib64 $(NATIVE_LIBRARIES)
