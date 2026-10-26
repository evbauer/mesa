include ../../../make/defaults-module.mk

MODULE_NAME := eos-wdm-av-test
SRCS := test_interp.f90
# eos pulls in forum, auto_diff, num, interp_1d/2d transitively via pkg-config;
# auto_diff is listed explicitly because test_interp.f90 uses it directly.
INTERNAL_DEPENDS_ON := const utils math chem auto_diff eos
# the WDM_AV modules and helm_polynomials are private eos modules (not installed
# to the public include dir), so add the eos build modules directory
INCLUDE_DIRS := -I$(MESA_DIR)/build/eos/modules
BINTYPE := executable

include ../../../make/Makefile

# The WDM_AV load/interp modules are private eos modules that nothing inside eos
# uses yet, so MESA's build does not generate their .mod files.  Force their
# generation (via the eos sub-make) before this test's dependency scan runs.
.PHONY: force_wdm_mods
force_wdm_mods:
	@$(MAKE) --no-print-directory -C $(MESA_DIR)/eos \
	   ../build/eos/modules/eos_wdm_av_def.mod \
	   ../build/eos/modules/eos_wdm_av_load.mod \
	   ../build/eos/modules/eos_wdm_av_interp.mod

$(BUILD_DIR_MODULE)/depend: force_wdm_mods

# copy the built executable here as ./test_interp so it can be run directly
all: ./test_interp
./test_interp: $(OBJ_OUT)
	cp $(OBJ_OUT) ./test_interp

# build and run the interpolation test
test: ./test_interp
	./test_interp
