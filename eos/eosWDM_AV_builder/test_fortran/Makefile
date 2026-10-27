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

# The WDM_AV load/interp/mix modules are private eos modules that nothing inside
# eos uses yet, so MESA's build does not generate their .mod files.  Force their
# generation (via the eos sub-make) before this test's dependency scan runs.
.PHONY: force_wdm_mods
# Ask the eos build for the WDM_AV module files. This may report "no rule" on a
# second run when they already exist, so ignore the sub-make status and check
# for the files instead.
force_wdm_mods:
	-@$(MAKE) --no-print-directory -C $(MESA_DIR)/eos \
	   ../build/eos/modules/eos_wdm_av_def.mod \
	   ../build/eos/modules/eos_wdm_av_load.mod \
	   ../build/eos/modules/eos_wdm_av_interp.mod \
	   ../build/eos/modules/eos_wdm_av_mix.mod 2>/dev/null
	@test -f $(MESA_DIR)/build/eos/modules/eos_wdm_av_mix.mod || \
	   (echo "eos_wdm_av_*.mod not found; run make in $(MESA_DIR)/eos first" && exit 1)

$(BUILD_DIR_MODULE)/depend: force_wdm_mods

# copy the built executable here as ./test_interp so it can be run directly
all: ./test_interp ./test_mix
./test_interp: $(OBJ_OUT)
	cp $(OBJ_OUT) ./test_interp

# test_mix: a second standalone program built from the same flags/libs.  It is
# compiled and linked directly (the standard SRCS path builds only test_interp).
./test_mix: test_mix.f90 force_wdm_mods
	$(FCOMPILE) test_mix.f90 -o test_mix.o
	$(EXECUTABLE) -o ./test_mix test_mix.o $(LIB_DEP_ARGS)

# build and run the interpolation test
test: ./test_interp
	./test_interp

test2: ./test_mix
	./test_mix

