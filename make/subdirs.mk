include $(MAKE_DIR)/defaults-module.mk

SUBDIRS := const utils math mtx auto_diff num interp_1d interp_2d forum chem eos colors rates neu net gyre kap ionization atm turb star_data star astero binary

ifeq ($(WITH_ADIPLS),yes)
  SUBDIRS += adipls
endif
