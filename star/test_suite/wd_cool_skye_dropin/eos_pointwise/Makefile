# WP8 point-wise eos driver (scratch branch only). Build with
#   env MESA_DIR=/Users/bauer39/work/mesas/evbauer/mesa make -C <this dir>
include ../../../../make/defaults-module.mk

MODULE_NAME := wp8-dropin-eval
SRCS := dropin_eval.f90
INTERNAL_DEPENDS_ON := const utils math chem auto_diff eos
# eos_wdm_av_def (drop-in counters) is a private eos module
INCLUDE_DIRS := -I$(MESA_DIR)/build/eos/modules
BINTYPE := executable

include ../../../../make/Makefile

all: ./dropin_eval
./dropin_eval: $(OBJ_OUT)
	cp $(OBJ_OUT) ./dropin_eval
