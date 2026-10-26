include ../../../make/defaults-module.mk

MODULE_NAME := make-mesa-tables
SRCS := src/make_mesa_tables.f90
INTERNAL_DEPENDS_ON := const utils math chem eos
INCLUDE_DIRS :=
BINTYPE := executable

include ../../../make/Makefile

# Build then run the generator from this directory (reads inlist_mesa_tables,
# writes under $(MESA_DIR)/data/eosWDM_AV_data/).
tables: $(OBJ_OUT)
	$(OBJ_OUT)
