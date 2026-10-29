FC = gfortran
FFLAGS = -O0 -g -Wall -Wextra -fcheck=all
FREEEOS_LIBS =
REFERENCE_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
BUILD_DIR := $(REFERENCE_DIR)build
SOURCE := $(REFERENCE_DIR)direct_freeeos.f90

.PHONY: all check
all: $(BUILD_DIR)/direct_freeeos

$(BUILD_DIR):
	mkdir -p "$@"

$(BUILD_DIR)/direct_freeeos.o: $(SOURCE) | $(BUILD_DIR)
	$(FC) $(FFLAGS) -c "$<" -o "$@"

check: $(BUILD_DIR)/direct_freeeos.o
	@printf '%s\n' 'Compile-only scaffold check; no FreeEOS calls or accuracy validation.'

$(BUILD_DIR)/direct_freeeos: $(BUILD_DIR)/direct_freeeos.o
	@test -n "$(strip $(FREEEOS_LIBS))" || { printf '%s\n' 'Set FREEEOS_LIBS from verified FreeEOS metadata; dependency is not bundled.'; exit 1; }
	$(FC) $(FFLAGS) "$<" $(FREEEOS_LIBS) -o "$@"
