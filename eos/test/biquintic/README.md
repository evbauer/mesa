# FreeEOS biquintic development and test harness

The development and test tools live in `mesa/eos/test/biquintic`. Production
Fortran modules live in `mesa/eos/private`: `runtime_sha256`,
`arithmetic_weights`, `freeeos_biquintic` (the AD-control variant),
`freeeos_table`, `freeeos_legacy_nodes`, `freeeos_supplement` and
`freeeos_family`. Build recipes refer to those modules directly; historical
Fortran variants retained here are comparison inputs.

Preprocessing and its policy files remain in
[`eos/eosFreeEOS_builder/biquintic`](../../eosFreeEOS_builder/biquintic).
The Python symlinks in this test tree preserve imports and the actual policy
checksum. Do not edit their targets through this harness.

## Current production validation

Use the [supplemental runtime testing instructions](supplemental_hdf5/README.md)
for private build scripts, production adapter tests and FORUM fixture checks.
Production loading reuses the standard FreeEOS table reader, loaded nodal
arrays and existing binary caches. Only mixed derivatives, masks and metadata
come from HDF5 through FORUM. The inlist switch is now
`FreeEOS_use_HELM_interpolation` (default false); the supplement override
remains `FreeEOS_biquintic_table_dir`.

Run computational checks on the allocated compute node. Coordinate CPUs and
memory with the lead, save outputs in fresh directories, and leave enough
allocation time for a clean checkpoint. The lead owns production builds and
before/after numerical parity; layout checks alone establish no scientific
or numerical equivalence.

## Preserved historical tools

The other subdirectories preserve earlier standalone adapters, cell/table
controls, stellar harnesses and saved-data analyses. Source paths and live
inlist templates follow the current layout and switch name. Some tools retain
fixed historical hashes, binaries, allocation identities and deadlines. They
can reject current sources or renamed templates; do not replace expected
hashes to make a historical provenance check pass. Archived inlists retain
the earlier switch name. These tools are not a single current acceptance suite.

`test_fortran/freeeos_biquintic.f90` and
`arithmetic_control/freeeos_biquintic.f90` are historical variants, distinct
from the production AD-control module. `runtime_adapter/` is the historical
fixed-mixture standalone adapter; current production legacy-EP and family
checks are under `supplemental_hdf5/`. Historical `co_table_support/` standalone
family tests predate the required legacy request context; use the current
`test_legacy_family` for production acceptance.

README commands in historical directories describe the associated workflow.
They may point at frozen experiment outputs and require earlier binaries or
inputs. Do not run them unmodified or overwrite their evidence. `run_validation.py`
is an earlier aggregate regression driver, not the production adapter gate.
The current scientific assessment is owned by the lead in
[the main report](../../../../freeeos_biquintic_progress/freeeos_biquintic_progress.pdf).
