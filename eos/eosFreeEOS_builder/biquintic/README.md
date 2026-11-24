# FreeEOS biquintic table preprocessing

This directory builds, verifies and packages supplemental FreeEOS table data.
It is not part of the production EOS source list.

- `freeeos_potential.py` defines constraints, mixed-derivative stencils and
  the numerical policy. Its checksum is stored in the supplements.
- `supplemental_hdf5/produce.py`, `verify.py` and `package.py` generate and
  distribute the three extra mixed derivatives, support masks and metadata.
- `family_cache/produce.py` supplies the audited source reader and earlier
  text-payload producer; `family_cache/verify.py` and `export_cache.py` retain
  explicit preprocessing/comparison utilities.

Production interpolation, reconstruction, FORUM reads, composition selection
and runtime digest checks live in `eos/private/`. Enable the EOS control
`FreeEOS_use_HELM_interpolation`, default `.false.`. Normal users install
MESA with its distributed data; no preprocessing is required.

Development tests, evaluation drivers, figure scripts and historical prototypes
live in `eos/test/biquintic/`. Test-side links share these preprocessing sources
rather than duplicating the checksum-bearing policy implementation.

Run numerical preprocessing on the inherited compute node with fresh output
paths. Preserve distributed source tables and frozen project artifacts.
See [supplemental_hdf5/README.md](supplemental_hdf5/README.md) for tools and
[supplemental_hdf5/FORMAT.md](supplemental_hdf5/FORMAT.md) for the data contract.
