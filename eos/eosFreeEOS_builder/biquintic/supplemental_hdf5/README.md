# FreeEOS supplemental HDF5 preprocessing

These tools generate, verify and package the three extra mixed derivatives
needed by the biquintic Hermite interpolant. The family contains 89 supplements:
83 ordinary tables and six C/O variants. Production is a development step;
normal installation extracts the supplied archive beside the FreeEOS tables.

Runtime code lives in `eos/private/`: `freeeos_supplement.f90`,
`freeeos_legacy_nodes.f90`, `freeeos_table.f90`, `freeeos_biquintic.f90`,
`freeeos_family.f90`, `arithmetic_weights.f90` and `runtime_sha256.f90`.
It uses the standard legacy reader and binary caches, reconstructs six
constraints from nodal columns and reads only extra derivatives, masks and
metadata through FORUM. Enable it with
`FreeEOS_use_HELM_interpolation = .true.` in `&eos`; the default is `.false.`.

## Preprocessing tools

Run on the inherited compute node, using new output directories:

```sh
python3 -B produce.py --help
python3 -B verify.py --help
python3 -B package.py --help
```

The producer accepts source/output directories and a worker count. It uses
existing source tabulations rather than solving FreeEOS for new physics.
The verifier checks the family and stored identities. The packager creates
a compressed distribution and checks staged extraction. The numerical policy
remains in `../freeeos_potential.py`, with its original contents and checksum;
these file moves do not require rebuilding the tables.

See [FORMAT.md](FORMAT.md) for datasets, masks, coordinates and provenance.
Production retains the standard loader's known `02z60x`/`04z50x` Q-header
behavior: loaded count, origin and spacing define the grid, and shorter
loaded tables use the aligned supplement prefix.

## Validation and historical tools

Private builds, adapter tests, numerical comparison drivers and historical
raw-reader tests live in `eos/test/biquintic/supplemental_hdf5/`. Read its
README before running them. The historical source reader is test code and
is not linked into production. Scientific results and numbered experiments
are maintained in the separate `freeeos_biquintic_progress` repository.
