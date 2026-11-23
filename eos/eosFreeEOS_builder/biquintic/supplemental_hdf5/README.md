# Source-backed FreeEOS HDF5 supplements

The native biquintic backend reads the selected distributed `.data` file once
per table load. It reconstructs the six lower Helmholtz constraints from the
printed nodal Pgas, E, S, Cv, chiRho and chiT. The matching `.bq.h5` contributes
only three higher mixed derivatives, the accepted-node mask and numerical
metadata. The interpolation kernel continues to use nine constraints in memory.
There are no new FreeEOS physics calculations and no legacy interpolated values
used to reconstruct these six constraints.

For normal use, hydrate the checkout's Git LFS data, run the standard MESA
`./install`, and set `FreeEOS_use_biquintic_interpolation = .true.` in the EOS
controls. Installation overlays the supplied supplements beside the original
FreeEOS tables. No preprocessing command or directory override is required.
The producer commands below are for explicit regeneration and development.

`freeeos_source.f90` parses printed Q blocks, audits their rectangularity and
uniform Q/T coordinates, and retains the canonical printed axes. Two original
header anomalies are accepted only for their known filenames and declared values:
`02z60x` has a redundant Q maximum of 3.69; `04z50x` declares 445 Q points and a
maximum of 3.23. Both actually contain 527 points ending at 5.69. Original rows
remain untouched, including the `02z70x` block at zero-based Q index 520 whose
printed density is low by 1e-4. That block is authenticated and excluded.
Other printed-density discrepancies are masked by the established source policy.

The low constraints use the same stable analytic radiation subtraction and
physical-to-native transformation as `freeeos_potential.source_constraints`.
Complete support is a centered radius-two cross applied twice, matching the
mean-route mixed derivative estimator. The loader compares the reconstructed
mask to the HDF5 mask exactly before assembling a table. Invalid nodes have zero
placeholders. Accepted nodes must have finite constraints. `assemble_table` in
`freeeos_table` supplies the narrow API that keeps source and HDF5 I/O separate
from the evaluator. `get_node_constraints` supports reference validation.

`freeeos_supplement.f90` uses MESA FORUM's `hdf5io_t`, `OPEN_FILE_RO`, shape
queries and numerical reads. Python shape `(ntau,nu)` becomes Fortran
`(nu,ntau)`. Required attributes are scalar int32/float64; SHA digests are
32-element int32 datasets. Shape and rank checks precede typed reads. Recoverable
schema, identity, mask and finite-value errors explicitly close the handle.
FORUM does not expose a general recoverable I/O or datatype-validation API:
corrupt containers and unsupported type conversions can terminate inside FORUM.
The producer and external verifier check exact on-disk types and lossless filters.

The runtime hashes the selected source file and supplement with the existing
native Fortran SHA-256 module **at load time**. It compares the supplement digest
to `supplements.manifest`, the numeric source digest to the manifest, the source
file digest to that metadata, and the preprocessing digest to the declared
implementation policy. It also checks composition, family, grid and radiation
metadata. There are no external hash commands and no hash work per EOS query.
Hashes are integrity checks against the distributed manifest, not a signature
or a guarantee about a file modified after loading. Tables are cached immutably
by composition and suffix; a different source/supplement directory is an explicit
error in the same process. The legacy X/Z-only cache is not used for reconstruction.

The default native directory is `mesa_data_dir/eosFreeEOS_data/`. The existing
`FreeEOS_biquintic_table_dir` override remains supported. Sources always come
from the MESA data directory at the native call site; standalone family calls
may optionally supply `source_directory`, defaulting to their payload directory.
Existing argument lists remain valid. Explicit old `family.manifest` V1
directories remain comparison inputs. A new `supplements.manifest` has distinct
magic `FREEEOS_BIQUINTIC_SUPPLEMENTS_V1`; a present malformed supplement manifest
cannot select the V1 reader. Missing or mismatched native data produce explicit
errors without ordinary/C/O substitution. Toggle-off routing and mapped
abundance partials retain their prior implementation.

Run generation using a Python environment with numpy and h5py, on an allocated
compute node. Arguments are portable and no report repository is required:

```sh
python3 produce.py --source-dir /path/to/mesa/data/eosFreeEOS_data \
  --output /path/to/staging/eosFreeEOS_data --workers 8
python3 verify.py --source-dir /path/to/mesa/data/eosFreeEOS_data \
  --supplement-dir /path/to/staging/eosFreeEOS_data
python3 package.py --source-dir /path/to/mesa/data/eosFreeEOS_data \
  --supplement-dir /path/to/staging/eosFreeEOS_data \
  --output /path/to/mesa/eos/eosFreeEOS_biquintic_data.tar.xz \
  --audit-output /path/to/archive_audit.json \
  --extraction-dir /path/to/new_extraction_directory --threads 8
```

An optional `produce.py --reference-dir /path/to/V1/family` checks source-generated
constraints and masks against SHA-authenticated archived V1 payloads exactly.
This oracle is unnecessary for collaborator regeneration. Generation writes
float64 mixed arrays with lossless gzip/shuffle/Fletcher32, int8 validity, and
small provenance records. All 89 sources, including six C/O variants, must be
present; `--limit` creates developer fixtures only. Packaging verifies coverage,
uses deterministic tar member ordering and numeric headers, and includes only
89 HDF5 files, 89 small audits and two family metadata files. It never includes
original `.data`, caches, objects, logs or extraction markers. A staged overlay
retains source-file links and revalidates every source and supplement identity.

For private validation, `build_private.sh MESA SDK OUTPUT` builds EOS without
installing libraries; `build_tests.sh MESA SDK OUTPUT` compiles the FORUM fixture,
all-table comparator and native drivers. `test_reader.py` covers an asymmetric
7x5 fixture and recoverable malformed/mismatched cases. `test_tables` compares
all assembled masks and mixed arrays exactly, and reports finite-precision
differences in reconstructed constraints, potential partials, states and state
partials. Runtime reconstruction involves independent Fortran exponentiation
and cancellation in radiation-dominated regimes; exact source-generation/V1
agreement does not imply bitwise runtime reconstruction. Report those differences
separately, without relaxing reference bounds or changing the numerical policy.
