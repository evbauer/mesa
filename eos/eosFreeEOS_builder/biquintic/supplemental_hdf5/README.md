# Legacy-loaded FreeEOS HDF5 supplements

The native biquintic backend obtains nodal Pgas, E, S, Cv, chiRho and chiT from
slot 1 of the existing FreeEOS table loader's `EosDT_XZ_Info%tbl1`. The same
loader reads the same distributed source files and legacy `.bin` caches with
its usual header/version/grid/filename conventions. Slot 1 is nodal physics;
slots 2:4 are bicubic coefficients and are not used to construct the potential.
Only the missing three higher mixed derivatives, complete support mask and
metadata come from the existing `.bq.h5` supplement. The biquintic kernel,
composition weights, analytic derivatives and radiation split are retained.

The family API accepts a real EOS request as optional `legacy_request` at the
end of either argument list. Supplemental-family calls require that context,
including warmed calls; missing context returns `family_legacy_error`.
`eosdt_eval` supplies the real handle request. Family selection passes its ix/iz
to `load_single_eosDT_table_by_id`, whose immutable prefix/suffix identities
keep ordinary and C/O nodal storage independent. Native loading never reparses
or hashes original source bodies. No new native persistent cache or physics
source generation is required.

`freeeos_legacy_nodes.f90` converts the loader's ln(P/E/S) with `math_lib::exp`
and reads Cv/chi directly. Production physical node coordinates use the legacy
header affine log10 T/Q grid. The HDF5/kernel coordinate grid is retained, and
its agreement with the loaded header and actual axes is audited with the
existing 1e-10-dex source/header coordinate tolerance. Printed row coordinates
can differ from their header in the last digits. Independent raw-reader values
therefore need not be bitwise identical: rounding, exp and radiation cancellation
can change values and particularly near-flat response derivatives. TASK085
reports those differences separately from source-cold/binary-warm parity.

HDF5 `valid` remains authoritative for source flags, density-coordinate audits
and stencil support, which are not retained as legacy binary columns. Every
eligible loaded node must also have finite/positive physical columns and finite
constraints. Shapes, metadata, policy and supplement digest are checked through
FORUM. Original source SHA records are compared between manifest and HDF5 as
provenance; original source files are not rescanned. Missing/mismatched data and
unsupported cells remain explicit failures. See [FORMAT.md](FORMAT.md).

For normal use, hydrate Git LFS assets, install MESA normally, and enable
`FreeEOS_use_biquintic_interpolation`. The supplied 89 supplements (83 ordinary,
six C/O) overlay the existing FreeEOS data directory. The native directory
override controls supplements. The real EOS request controls source loading.
The version-51 headers of `mesa-FreeEOS_02z60x.data` and
`mesa-FreeEOS_04z50x.data` report incorrect logQ extents. The latter also
declares 445 Q blocks while its body contains 527. Both bodies extend to
logQ=5.69; the headers report 3.69 and 3.23, respectively. This refactor
deliberately leaves those headers and the standard loader's handling unchanged.
The adapter follows `Locate_logQ`: the loaded count, origin and spacing
define the uniform grid, rather than the erroneous reported maximum. A
shorter loaded table uses the aligned prefix of the authenticated supplement;
the 04z50x native table therefore has 445 Q nodes, while 02z60x retains 527.
Neither table is rejected as a whole for these header defects. Origin,
spacing, temperature-grid, interior-axis and physical-column checks remain.
No corrected source files, repaired-header binary caches, silent cropping or
fallback are distributed by this change.

## Private validation

Run computational commands on the inherited allocated node. Supply a writable
fresh output path; do not overwrite frozen evidence or installed libraries.
`build_private.sh MESA SDK OUTPUT` builds EOS against that tree's already
installed compatible dependencies. `build_tests.sh MESA SDK OUTPUT` uses its
private `OUTPUT/native_build` package tree. TASK085's isolated build script
instead uses the coherent TASK081 SDK26.6.1 dependencies and current EOS sources.

With a private EOS build, from this directory:

```sh
bash build_tests.sh /path/to/mesa /path/to/sdk /path/to/private-output
export MESA_CACHES_DIR=/path/to/private-cache
export MESA_TEMP_CACHES_DIR=/path/to/private-temp-cache
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
/path/to/private-output/tests/test_legacy_tables /path/to/mesa /path/to/eos-inlist \
  /path/to/mesa/eos/eosFreeEOS_builder/mass_frac_solar.txt /path/to/supplements
/path/to/private-output/tests/test_legacy_family /path/to/mesa /path/to/eos-inlist \
  /path/to/mesa/eos/eosFreeEOS_builder/mass_frac_solar.txt /path/to/supplements
```

`test_legacy_tables` attempts all 89 delivered tables through the actual standard
loader and requires matching shapes, finite nodes, zero invalid placeholders,
support and shutdown. The earlier all-89-table
pass and eleven-case paired stellar campaign included a temporary loader
repair; those results do not validate this version without that repair.
`test_legacy_family` uses real handles for 16-thread
simultaneous ordinary/C/O first use, 320 exact warmed-return comparisons,
linear/cubic composition stencils and context/support failures. Budget those
threads within the one allocation.

`test_reader.py --driver OUTPUT/tests/test_reader --output FRESH --sdk SDK`
checks an asymmetric FORUM fixture and malformed metadata. Historical
`test_loader`, `test_loader.py` and `test_tables` deliberately link only the
explicit `historical_source_test.f90` + `freeeos_source.f90` comparison utility.
They exercise the archived independent source-reader/V1 contract, including
source hashing; they are not production loader tests. Those modules are omitted
from `eos/Makefile` and never restore a raw-source fallback. The historical
all-table comparator takes source directory, supplement directory and V1
reference directory; its comparison numbers are not an oracle for bitwise
legacy-nodal arithmetic. The producer/verify/package Python utilities remain
explicit development tools. No producer or archive repack was run for TASK085.
