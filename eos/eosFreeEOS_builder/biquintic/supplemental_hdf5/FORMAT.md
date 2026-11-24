# FreeEOS supplemental HDF5 tables

Implemented format version 1 (TASK078), with the TASK085 legacy-loader adapter.
The default native runtime obtains nodal physics through the existing FreeEOS
loader and its unchanged legacy binary cache, plus these supplements. Explicit
`FREEEOS_BIQUINTIC_V1` text directories remain comparison inputs.

## Purpose and file identity

Keep the distributed FreeEOS `.data` files authoritative for the six nodal
potential constraints through second order. Store only the three estimated
higher mixed derivatives, their accepted-node mask and identity metadata here.

Use one supplemental file per source table, for example:

```text
mesa-FreeEOS_02z70x.data       -> mesa-FreeEOS_02z70x.bq.h5
mesa-FreeEOS_80z00x_CO.data    -> mesa-FreeEOS_80z00x_CO.bq.h5
```

Ordinary and C/O files remain distinct even at identical X/Z. Selection follows
the EOS handle's configured suffix at each participating Z layer. The default
payload directory is `mesa_data_dir/eosFreeEOS_data`; the existing native
directory override remains supported. The native call site supplies its real EOS request to the family API. That
request and the selected composition indices go to the existing
`load_single_eosDT_table_by_id` routine. The source-directory argument is retained
for API/context compatibility; it does not authorize an independent source reader.

## Numerical datasets

Use natural-log coordinates `u=ln(Q)` and `tau=ln(T)`, with physical density
coordinate `ln(rho)=u+2*tau-12*ln(10)`. All three datasets contain derivatives
of the same reference-subtracted specific Helmholtz potential, in erg/g:

| Dataset | Quantity | Stored type |
| --- | --- | --- |
| `f_uu_tau` | d^3 f / du^2 d tau | IEEE float64 |
| `f_u_tau_tau` | d^3 f / du d tau^2 | IEEE float64 |
| `f_uu_tau_tau` | d^4 f / du^2 d tau^2 | IEEE float64 |
| `valid` | Complete accepted nodal jet, 0 or 1 | int8 |

In h5py, write arrays with shape `(ntau, nu)`, temperature as the slow axis.
FORUM receives the corresponding Fortran shape `(nu, ntau)`, addressed as
`array(iu,itau)`. Read `valid` into default integers, validate values 0/1, then
convert to the runtime mask; FORUM need not expose an int8 overload.

Use zero derivative placeholders at invalid nodes. The mask is authoritative
and must agree with the source-floor, coordinate-audit and stencil policy.
Lossless chunked compression and shuffle are permitted; float precision must
not be reduced. Compression does not change the reader's numerical contract.

Do not duplicate P, e, s, Cv, chiRho, chiT or the six reconstructed constraints.
Dense arrays do not need repeated per-node index columns.

## Required metadata

Write scalar numerical root attributes using int32 or float64, matching FORUM
default integer and real(dp) reads:

| Attribute | Type | Meaning |
| --- | --- | --- |
| `format_version` | int32 | 1 |
| `source_version` | int32 | Distributed source format identifier |
| `X`, `Z` | float64 | Source composition |
| `family_id` | int32 | 0 ordinary, 1 C/O |
| `nu`, `ntau` | int32 | Grid dimensions |
| `u_min`, `du`, `tau_min`, `dtau` | float64 | Native grid origin and spacing |
| `radiation_split` | int32 | 1 for the current analytic split |
| `radiation_constant` | float64 | Constant used for the split and restoration |
| `stencil_radius` | int32 | 2 for the current centered estimator |
| `twist_route` | int32 | 1 for the current mean of the two derivative routes |
| `minimum_log10_T` | float64 | 3.2 for the current source-floor policy |
| `preprocessing_version` | int32 | Version of the complete numerical policy |

Store `source_sha256` and `preprocessing_sha256` as rank-one int32 datasets
of length 32, containing the digest bytes as integers in 0..255. These are
small identity records, not new thermodynamic columns. This avoids relying
on character-attribute conversion in the installed FORUM/HDF5 stack.
Optional human-readable provenance strings may be written for inspection;
the runtime contract must not depend on reading them through FORUM.

The preprocessing digest for the implemented policy is
`40b513e48b512e44ea2c020cd29c000b60128661260f38c0c0b7365c2ae10a7d`.
It identifies the complete existing `freeeos_potential.py` implementation.
Changing that implementation requires an explicit policy update, regenerated
supplements and corresponding reader validation, even if version 1 is retained.

## Family manifest and runtime integrity

Use `supplements.manifest` with the distinct magic
`FREEEOS_BIQUINTIC_SUPPLEMENTS_V1`, followed by an entry count and records:

```text
X Z available supplement_filename source_sha256 supplement_sha256
```

The delivered family contains 89 available entries: 83 ordinary and six C/O.
An entry's suffix is part of its identity. An unavailable or missing selected
variant returns an error and cannot select its ordinary counterpart.
Runtime supplement digests are calculated by the existing Fortran SHA-256
module once at table load and compared against the manifest. Original-source
digests remain distributed provenance records and are compared between manifest
and HDF5 metadata; the production native path does not rescan source files.
Original-source loading and cache acceptance follow the existing legacy header,
version, grid, filename and request conventions. No native disk cache is added.
Digest records do not provide cryptographic signing or detect changes after an
immutable cached table has been loaded.

## Reader and source reconstruction

Read through `forum_m`, using `hdf5io_t(filename, OPEN_FILE_RO)` and its
attribute/dataset accessors. Check required names and shapes before typed
reads; FORUM rank/shape assertions can terminate the process. Bound dimensions
before allocating. Validate grid, composition/family, numerical policy and
distributed source provenance before accepting the data. Close the HDF5 handle on completed
reads and on recoverable validation errors.

At table-load time, `freeeos_family` passes the real `EoS_General_Info` request
and selected ix/iz to `load_single_eosDT_table_by_id`. The returned
`EosDT_XZ_Info%tbl1` has shape `(sz_per_eos_point,nv,nQ,nT)`: slot 1 contains
untouched nodal physics and slots 2:4 contain legacy bicubic coefficients.
`freeeos_legacy_nodes` reads only slot 1. Use `math_lib::exp` on the loader's
ln(Pgas), ln(E) and ln(S), which were formed by multiplying its parsed log10
values by ln10. Cv, chiRho and chiT are direct nodal columns. Apply the declared
radiation split and coordinate transform, then insert the three read datasets.
Retain nine constraints in memory for the existing biquintic evaluator.

The loader's decimal conversion/ln transform and exp convention need not
produce bitwise equality with the historical independent raw reader. Quantify
such differences rather than assert old-reader parity. Source-cold and
binary-cache-warm uses of the same legacy loader must agree.

Require version 51, matching Q origin/spacing and temperature grid, aligned
interior Q axes, and finite/positive physical columns at HDF5-eligible nodes.
The runtime Q count is the smaller of the loaded and supplemental counts;
all lower and higher constraints and masks use that same leading region.
The legacy locator defines Q from count/origin/spacing, so the adapter does
not treat a wrong reported Q maximum as a whole-table mismatch.
The distributed `valid` mask remains authoritative for original source flags,
the density-coordinate audit and derivative-stencil support: these flags are
not columns of the legacy binary cache. Do not infer their absence from finite
cached values, widen support or silently crop a mismatched supplement.

A supplemental family call without a real request returns `family_legacy_error`,
including warm calls. Request suffix/prefix identity must match the selected
family. The standard legacy loader must resolve ordinary/C/O identities
correctly; an X/Z-only global cache is insufficient. The two known malformed
version-51 source headers (`02z60x` and `04z50x`) are documented, with no repair
in this refactor. A shorter legacy grid uses the supplement's aligned prefix,
without extending the loaded table or changing its reported Q maximum. The standard
loader and binary-cache header conventions remain in force; no corrected cache
is supplied as a workaround. No production independent ASCII
parser or original-source SHA scan remains. `freeeos_source.f90` is retained as
an off-production historical comparison utility and is omitted from eos/Makefile.

Malformed HDF5 types or a damaged container can still trigger FORUM's own
fatal error policy; preflight checks do not promise general recoverable I/O.
Missing files, mismatched metadata and invalid support remain explicit
errors without substitution of a different family.

## Validation before use

First round-trip a nonsquare, asymmetric analytic fixture through Python and
FORUM to check orientation, types, mask and metadata. Compare assembled nodal
jets and returned physical fields/partials against the current complete text
payloads, including ordinary/C/O selection and multiple handles. Preserve
the archived text files as a numerical reference. Then run bounded stellar
toggle comparisons and their declared restart checks on a compute allocation.
