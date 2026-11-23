# FreeEOS supplemental HDF5 tables

Implemented format version 1, 23 November 2026 (TASK078). The default native
runtime reads distributed source tables plus these supplements. Explicit
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
directory override remains supported. The native call site supplies that
original source directory independently of the supplement-directory override.

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
Runtime source and supplement digests are calculated by the native Fortran
SHA-256 module once at table load. They are compared against the manifest and
numeric metadata; this is separate from producer and archive verification.
The digest records provide distributed-file integrity, not cryptographic signing
or detection of changes after an immutable cached table has been loaded.

## Reader and source reconstruction

Read through `forum_m`, using `hdf5io_t(filename, OPEN_FILE_RO)` and its
attribute/dataset accessors. Check required names and shapes before typed
reads; FORUM rank/shape assertions can terminate the process. Bound dimensions
before allocating. Validate grid, composition/family, numerical policy and
source identity before accepting the data. Close the HDF5 handle on completed
reads and on recoverable validation errors.

At table-load time, reconstruct the six lower constraints from the selected
source's nodal values, apply the same radiation subtraction and coordinate
transform, and insert the three read datasets. Retain nine constraints in
memory for the existing interpolation evaluator. Do not reconstruct from
legacy interpolated query outputs or assume a global X/Z-only source cache
has the correct configured suffix.

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
