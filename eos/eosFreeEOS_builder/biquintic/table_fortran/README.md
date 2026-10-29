# Experimental standalone FreeEOS table adapter

This is an **EXPERIMENTAL standalone adapter**, not a MESA-integrated EOS
backend. It reads only the explicit text contract below, reuses the existing
biquintic cell kernel, and derives all thermodynamics from the resulting total
specific Helmholtz potential. It does not edit or invoke runtime EOS routing,
defaults, legacy caches, interpolation fallbacks, clipping, or source-table
preprocessing. Source accuracy and table conventions remain unresolved.

## Build and validation

Use the compiler compatible with the **installed** MESA `const`, `utils`,
`math`, and `auto_diff` modules/libraries. In this workspace that is the supplied
gfortran-15 SDK, not system gfortran 10. From the workspace root:

```bash
export MESASDK_ROOT="$PWD/mesasdk-x86_64-linux-26.6.1"
export PATH="$MESASDK_ROOT/bin:$PATH"
cd mesa/eos/eosFreeEOS_builder/biquintic/table_fortran
make -j2 test
make -j2 BUILD_DIR=build/debug \
  FFLAGS='-O0 -g -std=f2008 -Wall -Wextra -Wno-compare-reals -Wno-unused-dummy-argument -fcheck=all -fbacktrace -finit-real=snan -ffp-contract=off -fopenmp' test
```

`MESA_DIR`, `FC`, `FFLAGS`, and `BUILD_DIR` may be overridden on the make command
line. Use a fresh build directory or `make -B` when changing flags. HELM's
`eos/private/helm_polynomials.f90` and
`../test_fortran/freeeos_biquintic.f90` are read-only build inputs; all newly
compiled objects/modules, binaries, fixtures, and logs stay in this directory's
`build/` tree. No dependency rebuild or alteration of the earlier cell-test
binaries occurs. `make test` uses only Python's standard library for protocol
tests; it does not duplicate the lead's dense real-table comparison.

## Text contract

There are six header records, then exactly `nq*nt` node records:

```text
FREEEOS_BIQUINTIC_V1
nq nt
X Z crad split_integer
source_sha256
preprocessing_sha256 radius route
u_min du tau_min dtau
iq it valid f f_u f_uu f_tau f_u_tau f_uu_tau f_tau_tau f_u_tau_tau f_uu_tau_tau
...
```

The final node derivative is `(2,2)`, namely `d^4 f / du^2 dtau^2`.
Nine values are read in Fortran `nodal(0:2,0:2)` order: **u derivative order
fastest**, then tau order. They are natural-log-coordinate derivatives of a
specific potential in `erg/g`, not derivatives with respect to base-10 logs,
linear temperature/density, or cell fractions. Spacing factors are supplied by
the kernel, never pre-applied to the cache values.

Nodes are ordered with tau outermost and u innermost. `iq,it` are 1-based and
must match their exact sequential positions. `valid` is integer 0 or 1; every
invalid record must contain nine finite zeros. Masks already represent complete
node jets. A cell is usable only if **all four corners** are valid.

Both hashes must be exactly 64 hexadecimal characters, case-insensitive on
read. `radius` is integer 1, 2, or 3; `route` is exactly `mean`, `first`, or
`second`. Provenance is retained in the loaded metadata; no JSON companion is
needed or trusted instead of this header. The hashes are declarations, not
checksums of the cache contents: supplying the optional expected source hash
authenticates only its equality to the header, not the underlying source bytes.

`nq,nt` must each be in `[2,100000]`, with at most 2,000,000 total nodes.
`X,Z,crad`, grid headers, and all node values must be finite; `X,Z >= 0`,
`X+Z <= 1`, and `crad >= 0`. Only split flags 0 and 1 are accepted.
Both spacings must be positive and give finite, representably distinct axes.
Each actual interval must also satisfy the kernel's safe range
`8*huge(dp)^(-1/3) <= spacing <= sqrt(huge(dp))/2`, and exceed
`8*epsilon(dp)*max(abs(lower),abs(upper),tiny(dp))`.

Records are whitespace-separated, at most 4095 characters; individual tokens
are at most 128 characters. Decimal real tokens may use `E/e/D/d` exponents.
List-directed repetition, commas, slashes, comments, missing fields, extra
fields, extra records (including blank trailing records), malformed indices,
NaN/Inf cache tokens, and overflowed real literals are rejected. No old V1
without the provenance line or any legacy MESA cache is silently adopted.

### Declared reference split, not a pure material EOS

For `split=1`, the node potential has had the **declared analytic reference**
`F_rad = -crad*exp(4*t-x)/3` subtracted, and this adapter restores precisely that
same reference, including all derivatives. This does **not** establish that the
remaining table is a pure material EOS. The lead's source audit suggests an
unverified additional radiation-like remainder around +4.40684 ppm; this
adapter neither removes that remainder nor changes the pressure-strip
reconstruction or radiation constant. For `split=0`, nodes contain the total
potential and no analytic potential is added. In both cases `Pgas` is the
pressure derived from the **total** interpolated potential minus exactly
`crad*exp(4*t)/3` from the header. There is no compiled-in radiation default.

## Fortran API

```fortran
use const_def, only: dp
use auto_diff_real_2var_order3_module
use freeeos_table

type(biquintic_table) :: table
type(table_metadata) :: metadata
type(potential_partials) :: partials
type(thermodynamic_state) :: state
type(auto_diff_real_2var_order3) :: x_ad, t_ad, potential_ad
real(dp) :: x, t
integer :: status

call load_table(path, expected_x, expected_z, table, status, expected_sha)
metadata = get_metadata(table)
call evaluate_table_real(table, x, t, partials, status)
x_ad = x
x_ad%d1val2 = 1.0_dp
t_ad = t
t_ad%d1val1 = 1.0_dp
call evaluate_table_ad(table, x_ad, t_ad, potential_ad, status)
call derive_state(partials, x, t, metadata%crad, state, status)
```

Check each `status` before using its results. Omit `expected_sha` to disable
only the expected-source-hash equality check; the source header hash is still
mandatory and validated. Expected X and Z each match to an absolute tolerance
of `32*epsilon(dp)` (about 7.1e-15); there is no composition interpolation.
Failed loading clears the destination, including any earlier loaded table.
Allocation/open failure returns status 10. Metadata includes dimensions,
composition, `crad`, split, axes, both hashes, radius, and route; table arrays are
private. Evaluations perform no I/O or mutable caching and can share a loaded
read-only table between callers. Do not reload it concurrently.

Inputs are **physical natural logarithms** `x=ln(rho)`, `t=ln(T)` in cgs.
The adapter converts `u=x-2*t+12*ln(10)`, `tau=t`. Returned real fields and
`partial_vector` order are:

```text
f f_x f_t f_xx f_xt f_tt f_xxx f_xxt f_xtt f_ttt
```

The real transformation uses `D_x=D_u`, `D_t=D_tau-2*D_u`, through total
degree three. Restored radiation derivatives are
`D_x^i D_t^j F_rad = (-1)^i 4^j F_rad`.
The AD path independently transforms the caller's physical AD inputs, directly
evaluates the polynomial, and independently adds the exponential radiation.
With the seeds above, `ad_vector(potential_ad)` follows the physical real-field
order. Arbitrary seeds, including nonlinear second/third input derivatives,
are supported; their outputs are derivatives with respect to those seeds, not
necessarily x and t. AD variable 2 is x and variable 1 is t for the canonical
physical seeds.

### Strict bounds and face convention

Temperature bounds are exact, inclusive, and reject even a next-representable
out-of-range t. The computed native u has a narrowly bounded representational
endpoint reconciliation, needed because valid physical endpoints can round
outside after the shear conversion. Only a computed u outside an outer
endpoint by at most

```text
8*epsilon(dp)*max(1,abs(x),abs(2*t),12*ln(10),abs(u_min),abs(u_max))
```

is reset to that endpoint. No other out-of-range coordinate is accepted. The
same scalar reset is used for AD without altering its derivative seeds;
analytic radiation still uses the original physical x,t. This is not an
extrapolation region or a general query clamp. At an exact interior face,
search selects the cell to the right/above; at the outer upper endpoint it
selects the last cell. No search for an alternative valid cell occurs if this
selected cell is masked. Generic tables are C2, not C3: thermodynamic first
partials involving third potential derivatives may have one-sided jumps.

## Thermodynamics and numerical state partials

All closure quantities use the same total potential, with `rho=exp(x)`,
`T=exp(t)`, and the following definitions:

```text
P       = rho*f_x
Pgas    = P - crad*T^4/3
E       = f - f_t
S       = -f_t/T
Cv      = (f_t - f_tt)/T
chiRho  = 1 + f_xx/f_x
chiT    = f_xt/f_x
gm1     = f_xt/(T*Cv)
gamma3  = 1 + gm1
gamma1  = chiRho + chiT*gm1
Cp      = Cv + f_x*chiT^2/(T*chiRho)
grad_ad = gm1/gamma1
```

`thermodynamic_state%value`, `%dlnrho`, and `%dlnT` each have 11 entries in
`state_names` order:
`P Pgas E S Cv Cp chiRho chiT gamma1 gamma3 grad_ad`.
Public `state_*` integer indices avoid hard-coded offsets. `%dlnrho` and
`%dlnT` are ordinary numerical `d(value)/dx` and `d(value)/dt`, **not**
derivatives of `ln(value)` and not base-10 derivatives. They are analytically
differentiated closure expressions using the third partials of the same
interpolated potential, not finite-difference runtime queries. Tests compare
them to numerical central differences away from cell faces.

Pressure is `dyn/cm^2`, E and f are `erg/g`, S/Cv/Cp are `erg/g/K`, and the
remaining quantities are dimensionless. This deliberately strict adapter
requires every state value above, as well as `gm1`, to be positive and finite.
In particular, nonpositive E/S or chiT are rejected rather than repaired by
an entropy convention, clipping, or a legacy EOS. These requirements are an
experimental validation contract, not evidence that the unresolved source
conventions satisfy it everywhere.

| Status | Symbol | Meaning |
| --- | --- | --- |
| 0 | `table_ok` | Accepted, finite result |
| 10 | `table_io_error` | Open/allocation failure |
| 11 | `table_format_error` | Invalid header/record/axis/version |
| 12 | `table_mismatch` | Expected composition or source hash mismatch |
| 13 | `table_unloaded` | No successfully loaded table |
| 14 | `table_invalid_input` | Nonfinite query/AD seed/closure input, bad expected hash |
| 15 | `table_out_of_bounds` | Outside the strict grid envelope |
| 16 | `table_invalid_cell` | Selected cell lacks four complete valid corners |
| 17 | `table_arithmetic_failure` | Nonfinite arithmetic or unsafe kernel evaluation |
| 18 | `table_nonphysical` | Nonpositive required closure output |

Failed real/AD/state outputs are entirely quiet NaN; no useful-looking zero
state is returned. Potential evaluation status is independent of closure
status: a finite polynomial can still produce rejected thermodynamics.

## Stable streaming driver

```bash
build/table_query CACHE EXPECTED_X EXPECTED_Z [EXPECTED_SOURCE_SHA256] < queries.txt
```

Each input record contains exactly two whitespace-separated numbers `x t`.
The process loads once, prints and flushes one header, then prints and flushes
one numerical row per query without waiting for EOF. Output has **58 columns**:

1. `x t status_real status_ad status_state`
2. Ten real potential fields in the order above.
3. Ten independently evaluated AD fields with `ad_` prefixes.
4. For each `state_names` entry, its value followed by `<name>_x <name>_t`.

Header spelling includes lowercase `gamma1` and `gamma3` and `P` for total
pressure. Derivative suffixes `_x/_t` mean natural-log coordinate derivatives
of the numerical state value. Rejected queries still produce a row with
statuses and NaN result columns and processing continues. Closure uses the
real result; if real evaluation fails, its status is propagated to
`status_state`. Finite/overflowed decimal query literals are parsed, with
overflowed coordinates reported as status 14. Lexically malformed queries,
including NaN/Inf words, cause exit 2; load rejection causes exit 1 and no
stdout header. Exit 0 means the stream was processed, **not** that every query
was accepted. Empty input and a final record without newline are supported.

## Test coverage and performance

`test_table.f90` creates tiny 3-by-2, two-cell fixtures for an ideal gas with
and without radiation and with total versus reference-split node jets.
It checks analytic physical partials/state values, independent real/AD
agreement, nonlinear AD chain rules through order three, numerical state
partials in both cells, inclusive boundaries, invalid corners, malformed
expected hashes, unloaded/failed reload behavior, NaN outputs, and rejected
nonpositive closure. Interpolation of exponential total node data is not
exact; tolerances distinguish this truncation error from real/AD roundoff.
The protocol suite checks truncation, record shape/order, finite nodes and
headers, bounded/unsafe grids, hashes/radius/routes, rejection of old caches,
256-query batches, and a persistent pipe kept open between queries.

The analytic test also reports a 100,000-query real-potential loop with no
I/O; an optimized checked build here takes about 0.027 CPU seconds. That is
only a tiny-table smoke benchmark, excludes loading/AD/closure, and is not a
claim about real-table dense-audit throughput or EOS accuracy.
