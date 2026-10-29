# Standalone FreeEOS biquintic cell kernel

This directory evaluates one tensor-product quintic Hermite cell at a time in
`u = ln(Q)` and `tau = ln(T)`, without loading tables or invoking MESA's runtime
EOS. It compiles the existing `eos/private/helm_polynomials.f90` into a local
object and follows the standalone WDM_AV test harness's compiler/link pattern.
Only WDM_AV's `qwts` spacing/sign convention and `fi5` corner packing are reused;
there is no isotope rescaling, mixing, or other physics here.

## Build and test

Use a compiler compatible with the **existing** MESA module files and libraries
in `MESA_DIR/build`. Required packages are `const`, `utils`, `math`, and
`auto_diff`, plus their existing external link dependencies (CRlibm in this
workspace). `pkg-config` supplies the flags; the private order-3 AD module path
is added explicitly. No target builds, installs, or changes MESA dependencies.
`MESASDK_ROOT`, when set, also supplies an SDK library search path for relocated
SDK installations. The system gfortran 10 here is not compatible with the
existing gfortran-15 MESA module files; use the supplied SDK.

From the workspace root:

```bash
export MESASDK_ROOT="$PWD/mesasdk-x86_64-linux-26.6.1"
export PATH="$MESASDK_ROOT/bin:$PATH"
cd mesa/eos/eosFreeEOS_builder/biquintic/test_fortran
mkdir -p build/optimized build/debug
make -B -j2 test > build/optimized/test.log 2>&1
make -B -j2 BUILD_DIR=build/debug BIN_DIR=build/debug \
  FFLAGS='-O0 -g -std=f2008 -Wall -Wextra -Wno-compare-reals -Wno-unused-dummy-argument -fcheck=all -fbacktrace -finit-real=snan -ffp-contract=off -fopenmp' \
  test > build/debug/test.log 2>&1
```

The default optimized binaries are `./cell_driver` and `./test_biquintic`, with
objects/modules in `build/optimized`. The independent debug build keeps both
binaries and objects/modules in `build/debug`, leaving the optimized driver
available to the lead. `make test-kernel` and `make test-driver` run the suites
separately. `MESA_DIR`, `FC`, `BUILD_DIR`, `BIN_DIR`, and `FFLAGS` can be overridden
on the make command line. Use a fresh build directory or `make -B` when changing
compiler flags. All generated build products, smoke-test inputs/outputs, and
logs live under the ignored `build/` directory, except the two ignored default
binaries.

## Fortran API

```fortran
use const_def, only: dp
use auto_diff_real_2var_order3_module
use freeeos_biquintic

real(dp) :: nodal(0:2,0:2,2,2), u0, u1, t0, t1, u, tau
type(freeeos_partials) :: partials
type(auto_diff_real_2var_order3) :: u_ad, tau_ad, f_ad
integer :: status

call evaluate_cell_real(nodal, u0, u1, t0, t1, u, tau, partials, status)
u_ad = u
u_ad%d1val2 = 1.0_dp
tau_ad = tau
tau_ad%d1val1 = 1.0_dp
call evaluate_cell_ad(nodal, u0, u1, t0, t1, u_ad, tau_ad, f_ad, status)
```

`nodal(a,b,corner_u,corner_tau)` is the physical
`d**(a+b) f / du**a dtau**b`, for `a,b = 0,1,2`, at the four corners, in
`erg/g`. These are derivatives with respect to **natural-log coordinates**, not
linear `Q`/`T`, base-10 logs, or normalized cell coordinates. Do not multiply
them by grid spacings before calling the kernel. The kernel supplies the
Hermite basis factors `du`, `du**2`, `dtau`, and `dtau**2` itself. Corner indices
1 and 2 denote the corresponding lower and upper endpoints. `t0,t1` are
endpoints in `tau`, not temperatures. Different spacings on each axis and in
adjacent cells are supported. All nine tensor derivatives per corner are used,
including `(2,2)`, although outputs stop at total degree 3.

The real result fields, in the driver output order, are:

| Field | Derivative | AD field with the seeds above |
| --- | --- | --- |
| `f` | value | `val` |
| `f_u` | `(1,0)` | `d1val2` |
| `f_tau` | `(0,1)` | `d1val1` |
| `f_uu` | `(2,0)` | `d2val2` |
| `f_u_tau` | `(1,1)` | `d1val1_d1val2` |
| `f_tau_tau` | `(0,2)` | `d2val1` |
| `f_uuu` | `(3,0)` | `d3val2` |
| `f_uu_tau` | `(2,1)` | `d1val1_d2val2` |
| `f_u_tau_tau` | `(1,2)` | `d2val1_d1val2` |
| `f_tau_tau_tau` | `(0,3)` | `d3val1` |

The real path uses HELM's analytic derivatives and `h5`. The separate AD path
evaluates the potential polynomial directly and returns only `f_ad`, whose
derivative fields follow the caller's seeds, including nonlinear chain rules.
Endpoints and nodal data are held constant in the AD path. The table above
uses AD variable 1 = `tau` and variable 2 = `u`, matching WDM_AV's convention.
Other seeds are valid; their resulting derivatives are not necessarily native
`u,tau` partials.

### Validation and status

Check `status` before consuming outputs. Every returned field is quiet NaN on
failure, in both paths. Validation checks all 36 nodal entries, all endpoints,
both query values, and all AD input fields for finiteness.

| Symbol | Value | Meaning |
| --- | --- | --- |
| `cell_ok` | 0 | Accepted cell and finite outputs |
| `cell_invalid_input` | 1 | Nonfinite nodal, endpoint, query, or AD input |
| `cell_invalid_spacing` | 2 | Nonincreasing, overflowed, or numerically unsafe spacing |
| `cell_out_of_bounds` | 3 | Query outside the closed cell |
| `cell_arithmetic_failure` | 4 | Accepted inputs produce any nonfinite output |

Both spacings must be between `8*huge(1.0_dp)**(-1.0_dp/3.0_dp)` and
`0.5*sqrt(huge(1.0_dp))`, inclusive. This conservatively keeps basis spacing
factors and inverse-cube derivative factors representable; arbitrary finite
nodal magnitudes or AD seeds can still overflow during evaluation. Arithmetic
failure detection assumes IEEE exceptions are not trapped. The standalone
tests deliberately exercise overflow without enabling `-ffpe-trap`.

The domain is exactly `[u0,u1] x [t0,t1]`. There is no extrapolation, tolerance,
or boundary clamping, even for a query one representable number outside.
Coordinate conversions done by the caller must respect this domain; affine
round-trip floating-point arithmetic need not reproduce an endpoint exactly.

### Sheared physical coordinates

For `x = ln(rho) = u + 2*tau - 12*ln(10)`, physical operators are
`d/dx|tau = d/du` and `d/dtau|x = d/dtau|u - 2*d/du`. Thus, for example,
`F_x = f_u`, `F_tau = f_tau - 2*f_u`, and
`F_tau_tau = f_tau_tau - 4*f_u_tau + 4*f_uu`. The tests independently apply
this transformation through total degree 3. For physical AD partials, seed
`x_ad` as variable 2, `tau_ad` as variable 1, and form
`u_ad = x_ad - 2*tau_ad + 12*ln10` before calling the AD kernel. The kernel
itself does not apply this shear or compute thermodynamic response functions.

## Streaming driver format

`cell_driver` reads stdin until EOF and emits **three lines per accepted cell**.
One cell consists of exactly six newline-delimited records:

```text
u0 u1 t0 t1
u tau
f f_u f_uu f_tau f_u_tau f_uu_tau f_tau_tau f_u_tau_tau f_uu_tau_tau
f f_u f_uu f_tau f_u_tau f_uu_tau f_tau_tau f_u_tau_tau f_uu_tau_tau
f f_u f_uu f_tau f_u_tau f_uu_tau f_tau_tau f_u_tau_tau f_uu_tau_tau
f f_u f_uu f_tau f_u_tau f_uu_tau f_tau_tau f_u_tau_tau f_uu_tau_tau
```

The last four records are corners `(u0,t0)`, `(u1,t0)`, `(u0,t1)`, `(u1,t1)`,
respectively. Within each corner, `a` varies fastest in `nodal(a,b,...)`, as
shown. The final `(2,2)` field is an **input only**. Plain decimal numbers and
Fortran `E`/`D` exponent notation are accepted. Separate values with spaces or
tabs; CRLF and a final record without a newline are accepted. Each record must
contain exactly 4, 2, or 9 numeric tokens, as appropriate, and fewer than 4096
characters. Blank records, comments, commas, slash termination, repetition
syntax, and surplus tokens are not part of this format.

Concatenate six-record cells without separators to amortize process startup:

```bash
./cell_driver < constant_cell.txt
cat constant_cell.txt affine_cell.txt | ./cell_driver
./cell_driver < many_cells.txt > outputs.txt
```

Each accepted cell emits:

1. The unchanged header
   `f f_u f_tau f_uu f_u_tau f_tau_tau f_uuu f_uu_tau f_u_tau_tau f_tau_tau_tau`.
2. Ten real-path values in that order, in `ES25.16E3` format.
3. Ten AD-path values in the same order, with the native seeds above.

Stdout is flushed after every complete cell, so a persistent pipe can exchange
one cell at a time without closing stdin. Empty input and EOF between cells
exit successfully. Malformed or incomplete records, including a partial next
header at EOF, exit 2. A kernel rejection exits 1 with the cell index and
numeric kernel status on stderr. Rejected cells emit **no** output lines;
previously accepted cells remain in stdout. Processing stops at the first
error, rather than skipping a rejected cell. There is no startup banner,
summary, or extra stdout text, preserving the single-cell three-line behavior.

## Analytic coverage and limitations

`test_biquintic.f90` exercises all ten output fields:

- All 36 monomials `u**a*tau**b`, `a,b <= 5`, plus a dense `(5,5)` polynomial
  at `erg/g` scale. Probes include interior, near-face, exact-face, and all
  corners on unequal axis spacings.
- Arbitrary independent corner jets and low-order endpoint derivatives,
  including the supplied `(2,1)` and `(1,2)` derivatives; separate real/AD
  agreement for arbitrary nodal data and nonzero coordinate origins.
- Two-exponential analytic convergence on four halved grids, grouped by
  derivative degree. Expected orders are 6, 5, 4, and 3 for degrees 0..3.
- Two-cell faces in both orientations: C2 continuity and tangential/mixed
  third-partial continuity, while explicitly observing allowed jumps in the
  third **normal** derivative.
- Exact reproduction of an independently differentiated polynomial in the
  sheared `x,tau` coordinates, including physical real and AD partials, plus
  nonlinear AD input chain rules through degree 3.
- NaN and both infinities in every scalar, nodal slot, and AD field;
  nonincreasing/unsafe spacings; next-representable out-of-cell queries;
  finite-input arithmetic overflow; NaN poisoning of all failed outputs.

`test_driver.sh` checks constant/affine numerical outputs and field order,
multi-cell streams including a 256-cell batch and a persistent pipe kept open
between queries, whitespace/EOF behavior, every
truncated/partial record, malformed tokens, invalid cells and nonfinite data,
and valid output followed by a rejected or incomplete next cell.

This is a local interpolation kernel, not a table loader, cell search,
preprocessor, thermodynamic closure, MESA runtime integration, or validation
of archived FreeEOS data. Shared cells are C2 only when their complete face
jets agree on the same tangential interval; generic local biquintic Hermite
interpolation is not C3 across faces. Cancellation and roundoff remain,
especially for nearly zero derivatives of large potentials or physical shear
combinations. The tests report explicitly scaled errors rather than claiming
small relative errors when an exact derivative is zero. Debug compilation
covers only these standalone sources; existing MESA dependency libraries are
reused unchanged. Python preprocessing, real-table auditing, and the report
remain with the lead.
