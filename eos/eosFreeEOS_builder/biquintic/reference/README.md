# Direct FreeEOS reference feasibility

## Current scope: existing tables only

The user stopped download/build pursuit on 2026-10-28. **The compile-only
direct-call scaffold is deferred**, not an active dependency acquisition task.
No solver files were obtained and no download process remains active.
`TABLE_ONLY_VALIDATION.md` specifies printed-node precision, iq520 exclusion,
current-MESA radiation reconstruction, gauge metadata limitations and
leakage-free archived-field holdout recommendations for the lead.
`SEARCH_EVIDENCE.md` records saved local/network evidence and its limits.
No table-only experiment establishes fresh-source accuracy. The dependency
checklist below is retained for a separately authorized future task only.

`RADIATION_AUDIT.md` documents the executable table-implied audit and its pinned
input. It finds current-CRAD pressure stripping consistent at printed precision,
while high-temperature energy/linear-density pressure proxies retain a distinct
~4.4 ppm radiation-like excess. Neither result verifies upstream constants.

## Status and scope (2026-10-28)

**No usable direct FreeEOS source/library has been located in the searched
readable local locations. Fresh-source accuracy remains untested.** Installed
MESA tables, MESA's `eosDT_get`/FreeEOS interpolation component, consistency
identities, and table holdouts are not independent executions of FreeEOS.

The external call is reliably established by the existing builder, so this
directory supplies a **compile-checked, unlinked scaffold**, not a validated
FreeEOS executable. No solver release, convergence/error semantics, historical
radiation constant, or absolute energy/entropy reference has been inferred from
the mere presence of tabulations. No dependency archive was downloaded;
resumed HEAD-request failures are recorded in `SEARCH_EVIDENCE.md`.

All new artifacts are under this directory. No builder, installed tables,
caches, preprocessing, `test_fortran`, or report sources are changed. The lead
owns preprocessing/test coordination; the kernel/report agent owns the report.

## Dependency search and evidence

Searches used `find` with case-insensitive `*freeeos*`, `*free_eos*`,
`*free-eos*`, and `*mesasdk*` patterns as applicable, and `rg` for source calls
and build configuration. See `SEARCH_EVIDENCE.md` for commands and limitations.

- Workspace: the builder call, table archives, installed ASCII tables, binary
  caches, reports, and MESA interpolation code are present. No matching solver
  source or separately named solver library was found.
- `/g/g14/bauer39` outside this workspace: no matching dependency names found.
  `/g/g14/bauer39/mesa`, `/g/g14/bauer39/mesasdk`, and
  `/g/g14/bauer39/.local` do not exist.
- Conventional system roots and `/usr/tce/packages` were inspected; inaccessible
  directories and search completion limits are explicitly recorded, not treated
  as evidence of absence. No claim is made about every shared filesystem.
- `MESASDK_ROOT`, `MESA_DIR`, `FREE_EOS_DIR`, and `PKG_CONFIG_PATH` are unset in
  this session. `LD_LIBRARY_PATH` only names the TeX Live directory.
  `pkg-config --modversion free-eos` and `--libs free-eos` both fail because the
  package is not found; `ldconfig -p` has no matching solver entry.
- All `.a` files under `mesa/build` were inspected with `nm --defined-only`:
  no `free_eos`/`free_eos_` text definition. `libeos.a` supplies interpolation
  routines and `free_eos_handle`, not the external thermodynamic solver.
  `mesa/build/free-eos-builder` is absent.
- `../../Makefile:7` (relative to here) declares
  `EXTERNAL_DEPENDS_ON := free-eos`. `mesa/make/setup-depends.mk` resolves external
  dependencies through `pkg-config`. The builder README's older instruction to
  edit `make/makefile` and set `FREE_EOS_DIR` does not describe this checkout's
  actual build layout; that file is absent. Do not edit either build system for
  this standalone reference task.

Evidence base: MESA revision `7b19fe283fc5315a7aaeb02f353e5ee4e611f3ed`.
The parent workspace itself is not a Git repository. Pre-existing plotter and
plan changes are outside this task.

### Input fingerprints (SHA-256)

| Workspace-relative input | SHA-256 |
| --- | --- |
| `mesa/eos/eosFreeEOS_builder/src/free_eos_table.f90` | `d2e3243befa696ddc41a3c108cf850bc93ce0c64d498ac3846618746644c0137` |
| `mesa/eos/eosFreeEOS_builder/inlist_solar` | `d41cb4848c32bbb50feae430bb29ae1fb54d2e86c862ea39923e2707ebc67dc7` |
| `mesa/eos/eosFreeEOS_builder/mass_frac_solar.txt` | `bd62922804e073a2a410480994b8a28e85b0bc53e469729d6ed5136ea9da10a6` |
| `mesa/data/eosFreeEOS_data/mesa-FreeEOS_02z70x.data` | `e6ca7ce27bee8197bb257a8c4886a06366e6951e58958fbed3735c877ee68502` |

## Established release, options, API, and composition

`mesa/docs/source/eos/overview.rst:148` identifies the **table family** as
FreeEOS **2.2.1, EOS1**, with scaled-solar composition. The original report
`freeeos_report/freeeos_interpolation_report.tex:86` cites that documentation.
Neither is proof of an installed direct solver or a complete historical build
manifest. Table-format `version=51` is not a FreeEOS release number.

`inlist_solar` selects `eos_version=1`. Builder
`src/free_eos_table.f90:78` maps this to `ifopt=3, ifmod=1, ifion=-2`.
`kif=2` is declared at line 40 as the density/temperature query mode. These are
the scaffold's exact constants. EOS1a (`ifion=-1`) and EOS4 (`ifmod=101`) are
not interchangeable substitutes. Detailed meanings of individual flags and
any additional internal controls must come from the actual 2.2.1 source/docs;
they are not invented here.

Builder lines 105–129 establish default Fortran integers, double precision,
20 abundances, six length-three real arrays, and this positional external ABI:

```text
free_eos(ifopt, ifmod, ifion, kif, eps, Neps,
         lnRho_in, lnT_in, logf, T, Rho, lnRho_out, P, logP,
         Cf, Cp, Qf, Qp, Sf, St, grad_ad, RTP, Rmue, fh2,
         fhe2, fhe3, xmu1, xmu3, eta, gamma1, gamma2, gamma3,
         h2rat, h2plusrat, lambda, gamma_e, degeneracy,
         pressure, density, energy, enthalpy, entropy, iter)
```

The builder uses an implicit external interface and passes the same `lnRho`
variable for input/output; the scaffold preserves that call pattern rather
than guessing dummy-argument intents. `const_def` uses
`selected_real_kind(p=15)` (real64 here); the scaffold uses `real64`.
Compiler options changing default real/integer kinds would invalidate this
assumption. Output `iter` is retained raw, not assigned an invented success
test. Source-side input/output intent and convergence semantics remain to be
confirmed before trusting any returned row.

The exact mixture is `X=0.70`, `Y=1-X-Z=0.28`, `Z=0.02`.
The metal pattern in `mass_frac_solar.txt` sums to **1.000000 in decimal**;
the builder does not renormalize it. Its scaled mass fractions are:

| Element | Mass fraction | Element | Mass fraction |
| --- | ---: | --- | ---: |
| H | 0.70000000 | He | 0.28000000 |
| C | 0.00343672 | N | 0.00100670 |
| O | 0.00934712 | Ne | 0.00209662 |
| Na | 0.00004180 | Mg | 0.00079848 |
| Al | 0.00007206 | Si | 0.00088114 |
| P | 0.00000846 | S | 0.00047026 |
| Cl | 0.00000584 | Ar | 0.00008670 |
| Ca | 0.00007792 | Ti | 0.00000390 |
| Cr | 0.00002234 | Mn | 0.00001558 |
| Fe | 0.00152866 | Ni | 0.00009970 |

`eps(i)=mass_frac(i)/atom_wgt(i)`, using exactly these builder/FreeEOS weights,
not MESA isotope masses and not a number-fraction normalization:

```text
H He C N O Ne Na Mg Al Si P S Cl Ar Ca Ti Cr Mn Fe Ni
1.007825 4.0026 12.0111 14.0067 15.9994
20.179 22.98977 24.305 26.9815 28.086
30.9738 32.06 35.453 39.948 40.08
47.9 51.996 54.938 55.847 58.71
```

The scaffold reads the existing metal list, verifies order and normalization,
scales it exactly as the builder does, and records both mass fractions and
`eps` in each output. Preserve its checksum for every experiment.

## Radiation, reference gauges, and table fields

The native path at builder lines 131–165 takes `P`, `energy`, `entropy`, and
their derivatives directly from the solver, but stores
`logPgas=log10(P-radiation_pressure(T))`. The original report's radiation
discussion (line 943) interprets these solver outputs as **total** quantities.
There is no radiation stripping from energy/entropy in the builder. This
interpretation still needs a direct 2.2.1 source/constant check; a table-only
identity test cannot establish the solver's conventions.

| Native ASCII field | Builder meaning / convention |
| --- | --- |
| `logT`, `logQ`, `logRho` | Base-10 coordinates; `q=log10rho-2*log10T+12`. Solver inputs are natural logs. |
| `logPgas` | Base-10 gas pressure, cgs `dyn cm^-2`, after MESA radiation subtraction. |
| `logE`, `logS` | Base-10 solver specific energy (`erg g^-1`) and entropy (`erg g^-1 K^-1`); no radiation/gauge adjustment. |
| `chiRho`, `chiT` | `pressure(2:3)`, logarithmic derivatives of solver/total pressure, not gas-pressure susceptibilities. |
| `Cp`, `Cv` | Solver `Cp` and `energy(3)/T`; `erg g^-1 K^-1`, unstripped. |
| `dE_dRho`, `dS_dT`, `dS_dRho` | `energy(2)/rho`, `entropy(3)/T`, `entropy(2)/rho`; physical-coordinate derivatives of unstripped quantities. |
| `gamma1`, `gamma3`, `grad_ad` | Solver outputs, unstripped. |
| `mu`, `log10_free_e`, `eta` | `1/xmu1`, `log(xmu3)/ln10`, `eta`. No separate conversion using isotope weights. |
| `MESA` | Native/fallback flag; 0 marks the accepted native path, 1 the fallback path. |
| `dpe`, `dsp`, `dse` | Redundant identities computed from solver outputs; not source-accuracy errors. |

`energy(2:3)` and `entropy(2:3)` are derivatives with respect to **natural**
logarithms of density and temperature, as demonstrated by the builder's
physical derivative conversions. Array entries for degeneracy, density,
enthalpy, and `pressure(1)` are preserved raw by the scaffold; their precise
semantics are not assumed from names alone.

Current MESA `radiation_pressure(T)` is `crad*T*T*T*T/3`
(`mesa/eos/public/eos_lib.f90:88`). `crad` is derived from CODATA-2018
`k=1.380649e-16`, `h=6.62607015e-27`, `c=2.99792458e10`, and `pi`
in `mesa/const/public/const_def.f90:71`; evaluating that expression gives
approximately **7.56573325027999926e-15 erg cm^-3 K^-4**.
This value is **not established as the historical table-generation value or
FreeEOS 2.2.1's radiation constant**. The scaffold therefore exports raw
solver pressure and does not silently subtract today's MESA radiation EOS.

For a verified common constant `a`, analytic radiation is
`P_rad=a*T^4/3`, `e_rad=a*T^4/rho`, `s_rad=4*a*T^3/(3*rho)`,
`f_rad=-a*T^4/(3*rho)`, `Cv_rad=4*a*T^3/rho`.
Total-pressure reconstruction from ASCII uses the historical **builder**
subtraction constant; consistent material subtraction from solver outputs
uses the verified **solver** convention. These need not be presumed equal.
For material pressure, transform susceptibilities as
`chiRho_gas=P_total*chiRho_total/P_gas` and
`chiT_gas=(P_total*chiT_total-4*P_rad)/P_gas`, with the corresponding radiation
subtractions in energy/entropy and all derivatives. Recompute `Cp`, adiabatic
coefficients, and potential derivatives consistently; do not subtract a
guessed scalar `Cp_rad` or mix gas pressure with total energy/entropy.

Absolute energy/entropy zeros are carried through unchanged by the builder.
The actual solver gauges, binding-energy conventions, rest-mass treatment,
and constants are unavailable locally and must be audited in the release.
An energy offset and an entropy offset shift `f=e-T*s`; do not fit away
disagreements and label the fit independent accuracy. Record any verified
constant-gauge conversion separately and apply it coherently to the potential.

Fallback rows are not direct FreeEOS references. The original report (line
963) notes their zero diagnostics and a native/base-10 versus fallback/natural-log
conversion concern. Do not use them to calibrate gauges or radiation. Do not
repair or regenerate tables in this task.

## Deferred dependency and build/run checklist

**Not authorized in the current table-only scope.** Do not resume downloads,
dependency builds or direct queries from this checklist without a new scope
instruction. No table generation is required or permitted by this harness.

1. Supply a readable **FreeEOS 2.2.1 source tree and matching built library**,
   or request approval to acquire/build one under `reference/`. Network access,
   downloads, installing software elsewhere, and changing shared artifacts
   require approval. A table archive is insufficient.
2. Inspect the release/version evidence, README/build instructions, actual
   `free_eos` declaration, EOS1 control mapping, radiation inclusion and
   constants, energy/entropy gauges, convergence/`iter` behavior, and any runtime
   data files. Record hashes, build options, compiler, full link line, and
   runtime paths under `reference/`. Identify any documented initialization;
   the existing builder calls none before the external solver, but this is not
   permission to invent or omit a release-required initialization routine.
3. Follow that release's build instructions in a private `reference/dependency/`
   tree/build directory with a reference-owned install prefix if needed. Exact
   upstream configure/CMake targets and transitive libraries cannot responsibly
   be specified without its source. Do not run MESA `rn_all`, table exporters,
   EOS extraction targets, or the table-builder executable. No default-kind
   promotion flags. Verify `free_eos_` via `nm` for a GNU Fortran build and verify
   ABI compatibility with the established call before execution.
4. Build this scaffold independently of MESA. `make check` only compiles it;
   actual linking requires explicit, verified `FREEEOS_LIBS`:

   ```sh
   cd mesa/eos/eosFreeEOS_builder/biquintic/reference
   make check
   pkg-config --modversion free-eos
   pkg-config --libs free-eos
   make FREEEOS_LIBS="$(pkg-config --libs free-eos)"
   ```

   The last three commands are conditional on a verified package being made
   available. Otherwise use `make FREEEOS_LIBS='/absolute/path/to/verified/archive
   plus its documented transitive link flags'` as one shell argument, not a
   guessed `-lfreeeos` library name. Keep runtime search paths private. If flags
   or dependency builds change, force relinking with `make -B` and record them.
5. After the API/conventions audit and link, run from **this directory** only:

   ```sh
   mkdir -p runs
   ./build/direct_freeeos ../../mass_frac_solar.txt queries_example.txt smoke_eos1.dat
   ```

   Output is a new `runs/smoke_eos1.dat`; existing files are never overwritten.
   The illustrative coordinates are not certified native nodes and do not
   constitute validation. Raw nonfinite/sentinel results are not replaced by
   MESA/HELM, logged, clipped, or declared successful merely because a call
   returns. Establish source-defined success and physical validity first.
6. Obtain query lists from the lead's audited native masks and saved report
   samples. For **matched nodes**, use native-only, positive finite rows and
   initially `log10T>=3.2`, with complete clean stencils. Retain node IDs, flags,
   original printed precision, and requested versus returned density. Use the
   actual printed temperature and q block, with
   `log10rho=q+2*log10T-12`. The builder accumulates q upward and T downward;
   floating-point accumulation/printed rounding must not be mistaken for
   exact-coordinate source errors. The ASCII `logRho` field records the requested
   density, not necessarily the returned solver density. Save these two-column
   coordinates under `reference/` and call the driver to a dedicated
   `runs/matched_nodes_*.dat` basename.
7. For **off-node** direct accuracy, retain the report/lead's exact off-node
   physical coordinates and valid native-stencil selection, not new unrelated
   random samples. Query those identical `(rho,T)` through the direct driver
   and each interpolator, writing separate `runs/off_nodes_*.dat` outputs.
   Compare pressure, energy, entropy, heat capacities, and derivatives under
   verified common radiation/gauge conventions. H/He ionization and molecular
   regions need separate coverage; high-temperature radiation-dominated tests
   wait until the constant audit is complete. Do not count decimated-table
   holdouts, self-consistency identities, or agreement with MESA's interpolated
   `eosDT_get` as fresh-source accuracy.

## Validation performed here

- GNU Fortran 10.3.1: `make check` succeeds with `-O0 -g -Wall -Wextra
  -fcheck=all`; no compiler diagnostics. See `logs/compile_check.log`.
- `gfortran build/direct_freeeos.o -o build/dependency_probe` fails, exit 1,
  with **undefined reference to `free_eos_`**. See `logs/dependency_probe.log`.
  This demonstrates the scaffold's unresolved external dependency; by itself
  it does not prove no library exists anywhere.
- Default `make` refuses to link without explicit `FREEEOS_LIBS` (make exit 2).
  See `logs/dependency_guard.log`.
- **No direct calls, matched-node comparisons, or off-node source-accuracy
  tests ran.** The remaining dependency is a provenance-verified 2.2.1 solver
  and its convention/ABI audit, not a regenerated table or cache.
