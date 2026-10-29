# Streaming default-MESA active-fraction probe

This harness evaluates the **actual default MESA EOS blend**, not an isolated
component or a guessed component-domain cut. `use_WDM_AV=.false.` is the sole
inlist override. Initialization follows plotter `Setup_eos`: `const_init`,
`math_init`, `chem_init`, `eos_init`, and handle allocation run once per process.
No MESA/default/runtime/table sources are changed or rebuilt.

## Build and run

From the workspace root:

```sh
make -C mesa/eos/eosFreeEOS_builder/biquintic/mesa_probe
printf '3.53 6.51\n2.33 6.765\n' | python3 -B mesa/eos/eosFreeEOS_builder/biquintic/mesa_probe/run_probe.py
python3 -B mesa/eos/eosFreeEOS_builder/biquintic/mesa_probe/test_smoke.py
```

The matching existing SDK is `mesasdk-x86_64-linux-26.6.1`, GNU Fortran 15.2.0.
The Makefile only builds `mesa_probe/build/mesa_probe`. It obtains installed
module/link flags from pkg-config for `mesa-eos mesa-const mesa-math mesa-chem
mesa-utils`, including their existing transitive dependencies. No MESA target
is invoked. `014_active_domain/build.log` records the complete expanded command.

The equivalent exact compile command, with the same pkg-config resolution:

```sh
PROBE_DIR="$PWD/mesa/eos/eosFreeEOS_builder/biquintic/mesa_probe"
MESA_DIR="$PWD/mesa"
SDK_DIR="$PWD/mesasdk-x86_64-linux-26.6.1"
PKG_CONFIG_PATH="$(printf '%s:' "$MESA_DIR"/build/*/lib/pkgconfig)$SDK_DIR/lib/pkgconfig:${PKG_CONFIG_PATH:-}"
export PKG_CONFIG_PATH
MODULE_FLAGS="$(pkg-config --define-variable=math_slot=crmath --cflags mesa-eos mesa-const mesa-math mesa-chem mesa-utils)"
LINK_FLAGS="$(pkg-config --define-variable=math_slot=crmath --libs --static mesa-eos mesa-const mesa-math mesa-chem mesa-utils)"
mkdir -p "$PROBE_DIR/build"
MESASDK_ROOT="$SDK_DIR" LD_LIBRARY_PATH="$SDK_DIR/lib:$SDK_DIR/lib64:${LD_LIBRARY_PATH:-}" \
"$SDK_DIR/bin/gfortran" -O2 -g -fopenmp -ftrampoline-impl=heap -ffp-contract=off -ffree-line-length-none \
$MODULE_FLAGS "$PROBE_DIR/probe.f90" -o "$PROBE_DIR/build/mesa_probe" \
-L"$SDK_DIR/lib" -L"$SDK_DIR/lib64" -Wl,-rpath,"$SDK_DIR/lib" -Wl,-rpath,"$SDK_DIR/lib64" $LINK_FLAGS
```

Some existing `.pc` private link flags contain the original SDK installation
path. The explicit current SDK library search/runtime paths resolve those
libraries here; no installed `.pc` file is edited.

## Stable stream interface

Launch `run_probe.py` with pipes; it inherits stdin/stdout without buffering the
query batch. The driver emits and flushes one header after initialization, then
one row per nonblank/non-comment input record, **before EOF**. Send exactly two
numbers per line: base-10 density and temperature. Blank lines and lines starting
with `#` are skipped. MESA diagnostics and composition IDs go to stderr, not the
numeric stdout stream.

```text
logRho logT ierr frac_HELM frac_OPAL_SCVH frac_FreeEOS frac_PC frac_Skye frac_CMS frac_ideal frac_WDM_AV Cv dCv_dlnRho dCv_dlnT
```

`ierr=0` means successful MESA call and finite exported fields, with fractions
in [0,1] summing to one within 1e-10. Fractions use the actual MESA `i_frac`
indices. Cv is in erg g^-1 K^-1; its partials are with respect to **natural**
log density/temperature, despite the base-10 input coordinates. WDM_AV must be
disabled in the allocated handle. No fractions are invented for failed queries.

Nonzero MESA statuses are preserved. Harness statuses: -1001 malformed/extra
input tokens; -1002 nonfinite coordinates; -1003 unsafe floating-point exponent
range; -1004 nonfinite basic result or exported Cv partial; -1005 invalid fraction
range/sum. Every failed row returns NaN fractions/Cv/partials and processing
continues. Initialization/composition failures terminate nonzero before the header.
This is a finite-status guard, not a guarantee MESA supports every finite physical
input or every extreme EOS domain.

Wrapper optional arguments: `--mesa-dir`, `--sdk-dir`, `--inlist`, `--massfile`.
The executable itself takes exactly `MESA_DIR INLIST MASSFILE`; using the wrapper
is required for report-cache setup. The default inlist changes no other control.
Supplying another inlist changes the declared experiment; its hash is recorded.

## Composition and cache safety

H1=.70, He4=.28, Z=.02; metal fractions come from the pinned
`mass_frac_solar.txt` SHA-256
`bd62922804e073a2a410480994b8a28e85b0bc53e469729d6ed5136ea9da10a6`.
Its exact decimal sum is 1.000000; binary64 sum here is 0.99999999999999989.
The 18-element pattern is normalized by that sum, then multiplied by .02;
H/He are not rescaled. The final total is checked within 1e-12.

Representative isotopes, in order:
H1, He4, C12, N14, O16, Ne20, Na23, Mg24, Al27, Si28, P31, S32, Cl35, Ar40,
Ca40, Ti48, Cr52, Mn55, Fe56, Ni58. All 20 resolve to supported MESA chem IDs;
the stderr startup record lists actual IDs and mass fractions. Missing IDs are
fatal, not silently replaced. This preserves FreeEOS table X/Z, but approximates
the elemental distribution with one isotope per element. It does **not** use
FreeEOS's elemental-average atomic weights as isotope masses or alter isotope
mass fractions to mimic them. Mean ion/electron properties and Skye selection
use the actual MESA composition, unlike the old H/He/C/O=.70/.28/.01/.01 plotter.

All runtime products are owned by
`freeeos_biquintic_progress/experiments/014_active_domain/`. The wrapper sets
`MESA_CACHES_DIR` and `MESA_TEMP_CACHES_DIR` there, symlinks all 30 installed EOS
bins into `runtime_cache/eosDT_cache`, and requires temporary-cache safety to
remain enabled. A cache write would go to a private temporary file then replace
a private cache link, not truncate its installed target. The wrapper detects
and rejects any replaced link after execution; the smoke test regenerates none.
Installed cache size/mtime inventories are compared before/after (not full
content hashes); the observed inventory remains unchanged. No installed cache
regeneration or table generation is requested.

## Smoke evidence and lead handoff

Eight individually flushed queries pass: both requested points, NaN/Infinity,
unsafe exponent, malformed/extra-token input, and successful recovery. Outputs
are in `smoke.tsv`, `smoke.stderr`, `smoke_results.json`, and `probe_metadata.json`
under `014_active_domain`; all cache links remain installed-bin symlinks.

| q | logT | logRho=q+2logT-12 | ierr | FreeEOS | Skye | Cv |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 2.51 | 6.51 | 3.53 | 0 | 0 | 1 | 1.232518553e8 |
| .8 | 6.765 | 2.33 | 0 | .83692 | .16308 | 1.917896893e8 |

Other fractions are zero at these two points. Feed every lead face/query
coordinate into this persistent process and use the returned fraction/status,
not component-owned q/T cuts, to characterize active contributions. The probe
does not implement face analysis or assign a Skye-only defect to FreeEOS.
