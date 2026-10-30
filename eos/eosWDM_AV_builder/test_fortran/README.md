# WDM_AV Fortran interpolant tests (WP2)

Standalone test harness for the per-species Hermite interpolant of the WDM_AV
tabulated equation of state. It exercises the three private `eos` modules

- `eos/private/eos_wdm_av_def.f90`   (table / derivs types, module store, lookup)
- `eos/private/eos_wdm_av_load.f90`  (HDF5 reader via forum, lazy load)
- `eos/private/eos_wdm_av_interp.f90` (bicubic/biquintic eval, real + auto_diff,
  isotope rescaling)

## Building and running

```bash
export MESA_DIR=/path/to/this/mesa/checkout
cd $MESA_DIR/eos && make                      # build the eos library first
cd eosWDM_AV_builder/test_fortran && make && ./test_interp
```

`make` links the built MESA libraries (modeled on `eos/plotter/Makefile`) and
copies the executable here as `./test_interp`. Because nothing inside `eos` uses
the `load`/`interp` modules yet, MESA's build does not emit their `.mod` files;
the Makefile force-generates them (via the `eos` sub-make) before compiling, and
adds `$(MESA_DIR)/build/eos/modules` to the include path (those private modules
are not installed to the public include dir). `make test` builds and runs.

## What is tested

- (a) Basis cardinality of the HELM quintic/cubic Hermite polynomials.
- (b) A generic two-exponential function (nontrivial in both axes, exact
  derivatives to all orders): node reproduction, held-out accuracy, finite-
  difference consistency of every returned derivative through third order,
  auto_diff-vs-real agreement for the `f`, `fx`, `fxx` outputs, and edge
  continuity of second derivatives (continuous for order 5, jumps for order 3).
- (b) The physics model (fully ionized ideal ions + ideal nondegenerate
  electrons, no radiation) for he4, order 5 and 3: held-out errors in
  `f, P = rho*fx, s = -ft/T, cV = (ft-ftt)/T, chi_rho = 1 + fxx/fx,
  chi_T = fxt/fx`, and node reproduction.
- (c) Isotope rescaling (Decision D6): a ne20 table evaluated for ne22 vs the
  analytic ne22 model, including points that map exactly onto table nodes
  (exact to roundoff) and the auto_diff rescaling path.
- (d) HDF5 loader: reads a conforming file (the builder-produced
  `data/eosWDM_AV_data/ideal/eosWDM_AV_he4.h5` if present, else a synthetic file
  from `make_test_h5.py`), compares against the analytic model, verifies that a
  file with the wrong `format_version` is refused, and exercises the lazy
  `eos_wdm_av_load_isotope` store path.

## Files

- `test_interp.f90`  the test program (PASS/FAIL lines with measured numbers)
- `make_test_h5.py`  writes a minimal conforming processed HDF5 table from the
  same analytic model (`python3 -I make_test_h5.py <out_dir>`); the test invokes
  it automatically into `scratch/test_data/`
- `Makefile`
- `scratch/`  build/run scratch (generated HDF5 test files); not a deliverable

## Notes

- `f` crosses zero (`f = -kT P1/m_u`), so its held-out error is reported scaled
  by the thermal energy `kT/m_u`; `cV = (ft-ftt)/T` is a difference of two large
  nearly-equal quantities, so its error reflects that cancellation (Decision D16
  recommends carrying `s` explicitly to avoid it).
- String attributes (isotope, element) are intentionally not read by the loader:
  forum reads character attributes with a fixed-length native-char mem type and
  HDF5 has no conversion path from the string types h5py writes, and a failed
  forum read aborts (`error stop`). Isotope identity comes from the caller /
  file name (FORMAT.md section 3). All numeric attributes read fine.
