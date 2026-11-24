# Supplemental HDF5 runtime tests

This directory contains current production adapter tests and historical
raw-source comparison tools. Production modules live in `eos/private`;
these scripts link against a private EOS build. They do not install EOS
libraries or introduce another production source reader or disk cache.

## Current production checks

Run inside the inherited single-node allocation, on its compute node. From
the workspace root, choose a **fresh** output directory; do not reuse frozen
numbered experiment outputs. The lead owns the private EOS build and numerical
before/after comparison. Example commands for an authorized new validation:

```bash
mesa_directory="$PWD/mesa"
sdk_directory="$PWD/mesasdk-x86_64-linux-26.6.1"
output="$PWD/freeeos_biquintic_progress/experiments/089_runtime_layout/worker_new_validation"
tests="$PWD/mesa/eos/test/biquintic/supplemental_hdf5"
# Set MESA_BUILD_JOBS to the build budget agreed with the lead; default is 16.
export MESA_BUILD_JOBS=16
bash "$tests/build_private.sh" "$mesa_directory" "$sdk_directory" "$output" eos
bash "$tests/build_tests.sh" "$mesa_directory" "$sdk_directory" "$output"
```

`build_private.sh` checks installed archive hashes before and after building.
`build_tests.sh` compiles the drivers below against the supplied private EOS
archive and existing dependency libraries. Both scripts accept absolute or
relative directory arguments and resolve the workspace from this relocated
location. The `stellar` stage builds the private downstream stellar packages
only after the EOS stage; it is not needed for these adapter checks.

Prepare a fresh EOS inlist using `FreeEOS_use_HELM_interpolation = .true.` and
`FreeEOS_biquintic_table_dir` set to the desired supplements directory. Use the
same original FreeEOS data and existing binary caches as the production
loader. Use the SDK library paths in `LD_LIBRARY_PATH`. The family stress test
explicitly requests 16 OpenMP threads; reserve that budget and unset any
inherited `OMP_THREAD_LIMIT` below 16 (or set it to 16). Reusable build scripts
do not impose an OpenMP thread limit. A worker may cap its own local smoke
checks separately; those checks are not a 16-thread family result.

| Executable | Arguments | Checks |
| --- | --- | --- |
| `test_legacy_tables` | MESA directory, EOS inlist, mass-fraction file, supplements directory | Standard loader and its legacy EP arrays supply lower constraints; all 89 tables (83 ordinary, six C/O), masks, dimensions, finite values and shutdown are checked. |
| `test_legacy_family` | Same four arguments | Real public EOS handles, missing request context, ordinary/C/O suffix identity, cold/warm concurrency, composition modes and support errors. |
| `native_driver` | See `native_integration/native_driver.f90` | Public component/final EOS API query driver. |
| `dual_handle_driver` | See `co_table_support/dual_handle_driver.f90` | Independent handles with distinct C/O suffixes. |
| `test_control` | MESA directory, mass-fraction file | Renamed control defaults false, public getter/setter round trip, old control name unregistered. |
| `test_reader` | HDF5 filename, expected status | FORUM reader schema and identity validation; synthetic Python fixtures are supported. |

The mass-fraction file is
`mesa/eos/eosFreeEOS_builder/mass_frac_solar.txt`. For example, after preparing
`$output/inlist` and choosing `supplements`, run:

```bash
export MESASDK_ROOT="$sdk_directory"
export LD_LIBRARY_PATH="$sdk_directory/lib:$sdk_directory/lib64:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=16 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"$output/tests/test_legacy_tables" "$mesa_directory" "$output/inlist" \
  "$mesa_directory/eos/eosFreeEOS_builder/mass_frac_solar.txt" "$supplements"
"$output/tests/test_legacy_family" "$mesa_directory" "$output/inlist" \
  "$mesa_directory/eos/eosFreeEOS_builder/mass_frac_solar.txt" "$supplements"
python3 -B "$tests/test_reader.py" --driver "$output/tests/test_reader" \
  --sdk "$sdk_directory" --output "$output/reader_fixtures"
```

## Historical tools and preprocessing

`freeeos_source.f90` and `historical_source_test.f90` are linked only into
`test_tables` and `test_loader`. Those tests independently read raw ASCII
sources to compare constraints and authenticate inputs. Their Python wrappers
are `test_loader.py` and `summarize_validation.py`. They are historical
comparisons, not acceptance tests of the production legacy-EP loading path.
`native_validation.py` and `stellar_validation.py` retain earlier comparison
workflows; use fresh outputs and inspect their required inputs before any run.
Frozen evidence and earlier expected hashes are not updated by relocation.

`produce.py`, `verify.py` and `package.py` are symlinks into the preprocessing
directory. Do not edit them through this test tree. Generation and packaging
instructions are in
[the builder README](../../../eosFreeEOS_builder/biquintic/supplemental_hdf5/README.md);
the schema is in
[FORMAT.md](../../../eosFreeEOS_builder/biquintic/supplemental_hdf5/FORMAT.md).
