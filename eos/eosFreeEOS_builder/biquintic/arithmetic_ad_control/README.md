# TASK059: isolated third-order AD lift

This diagnostic keeps053's real evaluator and arithmetic helper unchanged,
and changes only `evaluate_cell_ad` in a new private kernel copy. Canonical
and053code/cache/binaries stay frozen. All original validation and failed-output
poisoning remain; no fallback, clipping or state/composition policy changes.

At the input AD values, the real evaluator supplies ten partials through total
degree3. Let `delta_u=u-u%val` and `delta_t=tau-tau%val`. The AD result is the
degree3 Taylor polynomial in these zero-valued AD deltas, with the usual
factorials, including both mixed third-order terms. Because terms of total
degree4 or higher in zero-valued deltas vanish in an order3 jet algebra, this
is the exact mathematical composition rule for arbitrary nonlinear input
third-order jets. It does not replace the source polynomial by a cubic over
the cell; the expansion point and stable real jets are recomputed per query.

Only third-order AD is supplied. Floating-point lifting adds ordinary rounding,
and the real polynomial's remaining arithmetic error is inherited, not fixed.
Finite inputs whose composed derivatives overflow still fail strictly.

The Makefile compiles small private drivers against the existing050private EOS
and installed dependencies, with build/library parallelism<=2 and8GiB memory.
It never rebuilds/installs EOS. The unchanged2059-pair cell regression and047's
118-check adapter suite are run in private outputs.

`jet_driver.f90` accepts arbitrary ten-component input AD jets. `validate.py`
replays160 actual cached points against053's frozen80-digit reference, checks
unchanged real output, and composes the original degree-five Hermite basis
directly in an independent80-digit two-variable series ring. Identity, linear
shear, nonlinear all-orders, zero-first/nonlinear and constant seeds exercise
all mixed orders. Separate exact-series composition of the stable real jets
isolates lifting roundoff from inherited polynomial-evaluation error.

The bounded400-scenario actual-API check compares all124 payload fields with
053, including explicit faults and pure/blended/inactive routing. That API
uses the real path: parity does not constitute an end-to-end AD stellar test.
No full consumer/stellar rerun or promotion follows from this diagnostic.

The frozen JSON labels `full_payload_bitwise_equal` and
`identity_real_unchanged_bitwise` are legacy names, not bit-pattern assertions.
The implementation uses `numpy.array_equal`; the API comparison sets
`equal_nan=True`. These flags mean exact numerical array equality, with
matching NaNs accepted in API records and signed zeros not distinguished.
The meaningful API evidence is202 finite successful-return arrays agreeing
and198 matching expected error statuses (99composition102,99input103).
Arrays attached to failed public calls are not physical EOS results. Output
files and flag names remain unchanged; no numerical rerun accompanies this
terminology clarification.

All sources, commands, raw seeds/responses and reports remain under the new059
scope. Scripts reject existing output datasets. The report documents the
measured remaining errors and exact baseline/variant identities.
