# TASK053: isolated same-polynomial arithmetic control

Layout note (24 November 2026): `arithmetic_weights.f90` is now built from
`mesa/eos/private/arithmetic_weights.f90`; the local kernel is the retained
historical arithmetic-control variant. Current production acceptance uses
`supplemental_hdf5/` real legacy-EP adapter tests. Earlier evidence below is
unchanged.

This is a diagnostic alternative, not an edit or replacement of the frozen
native kernel, runtime adapter, API bridge, installed EOS, or stellar build.
All source-table values, cache bytes, preprocessing, cell geometry, source
composition, radiation restoration, guards and27-slot mapping remain unchanged.

`arithmetic_weights.f90` evaluates the same tensor Hermite polynomial in
binary64 using Horner differentiation of explicit degree-five coefficients.
Before contraction it subtracts the lower-left corner's constant value from
all four value nodes and adds it back only to the potential. The tensor value
basis sums to unity mathematically, so this changes evaluation conditioning,
not the represented polynomial. It changes no derivative nodes. Each36-term
sum uses Neumaier compensation. This is neither smoothing nor an f/T refit,
and it does not change the interpolant's continuity class or source gauge.

`freeeos_biquintic.f90` is an isolated copy of the baseline module, with only
the real evaluator delegated to this helper after its original validation.
The public API and all original guard/output poisoning behavior remain.
The AD evaluator deliberately remains the baseline arithmetic; the private
runtime adapter calls the real path. Do not infer that the AD path has gained
the arithmetic improvement, or promote this copy to production on that basis.

The Makefile compiles only small private drivers, linking the lead's existing
experiment050 private EOS and installed dependencies. It never builds EOS or
writes into an installed prefix. Baseline and variant binaries have separate
names/directories; no frozen binary is overwritten.

`repeat.py` imports unchanged051 solver/reference recipes and changes their
in-process output location and client executable location only. It writes a
new full2544-trial variant dataset and independent protocol/forcing audit.
The current hardened client is unchanged and its actual hash is recorded.
The original051 dataset, including both tight-stress failures, is read-only.

`validate.py` compares cached native jets with80-digit decimal evaluation of
the exact binary64 corner data; samples96 declared-face interior offsets and
64 deterministic clean cells. It also replays051's neighboring binary64
states without changing forcing/tolerance and compares baseline/variant API
success, guards, routing, hybrid rows and opt-out results on matched queries.
No source accuracy or exact-face API representability claim follows.

All output is under experiment053. `run.sh` enforces inherited allocation
placement,2-library-thread limits and8GiB virtual memory. The complete original
kernel and047 adapter suites are separately run against private variant
drivers. Every scientific script refuses to overwrite its existing dataset.
The report distinguishes improvements, unchanged/default outcomes, residual
limitations and the frozen baseline failures; this control does not authorize
stellar or production replacement.
