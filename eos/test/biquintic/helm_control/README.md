# Native HELM control (TASK048)

This private benchmark measures algebraic one-sided heat-capacity traces
on the existing MESA HELM table, at fixed audited solar composition.
It does not change or install an EOS implementation.

`helm_probe.f90` links the installed `eos_lib:helmeos2_eval`, exports the
actual linear axes and composition, and returns electron/positron, ion,
radiation, Coulomb-by-subtraction, and total Cv and physical log partials.
`native.py` differentiates the same nine-slot Hermite polynomial, using
each linear cell width. Final endpoint traces use 80-digit arithmetic;
independent extended-precision evaluations check the result and offsets.
The table's logarithmic spacing must not be substituted for these widths.

`run_control.py` freezes 96 nearby native faces around the canonical 24
FreeEOS cases and 48 faces around a predefined 12-state hotter lattice.
At tangential native knots it uses the upper linear midpoint, avoiding
ambiguous intersections. It checks all used source/cache corners, runs
six shrinking normal distances on each side, and calls the unchanged
default-MESA mixture probe using private cache/output directories.
The two EOS families' faces are not at matched physical coordinates.
The 36 seed states also have exact-coordinate HELM queries.

`validate_results.py` checks endpoint precision and finite differences of
the installed kernel's returned total Cv inside cells. Finite differences
and direct double-precision derivatives can reach cancellation floors;
these are recorded, not silently treated as exact limits.
`summarize.py` creates CSV traces and a two-panel figure.

The runners deliberately guard the authorized node/job (rzwhippet14,
5528075), impose 8 GiB address-space limits, and set numerical threads to
two. For a later allocation, review and update the explicit guards first.
See the experiment's `commands.sh` for reproduction. Build only the
single private driver, with `make -j1`; no MESA rebuild is needed.

Tests: `python3 -m unittest discover -s PATH/helm_control -v`.
The artifact root is `freeeos_biquintic_progress/experiments/048_helm_control/`.
Its `results.json`, `independent_validation.json`, `face_traces.csv`, and
dated task report carry the quantitative results and limitations.
