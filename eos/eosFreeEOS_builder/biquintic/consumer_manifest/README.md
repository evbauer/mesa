# Frozen controlled-consumer sample

`build_manifest.py` selects 24 exact archived experiment014 faces before any
solver outcome is observed: six ordinary and six tail cases per native face
direction. Ordinary ranks are evenly spaced over the central 20–70% of the
absolute archived Cv-temperature-gradient jump distribution; tail cases are
the six largest. Only successful actual-MESA pure-FreeEOS seed queries qualify.
This is a stratified diagnostic sample, not an unbiased population estimate.

The script validates complete support and records one-sided values for the
declared radius2/mean, centered, sourcefloor3.2, radiation-split candidate.
Selection is based on the previously measured total-potential archive, not new
solver outcomes. Seed classification does not classify an entire solver path.

From the workspace root:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python3 -B mesa/eos/eosFreeEOS_builder/biquintic/consumer_manifest/build_manifest.py
python3 -B -m unittest discover -s mesa/eos/eosFreeEOS_builder/biquintic/consumer_manifest -v
```

The output is `freeeos_biquintic_progress/experiments/038_solver_manifest/queries_v2.json`.
V2 corrects reversed X/Z metadata labels in V1 without changing any case or
numerical state. The old file is retained as superseded provenance; a guard
checks that the case records remain identical. Header-order tests cover the correction.
Changes to an existing manifest are rejected; a changed study needs a new
output directory and separately declared sample, not overwrite-and-retune.

## Routing, independent accounting, and conditioning

`probe_paths.py --input CSV --output DIRECTORY` deduplicates solver queries,
invokes the unchanged installed MESA probe with experiment-private caches,
and records actual component fractions for each trial and stage. It checks
query order, statuses and normalized fractions. Finite samples are not a
proof of continuous-path selection, and the candidate is not installed.

`review_consumer.py` checks saved scalar 039 paired forcing and trial counts
against routing 041. `review_coupled.py` independently checks coupled 042
common-forcing pairs, own-reference errors and state/span routing 043.
Both scripts account for saved trials, rather than independently recomputing
the EOS or its bracketed roots. Their defaults use those numbered archives.

`conditioning.py` records cancellation indicators at the same 24 seeds in
experiment 044. It compares arithmetic for f and f/T under declared reference
zeros, not fitted interpolants or source accuracy. Its tests explicitly
check gauge dependence and invalid heat-capacity rejection.
