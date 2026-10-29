# Bounded offline composition holdout

This isolated next-stage risk-assessment sidecar writes only
`freeeos_biquintic_progress/experiments/012_composition_control`. It does not
implement or import a MESA runtime selection/composition route. No solver, new
tables, installed caches, additive-volume mixing, or changes to FreeEOS mixture
physics are involved. The C4 lead, C3 sidecar, and adapter work remain untouched.

## Experiment

At fixed Z=.02, withhold X=.70 completely from the potential mixtures:

- Linear: X=.60/.80, weights (1/2,1/2).
- Cubic Lagrange: X=.50/.60/.80/.90, weights (-1/6,2/3,2/3,-1/6).
  **Negative weights, not convex:** positive source potentials do not imply
  positive heat capacities or compressibility; source errors can be amplified.

Read all five existing format-51 native grids, require identical authoritative
header fields and coordinate arrays, check composition headers, record and reverify every
source SHA-256. Use the lead's `source_constraints` imported explicitly as
`leadsource_constraints`: six native low-order constraints, logT>=3.2,
finite supplied columns, native flags, positive source Cv/chiRho and matching
printed density coordinates. The all-source intersection is passed through
both centered-radius-2 mean twist stages before any mixture is made.
Strict mode also requires matching display maxima; the explicitly authorized
Z=.02 metadata mode below instead audits their redundancy independently.

This centered BQ prototype does **not** cover every inherited own-use point
with logT>=3.1. Its unchanged source cut is logT>=3.20; two radius-2 twist
stages consume four native .02 temperature steps, so complete nodes/cells begin
at logT=3.28. The interval **[3.1,3.28)** remains uncovered. An upstream
FreeEOS low-temperature/off-by-3.1 fix alone does not close this representation
gap. No source-cut relaxation, one-sided stencil, extrapolation, fallback or
runtime coverage extension is implemented by this composition sidecar.

Mix all nine **nodal** potential jets with composition-only weights, then
evaluate the biquintic and transform to physical rho/T derivatives. Default:
analytic radiation-reference split with the current lead CRAD. Restore the same
radiation coefficient once, because the weights sum to one. E/S/Cv/chi are total;
only printed P is gas-only, and total P is reconstructed with current CRAD.
The split remainder is not a historically verified pure-material EOS.

The main independent archived-composition reference is printed X=.70 at
12,000 deterministically sampled common native SOURCE nodes, never fresh
physical accuracy. Require complete common cells, common source/corner material
masks and inherited own-use cuts. Material masks intersect all five printed
sources: current Prad/Pgas<1e-3, -8<=logQ<=5, logT<=7.8. Own-use means
logT>=3.1 and (logT>=6.5 or logQ<=-.5 or logRho<=-6.5), with inherited tolerance.
The .70 table participates in reference/domain selection only, never X weights.
Compare total/gas P, E, S, Cv, Cp, chiRho, chiT, Gamma1, **Gamma3** (not
Gamma3-1), and grad_ad against printed .70 values. Direct .70 redundant-column
closure is recorded separately to distinguish composition effects from source
thermodynamic closure.

The control applies the **same X weights to the same physical printed columns**,
not to logarithms and not MESA's actual X scheme. At SOURCE nodes it has no
rho/T interpolation error. Its optional off-node context uses explicitly
bilinear native-grid interpolation, not MESA's rho/T scheme either.

Start from the archived `freeeos_report/sample_residuals.csv`; the direct .70
radius2 own-use selection reproduces 6,713 queries at Z=.02. Intersect with the
common all-source stencil domain, reporting excluded queries. Compare mixed
versus direct .70 biquintics: **not accuracy**. Report shared all-own-use and
material-corner subsets, algebraic and finite-difference rho/T identities,
positivity/nonfinite failures, extrema, and nodal-mix/evaluated-jet commutation.
Finite-difference perturbations only test the returned functions; they do not
sample a fresh physical reference. Optional +/-h X weights at fixed rho,T verify
analytic composition partials of this same potential, not source chemical physics.

## Gauges and provenance

Historical E/S gauge alignment and fixed metal pattern at each Z are assumptions
**not verified** by filenames, headers, or the builder's present configuration.
Helium is assumed to vary as Y=1-X-Z. No inference about chemical potentials,
mixture entropy, or historical metal patterns is justified here.

For an exact composition-only gauge shift f -> f+A(X)+B(X)T, at fixed X:
P and Cv are unchanged, E -> E+A(X), S -> S-B(X). Tests apply exact analytic
gauge jets to each source and combine them with the same weights, checking
both invariants and the expected gauge changes. E/S positivity is gauge dependent.

There is an additional representation risk: a quintic polynomial in lnT cannot
represent B*T exactly between knots. The same nodal gauge shift is exact at
SOURCE knots but generally introduces off-node Cv/E/S drift. Tests deliberately
expose this instead of claiming nodal-polynomial gauge invariance. The real-table
audit separately measures exact analytic-gauge invariance and this nodal drift.
An analytic gauge treatment or independently aligned sources would be needed
before drawing stronger composition conclusions. This sidecar does not alter
the lead representation to implement such a treatment.

## Reproduction

From the workspace root, using the existing Python/numpy installation:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s mesa/eos/eosFreeEOS_builder/biquintic/composition_control -p 'test_*.py' -v
PYTHONDONTWRITEBYTECODE=1 python3 mesa/eos/eosFreeEOS_builder/biquintic/composition_control/run_control.py --derivative-x --other-z
PYTHONDONTWRITEBYTECODE=1 python3 mesa/eos/eosFreeEOS_builder/biquintic/composition_control/run_control.py --allow-audited-endpoint-metadata --derivative-x --other-z
```

`--other-z` cheaply repeats the archived holdout at Z=.00/.04, still focusing
Z=.02. `--sample-nodes` can lower the bounded sample. `--total-potential` runs
the optional total-polynomial control and replaces results in this same output
directory; retain the default split results for the primary assessment.

Artifacts: `results.json`, `RESULTS.md`, `Z02_samples.npz` (and Z00/Z04 when
requested and their input gates pass). Failed cases instead produce
`Z02_input_gate.json` / `Z04_input_gate.json`, with no composition predictions.
They record commands, code/source/query fingerprints, grid headers,
shared masks, sampled source indices, saved-query membership and per-point
returned fields. Tiny identities are not accuracy, stability, verified historical
composition physics, or runtime readiness. Results inform a later decision;
they do not authorize a runtime implementation.

The explicit audited invocation writes `results_audited_endpoint.json`,
`RESULTS_audited_endpoint.md`, and `Z02_samples_audited_endpoint.npz`, with
similarly suffixed auxiliary/gate artifacts. It preserves the initial strict
`results.json`, `RESULTS.md`, input-gate JSONs and Z=.00 auxiliary samples
unchanged, recording their hashes before and after the run.

## Observed primary input risk

The existing Z=.02/X=.60 header declares logQmax=3.69, whereas its count,
start and spacing imply 5.69 and the other four headers declare 5.69. This
fails both internal header-endpoint consistency and the required matching
header-grid gate. The strict audit records all five SHA-256s and the independent
6,713 direct-.70 query selection, but does **not** evaluate a Z=.02 composition
holdout or silently reinterpret that header. Z=.04/X=.50 declares a shorter
445-node grid versus 527 for the other sources, but actually contains 527*306
printed records; the lead reader rejects its inconsistent record count before
preprocessing. That auxiliary gate also fails, with all source hashes retained.
Z=.00 is an auxiliary passed-grid control, not a replacement for the primary.

## Explicit audited endpoint-metadata resolution

The user authorizes `--allow-audited-endpoint-metadata` for Z=.02 **only**.
This is a reader metadata resolution, not a repair to any source or cache and
not resampling. The raw headers and discrepancies remain in the results.
Without that explicit flag, the original strict gate still fails.

Audited Z=.02 loading uses the independent sidecar function
`endpoint_audit.read_endpoint_table(path, allow_audited_endpoint_metadata=True)`
for **all five sources**, never the default core `freeeos_potential.read_table`
on a wrong-endpoint file. The independent loader is strict by default:
`read_endpoint_table(path)` validates finite format/count/min/positive-step
metadata, exact record counts, derived/display maxima and printed block grids.
Only its explicit audited option can admit redundant display maxima, and only
at Z=.02 after the full printed-grid audit. Source constraints still come from
the unchanged lead `source_constraints`. Core reader validation can therefore
be strengthened later without invalidating this explicit sidecar loading path.
Regression tests replace the default core reader with an exception and confirm
that the audited source loader never calls it.

At Z=.02, `endpoint_audit.py` independently streams the original text, parses
every `logQ = logRho - 2*logT + 12` block label and every record, and requires:

- Exactly the declared integer num_logQ*num_logT records, with num_logT records
  per block and exactly num_logQ strictly monotonic block labels.
- Every block label agrees with declared count/min/step; both reconstructed
  axes are the authoritative arrays passed to unchanged lead source constraints.
  Counts, minima, steps, format and all
  other non-composition header fields still match across all five sources.
- Every printed T agrees with the grid. Every printed rho agrees with
  logQ+2*logT-12, except the already-known X=.70/Z=.02 iq=520 block, whose
  discrepancy must be exactly -1e-4 within 1e-10 tolerance. That entire block
  must remain invalid under the **unchanged** lead source mask. No new
  coordinate exception, coordinate repair, or label exception is authorized.
- All five source SHA-256s are recorded and reverified after evaluation.

Only the **named** `logQ_max` and `logT_max` fields are then excluded from
the header comparison; the implementation does not discard the entire grid
header. Here Q labels cover -10.09 through 5.69 in 527 steps of .03, each
with 306 printed temperatures. These independent records establish that
X=.60's declared 3.69 display maximum is redundant, not an authoritative
physical cutoff. The audit records the declared maximum, derived endpoint,
per-source evidence, exact flag/command and exclusion separately.

Tests reject altered actual grids, counts, minima, steps, block labels,
nonmonotonic ordering and unexpected printed rho/T discrepancies. The flag
does **not** waive Z=.04's incorrect record-count header; that case remains
deferred. The primary Z=.02 holdout can now be scientifically useful without
claiming that inconsistent display metadata have been silently fixed.
