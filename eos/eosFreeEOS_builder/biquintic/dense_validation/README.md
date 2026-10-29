# Dense inherited-query validation

`../dense_audit.py` audits the archived `X=0.70, Z=0.02` table at exactly the
queries in `freeeos_report/plots/offgrid_sample_data.npz`. It does not sample new
points, rerun the legacy cache, call a fresh FreeEOS, rebuild MESA, or duplicate
the lead's original 12,000-point matched-query audit. The inherited archive has
16 independently stratified, strictly off-node points per cell, over 86,179
cells: 1,378,864 saved queries.

## Reproduce

From the workspace root, using the existing Python/NumPy/Matplotlib environment:

```bash
mkdir -p mesa/eos/eosFreeEOS_builder/biquintic/dense_validation/build
python3 -B -m unittest discover \
  -s mesa/eos/eosFreeEOS_builder/biquintic/dense_validation \
  -p 'test_*.py' -v \
  > mesa/eos/eosFreeEOS_builder/biquintic/dense_validation/build/tests.log 2>&1
mkdir -p freeeos_biquintic_progress/experiments/004_dense_validation
python3 -B mesa/eos/eosFreeEOS_builder/biquintic/dense_audit.py \
  --output freeeos_biquintic_progress/experiments/004_dense_validation \
  --chunk-size 16384 --suspicious-cells 32 \
  > freeeos_biquintic_progress/experiments/004_dense_validation/run.log 2>&1
```

The script imports the lead's `source_constraints`, `construct_twists`,
`evaluate`, and `thermodynamic_partials` without changing them. Only 1D query
chunks are passed to the core evaluator. `--table` and `--archive` override
read-only inputs. `--no-figures`, `--no-query-archive`, and `--suspicious-cells 0`
can disable the respective products. Python bytecode writing is disabled by
the script; Matplotlib's configuration cache is redirected inside the output
directory. Synthetic tests use only this directory's ignored `build/` tree.
The default output is the progress report's
`freeeos_biquintic_progress/experiments/004_dense_validation`, never the original
assessment. The inherited archive under `freeeos_report/plots` remains a
read-only input. `relocation_manifest.json` in the corrected experiment records
the original mistaken destination, the proper rerun command, verification, and
cleanup of the generated duplicate; no original assessment assets are changed.

## Candidates and domain

- `radius2_mean_total`: total potential constraints, centered radius-2 mean
  twists, interpolated total potential, no radiation addition.
- `radius2_mean_analytic_split`: subtract analytic radiation from low-order
  physical constraints before native transformation/twists, interpolate the
  remaining potential, then add exact analytic radiation to native derivatives
  before computing thermodynamic state and partials.

Both candidates use identical inherited queries on the intersection of their
complete-stencil domains. The core source mask excludes non-native nodes,
`logT < 3.2`, printed-density disagreement at absolute log10 tolerance `1e-10`,
nonfinite constraints, and invalid source Cv/pressure/chiRho. These exclusions
propagate through both centered twist stages, not just the cell's four direct
corners. Queries below `logT=3.2` are also excluded explicitly. There is no
extrapolation or incomplete-stencil fallback.

`results.json` reports ordered, mutually exclusive exclusions summing to all
saved queries, overlapping diagnostic exclusions, source masks, complete-cell
counts, per-candidate available counts, and the final common intersection.
Archive shape, cell indices, strictly off-node phases, boolean masks, and the
inherited eligibility formula are validated before evaluation.

Subset masks are frozen from the saved archive, not reclassified using the
candidate's pressure:

- `all_valid`: every accepted complete-stencil saved query.
- `broad_own_use`: `all_valid & saved eligible`; inherited full-use rule
  `logT>=3.1 and (logT>=6.5 or logQ<=-0.5 or logRho<=-6.5)`, using the original
  `1e-12` boundary tolerance, with the stricter candidate domain above.
- `material_own_use`: `broad_own_use & saved material`.
- `material_valid`: `all_valid & saved material`, without the own-use cut.

The inherited material flag means legacy `Prad/Pgas < 1e-3`. Broad own-use is an
offline geometrical selection rule, not a freshly evaluated runtime blend
fraction. Queries outside these subsets must not be described as own-use
failures.

## Diagnostics and interpretation

Each candidate/subset has separate finite, positive, negative, zero,
nonpositive, and nonfinite counts for `Cv`, `P_rho`, `Pgas`, `E`, and `S`.
The failure union and the count with all five fields finite and positive are
reported without double-counting a query. Minimum finite values, worst finite
failures, and coordinate examples include saved flat index, sample phase index,
source-cell indices, `logQ`, `logT`, and `logRho`. Every failure coordinate is
preserved in `positivity_failures.npz`, not just the JSON examples. E/S positivity
is a requested sign diagnostic, not a general thermodynamic stability theorem.

Returned-function first-law/Maxwell identities use physical natural-log
partials from `thermodynamic_partials`:

```text
dpe = (rho*E_x + P_t - P)/P
dse = T*S_t/E_t - 1
dsp = -rho*T*S_x/P_t - 1
```

The physical-column route is computed separately from `E_rho`, `P_T`, `S_T`,
`E_T`, and `S_rho`, at the same queries. Balanced returned-function residuals
normalize the additive cancellation by the largest participating term.
Condition numbers are sums of absolute terms divided by the original
denominator; the additional well-conditioned summary uses a limit of 1000 and
counts residuals exceeding `1e-10`. Every distribution records total, finite,
and nonfinite counts, signed quantiles, and absolute median/p90/p99/max. A zero
denominator is excluded from finite quantiles but explicitly counted, not
silently treated as an identity pass.

The archive contains legacy returned-function (`mathematical_*`) and
physical-column (`physical_*`) residuals plus `Cv_mismatch`, **not** legacy
`Cv`, pressure, energy, entropy, or response values. Their distributions are
reported on the same accepted subsets as historical closure context only.
`Cv_mismatch` is legacy returned `dE/dT` divided by the legacy Cv column minus
one; it is not a candidate error, candidate Cv change, or source-accuracy
measurement. Candidate-relative departures are produced only if actual
`legacy_<field>` or `<field>` state arrays exist in a supplied saved archive;
missing fields are explicitly listed. A zero/nonfinite legacy denominator is
counted as unavailable. Total-versus-split relative departures are explicitly
labeled representation sensitivity, not fresh-source error.

Suspicious cells prioritize requested positivity failures in material own-use,
then broad own-use, then all valid queries. The largest total/split Cv changes
supplement this selection. Centered radius-1 mean, radius-2 first/second routes,
and radius-3 mean variants are evaluated only at the inherited saved phases in
those cells, with incomplete variants explicitly excluded. Their sign counts
and departures from the matching radiation radius-2 mean appear in
`suspicious_cell_twist_sensitivity.json`; no new random queries are introduced.

## Products

- `results.json`: provenance/checksums/command, complete counts, positivity,
  identities, legacy context, optional departures, timing, and limitations.
- `query_diagnostics.npz`: original coordinate shape, frozen flags, common
  mask, candidate fields and identity/conditioning arrays. Excluded values are
  NaN. Saved positive fields are directly auditable at every accepted query.
- `positivity_failures.npz`: all failed indices, values, coordinates, source
  cells and inherited subset flags, separately for each candidate/field.
- `cell_maps.npz`: sparse per-source-cell median absolute identity residuals,
  selected counts, failure counts and maxima, with original grid/cell indices.
- Candidate identity-median maps, candidate positivity-failure maps, and saved
  legacy-context maps, each as PNG/PDF when figures are enabled.

Figures distinguish `all_valid` from `material_own_use`. They use physical
log-density/temperature axes and `coolwarm` centered at log10 `-3`, displaying
`[-7,+1]`; colors are clipped for display only. Identity pixels are per-cell
median absolute residuals over accepted subset queries. Positivity pixels are
maximum finite negative magnitude divided by that same cell/subset's median
absolute field value; cells without finite negative failures are grey, while
nonfinite/exact-zero failures have black markers. White contours mark the
inherited broad own-use boundary. Identical blue identity maps demonstrate
potential closure, not physical accuracy or stability outside the sampled
queries. Generic local Hermite C3 face jumps and the lead's table-only source
limitations remain unchanged.

## Completed dense run

The initial full audit completes in 45.7 seconds with Python 3.13.2 and NumPy
2.1.2; the corrected progress-report rerun records its own timing and exact
command in `results.json` and `run.log`. All 12 focused/synthetic tests pass,
including a progress-report default-output regression check. An artifact recount
confirms that saved coordinates, cell indices, eligible/material flags, and
every reported requested sign-failure count match the generated query archive.

Of the original 1,378,864 saved queries, 1,261,376 are accepted. Ordered
exclusions are 33,104 below `logT=3.2`, 1,504 with direct printed-density corner
mismatch, and 82,880 outside the common complete twist-stencil domain after
earlier cuts. No excluded query is extrapolated.

| Accepted subset | Queries | Total-potential failed queries | Analytic-split failed queries |
| --- | ---: | ---: | ---: |
| All valid | 1,261,376 | 55 | 55 |
| Broad own-use | 978,486 | 0 | 0 |
| Material own-use | 305,095 | 0 | 0 |
| Material valid, without own-use cut | 587,985 | 55 | 55 |

For each candidate the all-valid failures comprise 37 negative Cv values, 18
negative `P_rho` values, and 8 negative E values, with overlaps yielding 55
distinct failed queries. Pgas and S have no sign failures; no returned state or
partial field is nonfinite. Every failed query has inherited `eligible=False`.
The worst total-potential values occur around source cell `(450,118)`:
Cv `-2.3177282e10`, `P_rho -2.0684489e14`, and E `-2.3342812e12`, at distinct
saved queries with `logQ~3.42..3.44` and `logT~4.47..4.48`. Exact coordinates,
units inherited from the core, and all failures are preserved in the products.

Maximum all-valid returned-function identity residuals are `6.96e-15` for the
total candidate and `5.14e-15` for analytic split; broad own-use maxima are below
`1.0e-15`. There are no well-conditioned residuals above `1e-10`. Material
own-use split-versus-total Cv departures have median absolute `5.44e-12`, p99
`1.33e-10`, and maximum `2.75e-10`; these are representation comparisons only.

Twist sensitivity uses 512 inherited queries in 32 selected cells. Radius-1
mean and the two radius-2 single routes each have 48 failed selected queries,
all outside broad own-use. Radius-3 mean accepts only 448 of the 512 selected
queries, with no failures among that reduced domain: this must **not** be
interpreted as a stability repair, because its stencil excludes 64 queries.
The actual saved archive has no legacy state-value arrays, so no candidate
state-relative departures from legacy are claimed. These results describe
sampled table-only interpolation, not fresh-source accuracy, a runtime EOS
blend, or positivity between untested query phases.
