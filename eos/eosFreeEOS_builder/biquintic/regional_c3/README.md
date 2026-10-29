# Regional coarse regularized C3 gate — experiment 032

This sidecar extends 031's explicitly roughness-selected C3 potential to the
**same 015 ionization, helium_ionization, and dense_tail rectangles and coarse
holdout arrays**. It imports `c3_spline_control/regularized.py` and its stable
`solve_regularized`, `state_at`, and basis helpers read-only. It does not edit
that API, another worker's code/docs/tests, the core, a source table, an old
performance artifact, a cache, EOS dispatch/defaults, or a lead report/PDF.
No compile, installation, debugger, approval, commit, or push is needed.

## Reproduce and bounds

From the workspace root, using the existing Python/NumPy/SciPy/threadpoolctl:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s mesa/eos/eosFreeEOS_builder/biquintic/regional_c3 -p 'test_*.py' -v
timeout 280s env PYTHONDONTWRITEBYTECODE=1 python3 mesa/eos/eosFreeEOS_builder/biquintic/regional_c3/run_control.py
PYTHONDONTWRITEBYTECODE=1 python3 mesa/eos/eosFreeEOS_builder/biquintic/regional_c3/validate_results.py
```

Use `--overwrite` only to replace this experiment's own products. The runner
executes **18 coarse fits only**: three patches × source response importance
**.01/.1** × lambda **.001/.01/.1**. There is no holdout-selected weight or
lambda, no patch reselection, no regenerated query sample, and no full fit in
this first gate. Optional full fits are explicitly deferred rather than started
without evidence of a sub-30-second regional completion.

Each new fit subprocess has a **30-second wall cap including import/build/QR**,
with a 240-second shared fit-launch/process budget and an absolute
**2026-10-28 20:00 PDT** finish fence. The enclosing reproduction `timeout`
also caps whole-process preparation, postprocessing, and hashing at 280 seconds.
Caps/timeouts/statuses are logged, not silently treated as solved. Six-case
coarse geometries have 952, 1008, and 912 unknowns respectively, strictly below
2,000; the API's augmented dense-work bound is 128 MiB. Each fit uses the API's
dense rank-revealing `gelsy` QR and 16 BLAS threads, with no sparse verification
or separate sparse solver development. A timed-out fit is killed by the parent;
the imported API's own longer alarm is not overridden by editing its globals.

## Exact pairing and no target leak

Read the pinned printed source SHA
`e6ca7ce27bee8197bb257a8c4886a06366e6951e58958fbed3735c877ee68502`.
Assert original patch starts/dimensions and both patch axes equal 015 arrays.
Coarse rows retain original even/even indices; the exact 015 odd/odd holdout
order is verified, including original source indices and local-donor flags.
Counts are 180/192/170 holdouts, 542 total. All original patch rows must remain
native valid under the default 3.2 floor and strict printed-coordinate rules.

Training constraints are converted **independently from retained rows only**.
Each fit worker receives only `fit_inputs.npz`: six supplied native constraints,
coarse axes, and retained index arrays. Neither withheld printed values, C4
predictions, BQ donors, nor outside-patch values enter fitting. Tests perturb
every nonretained/cross/withheld/outside row and require identical training
constraints, and reject displaced saved-query order and invalid selections.
Other tests check printed gamma-column mapping and paired improvement signs.

The **same single printed reference** serves C3 and archived C4. Pgas/E/S are
decoded from printed logs; Cp/Cv/chiRho/chiT and gamma1/gamma3/grad_ad use their
actual printed columns (gamma columns 13/14/15, not 8/9/10). Total P restores
current-CRAD radiation exactly as 015. The native constraints are the same
unsplit total specific potential; no gauge/radiation refit is performed. The
015-declared and current core hashes are both recorded, not asserted identical.

## Selector, comparisons, and interpretation

Degree-five tensor splines with double interior knots have C3 potential.
The API solves the explicit objective

`mean(weighted relative six-source-combination residual^2) + lambda * normalized integrated normal-fourth roughness`.

The six source combinations are E, TS, P/rho, P_rho, (P/rho)chiT, TCv; declared
response importance multiplies the last three residuals and enters the objective
squared. `sqrt(lambda)` multiplies roughness rows. Roughness is the exactly
quadrature-integrated normal fourth-derivative penalty from 031, normalized by
native interval spacing, patch area, and median |f|. No new selector is invented.
It resolves 029's unconstrained source-fit null modes; the augmented rank and
stationarity guard must pass. This is **a declared mathematical selector, not
physical truth, a source-accuracy guarantee, or a universal winning weight**.

For each importance, reuse matching 015 C4 **coarse** coefficients and saved
holdout values, without refitting. Nodal C4 values are evaluations of archived
coefficients. Report all eleven fields' retained-nodal residuals, holdout
median/p99/max/worst coordinates, strict per-query improvement/worsening versus
that matching C4, and sampled nodal/holdout positivity. Individual paired
reductions are `abs(C4-reference)/abs(reference) - abs(C3-reference)/abs(reference)`;
positive means improvement, with no denominator change between methods. Nodal
comparisons are retained **training** nodes; odd/odd holdouts remain separate.
Six-combination soft-reconciliation diagnostics and solver objectives are also
retained. Sampled positivity is not a whole-patch positivity certificate, and
C3 is not C4 or a globally blended regional representation.

The BQ comparison is **only the archived 015 global-donor coarse control** on
the same holdouts/reference. It uses clean donors outside the local patch, a
different information context. Its JSON summaries/positivity are copied without
donor/table regeneration. Individual BQ holdout predictions were not saved in
015, so there is no fabricated pointwise BQ paired ranking. The archived
patch-only subset count is preserved as context, not silently substituted for
the full holdout set. No region is labeled **active**: this task supplies no
actual installed-EOS fraction measurement on these regional arrays.

## Artifacts

All products belong to
`freeeos_biquintic_progress/experiments/032_regional_c3/`: `results.json`, small
`RESULTS.md` table, and per-patch `results.json`, `fit_inputs.npz`, six local
coefficient NPZs, six solver JSONs/stdout/stderr logs, and
`paired_queries_and_predictions.npz`. These coefficient/diagnostic arrays are
not EOS source tables or runtime caches. Fingerprints of source/API/archives,
versions, command, actual importance/lambda, runtime caps, rank, stationarity,
objectives, and memory are recorded. API diagnostics are retained verbatim;
if their text mentions the .1 default, the wrapper separately states the actual
.01/.1 argument and squared weighting. No old 015/029/031/017 products are written.

Each candidate's `objective_normalization` explicitly records the **six-node-row
mean divisor**, operator/target `1/sqrt(source_rows)`, potential coefficient
scale, source importance and its square, relative floors, lambda and its square
root, separate data/roughness objectives, and their augmented sum. The local
`solve.source_objective` text states the actual .01/.1 value; any default-.1 API
text is preserved as `solve.API_source_objective_text`, without editing the API.
Final regression independently replays both objective terms from recorded
coefficients/retained constraints, verifying the normalization for all 18 cases.

`validate_results.py` is the final **recorded-output regression after code
freeze**: replay archived coefficients, verify all 18 ranks/importance/lambda,
exact saved-query order and C4 pairing, error quantiles, individual reductions,
positivity and retained-only payload identity. It performs **no new fits**, source
reprocessing or cache access; hashes must stay fixed throughout regression.
It writes `validation_results.json`, adds provenance to `results.json`, and
produces 032's small `README.md` table. All-lambda results remain in `RESULTS.md`.

## Observed coarse gate

18/18 cases solve in 13.86 seconds total; individual augmented QR times are
about .18–.30 seconds, all full rank, with no sampled nonfinite/nonpositive
eleven-field states. **None of the 18 improves Cv holdout p99 over matching
archived C4 importance.** Pointwise wins still occur (for example 152/192
helium-ionization holdouts at importance .1/lambda .001), showing that individual
improvement and tail improvement are distinct. The explicit selector repairs
the null-space ambiguity, not necessarily regional source/holdout fidelity.
No full fit, universal winning weight, measured active geometry, or runtime
representation is selected by this result.
