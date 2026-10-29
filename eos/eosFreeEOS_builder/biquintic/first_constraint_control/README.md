# 033 exact-first C3 control

One bounded coarse active-geometry experiment. Preserve retained-node printed
E/S/P exactly through hard `E`, `TS` and `P/rho` combinations, while soft-fitting
the other three response combinations with the explicit normal-fourth-derivative
roughness selector frozen in 031. **Not selected runtime or fresh-source accuracy.**
No edits to 031, Linnaeus/032, core/runtime/default/table/cache or other workers.

## Reproduce

```sh
python3 -B mesa/eos/eosFreeEOS_builder/biquintic/first_constraint_control/test_constrained.py
python3 -B mesa/eos/eosFreeEOS_builder/biquintic/first_constraint_control/run_control.py --overwrite
```

Existing Python/NumPy/SciPy/threadpoolctl only. Set 16 BLAS threads before
NumPy import; imported 031 helper is checked before/after against SHA-256
`ba30d62f864e094bb52ad0b83f87554b8fff199c4146f1c1e1a6fc297b35edfd`.
All outputs are exclusive to
`freeeos_biquintic_progress/experiments/033_first_constraint/`: README,
`results.json`, `paired_controls.npz`, tests and validation logs. Audit arrays
are not new EOS tables/caches. Source/import/table/archive hashes, versions,
commands, runtime and Linux peak RSS are recorded.

Only three 1,064-unknown fits, one active-geometry coarse patch, lambdas
.001/.01/.1, 90 seconds per fit, no new solve at/after local 20:00. The three
controls take about 1.3 seconds total here, peak process RSS about 267 MiB.
The dense orthogonal basis is about 8.6 MiB and reduced least-squares matrix
about 15.5 MiB; no dense KKT normal equations are formed.

## Same query/reference contract

Use the exact experiment-015 active patch bounds and saved selected indices,
35x25 full grid and 18x13 retained coarse grid. All 234 source nodes are
native/valid. Copy exactly the same **204 odd/odd archived holdouts**, query order,
printed reference, and 031 soft/C4 baselines. No patch reselection, random query
regeneration, fine-grid twist reuse, solver queries or new C4/031 fitting.
`active_chiT` is geometry, not measured MESA active fractions.

The existing table SHA-256 is pinned to
`e6ca7ce27bee8197bb257a8c4886a06366e6951e58958fbed3735c877ee68502`.
Total pressure is current-CRAD `Pgas+CRAD*T^4/3`, with printed total energy,
entropy, Cv and chi conventions unchanged. No gauge/radiation fit. Hard total
P at known T also preserves Pgas; hard TS at known T preserves S.

## Equality-constrained solve

Use the degree-five double-interior-knot C3 basis and the same six linear
combinations from 015/029/031. Hard rows are `E`, `TS`, `P/rho`, all source
nodes: 702 equations. Soft rows are `P_rho`, `(P/rho)*chiT`, `TCv`, relative
source weights with response importance .1. Hard rows retain the same relative
equation scaling for conditioning, but are **equalities**, not large-weight
penalties. No rank-deficient hard/equality or reduced solution is accepted.

The objective within the hard affine constraint set is

```text
sum(weighted_relative_soft_response_residual^2)/(6*source_nodes)
 + lambda * [h_u^8 integral(f_uuuu^2) + h_tau^8 integral(f_tautautautau^2)]
            /(native_patch_area * median(abs(f))^2)
```

The six-field denominator intentionally matches 031; the omitted first-field
residuals are zero. The same exact tensor Gauss-quadrature derivative penalty
is reused read-only. Roughness is an explicit algebraic selector, not physical
energy, an arbitrary coefficient minimum norm or production parameter choice.

Equilibrate columns using hard plus selected soft/rough operator norms. Pivoted
QR of the hard transpose gives a particular solution and an orthogonal nullspace.
Require independent hard-row rank 702/702. Solve the projected soft/rough
problem with rank-revealing `gelsy` QR; require reduced rank 362/362. Check hard
relative residual <=1e-9 and normalized reduced gradient <=1e-10. Independently
reconstruct physical retained E/S/Pgas and require relative error <=1e-9. Any
hard failure stops remaining fits explicitly, rather than accepting minimum norm.
Actual physical first-constraint errors are <=6.35e-13; reduced stationarity is
below 8e-16. All three hard matrices and reduced systems meet full rank.

## Results

Cv archived-field relative errors on the same 204 holdouts:

| Candidate | Cv p99 | Cv max |
| --- | ---: | ---: |
| Archived C4 thermo .1 | 2.541e-4 | 9.516e-4 |
| 031 soft C3 lambda .001 | 2.852e-4 | 3.932e-4 |
| 033 hard C3 lambda .001 | 3.798e-4 | 5.381e-4 |
| 031 soft C3 lambda .01 | 3.146e-4 | 3.712e-4 |
| 033 hard C3 lambda .01 | 3.806e-4 | 5.380e-4 |
| 031 soft C3 lambda .1 | 4.711e-4 | 7.261e-4 |
| 033 hard C3 lambda .1 | 3.993e-4 | 5.366e-4 |

Hard preservation is feasible but **does not guarantee better heldout fidelity**.
Cv p99 worsens versus soft C3 at .001/.01 and improves at .1; maxima remain
below this archived C4 baseline but are not universally better than 031. Paired
Cv improves at 65.7%, 56.9%, 51.0% of targets respectively; aggregate/pointwise
behavior is not interchangeable. Hard heldout E p99 is about 1.15e-5 versus
soft C3 about 6.2e-6: nodal exactness can worsen between-node reproduction.

At lambdas .001/.01/.1, heldout `P_rho` p99 errors are
1.137e-4/1.136e-4/1.122e-4 and `P_T` p99 errors are
2.620e-4/2.610e-4/2.527e-4. These reference physical derivatives are reconstructed
from the printed total-pressure/chi columns, not fresh independent source data.
The full field/response/node/paired summaries are in JSON and raw arrays in NPZ.

All sampled retained/heldout/face fields and responses are finite/positive.
One-sided physical value/response-Jacobian face limits agree within about
1.7e-10 symmetric scaled roundoff on tested faces. No global positivity,
patch-end blend, C4 continuity, stellar solver or production-integration claim.
No offnode source truth is available. Heldout errors are archived-field fidelity,
and continuity/gradient identities are only construction checks.

Three tests pass: manufactured bicubic exact first fields with zero roughness/
objective, incompatible soft responses with **nonzero objective but first fields
still exact**, and dependent inconsistent hard rows rejected. Final audit also
checks hard-selected objective >= unconstrained 031 objective at the same lambda,
all paired IDs, pinned hashes, rank, hard physical precision and positivity.
This isolates the node-preservation versus curvature-reconciliation tradeoff;
it establishes no universal model winner or fresh-source accuracy.
