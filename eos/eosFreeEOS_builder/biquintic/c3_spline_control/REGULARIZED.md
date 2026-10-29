# 031 explicit roughness-selected C3 control

Frozen reusable implementation: `regularized.py`; tests: `test_regularized.py`.
No selected runtime, source solver, new tables or installed caches. The
unregularized experiment-029 rank failure remains unchanged.

## Reproduce

```sh
python3 -B mesa/eos/eosFreeEOS_builder/biquintic/c3_spline_control/test_regularized.py
python3 -B mesa/eos/eosFreeEOS_builder/biquintic/c3_spline_control/regularized.py --overwrite
```

Uses only the archived experiment-015 active geometry: 35x25 full source nodes,
18x13 coarse retained nodes, the same 204 odd/odd heldouts and 5,000 offnodes.
No C4 refits, patch reselection or new random queries. All output is exclusive
to `freeeos_biquintic_progress/experiments/031_regularized_c3/`:
README/JSON, paired coefficients/query/state arrays, test logs and validation log.
Input/table/source/archive hashes and prior solver attempts are recorded.

## Selector and normalization

The six linear relative thermodynamic combinations are unchanged from 029/015.
Fix response importance .1 on P_rho, P/rho*chiT and TCv. Let source potential
scale be `s=median(abs(f))`, native coordinates `u=ln10*logQ`, `tau=lnT`,
mesh spacings `h_u,h_tau`, and native patch area A. The objective is

```text
J = mean(weighted_relative_source_residual^2) + lambda * R
R = [h_u^8 integral(f_uuuu^2 du dtau)
     + h_tau^8 integral(f_tautautautau^2 du dtau)] / (A*s^2)
lambda = .001, .01, .1
```

The mean averages over **all six combinations/all nodes**, not separately per
field. Response importance enters squared. The roughness rows multiply by
sqrt(lambda). Coefficients are scaled by s. This is an explicit derivative
selector, not a coefficient minimum norm, fitted gauge, physical energy or
rotationally invariant prior. Changing h or s changes the declared selector;
lambda is not a universally transferable production parameter.

Use 2 Gauss points per interval for fourth-derivative factors and 6 for mass/
value factors. Degree-five splines have degree-one fourth derivatives; their
squares are degree two and value products degree ten, so this piecewise tensor
quadrature integrates the objective exactly. Sparse tensor factors implement
`K4_u tensor M_tau` and `M_u tensor K4_tau` without dense Gram matrices.

The roughness kernel is tensor bicubic polynomials, and >=4 source values per
axis eliminate their common homogeneous kernel. Thus the augmented objective
removes 029's structural ambiguity. Coarse QR confirms rank 1,064/1,064 for
every lambda; this does not prove physical superiority or global positivity.

## Resource and acceptance bounds

16 BLAS threads, 4,500 unknowns, 90 seconds per fit; no new solve after local
20:00. Coarse augmented matrices are 6,300x1,064, about 51 MiB dense, solved by
column-equilibrated rank-revealing `gelsy` QR. Full matrices are 24,834x3,744,
about 709 MiB dense, exceeding the explicit 128 MiB dense-work-array bound:
they use sparse column-equilibrated LSMR instead, not dense normal equations.

Full initialization uses exact single-to-double-knot insertion of the saved C4
coefficients and solves a **shifted-residual correction**. This is conditioning/
initialization, not a selector; the augmented objective is unique. Accept sparse
stops 2/5 (least-squares stationarity), at most 40,000 iterations, plus normalized
augmented-gradient norm <=1e-10. Full rank is analytic, not a claimed numerical
full QR measurement. The recorded LSMR condition estimate is not a rigorous
matrix condition number. A coarse sparse-versus-QR check differs by 6.1e-14 in
scaled coefficients. Whole-process peak RSS is about 292 MiB.

An initial zero-start full run hit 20,000 iterations for .001/.01; .1 returned
a residual stop with gradient about 2.3e-8. That historical .1 result **is not
accepted by the final stationarity policy**. Prior-attempt metadata is preserved;
final paired arrays contain only accepted final fits. No failed coefficient
vector or offnode prediction is exported as a selected solution.

## Results

Coarse Cv archived-field errors on exactly the same 204 heldouts:

| Candidate | Cv p99 | Cv max |
| --- | ---: | ---: |
| Saved C4 thermo .1 | 2.541e-4 | 9.516e-4 |
| C3 lambda .001 | 2.852e-4 | 3.932e-4 |
| C3 lambda .01 | 3.146e-4 | 3.712e-4 |
| C3 lambda .1 | 4.711e-4 | 7.261e-4 |

No universal winner: C3 Cv maxima improve but p99 is worse; Pgas/S p99 improve,
E can worsen. All node constraints remain **soft**, not exact. Coarse TCv nodal
p99 discrepancies are 3.16e-6, 8.40e-6 and 5.24e-5 as lambda increases. Full
JSON includes all six nodal combinations and all physical-field errors.

Full lambda .001 remains a failure at 40,000 iterations. Full .01/.1 converge
with 28,982/17,576 iterations and stationarity about 3.3e-13/3.1e-13. These use
the same 5,000 saved offnodes. Cv p99 departure from BQ is 1.345e-3/1.327e-3;
scaled Cv_t p99 departure is .186/.180. These are **model differences, not
offnode/source accuracy**. Offnode values divide by abs(BQ value), partials by
max(1,abs(BQ field),abs(BQ partial)). Physical `_x` holds T fixed; `_t` holds
rho fixed. No gauge/radiation constant is refit.

All accepted sampled nodal, heldout, offnode and face fields/responses are finite
and positive. One-sided `nextafter` response/Jacobian face limits agree within
about 2.3e-9 symmetric scaled roundoff on tested faces; this is not an exact
rational zero or C4 claim. Double knots remain C3, and sampled positivity is
not a bound over whole cells or patch ends. Patch ends are not blended.

Four regularized tests pass: exact fourth-derivative Gauss integration, penalty
of a nonzero 029 source-null witness, bicubic zero-roughness reproduction with
full-rank stationarity, and exact C4-to-C3 knot insertion. The companion six
spline tests also pass on the final code. No further experiments are part of
this frozen control. These are table-only selector sensitivity/consistency
results, not fresh-source accuracy or production parameter selection.
