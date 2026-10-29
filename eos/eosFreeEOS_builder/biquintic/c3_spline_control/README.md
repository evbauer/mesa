# Double-interior-knot C3 tensor spline control

Bounded table-only alternative potential feasibility audit, **not chosen runtime**.
Degree-five B-splines with interior multiplicity two have C3 potential, not C4.
This control tries the same six relative thermodynamic source-fit combinations
as experiment 015 at response weights **.01 and .1**, reusing its exact four
patches and saved query arrays. No core/runtime/default/table/cache changes,
downloads, other-worker output changes, commits or solver development.

## Reproduce

From the workspace root with existing Python/NumPy/SciPy/threadpoolctl:

```sh
OPENBLAS_NUM_THREADS=16 OMP_NUM_THREADS=16 MKL_NUM_THREADS=16 python3 -B mesa/eos/eosFreeEOS_builder/biquintic/c3_spline_control/test_spline_control.py
python3 -B mesa/eos/eosFreeEOS_builder/biquintic/c3_spline_control/run_control.py --overwrite
```

The runner sets 16 BLAS/OpenMP threads before NumPy imports, and applies a
threadpool limit to each QR. Hard maximum 4,500 unknowns; at most 16 attempted
new fits (four patches, two weights, coarse/full), coarse first. No new QR starts
at/after local 20:00. `--coarse-only` is available. Eight coarse fits complete
in about two seconds total; full fits are deferred for a **structural rank
failure**, not time/memory exhaustion. Largest attempted dense QR matrix is
11,950,848 bytes; observed whole-process peak RSS is about 235 MiB, including
source/archive arrays and library workspaces. No installation/build is needed.

All products live in
`freeeos_biquintic_progress/experiments/029_c3_spline/`: top-level README/JSON,
test/validation logs, and four patch result JSONs plus
`inputs_queries_and_null_audit.npz`. These archive existing coefficients/queries
and homogeneous null witnesses; they are **not new EOS tables or runtime caches**.
No deficient C3 QR coefficients are exported as selected EOS predictions.

## Exactly paired inputs and baseline

The pinned source table SHA-256 is
`e6ca7ce27bee8197bb257a8c4886a06366e6951e58958fbed3735c877ee68502`.
Copy the actual experiment-015 rectangles/selected original indices; assert
patch coordinates equal archived `patch_logQ`/`patch_logT`. Require all patch
source rows native/valid under the original default 3.2 temperature floor and
its printed-coordinate validity rules. Original iq520 is outside these patches.
Archive hashes and local source/version/command fingerprints are recorded.

Original patch starts are even. Coarse source rows are even/even retained rows;
all target original indices are odd/odd and therefore truly withheld on stride 2.
Use **the exact archived order** of `withheld_logQ`/`withheld_logT`, local donor
subset flags, and all 5,000 saved offnode query arrays per patch. Do not regenerate
random samples or reselect patches. Counts are 180/192/170/204 holdouts, 746 total.
`active_chiT` is the 015 geometry label, not measured active MESA fractions.

Archived single-knot C4 coefficient/value arrays and fit diagnostics are reused
without refitting 015. C4 heldout errors are recomputed against the same single
printed reference; the archives also preserve full nodal reconciliation,
offnode/model differences, Jacobians and positivity diagnostics. No C3 rank-failed
candidate is silently assigned those C4 values or included in an error ranking.

`freeeos_potential.py` now has an optional `minimum_logT` parameter, so its current
hash differs from the archived 015 snapshot. Both hashes and the current function
defaults are recorded; the unchanged default `minimum_logT=3.2` is used. All other
015 provenance files match at execution. This is explicit provenance, not a
claim that the whole source snapshot is unchanged.

## Identical source-fit objective

Use native `u=ln(10)*logQ`, `tau=ln(T)` and the existing total-potential six-jet
transform. Total P is reconstructed with current MESA CRAD from printed Pgas;
printed total E/S/Cv and chi fields retain their gauges/conventions. No analytic
radiation split, fitted gauge or fitted radiation constant is introduced here.

Physical-temperature operators at fixed rho are

```text
ft  = f_tau - 2*f_u
ftt = f_tautau - 4*f_utau + 4*f_uu
```

Fit linear combinations `E=f-ft`, `TS=-ft`, `P/rho=f_u`,
`P_rho=f_u+f_uu`, `(P/rho)*chiT=f_utau-2*f_uu`, and `T*Cv=ft-ftt`.
Each equation has relative scale `max(abs(target), .01*median(abs(target)))`.
The last three responses have multiplier .01 or .1; the first three have 1.
Use the same potential scaling and column equilibration as 015, with dense
rank-revealing `gelsy` QR and `cond=1e-13`. No regularizer, boundary-condition
choice, gauge minimum norm or additional high-jet constraints are added.

## Rank failure evidence

For N nodes per axis, sixfold clamped ends and double interior knots yield
`2*N+2` one-dimensional coefficients. Thus more DOFs do **not** automatically
make the six-node-constraint system unique, even when equations outnumber DOFs.
Every attempted patch/weight exhibits the same structural deficiency:

| Patch | Coarse source shape | Unknowns | Equations | QR rank | Nullity |
| --- | --- | ---: | ---: | ---: | ---: |
| ionization | 16 x 13 | 952 | 1,248 | 890 | 62 |
| helium_ionization | 17 x 13 | 1,008 | 1,326 | 944 | 64 |
| dense_tail | 18 x 11 | 912 | 1,188 | 850 | 62 |
| active_chiT | 18 x 13 | 1,064 | 1,404 | 998 | 66 |

The predicted nodal-map rank is `4*Nq*Nt+2*Nq+2*Nt`, with nullity
`2*Nq+2*Nt+4`; all eight observed QR ranks match. This follows the nested
one-dimensional value/first/second evaluation flag dimensions `N,N,2`: the six
supplied tensor-jet combinations omit the `(1,2)`, `(2,1)`, `(2,2)` flag blocks.
The homogeneous witness provides a direct numerical/nonuniqueness check
independent of the full rank count.

Construct one-dimensional nonzero splines `a(u)` and `b(tau)` with **value and
first derivative zero at every retained node**. Each 1D kernel has dimension two.
Their tensor product has all six nodal jets zero:

```text
w       = a*b       = 0
w_u     = a'*b      = 0
w_tau   = a*b'      = 0
w_uu    = a''*b     = 0
w_utau  = a'*b'     = 0
w_tautau= a*b''     = 0
```

Yet the same mode is nonzero at archived heldout/offnode queries. The saved
outer-product witness is one representative, not a complete null-space basis.
Its maximum spacing-scaled retained-node residual is below 5e-15, while its
withheld `TCv` direction can change by hundreds of arbitrary potential units.
Multiplying this mode changes offnode values/Jacobians without changing any six
source constraints. Those directions have **no physical normalization** and are
not reported as EOS values, accuracy errors or positivity failures.

Full double-knot grids have 3,168–3,744 unknowns (all below the 4,500 cap), but
the same structural null modes persist; eight expensive full QRs are therefore
deferred instead of selecting arbitrary deficient solutions. Lower weighted
residuals than C4 would not resolve uniqueness or establish better truth/fidelity.

The optional exact-E/TS/P fit is also deferred: enforcing those three combinations
pointwise and fitting the second responses softly preserves the **same** null
modes. A constrained least-squares nullspace/KKT formulation cannot remove them
by itself. A specified additional penalty, boundary/high-jet condition or explicit
minimum-norm selector would be needed, with separate conditioning, gauge,
holdout and positivity validation. No such selector is silently introduced.

## Continuity and tests

Six tests pass: double-knot multiplicity/partition through fifth derivative,
tensor polynomial representation, third-but-not-fourth basis continuity,
physical response/Jacobian one-sided continuity away from singular denominators,
six-combination equivalence/QR bound, and a nonzero tensor null witness that
leaves all six nodal constraints unchanged.

Polynomial representation is constructed through uniquely determined **1D
value/first/second** collocation and tensorized coefficients. It demonstrates
the space can represent degree-five-by-five polynomials, **not** that the
rank-deficient six-combination fit chooses those coefficients. The analytic
six-constraint fit correctly fails rank rather than being declared successful.

One-sided checks use `nextafter(face,-inf)` and `nextafter(face,+inf)`:
floating limits, not exact rational zeros. Third derivative continuity holds
to spacing-scaled roundoff; fourth derivative basis jumps are about 180 in that
scale, demonstrating C3 rather than C4. Physical response-Jacobian continuity
is checked on the manufactured nonsingular potential. For homogeneous witnesses,
only basis/native-jet continuity is reported: pressure ratios can be singular,
so their physical chi/Cp derivatives are deliberately not interpreted.

**Conclusion:** the requested unregularized double-knot alternative is
underdetermined by these six source constraints. The output records the rank
failure and archived C4 baseline, but cannot claim an improved C3 heldout
Cv/other-field comparison, reconciled nodal physics, positivity, runtime choice,
or fresh-source accuracy without an explicit additional selection rule.
