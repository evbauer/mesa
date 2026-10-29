# Bounded fixed-density self-inversion control

Offline sidecar only. Each model supplies its own energy target at the same
saved off-node truth temperature. Recovered temperatures test **self-inversion**,
not independent source accuracy, agreement between models, MESA active fractions,
or a stellar solver. No core, plotter, runtime defaults, installed caches/tables,
or other worker sources are edited. Imported evaluators are read-only.

## Run

From the workspace root, using the existing NumPy environment:

```sh
python3 -B mesa/eos/eosFreeEOS_builder/biquintic/inversion_control/test_control.py
python3 -B mesa/eos/eosFreeEOS_builder/biquintic/inversion_control/run_control.py --overwrite
```

All products go to
`freeeos_biquintic_progress/experiments/016_inversion_control/`:
`results.json`, `trials.csv`, `targets.npz`, `SUMMARY.md`, and the saved
`tests.log`. No dependency install or solver build is required. `-B` avoids
writing bytecode into imported worker directories. Existing audit products
require explicit `--overwrite`.

Default: 100 saved coordinates, four model-specific targets each (**400 targets**),
two starting guesses and two step strategies (**1,600 trials**). `--targets`
allows 1..125 coordinates, never more than 500 model-specific targets.
`--max-iterations` defaults to 60, maximum 100; `--rtol` defaults to 1e-10.

## Models and selection

- Biquintic total potential, radius-2 mean-route twist estimates.
- Biquintic analytic radiation-reference split, with coherent restoration.
- Existing radius-2 C3 split potential as an optional-method control included
  here through its read-only API; no C3 source modification.
- Exact analytic ideal-gas/radiation energy:
  `E=1.5e8*T+CRAD*T^4/rho`. The chosen material Cv and current CRAD define an
  analytic benchmark, not an approximation asserted accurate for this mixture.

Table SHA-256 is pinned to
`e6ca7ce27bee8197bb257a8c4886a06366e6951e58958fbed3735c877ee68502`.
Coordinates come unchanged from
`experiments/001_native_prototype/matched_queries.npz`, with its saved common/
own-use flags; deterministic selection seed is 20261028. Saved-query and
evaluator-file hashes are recorded in `results.json`.

Keep fixed `log10rho=q_truth+2*log10T_truth-12`. Brackets are
`[truth_logT-.30, truth_logT+.08]`, starting at either endpoint. The lower
start differs from truth by 35 q/T faces; the upper start by 9–10 faces.
Every q/T grid-face intersection of the fixed-density path is enumerated;
each intersected cell must have complete native stencils common to all table
models, clear of iq520 and other exclusions. Material node masks require
Prad/Pgas<1e-3, -8<=q<=5, logT<=7.8, and the declared table own-use cuts
(with the clean lower-temperature restriction logT>=3.2).

Finite positive energy/Cv and 0.05<=T*Cv/E<=20 are screened at every face and
segment midpoint for every model. The minimum observed screening ratio is
0.05008 for the table models. This is not an analytic positivity/monotonicity
proof between all checked points; iteration-time unsafe-Cv guards remain active.
These are table-domain masks, not MESA component-contribution measurements.

## Equation, Jacobian, and safeguards

Solve `E(rho,T)-E_target=0` in **lnT at fixed physical density**.
Since q=log10rho-2log10T+12, the interpolation-coordinate temperature derivative
requires the native shear. The existing thermodynamics transformation supplies
the analytic Jacobian `dE/dlnT=T*Cv` of the same potential. Newton never uses a
finite-difference Jacobian or a different model's target.

Full-step Newton applies the entire analytic-Jacobian update, with no damping
or bisection. An attempt outside the original declared-valid bracket is counted
and terminated without evaluating unverified cells. Such a failure is **not
proof of global Newton divergence**; the function could be valid farther away.

The safeguarded strategy updates the sign bracket, accepts a finite positive-
Jacobian Newton proposal only inside it with step <=75% of its width, and
otherwise uses midpoint bisection. Nonpositive/nonfinite Cv/Jacobian evaluations
are counted, including terminal evaluations; nonfinite energy terminates.
Convergence is relative energy residual <=1e-10, not a stopping test using the
known truth temperature. Recovered T and all counters are measured afterward.

## Results (2026-10-28)

237 candidates yield the 100 common conditioned brackets.

| Model family | Lower full Newton | Upper full Newton | Safeguarded both starts | Safeguarded iterations |
| --- | ---: | ---: | ---: | --- |
| Each of total/split/C3 table models | 24/100 | 100/100 | 200/200 | Lower median 6, max 8; upper median 5, max 6 |
| Exact analytic benchmark | 0/100 | 100/100 | 200/200 | Lower 6; upper 5 |

All **800 safeguarded trials converge**; no unsafe Cv is observed in this
screened sample. The 76/100 lower-start table failures and 100/100 analytic
failures are declared-bracket escapes, illustrating starting-guess/step-control
risk even for an exact smooth analytic function. They cannot be attributed
solely to grid-face smoothness defects.

Table safeguarded maximum relative T error is 7.78e-11; maximum relative energy
residual is 5.80e-11. Each table model uses 96 bisections across its 100 lower
starts and none for upper starts. Total/reference-split/C3 counts are nearly
identical on this material subset; no global superiority conclusion follows.
`trials.csv` includes individual recovered T, residuals, statuses, iteration/
safeguard/unsafe-Cv counts, bracket attempts, and accepted face crossings.

Four small analytic tests pass: batch exact-radiation Jacobian/recovery,
full-step escape versus safe recovery, deliberately unsafe Jacobian/Cv fallback,
and an invalid bracket. Separate finite-difference sanity checks on 16 saved
truths are recorded, **not used by Newton**. Their discrepancies become noisier
at smaller steps, particularly C3 (up to 3.6e-6 at step 1e-6); they are not
hidden as source errors or evidence of high-precision C3 accuracy.

This bounded, conditioned inversion control does not validate radiation-
dominated inversion, all ionization/phase regions, active EOS blends, evolution
solver behavior, or independently accurate thermodynamics.
