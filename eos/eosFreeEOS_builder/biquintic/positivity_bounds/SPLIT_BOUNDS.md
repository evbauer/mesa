# Experiment028: declared-reference split extension

```bash
OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 -B mesa/eos/eosFreeEOS_builder/biquintic/positivity_bounds/test_split_bounds.py -v
OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 -B mesa/eos/eosFreeEOS_builder/biquintic/positivity_bounds/bound_split.py
OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 -B mesa/eos/eosFreeEOS_builder/biquintic/positivity_bounds/verify_split_artifacts.py
```

Writes exclusively `freeeos_biquintic_progress/experiments/028_split_bounds/`.
Experiment024 is read-only input: its exact complete cells and frozen material
and full-one-own-use-branch masks are reused and checked against the split core.
No core/runtime/cache/debugger/approval changes or source-truth claim.

The existing `source_constraints(..., split_radiation=True)` remainder jets and
radius2-mean twists retain the declared current CRAD subtraction, including the
known ~4.4 ppm historical radiation remainder. This is a **reference split, not
a verified material EOS or determination of the actual source constant**.

For `Frad=-CRAD/3*exp(2*tau-u+12*ln(10))`, the analytic contribution to
`P_rho=f_u+f_uu` is exactly zero. Its representation certificate therefore bounds
only the remainder polynomial. The analytic `T Cv` contribution is positive
`4*CRAD*exp(2*tau-u+12*ln(10))`. On every adaptive leaf it is bounded below at
`tau_min,u_max` and above at `tau_max,u_min`, and added to the corresponding
remainder polynomial bound. Even a negative remainder can yield positive Cv.

Two status modes are retained, not conflated:

- **Representation**: polynomial-plus-analytic bounds with the same conservative
  polynomial guards as024 and a longdouble exponential-value guard.
- **Restored-core guarded**: additionally budgets floating radiation restoration
  and cancellation in `evaluate(..., add_radiation=True)` and `thermodynamics`.
  The per-leaf restoration guard is `4096*eps_float64*max_rad_Cv*(20/12)` for
  Cv numerator or `*(2/12)` for P_rho, using absolute analytic native-derivative
  contributions. Guard-only unresolved cases are explicitly counted.

These are **numerical guards, not directed-rounded exponential interval
enclosures, rigorous IEEE certificates, physical stability or source accuracy**.
All cells receive four actual-coordinate restored-core sample comparisons;
numeric guard exceedances, zero/negative samples, and representation-positive
versus core-nonpositive samples are explicitly counted. Samples calibrate the
implementation, never replace full-cell Bernstein plus monotonicity bounds.

Each mode has depth <=4, <=341 visited nodes per cell/field and <=13 live DFS
nodes; <=150000 cells, batches <=4096. Both statuses and margins are archived,
with all actual core-negative witnesses and unresolved cells. JSON records
source/core/sidecar SHA, exact baseline024 hashes, worst bounds and guards,
calibration counts and explicit `floating_proof_claim: false`.
The independent artifact verifier recounts every cell in both modes, reconstructs
64 spread-out split root bounds, replays every witness with `add_radiation=True`,
and checks that experiment024 and source/core checksums remain unchanged.

Lead integration uses `028_split_bounds/SUMMARY.md` and `results.json`, keeping
the material/full-one-branch primary row separate from broad and partial domains,
and representation positivity separate from restored-core guarded positivity.
The extension does not claim Pgas/E/S positivity or change the main report.
