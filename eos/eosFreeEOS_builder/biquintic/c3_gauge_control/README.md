# 034 known-gauge sensitivity of the coarse C3 selector

One bounded **table-only** gauge-sensitivity sidecar. No inferred original
FreeEOS gauge, new source accuracy, selected runtime, new tables/caches or
changes to frozen031/033, core or other workers. The artificial gauge is known
exactly; it is not estimated by fitting interpolation errors.

## Reproduce

```sh
python3 -B mesa/eos/eosFreeEOS_builder/biquintic/c3_gauge_control/test_gauge.py
python3 -B mesa/eos/eosFreeEOS_builder/biquintic/c3_gauge_control/run_control.py --overwrite
```

Existing Python/NumPy/SciPy/threadpoolctl only; 16 BLAS threads and frozen031
`solve_regularized` API. Helper SHA-256 is checked before and after:
`ba30d62f864e094bb52ad0b83f87554b8fff199c4146f1c1e1a6fc297b35edfd`.
All writes are exclusive to this source directory and
`freeeos_biquintic_progress/experiments/034_c3_gauge/`: README, JSON,
`paired_gauges.npz`, tests and validation logs. No changes to031/033 outputs.
The existing solar table hash is pinned to
`e6ca7ce27bee8197bb257a8c4886a06366e6951e58958fbed3735c877ee68502`.
Source/archive/command/version hashes and process RSS/time are recorded.

Exactly five **coarse** fits: baseline; two known shifts fitted as total
potential; two same shifts subtracted analytically before fitting and restored
analytically afterward. No full fits, new solve after local19:55 or parameter
search. Each fit uses1,064 unknowns, full augmented QR, the existing90-second
guard, lambda .01 and response importance .1. Total runtime about2.5seconds;
whole-process peak RSS about290MiB. No network/build/permission request needed.

## Same input/query conventions

Use only the exact archived015 active geometry patch: 18x13 retained native
coarse source nodes, the same204 archived odd/odd holdouts in saved order, and
known current-CRAD total-pressure/printed E/S conventions. The `active_chiT`
label is geometry, not measured MESA fractions. Archive033 heldouts are checked
identical and its lambda .01 values are saved without refitting; archive031
baseline coefficients are likewise reused for the reproduction check.

Let `f_new=f+a+b*T`, where the known test constants are

```text
a = median(printed retained E) = 1.2181677095425765e15 erg/g
b = median(printed retained S) = 9.747893323784697e8 erg/g/K
```

The two shifts use `(a,b)` and `(a,100*b)`. The constant a is the same in both;
there is no independent a-only fit or isolation of a's scaling effect.
Physical consequences, exactly, are

```text
E_new = E+a
S_new = S-b
P, Pgas, Cv, Cp, chiRho, chiT and their response partials unchanged
```

With native coordinates `u=ln10*logQ`, `tau=lnT`, add `a+b*T` to the `(0,0)`
low jet and `b*T` to `(0,1)`/`(0,2)`; all u-containing jets are unchanged.
At queries the exact analytic term is also differentiated through third order.
Total-fit comparisons subtract this **known gauge** from the returned native
potential/partials before reconstructing thermodynamics. Remainder-fit comparisons
use the fitted original-gauge remainder; raw shifted output restores the term
analytically. The E-a/S+b fieldwise removal crosscheck is archived.

## Why the selector is gauge sensitive

Frozen031 penalizes normalized integrated native normal fourth derivatives:

```text
R = [h_u^8 integral(f_uuuu^2) + h_tau^8 integral(f_tautautautau^2)]
    / (patch_area*median(abs(f))^2)
```

But `d_tau^4(b*T)=b*T`, so this penalty is **not** invariant to an entropy
gauge term with no physical Cv/pressure effect. A constant a has zero numerator
roughness, but changes potential scaling and E weights. Relative E/TS target
scales also change under these offsets. In addition, exponential b*T is not
exactly in a degree-five piecewise polynomial spline space; fitting it as part
of the total introduces approximation absent from exact analytic restoration.
The audit demonstrates **combined selector/scaling/representation sensitivity**,
not a numerical attribution exclusively to the roughness penalty.

The analytic pure-gauge penalty at a given fit scale s is

```text
R(b*T) = h_tau^8*b^2*(T_high^2-T_low^2)/(2*delta_tau*s^2)
```

It is saved for context, not subtracted from the fit objective: squared
roughness also contains a cross term with the original potential. No other
physics is fitted away. Scales are recomputed by the frozen API: baseline
4.087e15, total b1 2.337e15, total b100 5.327e17. The reference-subtracted
routes retain baseline scaling to source-subtraction roundoff.

The reference route removes exact known `a+b*T` from source low jets **before
all fit scales and weights are computed**, then restores analytically. This
defines a controlled reference convention, not automatic gauge invariance of
the original objective or discovery of the upstream gauge. For100b the input
subtraction roundtrip reaches2.33e-14 relative; that is recorded rather than
claimed exact bitwise recovery.

## Results and gates

Every fit has augmented rank1,064/1,064 and passes the frozen normalized
stationarity guard <=1e-10. The100b total fit has6.35e-11 stationarity, closer
to the guard than baseline6.37e-13; no conditioning/accuracy claim beyond that
gate. Baseline coefficients reproduce archived031 **exactly** in this run.

Known-gauge-removed relative departures from the baseline on the same204 targets:

| Route | Pgas p99 | Cv p99 | Cv max | chiRho p99 | chiT p99 |
| --- | ---: | ---: | ---: | ---: | ---: |
| Total a+bT | 3.451e-6 | 6.155e-4 | 1.113e-3 | 5.258e-5 | 1.784e-4 |
| Total a+100bT | 3.280e-4 | 3.782e-2 | 5.921e-2 | 5.076e-3 | 1.732e-2 |
| Subtract/restore a+bT | 4.629e-13 | 3.203e-11 | 4.803e-11 | 5.069e-12 | 1.343e-11 |
| Subtract/restore a+100bT | 1.094e-12 | 4.298e-11 | 5.260e-11 | 1.057e-11 | 2.152e-11 |

These are **gauge sensitivity/model departures**, not fresh-source errors.
The JSON also preserves printed-reference holdout fidelity and scaled physical
response-Jacobian departures. Values divide by abs(baseline field); Jacobians
divide by max(1,abs(baseline field),abs(baseline partial)). All states/coefficients
and paired IDs are archived, with no interpretation of raw gauge units as truth.

All recovered original-gauge sampled node/holdout fields are finite and positive.
Raw shifted S is nonpositive at102 b1 targets and all204 b100 targets, as expected
from a gauge offset. These signs are not invariant physical positivity failures;
raw signs and recovered status are reported separately. Sampled positivity does
not bound unsampled interiors or patch ends.

Three tests pass: exact a+bT physical field/response transformation; pure-tau
low-jet shift and source roundtrip; positive analytic normal-fourth penalty of
bT and its b-squared scaling. This gauge control exposes a reference-choice
requirement for the experimental roughness selector. It **does not establish
original FreeEOS gauge metadata, offnode source accuracy or runtime suitability**.
