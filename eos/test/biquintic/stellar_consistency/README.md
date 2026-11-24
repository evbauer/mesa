# Matched stellar-state consistency replay (056)

This diagnostic evaluates both fixed-composition backends through the frozen
050 component/final APIs at identical rho/T coordinates from explicitly named
052 candidate profiles. It neither evolves a star nor compares two different
stellar profiles as if they shared the same state. All named profile zones are
included; counts/quantiles are unweighted, not stellar-population statistics.

The three returned-function identities use total pressure, the established
radiation convention, and numerical state partials returned by MESA. The fourth
metric is `e_t/(T*Cv)-1`. Signed sum-of-absolute-terms-scaled balances accompany
the ratios. These are internal closure checks, not off-node physical accuracy.

Run only on the inherited allocated node; the reused client checks placement,
private cache links, successful response finiteness, fractions and callbacks.
The audit pins input/code/executable hashes, requires profile/replay fraction
agreement and unchanged installed EOS/cache inventory, and rejects unsuccessful
queries. It refuses to overwrite a completed audit. The exact frozen source
table/potential and composition are those of 047/050; no new source computation.

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 -B -m unittest discover -s mesa/eos/test/biquintic/stellar_consistency -v
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 -B mesa/eos/test/biquintic/stellar_consistency/audit.py --profiles freeeos_biquintic_progress/experiments/052_stellar_api/runs/candidate_100/LOGS/profile1.data freeeos_biquintic_progress/experiments/052_stellar_api/runs/candidate_100/LOGS/profile8.data freeeos_biquintic_progress/experiments/052_stellar_api/runs/candidate_100/LOGS/profile100.data --output freeeos_biquintic_progress/experiments/056_stellar_consistency/reproduction
python3 -B mesa/eos/test/biquintic/stellar_consistency/plot.py freeeos_biquintic_progress/experiments/056_stellar_consistency/hundred_step/summary.json
```

Three tests verify analytic ideal-gas-plus-radiation closure, detection of an
independently perturbed energy gradient, and the MESA profile parser. The plot
reads saved results only. Dashed/solid lines mean component/final, orange/blue
mean legacy/potential; the lower panel supplies the final routing weight.
Connected lines guide the eye through discrete zone samples; no additional
EOS interpolation is used for the diagnostics.
