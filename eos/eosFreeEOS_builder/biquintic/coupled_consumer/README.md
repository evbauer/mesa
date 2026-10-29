# Coupled pressure plus standard eps-grav standalone control

This fixed-composition table-only consumer uses the frozen V2 24-face manifest.
It does not invoke the MESA EOS API, a stellar residual, a burner, or any
runtime bridge. The candidate is the verified radius-2 mean split-reference
potential. The legacy control independently interpolates `lnPgas`, `Cv`,
`Cp`, `chiRho`, `chiT`, and `grad_ad` from the installed cache; it differentiates the
`lnPgas` spline for **actual total-pressure partials** after analytic radiation
restoration. `grad_ad` is zero-based cache slot 6, corresponding to
`eos_def.f90`'s one-based `i_grad_ad=7`; `eosdt_load_tables.f90` loads it
from `jgrad_ad=15` among EOS values (zero-based raw-row index 15 after the
leading logT field), and `eosdt_eval.f90` returns its own
interpolated value and natural-log partials. The cache-slot and source-node
mapping have a direct regression test. Pressure derivatives are not replaced
by independently interpolated response columns. This table adapter still
does not claim to reproduce final MESA selection or blending.

The residual is

```text
R1 = ln(Ptotal/Ptarget)
R2 = (epsgrav + H)/energy_scale
epsgrav = -T Cp [(1-grad_ad chiT)(t-t0) - grad_ad chiRho(x-x0)]
x = ln(rho), t = ln(T), dt = 1 second.
```

`gravity_control.gravity_residual` supplies the exact analytic value and
state Jacobian of the eps-grav term. The common forcing uses the geometric
mean of model pressures and arithmetic mean of negative model eps-grav values
at the same symmetric nominal target. The initial state and energy scale
`T0*sqrt(Cv_candidate*Cv_legacy)` are shared. Model-specific manufactured
controls use each model's own nominal target. An ideal gas with constant `Cv`
checks the solver, reference, and all four Jacobian entries.

For each frozen case, deterministic thermal scales `1e-4,1e-3,1e-2`, both
signs, and density/temperature increment ratios `+1,-1` place nominal start
and target on opposite sides of the seed face. Both models use fixed
natural-log boxes of half-width `.08` in `x` and `.015` in `t`, residual
infinity-norm tolerance `1e-10`, at most 25 Newton iterations, and the same
Armijo-style damping and invalid-trial accounting. The independent reference
brackets pressure in `x` at each fixed `t`, then brackets the reduced eps-grav
residual in `t`. Eleven samples inspect sign changes and reduced slopes;
sampled monotonicity is not a global uniqueness proof. Every bracket failure
is retained rather than excluded or retuned.

Run on the inherited allocated compute node with at most two numerical
threads, from the workspace root:

```sh
env OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2 python3 -B -m unittest discover -s mesa/eos/eosFreeEOS_builder/biquintic/coupled_consumer -p 'test_*.py' -v
env OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2 MPLCONFIGDIR=/tmp/freeeos_coupled_mpl python3 -B mesa/eos/eosFreeEOS_builder/biquintic/coupled_consumer/run_control.py
```

`--pilot` runs the predeclared four-case subset in an experiment-owned `pilot`
directory. The full run writes `trials.jsonl`, `jacobian_checks.json`,
`results.json`, a focused PDF figure, and `path_queries.csv` and
`span_queries.csv`. The CSV schema matches the scalar consumer and is for
the lead's independent MESA path probe. Finite samples, including straight
spans, do not prove continuous path purity or final blended-EOS behavior.

The unit test file matches `test_*.py` in this dedicated directory and is
automatically discovered by `biquintic/run_validation.py`.
