# Fixed-density thermal consumer control

This standalone test solves a prescribed, temperature-independent heating or
cooling backward-Euler step at fixed density and the audited X=.70, Z=.02
mixture. The physical residual is

`R(y) = y - y0 - dt*H/(T*Cv)`, where `y=ln(T)`, `dt=1 s`, and
`R_y = 1 + dt*H/(T*Cv)*(1+Cv_t/Cv)`.

The candidate uses the unchanged `freeeos_potential.py` kernel with radius-2
mean twists, analytic radiation split/restoration, and centered source floor
3.2. The legacy control reads the installed, independent-Cv spline coefficients
with the existing `audit_freeeos.py` evaluator. A constant-Cv smooth model
checks the solver and isolates backward-Euler temporal error. Neither table
model is routed through MESA's runtime EOS API.

The 24 exact archived face cases come only from
`freeeos_biquintic_progress/experiments/038_solver_manifest/queries_v2.json`.
Its V2 SHA-256 is fixed in the runner. V1's X/Z labels were reversed; its cases
were unchanged, but V1 is rejected by the canonical runner. The MESA fraction
of FreeEOS is measured at the frozen seed only; neither native-table validity
nor these solver results classify the whole path.

From the workspace root on the allocated compute node:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 -B -m unittest discover -s mesa/eos/eosFreeEOS_builder/biquintic/thermal_consumer -p 'test_*.py' -v
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=freeeos_biquintic_progress/experiments/039_thermal_consumer/.mplconfig python3 -B mesa/eos/eosFreeEOS_builder/biquintic/thermal_consumer/run_control.py
```

The regression runner also discovers `test_solver.py` automatically. Outputs
are experiment-owned in `experiments/039_thermal_consumer/`: `results.json`,
`trials.jsonl` with full iteration/residual histories, `jacobian_checks.json`,
`path_queries.csv`, `span_queries.csv`, and figures. The CSV path schema is
`case_id,trial_id,model,stage,logRho,logT`; `stage` distinguishes `initial`,
`solver_start`, `accepted_iterate`, `reference_root`, and rejected proposal
reasons. `span_queries.csv` is a separate finite, nine-point straight-span
sample. Runtime probing these rows can classify sampled states, but cannot
prove the fraction along an entire continuous path.

`common` trials share the same physical `y0`, `H`, bracket and Newton start
across candidate, legacy and smooth models. `off_face_control` moves the
physical start to the same side of the selected face and reports any other
native crossings. `manufactured` trials choose model-specific forcing for
exact and near-face targets. `restart_guess` changes only the Newton start;
`tolerance_sweep` changes only the nonlinear stopping threshold. The latter
two are supplemental controls, not extra members of the frozen sample.

Brent's method solves each model's own residual independently. Positive sampled
Jacobian values and sign-change counts are diagnostics, not a proof of global
root uniqueness between checked points. Root agreement tests numerical
convergence, not source-table accuracy. Every failure appears with a status in
`trials.jsonl` and the aggregate status count.
