# Matched FreeEOS compatibility control

This standalone diagnostic reads the corrected frozen 24-case V2 solver
manifest. It does not modify the source table, numerical kernel, cache, MESA
runtime route, or query selection. The source composition is X=.70, Z=.02.
V2 corrects reversed X/Z metadata in the superseded V1 manifest; the runner
checks exact identity of all 24 case records across versions.

Run on the inherited allocated compute node, with at most two numerical threads:

```sh
env OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2 MPLCONFIGDIR=/tmp/freeeos_compat_mpl python3 -B mesa/eos/eosFreeEOS_builder/biquintic/compatibility_control/test_compatibility.py
env OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2 MPLCONFIGDIR=/tmp/freeeos_compat_mpl python3 -B mesa/eos/eosFreeEOS_builder/biquintic/compatibility_control/run_control.py
```

The runner writes only `freeeos_biquintic_progress/experiments/040_matched_compatibility/`.
It evaluates Appendix A source-knot-line compatibility and direct one-sided
off-knot face jumps. The physical derivative relation uses
`d/d ln T|rho = d/d tau - 2 d/d u`. It checks analytic radiation split versus
unsplit interpolation, centered radius 1/2/3 and radius-2 differentiation
routes. Radius-2 mean with analytic radiation restoration is the candidate.
All variants retain exact supplied low-order nodal data. A zero-twist algebraic
control only isolates the estimated mixed-nodal contribution; it is not a
candidate interpolant.

Local resolution control follows experiment 021: retained-only grids at
strides 2/4/8, with identical odd/odd archived source knots withheld from all
three grids. Up to four nearby common targets are chosen per frozen face. The
source `Cv` column is the reproduction reference, not independent FreeEOS
accuracy. Stride 1 retains the targets and is not a holdout. Target reuse
across cases is recorded and should not be counted as independent samples.

The four-case printed-rounding control uses deterministic uniform draws within
the source's displayed half-decimal intervals for six consumed fields at
conservative local stencil support. It measures formatting sensitivity only;
it is neither a source-noise model nor a rigorous interval bound.

The figure compares the direct off-node scaled `Cv_t` jump with the largest
source-knot-line jump at the same face, and shows sensitivity to the centered
stencil radius/route. Every selected face is pure FreeEOS at its seed only.
No inference about actual solver-path selection, final blended EOS behavior,
or independent native-resolution source accuracy follows from this control.
