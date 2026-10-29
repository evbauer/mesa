# Bounded multislice smoother validation

Read-only use of the lead degree-5 C4 `smooth_patch.py` soft-fit APIs, never a runtime EOS change or a source oracle.

From the workspace root:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 mesa/eos/eosFreeEOS_builder/biquintic/smooth_validation/test_validation.py -v
PYTHONDONTWRITEBYTECODE=1 python3 mesa/eos/eosFreeEOS_builder/biquintic/smooth_validation/validate_multislice.py
PYTHONDONTWRITEBYTECODE=1 python3 mesa/eos/eosFreeEOS_builder/biquintic/smooth_validation/active_exact_faces.py
```

Default: four spawned patch workers, eight OpenBLAS/OMP/MKL threads each, 5,000 strict off-node queries per patch, 128 sampled internal faces per orientation/fit. Outputs are exclusively `freeeos_biquintic_progress/experiments/015_smooth_validation/`. `--workers 1`, `--patches NAME`, `--queries N`, and `--faces N` support bounded reruns without changing fit/source APIs.

Selection searches within declared physical boxes for the largest all-native/valid/material patch with clean full/coarse biquintic donors, globally even starts, odd 17–35 dimensions, and at least 144 odd/odd holdouts. Any shrink or failure is explicit; no flagged/low-temperature/fallback rows enter a fit. The smoother is patch-local, not a global fallback-influenced spline.

Compare weights through archived-node reconciliation, coarse-only held-out fidelity, boundary/interior tails, positivity, and numerical C4 limit checks—not just smoother curves. The global coarse biquintic control has outside-patch donors; a separate common subset uses patch-only coarse donors and identical retained source rows. Gamma3 is the corrected full Gamma3.

Simple interior knots guarantee C4 potential continuity, not C5, source accuracy, stability, gauge alignment, or continuity across unblended patch ends. Finite one-sided steps measure an O(step) gap and roundoff-sensitive extrapolated intercept, not an exact zero-step jump. Smoothing is not necessarily truth; all original/source/kernel files remain untouched, and no commits or pushes are made.

The exact-face follow-up reuses saved coefficients at three runtime-selected faces from experiment `014_active_domain`; it does not refit, rerun runtime probes, or extrapolate active fractions to the smoother. The selected active patch contains q=.77/logT=6.75647698345, q=.8/logT=6.765258774, and q=.88188515665/logT=6.56. Region-level interpretations and unresolved gates are in experiment `015_smooth_validation/INTERPRETATION.md`.
