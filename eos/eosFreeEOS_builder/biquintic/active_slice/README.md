# Actual active-FreeEOS slice, no refits

This exclusive sidecar writes only experiment `023_active_slice`. The lead
integrates its PNG/PDF into the main report; this worker does not edit the
report, core EOS, other workers, installed sources/caches, or archived fits.

1001 queries cover logQ=.50..1.10 at fixed logT=6.756476983450105, strictly
off native T knots. The .77 pure and .80 blended Q faces are included. All
native Q knots are marked vertically. The endpoints and all queries are
strictly inside the archived 015 `active_chiT` patch, not patch boundaries.

Three panels show:

- Cv for radius2/mean exact-constraint BQ, independently interpolated legacy Cv
  from the existing baseline cache, and archived C4 thermo_0.01/thermo_0.1.
- Physical Cv_t = dCv/dlnT at fixed rho for each **offline component**, not the
  derivative along this fixed-T density sweep and not runtime-blended MESA Cv_t.
- Actual default-MESA fracFreeEOS/fracSkye from **one** existing `mesa_probe`
  executable invocation streaming all 1001 points. `use_WDM_AV=.false.` is the
  existing probe's sole inlist override. No guessed own-use mask replaces it.

The saved 015 full coefficient arrays are evaluated directly by `evaluate_fit`.
They are total potentials, so radiation is not added again. No fitting function
is invoked. Exact BQ uses the same total-potential source convention as 015.
Cv/Cv_t reproduction at 256 archived queries guards evaluator compatibility
after subsequent reader-validation changes. Source/cache/archive/code SHA-256s
are recorded and reverified; no installed cache is rebuilt.

Runtime shading/labels use measured fractions only. Inactive component curves
remain visible as explicitly inactive context, never as selected MESA values.
Region boundaries are estimated between neighboring samples, not exact blend
boundaries. If there is no inactive region, none is invented.

Gas pressure and Cv positivity and finite Cv_t guard every curve. Fixed-rho
finite differences shift logQ by -2h/ln10 and logT by h/ln10. Source-Q knot
centers are excluded from that check to avoid averaging BQ's discontinuous
normal-third derivatives. These checks are returned-function diagnostics,
not fresh physical accuracy. Local C4 inside this patch implies no global EOS
C4 continuity, correct derivative levels, mixture accuracy or runtime readiness.
Historical gauge/metal-pattern assumptions are still not independently verified.

The probe uses experiment-023-owned cache symlinks and temporary directories;
installed cache inventories and symlink integrity are checked before/after.
The fixed-output `mesa_probe/run_probe.py` wrapper is **not** called because it
writes experiment 014. Its already-built executable is reused without rebuild.
Repeated plotting reuses the integrity-checked saved probe, not a second call.

From the workspace root:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s mesa/eos/eosFreeEOS_builder/biquintic/active_slice -p 'test_*.py' -v
PYTHONDONTWRITEBYTECODE=1 python3 mesa/eos/eosFreeEOS_builder/biquintic/active_slice/run_slice.py
```

Outputs: `active_freeeos_slice.png`, `active_freeeos_slice.pdf`, `queries.npz`,
`results.json`, `README.md`, and single-call probe TSV/input/stderr/metadata.
Curve colors are unconstrained; no logarithmic consistency map is produced.
