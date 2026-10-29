# FreeEOS potential interpolation experiments

This is a table-only, standalone prototype on `biquintic_FreeEOS`, not a
production MESA EOS backend. No selection/blending defaults, distributed tables,
installed caches, WDM_AV mixing, or FreeEOS low-temperature bug are changed.

The staged roadmap is `../BIQUINTIC_PLAN.md`. The separate living report is
`../../../../freeeos_biquintic_progress/freeeos_biquintic_progress.pdf` relative
to this directory. Its experiment outputs record numerical evidence and input
fingerprints; the initial scientific assessment remains separate.

## Scientific conventions

- Source logs are decimal; interpolant derivatives are in `u=lnQ,tau=lnT`.
  Physical coordinates satisfy `x=lnrho=u+2tau-12ln10` and `t=lnT=tau`.
  Physical temperature derivatives therefore use `D_tau-2D_u`.
- Only printed pressure is gas-only; E/S/Cv/chi remain total. Reconstruct total
  P with the declared current-MESA radiation constant. Use `f=e-Ts` and
  six first/second constraints; three higher mixed twists are estimated locally.
- Require native rows, finite supplied/redundant columns, positive source Cv
  and chiRho, logT>=3.2, density-coordinate consistency, and complete accepted
  centered stencils through both twist stages. No one-sided or fallback path.
- Total-potential and analytic-radiation-reference-split representations are
  experimental alternatives. Split mode subtracts declared analytic Frad from source
  constraints and restores it analytically; it does not change mixture physics
  or establish the historical FreeEOS radiation constant. The table-implied
  pressure-strip constant matches current MESA CRAD, while tight radiation
  derivative proxies retain a ~4.4ppm extra radiation-like contribution. The
  split remainder is not a verified pure-material EOS; preserve that residual.
- `gamma3` means Gamma3; `gamma3_minus1` is a separate internal quantity.
  Thermodynamic partials suffix `_x`/`_t` mean physical natural-log derivatives.
- A common potential enforces identities, not physical accuracy. Held-out rows
  test archived-field fidelity at coarser resolution. Independent arbitrary-query
  source accuracy remains unmeasured in this table-only phase.
- Shared-node biquintics are C2, not generally C3. Response values are continuous
  away from singularities; response Jacobians can jump. Better twists alone
  cannot alter normal-third jumps at tangential source knots.
- Source headers validate format, composition, counts, positive spacing, and
  displayed endpoints against count/minimum/step. The composition sidecar has
  an explicit solar-Z-only, independently audited display-endpoint resolution;
  it never changes rows or authorizes trimming/resampling.
- Component own-use cuts are not final EOS fractions. The scaled-solar
  `mesa_probe/` measures actual default selection without installing the
  candidate. Pure FreeEOS still contains some large response-Jacobian jumps,
  although many dense tails are Skye-only. Centered stencil support also leaves
  a clean low-temperature margin up to logT=3.28 that is distinct from the bug.

## Reproduction

From the workspace root, using local Python with numpy/scipy/matplotlib:

```bash
python3 -m unittest discover -s mesa/eos/eosFreeEOS_builder/biquintic -p 'test_*.py' -v
python3 mesa/eos/eosFreeEOS_builder/biquintic/run_audit.py --output freeeos_biquintic_progress/experiments/001_native_prototype
python3 mesa/eos/eosFreeEOS_builder/biquintic/validate_table.py --output freeeos_biquintic_progress/experiments/002_table_validation
python3 mesa/eos/eosFreeEOS_builder/biquintic/crosscheck_fortran.py --output freeeos_biquintic_progress/experiments/003_fortran_crosscheck
python3 mesa/eos/eosFreeEOS_builder/biquintic/face_obstacles.py --output freeeos_biquintic_progress/experiments/007_face_obstacles
```

`crosscheck_fortran.py` requires the streaming driver built using the compatible
local SDK and existing MESA libraries; see `test_fortran/README.md`. No command
here regenerates installed MESA caches or rebuilds its EOS module.

`export_cache.py --output PATH.cache` emits a metadata-bearing, plain-text
`FREEEOS_BIQUINTIC_V1` cache for a standalone adapter. It is not interchangeable
with any MESA cache. The default uses radius2/mean and analytic radiation split;
use `--total-potential` for the total-polynomial control. A companion JSON records
source/preprocessing/cache fingerprints, derivative route, and exclusions.
Generated caches must live in experiment-owned, ignored output directories.

## Implementation map

- `freeeos_potential.py`: strict preprocessing, Hermite evaluation, physical
  thermodynamics and their log-coordinate Jacobians.
- `run_audit.py`: original saved queries, paired consistency, twist routes,
  finite-difference returned-function checks and redundant-column closure.
- `validate_table.py`: archived-node reproduction, coarsening holdouts,
  real-table faces and analytic-radiation conditioning controls.
- `crosscheck_fortran.py` and `test_fortran/`: independent real/AD Fortran
  versus Python cell checks, analytic tests and invalid inputs.
- `face_obstacles.py`: twist-dependent interior-face versus fixed knot-face
  normal-third jumps.
- `smooth_curvature.py` and `smooth_patch.py`: restricted one-dimensional
  curvature reconciliation and C4 constrained-potential feasibility controls,
  not production representations or source accuracy estimates.
- `reference/`: table metadata/convention audit; direct-source scaffold deferred.
- `review/`: independent review snapshot and rounding-sensitivity work.
- `table_fortran/`, `export_cache.py`, `validate_adapter.py`: strict standalone
  metadata cache, physical-coordinate real/AD evaluation and table crosschecks.
- `c3_control/`, `smooth_validation/`: biseptic and regional C4 controls with
  common heldouts, smoothness/fidelity tradeoffs, and positivity diagnostics.
- `composition_control/`: bounded archived-X holdout, gauge checks, and explicit
  audited endpoint metadata; no runtime composition support.
- `mesa_probe/`, `active_faces.py`, `legacy_face_control.py`, `coverage_gap.py`:
  actual scaled-solar EOS fractions, active face tails, shrinking legacy gaps,
  and rejected-stencil coverage. Private caches reuse installed binary links.
- `inversion_control/`: bounded fixed-density energy self-inversion, not stellar
  solver validation and not an independent EOS accuracy test.

Keep generated data out of code commits. Keep future MESA changes separate
from the report repository and pre-existing plotter edits. Never push.
