# FreeEOS potential interpolation

This existing-table project on `biquintic_FreeEOS` implements a default-off
native MESA backend deriving thermodynamics and numerical partials from one
specific Helmholtz free energy. It retains the existing composition selection
and weights, including configured C/O variants. The original FreeEOS source
rows and EOS selection/blending defaults are unchanged. This remains an
experimental backend with support and stellar-convergence limitations.

The staged roadmap is `../BIQUINTIC_PLAN.md`. The separate living report is
`../../../../freeeos_biquintic_progress/freeeos_biquintic_progress.pdf` relative
to this directory. Its experiment outputs record numerical evidence and input
fingerprints; the initial scientific assessment remains separate.

Current implementation (23 November 2026): all 89 sources have compact HDF5
supplements. At loading, the backend reconstructs six lower constraints from
the distributed source rows and reads three mixed derivatives through FORUM.
All nine constraints remain in memory for queries. The archived full-text
representation remains available for explicit numerical comparisons.
Broad stellar-suite evidence predates this storage refactor: current per-case
outcomes include 86 successful tightened-routing evolutions, two skips, two
minimum-timestep failures and one user-stopped case. This is not a fully passing
replacement under default routing. Current decisions and validation status live in
`../../../../freeeos_biquintic_progress/NEXT_PHASE_PLAN.md`.

## Enabling the native backend

Hydrate the branch's EOS data archives through the usual Git-LFS workflow and
build/install MESA normally. `eos/Makefile` extracts
`eosFreeEOS_biquintic_data.tar.xz` after the original FreeEOS archive, overlaying
the 89 supplements into `mesa_data_dir/eosFreeEOS_data/`. In an EOS inlist use:

```fortran
&eos
   FreeEOS_use_biquintic_interpolation = .true.
/
```

The control defaults to false. `FreeEOS_biquintic_table_dir = ''` uses the
installed data directory; an explicit override is available for development.
Sources at the native MESA call site always come from the MESA data directory.
An invoked unsupported query fails explicitly; it does not silently substitute
legacy FreeEOS. Output-level blends with other EOSs retain their existing form.
Legacy molecular-weight, electron and degeneracy auxiliary outputs are retained.

The supplement adds 188.845 MiB to the original 1321.928 MiB compressed FreeEOS
archive. See [supplemental_hdf5/README.md](supplemental_hdf5/README.md) and
[FORMAT.md](supplemental_hdf5/FORMAT.md) for source-based generation, portable
verification/packaging commands, metadata and FORUM reader limitations.

## Runtime behavior and validation

Native calls retain legacy auxiliary quantities and analytic partials, but
evaluate only the six legacy columns through `eta`, including all inputs to
existing legacy guards. Toggle-off calls retain full legacy interpolation.
Each thread observes a successfully loaded native table under the shared load
lock before recording its index. Subsequent lookups read that immutable table
without the shared lock, while checking directory, source, grid and suffix.
Misses and failed loads retain the locked path and original error semantics.

Matched state replays show about 9–10 times the pre-optimization native warm
final-EOS throughput at 28 threads, and 1.16–1.20 times at one thread. Complete
return/partial, ordinary/C/O, concurrent cold/warm and strict-error parity checks
pass. Short C/O stellar restarts at one and 28 threads preserve all recorded
histories and final-model hashes. Startup reconstruction/authentication remains
once per table per process; neither these warm timings nor the short restarts
measure a new full stellar-suite speed ratio. Detailed task evidence is in
`../../../../freeeos_biquintic_progress/progress_reports/2026-11-23_runtime_profile.pdf`.

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
- FreeEOS-only diagnostics are not final EOS fractions. The scaled-solar
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
- `c3_spline_control/`: double-knot quintic C3 rank audit and explicit normal-
  fourth roughness selector; unsuccessful full fits remain recorded/rejected.
- `regional_c3/`: eighteen declared coarse fits paired with archived C4;
  no uniform Cv tail improvement and no runtime representation selected.
- `first_constraint_control/`: one coarse C3 fit with hard E/TS/P constraints
  and soft second responses; exact retained nodes do not ensure better holdouts.
- `c3_gauge_control/`: known artificial a+bT offsets reveal smoother selector/
  scaling/representation sensitivity; analytic removal is not source gauge inference.
- `edge_bounds/`: guarded whole-cell signs for newly supported shifted-native
  cells, including explicit3.1-only cells; geometry is not actual EOS fraction.
- `composition_control/`: bounded archived-X holdout, gauge checks, and explicit
  audited endpoint metadata; no runtime composition support.
- `mesa_probe/`, `active_faces.py`, `legacy_face_control.py`, `coverage_gap.py`:
  actual scaled-solar EOS fractions, active face tails, shrinking legacy gaps,
  and rejected-stencil coverage. Private caches reuse installed binary links.
- `inversion_control/`: bounded fixed-density energy self-inversion, not stellar
  solver validation and not an independent EOS accuracy test.
- `gauge_control.py`: finite-grid constant/linear-temperature gauge covariance,
  not inferred historical gauge errors or a change to source E/S conventions.
- `gravity_control.py`: offline standard-eps-grav EOS-term Jacobian at exact
  pure-FreeEOS faces, with fixed previous-state increments; not a stellar solve.
- `runtime_adapter/`: authenticated fixed-mixture 27-slot return contract,
  physical state partials, complete hybrid abundance rows and strict guards.
- `runtime_bridge/`: default-off callback, private EOS/stellar builds, actual
  API checks and bounded serial client; no general production interface claim.
- `runtime_consumers/`: frozen scalar/coupled trials through the actual APIs,
  including retained tight-tolerance failures and separate precision controls.
- `helm_control/`: native-coordinate HELM face benchmark, not an acceptance
  threshold transferred from another EOS or routing population.
- `stellar_control/`, `stellar_candidate/`: original-backend hybrid controls
  and matched fixed-composition stellar integration with private dependencies.
- `stellar_consistency/`: matched saved stellar-state component/final identities,
  not another evolution run or independent physical accuracy.
- `arithmetic_control/`: isolated same-polynomial real-arithmetic improvement,
  with unchanged AD limitations and the original consumer failures preserved.
- `arithmetic_ad_control/`: third-order AD jets composed from the stable real
  derivatives; no change to the cell polynomial, continuity class or stellar backend.
- `stellar_timestep/`: same-seed matched-age coarse/half exploratory comparison;
  strict original gate failures and profile-retention limitations stay explicit.
- `run_validation.py`: installed-file protections, checksums and available
  Python/Fortran/probe/private-API regression; requires compatible saved binaries.

Keep experiment-generated caches, binaries and plots out of code commits.
The supplemental distribution archive is the intentional Git-LFS data artifact.
Keep future MESA changes separate from the report repository and pre-existing
plotter edits. Never push without user authorization.
