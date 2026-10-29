# FreeEOS biquintic free-energy interpolation: test plan

## Objective and scope

On branch `biquintic_FreeEOS`, test whether a potential-based FreeEOS backend
can improve thermodynamic consistency without degrading physical accuracy,
stability, response-function smoothness, or MESA solver performance.

The scientific roadmap is `../../../freeeos_report/freeeos_interpolation_report.pdf`.
Maintain a separate living progress/results PDF at
`../../../freeeos_biquintic_progress/freeeos_biquintic_progress.pdf`, rebuilding
it after substantive milestones, experiments (including failed tests), and
scope decisions. Keep the initial assessment as background; log revisions,
reproduction settings, quantitative results, limitations, and next decisions
in the progress report rather than treating plans as completed tests.
This is a FreeEOS interpolation experiment, not further WDM_AV mixture
development. OPAL reconstruction and potential blending between different EOS
models are out of scope. Existing defaults and distributed tables remain the
baseline until the staged gates below pass.

The low-temperature FreeEOS table issue is a known, separate MESA bug.
Its expected upstream mitigation is to move the blend so FreeEOS contributes
nothing by logT=3.1. Repairing those tables or changing that blend is not part
of this project and is not a prerequisite for the clean-domain prototype.
The planned upstream change is context, not assumed to be present on this
branch; retain revision-specific selection masks in comparisons.

Start with the installed `mesa-FreeEOS_02z70x.data` mixture at X=0.70, Z=0.02.
Do not substitute the plotter's H/He/C/O mixture for the table's scaled-solar
metal pattern when evaluating a direct FreeEOS reference.

Allocation-resume scope update (28 October 2026): limit the current phase to
existing FreeEOS tabulations. Do not build new tables or pursue a new solver
installation. Use native-node reproduction, held-out-node/coarsening tests,
matched-query identities, stability, twist sensitivity, and cell-face/Jacobian
diagnostics. These tests measure reconstruction fidelity and internal behavior,
not independent full-resolution off-node accuracy. Direct-source validation
below is a deferred recommendation, not a dependency to continue table-only
experiments or a claim that this phase will deliver production readiness.

## Reuse and separation

- Borrow the corner packing, HELM Hermite basis weights, real/automatic-
  differentiation evaluation, and derivative tests in
  `../private/eos_wdm_av_interp.f90` and
  `../eosWDM_AV_builder/test_fortran/test_interp.f90`.
- Use `../private/helm_polynomials.f90` directly. Initially keep the
  FreeEOS adapter separate; do not refactor the working WDM_AV implementation
  as a prerequisite or reinterpret FreeEOS mixtures as isotope tables.
- Adapt preprocessing patterns from
  `../eosWDM_AV_builder/wdm_av_builder/derivs.py`, but apply FreeEOS-specific
  units, radiation conventions, source flags, and sheared-coordinate transforms.
- Use `../private/skye_thermodynamics.f90` as a convention and differentiation
  reference. Its physical-coordinate AD inputs differ from WDM_AV's log inputs;
  do not connect the interfaces without an explicitly tested transformation.
- Derive P, e, s, Cv, Cp, susceptibilities, adiabatic coefficients, and their
  numerical partials from the same potential. Only auxiliary quantities such as
  mu, electron abundance, and eta may remain separately interpolated.

## Stage 0: preserve and reproduce the baseline

1. Record the branch/revision, source-table checksums, compiler/SDK, radiation
   constant, composition, and EOS controls for each experiment.
2. Preserve the existing uncommitted plotter source/documentation changes.
   Stage only this project's changes for future commits; do not bundle those
   pre-existing edits into an implementation commit.
3. Reuse the report's saved query coordinates, not fresh unrelated samples:
   the 12,000 interior queries, 16 stratified queries per native cell,
   actual selection/fraction grid, and dense low-temperature slice.
4. Keep generated data and caches in an experiment-owned directory. Check
   extraction targets before rebuilding: the current EOS Makefile can delete
   installed table/cache directories when extraction targets run. Do not
   regenerate distributed tables as a side effect of a test build.

Gate: the legacy evaluator reproduces saved baseline diagnostics on the
same queries, with differences explained before comparing a new method.

## Stage 1: audit the source and construct a local prototype

Create a preprocessing/audit tool under a new `biquintic/` directory here.
Keep processed data distinct from the existing bicubic interpolation cache.

1. Decode all physical columns and retain native/fallback flags. Exclude
   known defective low-temperature data and stencils crossing them rather
   than repairing archived rows. Start in the original audit's clean interior
   (logT>=3.2 with complete valid native stencils), then validate the intended
   active domain separately. Keep the logRho=-8, logT=3.10 artifact only as
   an exclusion/stencil-isolation diagnostic, not a requirement to repair
   or reproduce valid FreeEOS thermodynamics there.
2. Begin with native-only cells and complete native derivative stencils.
   Compute higher derivatives locally, without allowing global spline
   influence from fallback data. Expand toward boundaries only with audited
   one-sided stencils, fresh source data, or a documented validity exclusion.
3. Work on the existing uniform grid in u=lnQ and tau=lnT. Reconstruct
   f=e-Ts and the six supplied low-order constraints in physical log
   coordinates, then transform all constraints to (u,tau). Test the inverse
   transform as well as the factors of ln(10).
4. Check redundant routes for mixed derivatives before reconciliation.
   Compare neighboring potential values with supplied gradients/curvatures
   using resolution-controlled differences or path integrals; tiny stored
   nodal diagnostics alone do not establish printed-grid integrability.
   Estimate the three missing twists in the interpolation coordinates using
   multiple stencil orders/routes; quantify disagreement rather than silently
   averaging incompatible data or setting the twists to zero.
5. Establish radiation conventions through direct-source checks. First verify
   total-quantity reconstruction; test consistent material-only subtraction
   and analytic radiation restoration before radiation-dominated validation.
   Do not mix gas pressure with unstripped energy/entropy or substitute an
   unverified table-generation radiation constant.
6. Test plain biquintic as the primary candidate, bicubic potential Hermite
   as a diagnostic control, and reference-subtracted biquintic if cancellation
   warrants it. Reference subtraction must act on the whole potential and
   every needed derivative, with an appropriate reference for this mixture.

Gate: explicit provenance/validity masks, internally audited constraints,
and sensitivity results on the declared clean domain. No low-temperature
bug fix is required. Full active-domain/production claims still require tests
up to the actual selection boundaries; after an upstream blend change, update
and rerun those tests without silently reusing old selection masks.

## Stage 2: verify a standalone Fortran interpolator

Build a small test executable against the existing MESA libraries before
changing the runtime FreeEOS route. Separate its private table/derivative
types from isotope rescaling and additive-volume mixing.

Tests:

- A tensor polynomial of degree at most five in each interpolation coordinate:
  exact cardinality, corner ordering, grid-spacing factors, and reproduction.
- Smooth analytic ideal-gas-plus-radiation and nontrivial two-coordinate
  functions: independent exact derivatives, including the (u,tau) shear.
- Real and AD agreement through the derivatives needed for all exported
  quantities and their first numerical partials (third potential derivatives).
- Centered finite-difference step sweeps; use coordinate-correct perturbations
  at fixed physical rho or T, not fixed Q by mistake.
- Interior, near-face, exact-face, exact-node, and outer-boundary probes;
  explicit rejection of invalid cells and out-of-range queries.
- C2 potential continuity across shared faces. Quantify possible jumps in
  third derivatives and hence Jacobian derivatives of response functions;
  do not demand mathematically impossible C3 continuity from generic local
  biquintic Hermite or mistake that limitation for an implementation bug.
- Provisional identity gate: well-scaled first-law/Maxwell balances below
  1e-10 where numerically well conditioned, with conditioning/cancellation
  exclusions declared. Both report diagnostic routes must agree because all
  thermodynamic columns now come from this same potential.

Gate: analytic tests, implementation consistency, and boundary handling pass
independently of the accuracy of the archived FreeEOS data.

## Stage 3: validate against FreeEOS, not just against identities

1. Obtain or locate FreeEOS 2.2.1 and confirm EOS1 options, abundances, atomic
   masses, energy/entropy references, and radiation. Its library was not found
   in the inspected local SDK/build tree; direct-source validation is an
   environment prerequisite, not something to replace with a cache comparison.
2. Establish agreement at native nodes before interpreting off-node differences.
   Use a dedicated query driver that reports failed source calls; never hide
   failures by invoking MESA/HELM fallback in the reference.
3. Evaluate legacy bicubic and the candidates at identical rho/T queries.
   Use all inexpensive interpolant samples and a stratified direct-source
   subset, expanded around defects. Include molecular regions, H/He and metal
   ionization, radiation transitions, cell faces, and actual selected edges.
4. Compare P, e, s, Cv, Cp, chiRho, chiT, Gamma1, Gamma3-1, grad_ad, and
   their relevant numerical partials. Scale zero-crossing f and reference-
   dependent quantities sensibly; do not quote relative f errors near f=0.
5. Separate four reports: physical accuracy, thermodynamic identities,
   positivity/stability, and Jacobian fidelity/smoothness. Report medians,
   p99, maxima, failure counts, coordinates, and provenance.
6. Reuse coarsened held-out tests as a resolution check, not as ground truth
   for full-resolution between-node accuracy. Vary resolution/twist stencils
   and perturb inputs within justified precision bounds to assess sensitivity.

Gate: no unexplained nonpositive Cv or P_rho in the proposed single-phase
active domain; converged estimates and no material accuracy regression versus
legacy at matched source queries. Set physical-error tolerances from the
direct-source baseline and intended application before accepting a candidate;
the 1e-10 identity target is not a physical-accuracy target.

## Stage 4: add an opt-in MESA FreeEOS path

Only after the fixed-table tests pass:

- Add dedicated FreeEOS potential definitions/loading/evaluation modules and
  a default-off control, with build and EOS-control wiring.
- Keep the existing FreeEOS component identity and selection/blend logic.
  Use the existing API to compare standalone component results before testing
  the final blend. A fixed-composition experiment must reject unsupported
  compositions or explicitly report legacy use, not silently misapply X=0.70.
- Introduce a distinct versioned potential cache with composition, source
  checksum, radiation/reference metadata, derivative scheme, and validity mask.
  Do not interpret legacy bicubic cache records as potential data.
- Ensure caller-visible thermodynamic outputs and density/temperature
  derivatives come from the same potential; preserve MESA's gas/total and
  log-return conventions. Negative f is allowed; log-return positivity is a
  separate interface requirement. Do not independently clip Cv to hide failures.
- Extend to other X/Z tables and interpolate potentials with composition-only
  weights before evaluating outputs. Align gauges and metal patterns, and
  validate composition partials. Interpolating separate output logs again would
  undo the consistency gain. This is needed before evolving composition-changing
  stellar models with the experimental path.
- In development, strict tests fail on invalid requested potential cells.
  Any later legacy fallback must be explicit and counted; consistency claims
  apply only where the potential path actually runs.

Gate: default-off behavior matches legacy, cache/API tests pass, and maps
separate pure potential, legacy fallback, and blended contributions. An unchanged
state-dependent EOS blend is not expected to preserve component consistency.

## Stage 5: solver-facing A/B tests and decision

Start with short, matched runs/restarts rather than a broad expensive grid:

- A solar model: envelope ionization, Gamma1/grad_ad and Brunt-Vaisala profile
  smoothness, and asteroseismic sensitivity.
- A low-mass/pre-main-sequence model: molecular and low-temperature conditions.
- A helium- or metal-rich model, including a WD envelope case where FreeEOS
  actually contributes; do not attribute a Skye-dominated core to FreeEOS.

Keep other physics/controls fixed and WDM_AV disabled. Compare EOS wall time,
memory, Newton iterations, retries, accepted timesteps, component fractions,
and evolutionary/profile quantities. Repeat with tighter temporal/spatial
resolution and relevant energy-equation forms. Separate EOS/blend discrepancies
from solver residual and time-discretization errors; rel_run_E_err alone is
not an acceptance metric. Run the existing EOS and WDM_AV regression tests
for any shared infrastructure touched.

Decision: retain an experimental backend unless the source-accuracy,
stability, smoothness, and solver gates all pass. If local biquintic response
partials are too rough, test a smoother constrained potential fit rather than
returning unrelated legacy thermodynamic coefficients. Changing production
defaults is a separate decision.

## Deliverables and commit boundaries

Proposed sequence, with no pushes:

1. Baseline/provenance tools, clean-domain masks, and stencil-isolation tests;
   no repair of the known low-temperature bug or its blend controls.
2. Standalone potential preprocessing and Hermite verification harness.
3. Fixed-composition source-accuracy/positivity/smoothness assessment and figures.
4. Default-off runtime route, metadata/cache and API regression tests.
5. Potential-level composition interpolation and composition tests.
6. Stellar A/B experiments and updated report with a go/no-go recommendation.

Small commits should contain only the changes for their milestone. Generated
large tables, raw exports, and caches should not be committed by default.
No implementation or default change is part of this planning document.

## Recommended first milestone

Start Stages 0-2 on the existing X=0.70, Z=0.02 data, while locating the direct
FreeEOS dependency for Stage 3. The first result should show: an audited clean
domain, a verified standalone biquintic kernel in the native sheared grid,
the same-query blue/orange identity comparison, and twist/positivity/cell-face
diagnostics. That is enough to judge whether the candidate merits runtime
integration without prematurely modifying EOS selection or mixture physics.

## Allocation progress — 28 October 2026

Stages 0–2 now have a standalone implementation under `biquintic/`: strict
native-only source preprocessing, local Hermite real/AD kernel, versioned
experiment cache and physical-coordinate table adapter. Optimized/debug
analytic and protocol tests pass; actual source-cell Python/Fortran comparisons
agree. Matched archived-query identities, redundant source columns, dense
positivity, coarsening holdouts, radiation-reference and printed-precision
controls are recorded in the separate progress report.

The next gate is not merely obtaining small identity residuals. Generic C2
biquintics retain some large response-Jacobian jumps in actual pure-FreeEOS
states, although many original dense tails are Skye-only. A streaming installed
MESA fraction probe and matched near-face legacy control establish that
distinction without installing the candidate. Biseptic C3 controls remove
jumps but worsen some archived-field holdouts; regional C4 soft fits improve
active-tail holdouts but do not uniformly outperform biquintic in ionization.
Keep exact biquintic, C3 and C4 alternatives as separate controls until this
fidelity/smoothness tradeoff is resolved or explicitly accepted for a bounded
default-off test. Fixed-density energy self-inversion does not exercise all
response partials and is not a substitute for stellar A/B validation.

Centered derivative masks also move the first complete potential nodes to
logT=3.28, above the clean-source threshold 3.2. The eventual upstream blend
change at 3.1 therefore does not by itself provide full candidate coverage.
Any future runtime path must fail or explicitly count fallback outside its
validated support; do not repair the cold tables or change EOS selection here.
Offline composition holdouts are feasible, with explicit audited metadata
handling where redundant display endpoints are wrong, but gauges, nonconvex
weights, source-family coverage and runtime composition remain separate gates.

### Latest bounded decision — table-only work remains the scope

Existing-table validation now includes native-only shifted stencils: the
standard source floor stays3.2, and a separate explicitly requested3.1 clean
control rejects every donor below the known-bug threshold. Both recover active
queries without changing installed source tables or EOS selection. Whole-cell
numerical Bernstein/reference-split bounds establish positive Cv and P_rho in
the declared full-own-use cells, subject to their numerical guard conventions;
they do not establish source accuracy or blend/whole-domain stability.

Double-knot quintic C3 source fits have structural null modes; a declared
integrated normal-fourth roughness selector removes that ambiguity. It improves
some active-geometry tails but none of18additional H/He/dense coarse cases
improves Cv p99 over matching C4. A hard E/TS/P-constrained C3 control is
feasible at one coarse patch but does not uniformly improve the same holdouts.
Do not treat a small identity residual, exact retained node, rank pass, or
smooth derivative as a physical-accuracy acceptance gate by itself.

The user explicitly limits this phase to existing tabulations; Stage3's
direct-source acquisition and new-table options are deferred, not prerequisites
to making useful table-only progress. Independent full-resolution offnode
accuracy remains unknown and must be reported as such. Before runtime work:

1. Evaluate the smoother's analytic-reference/gauge sensitivity, declared
   selector weights, patch-support and common-target refinement; avoid tuning
   and claiming performance on the same withheld targets.
   The bounded known a+bT audit034 now confirms combined gauge/scale/spline-
   representation sensitivity of the raw total-potential roughness objective;
   exact analytic reference removal preserves invariant fields. Declare a
   canonical reference/scaling convention or validate a gauge-invariant
   objective before selecting parameters. This is not a source-table gauge fix.
2. Make an explicit bounded experiment decision: exact nodal biquintic fidelity
   versus response-Jacobian continuity, or continue potential reconciliation.
   Keep all caller outputs consistent with the chosen potential.
3. If authorized after that decision, implement only a default-off, fixed-X/Z
   component comparison with strict metadata and counted invalid-cell fallback.
   Preserve selection/blend controls and exclude unsupported compositions.
4. Test actual solver consumers of Cp/chi and Cv partials, including standard
   gravitational heating and constant-density burning. The offline eps-grav
   control shows jumps are multiplied by thermodynamic increments; a large
   map residual alone does not demonstrate a large stellar Jacobian defect.

No runtime path is accepted or enabled by these table-only controls. The
separate progress PDF and RESUME.md carry the latest quantitative findings,
test status and next-allocation instructions.
