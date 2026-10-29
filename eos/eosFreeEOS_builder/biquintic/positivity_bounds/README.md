# Experiment024: numerical full-cell positivity bounds

Exclusive sidecar for the existing **centered radius2 mean TOTAL-potential**
biquintic. No core, runtime, cache, debugger, approval configuration, report
source, or installed EOS is edited. Lead report integration uses the generated
`freeeos_biquintic_progress/experiments/024_positivity_bounds/SUMMARY.md` and
`results.json`; the main report remains lead-owned.

## Reproduce

From the workspace root, with the existing Python/NumPy environment:

```bash
OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 -B mesa/eos/eosFreeEOS_builder/biquintic/positivity_bounds/test_bounds.py -v
OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 -B mesa/eos/eosFreeEOS_builder/biquintic/positivity_bounds/run_bounds.py
OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 -B mesa/eos/eosFreeEOS_builder/biquintic/positivity_bounds/verify_artifacts.py
```

The source is the small `mesa-FreeEOS_02z70x.data` table, SHA256
`e6ca7ce27bee8197bb257a8c4886a06366e6951e58958fbed3735c877ee68502`.
Only complete native cells from the unmodified core's radius2-mean construction
are used. Source constraints and twists are existing float64 values, treated
as fixed polynomial inputs, not recomputed with a new estimator.

## Polynomial and bounds

Let `u=ln(Q)`, `tau=ln(T)`, `x=(u-u0)/du`, `y=(tau-tau0)/dtau`.
All nine `f_(u^i tau^j)`, `i,j=0,1,2`, at each of four corners multiply the
core's quintic Hermite basis with native spacing powers. Their tensor product
forms the full **6x6 power coefficient** array in normalized `x,y`.
The core uses one global dq/dt in its H basis. Each sign polynomial is then
rescaled by actual stored endpoint width/global spacing, so its Bernstein
`[0,1]^2` rectangle covers the entire closed cell even when binary grid widths
differ very slightly. This preserves the core polynomial, not a new local H fit.
Physical density is `ln(rho)=u+2*tau-12*ln(10)`, so fixed-density derivatives
are `d/dln(T)=d/dtau-2*d/du`. In natural-log units:

- `T Cv = f_tau - 2 f_u - f_tautau + 4 f_utau - 4 f_uu`.
- `P_rho = f_u + f_uu` (total-pressure derivative at fixed temperature).

Differentiate power coefficients with the native `du,dtau` factors, then
convert/elevate both sign polynomials to degree `(5,5)` Bernstein coefficients.
For one axis, `b_i = sum_(k<=i) a_k * C(i,k)/C(5,k)`.
The exact-polynomial convex-hull property bounds its values by the minimum
and maximum Bernstein coefficients over the **entire closed rectangle**.
Strict `minimum > numerical_margin` certifies numerical positivity; samples
do not create certificates. Zero/nonpositive minima are ambiguous, not failures.

Each ambiguous polynomial is split dyadically into four subrectangles with
de Casteljau, depth-first. All leaves must pass for a parent positivity result.
Default depth is 4, with at most 341 visited nodes per cell/field, 13 live stack
nodes, and batches of 2048 cells. A cell cap of 150000 refuses larger inputs
rather than silently truncating. `--max-depth 3 --max-nodes 85` is a smaller
bound. Output overrides must stay under the experiment024 directory.

Negative witness candidates are evaluated by the **unmodified core at actual
logQ/logT coordinates**, forced into the original cell. Only a finite negative
returned Cv or P_rho is recorded as a failure. A negative coefficient alone
never fails a cell. Reaching the depth/node budget leaves an explicit unresolved
status and unresolved covering-leaf count.

## Floating-point scope

Assembly, derivatives, conversion, bounds and subdivision use `np.longdouble`.
The root guard is `4096*(eps_longdouble+eps_float64)*M`, where `M` is the
maximum Bernstein coefficient obtained from **absolute Hermite assembly and
derivative contributions**, before cancellation. Its longdouble component
guards sidecar calculation; its float64 component conservatively budgets core
evaluation disagreement. Children inherit the guard and add
`4096*eps_longdouble*max(abs(parent_coefficients))` per subdivision level.

These are conservative **numerical guards, not watertight IEEE certificates**:
there is no directed rounding, interval enclosure, or proven operation-count
error theorem. No uncertainty in source jets, source accuracy, global stability,
blended/runtime EOS, or independent radiation convention is certified.
JSON explicitly records `floating_proof_claim: false`, machine precision,
margin policy, source/core/sidecar checksums and calibration results.
The analytic radiation split is excluded because its exponential term is not
a biquintic polynomial. No Pgas, E or S full-cell certificates are claimed.

## Domain and calibration

The primary subset requires an **entire cell in one own-use branch**:
`tmin>=3.2` and `(tmin>=6.5 OR qmax<=-0.5 OR maxlogrho<=-6.5)`.
No corner-wise union or boundary tolerance substitutes for that condition.
Material screening requires all four **source** corners `Prad/Pgas<1e-3`,
the full q interval in `[-8,5]`, and the full t interval in `[3.2,7.8]`.
The radiation fraction screen is not a continuous material-fraction bound.
Fully valid cells outside this full-one-branch subset, including partial-domain
and union-only corner cells, are counted separately and never included in the
primary claim.

Self-contained tests cover analytic natural-coordinate signs, all nine jets,
full degree5x5 Hermite reproduction, equality with the core H basis,
degree-five Bernstein elevation, subdivision coordinate mapping, adaptive
positivity, zero ambiguity, witnessed negatives, bounded unresolved covers,
full-cell domain cuts, and synthetic sample agreement against the core.
The data run reproduces all nine jets on four corners of 128 cells spread across
the complete table domain and compares
the two sign polynomials against four deterministic **actual-coordinate** core
queries in every complete cell. Errors outside the chosen guards invalidate
the calibration. This is polynomial implementation calibration, not source
truth, and the finite queries are not the full-cell certificate mechanism.

## Products and lead integration

- `results.json`: complete subset/field counts, root versus adaptive passes,
  worst coordinates/bounds/margins, budgets, calibration and SHA256 provenance.
- `cell_bounds.npz`: original cell indices, frozen full-one-branch/material flags,
  every bound/margin, status, recursion depth/node/leaf counts; no raw query dump.
- `negative_witnesses.json`: all core-evaluated negative coordinates and values.
- `unresolved_cells.json`: every depth/budget ambiguity, never reclassified positive.
- `SUMMARY.md`: lead-ready counts with numerical-only and domain qualifiers.

The verifier recounts every archived status/subset, checks all recorded SHA256
values, reproduces domain masks and 64 spread-out root bounds, and independently
re-evaluates every negative witness against the core at its saved coordinates.

Lead wording: “Numerical Bernstein bounds with conservative roundoff guards
certify positivity of the implemented total-potential polynomial on the reported
material-screened/full-one-branch cells.” Keep outside-domain diagnostics
separate and retain unresolved counts. Do not replace “numerical” with “rigorous,”
infer source stability, or claim a cellwide material fraction from four corners.
