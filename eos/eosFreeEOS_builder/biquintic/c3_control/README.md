# Offline C3 biseptic Hermite feasibility control

Sidecar only: no runtime integration, primary-candidate replacement, source-core
edits, EOS/cache exports, tables, or solver builds. This tests the biquintic
roadmap's smoothness-versus-estimated-derivatives tradeoff; it does not compete
with the separately owned constrained C4 B-spline patch or table adapter.

## Construction

Import `freeeos_potential.source_constraints` unchanged. Coordinates are
`u=ln Q`, `tau=ln T`, with `x=ln rho=u+2 tau-12 ln 10` and physical
`D_t=D_tau-2 D_u`. Six supplied jets are `00,10,01,20,11,02`.
All `J_ab=D_u^a D_tau^b f`, `0<=a,b<=3`, are stored: **16 nodal jets**;
the default evaluator returns the **10 partials of total order <=3**.
`30` and `03` are estimated, not supplied; jets of total orders 4–6 are
also estimates, not additional FreeEOS observations.

In increasing total order, use centered native-grid first-difference operators
`delta_u,delta_tau` of radius 2 or 3 (formal orders 4 or 6):

```text
J30 = delta_u J20                  J03 = delta_tau J02
J21 = (delta_u J11 + delta_tau J20)/2
J12 = (delta_u J02 + delta_tau J11)/2
J31 = (delta_u J21 + delta_tau J30)/2
J22 = (delta_u J12 + delta_tau J21)/2
J13 = (delta_u J03 + delta_tau J12)/2
J32 = (delta_u J22 + delta_tau J31)/2
J23 = (delta_u J13 + delta_tau J22)/2
J33 = (delta_u J23 + delta_tau J32)/2
```

Redundant routes are explicit, equally averaged, and their disagreements are
reported. Each route requires **every** donor, including its zero-weight center,
valid in the parent jet mask; each average requires **both** routes. Final
complete nodes require all 16 masks, cells require all four complete corners.
No one-sided, fallback, cold, iq520, or zero-twist substitution. Conservative
source support has total index reach `<=4r`, each-axis reach `<=3r`; actual
per-jet masks are archived. Donors obey the source native/finite/positive/
density-consistent/logT>=3.2 rules. Own-use/material cuts are query cuts, not
donor cuts: accepted native donors outside those cuts remain allowed.

On `z=(u-u_left)/h_u`, the eight 1D cardinal polynomials solve exactly
`M_(endpoint,k),p = p!/(p-k)! endpoint^(p-k)`, `p>=k`, otherwise zero,
for endpoints 0,1 and derivative orders `k=0..3`. `M inverse` is obtained
by exact rational elimination, not a numerical inverse. Equivalently,

```text
H_left,k(z) = z^k/k! (1-z)^4 sum[j=0..3-k] binomial(3+j,j) z^j
H_right,k(z) = (-1)^k H_left,k(1-z)
f_cell = sum[i,j=left,right; a,b=0..3]
         J_ab(i,j) h_u^a h_tau^b H_i,a(z_u) H_j,b(z_tau)
```

Derivative bases carry `h^(k-order)`. Exact endpoint cardinality is used at
exact faces; tests also reproduce arbitrary degree-7-by-7 polynomials off-node
through all 16 partials. Shared jets give C3 potential continuity, hence
continuous response Jacobians where thermodynamic denominators are nonsingular.
No C4, monotonicity, convexity, or global positivity guarantee follows.

## Reproduce

From the workspace root, with existing Python and NumPy:

```bash
python3 -m unittest discover -s mesa/eos/eosFreeEOS_builder/biquintic/c3_control -p 'test_*.py' -v
python3 mesa/eos/eosFreeEOS_builder/biquintic/c3_control/run_control.py
```

Defaults read the existing solar printed table and
`freeeos_report/plots/offgrid_sample_data.npz`. Outputs live only in
`freeeos_biquintic_progress/experiments/011_c3_control/`:
`results.json`, `SUMMARY.md`, `tests.txt`, `domains.npz`, `face_queries.npz`.
Fingerprints, original flattened query IDs, coordinates, source/jet/coarse masks,
and forced left/right face cells are saved. Binary products are audit masks and
queries, **not new EOS tables**. Evaluation is chunked; no solver is built.

Both radii, total potential, and declared analytic-radiation-reference-split
controls are paired with the unchanged existing radius-matched mean-route
biquintic on the **same complete intersection of all eight models**. Saved
eligibility and own-use cuts are preserved; the material subset additionally
requires saved material and material source corners. P/Pgas/E/S/Cv/Cp/chi/Gamma/
grad_ad values and physical-log-coordinate response partials are compared.
Value departures divide by `abs(biquintic value)`; partial departures divide
by `max(1,abs(biquintic response),abs(biquintic partial))`. Faces use symmetric
scales including both traces. Nonfinite and positivity counts remain explicit.

Stride-2 tests construct the **entire 16-jet estimator using only even-even
coarse source rows**, then evaluate odd-odd withheld printed nodes in the common
complete domain, with identical queries for biquintic. Printed fields cannot
validate the ten estimated high jets directly; this measures the estimator's
effect on returned printed-field fidelity. It is not independent source
accuracy at full-resolution off-grid queries. Total/split modes preserve the
existing radiation convention; the split remainder is not verified pure material.

## Results And Limits

The saved-query comparison retains **964,659 own-use** queries, including
**293,331 material/own-use** queries. Radius-2 total C3 material p99 departures
from biquintic are `P=4.78e-7`, `Cv=1.83e-4`, and scaled `Cv_t=4.76e-2`.
Radius-3 gives `2.45e-7`, `1.56e-4`, and `4.75e-2`, respectively.
All sampled exact-face C3 value and response-Jacobian jumps are zero in this
implementation; material biquintic response-Jacobian maxima reach about 0.40.
That is a construction check, not proof of accuracy or superiority.

On **5,437 common material stride-2 withheld nodes**, C3 Cv p99 errors are
`1.35e-3`/`1.37e-3` for radii 2/3, versus biquintic
`6.99e-4`/`7.04e-4`. C3 maxima are about **2.05%–2.06%**, versus
biquintic **0.884%–0.886%**. No tested model has sampled nonpositive/nonfinite
P/Pgas/Cv/P_rho/chiRho/E/S on either saved or heldout subset. See experiment
`SUMMARY.md` and `results.json` for broad own-use tails, all fields, positivity,
per-jet route disagreement, domain loss, worst coordinates, and convergence.

The analytic exponential full-estimator maximum nodal-jet errors fall from
`5.73e-5` to `3.58e-6` (radius 2) and `2.93e-7` to `4.58e-9` (radius 3)
when natural spacing halves from 0.125 to 0.0625. This verifies controlled
smooth-input convergence, not printed-data integrability. More estimated
derivatives reduce donor coverage and can amplify noise/incompatibility.
**C3 smoothness is feasible, but this heldout result does not support replacing
the primary biquintic or claiming superior fidelity.**
