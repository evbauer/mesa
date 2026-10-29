# Variational normal-jet C3 experiment

Bounded table-only experimental high-jet selection, **not selected runtime**.
No core/C3-worker/runtime/default/cache changes, solver download/build, new EOS
tables, or fresh-source accuracy claim. Binary outputs contain audit masks and
paired queries/values, not interpolation coefficient caches.

## Reproduce

From the workspace root, using existing Python/NumPy/SciPy:

```sh
python3 -B mesa/eos/eosFreeEOS_builder/biquintic/variational_c3/test_variational.py
python3 -B mesa/eos/eosFreeEOS_builder/biquintic/variational_c3/run_control.py --overwrite
```

Only this source directory and
`freeeos_biquintic_progress/experiments/025_variational_c3/` receive writes.
The latter contains `README.md`, `results.json`, `audit_arrays.npz` and test/
validation logs. `-B` prevents bytecode writes in imported worker directories.
The runner pins the existing solar table SHA-256 to
`e6ca7ce27bee8197bb257a8c4886a06366e6951e58958fbed3735c877ee68502`.
It records imported/local source hashes, experiment-015 region hash, versions,
arguments, elapsed time and peak Linux process RSS. The full coarse-donor run
takes about one second here and peaks near 147 MiB; evaluations are limited to
1,000 queries by default, 2,000 maximum. No large saved offgrid archive is loaded.

## Objective and solve

Native interpolation coordinates are `u=ln(10)*logQ`, `tau=ln(T)`;
`ln(rho)=u+2*tau-12*ln(10)`. Preserve all six supplied native low jets, obtained
with the existing six-constraint transform. Replace only pure `(3,0)` and
`(0,3)` centered-FD estimates. On each contiguous valid retained line, use the
known endpoint `f`, normal `f'`, and normal `f''`, and find shared nodal normal
`f'''` minimizing the sum of `integral (f'''')^2 d(normal coordinate)` over its
septic segments. This objective is algebraic smoothness, not physical EOS energy.

For uniform spacing `h`, write endpoint data as
`v=(f_L,h*f'_L,h^2*f''_L,h^3*f'''_L,f_R,...,h^3*f'''_R)`.
The existing exact-rational cardinal septic polynomials give
`K_ij=integral H_i'''' H_j'''' dz`; segment energy is `h^-7*v^T*K*v`.
For the scaled unknowns `g=h^3*f'''`, their segment block is exactly

```text
K[(3,7),(3,7)] = [[16,4],[4,16]]
```

Assembling shared nodal unknowns gives a symmetric positive-definite
tridiagonal matrix: diagonal 16 at line ends, 32 at interior nodes, offdiagonal
4. Solve directly with SciPy `solve_banded`; no iterative/nonlinear solver.
Known-data contributions use a segmentwise subtracted quadratic Taylor
polynomial to reduce gauge/large-value cancellation. Its fourth derivative is
identically zero; this conditioning transformation changes neither the objective
nor any preserved source constraint/gauge. Four-point Gauss integration exactly
integrates the squared degree-three fourth-derivative polynomial and keeps the
reported energies nonnegative.

Free third jets at segment endpoints give natural zero normal fourth derivative;
interior fourth derivatives agree in that **one-dimensional normal solve**.
This is not full tensor C4 continuity. Quadratic/cubic polynomials reproduce
their exact third jets with zero minimum energy. Arbitrary degree-seven
reproduction is **not** expected with these natural boundary conditions.

All other estimated jets use the unchanged C3 route graph, strict centered
radius-2 donor masks, and both-route means. The `(2,1)`/`(1,2)` routes start
from supplied seconds; subsequent mixed jets incorporate the new normal thirds
where the existing graph does. All 16 shared nodal jets must be complete and
all four cell corners complete. No one-sided FD, zero substitute or invalid
bridge is introduced. Existing C3 evaluation/cardinal bases are imported read-only.

## Dependencies and masks

Unlike the FD baseline, each variational normal jet has a nonlocal dependency
on its **entire contiguous coarse native line segment**. Subsequent mixed
routes propagate those segment dependencies through their centered donor windows.
There is no finite global donor-radius claim; line end indices/lengths and all
per-jet masks are saved. Solves do not cross invalid, cold or original iq520
nodes. Source-constraint preservation is checked exactly on complete nodes.

Stride 2 uses original even/even rows only. Native masks also inspect original
flags, finite columns and printed-coordinate agreement, without fitting withheld
constraints. A separate conservative bridge guard removes both retained
endpoints of any accepted coarse edge whose intervening original row is invalid.
The pinned table has 148 density-direction and 8 temperature-direction such
edges; 275 otherwise valid retained nodes are excluded from **every** model.
This prevents long-line and mixed-FD bridging even over discarded odd rows.
It does not change the common patch target set. This stricter masked FD baseline
is an isolated sidecar control, not a rewrite of experiment 011/C3-worker results.

Global valid coarse donors are used, not solves truncated to the four patches.
The guarded run solves 348 density-normal segments (maximum 260 nodes) and
349 temperature-normal segments (maximum 125). Maximum normalized linear-system
stationarity residual is below 2.7e-16. All segment energy ratios are <=1 relative
to a feasible FD-interior competitor that uses minimizing thirds where FD lacks
endpoint support. This comparator is explicitly **not** the complete FD model
on an unsupported boundary. Median energy ratios are about 0.99988/0.99975;
small objective reduction does not establish physical improvement or integrability
of the printed constraints.

## Paired domains and conventions

Compare radius-2 mean-route BQ, radius-2 FD C3, and variational C3, all using the
same coherent **current-CRAD analytic-reference split/restoration**. Do not fit a
radiation constant or entropy/energy gauge; the remainder is not verified pure
material. Decode base-10 printed Pgas/E/S and total Cv/Cp/chi fields. Reconstructed
total P uses `Pgas+current_crad*T^4/3`; it is not an independent source field.

Use the exact experiment-015 rectangles: ionization, helium_ionization,
dense_tail and active_chiT. The latter is a geometry label, not measured MESA
component fractions. Require clean native targets, material `Prad/Pgas<1e-3`,
declared own-use, q in [-8,5], logT in [3.2,7.8], and material/own-use cell
corners. Intersect all three complete-cell masks **before comparison**.
All 746 eligible odd/odd patch targets are common and used: 180/192/170/204
by region. None is retained by stride 2. All constraints and estimated jets are
constructed from retained rows only, not fine-grid estimates.

One paired deterministic offnode query at coarse-cell coordinates (0.37,0.63)
is evaluated in each selected holdout cell. Raw values and physical-log-coordinate
Jacobians are archived. Departures from BQ are **model differences**, not offnode
accuracy. Values divide by abs(BQ value); Jacobian departures divide by
max(1,abs(BQ field),abs(BQ partial)). Physical `_x` partials hold T fixed; `_t`
partials hold rho fixed and therefore shift `u` by `-2*delta tau`.

## Findings and limits

On these same 746 printed holdouts, Cv p99/max relative errors are:

| Model | Cv p99 | Cv max |
| --- | ---: | ---: |
| BQ reference split | 2.872e-3 | 8.845e-3 |
| FD C3 | 6.550e-3 | 2.054e-2 |
| Variational C3 | 5.939e-3 | 1.846e-2 |

Variational Cv improves at 85.0% of targets versus FD C3, but aggregate Cv
remains worse than BQ. Improvement is not universal across fields: Pgas p99
increases from 3.013e-5 to 3.424e-5, while E/S p99 roughly halve.
Ionization Cv p99 improves strongly (3.320e-3 to 8.711e-4); dense-tail Cv remains
1.191e-2 versus BQ 6.165e-3. Full field/max/region results and paired raw errors
are in JSON/NPZ; the output README includes the full field summary.

Offnode Cv p99 departure from BQ decreases from 2.344e-3 to 1.976e-3, but scaled
Cv_t p99 departure increases from 0.545 to 0.570. Those changes are not evidence
of derivative accuracy. All sampled heldout/offnode/face physical fields and
positive responses are finite/positive; this is not a global positivity bound.

At 64 faces per native axis, both C3 controls give zero all-16-jet native trace
jumps and zero physical response-Jacobian jumps using the shared cardinal
construction. BQ response Jacobians can jump, as expected for C2 potential.
Central-difference checks at 32 common offnodes have maximum scaled Jacobian
disagreement below 8e-6 for these controls, with physical-log step 1e-5.
These are consistency/construction checks, **not independent fresh-source
accuracy**, and no formal convergence order is inferred.

Six tests cover known quadratic/cubic jets, SPD/minimum energy and stationarity,
nonlocal influence with invalid-line isolation, original midpoint bridge guards,
six-constraint preservation/shared C3 faces, and analytic physical Jacobians of
the returned smooth potential. The normal-jet selection is feasible and preserves
C3 smoothness, but these results do **not** justify replacing primary BQ or
selecting this nonlocal experimental construction for runtime use.
