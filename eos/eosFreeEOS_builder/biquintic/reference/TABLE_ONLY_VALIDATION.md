# Existing-tabulation validation contract

## Permitted claims and input

A withheld native row tests **reproduction of archived fields at coarser
resolution**, not independent solver/model accuracy at that node or a new
off-node point. Consistency identities, potential polynomial tests, positivity,
and agreement between interpolators are separate diagnostics, not fresh-source
accuracy. The unlinked direct-call scaffold is deferred.

Input: `mesa/data/eosFreeEOS_data/mesa-FreeEOS_02z70x.data`, SHA-256
`e6ca7ce27bee8197bb257a8c4886a06366e6951e58958fbed3735c877ee68502`.
Header: format version 51; X=.70, Y=.28, Z=.02; 306 temperatures and 527
q blocks; t=2.10..8.20, step .02; q=-10.09..5.69, step .03.
MESA documentation attributes this family to FreeEOS 2.2.1 EOS1. The header
does **not** identify release/build hashes, solver constants, gauges, or the
complete metal pattern. README records the exact current builder pattern,
FreeEOS weights and EOS1 options with hashes, not a historical build manifest.

## Printed-node precision and coordinates

Builder `src/free_eos_table.f90` prints header values with `f14.4` (line 339),
q-block labels with `f14.6` (351), and physical rows with
`(1p,99(e22.14))` (402): **15 significant digits**, not 22-digit physical
accuracy. Preserve strings and declared precision. Base-10 log rounding
propagates approximately as `ln(10)*delta_log10` into relative physical-field
error; potential cancellation can amplify it further.

q accumulates upward from its minimum (343/426); temperature accumulates
downward from its maximum in each block (346/394), then prints ascending.
For example the actual printed minimum is `2.10000000000009E+00`, not exact
nominal 2.1. Requested `log10rho=q+2*t-12` is computed before the call (363),
stored in `logRhos` (392), and printed (421): not necessarily solver-returned
density. Integer `(iq,it)` IDs and printed coordinates are distinct metadata.
Use a uniform-header grid only after coordinate audit and explicit exclusion;
do not silently overwrite inconsistent printed values. Held-out evaluations
must share the retained grid's declared coordinate convention.

### iq520 discrepancy (zero-based)

Read-only reference audit confirms all **306 rows** at q=5.510000 have
`printed_logRho-(q+2*printed_logT-12)` approximately **-1e-4**;
171 are native and 135 fallback. The delta range is
`[-0.00010000000001220144, -0.00009999999998555609]`.

| iq | printed t | printed logRho | Offset from header-uniform mapping |
| ---: | --- | --- | ---: |
| 519 | 3.20000000000009E+00 | -1.19999999999823E-01 | -3.19e-15 |
| 520 | 3.20000000000009E+00 | -9.00999999998238E-02 | -1.00000000003347e-4 |
| 520 | 8.20000000000000E+00 | 9.90990000000000E+00 | -9.99999999962142e-5 |
| 521 | 3.20000000000009E+00 | -5.99999999998246E-02 | -3.52e-15 |

The inspected loop/format does not explain the isolated offset as ordinary
accumulation/rounding: a -1e-4 q shift would appear in the six-decimal label.
Historical cause is **unknown**; no repair is made. Exclude all iq520 nodes
and interpolation/derivative stencils touching them, including expanded coarse
stencils. Outside saved q<=5 queries does not imply general stencil validity.

## Input and potential conventions

For native accepted rows (`MESA=0`, finite/positive required fields):

```text
t=log10T; T=10**t; rho=10**log10rho
Pgas=10**logPgas; e=10**logE; s=10**logS
P=Pgas + a_reconstruct*T**4/3
chiRho/chiT: unstripped total-pressure logarithmic susceptibilities
Cv/Cp: unstripped heat capacities
dE_dRho/dS_dT/dS_dRho: physical-coordinate derivatives
free_e=10**log10_free_e; mu is dimensionless
```

cgs: P `dyn cm^-2`, rho `g cm^-3`, e `erg g^-1`, s/Cv/Cp
`erg g^-1 K^-1`. Native log fields are **base 10**; physical-log derivative
coordinates below are **natural logarithms**. Exclude fallback rows/stencils:
their documented log-base concern and zero diagnostics are not valid FreeEOS
reference or gauge calibration.

The lead reconstructs total pressure with current MESA `crad`, approximately
`7.56573325027999926e-15 erg cm^-3 K^-4`, and total e/s/Cv.
Declare **current-MESA radiation reconstruction**, not verified historical
generation metadata. No solver was retrieved to establish its constant/gauges.
Do not mix Pgas with total e/s, fit away energy/entropy offsets, or claim
material subtraction is independently verified. Separate radiation-dominated
and strongly cancelling cases. A declared sensitivity bracket in `a` studies
convention sensitivity, not the true historical value.

For fixed composition, x=ln(rho), tau=ln(T), the six constraints are:

```text
f=e-T*s
f_x=P/rho                  f_tau=-T*s
f_xx=(P/rho)*(chiRho-1)     f_x_tau=(P/rho)*chiT
f_tau_tau=-T*(s+Cv)
```

For u=ln(Q)=x-2*tau+12*ln(10), F(u,tau)=f(x,tau):
`F_u=f_x`, `F_tau=2*f_x+f_tau`, `F_uu=f_xx`,
`F_u_tau=2*f_xx+f_x_tau`, and
`F_tau_tau=4*f_xx+4*f_x_tau+f_tau_tau`.
Redundant `f_x_tau=-T*rho*dS_dRho` and `Cv=T*dS_dT` remain **closure
diagnostics**, not independent truth or quantities to force into agreement.

## Holdout/coarsening recommendations (lead implementation)

1. Freeze hashes, grid convention, masks, low-T and iq520 exclusions, stencil
   width and reconstruction metadata. Start in declared clean t>=3.2 interior.
2. Use nested stride-2 and, where viable, stride-4 retained grids in q and t.
   Specify anchors/phases/edges; compare methods/resolutions on a common valid
   cell/withheld-node set and disclose sample/exclusion counts. This changes
   resolution deliberately, rather than establishing production-grid accuracy.
3. Withhold **all fields/constraints** at targets from fitting. Estimate
   missing twists/higher derivatives from retained data only. Retained twists
   precomputed using omitted fine-grid rows leak reference information; do not
   reuse them. Disclose use of identities imposed on retained rows.
4. Require every coarse interpolation/derivative footprint to avoid excluded,
   nonnative/nonfinite rows. Do not fill failing coarse stencils with fine-grid
   estimates. Apply iq520 exclusion to the widened footprints too.
5. Predict omitted archived nodes under the same coordinate/radiation/gauge
   convention. Compare P/e/s/Cv/Cp/chi/available derivatives separately from
   identities, polynomial tests, positivity and column closure. Report scaled
   absolute/relative errors, tails/maxima, negative heat capacities and regime
   counts; retain failures, do not clip or floor them away.
6. Label results **coarse-grid archived-field holdout errors**. Report common-
   sample resolution trends and ionization/radiation subsets. Tiny identities
   or nodal reproduction alone cannot set physical-accuracy tolerances.

No holdout script is implemented here, no tables/caches are regenerated, and
no fresh solver, official release hash or convergence/gauge audit is available.
