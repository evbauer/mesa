# Table-implied radiation audit

## Reproduce and scope

```sh
python3 -B mesa/eos/eosFreeEOS_builder/biquintic/reference/test_infer_radiation.py
python3 -B mesa/eos/eosFreeEOS_builder/biquintic/reference/infer_radiation.py --overwrite
```

The script pins the table SHA-256 to
`e6ca7ce27bee8197bb257a8c4886a06366e6951e58958fbed3735c877ee68502`
and refuses any other input. It writes only
`freeeos_biquintic_progress/experiments/005_radiation_audit/`:
`results.json` (distributions/selection/rounding/rho-T trends), `nodes.csv`
(selected node IDs and diagnostics), `trends.csv`, and `SUMMARY.md`.
Existing products require explicit `--overwrite`. No solver, new tables,
interpolation cache, constant fit, or material-physics fit is involved.

## Identities and selection

Using printed physical density/temperature and independent existing derivative
columns, with **total** chiT:

```text
P_Maxwell = -rho**2*T*dS_dRho/chiT
P_energy  = rho**2*dE_dRho/(1-chiT)
a_route   = 3*(P_route-Pgas)/T**4
Cp_check  = Cv + P_route*chiT**2/(rho*T*chiRho)
```

The two pressure routes share rho, chiT, and the originating table. Agreement
is table closure, not independent source accuracy. Their `a` estimates identify
the pressure-strip reconstruction, not automatically the energy radiation EOS.

Default selection: native flags; logT>=3.2; positive chiRho/chiT/Cp/Cv;
all finite printed fields; coordinate agreement within 1e-10 dex; iq520
explicitly excluded; complete 5x5 clean native neighborhood; current-CRAD
Prad/Pgas>=100; |1-chiT|>=0.25; positive finite identity-implied pressure.
Current CRAD is `7.565733250280004e-15` in cgs, calculated with the same
CODATA-2018 expression as the lead (algebraically equivalent to current MESA).
Selection is conditioned on that declared constant, not a fitted one.

## Results (2026-10-28)

Of 161,262 rows, 122,496 meet native/coordinate/finite checks; 116,047 have
clean neighborhoods; **25,195** meet the conditioned radiation selection.
Current Prad/Pgas ranges 100.06..225,168.87; |1-chiT| ranges
2.97064..3.00002, so the energy denominator is well separated from zero.

| Prad/Pgas cut | Nodes | Maxwell median a/CRAD-1 | Energy median | Route difference abs p99 | Cp current abs p99 |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 100 | 25,195 | -7.33e-15 | -8.10e-15 | 1.50e-14 | 1.26e-13 |
| 10,000 | 10,034 | -7.88e-15 | -9.21e-15 | 1.54e-14 | 9.04e-14 |
| 100,000 | 4,181 | -6.44e-15 | -6.66e-15 | 1.56e-14 | 5.09e-14 |
| 200,000 | 1,323 | -4.77e-15 | -4.22e-15 | 1.63e-14 | 5.11e-14 |

Across the main sample, relative deviations from current CRAD have abs p99
2.33e-13 (Maxwell) and 2.40e-13 (energy), maxima 2.68e-13 and 2.80e-13.
Thus **current MESA CRAD is sufficient for the table-implied total-pressure
reconstruction and Cp closure at approximately printed precision**.
Do not replace this pressure-strip constant with the energy proxy below.

Each printed scientific token supplies a half-last-decimal-unit estimate.
The script propagates rho/T/logPgas/derivative/chiT rounding through both
routes, plus a 20-epsilon floating-point allowance. Relative bounds have
median 6.58e-14 / 7.77e-14 and p99 2.76e-13 / 2.88e-13.
All route differences fall within the combined estimate (maximum ratio 0.131).
Current CRAD exceeds the per-node estimate at 45/62 nodes, by at most
1.19/1.21 times the estimated bound. These are printing/arithmetic estimates,
**not guaranteed historical input bounds, solver error bars, or confidence
intervals**. Small rho/T trends remain at approximately 1e-13; the saved bins
make them explicit rather than interpreting them as a changed physical constant.

## Pressure stripping versus energy radiation-like remainder

No other physics is fitted away. Two supplementary proxies are:

```text
a_erho_proxy = -rho**2*dE_dRho/T**4
a_linear_density_pressure_proxy = 3*P_current*(1-chiRho)/T**4
```

The first includes any material density dependence of energy. The second
assumes material pressure is locally linear in density; it is not a universal
identity. Neither should be substituted unconditionally for a source constant.

In a fixed diagnostic subset, logT>=7.5 and Prad/Pgas>=200,000 (**1,043 nodes**):

- Energy-density proxy median/CRAD-1: **+4.4068383811e-6**; range
  +4.4068382847e-6..+4.4068384746e-6, dimensional median ~7.56576659e-15.
- Linear-density pressure proxy median/CRAD-1: **+4.4068383782e-6**.
  Their difference abs p99 is 7.08e-14 CRAD; energy-proxy rounding p99 is
  7.83e-14 CRAD. The common ~4.4 ppm excess is far above printing estimates.
- The equivalent excess pressure is 88.3%..99.2% of **printed Pgas**, median
  95.5%. This explains why strong radiation dominance need not mean the printed
  gas field is a pure material-pressure reference.

The erho proxy flattens with increasing temperature/dominance; lower-temperature
or less-dominant nodes show material contamination. Entropy-density, Cv,
absolute E, and absolute S proxies are also retained; they contain material
terms, and E/S additionally carry reference gauges. They are not forced into
agreement. The table supports a **distinct radiation-like energy/remainder
contribution**, consistent with an internal radiation coefficient differing
from the pressure-strip coefficient. It does **not verify** the upstream
constant, distinguish all possible non-radiation physics, or establish gauges.
It is not correct to identify the shared pressure-route a with the energy
radiation coefficient simply because both pressure identities agree.

## Conditioning and implications for the lead

Inferring a is well conditioned here: subtraction-condition maxima ~1.02
(~1.00001 in the tight subset), because radiation dominates Pgas.
By contrast, total/Pgas reaches 225,170 and 1/chiRho reaches **2.92e7**:
recovering a small material remainder or susceptibility can strongly amplify
small total-field/convention errors. Good total-pressure/Cp identity closure
does not validate coarse interpolation of these remainders.

Keep current CRAD for pressure reconstruction. A coherent whole-potential
current-CRAD subtraction/restoration remains an algebraic reference split,
but its remainder includes the table-implied radiation-like excess; do not
label it verified pure material thermodynamics. A different analytic reference
would require coherent potential/derivative treatment while preserving the
declared printed-Pgas stripping convention. Do not tune a to minimize coarsening
errors or silently remove the observed residual as an alleged source correction.

Five small analytic tests pass, including unequal internal-radiation/strip
coefficients, pressure/Cp identities, print rounding, and iq520 stencil isolation.
These validate audit algebra only. No fresh-source accuracy is claimed.
