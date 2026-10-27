# Mixture-level WDM_AV validation vs Skye / HELM

Relative error of the additive-volume WDM_AV mixture (integrated eosDT path, `eosDT_get_component(i_eos_WDM_AV)`) against Skye / HELM run directly on the same five reference mixtures, at the reference grid nodes where both sides are finite. `S` is the F1-corrected entropy for Skye (reference minus `Delta s_ideal`); for HELM, radiation is stripped from WDM_AV's e/s/Cv and gas-only chiRho/chiT come from WDM_AV's dlnPgas columns. HELM grad_ad/gamma1 mix a radiation-definition difference and are taken only where Prad/Pgas < 1e-2. Full per-mixture, per-region data (p90, p99, max, all quantities) is in the CSV/JSON.

## Aggregate median relative error (median over the 5 mixtures)

### (a) Skye

| region | P | E | S | chiT | Cv | grad_ad |
|---|---|---|---|---|---|---|
| whole_valid | 2.5e-05 | 1.0e-03 | 1.3e-04 | 1.4e-04 | 1.6e-04 | 6.5e-05 |
| skye_valid | 1.6e-05 | 4.3e-04 | 1.7e-04 | 1.4e-04 | 1.5e-04 | 5.9e-05 |
| skye_valid_no_pair | 2.0e-05 | 5.5e-04 | 1.9e-04 | 1.5e-04 | 1.7e-04 | 7.6e-05 |
| target | 2.8e-04 | 4.8e-02 | 2.5e-04 | 9.2e-04 | 5.0e-04 | 7.1e-04 |

### (b) HELM

| region | P | E | S | chiT | Cv | grad_ad |
|---|---|---|---|---|---|---|
| whole_valid | 2.2e-04 | 3.1e-04 | 7.4e-03 | 2.1e-04 | 6.3e-05 | 7.7e-04 |
| skye_valid | 1.9e-04 | 2.6e-04 | 1.1e-02 | 1.7e-04 | 5.8e-05 | 5.5e-04 |
| skye_valid_no_pair | 2.2e-04 | 2.8e-04 | 1.3e-02 | 1.8e-04 | 6.1e-05 | 5.5e-04 |
| target | 1.9e-03 | 4.2e-03 | 1.2e-02 | 1.6e-03 | 3.5e-04 | 1.1e-03 |

## Per-mixture, target region (logRho in [-2,4], logT in [3,7])

| set | mixture | stat | P | E | S | chiT | Cv | grad_ad |
|---|---|---|---|---|---|---|---|---|
| skye | mix_c12_o16 | median | 2.8e-04 | 4.3e-02 | 2.5e-04 | 9.2e-04 | 4.5e-04 | 6.1e-04 |
| skye | mix_c12_o16 | p99 | 2.2e-02 | 5.2e-02 | 4.9e-02 | 5.0e-01 | 4.3e-01 | 4.7e-01 |
| skye | mix_c12_o16 | max | 5.0e-01 | 5.9e-02 | 1.5e-01 | 2.9e+00 | 2.8e+00 | 6.5e+00 |
| skye | mix_o16_ne20_mg24 | median | 2.8e-04 | 4.8e-02 | 2.3e-04 | 1.2e-03 | 5.0e-04 | 7.1e-04 |
| skye | mix_o16_ne20_mg24 | p99 | 1.5e-02 | 5.3e-02 | 5.0e-02 | 1.6e-01 | 3.5e-01 | 2.0e-01 |
| skye | mix_o16_ne20_mg24 | max | 1.3e-01 | 5.7e-02 | 9.3e-02 | 4.0e+00 | 1.8e+00 | 1.9e+00 |
| skye | mix_he4_c12 | median | 7.7e-04 | 3.9e-01 | 2.2e-03 | 5.1e-03 | 7.0e-03 | 5.7e-03 |
| skye | mix_he4_c12 | p99 | 1.1e-01 | 8.1e-01 | 2.0e-01 | 8.1e+00 | 4.2e-01 | 8.0e+00 |
| skye | mix_he4_c12 | max | 2.6e-01 | 8.3e-01 | 2.5e-01 | 2.2e+01 | 8.2e-01 | 2.2e+01 |
| skye | mix_h1_he4 | median | 1.1e-04 | 2.6e-02 | 3.3e-04 | 9.1e-04 | 1.1e-03 | 8.7e-04 |
| skye | mix_h1_he4 | p99 | 2.8e-02 | 1.1e+00 | 4.0e-02 | 4.3e+00 | 2.9e-01 | 4.1e+00 |
| skye | mix_h1_he4 | max | 8.8e-02 | 1.5e+00 | 2.5e-01 | 8.8e+00 | 2.4e+00 | 8.8e+00 |
| skye | mix_c12_o16_ne22 | median | 2.6e-04 | 4.9e-02 | 2.2e-04 | 8.4e-04 | 4.6e-04 | 5.5e-04 |
| skye | mix_c12_o16_ne22 | p99 | 1.6e-02 | 6.0e-02 | 5.0e-02 | 3.0e-01 | 3.1e-01 | 2.5e-01 |
| skye | mix_c12_o16_ne22 | max | 2.2e-01 | 6.6e-02 | 1.1e-01 | 2.1e+00 | 1.4e+00 | 1.6e+00 |
| helm | mix_c12_o16 | median | 1.9e-03 | 4.3e-03 | 1.3e-02 | 1.7e-03 | 3.5e-04 | 1.1e-03 |
| helm | mix_c12_o16 | p99 | 6.1e-03 | 2.7e-02 | 1.5e-01 | 8.1e-03 | 5.1e-03 | 2.3e-02 |
| helm | mix_c12_o16 | max | 6.8e-03 | 3.9e-02 | 1.9e-01 | 9.0e-03 | 5.9e-03 | 2.3e-02 |
| helm | mix_o16_ne20_mg24 | median | 1.8e-03 | 4.2e-03 | 1.1e-02 | 1.6e-03 | 2.8e-04 | 1.4e-03 |
| helm | mix_o16_ne20_mg24 | p99 | 5.0e-03 | 2.0e-02 | 1.4e-01 | 5.7e-03 | 3.1e-03 | 2.4e-02 |
| helm | mix_o16_ne20_mg24 | max | 5.3e-03 | 2.3e-02 | 2.0e-01 | 7.1e-03 | 4.5e-03 | 2.4e-02 |
| helm | mix_he4_c12 | median | 1.5e-02 | 3.3e-02 | 1.2e-02 | 1.2e-02 | 2.0e-03 | 3.6e-03 |
| helm | mix_he4_c12 | p99 | 5.2e-02 | 1.4e-01 | 1.2e-01 | 5.9e-02 | 2.4e-02 | 3.5e-02 |
| helm | mix_he4_c12 | max | 6.5e-02 | 1.7e-01 | 1.4e-01 | 8.4e-01 | 9.0e-02 | 5.3e+00 |
| helm | mix_h1_he4 | median | 1.9e-03 | 4.1e-03 | 8.7e-03 | 1.5e-03 | 6.9e-04 | 7.0e-04 |
| helm | mix_h1_he4 | p99 | 1.2e-02 | 3.4e-02 | 6.6e-02 | 1.4e-02 | 4.5e-03 | 2.6e-02 |
| helm | mix_h1_he4 | max | 3.6e-02 | 4.1e-02 | 7.7e-02 | 8.2e-02 | 1.4e-02 | 1.5e-01 |
| helm | mix_c12_o16_ne22 | median | 1.9e-03 | 4.1e-03 | 1.3e-02 | 1.6e-03 | 2.9e-04 | 1.1e-03 |
| helm | mix_c12_o16_ne22 | p99 | 5.2e-03 | 1.9e-02 | 1.3e-01 | 5.9e-03 | 2.8e-03 | 2.3e-02 |
| helm | mix_c12_o16_ne22 | max | 5.6e-03 | 2.3e-02 | 1.8e-01 | 6.5e-03 | 3.3e-03 | 2.3e-02 |

## Per-mixture, Skye-valid region (logT >= 6.5 OR logRho >= 3.7)

| set | mixture | stat | P | E | S | chiT | Cv | grad_ad |
|---|---|---|---|---|---|---|---|---|
| skye | mix_c12_o16 | median | 1.5e-05 | 3.5e-04 | 3.6e-05 | 1.2e-04 | 1.4e-04 | 5.2e-05 |
| skye | mix_c12_o16 | p99 | 2.8e-04 | 4.2e-02 | 7.1e-02 | 1.7e+00 | 4.4e-01 | 1.7e+00 |
| skye | mix_c12_o16 | max | 1.6e-03 | 4.5e-02 | 3.2e-01 | 3.3e+00 | 8.9e-01 | 3.3e+00 |
| skye | mix_o16_ne20_mg24 | median | 2.5e-05 | 7.3e-04 | 1.7e-04 | 1.4e-04 | 1.5e-04 | 6.2e-05 |
| skye | mix_o16_ne20_mg24 | p99 | 4.0e-04 | 4.8e-02 | 4.6e-02 | 1.3e+00 | 4.9e-01 | 1.3e+00 |
| skye | mix_o16_ne20_mg24 | max | 1.5e-03 | 5.0e-02 | 1.5e-01 | 2.6e+00 | 9.2e-01 | 2.6e+00 |
| skye | mix_he4_c12 | median | 3.5e-05 | 1.2e-03 | 3.1e-04 | 1.1e-03 | 1.3e-03 | 7.1e-04 |
| skye | mix_he4_c12 | p99 | 6.6e-03 | 3.7e-01 | 2.6e-01 | 1.0e+01 | 3.0e-01 | 1.0e+01 |
| skye | mix_he4_c12 | max | 3.8e-02 | 4.4e-01 | 6.0e-01 | 5.1e+01 | 8.7e-01 | 5.1e+01 |
| skye | mix_h1_he4 | median | 5.0e-06 | 8.3e-05 | 3.4e-04 | 1.5e-04 | 6.7e-05 | 2.3e-05 |
| skye | mix_h1_he4 | p99 | 5.2e-02 | 1.5e-02 | 3.7e-02 | 5.9e+00 | 5.5e-01 | 4.1e+00 |
| skye | mix_h1_he4 | max | 2.1e-01 | 1.9e-02 | 4.2e-02 | 5.0e+02 | 1.0e+03 | 3.4e+01 |
| skye | mix_c12_o16_ne22 | median | 1.6e-05 | 4.3e-04 | 5.2e-05 | 1.3e-04 | 1.5e-04 | 5.9e-05 |
| skye | mix_c12_o16_ne22 | p99 | 3.1e-04 | 5.0e-02 | 7.3e-02 | 1.8e+00 | 3.5e-01 | 1.8e+00 |
| skye | mix_c12_o16_ne22 | max | 1.5e-03 | 5.3e-02 | 3.2e-01 | 3.5e+00 | 6.1e-01 | 3.6e+00 |
| helm | mix_c12_o16 | median | 1.9e-04 | 2.4e-04 | 1.1e-02 | 1.5e-04 | 5.7e-05 | 5.1e-04 |
| helm | mix_c12_o16 | p99 | 2.8e-03 | 6.1e-03 | 2.0e-01 | 2.5e-03 | 1.9e-03 | 2.3e-02 |
| helm | mix_c12_o16 | max | 5.3e-03 | 1.1e-02 | 3.2e-01 | 1.5e-02 | 5.7e-02 | 1.3e-01 |
| helm | mix_o16_ne20_mg24 | median | 2.2e-04 | 2.8e-04 | 1.2e-02 | 1.7e-04 | 5.8e-05 | 4.7e-04 |
| helm | mix_o16_ne20_mg24 | p99 | 3.1e-03 | 8.7e-03 | 2.3e-01 | 3.6e-03 | 3.3e-03 | 2.4e-02 |
| helm | mix_o16_ne20_mg24 | max | 4.3e-03 | 1.6e-02 | 3.1e-01 | 3.5e-02 | 6.0e-02 | 3.3e-02 |
| helm | mix_he4_c12 | median | 7.8e-04 | 1.6e-03 | 8.6e-03 | 7.8e-04 | 6.0e-04 | 4.5e-03 |
| helm | mix_he4_c12 | p99 | 1.8e-02 | 4.0e-02 | 1.5e-01 | 2.0e-02 | 2.3e-02 | 3.1e-02 |
| helm | mix_he4_c12 | max | 2.5e-02 | 6.0e-02 | 2.1e-01 | 1.5e-01 | 2.5e-01 | 1.2e-01 |
| helm | mix_h1_he4 | median | 3.4e-05 | 7.0e-05 | 4.8e-03 | 6.7e-05 | 4.3e-05 | 7.9e-04 |
| helm | mix_h1_he4 | p99 | 1.9e-03 | 4.3e-03 | 7.0e-02 | 3.5e-02 | 2.1e-02 | 1.1e-01 |
| helm | mix_h1_he4 | max | 1.4e-02 | 6.2e-03 | 9.2e-02 | 2.1e-01 | 5.5e-01 | 1.5e-01 |
| helm | mix_c12_o16_ne22 | median | 1.9e-04 | 2.6e-04 | 1.1e-02 | 1.7e-04 | 6.1e-05 | 5.5e-04 |
| helm | mix_c12_o16_ne22 | p99 | 3.0e-03 | 6.7e-03 | 2.1e-01 | 2.9e-03 | 3.0e-03 | 2.3e-02 |
| helm | mix_c12_o16_ne22 | max | 4.4e-03 | 1.2e-02 | 3.4e-01 | 1.9e-02 | 6.6e-02 | 1.2e-01 |

## Entropy: F1 correction evidence (Skye)

Median relative error of S without and with the F1 correction (subtracting `Delta s_ideal` from the Skye reference). The correction removes the systematic entropy offset from Skye's double-counted ideal mixing entropy.

| mixture | Delta_s_ideal [erg/g/K] | S uncorr (skye_valid) | S corr (skye_valid) | S uncorr (target) | S corr (target) |
|---|---|---|---|---|---|
| mix_c12_o16 | 4.141e+06 | 4.3e-02 | 3.6e-05 | 2.9e-02 | 2.5e-04 |
| mix_o16_ne20_mg24 | 3.992e+06 | 5.8e-02 | 1.7e-04 | 2.6e-02 | 2.3e-04 |
| mix_he4_c12 | 7.791e+06 | 2.6e-02 | 3.1e-04 | 2.9e-02 | 2.2e-03 |
| mix_h1_he4 | 2.043e+07 | 1.3e-02 | 3.4e-04 | 2.9e-02 | 3.3e-04 |
| mix_c12_o16_ne22 | 4.464e+06 | 4.6e-02 | 5.2e-05 | 2.8e-02 | 2.2e-04 |

## HELM: interpolation floor vs observed error

Pure-species held-out interpolation median (mass-weighted over the mixture, a_o5 variant) vs the observed AV mixture median, so the AV model error can be isolated. A `*` marks a floor taken from the Skye held-out as a proxy (no HELM held-out table for that isotope).

### target region (logRho in [-2,4], logT in [3,7])

| mixture | source | P | E | S | Cv | chiRho | chiT |
|---|---|---|---|---|---|---|---|
| mix_c12_o16 | interp floor | 1.2e-09 | 2.6e-08 | 1.7e-09 | 1.3e-05 | 5.3e-08 | 3.2e-08 |
|  | observed AV | 1.9e-03 | 4.3e-03 | 1.3e-02 | 3.5e-04 | 3.7e-04 | 1.7e-03 |
| mix_o16_ne20_mg24 | interp floor | 8.6e-10* | 1.5e-08* | 1.8e-08* | 1.5e-05* | 8.0e-08* | 2.7e-07* |
|  | observed AV | 1.8e-03 | 4.2e-03 | 1.1e-02 | 2.8e-04 | 3.6e-04 | 1.6e-03 |
| mix_he4_c12 | interp floor | 1.3e-09* | 1.7e-08* | 4.7e-09* | 1.5e-05* | 9.1e-08* | 9.6e-08* |
|  | observed AV | 1.5e-02 | 3.3e-02 | 1.2e-02 | 2.0e-03 | 3.0e-03 | 1.2e-02 |
| mix_h1_he4 | interp floor | 1.3e-09* | 2.0e-08* | 8.1e-09* | 1.7e-05* | 1.2e-07* | 1.6e-07* |
|  | observed AV | 1.9e-03 | 4.1e-03 | 8.7e-03 | 6.9e-04 | 3.6e-04 | 1.5e-03 |
| mix_c12_o16_ne22 | interp floor | 1.2e-09* | 2.5e-08* | 2.2e-09* | 1.3e-05* | 5.5e-08* | 4.0e-08* |
|  | observed AV | 1.9e-03 | 4.1e-03 | 1.3e-02 | 2.9e-04 | 3.7e-04 | 1.6e-03 |

### Skye-valid region (logT >= 6.5 OR logRho >= 3.7)

| mixture | source | P | E | S | Cv | chiRho | chiT |
|---|---|---|---|---|---|---|---|
| mix_c12_o16 | interp floor | 2.7e-09 | 4.3e-08 | 9.1e-09 | 2.3e-05 | 2.4e-08 | 1.0e-07 |
|  | observed AV | 1.9e-04 | 2.4e-04 | 1.1e-02 | 5.7e-05 | 1.3e-05 | 1.5e-04 |
| mix_o16_ne20_mg24 | interp floor | 1.5e-09* | 2.5e-08* | 9.0e-08* | 2.7e-05* | 2.3e-08* | 4.8e-06* |
|  | observed AV | 2.2e-04 | 2.8e-04 | 1.2e-02 | 5.8e-05 | 1.6e-05 | 1.7e-04 |
| mix_he4_c12 | interp floor | 1.6e-09* | 2.6e-08* | 1.6e-08* | 2.5e-05* | 1.8e-08* | 9.5e-07* |
|  | observed AV | 7.8e-04 | 1.6e-03 | 8.6e-03 | 6.0e-04 | 7.3e-05 | 7.8e-04 |
| mix_h1_he4 | interp floor | 4.3e-10* | 8.8e-09* | 1.8e-08* | 2.6e-05* | 1.1e-08* | 2.1e-06* |
|  | observed AV | 3.4e-05 | 7.0e-05 | 4.8e-03 | 4.3e-05 | 5.8e-06 | 6.7e-05 |
| mix_c12_o16_ne22 | interp floor | 2.6e-09* | 4.2e-08* | 1.3e-08* | 2.3e-05* | 2.4e-08* | 3.3e-07* |
|  | observed AV | 1.9e-04 | 2.6e-04 | 1.1e-02 | 6.1e-05 | 1.5e-05 | 1.7e-04 |

## Coverage (integrated-path declines)

`eosDT_get_component(i_eos_WDM_AV)` calls the component directly (no blend box); a node is declined when the AV solver fails (no-converge / bad-domain / out-of-table), in which case WDM_AV returns the all-zero sentinel that `mix_vs_reference` records as flag = -1.

| set | mixture | both-valid / total | both-valid in target | ref-valid but WDM declined |
|---|---|---|---|---|
| skye | mix_c12_o16 | 5912/11946 | 934/2501 | 4496 |
| skye | mix_o16_ne20_mg24 | 5657/11946 | 761/2501 | 4776 |
| skye | mix_he4_c12 | 5799/11946 | 1021/2501 | 5108 |
| skye | mix_h1_he4 | 7179/11946 | 1534/2501 | 4271 |
| skye | mix_c12_o16_ne22 | 5837/11946 | 850/2501 | 4572 |
| helm | mix_c12_o16 | 4585/11946 | 630/2501 | 3499 |
| helm | mix_o16_ne20_mg24 | 4224/11946 | 504/2501 | 3669 |
| helm | mix_he4_c12 | 4843/11946 | 673/2501 | 3795 |
| helm | mix_h1_he4 | 5718/11946 | 1088/2501 | 3209 |
| helm | mix_c12_o16_ne22 | 4424/11946 | 548/2501 | 3657 |
