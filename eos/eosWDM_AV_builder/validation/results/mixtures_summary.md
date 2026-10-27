# Mixture-level WDM_AV validation vs Skye / HELM

Relative error of the additive-volume WDM_AV mixture (integrated eosDT path, `eosDT_get_component(i_eos_WDM_AV)`) against Skye / HELM run directly on the same five reference mixtures, at the reference grid nodes where both sides are finite. `S` is the F1-corrected entropy for Skye (reference minus `Delta s_ideal`); for HELM, radiation is stripped from WDM_AV's e/s/Cv and gas-only chiRho/chiT come from WDM_AV's dlnPgas columns. HELM grad_ad/gamma1 mix a radiation-definition difference and are taken only where Prad/Pgas < 1e-2. Full per-mixture, per-region data (p90, p99, max, all quantities) is in the CSV/JSON.

## Aggregate median relative error (median over the 5 mixtures)

### (a) Skye

| region | P | E | S | chiT | Cv | grad_ad |
|---|---|---|---|---|---|---|
| whole_valid | 2.5e-05 | 2.7e-05 | 6.7e-05 | 1.3e-04 | 1.4e-04 | 5.3e-05 |
| skye_valid | 1.3e-05 | 2.7e-05 | 1.2e-04 | 1.2e-04 | 1.4e-04 | 4.7e-05 |
| skye_valid_no_pair | 1.6e-05 | 2.8e-05 | 1.5e-04 | 1.4e-04 | 1.7e-04 | 6.7e-05 |
| target | 2.7e-04 | 6.2e-04 | 2.4e-04 | 9.1e-04 | 4.7e-04 | 6.4e-04 |

### (b) HELM

| region | P | E | S | chiT | Cv | grad_ad |
|---|---|---|---|---|---|---|
| whole_valid | 2.2e-04 | 3.1e-04 | 2.0e-04 | 2.2e-04 | 6.2e-05 | 8.5e-04 |
| skye_valid | 2.1e-04 | 2.6e-04 | 2.7e-04 | 1.6e-04 | 6.1e-05 | 5.5e-04 |
| skye_valid_no_pair | 2.2e-04 | 2.8e-04 | 3.2e-04 | 1.7e-04 | 6.2e-05 | 5.5e-04 |
| target | 1.9e-03 | 4.1e-03 | 3.4e-04 | 1.6e-03 | 3.0e-04 | 1.1e-03 |

### (c) Skye reference, tables built from P, E, F only (FD second derivatives)

| region | P | E | S | chiT | Cv | grad_ad |
|---|---|---|---|---|---|---|
| whole_valid | 2.5e-05 | 2.7e-05 | 7.0e-05 | 1.6e-04 | 1.8e-04 | 1.3e-04 |
| skye_valid | 1.3e-05 | 2.7e-05 | 1.3e-04 | 1.6e-04 | 2.1e-04 | 1.3e-04 |
| skye_valid_no_pair | 1.5e-05 | 2.8e-05 | 1.5e-04 | 1.8e-04 | 2.5e-04 | 1.4e-04 |
| target | 2.6e-04 | 6.1e-04 | 2.3e-04 | 9.2e-04 | 4.2e-04 | 6.9e-04 |

## Per-mixture, target region (logRho in [-2,4], logT in [3,7])

| set | mixture | stat | P | E | S | chiT | Cv | grad_ad |
|---|---|---|---|---|---|---|---|---|
| skye | mix_c12_o16 | median | 2.7e-04 | 6.2e-04 | 2.4e-04 | 8.7e-04 | 4.4e-04 | 5.6e-04 |
| skye | mix_c12_o16 | p99 | 2.2e-02 | 4.6e-02 | 4.7e-02 | 4.4e-01 | 3.1e-01 | 4.2e-01 |
| skye | mix_c12_o16 | max | 5.0e-01 | 2.0e-01 | 1.5e-01 | 2.9e+00 | 2.8e+00 | 6.5e+00 |
| skye | mix_o16_ne20_mg24 | median | 2.7e-04 | 6.8e-04 | 2.1e-04 | 1.0e-03 | 4.7e-04 | 6.4e-04 |
| skye | mix_o16_ne20_mg24 | p99 | 1.3e-02 | 3.6e-02 | 5.0e-02 | 1.6e-01 | 3.2e-01 | 2.0e-01 |
| skye | mix_o16_ne20_mg24 | max | 1.3e-01 | 1.2e-01 | 9.3e-02 | 4.0e+00 | 1.8e+00 | 1.9e+00 |
| skye | mix_he4_c12 | median | 7.7e-04 | 2.5e-03 | 1.8e-03 | 5.4e-03 | 6.9e-03 | 5.3e-03 |
| skye | mix_he4_c12 | p99 | 1.1e-01 | 8.6e-02 | 2.1e-01 | 8.3e+00 | 4.3e-01 | 8.2e+00 |
| skye | mix_he4_c12 | max | 2.6e-01 | 1.3e-01 | 2.5e-01 | 2.2e+01 | 8.2e-01 | 2.2e+01 |
| skye | mix_h1_he4 | median | 1.1e-04 | 3.7e-04 | 3.2e-04 | 9.1e-04 | 1.1e-03 | 8.8e-04 |
| skye | mix_h1_he4 | p99 | 2.8e-02 | 4.1e-02 | 3.9e-02 | 3.6e+00 | 3.2e-01 | 3.4e+00 |
| skye | mix_h1_he4 | max | 8.8e-02 | 8.6e-02 | 2.5e-01 | 8.8e+00 | 2.4e+00 | 8.8e+00 |
| skye | mix_c12_o16_ne22 | median | 2.5e-04 | 6.0e-04 | 2.0e-04 | 7.8e-04 | 4.4e-04 | 5.1e-04 |
| skye | mix_c12_o16_ne22 | p99 | 1.5e-02 | 4.1e-02 | 4.8e-02 | 2.9e-01 | 2.8e-01 | 2.5e-01 |
| skye | mix_c12_o16_ne22 | max | 2.2e-01 | 1.5e-01 | 1.1e-01 | 2.1e+00 | 1.4e+00 | 1.6e+00 |
| helm | mix_c12_o16 | median | 1.8e-03 | 4.1e-03 | 3.1e-04 | 1.6e-03 | 3.0e-04 | 1.0e-03 |
| helm | mix_c12_o16 | p99 | 6.0e-03 | 2.4e-02 | 9.0e-04 | 7.8e-03 | 4.7e-03 | 2.3e-02 |
| helm | mix_c12_o16 | max | 6.8e-03 | 3.9e-02 | 1.8e-03 | 9.0e-03 | 5.9e-03 | 2.3e-02 |
| helm | mix_o16_ne20_mg24 | median | 1.8e-03 | 4.1e-03 | 2.1e-04 | 1.6e-03 | 2.4e-04 | 1.3e-03 |
| helm | mix_o16_ne20_mg24 | p99 | 5.0e-03 | 2.0e-02 | 5.8e-04 | 5.6e-03 | 3.1e-03 | 2.4e-02 |
| helm | mix_o16_ne20_mg24 | max | 5.3e-03 | 2.3e-02 | 1.5e-03 | 7.1e-03 | 4.5e-03 | 2.4e-02 |
| helm | mix_he4_c12 | median | 1.5e-02 | 3.3e-02 | 6.7e-03 | 1.2e-02 | 2.1e-03 | 3.9e-03 |
| helm | mix_he4_c12 | p99 | 5.2e-02 | 1.5e-01 | 1.4e-02 | 5.8e-02 | 2.4e-02 | 2.6e-02 |
| helm | mix_he4_c12 | max | 6.5e-02 | 1.8e-01 | 1.4e-02 | 8.4e-01 | 9.0e-02 | 5.3e+00 |
| helm | mix_h1_he4 | median | 2.0e-03 | 4.3e-03 | 1.2e-02 | 1.6e-03 | 6.9e-04 | 7.0e-04 |
| helm | mix_h1_he4 | p99 | 1.2e-02 | 3.3e-02 | 6.1e-02 | 1.4e-02 | 5.3e-03 | 2.6e-02 |
| helm | mix_h1_he4 | max | 4.0e-02 | 4.2e-02 | 7.2e-02 | 8.2e-02 | 1.4e-02 | 2.5e-01 |
| helm | mix_c12_o16_ne22 | median | 1.9e-03 | 4.1e-03 | 3.4e-04 | 1.6e-03 | 2.6e-04 | 1.1e-03 |
| helm | mix_c12_o16_ne22 | p99 | 5.2e-03 | 1.9e-02 | 8.0e-04 | 5.8e-03 | 2.7e-03 | 2.3e-02 |
| helm | mix_c12_o16_ne22 | max | 5.6e-03 | 2.3e-02 | 1.0e-03 | 6.5e-03 | 3.3e-03 | 2.3e-02 |
| skye_fd | mix_c12_o16 | median | 2.6e-04 | 6.1e-04 | 2.3e-04 | 9.0e-04 | 3.5e-04 | 6.1e-04 |
| skye_fd | mix_c12_o16 | p99 | 2.2e-02 | 3.9e-02 | 4.7e-02 | 4.7e-01 | 1.2e+00 | 4.9e-01 |
| skye_fd | mix_c12_o16 | max | 5.0e-01 | 2.1e-01 | 1.5e-01 | 2.6e+00 | 3.1e+00 | 1.5e+00 |
| skye_fd | mix_o16_ne20_mg24 | median | 2.7e-04 | 6.8e-04 | 2.1e-04 | 1.0e-03 | 4.2e-04 | 6.9e-04 |
| skye_fd | mix_o16_ne20_mg24 | p99 | 1.3e-02 | 3.6e-02 | 5.0e-02 | 3.2e-01 | 1.2e+00 | 5.2e-01 |
| skye_fd | mix_o16_ne20_mg24 | max | 1.3e-01 | 1.3e-01 | 9.6e-02 | 3.7e+00 | 2.5e+00 | 2.4e+00 |
| skye_fd | mix_he4_c12 | median | 7.7e-04 | 2.5e-03 | 1.8e-03 | 5.8e-03 | 7.2e-03 | 5.8e-03 |
| skye_fd | mix_he4_c12 | p99 | 1.1e-01 | 8.3e-02 | 2.1e-01 | 8.3e+00 | 1.5e+00 | 8.2e+00 |
| skye_fd | mix_he4_c12 | max | 2.6e-01 | 1.2e-01 | 2.5e-01 | 2.2e+01 | 2.1e+00 | 2.2e+01 |
| skye_fd | mix_h1_he4 | median | 1.1e-04 | 3.7e-04 | 3.2e-04 | 9.2e-04 | 1.4e-03 | 9.5e-04 |
| skye_fd | mix_h1_he4 | p99 | 2.6e-02 | 3.9e-02 | 3.9e-02 | 3.6e+00 | 9.4e-01 | 3.1e+00 |
| skye_fd | mix_h1_he4 | max | 8.8e-02 | 8.6e-02 | 2.5e-01 | 8.7e+00 | 3.1e+00 | 8.8e+00 |
| skye_fd | mix_c12_o16_ne22 | median | 2.5e-04 | 5.9e-04 | 2.0e-04 | 8.2e-04 | 3.8e-04 | 5.6e-04 |
| skye_fd | mix_c12_o16_ne22 | p99 | 1.5e-02 | 3.9e-02 | 4.8e-02 | 3.7e-01 | 1.2e+00 | 5.0e-01 |
| skye_fd | mix_c12_o16_ne22 | max | 2.2e-01 | 1.5e-01 | 1.1e-01 | 2.1e+00 | 3.0e+00 | 1.4e+00 |

## Per-mixture, Skye-valid region (logT >= 6.5 OR logRho >= 3.7)

| set | mixture | stat | P | E | S | chiT | Cv | grad_ad |
|---|---|---|---|---|---|---|---|---|
| skye | mix_c12_o16 | median | 1.2e-05 | 2.6e-05 | 3.3e-05 | 1.1e-04 | 1.3e-04 | 4.2e-05 |
| skye | mix_c12_o16 | p99 | 2.7e-04 | 6.7e-04 | 7.4e-02 | 1.6e+00 | 4.4e-01 | 1.6e+00 |
| skye | mix_c12_o16 | max | 2.0e-03 | 1.9e-03 | 3.2e-01 | 3.3e+00 | 9.0e-01 | 3.3e+00 |
| skye | mix_o16_ne20_mg24 | median | 1.9e-05 | 2.8e-05 | 1.2e-04 | 1.1e-04 | 1.2e-04 | 4.3e-05 |
| skye | mix_o16_ne20_mg24 | p99 | 3.6e-04 | 9.0e-04 | 4.1e-02 | 1.3e+00 | 4.9e-01 | 1.3e+00 |
| skye | mix_o16_ne20_mg24 | max | 1.5e-03 | 3.6e-03 | 1.6e-01 | 2.6e+00 | 9.2e-01 | 2.6e+00 |
| skye | mix_he4_c12 | median | 4.1e-05 | 4.9e-05 | 2.9e-04 | 1.1e-03 | 1.2e-03 | 5.9e-04 |
| skye | mix_he4_c12 | p99 | 7.7e-03 | 4.3e-03 | 2.6e-01 | 1.0e+01 | 2.9e-01 | 1.0e+01 |
| skye | mix_he4_c12 | max | 3.8e-02 | 9.9e-03 | 6.0e-01 | 5.1e+01 | 8.7e-01 | 5.1e+01 |
| skye | mix_h1_he4 | median | 6.7e-06 | 6.9e-06 | 3.9e-04 | 2.5e-04 | 1.4e-04 | 4.9e-05 |
| skye | mix_h1_he4 | p99 | 4.3e-02 | 4.4e-03 | 3.7e-02 | 6.1e+00 | 6.9e-01 | 4.6e+00 |
| skye | mix_h1_he4 | max | 2.1e-01 | 1.8e-02 | 4.2e-02 | 5.0e+02 | 1.0e+03 | 3.4e+01 |
| skye | mix_c12_o16_ne22 | median | 1.3e-05 | 2.7e-05 | 4.3e-05 | 1.2e-04 | 1.4e-04 | 4.7e-05 |
| skye | mix_c12_o16_ne22 | p99 | 3.0e-04 | 7.5e-04 | 7.3e-02 | 1.8e+00 | 3.5e-01 | 1.8e+00 |
| skye | mix_c12_o16_ne22 | max | 2.4e-03 | 2.2e-03 | 3.2e-01 | 3.5e+00 | 6.0e-01 | 3.6e+00 |
| helm | mix_c12_o16 | median | 1.9e-04 | 2.4e-04 | 2.4e-04 | 1.6e-04 | 5.7e-05 | 5.1e-04 |
| helm | mix_c12_o16 | p99 | 2.8e-03 | 6.4e-03 | 1.4e-03 | 2.7e-03 | 1.9e-03 | 2.3e-02 |
| helm | mix_c12_o16 | max | 5.3e-03 | 1.1e-02 | 1.6e-03 | 1.5e-02 | 5.7e-02 | 1.3e-01 |
| helm | mix_o16_ne20_mg24 | median | 2.2e-04 | 2.7e-04 | 1.7e-04 | 1.6e-04 | 5.0e-05 | 4.7e-04 |
| helm | mix_o16_ne20_mg24 | p99 | 3.0e-03 | 8.5e-03 | 1.0e-03 | 3.6e-03 | 3.8e-03 | 2.4e-02 |
| helm | mix_o16_ne20_mg24 | max | 4.3e-03 | 1.6e-02 | 1.6e-03 | 3.5e-02 | 6.0e-02 | 3.3e-02 |
| helm | mix_he4_c12 | median | 1.2e-03 | 2.0e-03 | 5.2e-03 | 1.1e-03 | 8.1e-04 | 4.6e-03 |
| helm | mix_he4_c12 | p99 | 1.8e-02 | 4.1e-02 | 3.5e-02 | 2.1e-02 | 2.3e-02 | 3.1e-02 |
| helm | mix_he4_c12 | max | 3.1e-02 | 6.0e-02 | 5.0e-02 | 2.1e-01 | 2.5e-01 | 2.1e-01 |
| helm | mix_h1_he4 | median | 8.0e-05 | 1.5e-04 | 7.4e-03 | 1.3e-04 | 7.9e-05 | 7.9e-04 |
| helm | mix_h1_he4 | p99 | 1.9e-03 | 4.3e-03 | 7.6e-02 | 4.5e-02 | 2.7e-02 | 1.1e-01 |
| helm | mix_h1_he4 | max | 1.4e-02 | 6.2e-03 | 1.0e-01 | 2.1e-01 | 5.5e-01 | 2.2e-01 |
| helm | mix_c12_o16_ne22 | median | 2.1e-04 | 2.6e-04 | 2.7e-04 | 1.7e-04 | 6.1e-05 | 5.5e-04 |
| helm | mix_c12_o16_ne22 | p99 | 3.1e-03 | 7.0e-03 | 1.7e-03 | 3.1e-03 | 3.7e-03 | 2.3e-02 |
| helm | mix_c12_o16_ne22 | max | 4.4e-03 | 1.2e-02 | 1.9e-03 | 1.9e-02 | 6.6e-02 | 1.2e-01 |
| skye_fd | mix_c12_o16 | median | 1.2e-05 | 2.6e-05 | 3.4e-05 | 1.4e-04 | 2.0e-04 | 1.3e-04 |
| skye_fd | mix_c12_o16 | p99 | 2.6e-04 | 6.7e-04 | 7.6e-02 | 1.8e+00 | 1.1e+00 | 2.1e+00 |
| skye_fd | mix_c12_o16 | max | 7.2e-04 | 1.9e-03 | 3.2e-01 | 3.4e+00 | 1.5e+00 | 1.2e+02 |
| skye_fd | mix_o16_ne20_mg24 | median | 1.9e-05 | 2.8e-05 | 1.3e-04 | 1.4e-04 | 2.0e-04 | 1.3e-04 |
| skye_fd | mix_o16_ne20_mg24 | p99 | 3.5e-04 | 9.2e-04 | 4.3e-02 | 1.5e+00 | 1.1e+00 | 1.7e+00 |
| skye_fd | mix_o16_ne20_mg24 | max | 1.5e-03 | 3.6e-03 | 1.6e-01 | 2.8e+00 | 1.6e+00 | 2.2e+01 |
| skye_fd | mix_he4_c12 | median | 4.1e-05 | 5.5e-05 | 3.3e-04 | 1.3e-03 | 8.1e-04 | 8.8e-04 |
| skye_fd | mix_he4_c12 | p99 | 7.7e-03 | 4.3e-03 | 2.6e-01 | 1.0e+01 | 6.6e-01 | 1.0e+01 |
| skye_fd | mix_he4_c12 | max | 3.9e-02 | 9.9e-03 | 6.0e-01 | 5.1e+01 | 3.1e+00 | 5.1e+01 |
| skye_fd | mix_h1_he4 | median | 6.9e-06 | 8.2e-06 | 3.9e-04 | 3.0e-04 | 2.8e-04 | 9.6e-05 |
| skye_fd | mix_h1_he4 | p99 | 4.1e-02 | 4.1e-03 | 3.7e-02 | 5.4e+00 | 2.9e-01 | 4.7e+00 |
| skye_fd | mix_h1_he4 | max | 2.1e-01 | 1.7e-02 | 4.2e-02 | 3.4e+02 | 7.1e+02 | 1.6e+01 |
| skye_fd | mix_c12_o16_ne22 | median | 1.3e-05 | 2.7e-05 | 4.5e-05 | 1.6e-04 | 2.1e-04 | 1.4e-04 |
| skye_fd | mix_c12_o16_ne22 | p99 | 2.9e-04 | 7.5e-04 | 7.3e-02 | 1.9e+00 | 1.0e+00 | 2.3e+00 |
| skye_fd | mix_c12_o16_ne22 | max | 8.3e-04 | 2.2e-03 | 3.1e-01 | 3.6e+00 | 1.5e+00 | 2.9e+01 |

## Reference corrections: primary vs alternative

Corrections before comparison. (1) Skye S: subtract `Delta s_ideal` (Finding F1, Skye double-counts the ideal mixing entropy). (2) HELM S: ADD `Delta s_ideal` (HELM's mean-ion treatment carries no ideal mixing entropy, while WDM_AV adds it once). (3) Skye E: the ion-offset bug is fixed upstream (PR #1097), so the Skye mixture E is now additive and the comparison uses the RAW Skye E by default; the old `dE` offset correction is retained but OFF. The `primary` column is the reported comparison; `alt` is the other convention (S: uncorrected; Skye E: dE-offset applied; HELM S: raw). Median relative error.

| set | mixture | quantity | Delta_s / dE | primary (skye_valid) | alt (skye_valid) | primary (target) | alt (target) |
|---|---|---|---|---|---|---|---|
| skye | mix_c12_o16 | S | 4.141e+06 | 3.3e-05 | 3.9e-02 | 2.4e-04 | 2.5e-02 |
| skye | mix_c12_o16 | E | 6.979e+13 | 2.6e-05 | 4.0e-04 | 6.2e-04 | 1.8e-01 |
| skye | mix_o16_ne20_mg24 | S | 3.992e+06 | 1.2e-04 | 4.8e-02 | 2.1e-04 | 2.3e-02 |
| skye | mix_o16_ne20_mg24 | E | 1.377e+14 | 2.8e-05 | 6.7e-04 | 6.8e-04 | 2.7e-01 |
| skye | mix_he4_c12 | S | 7.791e+06 | 2.9e-04 | 2.4e-02 | 1.8e-03 | 2.7e-02 |
| skye | mix_he4_c12 | E | 2.293e+14 | 4.9e-05 | 1.5e-03 | 2.5e-03 | 4.2e-01 |
| skye | mix_h1_he4 | S | 2.043e+07 | 3.9e-04 | 1.5e-02 | 3.2e-04 | 2.9e-02 |
| skye | mix_h1_he4 | E | 1.279e+13 | 6.9e-06 | 7.3e-05 | 3.7e-04 | 2.9e-02 |
| skye | mix_c12_o16_ne22 | S | 4.464e+06 | 4.3e-05 | 4.2e-02 | 2.0e-04 | 2.4e-02 |
| skye | mix_c12_o16_ne22 | E | 8.308e+13 | 2.7e-05 | 4.7e-04 | 6.0e-04 | 1.9e-01 |
| helm | mix_c12_o16 | S | 4.141e+06 | 1.1e-02 | 2.4e-04 | 1.2e-02 | 3.1e-04 |
| helm | mix_o16_ne20_mg24 | S | 3.992e+06 | 1.1e-02 | 1.7e-04 | 1.1e-02 | 2.1e-04 |
| helm | mix_he4_c12 | S | 7.791e+06 | 9.7e-03 | 5.2e-03 | 1.1e-02 | 6.7e-03 |
| helm | mix_h1_he4 | S | 2.043e+07 | 5.7e-03 | 7.4e-03 | 8.8e-03 | 1.2e-02 |
| helm | mix_c12_o16_ne22 | S | 4.464e+06 | 1.1e-02 | 2.7e-04 | 1.2e-02 | 3.4e-04 |
| skye | mix_c12_o16 | S | 4.141e+06 | 3.4e-05 | 4.1e-02 | 2.3e-04 | 2.5e-02 |
| skye | mix_c12_o16 | E | 6.979e+13 | 2.6e-05 | 4.4e-04 | 6.1e-04 | 1.8e-01 |
| skye | mix_o16_ne20_mg24 | S | 3.992e+06 | 1.3e-04 | 5.0e-02 | 2.1e-04 | 2.2e-02 |
| skye | mix_o16_ne20_mg24 | E | 1.377e+14 | 2.8e-05 | 7.5e-04 | 6.8e-04 | 2.7e-01 |
| skye | mix_he4_c12 | S | 7.791e+06 | 3.3e-04 | 2.5e-02 | 1.8e-03 | 2.7e-02 |
| skye | mix_he4_c12 | E | 2.293e+14 | 5.5e-05 | 1.6e-03 | 2.5e-03 | 4.2e-01 |
| skye | mix_h1_he4 | S | 2.043e+07 | 3.9e-04 | 1.5e-02 | 3.2e-04 | 2.9e-02 |
| skye | mix_h1_he4 | E | 1.279e+13 | 8.2e-06 | 8.3e-05 | 3.7e-04 | 2.9e-02 |
| skye | mix_c12_o16_ne22 | S | 4.464e+06 | 4.5e-05 | 4.3e-02 | 2.0e-04 | 2.4e-02 |
| skye | mix_c12_o16_ne22 | E | 8.308e+13 | 2.7e-05 | 5.2e-04 | 5.9e-04 | 1.9e-01 |

## HELM: interpolation floor vs observed error

Pure-species held-out interpolation median (mass-weighted over the mixture, a_o5 variant) vs the observed AV mixture median, so the AV model error can be isolated. A `*` marks a floor taken from the Skye held-out as a proxy (no HELM held-out table for that isotope).

### target region (logRho in [-2,4], logT in [3,7])

| mixture | source | P | E | S | Cv | chiRho | chiT |
|---|---|---|---|---|---|---|---|
| mix_c12_o16 | interp floor | 1.2e-09 | 2.6e-08 | 1.7e-09 | 1.3e-05 | 5.3e-08 | 3.2e-08 |
|  | observed AV | 1.8e-03 | 4.1e-03 | 3.1e-04 | 3.0e-04 | 3.8e-04 | 1.6e-03 |
| mix_o16_ne20_mg24 | interp floor | 8.6e-10* | 1.7e-08* | 1.9e-08* | 1.5e-05* | 8.0e-08* | 2.7e-07* |
|  | observed AV | 1.8e-03 | 4.1e-03 | 2.1e-04 | 2.4e-04 | 4.0e-04 | 1.6e-03 |
| mix_he4_c12 | interp floor | 1.3e-09* | 2.2e-08* | 4.7e-09* | 1.5e-05* | 9.1e-08* | 9.6e-08* |
|  | observed AV | 1.5e-02 | 3.3e-02 | 6.7e-03 | 2.1e-03 | 2.9e-03 | 1.2e-02 |
| mix_h1_he4 | interp floor | 1.3e-09* | 2.3e-08* | 8.1e-09* | 1.7e-05* | 1.2e-07* | 1.6e-07* |
|  | observed AV | 2.0e-03 | 4.3e-03 | 1.2e-02 | 6.9e-04 | 3.6e-04 | 1.6e-03 |
| mix_c12_o16_ne22 | interp floor | 1.2e-09* | 2.5e-08* | 2.3e-09* | 1.3e-05* | 5.5e-08* | 4.0e-08* |
|  | observed AV | 1.9e-03 | 4.1e-03 | 3.4e-04 | 2.6e-04 | 3.9e-04 | 1.6e-03 |

### Skye-valid region (logT >= 6.5 OR logRho >= 3.7)

| mixture | source | P | E | S | Cv | chiRho | chiT |
|---|---|---|---|---|---|---|---|
| mix_c12_o16 | interp floor | 2.7e-09 | 4.3e-08 | 9.1e-09 | 2.3e-05 | 2.4e-08 | 1.0e-07 |
|  | observed AV | 1.9e-04 | 2.4e-04 | 2.4e-04 | 5.7e-05 | 1.4e-05 | 1.6e-04 |
| mix_o16_ne20_mg24 | interp floor | 1.5e-09* | 2.5e-08* | 9.0e-08* | 2.7e-05* | 2.3e-08* | 4.8e-06* |
|  | observed AV | 2.2e-04 | 2.7e-04 | 1.7e-04 | 5.0e-05 | 1.7e-05 | 1.6e-04 |
| mix_he4_c12 | interp floor | 1.6e-09* | 2.6e-08* | 1.5e-08* | 2.5e-05* | 1.8e-08* | 9.3e-07* |
|  | observed AV | 1.2e-03 | 2.0e-03 | 5.2e-03 | 8.1e-04 | 1.1e-04 | 1.1e-03 |
| mix_h1_he4 | interp floor | 4.3e-10* | 8.8e-09* | 1.8e-08* | 2.6e-05* | 1.1e-08* | 2.1e-06* |
|  | observed AV | 8.0e-05 | 1.5e-04 | 7.4e-03 | 7.9e-05 | 9.5e-06 | 1.3e-04 |
| mix_c12_o16_ne22 | interp floor | 2.6e-09* | 4.2e-08* | 1.3e-08* | 2.3e-05* | 2.4e-08* | 3.3e-07* |
|  | observed AV | 2.1e-04 | 2.6e-04 | 2.7e-04 | 6.1e-05 | 1.6e-05 | 1.7e-04 |

## Coverage (integrated-path declines)

`eosDT_get_component(i_eos_WDM_AV)` calls the component directly (no blend box); a node is declined when the AV solver fails (no-converge / bad-domain / out-of-table), in which case WDM_AV returns the all-zero sentinel that `mix_vs_reference` records as flag = -1.

| set | mixture | both-valid / total | both-valid in target | ref-valid but WDM declined |
|---|---|---|---|---|
| skye | mix_c12_o16 | 8396/11946 | 1232/2501 | 2012 |
| skye | mix_o16_ne20_mg24 | 8251/11946 | 1017/2501 | 2182 |
| skye | mix_he4_c12 | 8030/11946 | 1390/2501 | 2877 |
| skye | mix_h1_he4 | 8721/11946 | 1921/2501 | 2729 |
| skye | mix_c12_o16_ne22 | 8230/11946 | 1109/2501 | 2179 |
| helm | mix_c12_o16 | 6075/11946 | 840/2501 | 2009 |
| helm | mix_o16_ne20_mg24 | 5778/11946 | 681/2501 | 2115 |
| helm | mix_he4_c12 | 6267/11946 | 945/2501 | 2371 |
| helm | mix_h1_he4 | 6754/11946 | 1359/2501 | 2173 |
| helm | mix_c12_o16_ne22 | 5925/11946 | 753/2501 | 2156 |
| skye_fd | mix_c12_o16 | 8257/11946 | 1228/2501 | 2151 |
| skye_fd | mix_o16_ne20_mg24 | 8088/11946 | 1016/2501 | 2345 |
| skye_fd | mix_he4_c12 | 7883/11946 | 1387/2501 | 3024 |
| skye_fd | mix_h1_he4 | 8511/11946 | 1920/2501 | 2939 |
| skye_fd | mix_c12_o16_ne22 | 8070/11946 | 1106/2501 | 2339 |
