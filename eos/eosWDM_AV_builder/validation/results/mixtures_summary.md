# Mixture-level WDM_AV validation vs Skye / HELM

Relative error of the additive-volume WDM_AV mixture (integrated eosDT path, `eosDT_get_component(i_eos_WDM_AV)`) against Skye / HELM run directly on the same eleven reference mixtures, at the reference grid nodes where both sides are finite. `S` is the F1-corrected entropy for Skye (reference minus `Delta s_ideal`); for HELM, radiation is stripped from WDM_AV's e/s/Cv and gas-only chiRho/chiT come from WDM_AV's dlnPgas columns. HELM grad_ad/gamma1 mix a radiation-definition difference and are taken only where Prad/Pgas < 1e-2. Full per-mixture, per-region data (p90, p99, max, all quantities) is in the CSV/JSON.

## Aggregate median relative error (median over all mixtures of the set)

### (a) Skye

| region | P | E | S | chiT | Cv | grad_ad |
|---|---|---|---|---|---|---|
| whole_valid | 2.5e-05 | 2.7e-05 | 4.0e-05 | 1.3e-04 | 1.4e-04 | 5.3e-05 |
| skye_valid | 1.3e-05 | 2.7e-05 | 5.8e-05 | 1.2e-04 | 1.4e-04 | 4.7e-05 |
| skye_valid_no_pair | 1.6e-05 | 2.8e-05 | 7.0e-05 | 1.4e-04 | 1.7e-04 | 6.7e-05 |
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
| whole_valid | 2.5e-05 | 2.7e-05 | 4.4e-05 | 1.6e-04 | 1.8e-04 | 1.3e-04 |
| skye_valid | 1.3e-05 | 2.7e-05 | 5.9e-05 | 1.6e-04 | 2.1e-04 | 1.3e-04 |
| skye_valid_no_pair | 1.5e-05 | 2.8e-05 | 7.0e-05 | 1.8e-04 | 2.5e-04 | 1.4e-04 |
| target | 2.6e-04 | 6.1e-04 | 2.3e-04 | 9.2e-04 | 4.2e-04 | 6.9e-04 |

### (d) Skye, common-electron-density rule (D27)

| region | P | E | S | chiT | Cv | grad_ad |
|---|---|---|---|---|---|---|
| whole_valid | 1.3e-05 | 9.4e-06 | 2.0e-05 | 9.7e-06 | 2.2e-05 | 6.9e-06 |
| skye_valid | 8.7e-06 | 1.1e-05 | 3.3e-05 | 1.2e-05 | 1.8e-05 | 6.3e-06 |
| skye_valid_no_pair | 1.2e-05 | 1.8e-05 | 8.3e-05 | 1.3e-05 | 3.0e-05 | 1.0e-05 |
| target | 7.2e-05 | 1.6e-04 | 5.3e-05 | 1.9e-04 | 3.3e-05 | 6.2e-05 |

### (e) HELM, common-electron-density rule (D27)

| region | P | E | S | chiT | Cv | grad_ad |
|---|---|---|---|---|---|---|
| whole_valid | 1.2e-04 | 2.3e-04 | 1.5e-04 | 1.4e-04 | 4.6e-05 | 1.5e-01 |
| skye_valid | 3.1e-05 | 6.6e-05 | 1.6e-04 | 4.7e-05 | 2.3e-05 | 2.7e-01 |
| skye_valid_no_pair | 2.0e-04 | 2.5e-04 | 2.7e-04 | 1.3e-04 | 4.8e-05 | 9.4e-02 |
| target | 1.9e-03 | 4.4e-03 | 3.1e-04 | 1.5e-03 | 1.5e-04 | 8.9e-04 |

## Two mixing rules: additive volume vs common electron density (D27)

Median relative error against the same reference, additive volume (`pressure`) vs the common-electron-density rule (`electron_density`), for the two regions that matter. Lower is better; the common-n_e rule is expected to help against Skye (linear mixing reference) and to leave the HELM mean-ion Coulomb difference essentially unchanged.

### Skye: target region (logRho in [-2,4], logT in [3,7])

| rule | P | E | S | chiT | Cv |
|---|---|---|---|---|---|
| additive volume | 2.7e-04 | 6.2e-04 | 2.4e-04 | 9.1e-04 | 4.7e-04 |
| common n_e | 7.2e-05 | 1.6e-04 | 5.3e-05 | 1.9e-04 | 3.3e-05 |

### Skye: Skye-valid region (logT >= 6.5 OR logRho >= 3.7)

| rule | P | E | S | chiT | Cv |
|---|---|---|---|---|---|
| additive volume | 1.3e-05 | 2.7e-05 | 5.8e-05 | 1.2e-04 | 1.4e-04 |
| common n_e | 8.7e-06 | 1.1e-05 | 3.3e-05 | 1.2e-05 | 1.8e-05 |

### Skye: six worst mixtures (ranked by additive-volume target-gap P median), both rules

| mixture | region | rule | P | E | S | chiT | Cv |
|---|---|---|---|---|---|---|---|
| mix_he4_c12_o16 | target | AV | 8.6e-04 | 2.5e-03 | 1.6e-03 | 6.3e-03 | 7.3e-03 |
| mix_he4_c12_o16 | target | n_e | 5.5e-04 | 8.9e-04 | 2.6e-04 | 4.4e-04 | 6.8e-04 |
| mix_he4_c12_o16 | skye_valid | AV | 5.1e-05 | 6.3e-05 | 3.1e-04 | 1.1e-03 | 1.1e-03 |
| mix_he4_c12_o16 | skye_valid | n_e | 2.9e-05 | 3.7e-05 | 6.2e-05 | 1.2e-04 | 1.8e-04 |
| mix_he4_c12 | target | AV | 7.7e-04 | 2.5e-03 | 1.8e-03 | 5.4e-03 | 6.9e-03 |
| mix_he4_c12 | target | n_e | 4.8e-04 | 8.9e-04 | 3.8e-04 | 3.9e-04 | 7.7e-04 |
| mix_he4_c12 | skye_valid | AV | 4.1e-05 | 4.9e-05 | 2.9e-04 | 1.1e-03 | 1.2e-03 |
| mix_he4_c12 | skye_valid | n_e | 1.5e-05 | 1.9e-05 | 3.5e-05 | 1.3e-04 | 1.4e-04 |
| mix_he4_c12_10_90 | target | AV | 5.9e-04 | 1.5e-03 | 9.0e-04 | 2.8e-03 | 3.1e-03 |
| mix_he4_c12_10_90 | target | n_e | 3.6e-04 | 5.5e-04 | 2.8e-04 | 4.1e-04 | 2.7e-04 |
| mix_he4_c12_10_90 | skye_valid | AV | 4.7e-05 | 5.6e-05 | 1.5e-04 | 5.7e-04 | 5.6e-04 |
| mix_he4_c12_10_90 | skye_valid | n_e | 4.2e-05 | 4.8e-05 | 2.3e-05 | 5.1e-05 | 8.4e-05 |
| mix_c12_o16_ne20 | target | AV | 3.0e-04 | 7.8e-04 | 2.9e-04 | 1.2e-03 | 6.5e-04 |
| mix_c12_o16_ne20 | target | n_e | 7.2e-05 | 2.6e-04 | 5.3e-05 | 1.9e-04 | 3.7e-05 |
| mix_c12_o16_ne20 | skye_valid | AV | 2.0e-05 | 3.1e-05 | 5.8e-05 | 1.6e-04 | 1.9e-04 |
| mix_c12_o16_ne20 | skye_valid | n_e | 1.3e-05 | 1.6e-05 | 4.6e-05 | 1.2e-05 | 2.5e-05 |
| mix_o16_ne20_mg24 | target | AV | 2.7e-04 | 6.8e-04 | 2.1e-04 | 1.0e-03 | 4.7e-04 |
| mix_o16_ne20_mg24 | target | n_e | 2.1e-04 | 5.3e-04 | 6.9e-05 | 2.2e-04 | 1.6e-05 |
| mix_o16_ne20_mg24 | skye_valid | AV | 1.9e-05 | 2.8e-05 | 1.2e-04 | 1.1e-04 | 1.2e-04 |
| mix_o16_ne20_mg24 | skye_valid | n_e | 9.6e-06 | 1.3e-05 | 6.1e-05 | 1.2e-05 | 1.1e-05 |
| mix_c12_o16 | target | AV | 2.7e-04 | 6.2e-04 | 2.4e-04 | 8.7e-04 | 4.4e-04 |
| mix_c12_o16 | target | n_e | 5.0e-05 | 6.3e-05 | 3.3e-05 | 1.2e-04 | 3.0e-05 |
| mix_c12_o16 | skye_valid | AV | 1.2e-05 | 2.6e-05 | 3.3e-05 | 1.1e-04 | 1.3e-04 |
| mix_c12_o16 | skye_valid | n_e | 7.7e-06 | 9.4e-06 | 2.9e-05 | 8.0e-06 | 1.8e-05 |

### HELM: target region (logRho in [-2,4], logT in [3,7])

| rule | P | E | S | chiT | Cv |
|---|---|---|---|---|---|
| additive volume | 1.9e-03 | 4.1e-03 | 3.4e-04 | 1.6e-03 | 3.0e-04 |
| common n_e | 1.9e-03 | 4.4e-03 | 3.1e-04 | 1.5e-03 | 1.5e-04 |

### HELM: Skye-valid region (logT >= 6.5 OR logRho >= 3.7)

| rule | P | E | S | chiT | Cv |
|---|---|---|---|---|---|
| additive volume | 2.1e-04 | 2.6e-04 | 2.7e-04 | 1.6e-04 | 6.1e-05 |
| common n_e | 3.1e-05 | 6.6e-05 | 1.6e-04 | 4.7e-05 | 2.3e-05 |

### HELM: six worst mixtures (ranked by additive-volume target-gap P median), both rules

| mixture | region | rule | P | E | S | chiT | Cv |
|---|---|---|---|---|---|---|---|
| mix_he4_c12_o16 | target | AV | 1.9e-02 | 4.2e-02 | 6.6e-03 | 1.5e-02 | 2.0e-03 |
| mix_he4_c12_o16 | target | n_e | 2.0e-02 | 4.4e-02 | 6.4e-03 | 1.6e-02 | 1.2e-03 |
| mix_he4_c12_o16 | skye_valid | AV | 1.9e-03 | 2.6e-03 | 5.5e-03 | 1.6e-03 | 7.8e-04 |
| mix_he4_c12_o16 | skye_valid | n_e | 2.4e-04 | 4.3e-04 | 3.4e-03 | 2.7e-04 | 1.5e-04 |
| mix_he4_c12 | target | AV | 1.5e-02 | 3.3e-02 | 6.7e-03 | 1.2e-02 | 2.1e-03 |
| mix_he4_c12 | target | n_e | 1.6e-02 | 3.4e-02 | 6.4e-03 | 1.3e-02 | 1.1e-03 |
| mix_he4_c12 | skye_valid | AV | 1.2e-03 | 2.0e-03 | 5.2e-03 | 1.1e-03 | 8.1e-04 |
| mix_he4_c12 | skye_valid | n_e | 1.5e-04 | 2.9e-04 | 3.3e-03 | 1.9e-04 | 1.0e-04 |
| mix_he4_c12_10_90 | target | AV | 9.6e-03 | 2.1e-02 | 3.5e-03 | 7.2e-03 | 1.0e-03 |
| mix_he4_c12_10_90 | target | n_e | 1.0e-02 | 2.3e-02 | 3.4e-03 | 7.7e-03 | 7.6e-04 |
| mix_he4_c12_10_90 | skye_valid | AV | 8.4e-04 | 1.2e-03 | 2.7e-03 | 8.2e-04 | 4.3e-04 |
| mix_he4_c12_10_90 | skye_valid | n_e | 1.0e-04 | 2.1e-04 | 1.7e-03 | 1.3e-04 | 6.9e-05 |
| mix_c12_o16_ne20 | target | AV | 2.6e-03 | 5.7e-03 | 4.1e-04 | 2.1e-03 | 3.7e-04 |
| mix_c12_o16_ne20 | target | n_e | 2.6e-03 | 6.1e-03 | 3.7e-04 | 2.1e-03 | 1.8e-04 |
| mix_c12_o16_ne20 | skye_valid | AV | 2.9e-04 | 3.6e-04 | 3.3e-04 | 2.2e-04 | 7.7e-05 |
| mix_c12_o16_ne20 | skye_valid | n_e | 4.3e-05 | 8.8e-05 | 2.0e-04 | 4.7e-05 | 2.3e-05 |
| mix_h1_he4 | target | AV | 2.0e-03 | 4.3e-03 | 1.2e-02 | 1.6e-03 | 6.9e-04 |
| mix_h1_he4 | target | n_e | 2.2e-03 | 4.3e-03 | 1.2e-02 | 1.7e-03 | 4.9e-04 |
| mix_h1_he4 | skye_valid | AV | 8.0e-05 | 1.5e-04 | 7.4e-03 | 1.3e-04 | 7.9e-05 |
| mix_h1_he4 | skye_valid | n_e | 4.9e-06 | 3.4e-05 | 5.2e-03 | 7.4e-05 | 3.9e-05 |
| mix_c12_o16_ne22 | target | AV | 1.9e-03 | 4.1e-03 | 3.4e-04 | 1.6e-03 | 2.6e-04 |
| mix_c12_o16_ne22 | target | n_e | 1.9e-03 | 4.4e-03 | 3.1e-04 | 1.5e-03 | 1.5e-04 |
| mix_c12_o16_ne22 | skye_valid | AV | 2.1e-04 | 2.6e-04 | 2.7e-04 | 1.7e-04 | 6.1e-05 |
| mix_c12_o16_ne22 | skye_valid | n_e | 3.1e-05 | 6.6e-05 | 1.6e-04 | 1.0e-04 | 3.7e-05 |

### Composition series, both rules (target-gap median P and E)

#### C/O (minority O)

| mixture | X_min | skye P (AV) | skye P (n_e) | skye E (AV) | skye E (n_e) | helm P (AV) | helm P (n_e) |
|---|---|---|---|---|---|---|---|
| mix_c12_o16_99_01 | 0.01 | 6.6e-06 | 1.8e-06 | 1.8e-05 | 2.5e-06 | 6.4e-05 | 6.5e-05 |
| mix_c12_o16_90_10 | 0.10 | 7.0e-05 | 1.6e-05 | 1.8e-04 | 2.3e-05 | 5.9e-04 | 6.0e-04 |
| mix_c12_o16 | 0.50 | 2.7e-04 | 5.0e-05 | 6.2e-04 | 6.3e-05 | 1.8e-03 | 1.9e-03 |

#### He/C

| mixture | X_min | skye P (AV) | skye P (n_e) | skye E (AV) | skye E (n_e) | helm P (AV) | helm P (n_e) |
|---|---|---|---|---|---|---|---|
| mix_he4_c12_99_01 | 0.01 | 2.2e-05 | 1.2e-05 | 8.9e-05 | 2.6e-05 | 3.9e-04 | 4.1e-04 |
| mix_he4_c12_10_90 | 0.10 | 5.9e-04 | 3.6e-04 | 1.5e-03 | 5.5e-04 | 9.6e-03 | 1.0e-02 |
| mix_he4_c12 | 0.50 | 7.7e-04 | 4.8e-04 | 2.5e-03 | 8.9e-04 | 1.5e-02 | 1.6e-02 |

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
| skye | mix_c12_o16_90_10 | median | 7.0e-05 | 1.8e-04 | 8.6e-05 | 2.7e-04 | 1.5e-04 | 1.9e-04 |
| skye | mix_c12_o16_90_10 | p99 | 6.8e-03 | 1.0e-02 | 5.1e-03 | 1.2e-01 | 5.1e-02 | 1.4e-01 |
| skye | mix_c12_o16_90_10 | max | 3.4e-02 | 5.7e-02 | 6.7e-02 | 7.8e-01 | 7.5e-01 | 7.9e-01 |
| skye | mix_c12_o16_99_01 | median | 6.6e-06 | 1.8e-05 | 1.0e-05 | 2.9e-05 | 1.6e-05 | 2.1e-05 |
| skye | mix_c12_o16_99_01 | p99 | 6.3e-04 | 9.5e-04 | 4.9e-04 | 1.2e-02 | 4.8e-03 | 1.4e-02 |
| skye | mix_c12_o16_99_01 | max | 1.4e-03 | 2.2e-03 | 6.7e-04 | 3.3e-02 | 1.4e-01 | 9.1e-02 |
| skye | mix_he4_c12_99_01 | median | 2.2e-05 | 8.9e-05 | 4.7e-05 | 1.6e-04 | 1.9e-04 | 1.5e-04 |
| skye | mix_he4_c12_99_01 | p99 | 2.0e-03 | 3.0e-03 | 1.4e-02 | 2.6e-01 | 1.2e-02 | 2.6e-01 |
| skye | mix_he4_c12_99_01 | max | 3.3e-03 | 5.4e-03 | 2.7e-01 | 7.7e-01 | 5.0e-01 | 7.7e-01 |
| skye | mix_he4_c12_10_90 | median | 5.9e-04 | 1.5e-03 | 9.0e-04 | 2.8e-03 | 3.1e-03 | 2.6e-03 |
| skye | mix_he4_c12_10_90 | p99 | 7.6e-02 | 5.3e-02 | 8.2e-02 | 4.2e+00 | 7.5e-01 | 4.0e+00 |
| skye | mix_he4_c12_10_90 | max | 8.2e-01 | 1.9e-01 | 1.4e-01 | 3.5e+02 | 8.3e+00 | 1.0e+02 |
| skye | mix_c12_o16_ne20 | median | 3.0e-04 | 7.8e-04 | 2.9e-04 | 1.2e-03 | 6.5e-04 | 7.8e-04 |
| skye | mix_c12_o16_ne20 | p99 | 1.9e-02 | 4.4e-02 | 4.9e-02 | 3.7e-01 | 2.8e-01 | 3.1e-01 |
| skye | mix_c12_o16_ne20 | max | 2.5e-01 | 1.4e-01 | 1.1e-01 | 1.8e+00 | 1.3e+00 | 9.5e-01 |
| skye | mix_he4_c12_o16 | median | 8.6e-04 | 2.5e-03 | 1.6e-03 | 6.3e-03 | 7.3e-03 | 5.8e-03 |
| skye | mix_he4_c12_o16 | p99 | 8.8e-02 | 7.2e-02 | 1.3e-01 | 5.8e+00 | 4.5e-01 | 5.8e+00 |
| skye | mix_he4_c12_o16 | max | 2.3e-01 | 1.3e-01 | 1.4e-01 | 2.0e+01 | 6.4e-01 | 2.0e+01 |
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
| helm | mix_c12_o16_90_10 | median | 5.9e-04 | 1.3e-03 | 1.1e-04 | 5.0e-04 | 1.1e-04 | 5.2e-04 |
| helm | mix_c12_o16_90_10 | p99 | 1.9e-03 | 6.9e-03 | 2.7e-04 | 2.4e-03 | 1.5e-03 | 2.3e-02 |
| helm | mix_c12_o16_90_10 | max | 2.2e-03 | 1.0e-02 | 5.5e-04 | 2.7e-03 | 1.9e-03 | 2.3e-02 |
| helm | mix_c12_o16_99_01 | median | 6.4e-05 | 1.4e-04 | 1.1e-05 | 5.4e-05 | 1.2e-05 | 1.6e-04 |
| helm | mix_c12_o16_99_01 | p99 | 2.0e-04 | 7.3e-04 | 3.0e-05 | 2.6e-04 | 1.7e-04 | 2.3e-02 |
| helm | mix_c12_o16_99_01 | max | 2.3e-04 | 1.1e-03 | 6.0e-05 | 2.9e-04 | 2.1e-04 | 2.3e-02 |
| helm | mix_he4_c12_99_01 | median | 3.9e-04 | 8.4e-04 | 1.9e-04 | 3.2e-04 | 6.6e-05 | 2.5e-04 |
| helm | mix_he4_c12_99_01 | p99 | 1.2e-03 | 3.3e-03 | 3.8e-04 | 1.5e-03 | 8.7e-04 | 2.3e-02 |
| helm | mix_he4_c12_99_01 | max | 1.4e-03 | 3.8e-03 | 3.9e-04 | 1.8e-03 | 1.3e-03 | 2.3e-02 |
| helm | mix_he4_c12_10_90 | median | 9.6e-03 | 2.1e-02 | 3.5e-03 | 7.2e-03 | 1.0e-03 | 2.7e-03 |
| helm | mix_he4_c12_10_90 | p99 | 3.6e-02 | 1.3e-01 | 8.6e-03 | 3.7e-02 | 1.5e-02 | 3.0e-01 |
| helm | mix_he4_c12_10_90 | max | 4.8e-02 | 1.8e-01 | 8.8e-03 | 5.5e-02 | 1.6e-02 | 3.0e-01 |
| helm | mix_c12_o16_ne20 | median | 2.6e-03 | 5.7e-03 | 4.1e-04 | 2.1e-03 | 3.7e-04 | 1.4e-03 |
| helm | mix_c12_o16_ne20 | p99 | 7.2e-03 | 2.7e-02 | 9.7e-04 | 8.2e-03 | 4.1e-03 | 2.3e-02 |
| helm | mix_c12_o16_ne20 | max | 7.8e-03 | 3.3e-02 | 2.0e-03 | 9.6e-03 | 5.4e-03 | 2.3e-02 |
| helm | mix_he4_c12_o16 | median | 1.9e-02 | 4.2e-02 | 6.6e-03 | 1.5e-02 | 2.0e-03 | 5.1e-03 |
| helm | mix_he4_c12_o16 | p99 | 5.7e-02 | 1.8e-01 | 1.6e-02 | 6.1e-02 | 3.0e-02 | 6.8e-02 |
| helm | mix_he4_c12_o16 | max | 6.6e-02 | 2.3e-01 | 1.6e-02 | 7.4e-02 | 3.1e-02 | 1.3e-01 |
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
| skye_fd | mix_c12_o16_90_10 | median | 7.0e-05 | 1.8e-04 | 8.6e-05 | 2.8e-04 | 1.6e-04 | 2.3e-04 |
| skye_fd | mix_c12_o16_90_10 | p99 | 6.2e-03 | 8.9e-03 | 5.0e-03 | 3.1e-01 | 1.6e+00 | 6.5e-01 |
| skye_fd | mix_c12_o16_90_10 | max | 3.2e-02 | 5.6e-02 | 6.7e-02 | 1.3e+00 | 2.7e+00 | 1.5e+00 |
| skye_fd | mix_c12_o16_99_01 | median | 6.6e-06 | 1.8e-05 | 1.0e-05 | 3.4e-05 | 7.0e-05 | 6.4e-05 |
| skye_fd | mix_c12_o16_99_01 | p99 | 6.1e-04 | 8.4e-04 | 4.7e-04 | 3.8e-01 | 1.7e+00 | 7.4e-01 |
| skye_fd | mix_c12_o16_99_01 | max | 1.4e-03 | 1.9e-03 | 6.1e-04 | 9.7e-01 | 3.1e+00 | 1.6e+00 |
| skye_fd | mix_he4_c12_99_01 | median | 2.2e-05 | 8.9e-05 | 4.9e-05 | 1.7e-04 | 2.5e-04 | 1.8e-04 |
| skye_fd | mix_he4_c12_99_01 | p99 | 2.0e-03 | 3.0e-03 | 1.4e-02 | 7.2e-01 | 1.4e+00 | 4.7e-01 |
| skye_fd | mix_he4_c12_99_01 | max | 3.3e-03 | 5.4e-03 | 2.7e-01 | 1.4e+00 | 3.0e+00 | 8.6e-01 |
| skye_fd | mix_he4_c12_10_90 | median | 5.9e-04 | 1.5e-03 | 9.0e-04 | 3.1e-03 | 3.2e-03 | 2.7e-03 |
| skye_fd | mix_he4_c12_10_90 | p99 | 6.8e-02 | 4.1e-02 | 8.2e-02 | 4.1e+00 | 9.5e-01 | 4.1e+00 |
| skye_fd | mix_he4_c12_10_90 | max | 8.2e-01 | 1.5e-01 | 1.4e-01 | 3.5e+02 | 8.3e+00 | 1.0e+02 |
| skye_fd | mix_c12_o16_ne20 | median | 3.0e-04 | 7.8e-04 | 2.9e-04 | 1.2e-03 | 5.7e-04 | 8.2e-04 |
| skye_fd | mix_c12_o16_ne20 | p99 | 1.9e-02 | 4.1e-02 | 4.9e-02 | 4.0e-01 | 1.1e+00 | 5.1e-01 |
| skye_fd | mix_c12_o16_ne20 | max | 2.5e-01 | 1.4e-01 | 1.1e-01 | 1.6e+00 | 2.9e+00 | 5.9e+00 |
| skye_fd | mix_he4_c12_o16 | median | 8.5e-04 | 2.5e-03 | 1.7e-03 | 6.8e-03 | 7.7e-03 | 6.2e-03 |
| skye_fd | mix_he4_c12_o16 | p99 | 8.1e-02 | 6.3e-02 | 1.3e-01 | 5.9e+00 | 8.6e-01 | 5.8e+00 |
| skye_fd | mix_he4_c12_o16 | max | 2.3e-01 | 1.1e-01 | 1.4e-01 | 2.0e+01 | 1.7e+00 | 2.0e+01 |
| skye_ne | mix_c12_o16 | median | 5.0e-05 | 6.3e-05 | 3.3e-05 | 1.2e-04 | 3.0e-05 | 5.9e-05 |
| skye_ne | mix_c12_o16 | p99 | 1.5e-02 | 3.1e-02 | 4.5e-02 | 1.0e-01 | 2.4e-01 | 1.8e-01 |
| skye_ne | mix_c12_o16 | max | 2.1e-01 | 1.6e-01 | 1.0e-01 | 1.1e+00 | 1.1e+00 | 4.4e-01 |
| skye_ne | mix_o16_ne20_mg24 | median | 2.1e-04 | 5.3e-04 | 6.9e-05 | 2.2e-04 | 1.6e-05 | 6.0e-05 |
| skye_ne | mix_o16_ne20_mg24 | p99 | 1.1e-02 | 2.3e-02 | 4.5e-02 | 7.0e-02 | 2.1e-01 | 1.4e-01 |
| skye_ne | mix_o16_ne20_mg24 | max | 8.7e-02 | 1.1e-01 | 7.4e-02 | 6.5e-01 | 6.6e-01 | 2.1e-01 |
| skye_ne | mix_he4_c12 | median | 4.8e-04 | 8.9e-04 | 3.8e-04 | 3.9e-04 | 7.7e-04 | 1.6e-04 |
| skye_ne | mix_he4_c12 | p99 | 3.9e-02 | 3.9e-02 | 2.0e-01 | 2.2e-01 | 4.2e-01 | 2.6e-01 |
| skye_ne | mix_he4_c12 | max | 9.4e-02 | 1.3e-01 | 2.5e-01 | 3.3e-01 | 4.4e-01 | 3.2e-01 |
| skye_ne | mix_h1_he4 | median | 7.3e-05 | 1.6e-04 | 1.4e-04 | 2.2e-04 | 4.5e-04 | 9.7e-05 |
| skye_ne | mix_h1_he4 | p99 | 1.2e-02 | 1.2e-02 | 4.1e-02 | 7.7e-01 | 3.5e-01 | 8.1e-01 |
| skye_ne | mix_h1_he4 | max | 4.6e-02 | 3.6e-02 | 2.5e-01 | 7.8e+00 | 2.1e+00 | 6.1e+00 |
| skye_ne | mix_c12_o16_ne22 | median | 5.1e-05 | 7.3e-05 | 3.9e-05 | 1.7e-04 | 3.3e-05 | 6.2e-05 |
| skye_ne | mix_c12_o16_ne22 | p99 | 1.3e-02 | 2.7e-02 | 4.5e-02 | 9.4e-02 | 2.4e-01 | 1.9e-01 |
| skye_ne | mix_c12_o16_ne22 | max | 7.6e-02 | 1.1e-01 | 6.6e-02 | 6.9e-01 | 8.7e-01 | 2.2e-01 |
| skye_ne | mix_c12_o16_90_10 | median | 1.6e-05 | 2.3e-05 | 1.7e-05 | 3.9e-05 | 1.0e-05 | 1.6e-05 |
| skye_ne | mix_c12_o16_90_10 | p99 | 3.2e-03 | 7.7e-03 | 5.1e-03 | 2.1e-02 | 4.1e-02 | 3.1e-02 |
| skye_ne | mix_c12_o16_90_10 | max | 3.7e-02 | 5.8e-02 | 6.7e-02 | 8.6e-02 | 3.7e-01 | 3.0e-01 |
| skye_ne | mix_c12_o16_99_01 | median | 1.8e-06 | 2.5e-06 | 2.9e-06 | 4.0e-06 | 1.2e-06 | 1.6e-06 |
| skye_ne | mix_c12_o16_99_01 | p99 | 2.8e-04 | 6.0e-04 | 5.2e-04 | 1.9e-03 | 4.0e-03 | 3.0e-03 |
| skye_ne | mix_c12_o16_99_01 | max | 2.3e-03 | 2.8e-03 | 1.3e-03 | 6.0e-03 | 6.7e-03 | 3.8e-03 |
| skye_ne | mix_he4_c12_99_01 | median | 1.2e-05 | 2.6e-05 | 1.3e-05 | 9.1e-06 | 2.9e-05 | 6.8e-06 |
| skye_ne | mix_he4_c12_99_01 | p99 | 5.8e-04 | 1.0e-03 | 1.3e-03 | 3.5e-03 | 3.5e-03 | 3.4e-03 |
| skye_ne | mix_he4_c12_99_01 | max | 1.1e-03 | 1.8e-02 | 2.7e-01 | 2.1e-01 | 6.9e-01 | 4.0e-01 |
| skye_ne | mix_he4_c12_10_90 | median | 3.6e-04 | 5.5e-04 | 2.8e-04 | 4.1e-04 | 2.7e-04 | 1.7e-04 |
| skye_ne | mix_he4_c12_10_90 | p99 | 2.5e-02 | 3.4e-02 | 7.7e-02 | 4.5e-01 | 7.2e-01 | 4.5e-01 |
| skye_ne | mix_he4_c12_10_90 | max | 2.0e-01 | 2.3e-01 | 1.4e-01 | 7.5e-01 | 7.8e-01 | 5.4e-01 |
| skye_ne | mix_c12_o16_ne20 | median | 7.2e-05 | 2.6e-04 | 5.3e-05 | 1.9e-04 | 3.7e-05 | 6.8e-05 |
| skye_ne | mix_c12_o16_ne20 | p99 | 1.5e-02 | 2.7e-02 | 4.5e-02 | 7.7e-02 | 2.6e-01 | 1.9e-01 |
| skye_ne | mix_c12_o16_ne20 | max | 8.5e-02 | 1.2e-01 | 6.1e-02 | 6.5e-01 | 8.1e-01 | 2.2e-01 |
| skye_ne | mix_he4_c12_o16 | median | 5.5e-04 | 8.9e-04 | 2.6e-04 | 4.4e-04 | 6.8e-04 | 1.8e-04 |
| skye_ne | mix_he4_c12_o16 | p99 | 4.1e-02 | 5.2e-02 | 1.2e-01 | 2.9e-01 | 4.5e-01 | 3.5e-01 |
| skye_ne | mix_he4_c12_o16 | max | 9.4e-02 | 1.5e-01 | 1.4e-01 | 5.4e-01 | 5.6e-01 | 4.6e-01 |
| helm_ne | mix_c12_o16 | median | 1.9e-03 | 4.4e-03 | 2.8e-04 | 1.5e-03 | 1.5e-04 | 8.8e-04 |
| helm_ne | mix_c12_o16 | p99 | 7.2e-03 | 3.6e-02 | 1.6e-03 | 6.8e-03 | 2.0e-03 | 3.6e-01 |
| helm_ne | mix_c12_o16 | max | 7.8e-03 | 4.7e-02 | 2.2e-03 | 7.3e-03 | 2.0e-03 | 3.7e-01 |
| helm_ne | mix_o16_ne20_mg24 | median | 1.9e-03 | 4.3e-03 | 1.9e-04 | 1.5e-03 | 1.1e-04 | 9.8e-04 |
| helm_ne | mix_o16_ne20_mg24 | p99 | 6.0e-03 | 2.5e-02 | 1.4e-03 | 5.3e-03 | 1.8e-03 | 3.6e-01 |
| helm_ne | mix_o16_ne20_mg24 | max | 6.3e-03 | 2.7e-02 | 1.8e-03 | 5.9e-03 | 1.8e-03 | 3.7e-01 |
| helm_ne | mix_he4_c12 | median | 1.6e-02 | 3.4e-02 | 6.4e-03 | 1.3e-02 | 1.1e-03 | 4.2e-03 |
| helm_ne | mix_he4_c12 | p99 | 6.3e-02 | 1.9e-01 | 1.4e-02 | 5.6e-02 | 2.4e-02 | 3.5e-01 |
| helm_ne | mix_he4_c12 | max | 6.7e-02 | 2.1e-01 | 1.4e-02 | 6.4e-02 | 2.5e-02 | 3.7e-01 |
| helm_ne | mix_h1_he4 | median | 2.2e-03 | 4.3e-03 | 1.2e-02 | 1.7e-03 | 4.9e-04 | 7.2e-04 |
| helm_ne | mix_h1_he4 | p99 | 2.5e-02 | 4.3e-02 | 6.1e-02 | 1.4e-02 | 4.3e-03 | 3.3e-01 |
| helm_ne | mix_h1_he4 | max | 6.2e-02 | 5.2e-02 | 7.3e-02 | 7.6e-02 | 4.6e-03 | 3.6e-01 |
| helm_ne | mix_c12_o16_ne22 | median | 1.9e-03 | 4.4e-03 | 3.1e-04 | 1.5e-03 | 1.5e-04 | 8.9e-04 |
| helm_ne | mix_c12_o16_ne22 | p99 | 6.6e-03 | 2.7e-02 | 1.1e-03 | 5.7e-03 | 2.1e-03 | 3.6e-01 |
| helm_ne | mix_c12_o16_ne22 | max | 6.8e-03 | 2.9e-02 | 1.6e-03 | 6.3e-03 | 2.1e-03 | 3.7e-01 |
| helm_ne | mix_c12_o16_90_10 | median | 6.0e-04 | 1.4e-03 | 9.3e-05 | 4.9e-04 | 4.9e-05 | 3.5e-04 |
| helm_ne | mix_c12_o16_90_10 | p99 | 2.2e-03 | 9.5e-03 | 4.7e-04 | 2.1e-03 | 6.2e-04 | 3.6e-01 |
| helm_ne | mix_c12_o16_90_10 | max | 2.4e-03 | 1.2e-02 | 6.6e-04 | 2.2e-03 | 6.3e-04 | 3.7e-01 |
| helm_ne | mix_c12_o16_99_01 | median | 6.5e-05 | 1.5e-04 | 9.9e-06 | 5.2e-05 | 5.3e-06 | 1.3e-04 |
| helm_ne | mix_c12_o16_99_01 | p99 | 2.4e-04 | 9.8e-04 | 5.2e-05 | 2.2e-04 | 6.7e-05 | 3.6e-01 |
| helm_ne | mix_c12_o16_99_01 | max | 2.6e-04 | 1.2e-03 | 7.3e-05 | 2.4e-04 | 6.8e-05 | 3.7e-01 |
| helm_ne | mix_he4_c12_99_01 | median | 4.1e-04 | 8.5e-04 | 1.8e-04 | 3.2e-04 | 3.7e-05 | 1.8e-04 |
| helm_ne | mix_he4_c12_99_01 | p99 | 1.6e-03 | 4.1e-03 | 3.8e-04 | 1.3e-03 | 5.9e-04 | 3.5e-01 |
| helm_ne | mix_he4_c12_99_01 | max | 1.7e-03 | 4.4e-03 | 3.9e-04 | 1.5e-03 | 6.0e-04 | 3.7e-01 |
| helm_ne | mix_he4_c12_10_90 | median | 1.0e-02 | 2.3e-02 | 3.4e-03 | 7.7e-03 | 7.6e-04 | 2.8e-03 |
| helm_ne | mix_he4_c12_10_90 | p99 | 4.3e-02 | 2.1e-01 | 8.5e-03 | 3.8e-02 | 1.5e-02 | 3.6e-01 |
| helm_ne | mix_he4_c12_10_90 | max | 4.7e-02 | 2.5e-01 | 8.8e-03 | 4.4e-02 | 1.6e-02 | 3.7e-01 |
| helm_ne | mix_c12_o16_ne20 | median | 2.6e-03 | 6.1e-03 | 3.7e-04 | 2.1e-03 | 1.8e-04 | 1.1e-03 |
| helm_ne | mix_c12_o16_ne20 | p99 | 8.9e-03 | 3.6e-02 | 1.6e-03 | 8.0e-03 | 2.8e-03 | 3.6e-01 |
| helm_ne | mix_c12_o16_ne20 | max | 9.4e-03 | 4.1e-02 | 2.2e-03 | 8.8e-03 | 2.8e-03 | 3.7e-01 |
| helm_ne | mix_he4_c12_o16 | median | 2.0e-02 | 4.4e-02 | 6.4e-03 | 1.6e-02 | 1.2e-03 | 5.0e-03 |
| helm_ne | mix_he4_c12_o16 | p99 | 7.0e-02 | 2.3e-01 | 1.6e-02 | 6.6e-02 | 3.0e-02 | 3.6e-01 |
| helm_ne | mix_he4_c12_o16 | max | 7.5e-02 | 2.7e-01 | 1.6e-02 | 7.2e-02 | 3.1e-02 | 3.7e-01 |

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
| skye | mix_c12_o16_90_10 | median | 4.3e-06 | 6.6e-06 | 2.6e-05 | 3.6e-05 | 4.4e-05 | 1.5e-05 |
| skye | mix_c12_o16_90_10 | p99 | 8.7e-05 | 2.0e-04 | 1.3e-02 | 5.9e-01 | 4.4e-01 | 6.1e-01 |
| skye | mix_c12_o16_90_10 | max | 5.9e-04 | 1.6e-03 | 2.4e-01 | 1.2e+00 | 9.0e-01 | 1.2e+00 |
| skye | mix_c12_o16_99_01 | median | 4.7e-07 | 6.4e-07 | 3.8e-06 | 3.8e-06 | 4.9e-06 | 1.7e-06 |
| skye | mix_c12_o16_99_01 | p99 | 8.4e-06 | 2.2e-05 | 1.2e-03 | 6.4e-02 | 4.4e-01 | 3.3e-01 |
| skye | mix_c12_o16_99_01 | max | 2.6e-05 | 7.3e-05 | 3.0e-03 | 1.2e-01 | 9.0e-01 | 8.6e-01 |
| skye | mix_he4_c12_99_01 | median | 1.0e-06 | 1.6e-06 | 8.2e-06 | 3.1e-05 | 3.6e-05 | 1.9e-05 |
| skye | mix_he4_c12_99_01 | p99 | 9.3e-05 | 1.4e-04 | 7.2e-03 | 3.6e-01 | 5.0e-02 | 3.5e-01 |
| skye | mix_he4_c12_99_01 | max | 2.8e-04 | 6.2e-04 | 2.0e+00 | 2.0e+00 | 8.7e-01 | 2.0e+00 |
| skye | mix_he4_c12_10_90 | median | 4.7e-05 | 5.6e-05 | 1.5e-04 | 5.7e-04 | 5.6e-04 | 3.0e-04 |
| skye | mix_he4_c12_10_90 | p99 | 2.6e-03 | 1.9e-03 | 1.9e-01 | 4.2e+00 | 5.2e-01 | 4.2e+00 |
| skye | mix_he4_c12_10_90 | max | 1.3e-02 | 4.2e-03 | 2.8e-01 | 1.9e+01 | 8.1e-01 | 1.9e+01 |
| skye | mix_c12_o16_ne20 | median | 2.0e-05 | 3.1e-05 | 5.8e-05 | 1.6e-04 | 1.9e-04 | 6.4e-05 |
| skye | mix_c12_o16_ne20 | p99 | 4.0e-04 | 9.8e-04 | 7.1e-02 | 2.4e+00 | 4.4e-01 | 2.3e+00 |
| skye | mix_c12_o16_ne20 | max | 2.2e-03 | 3.2e-03 | 1.5e-01 | 4.7e+00 | 8.9e-01 | 4.7e+00 |
| skye | mix_he4_c12_o16 | median | 5.1e-05 | 6.3e-05 | 3.1e-04 | 1.1e-03 | 1.1e-03 | 6.2e-04 |
| skye | mix_he4_c12_o16 | p99 | 5.1e-03 | 4.8e-03 | 1.5e-01 | 9.6e+00 | 3.6e-01 | 9.6e+00 |
| skye | mix_he4_c12_o16 | max | 2.2e-02 | 1.2e-02 | 2.5e-01 | 4.4e+01 | 8.1e-01 | 4.4e+01 |
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
| helm | mix_c12_o16_90_10 | median | 6.1e-05 | 7.8e-05 | 7.9e-05 | 5.1e-05 | 2.1e-05 | 1.7e-04 |
| helm | mix_c12_o16_90_10 | p99 | 9.2e-04 | 2.0e-03 | 5.1e-04 | 4.6e-03 | 6.6e-03 | 2.3e-02 |
| helm | mix_c12_o16_90_10 | max | 1.7e-03 | 3.3e-03 | 6.7e-04 | 7.6e-02 | 9.3e-02 | 3.2e-02 |
| helm | mix_c12_o16_99_01 | median | 6.6e-06 | 8.6e-06 | 8.7e-06 | 5.6e-06 | 2.4e-06 | 5.8e-05 |
| helm | mix_c12_o16_99_01 | p99 | 9.8e-05 | 2.3e-04 | 5.6e-05 | 6.1e-04 | 7.2e-04 | 2.3e-02 |
| helm | mix_c12_o16_99_01 | max | 1.5e-04 | 5.6e-04 | 7.3e-05 | 9.9e-03 | 1.2e-02 | 3.2e-02 |
| helm | mix_he4_c12_99_01 | median | 2.9e-05 | 5.4e-05 | 1.5e-04 | 3.0e-05 | 2.5e-05 | 1.7e-04 |
| helm | mix_he4_c12_99_01 | p99 | 4.5e-04 | 1.0e-03 | 1.2e-03 | 3.7e-03 | 6.8e-03 | 2.5e-02 |
| helm | mix_he4_c12_99_01 | max | 6.1e-04 | 1.4e-03 | 2.0e-03 | 6.2e-02 | 1.1e-01 | 3.2e-02 |
| helm | mix_he4_c12_10_90 | median | 8.4e-04 | 1.2e-03 | 2.7e-03 | 8.2e-04 | 4.3e-04 | 3.2e-03 |
| helm | mix_he4_c12_10_90 | p99 | 1.2e-02 | 2.7e-02 | 2.0e-02 | 1.9e-02 | 5.3e-02 | 1.1e-01 |
| helm | mix_he4_c12_10_90 | max | 1.6e-02 | 4.3e-02 | 2.8e-02 | 3.1e-01 | 9.2e-01 | 3.3e-01 |
| helm | mix_c12_o16_ne20 | median | 2.9e-04 | 3.6e-04 | 3.3e-04 | 2.2e-04 | 7.7e-05 | 7.3e-04 |
| helm | mix_c12_o16_ne20 | p99 | 3.8e-03 | 9.9e-03 | 1.9e-03 | 3.8e-03 | 2.8e-03 | 2.3e-02 |
| helm | mix_c12_o16_ne20 | max | 5.4e-03 | 1.7e-02 | 2.1e-03 | 5.3e-03 | 4.8e-02 | 1.0e-01 |
| helm | mix_he4_c12_o16 | median | 1.9e-03 | 2.6e-03 | 5.5e-03 | 1.6e-03 | 7.8e-04 | 5.7e-03 |
| helm | mix_he4_c12_o16 | p99 | 2.5e-02 | 6.0e-02 | 3.5e-02 | 2.7e-02 | 3.0e-02 | 6.6e-02 |
| helm | mix_he4_c12_o16 | max | 3.6e-02 | 9.3e-02 | 4.5e-02 | 2.5e-01 | 2.2e-01 | 1.3e-01 |
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
| skye_fd | mix_c12_o16_90_10 | median | 4.3e-06 | 6.7e-06 | 2.8e-05 | 6.4e-05 | 1.3e-04 | 8.3e-05 |
| skye_fd | mix_c12_o16_90_10 | p99 | 8.1e-05 | 2.1e-04 | 1.3e-02 | 1.0e+00 | 1.1e+00 | 1.1e+00 |
| skye_fd | mix_c12_o16_90_10 | max | 5.8e-04 | 1.6e-03 | 2.4e-01 | 1.7e+00 | 3.0e+00 | 5.7e+01 |
| skye_fd | mix_c12_o16_99_01 | median | 4.7e-07 | 6.6e-07 | 4.3e-06 | 2.0e-05 | 7.5e-05 | 5.6e-05 |
| skye_fd | mix_c12_o16_99_01 | p99 | 8.2e-06 | 2.2e-05 | 1.2e-03 | 9.9e-01 | 1.2e+00 | 9.9e-01 |
| skye_fd | mix_c12_o16_99_01 | max | 2.6e-05 | 7.3e-05 | 3.0e-03 | 1.9e+00 | 3.5e+00 | 4.8e+01 |
| skye_fd | mix_he4_c12_99_01 | median | 1.0e-06 | 1.8e-06 | 9.1e-06 | 5.4e-05 | 1.3e-04 | 9.1e-05 |
| skye_fd | mix_he4_c12_99_01 | p99 | 9.0e-05 | 1.4e-04 | 7.4e-03 | 1.2e+00 | 7.3e-01 | 6.4e-01 |
| skye_fd | mix_he4_c12_99_01 | max | 2.8e-04 | 6.2e-04 | 2.0e+00 | 4.3e+01 | 2.8e+01 | 1.1e+01 |
| skye_fd | mix_he4_c12_10_90 | median | 4.7e-05 | 5.8e-05 | 1.7e-04 | 6.0e-04 | 4.6e-04 | 4.7e-04 |
| skye_fd | mix_he4_c12_10_90 | p99 | 2.6e-03 | 1.9e-03 | 1.9e-01 | 4.3e+00 | 8.8e-01 | 4.3e+00 |
| skye_fd | mix_he4_c12_10_90 | max | 1.2e-02 | 4.2e-03 | 2.8e-01 | 1.9e+01 | 2.2e+00 | 1.9e+01 |
| skye_fd | mix_c12_o16_ne20 | median | 2.0e-05 | 3.1e-05 | 5.9e-05 | 2.0e-04 | 2.4e-04 | 1.6e-04 |
| skye_fd | mix_c12_o16_ne20 | p99 | 3.7e-04 | 9.8e-04 | 7.3e-02 | 2.6e+00 | 1.1e+00 | 3.0e+00 |
| skye_fd | mix_c12_o16_ne20 | max | 1.2e-03 | 3.2e-03 | 1.5e-01 | 4.9e+00 | 1.5e+00 | 6.8e+01 |
| skye_fd | mix_he4_c12_o16 | median | 5.2e-05 | 6.9e-05 | 3.6e-04 | 1.3e-03 | 1.0e-03 | 9.0e-04 |
| skye_fd | mix_he4_c12_o16 | p99 | 5.0e-03 | 4.8e-03 | 1.5e-01 | 9.8e+00 | 7.7e-01 | 9.7e+00 |
| skye_fd | mix_he4_c12_o16 | max | 2.0e-02 | 1.2e-02 | 2.5e-01 | 4.4e+01 | 2.7e+00 | 4.4e+01 |
| skye_ne | mix_c12_o16 | median | 7.7e-06 | 9.4e-06 | 2.9e-05 | 8.0e-06 | 1.8e-05 | 5.0e-06 |
| skye_ne | mix_c12_o16 | p99 | 2.2e-04 | 4.9e-04 | 6.8e-02 | 1.1e-01 | 4.2e-01 | 3.1e-01 |
| skye_ne | mix_c12_o16 | max | 6.1e-04 | 1.6e-03 | 3.2e-01 | 2.0e-01 | 9.0e-01 | 9.2e-01 |
| skye_ne | mix_o16_ne20_mg24 | median | 9.6e-06 | 1.3e-05 | 6.1e-05 | 1.2e-05 | 1.1e-05 | 5.4e-06 |
| skye_ne | mix_o16_ne20_mg24 | p99 | 2.3e-04 | 5.3e-04 | 3.7e-02 | 4.6e-02 | 4.7e-01 | 3.5e-01 |
| skye_ne | mix_o16_ne20_mg24 | max | 5.5e-04 | 1.7e-03 | 1.6e-01 | 1.6e-01 | 9.2e-01 | 9.1e-01 |
| skye_ne | mix_he4_c12 | median | 1.5e-05 | 1.9e-05 | 3.5e-05 | 1.3e-04 | 1.4e-04 | 1.8e-05 |
| skye_ne | mix_he4_c12 | p99 | 5.6e-03 | 1.6e-03 | 2.5e-01 | 3.3e-01 | 4.5e-01 | 3.4e-01 |
| skye_ne | mix_he4_c12 | max | 3.8e-02 | 2.4e-03 | 6.0e-01 | 7.8e-01 | 9.0e-01 | 9.4e-01 |
| skye_ne | mix_h1_he4 | median | 1.3e-06 | 1.1e-06 | 1.7e-04 | 2.2e-05 | 8.9e-06 | 9.5e-06 |
| skye_ne | mix_h1_he4 | p99 | 6.1e-02 | 4.1e-03 | 3.7e-02 | 1.2e+00 | 8.5e-01 | 1.7e+00 |
| skye_ne | mix_h1_he4 | max | 2.8e-01 | 1.8e-02 | 4.2e-02 | 4.6e+02 | 9.6e+02 | 1.1e+02 |
| skye_ne | mix_c12_o16_ne22 | median | 8.7e-06 | 1.1e-05 | 3.3e-05 | 8.9e-06 | 2.1e-05 | 6.3e-06 |
| skye_ne | mix_c12_o16_ne22 | p99 | 2.2e-04 | 5.0e-04 | 6.8e-02 | 2.3e-01 | 3.4e-01 | 2.7e-01 |
| skye_ne | mix_c12_o16_ne22 | max | 6.2e-04 | 1.6e-03 | 3.2e-01 | 2.5e-01 | 6.0e-01 | 1.1e+00 |
| skye_ne | mix_c12_o16_90_10 | median | 2.8e-06 | 3.4e-06 | 1.2e-05 | 2.9e-06 | 6.4e-06 | 1.8e-06 |
| skye_ne | mix_c12_o16_90_10 | p99 | 4.7e-05 | 1.1e-04 | 1.2e-02 | 2.7e-02 | 4.3e-01 | 3.2e-01 |
| skye_ne | mix_c12_o16_90_10 | max | 5.9e-04 | 1.6e-03 | 2.4e-01 | 2.9e-01 | 9.0e-01 | 9.2e-01 |
| skye_ne | mix_c12_o16_99_01 | median | 3.0e-07 | 3.7e-07 | 1.8e-06 | 3.2e-07 | 6.8e-07 | 2.4e-07 |
| skye_ne | mix_c12_o16_99_01 | p99 | 4.4e-06 | 1.1e-05 | 1.1e-03 | 2.0e-02 | 4.2e-01 | 3.1e-01 |
| skye_ne | mix_c12_o16_99_01 | max | 8.1e-06 | 2.9e-05 | 3.0e-03 | 4.7e-02 | 9.0e-01 | 9.2e-01 |
| skye_ne | mix_he4_c12_99_01 | median | 3.1e-07 | 2.1e-07 | 1.7e-06 | 4.8e-06 | 5.2e-06 | 8.0e-07 |
| skye_ne | mix_he4_c12_99_01 | p99 | 5.5e-05 | 4.7e-05 | 1.6e-03 | 2.2e-02 | 4.4e-01 | 3.1e-01 |
| skye_ne | mix_he4_c12_99_01 | max | 2.8e-04 | 6.3e-04 | 2.0e+00 | 9.1e-01 | 9.1e-01 | 9.5e-01 |
| skye_ne | mix_he4_c12_10_90 | median | 4.2e-05 | 4.8e-05 | 2.3e-05 | 5.1e-05 | 8.4e-05 | 1.4e-05 |
| skye_ne | mix_he4_c12_10_90 | p99 | 2.1e-03 | 8.5e-04 | 1.7e-01 | 2.6e-01 | 5.7e-01 | 3.2e-01 |
| skye_ne | mix_he4_c12_10_90 | max | 7.9e-03 | 2.3e-03 | 2.8e-01 | 4.1e-01 | 9.0e-01 | 9.3e-01 |
| skye_ne | mix_c12_o16_ne20 | median | 1.3e-05 | 1.6e-05 | 4.6e-05 | 1.2e-05 | 2.5e-05 | 6.8e-06 |
| skye_ne | mix_c12_o16_ne20 | p99 | 2.4e-04 | 5.4e-04 | 6.6e-02 | 9.2e-02 | 4.2e-01 | 3.1e-01 |
| skye_ne | mix_c12_o16_ne20 | max | 4.1e-04 | 1.6e-03 | 1.5e-01 | 2.1e-01 | 8.9e-01 | 9.2e-01 |
| skye_ne | mix_he4_c12_o16 | median | 2.9e-05 | 3.7e-05 | 6.2e-05 | 1.2e-04 | 1.8e-04 | 1.9e-05 |
| skye_ne | mix_he4_c12_o16 | p99 | 4.1e-03 | 1.3e-03 | 1.4e-01 | 2.6e-01 | 4.4e-01 | 3.1e-01 |
| skye_ne | mix_he4_c12_o16 | max | 1.6e-02 | 2.3e-03 | 2.5e-01 | 5.7e-01 | 9.0e-01 | 9.2e-01 |
| helm_ne | mix_c12_o16 | median | 3.0e-05 | 6.0e-05 | 1.4e-04 | 3.0e-05 | 1.8e-05 | 2.7e-01 |
| helm_ne | mix_c12_o16 | p99 | 2.3e-03 | 5.5e-03 | 1.4e-03 | 2.2e-03 | 1.9e-03 | 4.3e+00 |
| helm_ne | mix_c12_o16 | max | 3.7e-03 | 1.0e-02 | 1.7e-03 | 3.6e-03 | 2.0e-03 | 6.8e+00 |
| helm_ne | mix_o16_ne20_mg24 | median | 3.7e-05 | 7.5e-05 | 1.0e-04 | 3.7e-05 | 1.6e-05 | 2.7e-01 |
| helm_ne | mix_o16_ne20_mg24 | p99 | 2.7e-03 | 6.9e-03 | 1.0e-03 | 2.6e-03 | 1.7e-03 | 4.3e+00 |
| helm_ne | mix_o16_ne20_mg24 | max | 4.5e-03 | 1.4e-02 | 1.6e-03 | 4.5e-03 | 1.8e-03 | 6.8e+00 |
| helm_ne | mix_he4_c12 | median | 1.5e-04 | 2.9e-04 | 3.3e-03 | 1.9e-04 | 1.0e-04 | 2.6e-01 |
| helm_ne | mix_he4_c12 | p99 | 1.7e-02 | 3.7e-02 | 3.4e-02 | 1.6e-02 | 2.2e-02 | 4.3e+00 |
| helm_ne | mix_he4_c12 | max | 2.6e-02 | 5.9e-02 | 5.0e-02 | 2.6e-02 | 2.4e-02 | 6.7e+00 |
| helm_ne | mix_h1_he4 | median | 4.9e-06 | 3.4e-05 | 5.2e-03 | 7.4e-05 | 3.9e-05 | 2.6e-01 |
| helm_ne | mix_h1_he4 | p99 | 1.7e-03 | 3.2e-03 | 7.3e-02 | 2.0e-02 | 1.2e-02 | 4.3e+00 |
| helm_ne | mix_h1_he4 | max | 1.3e-02 | 1.2e-02 | 1.0e-01 | 1.5e-01 | 3.4e-01 | 5.8e+00 |
| helm_ne | mix_c12_o16_ne22 | median | 3.1e-05 | 6.6e-05 | 1.6e-04 | 1.0e-04 | 3.7e-05 | 2.7e-01 |
| helm_ne | mix_c12_o16_ne22 | p99 | 2.5e-03 | 6.1e-03 | 1.6e-03 | 6.8e-03 | 4.9e-03 | 4.3e+00 |
| helm_ne | mix_c12_o16_ne22 | max | 4.1e-03 | 1.1e-02 | 2.0e-03 | 6.5e-02 | 1.2e-01 | 6.8e+00 |
| helm_ne | mix_c12_o16_90_10 | median | 9.5e-06 | 2.1e-05 | 4.9e-05 | 9.6e-06 | 6.6e-06 | 2.7e-01 |
| helm_ne | mix_c12_o16_90_10 | p99 | 7.3e-04 | 1.7e-03 | 5.0e-04 | 6.9e-04 | 6.0e-04 | 4.3e+00 |
| helm_ne | mix_c12_o16_90_10 | max | 1.2e-03 | 3.1e-03 | 6.6e-04 | 1.1e-03 | 1.1e-03 | 6.8e+00 |
| helm_ne | mix_c12_o16_99_01 | median | 1.0e-06 | 2.8e-06 | 5.4e-06 | 1.0e-06 | 8.4e-07 | 2.7e-01 |
| helm_ne | mix_c12_o16_99_01 | p99 | 7.8e-05 | 1.9e-04 | 5.5e-05 | 7.4e-05 | 6.5e-05 | 4.3e+00 |
| helm_ne | mix_c12_o16_99_01 | max | 1.3e-04 | 5.6e-04 | 7.4e-05 | 1.2e-04 | 1.1e-03 | 6.8e+00 |
| helm_ne | mix_he4_c12_99_01 | median | 3.3e-06 | 8.5e-06 | 9.7e-05 | 4.2e-06 | 3.2e-06 | 2.6e-01 |
| helm_ne | mix_he4_c12_99_01 | p99 | 4.3e-04 | 9.0e-04 | 1.1e-03 | 4.0e-04 | 5.3e-04 | 4.3e+00 |
| helm_ne | mix_he4_c12_99_01 | max | 6.5e-04 | 1.4e-03 | 2.0e-03 | 6.2e-04 | 8.3e-04 | 6.7e+00 |
| helm_ne | mix_he4_c12_10_90 | median | 1.0e-04 | 2.1e-04 | 1.7e-03 | 1.3e-04 | 6.9e-05 | 2.6e-01 |
| helm_ne | mix_he4_c12_10_90 | p99 | 1.1e-02 | 2.4e-02 | 2.0e-02 | 1.0e-02 | 1.4e-02 | 4.3e+00 |
| helm_ne | mix_he4_c12_10_90 | max | 1.7e-02 | 4.2e-02 | 2.8e-02 | 1.7e-02 | 1.5e-02 | 6.8e+00 |
| helm_ne | mix_c12_o16_ne20 | median | 4.3e-05 | 8.8e-05 | 2.0e-04 | 4.7e-05 | 2.3e-05 | 2.7e-01 |
| helm_ne | mix_c12_o16_ne20 | p99 | 3.5e-03 | 8.4e-03 | 1.9e-03 | 3.3e-03 | 2.7e-03 | 4.3e+00 |
| helm_ne | mix_c12_o16_ne20 | max | 5.7e-03 | 1.6e-02 | 2.1e-03 | 5.5e-03 | 2.8e-03 | 6.8e+00 |
| helm_ne | mix_he4_c12_o16 | median | 2.4e-04 | 4.3e-04 | 3.4e-03 | 2.7e-04 | 1.5e-04 | 2.7e-01 |
| helm_ne | mix_he4_c12_o16 | p99 | 2.4e-02 | 5.4e-02 | 3.4e-02 | 2.3e-02 | 2.9e-02 | 4.3e+00 |
| helm_ne | mix_he4_c12_o16 | max | 3.8e-02 | 9.1e-02 | 4.5e-02 | 3.7e-02 | 3.1e-02 | 6.8e+00 |

## Composition trend (target gap): error vs minority mass fraction

Target-gap median relative error in P and E for the C/O and He/C series, ordered by the minority (smaller) mass fraction X_min. Additive volume is exact for a pure species, and the electron-fragmentation error is quadratic in the electron-density split, so for small X_min the error is expected to fall roughly in proportion to X_min. Compare P/E against `validation/fragmentation_estimate.py`.

### C/O (minority O)

| mixture | X_min | skye P | skye E | helm P | helm E |
|---|---|---|---|---|---|
| mix_c12_o16_99_01 | 0.01 | 6.6e-06 | 1.8e-05 | 6.4e-05 | 1.4e-04 |
| mix_c12_o16_90_10 | 0.10 | 7.0e-05 | 1.8e-04 | 5.9e-04 | 1.3e-03 |
| mix_c12_o16 | 0.50 | 2.7e-04 | 6.2e-04 | 1.8e-03 | 4.1e-03 |

### He/C (minority = smaller of He,C)

| mixture | X_min | skye P | skye E | helm P | helm E |
|---|---|---|---|---|---|
| mix_he4_c12_99_01 | 0.01 | 2.2e-05 | 8.9e-05 | 3.9e-04 | 8.4e-04 |
| mix_he4_c12_10_90 | 0.10 | 5.9e-04 | 1.5e-03 | 9.6e-03 | 2.1e-02 |
| mix_he4_c12 | 0.50 | 7.7e-04 | 2.5e-03 | 1.5e-02 | 3.3e-02 |

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
| skye | mix_c12_o16_90_10 | S | 1.832e+06 | 2.6e-05 | 1.7e-02 | 8.6e-05 | 1.1e-02 |
| skye | mix_c12_o16_90_10 | E | 2.255e+13 | 6.6e-06 | 1.2e-04 | 1.8e-04 | 6.4e-02 |
| skye | mix_c12_o16_99_01 | S | 3.060e+05 | 3.8e-06 | 2.8e-03 | 1.0e-05 | 1.9e-03 |
| skye | mix_c12_o16_99_01 | E | 2.425e+12 | 6.4e-07 | 1.3e-05 | 1.8e-05 | 7.4e-03 |
| skye | mix_he4_c12_99_01 | S | 4.639e+05 | 8.2e-06 | 1.3e-03 | 4.7e-05 | 1.4e-03 |
| skye | mix_he4_c12_99_01 | E | 6.095e+12 | 1.6e-06 | 3.7e-05 | 8.9e-05 | 1.7e-02 |
| skye | mix_he4_c12_10_90 | S | 4.674e+06 | 1.5e-04 | 2.5e-02 | 9.0e-04 | 2.1e-02 |
| skye | mix_he4_c12_10_90 | E | 1.375e+14 | 5.6e-05 | 8.1e-04 | 1.5e-03 | 3.2e-01 |
| skye | mix_c12_o16_ne20 | S | 5.229e+06 | 5.8e-05 | 4.9e-02 | 2.9e-04 | 2.8e-02 |
| skye | mix_c12_o16_ne20 | E | 1.205e+14 | 3.1e-05 | 6.6e-04 | 7.8e-04 | 2.5e-01 |
| skye | mix_he4_c12_o16 | S | 9.528e+06 | 3.1e-04 | 4.5e-02 | 1.6e-03 | 3.7e-02 |
| skye | mix_he4_c12_o16 | E | 4.052e+14 | 6.3e-05 | 2.4e-03 | 2.5e-03 | 5.4e-01 |
| helm | mix_c12_o16 | S | 4.141e+06 | 1.1e-02 | 2.4e-04 | 1.2e-02 | 3.1e-04 |
| helm | mix_o16_ne20_mg24 | S | 3.992e+06 | 1.1e-02 | 1.7e-04 | 1.1e-02 | 2.1e-04 |
| helm | mix_he4_c12 | S | 7.791e+06 | 9.7e-03 | 5.2e-03 | 1.1e-02 | 6.7e-03 |
| helm | mix_h1_he4 | S | 2.043e+07 | 5.7e-03 | 7.4e-03 | 8.8e-03 | 1.2e-02 |
| helm | mix_c12_o16_ne22 | S | 4.464e+06 | 1.1e-02 | 2.7e-04 | 1.2e-02 | 3.4e-04 |
| helm | mix_c12_o16_90_10 | S | 1.832e+06 | 4.5e-03 | 7.9e-05 | 5.2e-03 | 1.1e-04 |
| helm | mix_c12_o16_99_01 | S | 3.060e+05 | 7.6e-04 | 8.7e-06 | 8.7e-04 | 1.1e-05 |
| helm | mix_he4_c12_99_01 | S | 4.639e+05 | 5.8e-04 | 1.5e-04 | 6.7e-04 | 1.9e-04 |
| helm | mix_he4_c12_10_90 | S | 4.674e+06 | 8.4e-03 | 2.7e-03 | 9.9e-03 | 3.5e-03 |
| helm | mix_c12_o16_ne20 | S | 5.229e+06 | 1.3e-02 | 3.3e-04 | 1.4e-02 | 4.1e-04 |
| helm | mix_he4_c12_o16 | S | 9.528e+06 | 1.6e-02 | 5.5e-03 | 1.8e-02 | 6.6e-03 |
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
| skye | mix_c12_o16_90_10 | S | 1.832e+06 | 2.8e-05 | 1.7e-02 | 8.6e-05 | 1.1e-02 |
| skye | mix_c12_o16_90_10 | E | 2.255e+13 | 6.7e-06 | 1.4e-04 | 1.8e-04 | 6.4e-02 |
| skye | mix_c12_o16_99_01 | S | 3.060e+05 | 4.3e-06 | 2.9e-03 | 1.0e-05 | 1.9e-03 |
| skye | mix_c12_o16_99_01 | E | 2.425e+12 | 6.6e-07 | 1.5e-05 | 1.8e-05 | 7.2e-03 |
| skye | mix_he4_c12_99_01 | S | 4.639e+05 | 9.1e-06 | 1.4e-03 | 4.9e-05 | 1.4e-03 |
| skye | mix_he4_c12_99_01 | E | 6.095e+12 | 1.8e-06 | 4.0e-05 | 8.9e-05 | 1.7e-02 |
| skye | mix_he4_c12_10_90 | S | 4.674e+06 | 1.7e-04 | 2.6e-02 | 9.0e-04 | 2.1e-02 |
| skye | mix_he4_c12_10_90 | E | 1.375e+14 | 5.8e-05 | 9.0e-04 | 1.5e-03 | 3.2e-01 |
| skye | mix_c12_o16_ne20 | S | 5.229e+06 | 5.9e-05 | 5.1e-02 | 2.9e-04 | 2.8e-02 |
| skye | mix_c12_o16_ne20 | E | 1.205e+14 | 3.1e-05 | 7.4e-04 | 7.8e-04 | 2.5e-01 |
| skye | mix_he4_c12_o16 | S | 9.528e+06 | 3.6e-04 | 4.7e-02 | 1.7e-03 | 3.7e-02 |
| skye | mix_he4_c12_o16 | E | 4.052e+14 | 6.9e-05 | 2.6e-03 | 2.5e-03 | 5.4e-01 |
| skye | mix_c12_o16 | S | 4.141e+06 | 2.9e-05 | 2.6e-02 | 3.3e-05 | 3.8e-02 |
| skye | mix_c12_o16 | E | 6.979e+13 | 9.4e-06 | 2.6e-04 | 6.3e-05 | 2.1e-01 |
| skye | mix_o16_ne20_mg24 | S | 3.992e+06 | 6.1e-05 | 3.3e-02 | 6.9e-05 | 4.4e-02 |
| skye | mix_o16_ne20_mg24 | E | 1.377e+14 | 1.3e-05 | 4.4e-04 | 5.3e-04 | 3.1e-01 |
| skye | mix_he4_c12 | S | 7.791e+06 | 3.5e-05 | 2.1e-02 | 3.8e-04 | 3.3e-02 |
| skye | mix_he4_c12 | E | 2.293e+14 | 1.9e-05 | 8.3e-04 | 8.9e-04 | 4.8e-01 |
| skye | mix_h1_he4 | S | 2.043e+07 | 1.7e-04 | 1.3e-02 | 1.4e-04 | 3.0e-02 |
| skye | mix_h1_he4 | E | 1.279e+13 | 1.1e-06 | 2.8e-05 | 1.6e-04 | 3.1e-02 |
| skye | mix_c12_o16_ne22 | S | 4.464e+06 | 3.3e-05 | 2.8e-02 | 3.9e-05 | 4.0e-02 |
| skye | mix_c12_o16_ne22 | E | 8.308e+13 | 1.1e-05 | 3.0e-04 | 7.3e-05 | 2.3e-01 |
| skye | mix_c12_o16_90_10 | S | 1.832e+06 | 1.2e-05 | 1.1e-02 | 1.7e-05 | 1.7e-02 |
| skye | mix_c12_o16_90_10 | E | 2.255e+13 | 3.4e-06 | 8.1e-05 | 2.3e-05 | 8.2e-02 |
| skye | mix_c12_o16_99_01 | S | 3.060e+05 | 1.8e-06 | 1.8e-03 | 2.9e-06 | 2.8e-03 |
| skye | mix_c12_o16_99_01 | E | 2.425e+12 | 3.7e-07 | 8.6e-06 | 2.5e-06 | 9.6e-03 |
| skye | mix_he4_c12_99_01 | S | 4.639e+05 | 1.7e-06 | 1.2e-03 | 1.3e-05 | 1.7e-03 |
| skye | mix_he4_c12_99_01 | E | 6.095e+12 | 2.1e-07 | 2.1e-05 | 2.6e-05 | 2.3e-02 |
| skye | mix_he4_c12_10_90 | S | 4.674e+06 | 2.3e-05 | 2.1e-02 | 2.8e-04 | 2.7e-02 |
| skye | mix_he4_c12_10_90 | E | 1.375e+14 | 4.8e-05 | 4.6e-04 | 5.5e-04 | 3.7e-01 |
| skye | mix_c12_o16_ne20 | S | 5.229e+06 | 4.6e-05 | 3.3e-02 | 5.3e-05 | 4.6e-02 |
| skye | mix_c12_o16_ne20 | E | 1.205e+14 | 1.6e-05 | 4.3e-04 | 2.6e-04 | 3.0e-01 |
| skye | mix_he4_c12_o16 | S | 9.528e+06 | 6.2e-05 | 3.8e-02 | 2.6e-04 | 5.1e-02 |
| skye | mix_he4_c12_o16 | E | 4.052e+14 | 3.7e-05 | 1.3e-03 | 8.9e-04 | 6.1e-01 |
| helm | mix_c12_o16 | S | 4.141e+06 | 6.4e-03 | 1.4e-04 | 1.4e-02 | 2.8e-04 |
| helm | mix_o16_ne20_mg24 | S | 3.992e+06 | 6.4e-03 | 1.0e-04 | 1.3e-02 | 1.9e-04 |
| helm | mix_he4_c12 | S | 7.791e+06 | 6.2e-03 | 3.3e-03 | 1.2e-02 | 6.4e-03 |
| helm | mix_h1_he4 | S | 2.043e+07 | 3.9e-03 | 5.2e-03 | 9.1e-03 | 1.2e-02 |
| helm | mix_c12_o16_ne22 | S | 4.464e+06 | 6.8e-03 | 1.6e-04 | 1.4e-02 | 3.1e-04 |
| helm | mix_c12_o16_90_10 | S | 1.832e+06 | 2.8e-03 | 4.9e-05 | 5.9e-03 | 9.3e-05 |
| helm | mix_c12_o16_99_01 | S | 3.060e+05 | 4.6e-04 | 5.4e-06 | 9.7e-04 | 9.9e-06 |
| helm | mix_he4_c12_99_01 | S | 4.639e+05 | 3.8e-04 | 9.7e-05 | 7.4e-04 | 1.8e-04 |
| helm | mix_he4_c12_10_90 | S | 4.674e+06 | 5.1e-03 | 1.7e-03 | 1.1e-02 | 3.4e-03 |
| helm | mix_c12_o16_ne20 | S | 5.229e+06 | 8.0e-03 | 2.0e-04 | 1.6e-02 | 3.7e-04 |
| helm | mix_he4_c12_o16 | S | 9.528e+06 | 1.0e-02 | 3.4e-03 | 2.0e-02 | 6.4e-03 |

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
| mix_c12_o16_90_10 | interp floor | 1.2e-09 | 2.5e-08 | 1.8e-09 | 1.3e-05 | 5.9e-08 | 3.9e-08 |
|  | observed AV | 5.9e-04 | 1.3e-03 | 1.1e-04 | 1.1e-04 | 1.2e-04 | 5.0e-04 |
| mix_c12_o16_99_01 | interp floor | 1.2e-09 | 2.5e-08 | 1.8e-09 | 1.3e-05 | 6.1e-08 | 4.1e-08 |
|  | observed AV | 6.4e-05 | 1.4e-04 | 1.1e-05 | 1.2e-05 | 1.3e-05 | 5.4e-05 |
| mix_he4_c12_99_01 | interp floor | 1.4e-09* | 1.9e-08* | 7.5e-09* | 1.7e-05* | 1.2e-07* | 1.5e-07* |
|  | observed AV | 3.9e-04 | 8.4e-04 | 1.9e-04 | 6.6e-05 | 7.1e-05 | 3.2e-04 |
| mix_he4_c12_10_90 | interp floor | 1.2e-09* | 2.4e-08* | 2.4e-09* | 1.3e-05* | 6.7e-08* | 5.2e-08* |
|  | observed AV | 9.6e-03 | 2.1e-02 | 3.5e-03 | 1.0e-03 | 2.0e-03 | 7.2e-03 |
| mix_c12_o16_ne20 | interp floor | 1.2e-09* | 2.4e-08* | 5.4e-09* | 1.4e-05* | 6.0e-08* | 8.8e-08* |
|  | observed AV | 2.6e-03 | 5.7e-03 | 4.1e-04 | 3.7e-04 | 5.5e-04 | 2.1e-03 |
| mix_he4_c12_o16 | interp floor | 1.3e-09* | 2.5e-08* | 2.8e-09* | 1.4e-05* | 6.7e-08* | 5.6e-08* |
|  | observed AV | 1.9e-02 | 4.2e-02 | 6.6e-03 | 2.0e-03 | 4.0e-03 | 1.5e-02 |

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
| mix_c12_o16_90_10 | interp floor | 2.8e-09 | 4.3e-08 | 8.4e-09 | 2.3e-05 | 2.3e-08 | 9.5e-08 |
|  | observed AV | 6.1e-05 | 7.8e-05 | 7.9e-05 | 2.1e-05 | 4.6e-06 | 5.1e-05 |
| mix_c12_o16_99_01 | interp floor | 2.8e-09 | 4.3e-08 | 8.2e-09 | 2.3e-05 | 2.3e-08 | 9.3e-08 |
|  | observed AV | 6.6e-06 | 8.6e-06 | 8.7e-06 | 2.4e-06 | 5.0e-07 | 5.6e-06 |
| mix_he4_c12_99_01 | interp floor | 3.8e-10* | 8.3e-09* | 2.3e-08* | 2.8e-05* | 1.4e-08* | 1.7e-06* |
|  | observed AV | 2.9e-05 | 5.4e-05 | 1.5e-04 | 2.5e-05 | 2.6e-06 | 3.0e-05 |
| mix_he4_c12_10_90 | interp floor | 2.6e-09* | 4.0e-08* | 9.7e-09* | 2.3e-05* | 2.2e-08* | 2.6e-07* |
|  | observed AV | 8.4e-04 | 1.2e-03 | 2.7e-03 | 4.3e-04 | 7.1e-05 | 8.2e-04 |
| mix_c12_o16_ne20 | interp floor | 2.4e-09* | 3.9e-08* | 2.6e-08* | 2.4e-05* | 2.4e-08* | 1.1e-06* |
|  | observed AV | 2.9e-04 | 3.6e-04 | 3.3e-04 | 7.7e-05 | 2.2e-05 | 2.2e-04 |
| mix_he4_c12_o16 | interp floor | 2.2e-09* | 3.6e-08* | 1.2e-08* | 2.4e-05* | 2.2e-08* | 4.3e-07* |
|  | observed AV | 1.9e-03 | 2.6e-03 | 5.5e-03 | 7.8e-04 | 1.6e-04 | 1.6e-03 |

## Coverage (integrated-path declines)

`eosDT_get_component(i_eos_WDM_AV)` calls the component directly (no blend box); a node is declined when the AV solver fails (no-converge / bad-domain / out-of-table), in which case WDM_AV returns the all-zero sentinel that `mix_vs_reference` records as flag = -1.

| set | mixture | both-valid / total | both-valid in target | ref-valid but WDM declined |
|---|---|---|---|---|
| skye | mix_c12_o16 | 8396/11946 | 1232/2501 | 2012 |
| skye | mix_o16_ne20_mg24 | 8251/11946 | 1017/2501 | 2182 |
| skye | mix_he4_c12 | 8030/11946 | 1390/2501 | 2877 |
| skye | mix_h1_he4 | 8721/11946 | 1921/2501 | 2729 |
| skye | mix_c12_o16_ne22 | 8230/11946 | 1109/2501 | 2179 |
| skye | mix_c12_o16_90_10 | 8397/11946 | 1232/2501 | 2097 |
| skye | mix_c12_o16_99_01 | 8398/11946 | 1232/2501 | 2124 |
| skye | mix_he4_c12_99_01 | 8034/11946 | 1404/2501 | 3073 |
| skye | mix_he4_c12_10_90 | 8026/11946 | 1382/2501 | 2688 |
| skye | mix_c12_o16_ne20 | 8237/11946 | 1115/2501 | 2192 |
| skye | mix_he4_c12_o16 | 7824/11946 | 1235/2501 | 2865 |
| helm | mix_c12_o16 | 6075/11946 | 840/2501 | 2009 |
| helm | mix_o16_ne20_mg24 | 5778/11946 | 681/2501 | 2115 |
| helm | mix_he4_c12 | 6267/11946 | 945/2501 | 2371 |
| helm | mix_h1_he4 | 6754/11946 | 1359/2501 | 2173 |
| helm | mix_c12_o16_ne22 | 5925/11946 | 753/2501 | 2156 |
| helm | mix_c12_o16_90_10 | 6070/11946 | 835/2501 | 2087 |
| helm | mix_c12_o16_99_01 | 6071/11946 | 835/2501 | 2097 |
| helm | mix_he4_c12_99_01 | 6258/11946 | 943/2501 | 2619 |
| helm | mix_he4_c12_10_90 | 6264/11946 | 948/2501 | 2034 |
| helm | mix_c12_o16_ne20 | 5919/11946 | 754/2501 | 2132 |
| helm | mix_he4_c12_o16 | 6062/11946 | 832/2501 | 2293 |
| skye_fd | mix_c12_o16 | 8257/11946 | 1228/2501 | 2151 |
| skye_fd | mix_o16_ne20_mg24 | 8088/11946 | 1016/2501 | 2345 |
| skye_fd | mix_he4_c12 | 7883/11946 | 1387/2501 | 3024 |
| skye_fd | mix_h1_he4 | 8511/11946 | 1920/2501 | 2939 |
| skye_fd | mix_c12_o16_ne22 | 8070/11946 | 1106/2501 | 2339 |
| skye_fd | mix_c12_o16_90_10 | 8258/11946 | 1228/2501 | 2236 |
| skye_fd | mix_c12_o16_99_01 | 8259/11946 | 1228/2501 | 2263 |
| skye_fd | mix_he4_c12_99_01 | 7888/11946 | 1400/2501 | 3219 |
| skye_fd | mix_he4_c12_10_90 | 7879/11946 | 1375/2501 | 2835 |
| skye_fd | mix_c12_o16_ne20 | 8098/11946 | 1111/2501 | 2331 |
| skye_fd | mix_he4_c12_o16 | 7678/11946 | 1231/2501 | 3011 |
| skye_ne | mix_c12_o16 | 9409/11946 | 1525/2501 | 999 |
| skye_ne | mix_o16_ne20_mg24 | 9319/11946 | 1361/2501 | 1114 |
| skye_ne | mix_he4_c12 | 9404/11946 | 1634/2501 | 1503 |
| skye_ne | mix_h1_he4 | 9906/11946 | 2041/2501 | 1544 |
| skye_ne | mix_c12_o16_ne22 | 9282/11946 | 1428/2501 | 1127 |
| skye_ne | mix_c12_o16_90_10 | 9409/11946 | 1525/2501 | 1085 |
| skye_ne | mix_c12_o16_99_01 | 9409/11946 | 1525/2501 | 1113 |
| skye_ne | mix_he4_c12_99_01 | 9404/11946 | 1634/2501 | 1703 |
| skye_ne | mix_he4_c12_10_90 | 9404/11946 | 1634/2501 | 1310 |
| skye_ne | mix_c12_o16_ne20 | 9263/11946 | 1429/2501 | 1166 |
| skye_ne | mix_he4_c12_o16 | 9223/11946 | 1525/2501 | 1466 |
| helm_ne | mix_c12_o16 | 7444/11946 | 943/2501 | 640 |
| helm_ne | mix_o16_ne20_mg24 | 7150/11946 | 788/2501 | 743 |
| helm_ne | mix_he4_c12 | 7643/11946 | 1053/2501 | 995 |
| helm_ne | mix_h1_he4 | 8043/11946 | 1453/2501 | 884 |
| helm_ne | mix_c12_o16_ne22 | 7299/11946 | 864/2501 | 782 |
| helm_ne | mix_c12_o16_90_10 | 7444/11946 | 943/2501 | 713 |
| helm_ne | mix_c12_o16_99_01 | 7444/11946 | 943/2501 | 724 |
| helm_ne | mix_he4_c12_99_01 | 7643/11946 | 1053/2501 | 1234 |
| helm_ne | mix_he4_c12_10_90 | 7643/11946 | 1053/2501 | 655 |
| helm_ne | mix_c12_o16_ne20 | 7286/11946 | 859/2501 | 765 |
| helm_ne | mix_he4_c12_o16 | 7441/11946 | 943/2501 | 914 |
