# Held-out-node interpolation error

Relative error of each quantity at fine nodes omitted from the 10-points/decade coarse grid, vs the source values. Variants: a_o5 = biquintic/source 2nd derivs, b_o5 = biquintic/FD 2nd derivs, a_o3 = bicubic. Full per-region data (incl. p90 and the whole valid plane) is in heldout_summary.csv / .json.

## Skye-valid region: logT >= 6.5 OR logRho >= 3.7

| table | variant | stat | P | s | c_V | chi_rho | chi_T |
|---|---|---|---|---|---|---|---|
| skye_he4 | a_o5 | median | 1.3e-11 | 1.3e-10 | 1.3e-06 | 3.4e-09 | 1.4e-09 |
| skye_he4 | a_o5 | p99 | 1.2e-03 | 7.0e-02 | 2.1e+00 | 1.7e+02 | 9.6e-01 |
| skye_he4 | a_o5 | max | 4.7e-02 | 1.6e+00 | 1.1e+03 | 1.4e+05 | 6.3e+03 |
| skye_he4 | b_o5 | median | 3.9e-10 | 7.2e-09 | 6.3e-05 | 2.4e-08 | 7.0e-06 |
| skye_he4 | b_o5 | p99 | 4.0e-02 | 8.7e-02 | 2.1e+00 | 4.5e+01 | 2.3e+00 |
| skye_he4 | b_o5 | max | 7.8e-01 | 1.8e+00 | 1.5e+03 | 1.6e+05 | 1.7e+04 |
| skye_he4 | a_o3 | median | 1.2e-08 | 1.8e-07 | 8.7e-03 | 2.4e-05 | 3.3e-07 |
| skye_he4 | a_o3 | p99 | 1.7e-02 | 5.7e-02 | 2.0e+00 | 5.0e+01 | 1.0e+00 |
| skye_he4 | a_o3 | max | 8.4e-02 | 1.8e+00 | 1.4e+03 | 1.5e+04 | 5.0e+03 |
| skye_c12 | a_o5 | median | 1.6e-11 | 4.5e-10 | 1.3e-06 | 2.8e-09 | 5.7e-09 |
| skye_c12 | a_o5 | p99 | 1.2e-03 | 4.8e-02 | 2.3e+00 | 3.0e+02 | 9.4e-01 |
| skye_c12 | a_o5 | max | 5.0e-02 | 2.6e+00 | 1.7e+03 | 1.2e+05 | 8.5e+03 |
| skye_c12 | b_o5 | median | 2.3e-10 | 1.2e-08 | 6.8e-05 | 2.5e-08 | 1.4e-05 |
| skye_c12 | b_o5 | p99 | 3.4e-02 | 5.4e-02 | 1.8e+00 | 1.0e+02 | 1.9e+00 |
| skye_c12 | b_o5 | max | 7.9e-01 | 2.7e+00 | 1.5e+03 | 2.0e+05 | 3.1e+04 |
| skye_c12 | a_o3 | median | 1.1e-08 | 2.2e-07 | 7.8e-03 | 2.7e-05 | 3.7e-07 |
| skye_c12 | a_o3 | p99 | 1.7e-02 | 4.0e-02 | 2.0e+00 | 8.0e+01 | 9.7e-01 |
| skye_c12 | a_o3 | max | 9.0e-02 | 2.5e+00 | 6.9e+02 | 1.5e+04 | 7.8e+03 |
| skye_o16 | a_o5 | median | 1.7e-11 | 6.4e-10 | 1.3e-06 | 2.8e-09 | 6.9e-09 |
| skye_o16 | a_o5 | p99 | 1.2e-03 | 5.1e-02 | 2.3e+00 | 3.5e+02 | 9.6e-01 |
| skye_o16 | a_o5 | max | 5.1e-02 | 2.3e+00 | 2.5e+03 | 1.4e+05 | 6.4e+02 |
| skye_o16 | b_o5 | median | 2.2e-10 | 1.5e-08 | 6.8e-05 | 2.5e-08 | 1.6e-05 |
| skye_o16 | b_o5 | p99 | 3.5e-02 | 5.8e-02 | 1.8e+00 | 1.3e+02 | 1.8e+00 |
| skye_o16 | b_o5 | max | 7.9e-01 | 2.3e+00 | 2.1e+03 | 2.1e+05 | 2.0e+03 |
| skye_o16 | a_o3 | median | 1.1e-08 | 2.3e-07 | 7.6e-03 | 2.7e-05 | 3.7e-07 |
| skye_o16 | a_o3 | p99 | 1.7e-02 | 4.3e-02 | 2.0e+00 | 9.0e+01 | 1.0e+00 |
| skye_o16 | a_o3 | max | 9.1e-02 | 2.1e+00 | 6.7e+02 | 1.6e+04 | 5.0e+02 |

## Skye-valid region excluding the pair-radiation corner (NOT (logT > 8.0 AND logRho < 0))

| table | variant | stat | P | s | c_V | chi_rho | chi_T |
|---|---|---|---|---|---|---|---|
| skye_he4 | a_o5 | median | 1.2e-11 | 9.7e-11 | 1.2e-06 | 2.2e-09 | 1.1e-09 |
| skye_he4 | a_o5 | p99 | 1.2e-04 | 8.0e-02 | 2.3e+00 | 1.9e-02 | 1.1e+00 |
| skye_he4 | a_o5 | max | 1.9e-03 | 1.6e+00 | 1.1e+03 | 4.4e+01 | 6.3e+03 |
| skye_he4 | b_o5 | median | 1.9e-10 | 7.1e-09 | 5.9e-05 | 2.4e-08 | 6.7e-06 |
| skye_he4 | b_o5 | p99 | 5.9e-03 | 9.8e-02 | 2.4e+00 | 9.3e-03 | 2.5e+00 |
| skye_he4 | b_o5 | max | 5.7e-02 | 1.8e+00 | 1.5e+03 | 7.6e+01 | 1.7e+04 |
| skye_he4 | a_o3 | median | 9.1e-09 | 1.8e-07 | 8.1e-03 | 2.3e-05 | 3.1e-07 |
| skye_he4 | a_o3 | p99 | 2.1e-03 | 6.8e-02 | 2.1e+00 | 1.3e-02 | 1.1e+00 |
| skye_he4 | a_o3 | max | 1.1e-02 | 1.8e+00 | 1.4e+03 | 4.8e+00 | 5.0e+03 |
| skye_c12 | a_o5 | median | 1.2e-11 | 3.4e-10 | 1.3e-06 | 1.8e-09 | 3.2e-09 |
| skye_c12 | a_o5 | p99 | 1.4e-04 | 5.3e-02 | 2.5e+00 | 2.2e-02 | 1.1e+00 |
| skye_c12 | a_o5 | max | 2.0e-03 | 2.6e+00 | 1.7e+03 | 1.3e+02 | 8.5e+03 |
| skye_c12 | b_o5 | median | 1.3e-10 | 1.1e-08 | 6.3e-05 | 2.4e-08 | 1.3e-05 |
| skye_c12 | b_o5 | p99 | 5.9e-03 | 5.5e-02 | 2.0e+00 | 1.3e-02 | 2.0e+00 |
| skye_c12 | b_o5 | max | 6.7e-02 | 2.7e+00 | 1.5e+03 | 2.3e+02 | 3.1e+04 |
| skye_c12 | a_o3 | median | 7.2e-09 | 2.1e-07 | 7.3e-03 | 2.7e-05 | 3.7e-07 |
| skye_c12 | a_o3 | p99 | 2.3e-03 | 4.6e-02 | 2.1e+00 | 2.2e-02 | 1.1e+00 |
| skye_c12 | a_o3 | max | 1.1e-02 | 2.5e+00 | 6.9e+02 | 1.5e+01 | 7.8e+03 |
| skye_o16 | a_o5 | median | 1.2e-11 | 4.9e-10 | 1.3e-06 | 1.8e-09 | 4.0e-09 |
| skye_o16 | a_o5 | p99 | 1.5e-04 | 5.7e-02 | 2.5e+00 | 2.3e-02 | 1.1e+00 |
| skye_o16 | a_o5 | max | 2.1e-03 | 2.3e+00 | 2.5e+03 | 1.7e+02 | 6.4e+02 |
| skye_o16 | b_o5 | median | 1.3e-10 | 1.3e-08 | 6.4e-05 | 2.5e-08 | 1.5e-05 |
| skye_o16 | b_o5 | p99 | 6.0e-03 | 5.9e-02 | 2.0e+00 | 1.5e-02 | 2.0e+00 |
| skye_o16 | b_o5 | max | 6.9e-02 | 2.3e+00 | 2.1e+03 | 3.0e+02 | 2.0e+03 |
| skye_o16 | a_o3 | median | 7.3e-09 | 2.3e-07 | 7.1e-03 | 2.7e-05 | 3.7e-07 |
| skye_o16 | a_o3 | p99 | 2.3e-03 | 4.8e-02 | 2.0e+00 | 2.4e-02 | 1.1e+00 |
| skye_o16 | a_o3 | max | 1.1e-02 | 2.1e+00 | 6.7e+02 | 1.9e+01 | 5.0e+02 |

## WDM_AV target region: logRho in [-2, 4] AND logT in [3, 7]

| table | variant | stat | P | s | c_V | chi_rho | chi_T |
|---|---|---|---|---|---|---|---|
| skye_he4 | a_o5 | median | 2.0e-11 | 1.2e-10 | 9.6e-07 | 7.7e-09 | 1.7e-09 |
| skye_he4 | a_o5 | p99 | 2.3e-04 | 7.0e-04 | 2.1e-01 | 1.4e-03 | 9.3e-02 |
| skye_he4 | a_o5 | max | 1.4e-01 | 1.7e-01 | 4.1e+01 | 1.7e+00 | 1.3e+01 |
| skye_he4 | b_o5 | median | 4.0e-09 | 5.6e-09 | 4.4e-05 | 2.2e-07 | 6.6e-06 |
| skye_he4 | b_o5 | p99 | 4.9e-04 | 5.5e-03 | 1.4e+00 | 2.6e-03 | 5.3e-01 |
| skye_he4 | b_o5 | max | 1.5e-01 | 1.7e-01 | 3.4e+01 | 1.4e+00 | 1.2e+01 |
| skye_he4 | a_o3 | median | 3.8e-08 | 1.6e-07 | 5.6e-03 | 1.0e-04 | 2.5e-07 |
| skye_he4 | a_o3 | p99 | 2.5e-05 | 6.0e-04 | 8.7e-01 | 8.9e-04 | 4.9e-02 |
| skye_he4 | a_o3 | max | 1.4e-01 | 1.7e-01 | 2.9e+01 | 2.3e+00 | 1.0e+01 |
| skye_c12 | a_o5 | median | 2.0e-11 | 2.6e-10 | 9.1e-07 | 8.2e-09 | 3.2e-09 |
| skye_c12 | a_o5 | p99 | 4.3e-04 | 3.8e-04 | 1.9e-01 | 3.1e-03 | 9.2e-02 |
| skye_c12 | a_o5 | max | 8.8e-01 | 1.2e-01 | 2.6e+01 | 9.4e+00 | 5.2e+02 |
| skye_c12 | b_o5 | median | 2.1e-09 | 9.1e-09 | 4.9e-05 | 2.4e-07 | 8.8e-06 |
| skye_c12 | b_o5 | p99 | 5.2e-04 | 2.3e-03 | 1.1e+00 | 3.0e-03 | 2.9e-01 |
| skye_c12 | b_o5 | max | 3.8e-01 | 8.4e-02 | 2.0e+01 | 9.6e-01 | 2.8e+01 |
| skye_c12 | a_o3 | median | 3.6e-08 | 1.9e-07 | 5.4e-03 | 1.1e-04 | 3.7e-07 |
| skye_c12 | a_o3 | p99 | 4.9e-05 | 3.7e-04 | 7.7e-01 | 1.0e-03 | 6.1e-02 |
| skye_c12 | a_o3 | max | 1.4e+00 | 1.3e-01 | 2.9e+01 | 2.4e+01 | 5.3e+01 |
| skye_o16 | a_o5 | median | 2.0e-11 | 2.9e-10 | 9.1e-07 | 8.4e-09 | 3.2e-09 |
| skye_o16 | a_o5 | p99 | 4.5e-04 | 4.2e-04 | 1.5e-01 | 4.1e-03 | 7.8e-02 |
| skye_o16 | a_o5 | max | 4.9e-01 | 1.0e-01 | 1.2e+01 | 3.1e+00 | 3.2e+01 |
| skye_o16 | b_o5 | median | 1.6e-09 | 1.1e-08 | 4.9e-05 | 2.5e-07 | 9.6e-06 |
| skye_o16 | b_o5 | p99 | 7.1e-04 | 2.0e-03 | 9.8e-01 | 3.9e-03 | 2.5e-01 |
| skye_o16 | b_o5 | max | 2.5e-01 | 8.6e-02 | 9.4e+00 | 1.4e+00 | 2.4e+01 |
| skye_o16 | a_o3 | median | 3.7e-08 | 2.0e-07 | 5.3e-03 | 1.2e-04 | 3.7e-07 |
| skye_o16 | a_o3 | p99 | 6.3e-05 | 2.4e-04 | 3.1e-01 | 1.3e-03 | 4.2e-02 |
| skye_o16 | a_o3 | max | 5.8e+00 | 1.2e-01 | 1.3e+01 | 7.6e+00 | 2.7e+01 |
