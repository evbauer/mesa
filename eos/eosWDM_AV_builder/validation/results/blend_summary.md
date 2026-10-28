# WDM_AV blend-zone analysis (full eosDT cascade)

Internal consistency of `eosDT_get` across the WDM_AV blend box: returned partials vs central finite differences of the returned values, and the value step at each edge (pure-WDM_AV trend minus pure-neighbour trend, both extrapolated to the ramp centre).

## (a) FD consistency: max relative mismatch (returned vs central FD)

| set | mixture | cut | quantity | max in-ramp | max out-ramp |
|---|---|---|---|---|---|
| skye | c12_o16 | rho_right_T6.0 | lnPgas | 6.45e-03 | 1.86e-04 |
| skye | c12_o16 | rho_right_T6.0 | lnE | 3.50e-02 | 8.74e-04 |
| skye | c12_o16 | rho_right_T6.0 | lnS | 6.49e+00 | 2.81e-02 |
| skye | c12_o16 | rho_right_T6.0 | Cv | 8.37e-01 | 2.72e+00 |
| skye | c12_o16 | rho_right_T6.0 | grad_ad | 4.08e-02 | 1.13e-01 |
| skye | c12_o16 | rho_right_T6.8 | lnPgas | 2.24e-07 | 1.60e-06 |
| skye | c12_o16 | rho_right_T6.8 | lnE | 5.27e-07 | 4.87e-06 |
| skye | c12_o16 | rho_right_T6.8 | lnS | 1.42e-04 | 1.92e-04 |
| skye | c12_o16 | rho_right_T6.8 | Cv | 2.22e-06 | 1.55e-05 |
| skye | c12_o16 | rho_right_T6.8 | grad_ad | 3.01e-05 | 1.37e-04 |
| skye | c12_o16 | T_top_R2.0 | lnPgas | 7.34e-04 | 1.65e-03 |
| skye | c12_o16 | T_top_R2.0 | lnE | 6.09e-04 | 1.61e-02 |
| skye | c12_o16 | T_top_R2.0 | lnS | 1.08e-03 | 2.60e-02 |
| skye | c12_o16 | T_top_R2.0 | Cv | 3.39e-05 | 4.62e-02 |
| skye | c12_o16 | T_top_R2.0 | grad_ad | 5.54e-05 | 5.39e-01 |
| skye | c12_o16 | T_top_R4.0 | lnPgas | 7.78e-06 | 4.16e-05 |
| skye | c12_o16 | T_top_R4.0 | lnE | 7.59e-06 | 2.31e-05 |
| skye | c12_o16 | T_top_R4.0 | lnS | 5.42e-05 | 1.35e-03 |
| skye | c12_o16 | T_top_R4.0 | Cv | 2.38e-04 | 3.81e-01 |
| skye | c12_o16 | T_top_R4.0 | grad_ad | 1.13e-04 | 5.10e-01 |
| skye | c12_o16 | rho_left_T6.8 | lnPgas | 3.03e-07 | 4.15e-07 |
| skye | c12_o16 | rho_left_T6.8 | lnE | 1.81e-06 | 2.12e-06 |
| skye | c12_o16 | rho_left_T6.8 | lnS | 1.83e-06 | 1.84e-06 |
| skye | c12_o16 | rho_left_T6.8 | Cv | 2.16e-06 | 2.26e-05 |
| skye | c12_o16 | rho_left_T6.8 | grad_ad | 3.34e-06 | 8.38e-06 |
| skye | c12 | rho_right_T6.0 | lnPgas | 7.46e-03 | 4.78e-05 |
| skye | c12 | rho_right_T6.0 | lnE | 4.02e-02 | 2.24e-04 |
| skye | c12 | rho_right_T6.0 | lnS | 9.83e+00 | 7.93e-02 |
| skye | c12 | rho_right_T6.0 | Cv | 6.01e+00 | 6.52e+00 |
| skye | c12 | rho_right_T6.0 | grad_ad | 7.70e-01 | 7.26e-01 |
| skye | c12 | rho_right_T6.8 | lnPgas | 2.35e-07 | 1.54e-06 |
| skye | c12 | rho_right_T6.8 | lnE | 5.55e-07 | 4.63e-06 |
| skye | c12 | rho_right_T6.8 | lnS | 1.34e-06 | 3.67e-06 |
| skye | c12 | rho_right_T6.8 | Cv | 2.45e-06 | 1.55e-05 |
| skye | c12 | rho_right_T6.8 | grad_ad | 1.59e-04 | 3.77e-04 |
| skye | c12 | T_top_R2.0 | lnPgas | 1.85e-03 | 5.20e-05 |
| skye | c12 | T_top_R2.0 | lnE | 9.43e-04 | 3.45e-05 |
| skye | c12 | T_top_R2.0 | lnS | 4.07e-03 | 7.53e-05 |
| skye | c12 | T_top_R2.0 | Cv | 5.93e-05 | 1.98e-04 |
| skye | c12 | T_top_R2.0 | grad_ad | 2.60e-03 | 2.59e-04 |
| skye | c12 | T_top_R4.0 | lnPgas | 7.65e-06 | 7.22e-06 |
| skye | c12 | T_top_R4.0 | lnE | 7.44e-06 | 7.10e-06 |
| skye | c12 | T_top_R4.0 | lnS | 7.59e-06 | 7.96e-06 |
| skye | c12 | T_top_R4.0 | Cv | 2.46e-04 | 3.11e-04 |
| skye | c12 | T_top_R4.0 | grad_ad | 2.34e-03 | 3.39e-03 |
| skye | c12 | rho_left_T6.8 | lnPgas | 1.13e-05 | 1.49e-05 |
| skye | c12 | rho_left_T6.8 | lnE | 2.29e-06 | 2.12e-06 |
| skye | c12 | rho_left_T6.8 | lnS | 7.29e-06 | 1.06e-05 |
| skye | c12 | rho_left_T6.8 | Cv | 2.15e-06 | 2.26e-05 |
| skye | c12 | rho_left_T6.8 | grad_ad | 3.82e-06 | 9.47e-06 |
| helm | c12_o16 | rho_right_T6.0 | lnPgas | 6.53e-03 | 2.41e-05 |
| helm | c12_o16 | rho_right_T6.0 | lnE | 2.93e-02 | 9.44e-05 |
| helm | c12_o16 | rho_right_T6.0 | lnS | 1.51e+01 | 5.90e-02 |
| helm | c12_o16 | rho_right_T6.0 | Cv | 8.72e+00 | 2.84e+01 |
| helm | c12_o16 | rho_right_T6.0 | grad_ad | 1.66e+01 | 3.03e+01 |
| helm | c12_o16 | rho_right_T6.8 | lnPgas | 4.45e-06 | 6.20e-06 |
| helm | c12_o16 | rho_right_T6.8 | lnE | 8.82e-06 | 1.14e-05 |
| helm | c12_o16 | rho_right_T6.8 | lnS | 1.69e-04 | 2.29e-04 |
| helm | c12_o16 | rho_right_T6.8 | Cv | 2.01e-05 | 2.82e-05 |
| helm | c12_o16 | rho_right_T6.8 | grad_ad | 3.94e-04 | 5.09e-04 |
| helm | c12_o16 | T_top_R2.0 | lnPgas | 7.41e-04 | 2.60e-01 |
| helm | c12_o16 | T_top_R2.0 | lnE | 9.92e-03 | 7.42e-01 |
| helm | c12_o16 | T_top_R2.0 | lnS | 7.88e-03 | 1.84e+00 |
| helm | c12_o16 | T_top_R2.0 | Cv | 3.40e-01 | 4.54e-03 |
| helm | c12_o16 | T_top_R2.0 | grad_ad | 6.49e-01 | 1.26e-02 |
| helm | c12_o16 | T_top_R4.0 | lnPgas | 1.55e-05 | 1.81e-05 |
| helm | c12_o16 | T_top_R4.0 | lnE | 7.77e-06 | 8.40e-06 |
| helm | c12_o16 | T_top_R4.0 | lnS | 5.48e-05 | 8.04e-05 |
| helm | c12_o16 | T_top_R4.0 | Cv | 2.31e-04 | 3.11e-04 |
| helm | c12_o16 | T_top_R4.0 | grad_ad | 1.70e-03 | 2.45e-03 |
| helm | c12_o16 | rho_left_T6.8 | lnPgas | 4.89e-07 | 6.89e-07 |
| helm | c12_o16 | rho_left_T6.8 | lnE | 1.00e-05 | 1.71e-05 |
| helm | c12_o16 | rho_left_T6.8 | lnS | 1.83e-06 | 1.84e-06 |
| helm | c12_o16 | rho_left_T6.8 | Cv | 2.16e-06 | 2.26e-05 |
| helm | c12_o16 | rho_left_T6.8 | grad_ad | 3.34e-06 | 8.38e-06 |
| helm | c12 | rho_right_T6.0 | lnPgas | 7.55e-03 | 4.83e-05 |
| helm | c12 | rho_right_T6.0 | lnE | 3.55e-02 | 1.99e-04 |
| helm | c12 | rho_right_T6.0 | lnS | 9.83e+00 | 7.93e-02 |
| helm | c12 | rho_right_T6.0 | Cv | 4.79e+01 | 5.19e+01 |
| helm | c12 | rho_right_T6.0 | grad_ad | 5.14e+01 | 5.20e+01 |
| helm | c12 | rho_right_T6.8 | lnPgas | 4.43e-06 | 6.16e-06 |
| helm | c12 | rho_right_T6.8 | lnE | 1.14e-05 | 1.48e-05 |
| helm | c12 | rho_right_T6.8 | lnS | 2.69e-05 | 3.62e-05 |
| helm | c12 | rho_right_T6.8 | Cv | 2.05e-05 | 2.87e-05 |
| helm | c12 | rho_right_T6.8 | grad_ad | 3.87e-04 | 4.99e-04 |
| helm | c12 | T_top_R2.0 | lnPgas | 1.86e-03 | 6.68e-05 |
| helm | c12 | T_top_R2.0 | lnE | 1.94e-02 | 2.28e-05 |
| helm | c12 | T_top_R2.0 | lnS | 1.54e-02 | 7.53e-05 |
| helm | c12 | T_top_R2.0 | Cv | 6.25e-01 | 1.82e-04 |
| helm | c12 | T_top_R2.0 | grad_ad | 6.67e-01 | 1.13e-05 |
| helm | c12 | T_top_R4.0 | lnPgas | 1.51e-05 | 1.75e-05 |
| helm | c12 | T_top_R4.0 | lnE | 1.08e-05 | 1.23e-05 |
| helm | c12 | T_top_R4.0 | lnS | 1.01e-05 | 1.13e-05 |
| helm | c12 | T_top_R4.0 | Cv | 2.40e-04 | 3.25e-04 |
| helm | c12 | T_top_R4.0 | grad_ad | 1.75e-03 | 2.51e-03 |
| helm | c12 | rho_left_T6.8 | lnPgas | 1.15e-05 | 1.51e-05 |
| helm | c12 | rho_left_T6.8 | lnE | 9.05e-06 | 1.56e-05 |
| helm | c12 | rho_left_T6.8 | lnS | 7.29e-06 | 1.06e-05 |
| helm | c12 | rho_left_T6.8 | Cv | 2.15e-06 | 2.26e-05 |
| helm | c12 | rho_left_T6.8 | grad_ad | 3.82e-06 | 9.46e-06 |

## (b) value steps at the ramp centre

Jump = (WDM_AV-side trend) - (neighbour-side trend) at the ramp centre. For lnPgas/lnE/lnS the jump is d(lnX), i.e. the relative step in X; for Cv/grad_ad the relative step is jump/neighbour.

| set | mixture | cut | neighbour | lnPgas | lnE | lnS | lnS(exp) | Cv(rel) | grad_ad(rel) |
|---|---|---|---|---|---|---|---|---|---|
| skye | c12_o16 | rho_right_T6.0 | Skye | 5.57e-04 | -7.43e-05 | -2.73e-02 | -9.07e-02 | -1.82e+00 | -8.69e-01 |
| skye | c12_o16 | rho_right_T6.8 | Skye | 1.99e-03 | 1.73e-03 | -4.74e-02 | -4.82e-02 | 9.03e-03 | -1.04e-03 |
| skye | c12_o16 | T_top_R2.0 | Skye | -1.33e-03 | -5.83e-02 | -2.83e-02 | -1.18e-02 | -2.16e-01 | 2.56e-02 |
| skye | c12_o16 | T_top_R4.0 | Skye | -1.13e-02 | -1.01e-02 | -2.37e-02 | -2.85e-02 | 1.79e-02 | 1.06e-03 |
| skye | c12_o16 | rho_left_T6.8 | FreeEOS | 4.45e-04 | 8.17e-03 | 5.45e-03 |   -   | -1.40e-01 | 1.39e-03 |
| skye | c12 | rho_right_T6.0 | Skye | 1.12e-03 | 1.12e-03 | 1.77e-01 | 0.00e+00 | -1.27e+00 | -7.87e-01 |
| skye | c12 | rho_right_T6.8 | Skye | 2.09e-03 | 1.79e-03 | 1.62e-03 | 0.00e+00 | 8.65e-03 | -1.08e-03 |
| skye | c12 | T_top_R2.0 | Skye | -2.08e-04 | -5.77e-02 | -1.58e-02 | 0.00e+00 | -2.13e-01 | 2.63e-02 |
| skye | c12 | T_top_R4.0 | Skye | -1.10e-02 | -1.00e-02 | 4.41e-03 | 0.00e+00 | 1.74e-02 | 9.82e-04 |
| skye | c12 | rho_left_T6.8 | FreeEOS | 1.86e-02 | 6.50e-03 | 1.44e-02 |   -   | -1.40e-01 | 1.46e-03 |
| helm | c12_o16 | rho_right_T6.0 | Skye | 1.14e-02 | 9.74e-03 | -7.34e-03 | -9.07e-02 | -1.62e+00 | -6.15e-01 |
| helm | c12_o16 | rho_right_T6.8 | Skye | 1.23e-02 | 1.10e-02 | -5.79e-02 | -4.82e-02 | 2.91e-02 | -2.33e-02 |
| helm | c12_o16 | T_top_R2.0 | HELM | -2.20e-02 | -1.64e-02 | 2.82e-02 | 2.80e-02 | 4.65e-02 | -7.48e-03 |
| helm | c12_o16 | T_top_R4.0 | Skye | 7.48e-03 | -4.59e-03 | -2.22e-02 | -2.85e-02 | 3.42e-02 | 6.49e-03 |
| helm | c12_o16 | rho_left_T6.8 | FreeEOS | 7.05e-04 | -6.49e-03 | 5.46e-03 |   -   | -1.40e-01 | 1.39e-03 |
| helm | c12 | rho_right_T6.0 | Skye | 1.16e-02 | 1.25e-02 | 1.42e-01 | 0.00e+00 | -1.27e+00 | -7.85e-01 |
| helm | c12 | rho_right_T6.8 | Skye | 1.22e-02 | 1.34e-02 | -7.98e-03 | 0.00e+00 | 2.75e-02 | -2.16e-02 |
| helm | c12 | T_top_R2.0 | Skye | 1.11e-02 | -1.07e-01 | -1.44e-02 | 0.00e+00 | -2.17e-01 | 3.20e-02 |
| helm | c12 | T_top_R4.0 | Skye | 7.09e-03 | 1.74e-03 | 5.71e-03 | 0.00e+00 | 3.30e-02 | 5.81e-03 |
| helm | c12 | rho_left_T6.8 | FreeEOS | 1.88e-02 | -5.28e-03 | 1.44e-02 |   -   | -1.40e-01 | 1.46e-03 |
