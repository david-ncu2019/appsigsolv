# Timeseries Decomposition Report: LNJS_neu — dU

## Accepted Model

| Parameter | Value |
|---|---|
| Component | dU |
| Sigma_0 assumed (mm) | 6.0 |
| Sigma_hat a-posteriori (mm) | 5.35 |
| Polynomial degree | 1 |
| Seasonal periods (yr) | [0.25, 0.5, 1.0, 2.0, 3.78, 5.27, 10.58] |
| Jump dates | ['20221007'] |
| Polyline breaks | [] |
| Exp trend (b/day) | None |
| Exp relaxation | {} |
| Log relaxation | {} |
| n_params | 17 |
| Degrees of freedom (r) | 5028 |
| T_stat (SSR/σ²) | 4003.3807 |
| χ²_critical (K) | 5194.0745 |
| Unit variance factor (T/r) | 0.7962 |
| p-value | 1.000000 |
| DIA iterations | 1 |

## Variance Explained per Component (%)

| Component | Std (mm) | Variance Explained (%) |
|---|---|---|
| dU_trend | 11.86 | 17.27 |
| dU_0.25yr | 0.08 | 0.00 |
| dU_0.5yr | 0.76 | 0.07 |
| dU_1yr | 0.97 | 0.12 |
| dU_2yr | 0.23 | 0.01 |
| dU_3.78yr | 0.37 | 0.02 |
| dU_5.27yr | 0.57 | 0.04 |
| dU_10.58yr | 0.72 | 0.06 |
| dU_jump | 31.74 | 123.76 |
| dU_noise | 5.34 | 3.51 |

## Sigma Scan Summary

| Sigma (mm) | Accepted | T_stat | p-value | n_params | n_periods | n_polylines |
|---|---|---|---|---|---|---|
| 2.0 | ✗ | — | — | — | — | — |
| 3.0 | ✗ | — | — | — | — | — |
| 4.0 | ✗ | — | — | — | — | — |
| 5.0 | ✗ | — | — | — | — | — |
| 6.0 | ✓ | 4003.38 | 1.0000 | 17 | 7 | 0 |
| 7.0 | ✓ | 2941.26 | 1.0000 | 17 | 7 | 0 |
| 8.0 | ✓ | 2251.90 | 1.0000 | 17 | 7 | 0 |
| 9.0 | ✓ | 1779.28 | 1.0000 | 17 | 7 | 0 |
| 10.0 | ✓ | 1441.22 | 1.0000 | 17 | 7 | 0 |
| 11.0 | ✓ | 1191.09 | 1.0000 | 17 | 7 | 0 |
| 12.0 | ✓ | 1000.85 | 1.0000 | 17 | 7 | 0 |
| 13.0 | ✓ | 852.79 | 1.0000 | 17 | 7 | 0 |
| 14.0 | ✓ | 735.31 | 1.0000 | 17 | 7 | 0 |
| 15.0 | ✓ | 640.54 | 1.0000 | 17 | 7 | 0 |
| 16.0 | ✓ | 562.98 | 1.0000 | 17 | 7 | 0 |
| 17.0 | ✓ | 498.69 | 1.0000 | 17 | 7 | 0 |
| 18.0 | ✓ | 444.82 | 1.0000 | 17 | 7 | 0 |
| 19.0 | ✓ | 399.23 | 1.0000 | 17 | 7 | 0 |
| 20.0 | ✓ | 360.30 | 1.0000 | 17 | 7 | 0 |
