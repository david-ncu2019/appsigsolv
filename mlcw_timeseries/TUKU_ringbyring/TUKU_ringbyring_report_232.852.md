# Timeseries Decomposition Report: TUKU_ringbyring — 232.852

## Accepted Model

| Parameter | Value |
|---|---|
| Component | 232.852 |
| Sigma_0 assumed (mm) | 2.0 |
| Sigma_hat a-posteriori (mm) | 2.10 |
| Polynomial degree | 2 |
| Seasonal periods (yr) | [0.5, 1.0, 1.97, 2.5] |
| Jump dates | [] |
| Polyline breaks | [] |
| Exp trend (b/day) | None |
| Exp relaxation | {} |
| Log relaxation | {} |
| n_params | 11 |
| Degrees of freedom (r) | 241 |
| T_stat (SSR/σ²) | 265.2986 |
| χ²_critical (K) | 278.2127 |
| Unit variance factor (T/r) | 1.1008 |
| p-value | 0.135287 |
| DIA iterations | 0 |

## Variance Explained per Component (%)

| Component | Std (mm) | Variance Explained (%) |
|---|---|---|
| 232.852_trend | 9.68 | 92.60 |
| 232.852_0.5yr | 0.04 | 0.00 |
| 232.852_1yr | 0.33 | 0.11 |
| 232.852_1.97yr | 0.69 | 0.47 |
| 232.852_2.5yr | 0.96 | 0.91 |
| 232.852_noise | 2.05 | 4.16 |

## Sigma Scan Summary

| Sigma (mm) | Accepted | T_stat | p-value | n_params | n_periods | n_polylines |
|---|---|---|---|---|---|---|
| 2.0 | ✓ | 265.30 | 0.1353 | 11 | 4 | 0 |
| 3.0 | ✓ | 117.91 | 1.0000 | 11 | 4 | 0 |
| 4.0 | ✓ | 66.32 | 1.0000 | 11 | 4 | 0 |
| 5.0 | ✓ | 42.45 | 1.0000 | 11 | 4 | 0 |
| 6.0 | ✓ | 29.48 | 1.0000 | 11 | 4 | 0 |
| 7.0 | ✓ | 21.66 | 1.0000 | 11 | 4 | 0 |
| 8.0 | ✓ | 16.58 | 1.0000 | 11 | 4 | 0 |
| 9.0 | ✓ | 13.10 | 1.0000 | 11 | 4 | 0 |
| 10.0 | ✓ | 10.61 | 1.0000 | 11 | 4 | 0 |
| 11.0 | ✓ | 8.77 | 1.0000 | 11 | 4 | 0 |
| 12.0 | ✓ | 7.37 | 1.0000 | 11 | 4 | 0 |
| 13.0 | ✓ | 6.28 | 1.0000 | 11 | 4 | 0 |
| 14.0 | ✓ | 5.41 | 1.0000 | 11 | 4 | 0 |
| 15.0 | ✓ | 4.72 | 1.0000 | 11 | 4 | 0 |

## Anomalous Observations (w-test > 3.29)

| Date | w-stat | Residual (mm) |
|---|---|---|
| 2013-03-01 | -3.75 | -7.50 |
