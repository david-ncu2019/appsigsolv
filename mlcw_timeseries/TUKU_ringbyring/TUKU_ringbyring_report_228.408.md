# Timeseries Decomposition Report: TUKU_ringbyring — 228.408

## Accepted Model

| Parameter | Value |
|---|---|
| Component | 228.408 |
| Sigma_0 assumed (mm) | 2.0 |
| Sigma_hat a-posteriori (mm) | 1.95 |
| Polynomial degree | 1 |
| Seasonal periods (yr) | [0.5, 1.0, 5.82, 16.16] |
| Jump dates | ['20060223'] |
| Polyline breaks | [] |
| Exp trend (b/day) | None |
| Exp relaxation | {'20060223': [180]} |
| Log relaxation | {} |
| n_params | 12 |
| Degrees of freedom (r) | 251 |
| T_stat (SSR/σ²) | 238.2596 |
| χ²_critical (K) | 288.9551 |
| Unit variance factor (T/r) | 0.9492 |
| p-value | 0.708381 |
| DIA iterations | 0 |

## Variance Explained per Component (%)

| Component | Std (mm) | Variance Explained (%) |
|---|---|---|
| 228.408_trend | 6.54 | 59.75 |
| 228.408_0.5yr | 0.05 | 0.00 |
| 228.408_1yr | 0.21 | 0.06 |
| 228.408_5.82yr | 1.08 | 1.63 |
| 228.408_16.16yr | 1.84 | 4.71 |
| 228.408_jump | 0.30 | 0.13 |
| 228.408_exp | 1.47 | 3.04 |
| 228.408_noise | 1.90 | 5.07 |

## Sigma Scan Summary

| Sigma (mm) | Accepted | T_stat | p-value | n_params | n_periods | n_polylines |
|---|---|---|---|---|---|---|
| 2.0 | ✓ | 246.20 | 0.5911 | 11 | 4 | 0 |
| 3.0 | ✓ | 109.42 | 1.0000 | 11 | 4 | 0 |
| 4.0 | ✓ | 61.55 | 1.0000 | 11 | 4 | 0 |
| 5.0 | ✓ | 39.39 | 1.0000 | 11 | 4 | 0 |
| 6.0 | ✓ | 27.36 | 1.0000 | 11 | 4 | 0 |
| 7.0 | ✓ | 20.10 | 1.0000 | 11 | 4 | 0 |
| 8.0 | ✓ | 15.39 | 1.0000 | 11 | 4 | 0 |
| 9.0 | ✓ | 12.16 | 1.0000 | 11 | 4 | 0 |
| 10.0 | ✓ | 9.85 | 1.0000 | 11 | 4 | 0 |
| 11.0 | ✓ | 8.14 | 1.0000 | 11 | 4 | 0 |
| 12.0 | ✓ | 6.84 | 1.0000 | 11 | 4 | 0 |
| 13.0 | ✓ | 5.83 | 1.0000 | 11 | 4 | 0 |
| 14.0 | ✓ | 5.02 | 1.0000 | 11 | 4 | 0 |
| 15.0 | ✓ | 4.38 | 1.0000 | 11 | 4 | 0 |

## Anomalous Observations (w-test > 3.29)

| Date | w-stat | Residual (mm) |
|---|---|---|
| 2006-02-23 | -6.43 | -12.86 |
