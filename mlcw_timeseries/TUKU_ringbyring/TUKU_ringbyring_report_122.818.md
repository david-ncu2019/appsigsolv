# Timeseries Decomposition Report: TUKU_ringbyring — 122.818

## Accepted Model

| Parameter | Value |
|---|---|
| Component | 122.818 |
| Sigma_0 assumed (mm) | 2.0 |
| Sigma_hat a-posteriori (mm) | 1.82 |
| Polynomial degree | 1 |
| Seasonal periods (yr) | [0.5, 1.0, 4.78] |
| Jump dates | ['20080520', '20210702'] |
| Polyline breaks | [] |
| Exp trend (b/day) | None |
| Exp relaxation | {'20080520': [180]} |
| Log relaxation | {} |
| n_params | 11 |
| Degrees of freedom (r) | 253 |
| T_stat (SSR/σ²) | 208.9936 |
| χ²_critical (K) | 291.1017 |
| Unit variance factor (T/r) | 0.8261 |
| p-value | 0.979952 |
| DIA iterations | 0 |

## Variance Explained per Component (%)

| Component | Std (mm) | Variance Explained (%) |
|---|---|---|
| 122.818_trend | 17.95 | 86.47 |
| 122.818_0.5yr | 0.26 | 0.02 |
| 122.818_1yr | 0.89 | 0.21 |
| 122.818_4.78yr | 0.51 | 0.07 |
| 122.818_jump | 1.42 | 0.54 |
| 122.818_exp | 3.07 | 2.53 |
| 122.818_noise | 1.78 | 0.85 |

## Sigma Scan Summary

| Sigma (mm) | Accepted | T_stat | p-value | n_params | n_periods | n_polylines |
|---|---|---|---|---|---|---|
| 2.0 | ✓ | 234.94 | 0.7990 | 10 | 3 | 0 |
| 3.0 | ✓ | 104.42 | 1.0000 | 10 | 3 | 0 |
| 4.0 | ✓ | 58.73 | 1.0000 | 10 | 3 | 0 |
| 5.0 | ✓ | 37.59 | 1.0000 | 10 | 3 | 0 |
| 6.0 | ✓ | 26.10 | 1.0000 | 10 | 3 | 0 |
| 7.0 | ✓ | 19.18 | 1.0000 | 10 | 3 | 0 |
| 8.0 | ✓ | 14.68 | 1.0000 | 10 | 3 | 0 |
| 9.0 | ✓ | 11.60 | 1.0000 | 10 | 3 | 0 |
| 10.0 | ✓ | 9.40 | 1.0000 | 10 | 3 | 0 |
| 11.0 | ✓ | 7.77 | 1.0000 | 10 | 3 | 0 |
| 12.0 | ✓ | 6.53 | 1.0000 | 10 | 3 | 0 |
| 13.0 | ✓ | 5.56 | 1.0000 | 10 | 3 | 0 |
| 14.0 | ✓ | 4.79 | 1.0000 | 10 | 3 | 0 |
| 15.0 | ✓ | 4.18 | 1.0000 | 10 | 3 | 0 |
