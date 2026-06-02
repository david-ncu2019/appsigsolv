# Timeseries Decomposition Report: TUKU_ringbyring — 86.914

## Accepted Model

| Parameter | Value |
|---|---|
| Component | 86.914 |
| Sigma_0 assumed (mm) | 2.0 |
| Sigma_hat a-posteriori (mm) | 1.90 |
| Polynomial degree | 2 |
| Seasonal periods (yr) | [0.5, 1.0, 12.87] |
| Jump dates | ['20060223', '20090915', '20140617', '20210407'] |
| Polyline breaks | [] |
| Exp trend (b/day) | None |
| Exp relaxation | {'20140617': [180], '20210407': [90]} |
| Log relaxation | {} |
| n_params | 15 |
| Degrees of freedom (r) | 249 |
| T_stat (SSR/σ²) | 225.5298 |
| χ²_critical (K) | 286.8078 |
| Unit variance factor (T/r) | 0.9057 |
| p-value | 0.854652 |
| DIA iterations | 0 |

## Variance Explained per Component (%)

| Component | Std (mm) | Variance Explained (%) |
|---|---|---|
| 86.914_trend | 12.98 | 60.28 |
| 86.914_0.5yr | 0.85 | 0.26 |
| 86.914_1yr | 1.16 | 0.48 |
| 86.914_12.87yr | 1.37 | 0.67 |
| 86.914_jump | 0.92 | 0.30 |
| 86.914_exp | 4.29 | 6.58 |
| 86.914_noise | 1.85 | 1.22 |

## Sigma Scan Summary

| Sigma (mm) | Accepted | T_stat | p-value | n_params | n_periods | n_polylines |
|---|---|---|---|---|---|---|
| 2.0 | ✓ | 241.91 | 0.6484 | 13 | 3 | 0 |
| 3.0 | ✓ | 107.51 | 1.0000 | 13 | 3 | 0 |
| 4.0 | ✓ | 60.48 | 1.0000 | 13 | 3 | 0 |
| 5.0 | ✓ | 38.71 | 1.0000 | 13 | 3 | 0 |
| 6.0 | ✓ | 26.88 | 1.0000 | 13 | 3 | 0 |
| 7.0 | ✓ | 19.75 | 1.0000 | 13 | 3 | 0 |
| 8.0 | ✓ | 15.12 | 1.0000 | 13 | 3 | 0 |
| 9.0 | ✓ | 11.95 | 1.0000 | 13 | 3 | 0 |
| 10.0 | ✓ | 9.68 | 1.0000 | 13 | 3 | 0 |
| 11.0 | ✓ | 8.00 | 1.0000 | 13 | 3 | 0 |
| 12.0 | ✓ | 6.72 | 1.0000 | 13 | 3 | 0 |
| 13.0 | ✓ | 5.73 | 1.0000 | 13 | 3 | 0 |
| 14.0 | ✓ | 4.94 | 1.0000 | 13 | 3 | 0 |
| 15.0 | ✓ | 4.30 | 1.0000 | 13 | 3 | 0 |
