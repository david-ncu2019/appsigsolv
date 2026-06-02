# Timeseries Decomposition Report: TUKU_ringbyring — 252.265

## Accepted Model

| Parameter | Value |
|---|---|
| Component | 252.265 |
| Sigma_0 assumed (mm) | 3.0 |
| Sigma_hat a-posteriori (mm) | 2.74 |
| Polynomial degree | 2 |
| Seasonal periods (yr) | [0.5, 1.0, 8.24, 16.69] |
| Jump dates | ['20140423', '20160501', '20211203'] |
| Polyline breaks | [] |
| Exp trend (b/day) | None |
| Exp relaxation | {'20211203': [30]} |
| Log relaxation | {} |
| n_params | 15 |
| Degrees of freedom (r) | 249 |
| T_stat (SSR/σ²) | 207.3090 |
| χ²_critical (K) | 286.8078 |
| Unit variance factor (T/r) | 0.8326 |
| p-value | 0.974632 |
| DIA iterations | 0 |

## Variance Explained per Component (%)

| Component | Std (mm) | Variance Explained (%) |
|---|---|---|
| 252.265_trend | 40.99 | 111.71 |
| 252.265_0.5yr | 0.05 | 0.00 |
| 252.265_1yr | 0.70 | 0.03 |
| 252.265_8.24yr | 3.09 | 0.64 |
| 252.265_16.69yr | 3.33 | 0.74 |
| 252.265_jump | 3.51 | 0.82 |
| 252.265_exp | 5.55 | 2.05 |
| 252.265_noise | 2.66 | 0.47 |

## Sigma Scan Summary

| Sigma (mm) | Accepted | T_stat | p-value | n_params | n_periods | n_polylines |
|---|---|---|---|---|---|---|
| 2.0 | ✗ | — | — | — | — | — |
| 3.0 | ✓ | 210.62 | 0.9666 | 14 | 4 | 0 |
| 4.0 | ✓ | 118.47 | 1.0000 | 14 | 4 | 0 |
| 5.0 | ✓ | 75.82 | 1.0000 | 14 | 4 | 0 |
| 6.0 | ✓ | 52.66 | 1.0000 | 14 | 4 | 0 |
| 7.0 | ✓ | 38.69 | 1.0000 | 14 | 4 | 0 |
| 8.0 | ✓ | 29.62 | 1.0000 | 14 | 4 | 0 |
| 9.0 | ✓ | 23.40 | 1.0000 | 14 | 4 | 0 |
| 10.0 | ✓ | 18.96 | 1.0000 | 14 | 4 | 0 |
| 11.0 | ✓ | 15.67 | 1.0000 | 14 | 4 | 0 |
| 12.0 | ✓ | 13.16 | 1.0000 | 14 | 4 | 0 |
| 13.0 | ✓ | 11.22 | 1.0000 | 14 | 4 | 0 |
| 14.0 | ✓ | 9.67 | 1.0000 | 14 | 4 | 0 |
| 15.0 | ✓ | 8.42 | 1.0000 | 14 | 4 | 0 |
