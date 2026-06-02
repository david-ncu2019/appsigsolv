# Timeseries Decomposition Report: TUKU_ringbyring — 220.335

## Accepted Model

| Parameter | Value |
|---|---|
| Component | 220.335 |
| Sigma_0 assumed (mm) | 2.0 |
| Sigma_hat a-posteriori (mm) | 1.87 |
| Polynomial degree | 0 |
| Seasonal periods (yr) | [0.5, 1.0, 6.99, 12.39] |
| Jump dates | ['20140617'] |
| Polyline breaks | [] |
| Exp trend (b/day) | None |
| Exp relaxation | {'20140617': [90]} |
| Log relaxation | {} |
| n_params | 11 |
| Degrees of freedom (r) | 250 |
| T_stat (SSR/σ²) | 218.1986 |
| χ²_critical (K) | 287.8815 |
| Unit variance factor (T/r) | 0.8728 |
| p-value | 0.927479 |
| DIA iterations | 0 |

## Variance Explained per Component (%)

| Component | Std (mm) | Variance Explained (%) |
|---|---|---|
| 220.335_trend | 0.00 | 0.00 |
| 220.335_0.5yr | 0.10 | 0.04 |
| 220.335_1yr | 0.41 | 0.75 |
| 220.335_6.99yr | 0.53 | 1.25 |
| 220.335_12.39yr | 2.18 | 21.09 |
| 220.335_jump | 1.23 | 6.71 |
| 220.335_exp | 2.70 | 32.29 |
| 220.335_noise | 1.83 | 14.81 |

## Sigma Scan Summary

| Sigma (mm) | Accepted | T_stat | p-value | n_params | n_periods | n_polylines |
|---|---|---|---|---|---|---|
| 2.0 | ✓ | 225.03 | 0.8793 | 10 | 4 | 0 |
| 3.0 | ✓ | 100.01 | 1.0000 | 10 | 4 | 0 |
| 4.0 | ✓ | 56.26 | 1.0000 | 10 | 4 | 0 |
| 5.0 | ✓ | 36.01 | 1.0000 | 10 | 4 | 0 |
| 6.0 | ✓ | 25.00 | 1.0000 | 10 | 4 | 0 |
| 7.0 | ✓ | 18.37 | 1.0000 | 10 | 4 | 0 |
| 8.0 | ✓ | 14.06 | 1.0000 | 10 | 4 | 0 |
| 9.0 | ✓ | 11.11 | 1.0000 | 10 | 4 | 0 |
| 10.0 | ✓ | 9.00 | 1.0000 | 10 | 4 | 0 |
| 11.0 | ✓ | 7.44 | 1.0000 | 10 | 4 | 0 |
| 12.0 | ✓ | 6.25 | 1.0000 | 10 | 4 | 0 |
| 13.0 | ✓ | 5.33 | 1.0000 | 10 | 4 | 0 |
| 14.0 | ✓ | 4.59 | 1.0000 | 10 | 4 | 0 |
| 15.0 | ✓ | 4.00 | 1.0000 | 10 | 4 | 0 |

## Anomalous Observations (w-test > 3.29)

| Date | w-stat | Residual (mm) |
|---|---|---|
| 2014-06-17 | -5.08 | -10.15 |
