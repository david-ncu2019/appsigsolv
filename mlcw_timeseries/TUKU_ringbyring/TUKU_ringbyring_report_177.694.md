# Timeseries Decomposition Report: TUKU_ringbyring — 177.694

## Accepted Model

| Parameter | Value |
|---|---|
| Component | 177.694 |
| Sigma_0 assumed (mm) | 3.0 |
| Sigma_hat a-posteriori (mm) | 2.34 |
| Polynomial degree | 0 |
| Seasonal periods (yr) | [0.5, 1.0, 5.96, 12.87] |
| Jump dates | ['20050422', '20110113', '20120221'] |
| Polyline breaks | [] |
| Exp trend (b/day) | None |
| Exp relaxation | {'20050422': [180], '20110113': [180], '20120221': [180]} |
| Log relaxation | {} |
| n_params | 15 |
| Degrees of freedom (r) | 247 |
| T_stat (SSR/σ²) | 150.6791 |
| χ²_critical (K) | 284.6599 |
| Unit variance factor (T/r) | 0.6100 |
| p-value | 1.000000 |
| DIA iterations | 0 |

## Variance Explained per Component (%)

| Component | Std (mm) | Variance Explained (%) |
|---|---|---|
| 177.694_trend | 0.00 | 0.00 |
| 177.694_0.5yr | 0.23 | 0.06 |
| 177.694_1yr | 0.27 | 0.08 |
| 177.694_5.96yr | 1.19 | 1.52 |
| 177.694_12.87yr | 2.90 | 9.05 |
| 177.694_jump | 3.80 | 15.60 |
| 177.694_exp | 12.26 | 162.27 |
| 177.694_noise | 2.28 | 5.59 |

## Sigma Scan Summary

| Sigma (mm) | Accepted | T_stat | p-value | n_params | n_periods | n_polylines |
|---|---|---|---|---|---|---|
| 2.0 | ✓ | 255.25 | 0.2974 | 18 | 6 | 2 |
| 3.0 | ✓ | 243.71 | 0.6002 | 12 | 4 | 0 |
| 4.0 | ✓ | 137.08 | 1.0000 | 12 | 4 | 0 |
| 5.0 | ✓ | 87.73 | 1.0000 | 12 | 4 | 0 |
| 6.0 | ✓ | 60.93 | 1.0000 | 12 | 4 | 0 |
| 7.0 | ✓ | 44.76 | 1.0000 | 12 | 4 | 0 |
| 8.0 | ✓ | 34.27 | 1.0000 | 12 | 4 | 0 |
| 9.0 | ✓ | 27.08 | 1.0000 | 12 | 4 | 0 |
| 10.0 | ✓ | 21.93 | 1.0000 | 12 | 4 | 0 |
| 11.0 | ✓ | 18.13 | 1.0000 | 12 | 4 | 0 |
| 12.0 | ✓ | 15.23 | 1.0000 | 12 | 4 | 0 |
| 13.0 | ✓ | 12.98 | 1.0000 | 12 | 4 | 0 |
| 14.0 | ✓ | 11.19 | 1.0000 | 12 | 4 | 0 |
| 15.0 | ✓ | 9.75 | 1.0000 | 12 | 4 | 0 |
