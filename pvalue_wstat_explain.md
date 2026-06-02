# Understanding Statistical Metrics in Geodetic Timeseries

If you are looking at a Timeseries Decomposition Report and seeing terms like **T_stat**, **chi2_critical**, **unit_var_factor**, **sigma_hat_mm**, **p-value**, and **w-stat**, this document explains what each one means and how they connect to the formal OMT/DIA statistical framework (Teunissen, TU Delft).

---

## 1. The Foundation: Sigma — "Expectation of Messiness"

Before we talk about the stats, we have to talk about **Sigma**.

Sigma ($\sigma_0$, `sigma_mm` in the report) is your **assumed a-priori measurement noise**. It represents how noisy you expect the data to be *before* fitting the model.

Think of it like this: when measuring GPS displacement, you know the receiver introduces some noise. You might say "I expect errors of about 5 mm." That 5 mm is your Sigma.

- If your data jumps by **50 mm** overnight, you know something real happened (an earthquake, equipment change, etc.).
- If your data shifts by **0.1 mm**, you don't care — it's well within your noise expectation.

**The sigma scan:** Because the true noise level is unknown, `appsigsolv` sweeps a range of sigma values (e.g., 2 mm to 20 mm in 0.5 mm steps) and finds the smallest sigma for which the model statistically passes the Overall Model Test. The selected sigma is **a-posteriori** — it was chosen *after* observing the data, not set a priori.

---

## 2. The Overall Model Test (OMT): T_stat and chi2_critical

The **Overall Model Test (OMT)** is the formal statistical test that decides whether a model fits the data well enough. It comes directly from geodetic quality control theory (Teunissen, 1988).

### The Test Statistic: T_stat

$$T = \frac{\text{SSR}}{\sigma_0^2} = \frac{\sum \hat{e}_i^2}{\sigma_0^2}$$

- **SSR** = sum of squared residuals (how far the observed points are from the fitted model)
- **σ₀** = assumed a-priori sigma
- Under a correct model with correct σ₀, $T$ follows a **chi-squared distribution** with $r$ degrees of freedom: $T \sim \chi^2(r)$

### The Critical Value: chi2_critical (K)

$$K = \chi^2_{1-\alpha}(r)$$

This is the threshold from the chi-squared table at significance level $\alpha$ (default 0.05) and $r$ degrees of freedom.

**Decision rule:**
- $T_{\text{stat}} \leq K$ → model **accepted** (residuals are consistent with assumed noise)
- $T_{\text{stat}} > K$ → model **rejected** (residuals are too large for the assumed sigma)

### The p-value

$$p = 1 - \chi^2_{\text{CDF}}(T_{\text{stat}},\ r)$$

This is equivalent to the chi-squared tail probability. The acceptance condition `p-value ≥ α` is mathematically identical to `T_stat ≤ chi2_critical`.

**In your report:** When sigma = 10.5 mm, the $T_{\text{stat}}$ was 5359.12 and the $\chi^2_{\text{critical}}$ was ~5541. Since $T \leq K$, the model was accepted, and the p-value jumped above 0.05.

---

## 3. Unit Variance Factor and A-Posteriori Sigma

These two fields tell you **how well the assumed sigma matches the actual data variability**, after the model is fitted.

### Unit Variance Factor (unit_var_factor)

$$\hat{\sigma}^2 / \sigma_0^2 = \frac{T_{\text{stat}}}{r} = \frac{\text{SSR}/r}{\sigma_0^2}$$

- **= 1.0** means the assumed sigma exactly matches the actual noise level — perfect calibration
- **< 1.0** means the model overfits, or the assumed sigma was too large (the fit is better than expected)
- **> 1.0** means the residuals are larger than the assumed noise (sigma was underestimated)

A well-accepted model typically has a unit variance factor between **0.7 and 1.1**.

### A-Posteriori Sigma (sigma_hat_mm)

$$\hat{\sigma} = \sqrt{\text{SSR}/r}$$

This is the **estimated actual noise level** from the data, in mm. Compare it to the assumed `sigma_mm`:
- If they are close, your sigma assumption was well-calibrated.
- If `sigma_hat_mm` is much smaller than `sigma_mm`, the model is somewhat over-specified.
- If `sigma_hat_mm` is larger than `sigma_mm`, the model may be missing components.

---

## 4. Degrees of Freedom (r)

$$r = n_{\text{obs}} - n_{\text{param}}$$

- $n_{\text{obs}}$ = total number of observations in the timeseries
- $n_{\text{param}}$ = number of model parameters (trend coefficients, seasonal amplitudes, jump steps, relaxation terms, etc.)

The degrees of freedom set the chi-squared distribution shape. More parameters reduce $r$, which makes it harder for the OMT to accept the model. This penalizes overly complex models — a parsimonious model is preferred.

---

## 5. The w-stat: The "Troublemaker" Test

While the p-value looks at the *entire model at once*, the **w-stat** (w-test) looks at **one specific observation** and asks: "Is this single data point an anomaly that the model cannot explain?"

### Formal Baarda w-test

The w-test is the formal Baarda identification test from the DIA framework. For a single epoch $i$, the w-statistic is:

$$w_i = \frac{\hat{e}_i}{\sigma_0 \cdot \sqrt{1 - h_{ii}}}$$

- $\hat{e}_i$ = residual at epoch $i$ (observed − model)
- $h_{ii}$ = diagonal element of the hat matrix $H = G(G^TG)^{-1}G^T$, computed efficiently via thin QR decomposition
- $\sigma_0$ = assumed a-priori sigma
- Under $H_0$ (no outlier), $w_i \sim \mathcal{N}(0, 1)$

**Why the hat matrix term?** The term $(1 - h_{ii})$ accounts for leverage — observations near the boundaries of the design matrix have higher leverage ($h_{ii}$ closer to 1) and are harder to detect as outliers. This correction ensures the w-test is statistically valid even for high-leverage points.

**Threshold:** $|w_i| > 3.29$ corresponds to a two-sided test at the 0.1% significance level. This is the standard "Anomalous Observation" threshold in geodetic quality control.

**In your report:** On `2022-06-17`, the `w-stat` was **11.46** with a residual of **131.81 mm** — this was flagged as an anomalous observation and prompted the DIA loop to add a jump at that epoch.

> **Note:** The simplified formula `w = residual / sigma` (without the hat matrix term) is an approximation that can underestimate the w-statistic for high-leverage observations. The implementation uses the full Baarda formula above.

---

## 6. How They Work Together: The DIA Loop

The software implements the **DIA (Detection, Identification, Adaptation)** framework from Teunissen (1990, TU Delft).

### Detection

At each DIA iteration, the OMT is evaluated:
- $T_{\text{stat}} > K$ → model **rejected** → proceed to Identification
- $T_{\text{stat}} \leq K$ → model **accepted** → exit loop, proceed to relaxation testing

### Identification

The Identification phase determines *which* model change best explains the current misfit. It tests four hypothesis groups using **formal w-tests**, then selects the alternative with the highest $|w|$:

1. **Datasnooping** — per-epoch unit vector test (finds single-epoch outliers)
2. **Missing periodic signal** — top Lomb-Scargle peaks tested as candidate sinusoidal columns
3. **Velocity break / polyline** — velocity spike and CUSUM candidates tested with formal w-test (CUSUM boundary from Ploberger-Krämer 1992: 1.36)
4. **Exponential trend** — AIC-detected decay rate `b` tested as an exponential column

A significance gate enforces $|w_{\text{best}}| \geq z_{\alpha/2}$ — if no alternative hypothesis passes, Identification returns "None" and the loop stops to prevent spurious adaptation.

### Adaptation

The winning alternative hypothesis is incorporated into the model:
- **Outlier** → insert a unit-step at the flagged epoch (outlier absorption)
- **Period** → add the new periodic component
- **Polyline** → add a velocity break at the candidate epoch
- **Exp trend** → add the exponential decay column

Then the OMT is re-evaluated. This cycles until the model is accepted or `max-iter` is reached.

---

## 7. Reading the Sigma Scan Summary Table

The sigma scan table in the report shows, for each candidate sigma value:

| Column | Meaning |
|---|---|
| Sigma (mm) | Assumed a-priori noise level |
| Accepted | ✓ if T_stat ≤ chi2_critical, ✗ otherwise |
| T_stat | Sum of squared residuals / σ² |
| p-value | Tail probability: 1 − χ²_CDF(T, r) |
| n_params | Number of model parameters |
| n_periods | Number of periodic components accepted |
| n_polylines | Number of polyline breaks accepted |

The **selected model** (shown in the "Accepted Model" table) is the accepted model at the smallest accepted sigma (tightest noise = best signal fit). Among models that pass at the same sigma, the one with fewest parameters (parsimony) wins; further ties are broken by highest p-value.

---

## 8. Quick Reference: All Report Fields

| Field | Formula | What a good value looks like |
|---|---|---|
| `Sigma_0 assumed (mm)` | User-specified or scan-selected σ₀ | Smallest sigma passing OMT |
| `Sigma_hat a-posteriori (mm)` | $\sqrt{\text{SSR}/r}$ | Close to `Sigma_0 assumed` |
| `Degrees of freedom (r)` | $n_{\text{obs}} - n_{\text{param}}$ | Large (many obs, few params) |
| `T_stat (SSR/σ²)` | $\text{SSR}/\sigma_0^2$ | ≤ chi2_critical |
| `χ²_critical (K)` | $\chi^2_{1-\alpha}(r)$ | Reference threshold |
| `Unit variance factor (T/r)` | $T/r$ | Near 1.0 (0.7–1.1 ideal) |
| `p-value` | $1 - \chi^2_{\text{CDF}}(T, r)$ | ≥ 0.05 (accepted) |
| `DIA iterations` | Loop count | 0 = clean fit, high = complex signal |

---

*For implementation details, see `appsigsolv/core/dia.py` — functions `calculate_omt()`, `_compute_w_stats()`, `_identify_best_alternative()`, and `run_omt_dia_loop()`.*
