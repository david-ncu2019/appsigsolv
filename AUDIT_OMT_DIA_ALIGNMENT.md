# Audit: appsigsolv Implementation vs. OMT/DIA Framework

**Date**: 2026-05-07  
**Scope**: Verification that appsigsolv correctly implements the Overall Model Test (OMT) and DIA (Detection-Identification-Adaptation) methodology from the Teunissen framework  
**Status**: ⚠️ **CRITICAL ISSUES FOUND**

---

## 1. OVERALL MODEL TEST (OMT) IMPLEMENTATION

### 1.1 Theoretical Requirements (Teunissen Framework)

The OMT is a global hypothesis test on the model:
- **Null hypothesis (H₀)**: The model is correct; residuals follow N(0, σ²)
- **Alternative hypothesis (H₁)**: The model is incorrect; residuals have systematic component
- **Test statistic**: T = (e'e) / σ₀² ~ χ²(r), where:
  - e = residuals vector
  - σ₀² = a priori variance of unit weight
  - r = redundancy = n_obs - n_param (degrees of freedom)
- **Decision**: 
  - If T ≤ χ²(α, r), **accept H₀** (model is acceptable)
  - If T > χ²(α, r), **reject H₀** (model is deficient, needs improvement)
- **p-value**: P(T > T_observed | H₀) = 1 - F(T_observed; r)

### 1.2 Code Implementation (dia.py:9-22)

```python
def calculate_omt(residuals, m_obs, n_param, sigma_m, alpha=0.05):
    r = m_obs - n_param
    if r <= 0:
        return np.inf, np.inf, 0.0, 0.0
    
    ssr = np.sum(residuals**2)
    T_stat = ssr / (sigma_m**2)
    omt = T_stat / r  # ← PROBLEM #1
    
    p_value = 1.0 - chi2.cdf(T_stat, df=r)
    K = chi2.ppf(1.0 - alpha, df=r)
    K_norm = K / r  # ← PROBLEM #2
    
    return T_stat, omt, p_value, K_norm
```

### 1.3 CRITICAL ISSUES

#### ❌ **ISSUE #1: Incorrect OMT Metric Definition**
- **Line 16**: `omt = T_stat / r`
- **Problem**: Dividing T by r gives a normalized test statistic (mean-squared residual), NOT the OMT
- **Theory**: The OMT should be T_stat itself, compared against critical value K = χ²(α, r)
- **Impact**: The code creates a misleading "OMT" value that isn't the actual test statistic
- **Consequence**: Model acceptance/rejection logic may be incorrect

**Evidence**:
```python
# Line 301-308: Uses omt for decision making
if p_value >= alpha:  # ← Correct: compares p_value to alpha
    model["_omt_stats"] = {
        ...
        "omt": omt,  # ← This "omt" is T/r, not the actual test statistic
        ...
    }
    return model
```

#### ❌ **ISSUE #2: K_norm Definition is Non-Standard**
- **Line 20**: `K_norm = K / r`
- **Problem**: K_norm is the critical value divided by DOF, which normalizes it relative to degrees of freedom
- **Theory**: This is NOT a standard geodetic metric. The critical OMT threshold is simply K = χ²(α, r)
- **Question**: What is K_norm used for? Is it being compared to the normalized omt?

#### ⚠️ **ISSUE #3: p-value Calculation Appears Correct**
- **Line 18**: `p_value = 1.0 - chi2.cdf(T_stat, df=r)` ✓
- **Correct**: This properly calculates P(χ² > T_stat) when H₀ is true
- **Decision**: Model accepted if p_value ≥ alpha (correct)

---

## 2. DIA CYCLE: DETECTION, IDENTIFICATION, ADAPTATION

### 2.1 Theoretical Requirements

The DIA method extends hypothesis testing by:
1. **Detection**: Global test (OMT) to detect if model has errors
2. **Identification**: Which parameter or component is responsible?
3. **Adaptation**: Update model and re-estimate

Expected workflow:
```
WHILE p_value < alpha (model rejected):
    1. Estimate parameters & residuals
    2. Global test (OMT) → if pass, accept
    3. Analyze residuals to identify source of rejection
    4. Adapt model (add basis function, parameter, etc.)
    5. Re-estimate and test
```

### 2.2 Code Implementation (dia.py:277-343)

```python
def run_omt_dia_loop(...):
    # Line 297-308: Initial fit and test
    for iteration in range(max_iter):
        G, m, e2, d_hat = estimate_time_func(model, date_list, dis_ts)
        residuals = dis_ts - d_hat
        T_stat, omt, p_value, K_norm = calculate_omt(residuals, len(dis_ts), n_param, sigma_m, alpha)
        
        # Line 303-308: ACCEPTANCE CONDITION
        if p_value >= alpha:  # ✓ Correct decision rule
            model["_omt_stats"] = {...}
            return model
```

✓ **Correct**: Accepts model when p_value ≥ alpha (fails to reject H₀)

```python
        # Line 310-312: Iteration termination
        if omt >= last_omt and iteration > 0:
            break
        last_omt = omt
```

⚠️ **ISSUE #4: Convergence Criterion is Non-Standard**
- **Problem**: The code exits if `omt >= last_omt` (no improvement in omt)
- **Question**: Should this be checking if omt is improving, or if we're stuck in a cycle?
- **Theory**: The standard approach is to adapt only if rejection is still occurring
- **Implementation**: The convergence check seems pragmatic but unconventional

### 2.3 Identification Phase (dia.py:314-342)

```python
adapt_type, adapt_val = robust_analyze_residuals(residuals, date_list, model)
```

✓ **Implements DIA identification**: Analyzes residuals to determine what to add
- Detects periods (periodicity detection)
- Detects polylines (velocity changes)
- Detects exp_trend (exponential trends)

---

## 3. RESIDUAL ANALYSIS & ADAPTATION (dia.py:227-276)

### 3.1 robust_analyze_residuals() Function

This function implements the "Identification" step by checking residuals for:

#### A. Periodicity Detection
```python
pgram = lombscargle(days, residuals, freqs, normalize=True)
peaks, _ = find_peaks(pgram, height=max_pwr * 0.4)
```
✓ Uses Lomb-Scargle periodogram to detect periodic components
✓ Finds peaks in the power spectrum
⚠️ **ISSUE #5**: No formal statistical test for significance
- Missing: Does the detected period contribute significantly to explaining variance?
- Missing: Cross-validation of period detection

#### B. Velocity Change Detection
```python
velocity = np.diff(residuals) / np.diff(days)
mad_vel = np.median(np.abs(valid_vel - median_vel))
if (velocity[min_idx] < median_vel - 3.5 * mad_vel) or ...
```
✓ Uses robust statistics (median + MAD) instead of mean ± std
✓ 3.5σ threshold is reasonable for outlier detection
⚠️ **ISSUE #6**: The 3.5 threshold is heuristic, not derived from OMT framework

#### C. CUSUM-Based Break Detection
```python
def _cusum_break_date(...):
    cusum = np.cumsum(residuals - mu)
    cusum_range = (np.max(cusum) - np.min(cusum)) / (sigma * np.sqrt(n))
    if cusum_range < 1.36:  # ← 95% KS critical value
        return None
```
✓ Uses CUSUM for detecting cumulative shifts
✓ References KS critical value (1.36 for ~95% confidence)
⚠️ **ISSUE #7**: Mixing CUSUM (cumulative sum) with KS (Kolmogorov-Smirnov) critical value
- The KS critical value applies to the maximum distance between empirical and theoretical CDFs
- CUSUM uses different thresholds; 1.36 may not be optimal here
- Missing citation/justification for this specific threshold

---

## 4. STOCHASTIC MODEL & SIGMA ESTIMATION

### 4.1 Theoretical Requirement

**Critical in Teunissen OMT**: σ₀² (a priori variance of unit weight) must be:
1. **Known** (from instrument specs, or)
2. **Estimated** from reference data, or
3. **Scanned** to find the value that makes the model acceptable

### 4.2 Code Implementation (dia.py:403-430)

```python
def run_omt_sigma_scan(series, ..., sigma_min, sigma_max, sigma_step, ...):
    sigmas = np.arange(sigma_min, sigma_max + sigma_step * 0.5, sigma_step)
    for sigma_mm in sigmas:
        result, table_row = _run_single_sigma(sigma_mm, ...)
```

✓ **Correct approach**: Sigma scan is a standard practice in geodesy
- Varies σ from user-specified min/max
- For each σ, runs full OMT-DIA loop
- Selects model with lowest parameters, lowest σ, highest p-value

**Question**: How is sigma_min/max chosen? 
- If too small: May never achieve acceptance (all models rejected)
- If too large: May accept overfitted models

---

## 5. MODEL SELECTION CRITERION

### 5.1 Code (dia.py:427)

```python
best = min(scan_results, key=lambda m: (
    m["_omt_stats"]["n_param"], 
    m["_omt_stats"]["sigma_mm"], 
    -m["_omt_stats"]["p_value"]
))
```

✓ **Follows Occam's Razor**: 
1. First: Minimize parameters (simplicity)
2. Then: Minimize σ (goodness-of-fit)
3. Then: Maximize p-value (statistical confidence)

✓ **Correct implementation** of standard model selection

---

## 6. POST-SEISMIC RELAXATION (exp_trend, exp, log)

### 6.1 Design Matrices

#### Exponential Trend (modeling.py:94-98)
```python
def get_design_matrix4exp_trend(date_list, b_per_day: float):
    days = np.array([(d - t0).days ...], dtype=np.float64)
    return (np.exp(-b_per_day * days) - 1.0).reshape(-1, 1)
```
✓ **Correct**: Models long-term relaxation as exp(-b*t)
✓ **Starts at zero** (exp(0) - 1 = 0): Continuous at t=0
⚠️ **Note**: This is a global model component, NOT tied to a specific event

#### Exponential Relaxation Post-Event (modeling.py:63-74)
```python
def get_design_matrix4exp_func(date_list, exp_dict):
    for exp_onset, taus in exp_dict.items():
        exp_T = datetime2years([exp_onset])[0]
        exp_tau_yr = exp_tau / 365.25
        A[:, i] = np.array(t > exp_T).flatten() * (1 - np.exp(-1 * (t - exp_T) / exp_tau_yr))
```
✓ **Correct**: Only active after event (t > exp_T)
✓ **Starts at zero** at event time: (1 - exp(0)) = 0
✓ **Asymptotes to 1** as t → ∞: Natural relaxation model

#### Logarithmic Relaxation (modeling.py:76-92)
```python
def get_design_matrix4log_func(date_list, log_dict):
    ...
    A[:, i] = np.array(t > log_T).flatten() * np.nan_to_num(
        np.log(1 + (t - log_T) / log_tau_yr),
        nan=0, neginf=0
    )
```
✓ **Correct**: Models slower, long-term relaxation
⚠️ **ISSUE #8**: Replaces NaN/neginf with 0
- Pre-event times (t < log_T): (t - log_T) < 0, so log(1 + negative) → NaN or neginf
- Code maps these to 0, which is correct
- But relying on `nan_to_num` is fragile; explicit guard would be clearer:
```python
# Better:
arg = np.where(t > log_T, (t - log_T) / log_tau_yr, 0)
A[:, i] = np.log(1 + np.maximum(arg, 0))
```

---

## 7. CRITICAL ALGORITHMIC ISSUES

### 7.1 Empty Series Guard (dia.py:278-279)

```python
series_clean = series.dropna()
if series_clean.empty:
    return None
```

✓ **Correct**: Prevents fitting on empty data
⚠️ **ISSUE #9**: No error message returned to user
- Function silently returns None
- Caller doesn't know why estimation failed
- Recommendation: Raise exception or return structured error

### 7.2 Design Matrix Singularity (modeling.py:100-110)

```python
def get_design_matrix4time_func(date_list, model, ref_date=None):
    yr_diff = np.array(datetime2years(date_list))
    if ref_date is None:
        ref_date = date_list[0]
    ref_idx = date_list.index(ref_date) if ref_date in date_list else 0
    yr_diff -= yr_diff[ref_idx]  # ← Normalizes to reference date
```

✓ **Correct**: Centers polynomial on first observation
✓ **Reduces numerical issues** in least-squares
✓ **Proper normalization**: Polynomial conditioned well

### 7.3 Condition Number Guarding (dia.py:139-148)

```python
if np.linalg.cond(G_trial) > cond_threshold:
    continue
```

✓ **Excellent practice**: Avoids adding linearly dependent parameters
✓ Uses cond_threshold=1e8 (very loose guard)
⚠️ **ISSUE #10**: The condition number threshold may be too permissive
- Standard practice: cond < 1e4 (avoid ill-conditioning)
- 1e8 allows highly ill-conditioned matrices
- Recommendation: Tighten to 1e4 or 1e5

---

## 8. STATISTICAL SIGNIFICANCE TESTING

### 8.1 Missing: Formal Tests for Identified Components

When `robust_analyze_residuals()` detects a period, it doesn't test:
- **Null hypothesis**: The detected period is noise
- **Alternative hypothesis**: The period is a real signal
- **Test**: Could use F-test comparing models with/without period

**Current approach (heuristic)**:
```python
if not any(abs(period_yr - existing) < 0.05 for existing in model.get("periodic", [])):
    return "period", period_yr
```
✓ Avoids duplicate periods
⚠️ **ISSUE #11**: No statistical significance test
- Missing: Does this period reduce residual variance significantly?
- Missing: Is it above noise power floor?

### 8.2 Velocity Change Detection: Non-Standard Threshold

```python
if (velocity[min_idx] < median_vel - 3.5 * mad_vel) or ...
```

⚠️ **ISSUE #12**: The 3.5σ threshold is heuristic
- Not derived from OMT or DIA theory
- Teunissen framework typically uses α-level tests (e.g., 0.05)
- 3.5σ ≈ 0.0005 false-alarm rate (very conservative)
- Recommend: Derive from chi-square test instead

---

## 9. INITIALIZATION & USER INPUTS

### 9.1 Jump Date Detection (dia.py:24-84)

```python
def detect_jumps(series, ..., sigma_threshold=3.0, ...):
```

✓ **Implements jump detection** with rolling median + MAD
✓ Uses adaptive thresholding (99th percentile)
⚠️ **ISSUE #13**: Mixes statistical and heuristic thresholds
- Uses σ_threshold=3.0 (statistical)
- Then uses 99th percentile (heuristic)
- These thresholds are independent; interplay unclear

### 9.2 Period Detection: No Nyquist Check

```python
def auto_detect_periods(series, min_yr=0.2, max_yr=20.0):
    ...
    freqs = np.linspace(2*np.pi/(max_yr*365.25), 2*np.pi/(min_yr*365.25), 5000)
```

⚠️ **ISSUE #14**: No Nyquist frequency validation
- If data has observations spaced > 2*min_yr apart, periods < min_yr are aliased
- Example: Annual data (365 days apart) can only resolve periods ≥ 2 years
- Code doesn't check observation spacing
- Missing: Warning if time span < period × 3

---

## 10. SUMMARY TABLE OF ISSUES

| # | Issue | Severity | Impact | Location |
|---|-------|----------|--------|----------|
| 1 | Incorrect OMT = T/r definition | 🔴 CRITICAL | Misleading metric; decision logic unclear | dia.py:16 |
| 2 | K_norm = K/r non-standard | 🟡 MEDIUM | Unclear usage; non-standard metric | dia.py:20 |
| 3 | p-value calculation | ✓ OK | Correct | dia.py:18 |
| 4 | Convergence criterion (omt >= last_omt) | 🟡 MEDIUM | Non-standard but pragmatic | dia.py:310 |
| 5 | Period detection: No significance test | 🟡 MEDIUM | May add noise periods | dia.py:233-241 |
| 6 | Velocity threshold 3.5σ heuristic | 🟡 MEDIUM | Not from OMT framework | dia.py:256 |
| 7 | CUSUM uses KS critical value | 🟠 HIGH | Mixing incompatible tests | dia.py:164 |
| 8 | Log func uses nan_to_num fragile | 🟡 MEDIUM | Works but unclear | modeling.py:86-88 |
| 9 | Empty series returns None silently | 🟡 MEDIUM | No error feedback | dia.py:279 |
| 10 | Condition number 1e8 permissive | 🟡 MEDIUM | Allows ill-conditioning | dia.py:147 |
| 11 | No F-test for period significance | 🟡 MEDIUM | May overfit | dia.py:241 |
| 12 | 3.5σ velocity threshold heuristic | 🟡 MEDIUM | Non-standard decision rule | dia.py:256 |
| 13 | Jump detection mixes σ & percentile | 🟡 MEDIUM | Threshold logic unclear | dia.py:39-43 |
| 14 | No Nyquist frequency validation | 🟡 MEDIUM | May detect aliases | dia.py:94-96 |

---

## 11. RECOMMENDATIONS

### Immediate Fixes (Correctness)
1. **Fix Issue #1**: Change `omt = T_stat / r` to document that this is "normalized OMT" or use T_stat directly
2. **Fix Issue #7**: Replace KS critical value with proper CUSUM threshold or cite source
3. **Fix Issue #8**: Replace `nan_to_num` with explicit guards

### Short-term Improvements (Clarity)
4. Add docstrings explaining non-standard choices (K_norm, convergence criterion)
5. Add statistical tests for period/polyline detection
6. Tighten condition number threshold to 1e4 or document why 1e8 is acceptable

### Long-term (Theory Alignment)
7. Implement F-tests for component significance
8. Add Nyquist frequency validation
9. Replace heuristic thresholds with α-level tests from Teunissen framework
10. Add structured error returns instead of silent failures

---

## 12. WHAT IS "FAKING" AND WHAT IS NOT

### 12.1 Not Faking (Correct Implementation)

These components genuinely follow Teunissen methodology:

1. **P-value calculation** ✓
   - Line 18: `p_value = 1.0 - chi2.cdf(T_stat, df=r)`
   - This is the correct test statistic CDF

2. **DIA loop structure** ✓
   - Detection: Global OMT test (lines 301)
   - Identification: robust_analyze_residuals() (line 314)
   - Adaptation: Model updated based on residuals (lines 316-340)
   - Iteration: Repeats until acceptance (line 303)

3. **Design matrices** ✓
   - Polynomial, periodic, step, polyline, exponential, logarithmic all correctly implemented
   - Proper bases for all basis functions

4. **Model selection** ✓
   - Uses Occam's Razor: fewer parameters preferred
   - Sigma scan is standard practice

### 12.2 Is Faking (Problematic)

These ARE concerning and should be fixed:

1. **The "omt" variable is misleading** 🔴
   ```python
   omt = T_stat / r  # This is NOT the Overall Model Test!
   ```
   - Should be called something like `mean_squared_residual` or `test_statistic_normalized`
   - The actual OMT is just `T_stat`
   - This creates confusion for users reading the code or outputs

2. **CUSUM threshold mixes incompatible tests** 🔴
   ```python
   if cusum_range < 1.36:  # This is a KS critical value, not CUSUM threshold
   ```
   - This is technically wrong; should have proper CUSUM threshold
   - However, 1.36 might work pragmatically by coincidence

3. **Heuristic thresholds without justification** 🟡
   - 3.5σ for velocity outliers
   - 0.4 * max_pwr for period peaks
   - These work but lack theoretical foundation

### 12.3 Not Faking (Just Unclear Documentation)

1. **K_norm = K/r** — Non-standard but not wrong, just unclear purpose
2. **Convergence criterion (omt >= last_omt)** — Non-standard but pragmatic
3. **Condition number 1e8** — Loose but documented choice

---

## 13. PRACTICAL VERIFICATION

To verify if the code "fakes" results, test these scenarios:

### Test 1: Pure Noise (No Real Signal)
```
Input: White noise series with known σ
Expected: OMT test should reject for small σ, accept for large σ
If code fakes: Would accept unrealistically small σ or add spurious components
```

### Test 2: Known Signal + Noise
```
Input: Synthetic series: [constant] + [annual period] + [linear trend] + [noise]
Expected: Model should identify exactly these components
If code fakes: Would over-fit (add fake periods) or under-fit (miss real ones)
```

### Test 3: Jump Detection
```
Input: Series with known jump at specific date
Expected: detect_jumps() should find jump with ±few days accuracy
If code fakes: Would miss jumps or detect false ones
```

---

## 14. CONCLUSION

**Overall Assessment**: 🟡 **PARTIALLY ALIGNED WITH THEORY**

### What Works ✓
- P-value calculation and hypothesis testing are correct
- DIA cycle (detection-identification-adaptation) properly implemented
- Model selection logic is sound (Occam's Razor + sigma scan)
- Design matrices are mathematically correct
- Component extraction and reporting work as intended

### What Needs Fixing 🔴
1. **Rename `omt` variable** to clarify it's normalized T/r, not the actual test statistic
2. **Fix CUSUM threshold** to use proper CUSUM-specific critical value
3. **Document heuristic choices** (3.5σ, 0.4*max_pwr, etc.) with justification or cite sources

### What's Misleading 🟡
1. Calls the metric "omt" when it's T_stat/r (normalized form)
2. Uses thresholds that feel statistical but are actually heuristic
3. Missing significance tests for detected components

### Is It "Faking"?
**No, it's not deliberately fraudulent.** But it does obscure some non-standard choices:
- The code works mathematically
- But the terminology is confusing ("omt" should not be T/r)
- Some design choices lack proper statistical justification
- Risk: Users may think it strictly follows Teunissen when it doesn't entirely

**Recommendation**: Fix the terminology, add comments explaining non-standard choices, and add formal statistical tests for component identification.


