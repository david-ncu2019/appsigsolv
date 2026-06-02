# Fix Outer-Ring Decomposition Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix four root causes that make the DIA/OMT decomposition algorithm overfit, produce collinear long-period sinusoids, miss structural breaks, and report wtest as a signal component, for large-displacement outer-ring InSAR columns.

**Architecture:** All changes are confined to two files: `appsigsolv/core/dia.py` (sigma selection key, collinearity guard, CUSUM-based structural break) and `appsigsolv/utils/visualization.py` (exclude wtest from variance table). No CLI arguments, no CSV format, no batch script changes.

**Tech Stack:** Python, NumPy, SciPy (`linalg.lstsq`, `chi2`, `lombscargle`), pandas

---

## File Map

| File | Change |
|---|---|
| `appsigsolv/core/dia.py` | Fix 1: sigma selection key in `run_omt_sigma_scan`; Fix 2: add CUSUM break detection in `robust_analyze_residuals`; Fix 3: add collinearity guard in `prescreen_periods` |
| `appsigsolv/utils/visualization.py` | Fix 4: skip wtest column in variance table in `save_report` |

Absolute paths:
- `D:\1000_SCRIPTS\004_Project003\20260501_timeseries_signal_solver\appsigsolv\core\dia.py`
- `D:\1000_SCRIPTS\004_Project003\20260501_timeseries_signal_solver\appsigsolv\utils\visualization.py`

---

## Task 1: Fix sigma selection key (parsimony-first, largest acceptable sigma)

**Files:**
- Modify: `appsigsolv/core/dia.py`, line 314

**Why:** The current key `(sigma_mm, n_param, -p_value)` selects the tightest-noise model first, which incentivises the DIA loop to overfit (add many spurious periods) to squeeze residuals below a tight sigma threshold. The correct priority is: fewest parameters (parsimony), then largest acceptable sigma (most conservative noise assumption), then highest p-value.

- [ ] **Step 1: Change the selection key in `run_omt_sigma_scan`**

In `appsigsolv/core/dia.py`, find line 314:
```python
    best = min(scan_results, key=lambda m: (m["_omt_stats"]["sigma_mm"], m["_omt_stats"]["n_param"], -m["_omt_stats"]["p_value"]))
```

Replace with:
```python
    best = min(scan_results, key=lambda m: (m["_omt_stats"]["n_param"], -m["_omt_stats"]["sigma_mm"], -m["_omt_stats"]["p_value"]))
```

- [ ] **Step 2: Verify syntax**

Run:
```powershell
python -c "import ast; ast.parse(open(r'D:\1000_SCRIPTS\004_Project003\20260501_timeseries_signal_solver\appsigsolv\core\dia.py').read()); print('OK')"
```
Expected: `OK`

- [ ] **Step 3: Commit**

```powershell
git -C "D:\1000_SCRIPTS\004_Project003\20260501_timeseries_signal_solver" add appsigsolv/core/dia.py
git -C "D:\1000_SCRIPTS\004_Project003\20260501_timeseries_signal_solver" commit -m "fix: prefer parsimonious model at largest acceptable sigma in sigma scan"
```

---

## Task 2: Add collinearity guard in `prescreen_periods`

**Files:**
- Modify: `appsigsolv/core/dia.py`, function `prescreen_periods` (~lines 108–133)

**Why:** Periods longer than ~40% of the observation span (e.g., 13-15 yr against a 24.5-yr record) produce sine/cosine columns in the design matrix that are nearly linearly dependent on the polynomial trend. This causes the lstsq solver to produce component magnitudes that partially cancel, with trend explaining >100% of variance. The fix rejects any candidate period whose addition raises the design matrix condition number above a threshold.

- [ ] **Step 1: Add the collinearity guard inside `prescreen_periods`**

Replace the entire `prescreen_periods` function in `appsigsolv/core/dia.py`:

```python
def prescreen_periods(series: pd.Series, candidates_yr: list, cond_threshold: float = 1e8) -> list:
    from appsigsolv.core.modeling import build_design_matrix_cols
    values = series.values
    valid_mask = ~np.isnan(values)
    days = np.array([(d - series.index[0]).days for d in series.index])
    if valid_mask.sum() < 20:
        return []

    coeffs = np.polyfit(days[valid_mask], values[valid_mask], 2)
    detrended = values[valid_mask] - np.polyval(coeffs, days[valid_mask])

    freqs = np.linspace(2*np.pi/(20*365.25), 2*np.pi/(0.2*365.25), 10000)
    pgram = lombscargle(days[valid_mask], detrended, freqs, normalize=True)

    accepted = []
    power_threshold = np.percentile(pgram, 90)
    
    # Build a baseline design matrix: quadratic polynomial only
    # Columns: [1, t, t^2] evaluated at valid observation days
    t = days[valid_mask].astype(float)
    t_max = t.max() if t.max() > 0 else 1.0
    G_base = np.column_stack([np.ones(len(t)), t / t_max, (t / t_max)**2])

    for period_yr in candidates_yr:
        period_days = period_yr * 365.25
        target_f = 2*np.pi / period_days
        idx = np.argmin(np.abs(freqs - target_f))
        window = pgram[max(0, idx-50):min(len(pgram), idx+51)]
        has_peak = len(window) > 0 and window.max() > power_threshold and window.max() > 0.05

        if not has_peak:
            continue

        # Collinearity guard: check condition number with this period added
        cos_col = np.cos(2 * np.pi / period_days * t)
        sin_col = np.sin(2 * np.pi / period_days * t)
        # Append already-accepted periods too to simulate the full G
        G_trial = np.column_stack([G_base] + [
            col for p in accepted
            for col in (np.cos(2*np.pi/(p*365.25)*t), np.sin(2*np.pi/(p*365.25)*t))
        ] + [cos_col, sin_col])
        cond = np.linalg.cond(G_trial)
        if cond > cond_threshold:
            continue  # reject: too collinear with existing columns

        accepted.append(period_yr)
    return accepted
```

- [ ] **Step 2: Verify syntax**

```powershell
python -c "import ast; ast.parse(open(r'D:\1000_SCRIPTS\004_Project003\20260501_timeseries_signal_solver\appsigsolv\core\dia.py').read()); print('OK')"
```
Expected: `OK`

- [ ] **Step 3: Verify the import of `build_design_matrix_cols` is not needed (we inline the matrix construction)**

The code above builds G_base directly from `days` without importing anything extra — remove the unused import line from the function body:

The function as written does NOT call `build_design_matrix_cols`. Delete that import line from inside `prescreen_periods` so the function starts with:
```python
def prescreen_periods(series: pd.Series, candidates_yr: list, cond_threshold: float = 1e8) -> list:
    values = series.values
    ...
```

- [ ] **Step 4: Quick smoke test**

```powershell
cd "D:\1000_SCRIPTS\004_Project003\20260501_timeseries_signal_solver"
python -c "
import pandas as pd, numpy as np
from appsigsolv.core.dia import prescreen_periods
idx = pd.date_range('1997-01-01', periods=298, freq='30D')
s = pd.Series(np.sin(2*np.pi/365.25*np.arange(298)*30), index=idx)
result = prescreen_periods(s, [0.5, 1.0, 15.0])
print('Accepted periods:', result)
print('15yr rejected (collinear):', 15.0 not in result)
"
```
Expected output: `15.0 not in result` should print `True` (15yr is nearly collinear with the polynomial over a 24-yr span).

- [ ] **Step 5: Commit**

```powershell
git -C "D:\1000_SCRIPTS\004_Project003\20260501_timeseries_signal_solver" add appsigsolv/core/dia.py
git -C "D:\1000_SCRIPTS\004_Project003\20260501_timeseries_signal_solver" commit -m "fix: reject long periods collinear with polynomial trend in prescreen_periods"
```

---

## Task 3: Add CUSUM-based structural break detection in `robust_analyze_residuals`

**Files:**
- Modify: `appsigsolv/core/dia.py`, function `robust_analyze_residuals` (~lines 135–170)

**Why:** The 2018-2021 acceleration in the outer-ring columns is a sustained, gradual trend change, not a velocity spike. The existing inter-epoch velocity check only triggers on a single extreme epoch. A CUSUM (cumulative sum) test detects where a sustained shift in residual level begins, and returns that date as a polyline break candidate — even when individual inter-epoch velocities don't exceed the MAD threshold.

- [ ] **Step 1: Add CUSUM break detection as a third branch in `robust_analyze_residuals`**

Replace the entire `robust_analyze_residuals` function:

```python
def robust_analyze_residuals(residuals: np.ndarray, dates: list, model: dict):
    days = np.array([(d - dates[0]).days for d in dates])
    freqs = np.linspace(2*np.pi/(20*365.25), 2*np.pi/(0.2*365.25), 5000)
    pgram = lombscargle(days, residuals, freqs, normalize=True)
    
    max_pwr = np.max(pgram)
    peaks, _ = find_peaks(pgram, height=max_pwr * 0.4)
    
    if len(peaks) > 0 and max_pwr > 0.1:
        peak_powers = pgram[peaks]
        sorted_peaks = peaks[np.argsort(peak_powers)[::-1]]
        for p_idx in sorted_peaks:
            f_rad = freqs[p_idx]
            period_yr = round((2*np.pi / f_rad) / 365.25, 2)
            if not any(abs(period_yr - existing) < 0.05 for existing in model.get("periodic", [])):
                return "period", period_yr

    velocity = np.diff(residuals) / np.diff(days)
    dt = np.diff(days)
    valid_vel_mask = dt < np.percentile(dt, 95) if len(dt) > 10 else np.ones_like(dt, dtype=bool)
    
    if valid_vel_mask.sum() > 5:
        valid_vel = velocity[valid_vel_mask]
        median_vel = np.median(valid_vel)
        mad_vel = np.median(np.abs(valid_vel - median_vel))
        
        if mad_vel > 1e-12:
            min_idx = np.argmin(velocity)
            max_idx = np.argmax(velocity)
            if (velocity[min_idx] < median_vel - 3.5 * mad_vel) or (velocity[max_idx] > median_vel + 3.5 * mad_vel):
                extreme_idx = min_idx if abs(velocity[min_idx]) > abs(velocity[max_idx]) else max_idx
                break_date = dates[extreme_idx].strftime("%Y%m%d")
                if break_date not in model.get("polyline", []):
                    return "polyline", break_date

    # CUSUM-based structural break detection for gradual accelerations
    # Detects a sustained level/trend shift not caught by the velocity spike test
    if len(residuals) >= 20:
        cusum_break = _cusum_break_date(residuals, dates)
        if cusum_break is not None and cusum_break not in model.get("polyline", []):
            return "polyline", cusum_break

    return "None", None


def _cusum_break_date(residuals: np.ndarray, dates: list, min_segment: int = 10) -> str | None:
    """Return the date string of the maximum-CUSUM structural break, or None if not significant."""
    n = len(residuals)
    if n < 2 * min_segment:
        return None
    mu = np.mean(residuals)
    cusum = np.cumsum(residuals - mu)
    # The CUSUM test statistic: max |S_k| / (sigma * sqrt(n))
    sigma = np.std(residuals, ddof=1)
    if sigma < 1e-12:
        return None
    # Normalised CUSUM range
    cusum_range = (np.max(cusum) - np.min(cusum)) / (sigma * np.sqrt(n))
    # Significance threshold at ~95% level for the Kolmogorov-Smirnov statistic is ~1.36
    if cusum_range < 1.36:
        return None
    # Break point = index of maximum absolute deviation of CUSUM from zero
    break_idx = int(np.argmax(np.abs(cusum)))
    # Enforce minimum segment length on both sides
    if break_idx < min_segment or break_idx > n - min_segment:
        return None
    return dates[break_idx].strftime("%Y%m%d")
```

- [ ] **Step 2: Verify syntax**

```powershell
python -c "import ast; ast.parse(open(r'D:\1000_SCRIPTS\004_Project003\20260501_timeseries_signal_solver\appsigsolv\core\dia.py').read()); print('OK')"
```
Expected: `OK`

- [ ] **Step 3: Smoke test the CUSUM function**

```powershell
cd "D:\1000_SCRIPTS\004_Project003\20260501_timeseries_signal_solver"
python -c "
import numpy as np, pandas as pd
from datetime import datetime, timedelta
from appsigsolv.core.dia import _cusum_break_date

# Simulate a structural break at index 150 (gradual acceleration)
np.random.seed(0)
n = 200
dates = [datetime(1997,1,1) + timedelta(days=30*i) for i in range(n)]
residuals = np.concatenate([np.random.normal(0, 1, 150), np.random.normal(-5, 1, 50)])
result = _cusum_break_date(residuals, dates)
print('Break detected at:', result)
print('Expected: around 1997-01-01 + 150*30 days =', dates[150].strftime('%Y%m%d'))

# No-break case
flat = np.random.normal(0, 1, 200)
print('No-break result (should be None):', _cusum_break_date(flat, dates))
"
```
Expected: break detected near index 150, None for flat noise.

- [ ] **Step 4: Commit**

```powershell
git -C "D:\1000_SCRIPTS\004_Project003\20260501_timeseries_signal_solver" add appsigsolv/core/dia.py
git -C "D:\1000_SCRIPTS\004_Project003\20260501_timeseries_signal_solver" commit -m "fix: add CUSUM structural break detection for gradual accelerations in residuals"
```

---

## Task 4: Exclude wtest from variance-explained table in report

**Files:**
- Modify: `appsigsolv/utils/visualization.py`, function `save_report` (~lines 143–148)

**Why:** `{comp}_wtest` is stored in the `components` dict (it's the normalised residual = noise/sigma_m, a dimensionless diagnostic). The variance loop includes it, producing physically meaningless entries like "44 million % variance explained". It should be excluded alongside `{comp}_model`.

- [ ] **Step 1: Add wtest to the skip list in `save_report`**

In `appsigsolv/utils/visualization.py`, find line 144:
```python
        if col in (f"{comp}_model", f"{comp}_noise"):
```

Replace with:
```python
        if col in (f"{comp}_model", f"{comp}_noise", f"{comp}_wtest"):
```

- [ ] **Step 2: Verify syntax**

```powershell
python -c "import ast; ast.parse(open(r'D:\1000_SCRIPTS\004_Project003\20260501_timeseries_signal_solver\appsigsolv\utils\visualization.py').read()); print('OK')"
```
Expected: `OK`

- [ ] **Step 3: Commit**

```powershell
git -C "D:\1000_SCRIPTS\004_Project003\20260501_timeseries_signal_solver" add appsigsolv/utils/visualization.py
git -C "D:\1000_SCRIPTS\004_Project003\20260501_timeseries_signal_solver" commit -m "fix: exclude wtest column from variance-explained table in report"
```

---

## Task 5: End-to-end verification on XIGANG outer-ring columns

**Goal:** Confirm all four fixes improve decomposition quality on the problem columns.

- [ ] **Step 1: Delete existing output for two problem columns so they are re-processed**

```powershell
Remove-Item "D:\1000_SCRIPTS\004_Project003\20260427_InSAR_MLCW_v2\MLCW_decomposition\XIGANG_ringbyring\XIGANG_ringbyring_model_69.45.json"
Remove-Item "D:\1000_SCRIPTS\004_Project003\20260427_InSAR_MLCW_v2\MLCW_decomposition\XIGANG_ringbyring\XIGANG_ringbyring_model_82.066.json"
Remove-Item "D:\1000_SCRIPTS\004_Project003\20260427_InSAR_MLCW_v2\MLCW_decomposition\XIGANG_ringbyring\XIGANG_ringbyring_report_69.45.md"
Remove-Item "D:\1000_SCRIPTS\004_Project003\20260427_InSAR_MLCW_v2\MLCW_decomposition\XIGANG_ringbyring\XIGANG_ringbyring_report_82.066.md"
Remove-Item "D:\1000_SCRIPTS\004_Project003\20260427_InSAR_MLCW_v2\MLCW_decomposition\XIGANG_ringbyring\XIGANG_ringbyring_decomposed_69.45.csv"
Remove-Item "D:\1000_SCRIPTS\004_Project003\20260427_InSAR_MLCW_v2\MLCW_decomposition\XIGANG_ringbyring\XIGANG_ringbyring_decomposed_82.066.csv"
```

- [ ] **Step 2: Run decomposition on the two problem columns**

```powershell
$env:PYTHONPATH = "D:\1000_SCRIPTS\004_Project003\20260501_timeseries_signal_solver"
cd "D:\1000_SCRIPTS\004_Project003\20260501_timeseries_signal_solver"
python -m appsigsolv decompose `
  "D:\1000_SCRIPTS\004_Project003\20260427_InSAR_MLCW_v2\MLCW_timeseries\XIGANG_ringbyring.csv" `
  --component "69.45,82.066" `
  --date-col datetime `
  --unit mm `
  --poly-deg -1 `
  --periods "0.5,1" `
  --auto-periods 5 `
  --sigma-min 1.0 `
  --sigma-max 30.0 `
  --sigma-step 0.5 `
  --alpha 0.05 `
  --max-iter 1000 `
  --irregular `
  --no-plot `
  --no-relax `
  --output-dir "D:\1000_SCRIPTS\004_Project003\20260427_InSAR_MLCW_v2\MLCW_decomposition" `
  --cores 1
```

- [ ] **Step 3: Check the new reports — verify improvements**

```powershell
Get-Content "D:\1000_SCRIPTS\004_Project003\20260427_InSAR_MLCW_v2\MLCW_decomposition\XIGANG_ringbyring\XIGANG_ringbyring_report_69.45.md"
Get-Content "D:\1000_SCRIPTS\004_Project003\20260427_InSAR_MLCW_v2\MLCW_decomposition\XIGANG_ringbyring\XIGANG_ringbyring_report_82.066.md"
```

**Pass criteria (all must hold):**
1. The `wtest` row is absent from the Variance Explained table
2. `noise variance%` is < 100% for both columns (was 229% and 373%)
3. `trend variance%` is < 100% for 82.066 (was 162%)
4. The model JSON has fewer periods than before (69.45 had 4, 82.066 had 6 — both should drop)
5. At least one of the two columns has a `polyline` entry in the model JSON (the CUSUM break from ~2018-2021)

```powershell
Get-Content "D:\1000_SCRIPTS\004_Project003\20260427_InSAR_MLCW_v2\MLCW_decomposition\XIGANG_ringbyring\XIGANG_ringbyring_model_69.45.json"
Get-Content "D:\1000_SCRIPTS\004_Project003\20260427_InSAR_MLCW_v2\MLCW_decomposition\XIGANG_ringbyring\XIGANG_ringbyring_model_82.066.json"
```

- [ ] **Step 4: If all pass criteria hold, run the full batch to update all 12 problem columns**

Delete all existing outputs for the problem columns then run `batch_process_MLCW.py`:
```powershell
# Delete all problem column outputs to force reprocessing
$cols = "51.038","60.042","69.45","82.066","91.988","102.096","122.911","134.076","142.674","152.077","172.045","182.089"
$dir = "D:\1000_SCRIPTS\004_Project003\20260427_InSAR_MLCW_v2\MLCW_decomposition\XIGANG_ringbyring"
foreach ($c in $cols) {
    Remove-Item "$dir\XIGANG_ringbyring_model_$c.json" -ErrorAction SilentlyContinue
    Remove-Item "$dir\XIGANG_ringbyring_report_$c.md" -ErrorAction SilentlyContinue
    Remove-Item "$dir\XIGANG_ringbyring_decomposed_$c.csv" -ErrorAction SilentlyContinue
}
python "D:\1000_SCRIPTS\004_Project003\20260427_InSAR_MLCW_v2\batch_process_MLCW.py"
```
