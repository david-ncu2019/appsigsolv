# appsigsolv User Guide

**appsigsolv** — Applied Signal Solver for geodetic and geophysical timeseries decomposition.

This guide covers installation, both CLI commands, input/output formats, interpretation of results, and worked examples.

---

## Table of Contents

1. [What appsigsolv Does](#1-what-appsigsolv-does)
2. [Installation & Requirements](#2-installation--requirements)
3. [Quick Start](#3-quick-start)
4. [Input Data Formats](#4-input-data-formats)
5. [Command: `decompose`](#5-command-decompose)
6. [Command: `reconstruct`](#6-command-reconstruct)
7. [Understanding the Outputs](#7-understanding-the-outputs) *(7.1 JSON · 7.2 CSV · 7.3 Report fields)*
8. [Statistical Concepts (Plain Language)](#8-statistical-concepts-plain-language) *(OMT · T_stat · chi2_critical · unit_var_factor · sigma_hat · w-test · DIA)*
9. [Worked Examples](#9-worked-examples)
10. [Tips, Gotchas & FAQ](#10-tips-gotchas--faq)
11. [Changelog](#11-changelog)

---

## 1. What appsigsolv Does

`appsigsolv` fits a **parametric signal model** to a geodetic or geophysical timeseries and decomposes it into physically meaningful components:

- **Long-term trend** — polynomial (offset, linear velocity, acceleration)
- **Exponential decay trend** — monotonic decay shape `a·(exp(-b·t)−1)` for compaction or relaxation signals that asymptote toward a plateau
- **Seasonal signals** — sinusoidal waves at user-specified or auto-detected periods
- **Jumps** — instantaneous step changes from earthquakes, equipment changes, etc.
- **Polyline breaks** — changes in velocity rate (piecewise linear)
- **Post-event relaxation** — exponential or logarithmic transients after a jump

The model is selected automatically using the **DIA (Detection, Identification, Adaptation)** framework combined with the **Overall Model Test (OMT)** — a rigorous statistical procedure from geodetic quality control theory (Teunissen, TU Delft). The tool sweeps a range of assumed noise levels (sigma scan) and selects the most parsimonious accepted model.

**Supported data types:**
- GPS/GNSS displacement timeseries (NEU components)
- InSAR displacement timeseries
- Groundwater level timeseries
- Multilayer compaction well (MLCW) layer-by-layer measurements

---

## 2. Installation & Requirements

### Dependencies

`appsigsolv` uses the standard scientific Python stack:

```
numpy
scipy
pandas
matplotlib
```

Install with pip if needed:

```bash
pip install numpy scipy pandas matplotlib
```

### Running the Package

There is no `setup.py`. Run the package directly from the parent directory of `appsigsolv/`:

```bash
cd D:\1000_SCRIPTS\004_Project003\20260501_timeseries_signal_solver

python -m appsigsolv --help
```

All examples in this guide assume this is your working directory.

---

## 3. Quick Start

**Decompose a GPS Up component with defaults:**

```bash
python -m appsigsolv decompose gps_timeseries/TKJS_neu.csv --component dU
```

**Decompose all components in a file (Batch Mode):**

```bash
python -m appsigsolv decompose gps_timeseries/TKJS_neu.csv --component all
```

**Re-process a specific component, overwriting an existing result:**

```bash
python -m appsigsolv decompose gps_timeseries/TKJS_neu.csv --component dU --force
```

**Reconstruct at specific days of every month:**

```bash
python -m appsigsolv reconstruct gps_timeseries/TKJS_neu.csv \
    --json gps_timeseries/TKJS_neu/TKJS_neu_model_dU.json \
    --target-col dU --sampling-rate custom --custom-dates 1,6,11,16,21,26
```

---

## 4. Input Data Formats

### 4.1 GPS / Regular Timeseries (CSV)

The tool expects a CSV with:
- One **date column** (auto-detected by name: `gpsdate`, `datetime`, `date`, `Date`, `time`, `Time`)
- One or more **displacement columns** in millimetres (default unit)

```csv
datetime,gpsdate,dN,dE,dU,sN,sE,sU
2010.00548,2010-01-02,-39.99,86.92,-103.59,...
2010.01096,2010-01-06,-39.37,85.40,-99.89,...
```

The date column must be parseable by `pandas.to_datetime` (ISO format `YYYY-MM-DD` recommended).

The tool resamples to a **daily grid** by default and fills short gaps (up to 7 days) by interpolation. If your data is already daily and gap-free, the resampling is a no-op.

### 4.2 MLCW / Irregular Wide-Format (CSV)

Multi-layer compaction well data is typically wide-format with irregular observation times:

```csv
datetime,8.775,11.938,25.605,39.545,...
2003-12-03,0.0,0.0,0.0,...
2004-01-09,4.0,-2.0,-3.0,...
2004-02-15,7.0,-4.0,-5.0,...
```

- First column: date
- Remaining columns: layer depths (numeric column names, values in mm)
- Use `--irregular` flag to skip daily resampling (essential for this format)

**Numeric column names** (e.g. `86.15899999999999`) are automatically normalised to 3 decimal places (`86.159`) in all output filenames, JSON keys, and CSV headers to avoid floating-point precision artefacts in file names.

### 4.3 Units

- Default input unit: `mm` (millimetres)
- Internal computation unit: `m` (metres) — all output CSVs are back-scaled to `mm`
- Use `--unit m` if your input is already in metres

---

## 5. Command: `decompose`

The `decompose` command runs the full automated pipeline: load data → detect jumps → scan sigma → DIA loop → save results.

### Syntax

```bash
python -m appsigsolv decompose <input_csv> [OPTIONS]
```

### All Options

| Option | Default | Description |
|---|---|---|
| `input_csv` | *(required)* | Path to input CSV |
| `--component` | `dU` | Column name(s) to process. Use comma-separated list (`dN,dE,dU`) or `all` |
| `--date-col` | *(auto)* | Override the auto-detected date column name |
| `--unit` | `mm` | Input unit: `mm` or `m` |
| `--jumps` | | Extra jump dates to force: `YYYY-MM-DD,YYYY-MM-DD,...` |
| `--polylines` | | Extra polyline break dates: `YYYY-MM-DD,YYYY-MM-DD,...` |
| `--logs` | | Log relaxation terms: `YYYY-MM-DD:tau_days,...` |
| `--poly-deg` | `1` | Polynomial degree: `0`=offset only, `1`=linear trend, `2`=acceleration, `3`=cubic. Use `-1` to auto-select the degree with the tightest passing sigma (0–3). |
| `--periods` | `0.25,0.5,1.0,2.0` | Candidate periods (years) always included in the model search |
| `--auto-periods` | `5` | Auto-detect up to N additional dominant periods via Lomb-Scargle |
| `--sigma-min` | `2.0` | Start of sigma scan in mm |
| `--sigma-max` | `15.0` | End of sigma scan in mm |
| `--sigma-step` | `0.5` | Sigma step size in mm |
| `--alpha` | `0.05` | OMT significance level (model accepted when p-value ≥ alpha) |
| `--max-iter` | `5` | Max DIA iterations per sigma value |
| `--irregular` | *(flag)* | Skip daily resampling — required for MLCW or non-daily data |
| `--no-plot` | *(flag)* | Skip PNG figure generation |
| `--no-relax` | *(flag)* | Skip exponential/log relaxation testing after model acceptance |
| `--force` | *(flag)* | **Overwrite existing results.** By default, components with an existing `_model_<comp>.json` are skipped. Use `--force` to re-process them (e.g. when you want to try different parameters on a component). |
| `--output-dir` | *(same folder as CSV)* | Parent directory for outputs |
| `--cores` | `1` | Number of CPU cores for parallel sigma scan |
| `--exp-trend` | *(off)* | Exponential decay trend `exp(-b·t)−1`. Pass `auto` to auto-detect the best decay rate via AIC, or a numeric `b` in 1/days (e.g. `0.001` ≈ 2.7-yr time constant). Omit to disable. |

### Resumption and Overwrite Behaviour

By default, `decompose` **skips** any component that already has a `{stem}_model_{comp}.json` in the output directory. This means:

- **Interrupted batch runs** can be safely re-started — already-processed components are skipped automatically.
- **Manual re-processing** requires `--force` to overwrite an existing result.

```bash
# Skip already-done components (default behaviour)
python -m appsigsolv decompose station.csv --component all ...

# Overwrite a specific unsatisfactory result
python -m appsigsolv decompose station.csv --component 86.159 --force --poly-deg 2 ...

# Re-run all components from scratch
python -m appsigsolv decompose station.csv --component all --force ...
```

### Outputs

Results are saved in `<output_dir>/<csv_stem>/` (e.g., `gps_timeseries/TKJS_neu/`):

| File | Description |
|---|---|
| `<stem>_model_<comp>.json` | Accepted model configuration (reusable with `reconstruct`) |
| `<stem>_decomposed_<comp>.csv` | All decomposed signal components (displacement in mm) |
| `<stem>_decomposed_<comp>.png` | Diagnostic figure |
| `<stem>_report_<comp>.md` | Statistical summary report |

---

## 6. Command: `reconstruct`

The `reconstruct` command fits a **user-specified** model (either from a JSON file or CLI flags) to a timeseries and saves the fitted curve. 

### Syntax

```bash
python -m appsigsolv reconstruct <input_file> [OPTIONS]
```

### All Options

| Option | Default | Description |
|---|---|---|
| `input_file` | *(required)* | CSV or Excel file |
| `--json` | | Load model from a JSON config file (produced by `decompose`) |
| `--date-col` | *(auto)* | Override auto-detected date column |
| `--target-col` | *(auto)* | Override auto-detected displacement column |
| `--unit` | `mm` | Input unit: `mm` or `m` |
| `--poly` | | Polynomial degree |
| `--period` | | Periodic component in years (repeatable) |
| `--step` | | Jump step date `YYYYMMDD` (repeatable) |
| `--polyline` | | Polyline break date `YYYYMMDD` (repeatable) |
| `--exp DATE TAU` | | Exponential term: onset date `YYYYMMDD` and tau in days (repeatable) |
| `--log DATE TAU` | | Logarithmic term: onset date `YYYYMMDD` and tau in days (repeatable) |
| `--exp-trend B` | *(off)* | Exponential decay trend `b` value in 1/days. Loaded automatically when using `--json` from `decompose`. |
| `-o` / `--output` | *(auto)* | Output file path |
| `--sampling-rate` | *(none)* | Set to `daily` for every day, or `custom` for specific days of each month |
| `--custom-dates` | *(none)* | Comma-separated day numbers (1–31) used with `--sampling-rate custom` |
| `--ref-date` | *(first date)* | Reference epoch for model `YYYYMMDD` |

### Custom Sampling

When using `--sampling-rate custom`, you must provide `--custom-dates`. For example, `--custom-dates 1,15` will produce modeled values for the 1st and 15th of every month within the timeseries range. Day values exceeding the length of a specific month (e.g., 31 in February) are automatically clamped to the last day of that month.

### Output

A single CSV (or Excel) containing:
- `date` column
- `modeled` column (if sampling at original observation times)
- `reconstructed` column (if using `daily` or `custom` sampling)

The output filename will automatically include a suffix like `_daily`, `_custom`, or `_modeled` if `--output` is not specified.

---

## 7. Understanding the Outputs

### 7.1 Model JSON

```json
{
    "polynomial": 1,
    "periodic": [0.5, 1.0, 3.2, 4.75, 9.14, 2.0],
    "stepDate": ["20220617"],
    "polyline": [],
    "exp": {},
    "log": {},
    "exp_trend": null
}
```

If an exponential decay trend was detected, `exp_trend` will hold the fitted `b` value (in 1/days) instead of `null`, for example `"exp_trend": 0.000616`. This value is automatically used by `reconstruct --json`.

This file is directly usable as `--json` input to `reconstruct`.

### 7.2 Decomposed CSV

All displacement values are in **millimetres** (mm). The `_wtest` and `flagged` columns are dimensionless. Each row is one date.

| Column | Description |
|---|---|
| `date` | ISO date |
| `<comp>` | Observed displacement |
| `<comp>_trend` | Polynomial trend component |
| `<comp>_exp_trend` | Exponential decay trend `a·(exp(-b·t)−1)` (only present when detected or specified) |
| `<comp>_<T>yr` | Seasonal component at period T years |
| `<comp>_jump` | Cumulative step function (if jumps present) |
| `<comp>_exp_<date>` | Exponential relaxation term (if present) |
| `<comp>_log_<date>` | Logarithmic relaxation term (if present) |
| `<comp>_model` | Sum of all model components (fitted values) |
| `<comp>_noise` | Residual = observed − model |
| `<comp>_wtest` | Baarda w-statistic per observation: `ê_i / (σ · sqrt(1 − h_ii))` where `h_ii` is the hat-matrix diagonal |
| `flagged` | `True` when `|w-test| > 3.29` (anomalous observation at 0.1% level) |

### 7.3 Report Fields (Accepted Model Table)

The Markdown report includes an **Accepted Model** table. Here is what each field means:

| Field | Explanation |
|---|---|
| `Sigma_0 assumed (mm)` | The a-priori noise level selected by the sigma scan. This is the smallest sigma at which the OMT accepted the model. |
| `Sigma_hat a-posteriori (mm)` | Estimated actual noise level from the data: `sqrt(SSR / r)`. Should be close to `Sigma_0 assumed`. |
| `Polynomial degree` | Trend degree: 0 = offset, 1 = linear velocity, 2 = acceleration, 3 = cubic. With `--poly-deg -1`, the degree whose model passes OMT at the smallest sigma is selected; ties broken by fewest parameters. |
| `Seasonal periods (yr)` | All periodic components in the accepted model, in years. |
| `Jump dates` | Auto-detected or user-forced step-change epochs. |
| `Exp relaxation` | Post-jump exponential relaxation terms (onset → tau in days). |
| `Log relaxation` | Post-jump logarithmic relaxation terms. |
| `Exp trend (b/day)` | Exponential decay trend rate (null if not used). |
| `n_params` | Total number of fitted model parameters. Fewer is better (parsimony). |
| `Degrees of freedom (r)` | `n_obs − n_params`. Determines the chi-squared distribution shape for the OMT. |
| `T_stat (SSR/σ²)` | Overall Model Test statistic: sum of squared residuals divided by σ₀². Accepted when `T_stat ≤ χ²_critical`. |
| `χ²_critical (K)` | Chi-squared critical value at significance level α: `χ²_{1−α}(r)`. The formal threshold for model acceptance. |
| `Unit variance factor (T/r)` | `T_stat / r = σ̂² / σ₀²`. Near 1.0 means well-calibrated; < 1 means slight overfit; > 1 means underfit. |
| `p-value` | Tail probability `1 − χ²_CDF(T, r)`. Accepted when `p-value ≥ α`. Equivalent to `T_stat ≤ χ²_critical`. |
| `DIA iterations` | Number of Detection–Identification–Adaptation cycles needed before acceptance. |

---

## 8. Statistical Concepts (Plain Language)

### Sigma (σ₀) — "Expectation of Messiness"

Sigma (`sigma_mm` in the report) is your **assumed a-priori measurement noise** — how noisy you expect the data to be before fitting any model. The sigma scan sweeps a range (e.g., 2–20 mm) and selects the smallest sigma at which the Overall Model Test passes. The selected sigma is therefore **a-posteriori**: it was determined by observing the data, not set independently in advance.

### Overall Model Test — T_stat, chi2_critical, and p-value

The **Overall Model Test (OMT)** is the formal chi-squared test that decides whether a model fits. The test statistic is:

```
T_stat = SSR / σ₀²
```

where SSR is the sum of squared residuals. Under a correctly specified model, `T_stat ~ χ²(r)` (chi-squared with `r` degrees of freedom).

**Acceptance criterion:** `T_stat ≤ χ²_critical` — equivalently, `p-value ≥ α` (default α = 0.05).

- **T_stat ≤ χ²_critical** → model **accepted** (residuals are consistent with assumed noise)
- **T_stat > χ²_critical** → model **rejected** (residuals are too large)

The p-value is just the tail probability `1 − χ²_CDF(T_stat, r)`, so checking `p ≥ 0.05` is mathematically identical to checking `T_stat ≤ K`.

### Unit Variance Factor and A-Posteriori Sigma

Two fields in the report characterise how well the assumed sigma matches the actual data variability after fitting:

**Unit variance factor** = `T_stat / r = σ̂² / σ₀²`
- Near 1.0 → sigma assumption was well-calibrated
- < 1.0 → model slightly over-specified, or sigma was too large
- > 1.0 → residuals larger than expected; sigma may be underestimated

**Sigma_hat a-posteriori** = `sqrt(SSR / r)` in mm — the estimated actual noise level from the data. Compare to `Sigma_0 assumed`: close agreement means good calibration.

### Degrees of Freedom (r)

```
r = n_obs − n_params
```

Degrees of freedom control the chi-squared distribution shape used for the OMT. Adding model parameters (more seasonal components, jump terms, relaxation terms) reduces `r`, which raises the `χ²_critical` threshold. The tool prefers **parsimonious models** — when multiple models pass at the same sigma, the one with fewest parameters (parsimony) is preferred; further ties are broken by p-value.

### w-statistic — "Is This Specific Point an Outlier?"

The w-statistic implements the **formal Baarda datasnooping test** per observation epoch:

```
w_i = ê_i / (σ₀ · sqrt(1 − h_ii))
```

- `ê_i` = residual at epoch i
- `h_ii` = hat-matrix diagonal element (leverage), computed via thin QR decomposition
- Under H₀ (no outlier), `w_i ~ N(0, 1)`

If `|w_i| > 3.29`, the observation is flagged as an **anomalous observation** (two-sided test at 0.1% significance level). The CSV `flagged` column marks these rows.

> The simplified formula `w = residual / sigma` (without the `sqrt(1 − h_ii)` correction) is an approximation. The implementation uses the full Baarda formula, which correctly accounts for the leverage of each observation.

### DIA Loop — Detection, Identification, Adaptation

The sigma scan runs a DIA loop at each candidate sigma value:

1. **Detection** — evaluate the OMT. If accepted, exit.
2. **Identification** — test four alternative hypothesis groups using formal w-tests:
   - Per-epoch datasnooping (single outlier)
   - Missing periodic signals (Lomb-Scargle candidates)
   - Velocity break / polyline (CUSUM + velocity spike candidates)
   - Exponential trend (AIC-detected decay rate)
   
   The alternative with the highest `|w|` above `z_{α/2}` wins. If none passes, adaptation stops.
3. **Adaptation** — incorporate the winning alternative into the model, then repeat from Detection.

This continues until the OMT accepts or `--max-iter` is reached.

### Exponential Decay Trend — "What is it?"

Some signals do not increase or decrease at a steady rate — they start fast and then slow down, asymptoting toward a final plateau. Groundwater compaction and post-pumping recovery are typical examples.

The model used is:

```
displacement(t) = a · (exp(-b · t) − 1)
```

- `a` — the amplitude (total displacement at infinite time, in metres)
- `b` — the decay rate in 1/days. Larger `b` means faster convergence. The time constant τ = 1/b.
- At `t = 0`: displacement = 0 (anchored to the first observation)

Because `b` is fixed before OLS fitting (amplitude `a` is solved linearly), this component integrates cleanly with the OMT statistical framework. When using `--exp-trend auto`, the tool scans 20 candidate `b` values and selects the one with the lowest AIC over a plain linear model. The threshold is **ΔAIC > 10** (strong evidence criterion) to avoid adding an exponential trend for weak signals.

---

## 9. Worked Examples

### Example 1 — Batch Processing with Resumption

```bash
python -m appsigsolv decompose gps_timeseries/TKJS_neu.csv --component all
```

Processes all columns. If interrupted, re-running this command will skip already completed components automatically.

### Example 2 — Re-processing an Unsatisfactory Component

After reviewing the diagnostic plot, you decide a component needs a different polynomial degree:

```bash
python -m appsigsolv decompose station.csv --component 86.159 \
    --poly-deg 2 --force --irregular --no-relax \
    --periods 0.5,1 --output-dir ./MLCW_decomposition
```

`--force` deletes any existing JSON/CSV/PNG/report for `86.159` and re-processes from scratch.

### Example 3 — Custom Reconstruction for Comparison

After `decompose` produces a model, reconstruct it on exactly the 1st and 15th of every month:

```bash
python -m appsigsolv reconstruct gps_timeseries/TKJS_neu.csv \
    --json gps_timeseries/TKJS_neu/TKJS_neu_model_dU.json \
    --target-col dU \
    --sampling-rate custom \
    --custom-dates 1,15 \
    -o results/TKJS_semi_monthly.csv
```

### Example 4 — MLCW Column with Exponential Decay Trend

For compaction-well data where a layer shows a monotonic exponential-decay shape (e.g., rapid early compaction tapering off over years), use `--exp-trend auto` to let the tool find the best decay rate:

```bash
python -m appsigsolv decompose XIGANG_ringbyring.csv \
    --component 102.096 \
    --date-col datetime \
    --unit mm \
    --poly-deg -1 \
    --periods 0.5,1 \
    --auto-periods 5 \
    --sigma-min 1.0 \
    --sigma-max 30.0 \
    --sigma-step 0.5 \
    --irregular \
    --no-relax \
    --exp-trend auto \
    --output-dir ./MLCW_decomposition
```

The exponential trend component (`<comp>_exp_trend`) will appear in the decomposed CSV and as a dashed purple line in the trend panel of the diagnostic plot. The fitted `b` value is recorded in the JSON and the report.

To then reconstruct at specific dates using the accepted model:

```bash
python -m appsigsolv reconstruct XIGANG_ringbyring.csv \
    --json MLCW_decomposition/XIGANG_ringbyring/XIGANG_ringbyring_model_102.096.json \
    --target-col 102.096 \
    --sampling-rate custom \
    --custom-dates 1,15
```

The `exp_trend` key in the JSON is picked up automatically — no extra flags needed.

---

## 10. Tips, Gotchas & FAQ

### When should I use `--irregular`?

Use `--irregular` whenever your observations are **not daily** (e.g., MLCW or monthly InSAR). It prevents the tool from creating excessive synthetic points via resampling.

### How does automatic polynomial selection work?

If you pass `--poly-deg -1`, the tool will test degrees 0, 1, 2, and 3. It selects the degree whose accepted model passes OMT at the **smallest sigma** (tighter noise = better signal fit). Among degrees that pass at the same sigma, the one with fewer parameters wins; further ties are broken by p-value.

### How do I re-run a component I am not happy with?

Use `--force` with the specific `--component` name and any adjusted parameters. Without `--force`, the existing JSON will be detected and skipped.

### Output CSV units vs. input units

Input is mm (or m if `--unit m`). **All output CSV columns are in millimetres (mm).**

### Why does my batch script stop silently mid-run?

On Windows, running many stations in a single Python process can exhaust GDI handles as matplotlib accumulates figure objects. Best practice for batch scripts:

1. Use `cores=1` — avoids `ProcessPoolExecutor` spawn issues on Windows.
2. Call `plt.close('all')` and `gc.collect()` between stations in the `finally` block.
3. Import `matplotlib` and call `matplotlib.use("Agg")` at the top of the batch script before importing `appsigsolv`.

### Why do MLCW column names have ugly decimal filenames?

Float column names in CSVs (e.g. `86.15899999999999`) are a pandas read artefact. `appsigsolv` automatically normalises all numeric column names to 3 decimal places (`86.159`) during loading, so all output filenames, JSON keys, and CSV headers use the clean form.

---

## 11. Changelog

### v0.3.0 (2026-05)

**New features:**
- `--force` flag for `decompose`: overwrite existing results instead of skipping. Useful when manually re-processing unsatisfactory components with adjusted parameters.

**Bug fixes & robustness:**
- **Batch crash fix (Windows GDI exhaustion):** `save_plot` now calls `plt.close(fig)` + `plt.close('all')` to guarantee figure handles are released after every plot. Batch scripts should additionally call `plt.close('all')` + `gc.collect()` between stations.
- **Numeric column name normalisation:** All floating-point CSV column name artefacts (e.g. `86.15899999999999`) are now rounded to 3 decimal places at load time in both `load_and_preprocess` and the `component="all"` path of `decompose`. Output filenames are clean (e.g. `86.159`).
- **Per-component exception handling:** The entire fitting + extraction + output pipeline for each component is now wrapped in `try/except Exception`, so a single failing component logs its error and continues — it no longer crashes the entire batch or station run.
- **Adaptive plot x-axis:** `YearLocator` interval is now data-driven (1 yr for spans < 4 yr, 2 yr for 4–10 yr, 5 yr for > 10 yr), preventing a matplotlib crash on short timeseries.

**Algorithm:**
- **Sigma-first cross-degree selection:** When `--poly-deg -1`, the winning degree is now chosen by the smallest accepted sigma (tightest noise), not best p-value. Ties are broken by parameter count, then p-value.
- **`cores` default changed to `1`:** Avoids Windows `ProcessPoolExecutor` / `spawn` issues when calling `run_decompose` directly from a batch script (not under `if __name__ == '__main__'`).

### v0.2.0 (prior)

- Exponential decay trend component (`--exp-trend auto` / numeric b)
- Auto polynomial degree selection (`--poly-deg -1`)
- Sigma scan with OMT-based model acceptance
- DIA loop with Lomb-Scargle period detection, jump detection, relaxation testing
- `reconstruct` command with custom sampling rates

---

*Package version: 0.3.0 | appsigsolv — Applied Signal Solver*
