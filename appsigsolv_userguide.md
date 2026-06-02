# appsigsolv User Guide

**appsigsolv** — Applied Signal Solver for geodetic and geophysical timeseries decomposition.

This guide covers installation, both CLI commands, input/output formats, interpretation of results, and worked examples.

---

## Table of Contents

1. [What appsigsolv Does](#1-what-appsigolv-does)
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

The `decompose` command runs the full automated pipeline: load data -> detect jumps -> scan sigma -> DIA loop -> save results.

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
| `--poly-deg` | `1` | Polynomial degree: `0`=offset only, `1`=linear trend, `2`=acceleration, `3`=cubic. Use `-1` to auto-select the degree with the tightest passing sigma (from `--poly-deg-min` up to 3). |
| `--poly-deg-min` | `0` | Minimum polynomial degree considered when `--poly-deg -1`. Set to `1` to exclude offset-only models, `2` to require at least acceleration. Has no effect when `--poly-deg` is a fixed value. |
| `--periods` | `0.25,0.5,1.0,2.0` | Candidate periods (years) always included in the model search |
| `--auto-periods` | `5` | Auto-detect up to N additional dominant periods via Lomb-Scargle |
| `--sigma-min` | `2.0` | Start of sigma scan in mm |
| `--sigma-max` | `15.0` | End of sigma scan in mm |
| `--sigma-step` | `0.5` | Sigma step size in mm |
| `--alpha` | `0.05` | OMT significance level (model accepted when p-value >= alpha) |
| `--max-iter` | `5` | Max DIA iterations per sigma value |
| `--irregular` | *(flag)* | Skip daily resampling — required for MLCW or non-daily data |
| `--no-plot` | *(flag)* | Skip PNG figure generation |
| `--no-relax` | *(flag)* | Skip exponential/log relaxation testing after model acceptance |
| `--no-jump` | *(flag)* | Disable automatic jump detection. Manually specified `--jumps` dates are still applied. Suppresses hypothesis group 1 (datasnooping/outlier) in the DIA loop. |
| `--no-seasonal` | *(flag)* | Disable all periodic/seasonal component detection. Overrides `--periods` and `--auto-periods`. Suppresses hypothesis group 2 (missing periodic) in the DIA loop. |
| `--no-exp-trend` | *(flag)* | Disable exponential trend auto-detection. Suppresses hypothesis group 4 (exponential trend) in the DIA loop. Unlike omitting `--exp-trend` (which allows automatic detection), this actively prevents detection. |
| `--auto-sigma` | *(flag)* | Select the sigma value that **maximises polyline breakpoints** rather than minimising parameter count. Designed for pure piecewise-linear fitting where the goal is the finest velocity segmentation that still passes the OMT. |
| `--start-date` | *(none)* | Filter data to start at this date (inclusive). Format: `YYYY-MM-DD`. Applied after preprocessing, before jump detection. |
| `--end-date` | *(none)* | Filter data to end at this date (inclusive). Format: `YYYY-MM-DD`. |
| `--exp-trend` | *(off)* | Exponential decay trend `exp(-b.t)-1`. Pass `auto` to auto-detect the best decay rate via AIC, or a numeric `b` in 1/days (e.g. `0.001` ~= 2.7-yr time constant). Omit to allow the DIA loop to detect it as hypothesis group 4. |
| `--output-dir` | *(same folder as CSV)* | Parent directory for outputs |
| `--cores` | `1` | Number of CPU cores for parallel sigma scan |
| `--force` | *(flag)* | **Overwrite existing results.** By default, components with an existing `_model_<comp>.json` are skipped. Use `--force` to re-process them. |

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
| `<stem>_skipped_<comp>.txt` | Written when a component is auto-skipped (see below) |

**JSON note:** The saved model JSON strips internal control fields. Specifically, `no_seasonal`, `no_jump`, and `allowed_periods` are **not written** to the JSON file — they are runtime-only flags that do not belong in the persisted model configuration.

### Auto-Skip Behaviour

When a component cannot be fit, `decompose` writes a `_skipped_<comp>.txt` file and moves on to the next component. This happens in two cases:

- **No accepted model** — the sigma scan exhausted all candidate sigma values and polynomial degrees without the OMT accepting any model. This typically indicates a malfunctioned sensor, a highly irregular signal, or a data span too short for the selected periods. The skip report records the degrees tried and the sigma range.
- **Timeout** — fitting exceeded **180 seconds** per component. The fitting thread is abandoned and a skip report is written with the reason `"timed out"`. This protects batch runs against pathologically slow timeseries on Windows (where `signal.SIGALRM` is unavailable).

In both cases the batch continues without interruption. Inspect the skip report and the diagnostic plot of the raw series to decide whether to adjust parameters or exclude the station entirely.

```
# Example skip report content:
Component : dU
Reason    : No accepted model found after 57s (tried polynomial degrees [1, 2, 3], sigma 2.0–20.0 mm). Likely a malfunctioned or highly irregular timeseries.
```

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

**Fields not persisted in JSON:** Three runtime-only fields are used inside the DIA loop but stripped before saving: `no_seasonal`, `no_jump`, and `allowed_periods`. These control how the DIA searches for candidate components and have no meaning outside the fitting process. The saved JSON contains only the accepted model structure, which is directly usable as `--json` input to `reconstruct`.

### 7.2 Decomposed CSV

All displacement values are in **millimetres** (mm). The `_wtest` and `flagged` columns are dimensionless. Each row is one date.

| Column | Description |
|---|---|
| `date` | ISO date |
| `<comp>` | Observed displacement |
| `<comp>_trend` | Polynomial trend component (includes polyline breaks) |
| `<comp>_exp_trend` | Exponential decay trend `a.(exp(-b.t)-1)` (only present when detected or specified) |
| `<comp>_<T>yr` | Seasonal component at period T years |
| `<comp>_jump` | Cumulative step function (if jumps present) |
| `<comp>_exp_<date>` | Exponential relaxation term (if present) |
| `<comp>_log_<date>` | Logarithmic relaxation term (if present) |
| `<comp>_model` | Sum of all model components (fitted values) |
| `<comp>_noise` | Residual = observed − model |
| `<comp>_wtest` | Baarda w-statistic per observation: `e_i / (sigma . sqrt(1 - h_ii))` where `h_ii` is the hat-matrix diagonal |
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
| `Polyline breaks` | Velocity-break dates identified automatically by the DIA or user-specified. |
| `Exp relaxation` | Post-jump exponential relaxation terms (onset -> tau in days). |
| `Log relaxation` | Post-jump logarithmic relaxation terms. |
| `Exp trend (b/day)` | Exponential decay trend rate (null if not used). |
| `n_params` | Total number of fitted model parameters. Fewer is better (parsimony). |
| `Degrees of freedom (r)` | `n_obs - n_params`. Determines the chi-squared distribution shape for the OMT. |
| `T_stat (SSR/sigma²)` | Overall Model Test statistic: sum of squared residuals divided by sigma_0². Accepted when `T_stat <= chi²_critical`. |
| `chi²_critical (K)` | Chi-squared critical value at significance level alpha: `chi²_{1-alpha}(r)`. The formal threshold for model acceptance. |
| `Unit variance factor (T/r)` | `T_stat / r = sigma_hat² / sigma_0²`. Near 1.0 means well-calibrated; < 1 means slight overfit; > 1 means underfit. |
| `p-value` | Tail probability `1 - chi²_CDF(T, r)`. Accepted when `p-value >= alpha`. Equivalent to `T_stat <= chi²_critical`. |
| `DIA iterations` | Number of Detection–Identification–Adaptation cycles needed before acceptance. |

---

## 8. Statistical Concepts (Plain Language)

### Sigma (sigma_0) — "Expectation of Messiness"

Sigma (`sigma_mm` in the report) is your **assumed a-priori measurement noise** — how noisy you expect the data to be before fitting any model. The sigma scan sweeps a range (e.g., 2–20 mm) and selects the smallest sigma at which the Overall Model Test passes. The selected sigma is therefore **a-posteriori**: it was determined by observing the data, not set independently in advance.

**Two selection strategies:**

- **Default (parsimony):** Among all sigma values that produced an accepted model, the one with the fewest parameters wins. This gives the simplest model that still describes the data. Ties are broken by sigma (tighter noise preferred), then by p-value.
- **`--auto-sigma` (finest segmentation):** Selects the sigma that produces the most polyline breakpoints. This is useful when you know the signal contains many velocity changes (e.g., a slow-slip event sequence) and you want the model to resolve the finest segmentation that still passes the OMT. Use this together with `--poly-deg 1` (linear trend) for a pure piecewise-linear fit.

### Overall Model Test — T_stat, chi2_critical, and p-value

The **Overall Model Test (OMT)** is the formal chi-squared test that decides whether a model fits. The test statistic is:

```
T_stat = SSR / sigma_0²
```

where SSR is the sum of squared residuals. Under a correctly specified model, `T_stat ~ chi²(r)` (chi-squared with `r` degrees of freedom).

**Acceptance criterion:** `T_stat <= chi²_critical` — equivalently, `p-value >= alpha` (default alpha = 0.05).

- **T_stat <= chi²_critical** -> model **accepted** (residuals are consistent with assumed noise)
- **T_stat > chi²_critical** -> model **rejected** (residuals are too large)

The p-value is just the tail probability `1 - chi²_CDF(T_stat, r)`, so checking `p >= 0.05` is mathematically identical to checking `T_stat <= K`.

### Unit Variance Factor and A-Posteriori Sigma

Two fields in the report characterise how well the assumed sigma matches the actual data variability after fitting:

**Unit variance factor** = `T_stat / r = sigma_hat² / sigma_0²`
- Near 1.0 -> sigma assumption was well-calibrated
- < 1.0 -> model slightly over-specified, or sigma was too large
- > 1.0 -> residuals larger than expected; sigma may be underestimated

**Sigma_hat a-posteriori** = `sqrt(SSR / r)` in mm — the estimated actual noise level from the data. Compare to `Sigma_0 assumed`: close agreement means good calibration.

### Degrees of Freedom (r)

```
r = n_obs - n_params
```

Degrees of freedom control the chi-squared distribution shape used for the OMT. Adding model parameters (more seasonal components, jump terms, relaxation terms) reduces `r`, which raises the `chi²_critical` threshold. The tool prefers **parsimonious models** — when multiple models pass at the same sigma, the one with fewest parameters (parsimony) is preferred; further ties are broken by p-value.

### w-statistic — "Is This Specific Point an Outlier?"

The w-statistic implements the **formal Baarda datasnooping test** per observation epoch:

```
w_i = e_i / (sigma_0 . sqrt(1 - h_ii))
```

- `e_i` = residual at epoch i
- `h_ii` = hat-matrix diagonal element (leverage), computed via thin QR decomposition
- Under H_0 (no outlier), `w_i ~ N(0, 1)`

If `|w_i| > 3.29`, the observation is flagged as an **anomalous observation** (two-sided test at 0.1% significance level). The CSV `flagged` column marks these rows.

> The simplified formula `w = residual / sigma` (without the `sqrt(1 - h_ii)` correction) is an approximation. The implementation uses the full Baarda formula, which correctly accounts for the leverage of each observation.

### Jump Detection — Pre-Loop Stage

Before the DIA loop starts, `detect_jumps` scans the timeseries for abrupt level shifts using two complementary methods:

1. **Rolling-diff stage** — computes day-to-day differences and flags spikes exceeding a local MAD-based threshold, then validates each spike against a 30-day median-filtered trend. Catches instantaneous co-seismic jumps or equipment resets with dense surrounding data.
2. **Gap-level-shift stage** — scans the *raw* (non-gap-filled) series for data gaps longer than 30 days and measures the level shift across each gap (median of 30 observations before vs. 30 after). Flags the resumption date when the shift exceeds the adaptive threshold. Catches post-outage instrument resets where the rolling-diff stage would fail because the gap itself is filled by interpolation.

User-supplied `--jumps` dates are merged in after both stages. The `--no-jump` flag disables *both* automatic stages while still respecting manually specified `--jumps` dates.

A 90-day minimum spacing is enforced between all detected jump dates.

### DIA Loop — Detection, Identification, Adaptation

The sigma scan runs a DIA loop at each candidate sigma value. The loop tests four distinct hypothesis groups, each controlled by its own flag:

**1. Hypothesis group 1: Datasnooping (single-epoch outlier)**

The w-statistic is computed for every observation epoch. The epoch with the highest `|w|` above the critical value `z_{alpha/2}` is flagged as an outlier and converted into a jump (unit step) at that date, effectively removing it from the estimation.

- **Controlled by:** `--no-jump` (disables this group entirely)
- **Physical motivation:** A spiked measurement from a sensor glitch, lightning strike, or transient disturbance that should not propagate into the trend or seasonal estimates.

**2. Hypothesis group 2: Missing periodic signal (Lomb-Scargle)**

The residuals are searched for dominant spectral peaks using the Lomb-Scargle periodogram. The top 5 peaks are tested as candidate sine/cosine pairs using the formal w-test. Only periods within 0.05 years of the allowed set (the merged list of `--periods` and auto-detected `--auto-periods` candidates) are accepted. This prevents the DIA from adding periods that were not in the original candidate library.

- **Controlled by:** `--no-seasonal` (disables this group entirely). The allowed set is further restricted by `--periods` and `--auto-periods` (which define what periods exist in `allowed_periods`).
- **Physical motivation:** An unmodeled annual or semi-annual groundwater cycle, a thermal expansion signal, or a multi-year drought/recharge pattern.

**3. Hypothesis group 3: Velocity break / polyline (piecewise-linear)**

Two sub-strategies identify candidate break dates:

- **Velocity spike:** The first difference of the residuals is examined. The epoch with the largest positive or negative velocity spike (relative to the MAD-scaled median velocity) is tested as a polyline break.
- **CUSUM break:** The Ploberger-Kramer (1992) normalised CUSUM range test identifies structural breaks in the residual mean. If the CUSUM range exceeds the 1.36 threshold (~95% significance), the break date is tested.

The polyline candidate with the highest `|w|` wins. This hypothesis is never disabled — it is always active.

- **Physical motivation:** A change in subsidence rate due to new pumping regulations, an earthquake-induced poroelastic response, or a construction load onset.

**4. Hypothesis group 4: Exponential trend**

The AIC-based `auto_detect_exp_trend` function tests 20 candidate `b` values (decay rates) against the current residuals. A candidate is accepted when two guards pass: delta-AIC > 10 (strong evidence) and the exponential component explains at least 10% of the series variance. The winning `b` value is tested with the formal w-test.

- **Controlled by:** `--no-exp-trend` (disables this group entirely). Also affected by `--exp-trend`: if a pre-specified `b` value or auto-detected `b` is already in the model (set before the DIA runs), this group is not evaluated because `exp_trend` is not `None`.
- **Physical motivation:** A compaction well layer showing monotonic consolidation — rapid initial deformation that tapers over years toward a final settlement level. The exponential shape `a.(exp(-b.t)-1)` describes this: zero at the first observation, asymptotic toward amplitude `a` at infinite time.

**Selection rule (Teunissen):** Among all candidates from all four groups, the one with the highest `|w|` is adapted, provided `|w_winner| > z_{alpha/2}` (the two-tailed normal critical value). If no candidate passes this gate, the DIA stops.

**Adaptation** incorporates the winning component into the model, then the next iteration begins with re-estimation.

**Convergence guard:** The DIA stops when the unit variance factor stops improving between iterations (i.e., is no longer decreasing). This prevents the loop from adding components that provide negligible improvement, even before `--max-iter` is reached.

### Exponential Decay Trend — "What is it?"

Some signals do not increase or decrease at a steady rate — they start fast and then slow down, asymptoting toward a final plateau. Groundwater compaction and post-pumping recovery are typical examples.

The model used is:

```
displacement(t) = a . (exp(-b . t) - 1)
```

- `a` — the amplitude (total displacement at infinite time, in metres)
- `b` — the decay rate in 1/days. Larger `b` means faster convergence. The time constant tau = 1/b.
- At `t = 0`: displacement = 0 (anchored to the first observation)

Because `b` is fixed before OLS fitting (amplitude `a` is solved linearly), this component integrates cleanly with the OMT statistical framework. When using `--exp-trend auto`, the tool scans 20 candidate `b` values and selects the one with the lowest AIC over a plain linear model. The threshold is **delta-AIC > 10** (strong evidence criterion) to avoid adding an exponential trend for weak signals.

**Three ways to use exponential trend:**

1. **Omit `--exp-trend` entirely** — the DIA loop may detect an exponential trend as hypothesis group 4 during the identification phase. This is data-driven detection.
2. **`--exp-trend auto`** — scans for the best `b` value *before* the DIA loop starts, and the detected `b` is baked into the initial model. The DIA group 4 is then skipped because `exp_trend` is already non-None.
3. **`--exp-trend 0.001`** (numeric) — you specify the decay rate directly. The DIA group 4 is skipped.
4. **`--no-exp-trend`** — explicitly prevents any exponential trend from being added, by either pre-scan or DIA detection.

**When to use each mode:**

- For MLCW compaction-well layers: start with `--exp-trend auto` to see if the signal has a clear exponential decay shape. If the auto-detection finds no significant trend (skip report or no `_exp_trend` column in the CSV), fall back to the default linear model.
- For GPS with post-seismic deformation: omit `--exp-trend` and let the DIA detect it as part of the identification phase, since exponential transients here typically follow a known jump.
- For InSAR with steady subsidence: use `--no-exp-trend` unless you have a physical reason to expect exponential decay (e.g., a site transitioning from elastic to inelastic compaction).

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

### Example 3 — MLCW Layer with Exponential Decay Trend, No Seasonal, No Jumps

For a compaction-well layer showing monotonic consolidation without seasonal pumping cycles or abrupt steps:

```bash
python -m appsigsolv decompose XIGANG_ringbyring.csv \
    --component 102.096 \
    --date-col datetime \
    --unit mm \
    --poly-deg -1 \
    --no-seasonal \
    --no-jump \
    --no-relax \
    --exp-trend auto \
    --periods "" \
    --auto-periods 0 \
    --irregular \
    --output-dir ./MLCW_decomposition
```

`--no-seasonal` suppresses all seasonal detection (no annual or semi-annual pumping cycles). `--no-jump` prevents step detection (useful when instrument resets are already corrected). `--exp-trend auto` lets the tool find the best consolidation decay rate. The exponential trend component (`<comp>_exp_trend`) appears in the decomposed CSV and as a dashed purple line in the trend panel of the diagnostic plot.

If the auto-detected b value is worth keeping, reconstruct at specific dates:

```bash
python -m appsigsolv reconstruct XIGANG_ringbyring.csv \
    --json MLCW_decomposition/XIGANG_ringbyring/XIGANG_ringbyring_model_102.096.json \
    --target-col 102.096 \
    --sampling-rate custom \
    --custom-dates 1,15
```

### Example 4 — Restrict Date Range to a Subset of the Record

For a GPS station where only the post-2015 period is relevant (e.g., after a known equipment upgrade):

```bash
python -m appsigsolv decompose gps_timeseries/CHIN_neu.csv \
    --component dU \
    --start-date 2015-01-01 \
    --force --no-relax
```

The timeseries is truncated before jump detection and DIA processing. Use `--end-date` together with `--start-date` to isolate a specific window of interest.

### Example 5 — Pure Piecewise-Linear Fit with Auto-Sigma

For a signal where you expect multiple velocity changes (e.g., a slow-slip sequence) and want the finest segmentation:

```bash
python -m appsigsolv decompose gps_timeseries/CHIN_neu.csv \
    --component dU \
    --poly-deg 1 \
    --no-seasonal \
    --no-jump \
    --auto-sigma \
    --sigma-min 1.0 --sigma-max 20.0 --sigma-step 0.5 \
    --force
```

`--poly-deg 1` restricts the trend to a linear trend per segment. `--auto-sigma` selects the sigma that produces the most polyline breakpoints rather than the fewest parameters.

### Example 6 — Custom Reconstruction for Comparison

After `decompose` produces a model, reconstruct it on exactly the 1st and 15th of every month:

```bash
python -m appsigsolv reconstruct gps_timeseries/TKJS_neu.csv \
    --json gps_timeseries/TKJS_neu/TKJS_neu_model_dU.json \
    --target-col dU \
    --sampling-rate custom \
    --custom-dates 1,15 \
    -o results/TKJS_semi_monthly.csv
```

---

## 10. Tips, Gotchas & FAQ

### When should I use `--irregular`?

Use `--irregular` whenever your observations are **not daily** (e.g., MLCW or monthly InSAR). It prevents the tool from creating excessive synthetic points via resampling.

### When should I use `--no-seasonal`?

Use `--no-seasonal` when you know the signal has no periodic component — for example, a deep aquitard layer where seasonal pumping does not reach, or a compaction well in a purely inelastic deformation regime. Without this flag, the DIA may waste iterations trying to fit harmonics to noise. It also speeds up processing because Lomb-Scargle period detection is skipped.

### When should I use `--no-jump`?

Use `--no-jump` when the timeseries has already been corrected for steps (instrument resets, co-seismic offsets) in preprocessing, or when the signal genuinely has no abrupt changes (e.g., a compaction layer measured by a well-calibrated MLCW). Without this flag, the DIA may identify outlier epochs as step functions, adding unnecessary parameters.

### When should I use `--no-exp-trend`?

Use `--no-exp-trend` when the signal is well described by a linear or polynomial trend and you want to prevent the algorithm from chasing an exponential shape in the residuals. This is especially useful for InSAR timeseries where steady, linear subsidence is expected and any curvature comes from atmospheric noise rather than physical decay.

### When is `--auto-sigma` useful?

Use `--auto-sigma` when the physical process you are studying produces velocity changes rather than steps or seasonal cycles. Typical targets: slow-slip events in GPS, stick-slip motion in creepmeters, staged construction settlement. The selection criterion flips from "fewest parameters" to "most polyline breakpoints."

### How does automatic polynomial selection work?

If you pass `--poly-deg -1`, the tool tests degrees from `--poly-deg-min` (default `0`) up to `3`. It selects the degree whose accepted model passes OMT at the **smallest sigma** (tighter noise = better signal fit). Among degrees that pass at the same sigma, the one with fewer parameters wins; further ties are broken by p-value.

### How do I prevent offset-only (degree 0) models?

Use `--poly-deg-min 1` together with `--poly-deg -1`. This restricts the auto-selection to degrees 1–3, guaranteeing at least a linear trend in the accepted model. In batch scripts, add `poly_deg_min=1` to the `Namespace`:

```python
args = Namespace(
    poly_deg=-1,
    poly_deg_min=1,   # exclude degree-0 offset-only models
    ...
)
```

### How do I re-run a component I am not happy with?

Use `--force` with the specific `--component` name and any adjusted parameters. Without `--force`, the existing JSON will be detected and skipped.

### Output CSV units vs. input units

Input is mm (or m if `--unit m`). **All output CSV columns are in millimetres (mm).**

### How do I know if a station was skipped?

After a batch run, check the output directory for `*_skipped_*.txt` files — one per skipped component. Each file records the component name, the reason (no model or timeout), elapsed time, and the sigma/degree range that was tried. The raw data plot is still written so you can visually inspect whether the signal is physically meaningful or the sensor malfunctioned.

### Why does my batch script stop silently mid-run?

On Windows, running many stations in a single Python process can exhaust GDI handles as matplotlib accumulates figure objects. Best practice for batch scripts:

1. Use `cores=1` — avoids `ProcessPoolExecutor` spawn issues on Windows.
2. Call `plt.close('all')` and `gc.collect()` between stations in the `finally` block.
3. Import `matplotlib` and call `matplotlib.use("Agg")` at the top of the batch script before importing `appsigsolv`.

### Why do MLCW column names have ugly decimal filenames?

Float column names in CSVs (e.g. `86.15899999999999`) are a pandas read artefact. `appsigsolv` automatically normalises all numeric column names to 3 decimal places (`86.159`) during loading, so all output filenames, JSON keys, and CSV headers use the clean form.

### What happens with `--start-date` / `--end-date`?

Date filtering is applied **after** preprocessing (gap-filling, outlier removal) and **before** jump detection. This means:
- Outliers near the date boundaries are still removed based on the full-series MAD.
- Jump detection sees only the filtered window, so jumps outside the window do not affect the model.
- The reference epoch (t=0) for the polynomial and exponential trend is the first date in the filtered window.

---

## 11. Changelog

### v0.5.0 (2026-06)

**New flags for controlling DIA hypothesis groups:**

- **`--no-jump`:** Disables automatic jump detection (hypothesis group 1 — datasnooping/outlier). Manually specified `--jumps` dates are still applied. Essential for MLCW data where instrument resets have already been corrected in preprocessing.
- **`--no-seasonal`:** Disables all periodic/seasonal detection (hypothesis group 2 — Lomb-Scargle period search). Overrides `--periods` and `--auto-periods`. Speeds up processing for non-seasonal signals like deep aquitard compaction.
- **`--no-exp-trend`:** Disables exponential trend auto-detection (hypothesis group 4). Sets a sentinel value that prevents the DIA from identifying an exponential decay trend, even if the residuals show one. Use this for InSAR or GPS signals where linear subsidence is expected and curvature is noise.
- **`--auto-sigma`:** Selects the sigma value that maximises polyline breakpoints rather than minimising parameter count. Designed for pure piecewise-linear fitting where the goal is the finest velocity segmentation that still honours the OMT.

**Date range filtering:**

- **`--start-date` / `--end-date`:** Filter the timeseries to a specific date window before jump detection and DIA processing. Input format: `YYYY-MM-DD` (inclusive). Applied after gap-filling and outlier removal.

**Documentation improvements:**

- DIA hypothesis groups now documented individually with their controlling flags and physical motivation.
- `allowed_periods` mechanism explained: the DIA can only add periods within 0.05 years of the candidate set (merged from `--periods` and `--auto-periods`).
- Saved JSON filtering documented: `no_seasonal`, `no_jump`, and `allowed_periods` are runtime-only and stripped from persisted JSON files.
- New worked examples for date-range filtering, auto-sigma piecewise-linear fitting, and combined suppression flags for MLCW processing.

**`--exp-trend` interaction clarified:**

Three modes now documented: (1) omit `--exp-trend` to let DIA detect it, (2) `--exp-trend auto` for pre-scan, (3) numeric value for user-specified decay rate, plus (4) `--no-exp-trend` to disable entirely.

### v0.4.1 (2026-05)

**Bug fixes:**
- **Timing regression fix:** A step-function hypothesis (Heaviside H(t-t0)) introduced in a prior commit caused an 8x slowdown and widespread timeout failures (e.g., GFES, KTES). Root cause: the iterative DIA loop found CUSUM structural breaks on every iteration for under-assumed sigma values, accumulating up to 30 spurious step dates and bloating the design matrix. The in-loop step hypothesis has been removed entirely.
- **GFES/KTES restored:** Both stations that were incorrectly producing `_skipped_dU.txt` (timeout) now complete and produce accepted models as before.

**Improvements:**
- **Gap-level-shift detection in `detect_jumps`:** A second detection stage now scans raw (non-gap-filled) series for data gaps longer than 30 days and measures the level shift across each gap (median of 30 raw observations before vs. after). If the shift exceeds the adaptive MAD threshold, the resumption date is flagged as a jump. This catches post-outage instrument resets — such as a GPS station that went offline for 235 days and resumed at a 113 mm different level — that the original rolling-diff stage missed because the gap was filled by time-interpolation before differencing.

### v0.4.0 (2026-05)

**New features:**
- **Auto-skip with skip report:** When the sigma scan finds no accepted model, or when fitting exceeds 180 s, `decompose` writes a `_skipped_<comp>.txt` report and continues to the next component. Previously the batch would hang or silently print a WARNING with no output file. The skip report records the component name, reason, elapsed time, and parameters tried.
- **3-minute per-component timeout:** Fitting runs in a daemon thread; the main thread waits up to 180 s. If it exceeds this, the thread is abandoned, a skip report is written, and the batch continues. This prevents malfunctioned or highly irregular timeseries from blocking an entire batch on Windows (where `signal.SIGALRM` is unavailable).

**Performance:**
- **18x Lomb-Scargle speedup in DIA loop:** The `days[]` and `freqs[]` arrays in `_identify_best_alternative` are now pre-computed once before the DIA iteration loop (they are loop-invariant) and passed in, instead of being recomputed every iteration. The frequency grid was also reduced from 5 000 to 500 points — sufficient to resolve all peaks between 0.2–20 yr, cutting per-call time ~10x. Combined effect: a 15-year GPS series (5 476 pts, CHIN station) drops from ~18 minutes to ~60 seconds.

**Output / UX:**
- **Reduced verbose output:** Removed per-component log lines that cluttered batch output: `[periods] Forcing...`, `Series loaded: N points`, `[extract] Components: [...]`, and all four `[output] ... saved` lines. What remains: preprocessing stats, sigma-scan result, auto-deg selection, skip/timeout notices, and errors.

### v0.3.0 (2026-05)

**New features:**
- `--force` flag for `decompose`: overwrite existing results instead of skipping. Useful when manually re-processing unsatisfactory components with adjusted parameters.

**Bug fixes & robustness:**
- **Batch crash fix (Windows GDI exhaustion):** `save_plot` now calls `plt.close(fig)` + `plt.close('all')` to guarantee figure handles are released after every plot. Batch scripts should additionally call `plt.close('all')` + `gc.collect()` between stations.
- **Numeric column name normalisation:** All floating-point CSV column name artefacts (e.g. `86.15899999999999`) are now rounded to 3 decimal places at load time in both `load_and_preprocess` and the `component="all"` path of `decompose`. Output filenames are clean (e.g. `86.159`).
- **Per-component exception handling:** The entire fitting + extraction + output pipeline for each component is now wrapped in `try/except Exception`, so a single failing component logs its error and continues — it no longer crashes the entire batch or station run.
- **Adaptive plot x-axis:** `YearLocator` interval is now data-driven (1 yr for spans < 4 yr, 2 yr for 4–10 yr, 5 yr for > 10 yr), preventing a matplotlib crash on short timeseries.

**Algorithm:**
- **`--poly-deg-min`:** New option to set the minimum polynomial degree when using `--poly-deg -1`. Default `0` (old behaviour). Set to `1` to exclude offset-only models — useful when data clearly has a trend and degree-0 results are not physically meaningful.
- **Sigma-first cross-degree selection:** When `--poly-deg -1`, the winning degree is now chosen by the smallest accepted sigma (tightest noise), not best p-value. Ties are broken by parameter count, then p-value.
- **`cores` default changed to `1`:** Avoids Windows `ProcessPoolExecutor` / `spawn` issues when calling `run_decompose` directly from a batch script (not under `if __name__ == '__main__'`).

### v0.2.0 (prior)

- Exponential decay trend component (`--exp-trend auto` / numeric b)
- Auto polynomial degree selection (`--poly-deg -1`)
- Sigma scan with OMT-based model acceptance
- DIA loop with Lomb-Scargle period detection, jump detection, relaxation testing
- `reconstruct` command with custom sampling rates

---

*Package version: 0.5.0 | appsigsolv — Applied Signal Solver*
