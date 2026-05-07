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
7. [Understanding the Outputs](#7-understanding-the-outputs)
8. [Statistical Concepts (Plain Language)](#8-statistical-concepts-plain-language)
9. [Worked Examples](#9-worked-examples)
10. [Tips, Gotchas & FAQ](#10-tips-gotchas--faq)

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

### 4.3 Units

- Default input unit: `mm` (millimetres)
- Internal computation unit: `m` (metres) — all output CSVs are in metres
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
| `--poly-deg` | `1` | Polynomial degree: `0`=offset only, `1`=linear trend, `2`=acceleration. Use `-1` for auto-selection. |
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
| `--output-dir` | *(same folder as CSV)* | Parent directory for outputs |
| `--cores` | `1` | Number of CPU cores for parallel sigma scan |
| `--exp-trend` | *(off)* | Exponential decay trend `exp(-b·t)−1`. Pass `auto` to auto-detect the best decay rate via AIC, or a numeric `b` in 1/days (e.g. `0.001` ≈ 2.7-yr time constant). Omit to disable. |

### Resumption Capability

`decompose` now supports automatic resumption. If a batch run is interrupted, re-running the same command will automatically skip any components that already have a corresponding `{stem}_model_{comp}.json` file in the output directory. This ensures that only new or incomplete components are processed, saving significant time.

### Outputs

Results are saved in `<output_dir>/<csv_stem>/` (e.g., `gps_timeseries/TKJS_neu/`):

| File | Description |
|---|---|
| `<stem>_model_<comp>.json` | Accepted model configuration (reusable with `reconstruct`) |
| `<stem>_decomposed_<comp>.csv` | All decomposed signal components (in metres) |
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

All values are in **metres**. Each row is one date.

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
| `<comp>_wtest` | Normalised w-statistic per observation |
| `flagged` | `True` when `|w-test| > 3.29` (anomalous observation) |

---

## 8. Statistical Concepts (Plain Language)

### Sigma (σ) — "Expectation of Messiness"

Sigma is your assumed measurement noise. The tool scans a range and finds the smallest sigma for which the model statistically fits the data.

### p-value — "Does the Overall Model Fit?"

- **p-value ≥ 0.05** → model **accepted** (misfit is within expectations)
- **p-value < 0.05** → model **rejected** (misfit is too large)

### w-statistic — "Is This Specific Point an Outlier?"

For each observation: `w = residual / sigma`. If `|w| > 3.29`, the point is flagged as an anomaly.

### Exponential Decay Trend — "What is it?"

Some signals do not increase or decrease at a steady rate — they start fast and then slow down, asymptoting toward a final plateau. Groundwater compaction and post-pumping recovery are typical examples.

The model used is:

```
displacement(t) = a · (exp(-b · t) − 1)
```

- `a` — the amplitude (total displacement at infinite time, in metres)
- `b` — the decay rate in 1/days. Larger `b` means faster convergence. The time constant τ = 1/b.
- At `t = 0`: displacement = 0 (anchored to the first observation)

Because `b` is fixed before OLS fitting (amplitude `a` is solved linearly), this component integrates cleanly with the OMT statistical framework. When using `--exp-trend auto`, the tool scans 20 candidate `b` values and selects the one with the lowest AIC improvement over a plain linear model (threshold ΔAIC > 2).

---

## 9. Worked Examples

### Example 1 — Batch Processing with Resumption

```bash
python -m appsigsolv decompose gps_timeseries/TKJS_neu.csv --component all --cores 4
```

Processes all columns in parallel. If interrupted, re-running this command will skip already completed components.

### Example 2 — Custom Reconstruction for Comparison

After `decompose` produces a model, reconstruct it on exactly the 1st and 15th of every month:

```bash
python -m appsigsolv reconstruct gps_timeseries/TKJS_neu.csv \
    --json gps_timeseries/TKJS_neu/TKJS_neu_model_dU.json \
    --target-col dU \
    --sampling-rate custom \
    --custom-dates 1,15 \
    -o results/TKJS_semi_monthly.csv
```

### Example 3 — MLCW Column with Exponential Decay Trend

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

If you pass `--poly-deg -1`, the tool will test degrees 0, 1, and 2. It selects the one that results in an accepted model with the fewest parameters.

### Output CSV units vs. input units

Input is mm (or m if `--unit m`). **All output CSV columns are in metres.**

---

*Package version: 0.2.0 | appsigsolv — Applied Signal Solver*
