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

**Decompose all components in a file:**

```bash
python -m appsigsolv decompose gps_timeseries/TKJS_neu.csv --component all
```

**Reconstruct a modeled timeseries from a saved JSON config:**

```bash
python -m appsigsolv reconstruct gps_timeseries/TKJS_neu.csv \
    --json gps_timeseries/TKJS_neu/TKJS_neu_model_dU.json \
    --target-col dU --daily
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
| `--poly-deg` | `1` | Polynomial degree: `0`=offset only, `1`=linear trend, `2`=acceleration |
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

The `reconstruct` command fits a **user-specified** model (either from a JSON file or CLI flags) to a timeseries and saves the fitted curve. Use this to:

- Apply a known model from `decompose` to new data
- Build a specific parametric model manually and fit it
- Forward-model on a dense daily grid

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
| `--period` | | Periodic component in years (repeatable: `--period 0.5 --period 1.0`) |
| `--step` | | Jump step date `YYYYMMDD` (repeatable) |
| `--polyline` | | Polyline break date `YYYYMMDD` (repeatable) |
| `--exp DATE TAU` | | Exponential term: onset date `YYYYMMDD` and tau in days (repeatable) |
| `--log DATE TAU` | | Logarithmic term: onset date `YYYYMMDD` and tau in days (repeatable) |
| `-o` / `--output` | `<input>_modeled.<ext>` | Output file path |
| `--daily` | *(flag)* | Output on a dense daily grid instead of at observation times |
| `--ref-date` | *(first date)* | Reference epoch for model `YYYYMMDD` |

### Output

A single CSV (or Excel) identical to the input plus:
- `modeled` column — fitted model values in metres
- (with `--daily`) `reconstructed` column — dense daily model

### CLI flags override JSON

If you load `--json` and also pass `--poly`, `--period`, etc., the CLI flags take precedence over the JSON values. This lets you load a base model and tweak it without editing the JSON file.

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
    "log": {}
}
```

- `polynomial`: degree of the trend (0, 1, or 2)
- `periodic`: list of accepted period lengths in years
- `stepDate`: list of jump dates in `YYYYMMDD` format
- `polyline`: list of velocity-change dates
- `exp` / `log`: dict of `{"YYYYMMDD": tau_days}` for relaxation terms

This file is directly usable as `--json` input to `reconstruct`.

### 7.2 Decomposed CSV

All values are in **metres**. Each row is one date.

| Column | Description |
|---|---|
| `date` | ISO date |
| `<comp>` | Observed displacement |
| `<comp>_trend` | Polynomial trend component |
| `<comp>_<T>yr` | Seasonal component at period T years |
| `<comp>_jump` | Cumulative step function (if jumps present) |
| `<comp>_exp_<date>` | Exponential relaxation term (if present) |
| `<comp>_log_<date>` | Logarithmic relaxation term (if present) |
| `<comp>_model` | Sum of all model components (fitted values) |
| `<comp>_noise` | Residual = observed − model |
| `<comp>_wtest` | Normalised w-statistic per observation |
| `flagged` | `True` when `\|w-test\| > 3.29` (anomalous observation) |

### 7.3 Diagnostic Figure

The PNG has a two-panel layout:

- **Left panel:**
  - Upper 2/3 — Observed (scatter) vs. Model (line) in mm
  - Lower 1/3 — Long-term trend only

- **Right panel** — dynamically stacked subplots showing whichever of these are present:
  - Seasonal components
  - Jumps & relaxation terms
  - Residual noise with flagged anomalies highlighted

The figure title shows the accepted sigma, p-value, number of parameters, and polynomial degree.

### 7.4 Statistical Report

A Markdown file with four sections:

1. **Accepted Model** — all model parameters, OMT statistic, p-value, DIA iterations used
2. **Variance Explained per Component** — std (mm) and % variance for each component
3. **Sigma Scan Summary** — table showing accepted/rejected sigma values with p-values
4. **Anomalous Observations** — dates where `|w-stat| > 3.29` with residual size in mm

---

## 8. Statistical Concepts (Plain Language)

### Sigma (σ) — "Expectation of Messiness"

Sigma is your assumed measurement noise — how much random scatter you expect in a single daily observation. The tool does not require you to know this value in advance; instead it **scans a range** of sigma values and finds the smallest one for which the model statistically fits the data.

A smaller accepted sigma means the data is clean and the model explains it well. A large accepted sigma means there is substantial scatter or un-modeled signal.

### p-value — "Does the Overall Model Fit?"

The Overall Model Test (OMT) computes a chi-squared statistic from the sum of squared residuals and asks: *given our assumed sigma, is this level of misfit plausible by random chance alone?*

- **p-value ≥ 0.05** → model **accepted** (misfit is within expectations)
- **p-value < 0.05** → model **rejected** (misfit is too large; the model is missing something)

The DIA loop then adds model terms (a new period or a velocity break) and retests until acceptance.

### w-statistic — "Is This Specific Point an Outlier?"

For each observation: `w = residual / sigma`

- `|w| > 3.29` → the observation is flagged as anomalous (0.1% probability under the null)
- Flagged observations appear in the `flagged` column of the CSV and in the report table

Flagged points are **not removed** from the fit — they are retained and reported so you can decide their cause (instrument spike, real geophysical event, data error).

### DIA Loop — How the Model Is Built Automatically

1. Fit model → compute OMT p-value
2. If p-value < 0.05 (rejected): call Lomb-Scargle on residuals to find the strongest unmodeled frequency, or check for a velocity-break pattern
3. Add the identified term to the model and re-fit
4. Repeat up to `--max-iter` times per sigma value
5. If still rejected after max iterations: try the next sigma value and restart

The first sigma (starting from `--sigma-min`) that achieves acceptance is selected as the final model. In case of ties, the model with the fewest parameters wins; if still tied, the highest p-value wins.

### Relaxation Testing

After acceptance, if `--no-relax` is not set, the tool tests whether adding an exponential relaxation term after each detected jump improves the model (tau = 30, 90, or 180 days). This captures post-seismic or hydrological aftereffects.

---

## 9. Worked Examples

### Example 1 — GPS Up Component (Standard)

```bash
python -m appsigsolv decompose gps_timeseries/TKJS_neu.csv --component dU
```

Default settings: linear trend (`--poly-deg 1`), candidate periods 0.25, 0.5, 1.0, 2.0 yr plus up to 5 auto-detected, sigma scan from 2 to 15 mm.

**Expected outputs in `gps_timeseries/TKJS_neu/`:**
- `TKJS_neu_model_dU.json`
- `TKJS_neu_decomposed_dU.csv`
- `TKJS_neu_decomposed_dU.png`
- `TKJS_neu_report_dU.md`

### Example 2 — All Three GPS Components in Parallel

```bash
python -m appsigsolv decompose gps_timeseries/TKJS_neu.csv \
    --component all \
    --cores 4
```

Processes `dN`, `dE`, `dU` (and any other non-date columns) simultaneously using 4 CPU cores for the sigma scan.

### Example 3 — GPS with a Known Earthquake Jump

If a Mw 6.2 earthquake occurred on 2022-06-17 and you want to force a jump regardless of auto-detection:

```bash
python -m appsigsolv decompose gps_timeseries/TKJS_neu.csv \
    --component dU \
    --jumps 2022-06-17
```

### Example 4 — GPS with Acceleration (Quadratic Trend)

For a station showing clear acceleration (e.g., accelerating subsidence):

```bash
python -m appsigsolv decompose gps_timeseries/TKJS_neu.csv \
    --component dU \
    --poly-deg 2
```

### Example 5 — MLCW Irregular Data

```bash
python -m appsigsolv decompose mlcw_timeseries/TUKU_ringbyring.csv \
    --component all \
    --irregular \
    --unit mm \
    --poly-deg 1 \
    --no-relax
```

`--irregular` is required because observations are approximately monthly, not daily. `--component all` processes every depth layer column.

### Example 6 — Widen the Sigma Scan

If your data has high noise (e.g., InSAR coherence issues) and the default scan (2–15 mm) never finds an accepted model, extend the range:

```bash
python -m appsigsolv decompose my_insar.csv \
    --component displacement \
    --sigma-min 5.0 \
    --sigma-max 30.0 \
    --sigma-step 1.0
```

### Example 7 — Reconstruct from a Saved JSON Model

After `decompose` produces `TKJS_neu_model_dU.json`, apply the same model to produce a dense daily reconstruction:

```bash
python -m appsigsolv reconstruct gps_timeseries/TKJS_neu.csv \
    --json gps_timeseries/TKJS_neu/TKJS_neu_model_dU.json \
    --target-col dU \
    --daily \
    -o results/TKJS_neu_daily_model.csv
```

### Example 8 — Manual Model with `reconstruct`

Manually specify a linear trend + annual + semi-annual + one jump, without running the full DIA pipeline:

```bash
python -m appsigsolv reconstruct gps_timeseries/TKJS_neu.csv \
    --target-col dU \
    --poly 1 \
    --period 1.0 \
    --period 0.5 \
    --step 20220617 \
    -o results/TKJS_manual_model.csv
```

---

## 10. Tips, Gotchas & FAQ

### When should I use `--irregular`?

Use `--irregular` whenever your observations are **not daily**: MLCW monthly surveys, InSAR with variable revisit time, or any timeseries where resampling to a daily grid would create too many synthetic points. Without this flag, the tool resamples to daily and fills gaps up to 7 days, which is appropriate for GPS but not for coarser data.

### The sigma scan never finds an accepted model — what do I do?

1. **Extend the sigma range**: use a larger `--sigma-max` (e.g., 30 or 50 mm).
2. **Add more iterations**: `--max-iter 10` gives the DIA loop more chances to add model terms.
3. **Force known jumps**: if you know an event date, add it with `--jumps YYYY-MM-DD` so the tool does not have to auto-detect it.
4. **Add a polyline**: if there is a known velocity change, use `--polylines YYYY-MM-DD`.
5. **Check the data**: very large data gaps or systematic seasonal patterns not well-represented by sinusoids may prevent acceptance.

### What does `--component all` do exactly?

It reads every column name from the CSV header (excluding the auto-detected date column) and runs the full decompose pipeline independently on each one. Output files are named per component. For MLCW data with dozens of depth layers, this is a convenient batch mode.

### Can I reload and re-apply a model JSON from a previous run?

Yes. The JSON produced by `decompose` is the direct input format for `reconstruct --json`. You can also edit the JSON manually (e.g., add or remove a period value) and then reload it.

### How do I interpret very high w-statistics (e.g., w > 10)?

A w-statistic of 10+ means that observation is roughly 10× farther from the model than expected by noise. Common causes:
- A data transcription error (wrong unit, sign flip)
- A real but unmodeled abrupt event (earthquake, instrument swap)
- The sigma is too small (the tool still accepted the model at this sigma, but individual points are extreme)

The `flagged` column in the CSV marks all such points (`|w| > 3.29`). They remain in the fit; removal is left to user judgment.

### Can I use `reconstruct` on a different station with the same model?

Yes — load the JSON from station A and point `reconstruct` at station B's CSV. The model structure (periods, polynomial degree, jump dates) will be applied to the new data. This is useful for comparing stations or forcing a reference model.

### Output CSV units vs. input units

Input data is assumed to be in mm (or m if `--unit m`). All output CSV columns (including components) are in **metres**. Multiply by 1000 to convert back to mm for plotting or comparison.

### Parallelisation (`--cores`)

The sigma scan is embarrassingly parallel (each sigma value is independent). Setting `--cores 4` can reduce runtime by ~4× for wide sigma ranges. On Windows, make sure your script is called from a `if __name__ == "__main__":` guard if embedding the call in a script (not required when using the CLI directly).

---

*Package version: 0.1.0 | appsigsolv — Applied Signal Solver*
