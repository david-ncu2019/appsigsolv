# appsigsolv — Applied Geology's Signal Solver

[![Python](https://img.shields.io/badge/Python-3.10+-blue.svg)]()
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

A CLI tool for parametric time series decomposition built on the **Overall Model Test (OMT)** framework of Teunissen (TU Delft). It decomposes a deformation time series into physically interpretable components — trend, seasonal cycles, velocity breaks, jumps, and relaxation signals — while using statistical hypothesis testing to decide what belongs in the model.

The name is short for **Applied Geology's Signal Solver**.

---

## What it does

Given a time series of displacement (GPS, InSAR, MLCW extensometer, or groundwater level), appsigsolv builds a parametric model of the form:

```
d(t)  =  poly(t)            — polynomial trend (velocity, acceleration)
       + Σ aₖ·sin(ωₖt+φₖ)   — periodic signals (annual, semi-annual cycles)
       + Σ cⱼ·H(t−tⱼ)       — step offsets (jumps, co-seismic offsets)
       + Σ mⱼ·(t−tⱼ)·H(t−tⱼ) — polyline velocity breaks
       + Σ p·exp(−t/τ)      — exponential relaxation
       + Σ q·ln(1+t/τ)      — logarithmic relaxation
       + b·(exp(−β·t)−1)    — exponential trend (e.g. poroelastic decay)
```

The **DIA procedure** (Detection-Identification-Adaptation, Teunissen 2017) decides which of these components are significant:

1. **Detect** — the Overall Model Test (χ² test) checks if the residuals match the assumed noise level σ.
2. **Identify** — if the OMT fails, four hypothesis groups compete to explain the misfit:
   - *Group 1* — single-epoch outlier (datasnooping via w-test)
   - *Group 2* — missing periodic signal (Lomb-Scargle periodogram)
   - *Group 3* — velocity break / polyline (CUSUM on residuals → spike detection)
   - *Group 4* — exponential trend (AIC-based exponential fit)
3. **Adapt** — the winning hypothesis is added to the model; the loop repeats until the OMT passes.

The result is a model that is **as simple as possible but as complex as the data demands** — a parsimonious decomposition driven by statistical evidence, not guesswork.

---

## Quick start

```bash
pip install appsigsolv
```

**Decompose a GPS vertical component:**

```bash
appsigsolv decompose station.csv --component dU --periods "0.5,1.0" --output-dir ./results
```

**Decompose a multilayer compaction record (MLCW), full timeline:**

```bash
appsigsolv decompose mlcw.csv --component all --periods "0.5,1.0" \
  --sigma-min 2.0 --sigma-max 10.0 --sigma-step 0.5 \
  --poly-deg 1 --no-jump --output-dir ./results
```

**Piecewise-linear only (no seasonal, fit the segmentation):**

```bash
appsigsolv decompose timeseries.csv --component all \
  --no-seasonal --no-jump --auto-sigma --poly-deg 1 \
  --sigma-min 1.0 --sigma-max 8.0 --sigma-step 0.5
```

**Filter to a specific time window:**

```bash
appsigsolv decompose timeseries.csv --component dU \
  --start-date 2014-01-01 --end-date 2021-12-31
```

---

## Installation

### PyPI

```bash
pip install appsigsolv
```

### From source

```bash
git clone https://github.com/david-ncu2019/appsigsolv.git
cd appsigsolv
pip install -e .
```

**Dependencies:** `numpy`, `scipy`, `pandas`, `matplotlib` (automatically installed).

---

## Command reference

### `decompose` — parametric signal decomposition

```
appsigsolv decompose <input_csv> [options]
```

| Argument | Default | Description |
|----------|---------|-------------|
| `--component` | `dU` | Component to process, or `all` for batch multi-column |
| `--date-col` | auto | Date column name |
| `--unit` | `mm` | Input unit (`mm` or `m`); converted internally to metres |
| `--periods` | `0.25,0.5,1.0,2.0` | Candidate periods in years |
| `--auto-periods` | `5` | Auto-detect up to N additional dominant periods (Lomb-Scargle) |
| `--poly-deg` | `1` | Polynomial degree: 0=offset, 1=velocity, 2=acceleration. `-1` = auto-select best from `[poly-deg-min, 3]` |
| `--poly-deg-min` | `0` | Minimum degree when auto-selecting |
| `--sigma-min` | `2.0` | Minimum a-priori sigma (mm) |
| `--sigma-max` | `15.0` | Maximum a-priori sigma (mm) |
| `--sigma-step` | `0.5` | Sigma scan step (mm) |
| `--alpha` | `0.05` | OMT significance level |
| `--max-iter` | `5` | Max DIA iterations per sigma |
| `--jumps` | — | Extra jump dates (comma-separated, YYYY-MM-DD) |
| `--polylines` | — | Extra polyline break dates (comma-separated) |
| `--logs` | — | Extra log relaxation dates with tau (DATE:tau) |
| `--output-dir` | *csv parent* | Parent directory for output folder |
| `--cores` | `1` | CPU cores for parallel sigma scan |
| `--force` | — | Overwrite existing results (default: skip processed components) |

**Hypothesis group control flags:**

| Flag | Effect | Physical use case |
|------|--------|-------------------|
| `--no-jump` | Disable automatic jump detection (hyp group 1). Manual `--jumps` still honoured. | MLCW data with no physical mechanism for instantaneous jumps — all steps are velocity changes, not offsets. |
| `--no-seasonal` | Disable all periodic/seasonal detection (hyp group 2). Overrides `--periods` and `--auto-periods`. | When you want only the trend/segmentation without seasonal cycles. |
| `--no-exp-trend` | Disable exponential trend detection (hyp group 4). | Prevent spurious exponential fits in short records. |
| `--auto-sigma` | Select sigma that maximises polyline breakpoints (finest velocity segmentation) instead of most parsimonious model. | Pure piecewise-linear fitting where resolving every velocity change matters more than model simplicity. |

**Date range filters:**

| Flag | Effect |
|------|--------|
| `--start-date YYYY-MM-DD` | Clip data to start at this date (inclusive) |
| `--end-date YYYY-MM-DD` | Clip data to end at this date (inclusive) |

**Exponential trend:**

| Flag | Effect |
|------|--------|
| `--exp-trend auto` | Auto-detect best exponential decay rate b |
| `--exp-trend 0.001` | Fixed b = 0.001 /day (~2.7 yr time constant) |

### `reconstruct` — forward model from parameters

```bash
appsigsolv reconstruct <input_file> --json model.json [options]
```

Builds the synthetic signal from a saved model JSON. Supports custom sampling rates (`daily` or `custom` days-of-month) for comparison with irregular observations.

---

## Outputs

Each component generates these files inside `{output-dir}/{csv-stem}/`:

| File | Contents |
|------|----------|
| `*_decomposed_{comp}.csv` | Full decomposition: original, periodic, trend, step, polyline, exp, log, residual, w-test values, flags |
| `*_model_{comp}.json` | Best-fit model parameters (JSON) — load with `reconstruct --json` |
| `*_report_{comp}.md` | Human-readable report: sigma scan table, OMT statistics, accepted parameters |
| `*_{comp}.png` | Decomposition plot (4-panel: data+model, periodic, trend, residuals) |
| `*_skipped_{comp}.txt` | Skip report (if fitting failed or timed out) |

All displacements in the CSV are written in **mm** for readability; internal computation uses metres.

---

## Architecture

```
appsigsolv/
├── cli/
│   ├── parser.py             — Argparse configuration
│   ├── cmd_decompose.py      — Decompose command logic
│   └── cmd_reconstruct.py    — Reconstruct command logic
├── core/
│   ├── dia.py                — DIA algorithm (hyp group competition, sigma scan)
│   └── modeling.py           — Design matrices, parameter estimation, component extraction
├── io/
│   └── data_manager.py       — Loading, preprocessing (outlier removal, gap-fill), saving
└── utils/
    └── visualization.py      — Plotting and report generation
```

### The DIA hypothesis groups in detail

**Group 1 — Outlier detection (datasnooping):** Computes the w-test statistic for every epoch. The epoch with the largest |w| > 3.29 (99.9 % confidence) is flagged as an outlier and its influence is removed by inserting a unit-step at that date. *Disabled by `--no-jump`.*

**Group 2 — Missing periodic signal:** Runs a Lomb-Scargle periodogram on the residuals. The top-5 spectral peaks are evaluated; the strongest peak not already in the model is added as a new sine-cosine pair. *Disabled by `--no-seasonal`; restricted to the candidate set via `allowed_periods`.*

**Group 3 — Velocity break (polyline):** A CUSUM accumulation of residuals detects sustained drift. When the cumulative sum exceeds a spike threshold, a hinge (change in slope) is inserted at the CUSUM pivot point.

**Group 4 — Exponential trend:** Fits an exponential model to the residuals and evaluates its AIC improvement. If the AIC reduction is significant, the exponential component is added.

After each addition, the OMT is re-evaluated. The loop terminates when the OMT accepts the model (p ≥ α) or `--max-iter` is exceeded.

---

## Typical workflows

### GPS time series (daily, 10+ years)

```bash
appsigsolv decompose gps.csv --component dU --periods "0.5,1.0" \
  --poly-deg 2 --sigma-min 1.0 --sigma-max 8.0 --cores 4 --output-dir ./gps_results
```

Includes annual loading cycles, allows acceleration, uses tighter sigma range.

### MLCW (multilevel compaction, irregular, many components)

```bash
appsigsolv decompose mlcw.csv --component all --periods "0.5,1.0" \
  --poly-deg 1 --no-jump --sigma-min 2.0 --sigma-max 12.0 \
  --sigma-step 1.0 --output-dir ./mlcw_results
```

Each ring is a separate aquifer layer; `--no-jump` prevents spurious step offsets, `--poly-deg 1` captures steady compaction velocity between breakpoints.

### Piecewise-linear segmentation

```bash
appsigsolv decompose mlcw.csv --component all --no-seasonal --no-jump \
  --auto-sigma --poly-deg 1 --sigma-min 1.0 --sigma-max 10.0 --sigma-step 0.5
```

Use `--auto-sigma` to resolve the finest velocity segmentation. Each polyline segment represents a period of constant compaction rate.

---

## References

- Teunissen, P. J. G. (2017). *Distributional theory for the DIA method*. Journal of Geodesy.
- Teunissen, P. J. G. (2018). *The DIA method for testing and estimation — a review*. Journal of Geodesy.
- Amiri-Simkooei, A. R., Tiberius, C. C. J. M., & Teunissen, P. J. G. (2007). *Assessment of noise in GPS coordinate time series: methodology and results*. Journal of Geophysical Research.

---

## License

MIT — see [LICENSE](LICENSE).

---

## Related repositories

- [twostoolspy](https://github.com/david-ncu2019/twostoolspy) — 2S-TOOL Python port for InSAR time series analysis
- [timeseries_signal_solver](https://github.com/david-ncu2019/timeseries_signal_solver) — Earlier versions of the solver
