"""Data Management (Loading, Preprocessing, Saving)."""
import os
import json
import pandas as pd
import numpy as np
from scipy.stats import median_abs_deviation
from pathlib import Path

def load_and_preprocess(csv_path: str, component: str, date_col: str = "", unit: str = "mm", mad_threshold: float = 4.5, irregular: bool = False) -> tuple:
    df = pd.read_csv(csv_path)
    # Normalise numeric column names to 3 d.p. so float-precision artefacts
    # like '86.15899999999999' match the rounded key '86.159'.
    rename_map = {}
    for c in df.columns:
        try:
            rename_map[c] = f"{round(float(c), 3):g}"
        except ValueError:
            pass
    if rename_map:
        df = df.rename(columns=rename_map)
    
    detected_date_col = date_col
    if not detected_date_col:
        for c in ['gpsdate', 'datetime', 'date', 'Date', 'time', 'Time']:
            if c in df.columns:
                detected_date_col = c
                break
        if not detected_date_col:
            detected_date_col = df.columns[0]
            print(f"  [preprocess] Warning: No standard date column found. Defaulting to first column '{detected_date_col}'")
            
    df["_internal_date"] = pd.to_datetime(df[detected_date_col].astype(str), errors='coerce')
    df = df.dropna(subset=["_internal_date", component])
    df = df.set_index("_internal_date").sort_index()

    series = df[component].copy()
    series = pd.to_numeric(series, errors='coerce').dropna()
    series = series.groupby(series.index).median()
    
    if unit == 'mm':
        series = series / 1000.0  

    # Detrend before MAD: avoids flagging the long-term trend itself as outliers
    s_clean = series.dropna()
    if len(s_clean) >= 2:
        days_num = np.array([(d - s_clean.index[0]).days for d in s_clean.index], dtype=float)
        linear_fit = np.polyfit(days_num, s_clean.values, 1)
        trend_vals = np.polyval(linear_fit, days_num)
        residuals_for_mad = s_clean.values - trend_vals
        mad = median_abs_deviation(residuals_for_mad, scale="normal")
        days_full = np.array([(d - s_clean.index[0]).days for d in series.index], dtype=float)
        trend_full = np.polyval(linear_fit, days_full)
        detrended_full = series.values - trend_full
        outliers = np.abs(detrended_full) > mad_threshold * mad
    else:
        median = series.median()
        mad = median_abs_deviation(s_clean, scale="normal")
        outliers = np.abs(series - median) > mad_threshold * mad
    series[outliers] = np.nan
    print(f"  [preprocess] {component}: {outliers.sum()} outliers removed")

    if not irregular:
        full_dates = pd.date_range(start=series.index.min(), end=series.index.max(), freq="D")
        series = series.reindex(full_dates)
        series = series.interpolate(method="time", limit=7).bfill(limit=7).ffill(limit=7)
        print(f"  [preprocess] {component}: {len(series)} daily points after gap-fill")
    else:
        print(f"  [preprocess] {component}: Kept {len(series)} original irregular points (no daily gap-fill)")

    return series, detected_date_col

def load_and_clean_data_for_reconstruct(file_path, date_col=None, target_col=None, unit='mm'):
    print(f"Loading data from {file_path}...")
    if file_path.endswith('.csv'):
        df = pd.read_csv(file_path)
    elif file_path.endswith(('.xls', '.xlsx')):
        df = pd.read_excel(file_path)
    else:
        raise ValueError("Unsupported file format. Must be .csv or .xls/.xlsx")

    if not date_col:
        for c in ['gpsdate', 'date', 'Date', 'time', 'Time', 'datetime']:
            if c in df.columns:
                date_col = c
                break
    if not target_col:
        for c in ['dU', 'dN', 'dE', 'displacement_mm', 'disp', 'value', 'target', 'displacement']:
            if c in df.columns:
                target_col = c
                break
                
    if not date_col or not target_col:
        raise ValueError(f"Could not auto-detect columns in {df.columns}. Specify with --date-col/--target-col.")

    print(f"  Using columns: date='{date_col}', target='{target_col}'")
    
    df[date_col] = pd.to_datetime(df[date_col].astype(str), errors='coerce')
    df = df.dropna(subset=[date_col, target_col]).sort_values(by=date_col)
    
    dates = [d.to_pydatetime() for d in df[date_col]]
    raw_disp = df[target_col].values.astype(np.float64)
    
    if unit == 'mm':
        raw_disp *= 0.001
        
    return df, date_col, target_col, dates, raw_disp

def load_json_config(json_path):
    if not os.path.exists(json_path):
        return {}
    with open(json_path, 'r') as f:
        config = json.load(f)
    return config

def save_json_config(best_model: dict, comp: str, out_root: Path, stem: str) -> Path:
    model_config = {k: v for k, v in best_model.items() if not k.startswith("_")}
    json_path = out_root / f"{stem}_model_{comp}.json"
    with open(json_path, "w") as f:
        json.dump(model_config, f, indent=4)
    print(f"  [output] JSON config saved: {json_path}")
    return json_path

def save_csv(series: pd.Series, components: dict, comp: str, out_root: Path, stem: str) -> Path:
    df = pd.DataFrame({"date": series.index, comp: series.values})
    df = df.set_index("date")
    for col, s in components.items():
        df[col] = s.values

    df["flagged"] = np.abs(df[f"{comp}_wtest"]) > 3.29
    # Convert back to mm for output (internal computation is in metres)
    skip_cols = {f"{comp}_wtest", "flagged"}
    mm_cols = [c for c in df.columns if c not in skip_cols]
    df[mm_cols] = df[mm_cols] * 1000.0
    csv_path = out_root / f"{stem}_decomposed_{comp}.csv"
    df.to_csv(csv_path, float_format="%.6f")
    print(f"  [output] CSV saved: {csv_path}")
    return csv_path
