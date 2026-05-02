"""DIA logic, anomaly detection, and OMT loops."""
from datetime import datetime
import numpy as np
import pandas as pd
from scipy.stats import chi2
from scipy.signal import lombscargle, find_peaks
from .modeling import estimate_time_func

def calculate_omt(residuals, m_obs, n_param, sigma_m, alpha=0.05):
    r = m_obs - n_param
    if r <= 0:
        return np.inf, np.inf, 0.0, 0.0
    
    ssr = np.sum(residuals**2)
    T_stat = ssr / (sigma_m**2)
    omt = T_stat / r
    
    p_value = 1.0 - chi2.cdf(T_stat, df=r)
    K = chi2.ppf(1.0 - alpha, df=r)
    K_norm = K / r
    
    return T_stat, omt, p_value, K_norm

def detect_jumps(
    series: pd.Series,
    extra_dates: list,
    window_days: int = 365,
    sigma_threshold: float = 3.0,
    smooth_window: int = 30,
    adaptive_percentile: float = 99,
    min_days_apart: int = 90,
) -> list:
    from scipy.ndimage import median_filter

    valid = series.interpolate(method="time")
    diffs = np.diff(valid.values)
    abs_diffs = np.abs(diffs)

    rolling_median = pd.Series(abs_diffs).rolling(window=window_days, center=True, min_periods=1).median().values
    rolling_mad = pd.Series(abs_diffs).rolling(window=window_days, center=True, min_periods=1).apply(
        lambda x: np.median(np.abs(x - np.median(x))), raw=True
    ).values
    threshold = rolling_median + sigma_threshold * 1.4826 * rolling_mad
    candidates = [(valid.index[i + 1], diffs[i]) for i in np.where(abs_diffs > threshold)[0]]

    raw_diffs = np.abs(np.diff(series.values))
    raw_diffs_clean = raw_diffs[~np.isnan(raw_diffs)]
    adaptive_thr = max(np.percentile(raw_diffs_clean, adaptive_percentile) * 1000, 3.0) if len(raw_diffs_clean) > 0 else 3.0

    trend = median_filter(valid.values, size=smooth_window)
    trend_s = pd.Series(trend, index=valid.index)

    validated = []
    for date, _ in candidates:
        try:
            idx = trend_s.index.get_loc(date)
            before_vals = trend_s.values[max(0, idx - 30):idx]
            after_vals = trend_s.values[idx + 1:min(len(trend_s), idx + 31)]
            if len(before_vals) == 0 or len(after_vals) == 0:
                continue
            before = np.median(before_vals)
            after = np.median(after_vals)
            if abs((after - before) * 1000) >= adaptive_thr:
                validated.append(date)
        except Exception:
            continue

    validated.sort()
    filtered = []
    for d in validated:
        if not any(abs((d - f).days) < min_days_apart for f in filtered):
            filtered.append(d)

    for ds in extra_dates:
        try:
            dt = pd.Timestamp(ds)
            if not any(abs((dt - f).days) < min_days_apart for f in filtered):
                filtered.append(dt)
        except Exception:
            pass

    filtered.sort()
    return [d.to_pydatetime() for d in filtered]

def auto_detect_periods(series: pd.Series, max_periods: int = 5, min_yr: float = 0.2, max_yr: float = 20.0) -> list:
    values = series.values
    valid_mask = ~np.isnan(values)
    if valid_mask.sum() < 20:
        return []
        
    days = np.array([(d - series.index[0]).days for d in series.index])
    coeffs = np.polyfit(days[valid_mask], values[valid_mask], 2)
    detrended = values[valid_mask] - np.polyval(coeffs, days[valid_mask])
    
    freqs = np.linspace(2*np.pi/(max_yr*365.25), 2*np.pi/(min_yr*365.25), 5000)
    pgram = lombscargle(days[valid_mask], detrended, freqs, normalize=True)
    peaks, _ = find_peaks(pgram, distance=50)
    
    valid_peaks = []
    for p in peaks:
        f_rad = freqs[p]
        period_yr = (2*np.pi / f_rad) / 365.25
        valid_peaks.append((period_yr, pgram[p]))
                
    valid_peaks.sort(key=lambda x: x[1], reverse=True)
    return [round(p[0], 2) for p in valid_peaks[:max_periods]]

def prescreen_periods(series: pd.Series, candidates_yr: list) -> list:
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
    
    for period_yr in candidates_yr:
        period_days = period_yr * 365.25
        target_f = 2*np.pi / period_days
        idx = np.argmin(np.abs(freqs - target_f))
        window = pgram[max(0, idx-50):min(len(pgram), idx+51)]
        has_peak = len(window) > 0 and window.max() > power_threshold and window.max() > 0.05

        if has_peak:
            accepted.append(period_yr)
    return accepted

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
            
    return "None", None

def run_omt_dia_loop(series, jump_dates, initial_periods, initial_polylines, initial_logs, poly_deg, sigma_mm, alpha, max_iter):
    series_clean = series.dropna()
    if series_clean.empty:
        return None
    date_list = [d.to_pydatetime() if hasattr(d, "to_pydatetime") else d for d in series_clean.index]
    dis_ts = series_clean.values.copy()
    sigma_m = sigma_mm / 1000.0

    step_dates = [d.strftime("%Y%m%d") for d in jump_dates]
    model = {
        "polynomial": poly_deg,
        "periodic": list(initial_periods),
        "stepDate": step_dates,
        "polyline": list(initial_polylines),
        "exp": {},
        "log": initial_logs,
    }

    last_omt = 9999.0
    for iteration in range(max_iter):
        G, m, e2, d_hat = estimate_time_func(model, date_list, dis_ts)
        residuals = dis_ts - d_hat
        n_param = G.shape[1]
        T_stat, omt, p_value, K_norm = calculate_omt(residuals, len(dis_ts), n_param, sigma_m, alpha)

        if p_value >= alpha:
            model["_omt_stats"] = {
                "sigma_mm": sigma_mm, "omt": omt, "p_value": p_value,
                "K_norm": K_norm, "n_param": n_param, "iterations": iteration,
            }
            return model
            
        if omt >= last_omt and iteration > 0:
            break
        last_omt = omt

        adapt_type, adapt_val = robust_analyze_residuals(residuals, date_list, model)

        if adapt_type == "period":
            if adapt_val not in model["periodic"]:
                model["periodic"].append(adapt_val)
            else:
                adapt_type = "polyline"
                adapt_val = None
                
        if adapt_type == "polyline":
            if adapt_val is None:
                days = np.array([(d - date_list[0]).days for d in date_list])
                velocity = np.diff(residuals) / np.diff(days)
                ext_idx = np.argmax(np.abs(velocity))
                adapt_val = date_list[ext_idx].strftime("%Y%m%d")
            if adapt_val not in model["polyline"]:
                model["polyline"].append(adapt_val)
            else:
                break
                
        if adapt_type == "None":
            break
    return None

def test_relaxation(series, jump_dates, accepted_model, alpha):
    if not jump_dates: return accepted_model
    series_clean = series.dropna()
    if series_clean.empty: return accepted_model

    date_list = [d.to_pydatetime() if hasattr(d, "to_pydatetime") else d for d in series_clean.index]
    dis_ts = series_clean.values.copy()
    sigma_mm = accepted_model["_omt_stats"]["sigma_mm"]
    sigma_m = sigma_mm / 1000.0

    model = {k: v for k, v in accepted_model.items() if not k.startswith("_")}
    model["exp"] = {}

    current_stats = accepted_model["_omt_stats"]
    best_omt = current_stats["omt"]

    for jump_dt in jump_dates:
        jump_str = jump_dt.strftime("%Y%m%d")
        best_tau = None
        for tau_days in [30, 90, 180]:
            trial_model = {k: (list(v) if isinstance(v, list) else v) for k, v in model.items()}
            trial_exp = {k: list(v) for k, v in model["exp"].items()}
            trial_exp[jump_str] = [tau_days]
            trial_model["exp"] = trial_exp

            G, m, e2, d_hat = estimate_time_func(trial_model, date_list, dis_ts)
            residuals = dis_ts - d_hat
            n_param = G.shape[1]
            _, omt, p_value, _ = calculate_omt(residuals, len(dis_ts), n_param, sigma_m, alpha)

            if p_value >= alpha and omt < best_omt:
                best_omt = omt
                best_tau = tau_days

        if best_tau is not None:
            model["exp"][jump_str] = [best_tau]

    G, m, e2, d_hat = estimate_time_func(model, date_list, dis_ts)
    residuals = dis_ts - d_hat
    n_param = G.shape[1]
    T_stat, omt, p_value, K_norm = calculate_omt(residuals, len(dis_ts), n_param, sigma_m, alpha)
    model["_omt_stats"] = {
        "sigma_mm": sigma_mm, "omt": omt, "p_value": p_value,
        "K_norm": K_norm, "n_param": n_param, "iterations": current_stats["iterations"],
    }
    return model

def _run_single_sigma(sigma_mm, series, jump_dates, candidate_periods, initial_polylines, initial_logs, poly_deg, alpha, max_iter):
    result = run_omt_dia_loop(series, jump_dates, candidate_periods, initial_polylines, initial_logs, poly_deg, sigma_mm, alpha, max_iter)
    if result is not None:
        s = result["_omt_stats"]
        return result, {"sigma_mm": sigma_mm, "accepted": True, "p_value": s["p_value"],
                        "n_param": s["n_param"], "n_periods": len(result.get("periodic", [])),
                        "n_polylines": len(result.get("polyline", []))}
    return None, {"sigma_mm": sigma_mm, "accepted": False,
                  "p_value": None, "n_param": None,
                  "n_periods": None, "n_polylines": None}

def run_omt_sigma_scan(series, jump_dates, candidate_periods, initial_polylines, initial_logs, poly_deg, sigma_min, sigma_max, sigma_step, alpha, max_iter, no_relax=False, cores=1):
    sigmas = np.arange(sigma_min, sigma_max + sigma_step * 0.5, sigma_step)
    scan_results = []
    scan_table = []

    if cores > 1:
        import concurrent.futures
        from functools import partial
        worker = partial(_run_single_sigma, series=series, jump_dates=jump_dates, candidate_periods=candidate_periods, initial_polylines=initial_polylines, initial_logs=initial_logs, poly_deg=poly_deg, alpha=alpha, max_iter=max_iter)
        with concurrent.futures.ProcessPoolExecutor(max_workers=cores) as executor:
            for result, table_row in executor.map(worker, sigmas):
                scan_table.append(table_row)
                if result is not None:
                    scan_results.append(result)
    else:
        for sigma_mm in sigmas:
            result, table_row = _run_single_sigma(sigma_mm, series, jump_dates, candidate_periods, initial_polylines, initial_logs, poly_deg, alpha, max_iter)
            scan_table.append(table_row)
            if result is not None:
                scan_results.append(result)

    if not scan_results:
        return None, scan_table

    best = min(scan_results, key=lambda m: (m["_omt_stats"]["sigma_mm"], m["_omt_stats"]["n_param"], -m["_omt_stats"]["p_value"]))
    if not no_relax:
        best = test_relaxation(series, jump_dates, best, alpha)
    return best, scan_table
