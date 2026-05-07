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

def prescreen_periods(series: pd.Series, candidates_yr: list, cond_threshold: float = 1e8) -> list:
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

    # Baseline design matrix: normalised quadratic polynomial at valid observation days
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

        # Collinearity guard: reject if adding this period raises condition number too high
        cos_col = np.cos(2 * np.pi / period_days * t)
        sin_col = np.sin(2 * np.pi / period_days * t)
        extra_cols = [
            col for p in accepted
            for col in (np.cos(2*np.pi/(p*365.25)*t), np.sin(2*np.pi/(p*365.25)*t))
        ]
        G_trial = np.column_stack([G_base] + extra_cols + [cos_col, sin_col])
        if np.linalg.cond(G_trial) > cond_threshold:
            continue

        accepted.append(period_yr)
    return accepted

def _cusum_break_date(residuals: np.ndarray, dates: list, min_segment: int = 10):
    """Return the date string of the max-CUSUM structural break, or None if not significant."""
    n = len(residuals)
    if n < 2 * min_segment:
        return None
    mu = np.mean(residuals)
    cusum = np.cumsum(residuals - mu)
    sigma = np.std(residuals, ddof=1)
    if sigma < 1e-12:
        return None
    # Normalised CUSUM range; 1.36 is the ~95% KS critical value
    cusum_range = (np.max(cusum) - np.min(cusum)) / (sigma * np.sqrt(n))
    if cusum_range < 1.36:
        return None
    break_idx = int(np.argmax(np.abs(cusum)))
    if break_idx < min_segment or break_idx > n - min_segment:
        return None
    return dates[break_idx].strftime("%Y%m%d")

def auto_detect_exp_trend(series: pd.Series, b_candidates=None, aic_improvement_threshold: float = 10.0):
    """
    Test whether series is better described by [1, t, exp(-b*t)-1] than [1, t].
    Returns best b_per_day (1/days) or None if no significant improvement.

    Two guards must both pass:
    - ΔAIC > aic_improvement_threshold (default 10.0 = strong evidence, Burnham & Anderson)
    - Exp-trend component explains >= 10% of series variance
    """
    if b_candidates is None:
        b_candidates = np.logspace(-4, -1, 20)

    s = series.dropna()
    if len(s) < 10:
        return None

    t0 = s.index[0]
    days = np.array([(d - t0).days for d in s.index], dtype=np.float64)
    y = s.values.astype(np.float64)
    n = len(y)

    G_base = np.column_stack([np.ones(n), days])
    m_base, ssr_base, _, _ = np.linalg.lstsq(G_base, y, rcond=None)
    ssr_base = float(ssr_base[0]) if len(ssr_base) > 0 else float(np.sum((y - G_base @ m_base) ** 2))
    aic_base = n * np.log(max(ssr_base / n, 1e-30)) + 2 * 2

    best_b = None
    best_aic = aic_base - aic_improvement_threshold
    best_m_trial = None
    best_exp_col = None

    for b in b_candidates:
        exp_col = np.exp(-b * days) - 1.0
        G_trial = np.column_stack([np.ones(n), days, exp_col])
        if np.linalg.cond(G_trial) > 1e10:
            continue
        m_trial, ssr_trial, _, _ = np.linalg.lstsq(G_trial, y, rcond=None)
        ssr_trial = float(ssr_trial[0]) if len(ssr_trial) > 0 else float(np.sum((y - G_trial @ m_trial) ** 2))
        aic_trial = n * np.log(max(ssr_trial / n, 1e-30)) + 2 * 3
        if aic_trial < best_aic:
            best_aic = aic_trial
            best_b = b
            best_m_trial = m_trial
            best_exp_col = exp_col

    if best_b is None:
        return None

    # Variance guard: exp-trend component must explain >= 10% of series variance
    exp_component = best_exp_col * best_m_trial[2]
    if np.var(exp_component) / max(np.var(y), 1e-30) < 0.10:
        return None

    return best_b

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

    # CUSUM-based detection for gradual accelerations not caught by velocity spike test
    if len(residuals) >= 20:
        cusum_break = _cusum_break_date(residuals, dates)
        if cusum_break is not None and cusum_break not in model.get("polyline", []):
            return "polyline", cusum_break

    # Exp-trend detection: monotonic residual with decreasing rate suggests missed exp trend
    if len(residuals) >= 20 and model.get("exp_trend") is None:
        residual_series = pd.Series(residuals, index=pd.DatetimeIndex(dates))
        best_b = auto_detect_exp_trend(residual_series)
        if best_b is not None:
            return "exp_trend", best_b

    return "None", None

def run_omt_dia_loop(series, jump_dates, initial_periods, initial_polylines, initial_logs, poly_deg, sigma_mm, alpha, max_iter, exp_trend_b=None):
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
        "exp_trend": exp_trend_b,
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

        if adapt_type == "exp_trend":
            if model.get("exp_trend") is None:
                model["exp_trend"] = adapt_val
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

def _run_single_sigma(sigma_mm, series, jump_dates, candidate_periods, initial_polylines, initial_logs, poly_deg, alpha, max_iter, exp_trend_b=None):
    result = run_omt_dia_loop(series, jump_dates, candidate_periods, initial_polylines, initial_logs, poly_deg, sigma_mm, alpha, max_iter, exp_trend_b=exp_trend_b)
    if result is not None:
        s = result["_omt_stats"]
        return result, {"sigma_mm": sigma_mm, "accepted": True, "p_value": s["p_value"],
                        "n_param": s["n_param"], "n_periods": len(result.get("periodic", [])),
                        "n_polylines": len(result.get("polyline", []))}
    return None, {"sigma_mm": sigma_mm, "accepted": False,
                  "p_value": None, "n_param": None,
                  "n_periods": None, "n_polylines": None}

def run_omt_sigma_scan(series, jump_dates, candidate_periods, initial_polylines, initial_logs, poly_deg, sigma_min, sigma_max, sigma_step, alpha, max_iter, no_relax=False, cores=1, exp_trend_b=None):
    sigmas = np.arange(sigma_min, sigma_max + sigma_step * 0.5, sigma_step)
    scan_results = []
    scan_table = []

    if cores > 1:
        import concurrent.futures
        from functools import partial
        worker = partial(_run_single_sigma, series=series, jump_dates=jump_dates, candidate_periods=candidate_periods, initial_polylines=initial_polylines, initial_logs=initial_logs, poly_deg=poly_deg, alpha=alpha, max_iter=max_iter, exp_trend_b=exp_trend_b)
        with concurrent.futures.ProcessPoolExecutor(max_workers=cores) as executor:
            for result, table_row in executor.map(worker, sigmas):
                scan_table.append(table_row)
                if result is not None:
                    scan_results.append(result)
    else:
        for sigma_mm in sigmas:
            result, table_row = _run_single_sigma(sigma_mm, series, jump_dates, candidate_periods, initial_polylines, initial_logs, poly_deg, alpha, max_iter, exp_trend_b=exp_trend_b)
            scan_table.append(table_row)
            if result is not None:
                scan_results.append(result)

    if not scan_results:
        return None, scan_table

    best = min(scan_results, key=lambda m: (m["_omt_stats"]["n_param"], -m["_omt_stats"]["sigma_mm"], -m["_omt_stats"]["p_value"]))
    if not no_relax:
        best = test_relaxation(series, jump_dates, best, alpha)
    return best, scan_table
