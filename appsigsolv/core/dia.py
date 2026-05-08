"""DIA logic, anomaly detection, and OMT loops.

Implements the Teunissen Detection-Identification-Adaptation (DIA) framework:
  - Detection:       Overall Model Test  T = SSR/σ² ~ χ²(r),  accept if p ≥ α
  - Identification:  Formal w-test per alternative hypothesis (Baarda); select max|w|
  - Adaptation:      Extend functional model with identified component; re-estimate
"""
from datetime import datetime
import numpy as np
import pandas as pd
from scipy.stats import chi2, norm
from scipy.signal import lombscargle, find_peaks
from .modeling import estimate_time_func


# ---------------------------------------------------------------------------
# Core OMT statistics
# ---------------------------------------------------------------------------

def calculate_omt(residuals, m_obs, n_param, sigma_m, alpha=0.05):
    """
    Compute Overall Model Test statistics.

    Returns
    -------
    T_stat        : chi-squared test statistic  T = SSR / σ²  ~ χ²(r)
    unit_var_factor : a-posteriori unit variance factor  T/r = σ̂²/σ₀²
    p_value       : Prob(χ²(r) ≥ T_stat);  accept model when p_value ≥ alpha
    chi2_critical : critical value K = χ²_{1-α}(r)  (accept when T_stat ≤ K)
    sigma_hat_m   : a-posteriori noise estimate  σ̂ = sqrt(SSR/r)  [metres]
    """
    r = m_obs - n_param
    if r <= 0:
        return np.inf, np.inf, 0.0, np.inf, sigma_m

    ssr = np.sum(residuals**2)
    T_stat = ssr / (sigma_m**2)
    unit_var_factor = T_stat / r          # σ̂²/σ₀²; equals 1 when model fits perfectly

    p_value = 1.0 - chi2.cdf(T_stat, df=r)
    chi2_critical = chi2.ppf(1.0 - alpha, df=r)  # K; accept when T_stat ≤ K
    sigma_hat_m = np.sqrt(ssr / r)               # a-posteriori noise [metres]

    return T_stat, unit_var_factor, p_value, chi2_critical, sigma_hat_m


# ---------------------------------------------------------------------------
# Jump / outlier detection (pre-processing, not formal DIA identification)
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Period detection helpers (used to build alternative-hypothesis library)
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# CUSUM structural-break helper (used to build hypothesis candidates)
# ---------------------------------------------------------------------------

def _cusum_break_date(residuals: np.ndarray, dates: list, min_segment: int = 10):
    """Return the date string of the max-CUSUM structural break, or None if not significant.

    Uses the Ploberger-Krämer (1992) normalised CUSUM range test.
    The critical boundary 1.36 corresponds to the ~95% level for this test
    (not to be confused with the Kolmogorov-Smirnov statistic despite the
    same numerical value at α=0.05).
    """
    n = len(residuals)
    if n < 2 * min_segment:
        return None
    mu = np.mean(residuals)
    cusum = np.cumsum(residuals - mu)
    sigma = np.std(residuals, ddof=1)
    if sigma < 1e-12:
        return None
    # Ploberger-Krämer (1992) boundary: reject stability if range > 1.36·σ·√n
    cusum_range = (np.max(cusum) - np.min(cusum)) / (sigma * np.sqrt(n))
    if cusum_range < 1.36:
        return None
    break_idx = int(np.argmax(np.abs(cusum)))
    if break_idx < min_segment or break_idx > n - min_segment:
        return None
    return dates[break_idx].strftime("%Y%m%d")


# ---------------------------------------------------------------------------
# Exponential-trend auto-detection (AIC-based, used to build hypothesis library)
# ---------------------------------------------------------------------------

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

    exp_component = best_exp_col * best_m_trial[2]
    if np.var(exp_component) / max(np.var(y), 1e-30) < 0.10:
        return None

    return best_b


# ---------------------------------------------------------------------------
# Formal w-test (Baarda datasnooping — Identification phase)
# ---------------------------------------------------------------------------

def _compute_w_stats(residuals: np.ndarray, G: np.ndarray, sigma_m: float):
    """
    Compute per-epoch Baarda w-statistics for datasnooping.

    Under Q_yy = σ²I the residual covariance is Q_ê = σ²(I - H) where H = G(G^TG)^{-1}G^T.
    The w-statistic for epoch i simplifies to:
        w_i = ê_i / (σ · sqrt(1 - h_ii))
    where h_ii is the i-th diagonal of the hat matrix H.

    Under H₀, w_i ~ N(0,1).  Flag epoch i if |w_i| > z_{α/2}.
    """
    n = len(residuals)
    try:
        # Hat matrix diagonal via thin QR (numerically stable, O(n·p²))
        Q, _ = np.linalg.qr(G, mode='reduced')
        h_diag = np.sum(Q**2, axis=1)
    except np.linalg.LinAlgError:
        h_diag = np.zeros(n)

    denom = sigma_m * np.sqrt(np.maximum(1.0 - h_diag, 1e-12))
    return residuals / denom


def _w_test_for_candidate(residuals: np.ndarray, G: np.ndarray, sigma_m: float,
                           c_col: np.ndarray) -> float:
    """
    Compute the w-statistic for a single alternative hypothesis defined by
    column vector c_col (the additional design-matrix column for the candidate).

    Formula (Teunissen):
        w = (c^T Q_y^{-1} ê) / sqrt(c^T Q_y^{-1} Q_ê Q_y^{-1} c)

    Under Q_yy = σ²I this simplifies to:
        numerator   = c^T ê / σ²
        Q_ê = σ²(I - H),  so  c^T Q_y^{-1} Q_ê Q_y^{-1} c = c^T(I-H)c / σ²
        w = (c^T ê) / (σ · sqrt(c^T(I-H)c))
    """
    try:
        Q, _ = np.linalg.qr(G, mode='reduced')
        Hc = Q @ (Q.T @ c_col)           # H·c  (projection of c onto column space of G)
        ImHc = c_col - Hc                # (I-H)·c
        denom_sq = np.dot(c_col, ImHc)   # c^T(I-H)c
        if denom_sq < 1e-20:
            return 0.0
        w = np.dot(c_col, residuals) / (sigma_m * np.sqrt(denom_sq))
    except (np.linalg.LinAlgError, ValueError):
        w = 0.0
    return float(w)


# ---------------------------------------------------------------------------
# Identification: build candidate hypotheses and select via formal w-test
# ---------------------------------------------------------------------------

def _identify_best_alternative(residuals: np.ndarray, dates: list, model: dict,
                                G: np.ndarray, sigma_m: float, alpha: float,
                                days: np.ndarray = None, freqs: np.ndarray = None):
    """
    Teunissen Identification phase: evaluate all candidate alternative
    hypotheses with formal w-tests and return the one with max|w|.

    Candidate hypotheses
    --------------------
    1. Datasnooping   — unit vector c_i for each epoch (single-outlier test)
    2. Missing period — cosine/sine columns for each Lomb-Scargle peak
    3. Velocity break — piecewise-linear column at the estimated break date
    4. Exp trend      — exponential-trend column (if not already in model)

    Selection rule (Teunissen): adapt on candidate i* = argmax_i |w_i|,
    provided |w_{i*}| > z_{α/2}  (two-tailed normal critical value).

    Parameters
    ----------
    days  : pre-computed day-offset array (loop-invariant); computed here if None
    freqs : pre-computed Lomb-Scargle frequency grid (loop-invariant); computed here if None

    Returns
    -------
    (adapt_type, adapt_val, max_w_abs)
        adapt_type : "period" | "polyline" | "exp_trend" | "outlier" | "None"
        adapt_val  : the value to add to the model (period_yr, date_str, b_val, or epoch_idx)
        max_w_abs  : absolute w-statistic of winning hypothesis
    """
    if days is None:
        days = np.array([(d - dates[0]).days for d in dates], dtype=float)
    n = len(residuals)
    z_crit = norm.ppf(1.0 - alpha / 2.0)   # two-tailed critical value

    best_w = 0.0
    best_type = "None"
    best_val = None

    # --- Hypothesis group 1: datasnooping (single-epoch outlier) ------------
    w_epoch = _compute_w_stats(residuals, G, sigma_m)
    max_epoch_idx = int(np.argmax(np.abs(w_epoch)))
    max_epoch_w = float(np.abs(w_epoch[max_epoch_idx]))
    if max_epoch_w > best_w:
        best_w = max_epoch_w
        best_type = "outlier"
        best_val = max_epoch_idx

    # --- Hypothesis group 2: missing periodic signal ------------------------
    if n >= 20:
        if freqs is None:
            freqs = np.linspace(2*np.pi/(20*365.25), 2*np.pi/(0.2*365.25), 500)
        pgram = lombscargle(days, residuals, freqs, normalize=True)
        max_pwr = np.max(pgram)
        peaks, _ = find_peaks(pgram, height=max_pwr * 0.4)

        if len(peaks) > 0 and max_pwr > 0.1:
            peak_powers = pgram[peaks]
            sorted_peaks = peaks[np.argsort(peak_powers)[::-1]]
            for p_idx in sorted_peaks[:5]:   # check top-5 spectral peaks
                f_rad = freqs[p_idx]
                period_yr = round((2*np.pi / f_rad) / 365.25, 2)
                if any(abs(period_yr - existing) < 0.05 for existing in model.get("periodic", [])):
                    continue
                period_days = period_yr * 365.25
                # Alternative hypothesis design: two columns (cos + sin); use cos for w-test
                c_cos = np.cos(2*np.pi / period_days * days)
                c_sin = np.sin(2*np.pi / period_days * days)
                w_cos = abs(_w_test_for_candidate(residuals, G, sigma_m, c_cos))
                w_sin = abs(_w_test_for_candidate(residuals, G, sigma_m, c_sin))
                w_period = max(w_cos, w_sin)
                if w_period > best_w:
                    best_w = w_period
                    best_type = "period"
                    best_val = period_yr

    # --- Hypothesis group 3: piecewise-linear velocity break ----------------
    if n >= 20:
        # Velocity spike approach to locate candidate break date
        valid_vel_mask = np.ones(n - 1, dtype=bool)
        dt_arr = np.diff(days)
        if len(dt_arr) > 10:
            valid_vel_mask = dt_arr < np.percentile(dt_arr, 95)
        velocity = np.diff(residuals) / np.maximum(dt_arr, 1e-6)
        valid_vel = velocity[valid_vel_mask]
        if len(valid_vel) > 5:
            median_vel = np.median(valid_vel)
            mad_vel = np.median(np.abs(valid_vel - median_vel))
            if mad_vel > 1e-12:
                min_idx = int(np.argmin(velocity))
                max_idx = int(np.argmax(velocity))
                for cand_idx in [min_idx, max_idx]:
                    break_date = dates[cand_idx].strftime("%Y%m%d")
                    if break_date in model.get("polyline", []):
                        continue
                    t_start = (dates[cand_idx] - dates[0]).days
                    c_pline = np.maximum(days - t_start, 0.0)
                    w_pline = abs(_w_test_for_candidate(residuals, G, sigma_m, c_pline))
                    if w_pline > best_w:
                        best_w = w_pline
                        best_type = "polyline"
                        best_val = break_date

        # CUSUM candidate (Ploberger-Krämer 1992) as additional break hypothesis
        cusum_break = _cusum_break_date(residuals, dates)
        if cusum_break is not None and cusum_break not in model.get("polyline", []):
            try:
                cusum_dt = datetime.strptime(cusum_break, "%Y%m%d")
                t_start_c = (cusum_dt - dates[0]).days
                c_cusum = np.maximum(days - t_start_c, 0.0)
                w_cusum = abs(_w_test_for_candidate(residuals, G, sigma_m, c_cusum))
                if w_cusum > best_w:
                    best_w = w_cusum
                    best_type = "polyline"
                    best_val = cusum_break
            except (ValueError, TypeError):
                pass

    # --- Hypothesis group 4: exponential trend ------------------------------
    if n >= 20 and model.get("exp_trend") is None:
        residual_series = pd.Series(residuals, index=pd.DatetimeIndex(dates))
        best_b = auto_detect_exp_trend(residual_series)
        if best_b is not None:
            c_exp = np.exp(-best_b * days) - 1.0
            w_exp = abs(_w_test_for_candidate(residuals, G, sigma_m, c_exp))
            if w_exp > best_w:
                best_w = w_exp
                best_type = "exp_trend"
                best_val = best_b

    # Significance gate: winning hypothesis must exceed z_{α/2}
    if best_w < z_crit:
        return "None", None, best_w

    return best_type, best_val, best_w


# ---------------------------------------------------------------------------
# Legacy heuristic wrapper (kept for reference; superseded by _identify_best_alternative)
# ---------------------------------------------------------------------------

def robust_analyze_residuals(residuals: np.ndarray, dates: list, model: dict):
    """Heuristic residual analyser (pre-DIA exploration only).

    This function is retained for backward compatibility.
    The formal DIA loop now calls _identify_best_alternative() instead,
    which uses proper w-tests per Teunissen.
    """
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

    if len(residuals) >= 20:
        cusum_break = _cusum_break_date(residuals, dates)
        if cusum_break is not None and cusum_break not in model.get("polyline", []):
            return "polyline", cusum_break

    if len(residuals) >= 20 and model.get("exp_trend") is None:
        residual_series = pd.Series(residuals, index=pd.DatetimeIndex(dates))
        best_b = auto_detect_exp_trend(residual_series)
        if best_b is not None:
            return "exp_trend", best_b

    return "None", None


# ---------------------------------------------------------------------------
# Main DIA loop
# ---------------------------------------------------------------------------

def run_omt_dia_loop(series, jump_dates, initial_periods, initial_polylines,
                     initial_logs, poly_deg, sigma_mm, alpha, max_iter,
                     exp_trend_b=None):
    """
    Run the Teunissen DIA loop for a single assumed a-priori sigma.

    Parameters
    ----------
    series           : pd.Series  time-indexed displacement (metres)
    jump_dates       : list of datetime  known/detected step offsets
    initial_periods  : list of float     candidate periodic periods (years)
    initial_polylines: list of str       pre-specified velocity-break dates (YYYYMMDD)
    initial_logs     : dict              pre-specified logarithmic relaxations
    poly_deg         : int               polynomial degree for functional model
    sigma_mm         : float             a-priori measurement noise assumption [mm]
                                         NOTE: sigma_mm is an assumed a-priori value.
                                         The sigma scan selects the sigma that satisfies
                                         the OMT — the selected value is therefore an
                                         a-posteriori estimate, not a pre-set instrument spec.
    alpha            : float             significance level (e.g. 0.05)
    max_iter         : int               maximum DIA iterations
    exp_trend_b      : float or None     pre-specified exponential trend decay [1/day]

    Returns
    -------
    Accepted model dict (with '_omt_stats' key) or None if no model accepted.

    Workflow per iteration
    ----------------------
    1. Estimate parameters via least squares.
    2. DETECT: compute T_stat = SSR/σ²; accept if p_value ≥ alpha.
    3. IDENTIFY: call _identify_best_alternative() to find max|w| candidate.
    4. ADAPT: extend functional model with identified component.
    5. Repeat until accepted or no more candidates.
    """
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

    # Pre-compute loop-invariant arrays for _identify_best_alternative
    _days_arr = np.array([(d - date_list[0]).days for d in date_list], dtype=float)
    _n = len(dis_ts)
    # Frequency grid: 500 points is sufficient to resolve peaks; 5000 was wasteful for long series
    _freqs_arr = np.linspace(2*np.pi/(20*365.25), 2*np.pi/(0.2*365.25), 500) if _n >= 20 else None

    last_unit_var = 9999.0
    for iteration in range(max_iter):

        # --- ESTIMATION ---
        G, m, e2, d_hat = estimate_time_func(model, date_list, dis_ts)
        residuals = dis_ts - d_hat
        n_param = G.shape[1]

        # --- DETECTION ---
        T_stat, unit_var_factor, p_value, chi2_critical, sigma_hat_m = calculate_omt(
            residuals, len(dis_ts), n_param, sigma_m, alpha
        )

        if p_value >= alpha:
            # Model accepted: store full OMT statistics
            model["_omt_stats"] = {
                "sigma_mm":        sigma_mm,
                "sigma_hat_mm":    sigma_hat_m * 1000.0,   # a-posteriori noise estimate [mm]
                "T_stat":          T_stat,
                "chi2_critical":   chi2_critical,
                "unit_var_factor": unit_var_factor,         # σ̂²/σ₀²; ~1 means good fit
                "p_value":         p_value,
                "n_param":         n_param,
                "r":               len(dis_ts) - n_param,  # degrees of freedom
                "iterations":      iteration,
            }
            return model

        # Convergence guard: stop if unit variance factor is no longer improving
        if unit_var_factor >= last_unit_var and iteration > 0:
            break
        last_unit_var = unit_var_factor

        # --- IDENTIFICATION (formal w-test per Teunissen) ---
        adapt_type, adapt_val, max_w = _identify_best_alternative(
            residuals, date_list, model, G, sigma_m, alpha,
            days=_days_arr, freqs=_freqs_arr
        )

        # --- ADAPTATION ---
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

        if adapt_type == "outlier":
            # Outlier identified: flag by inserting a unit-step at that epoch
            # (equivalent to removing the epoch from estimation)
            outlier_date = date_list[adapt_val].strftime("%Y%m%d")
            if outlier_date not in model["stepDate"]:
                model["stepDate"].append(outlier_date)
            else:
                break

        if adapt_type == "polyline":
            if adapt_val is None:
                # Fallback: use maximum velocity as break date
                days_arr = np.array([(d - date_list[0]).days for d in date_list])
                velocity = np.diff(residuals) / np.maximum(np.diff(days_arr), 1e-6)
                ext_idx = int(np.argmax(np.abs(velocity)))
                adapt_val = date_list[ext_idx].strftime("%Y%m%d")
            if adapt_val not in model["polyline"]:
                model["polyline"].append(adapt_val)
            else:
                break

        if adapt_type == "None":
            break

    return None


# ---------------------------------------------------------------------------
# Post-seismic relaxation test (called after DIA loop accepts a model)
# ---------------------------------------------------------------------------

def test_relaxation(series, jump_dates, accepted_model, alpha):
    """Test whether adding exponential relaxation after each jump improves the model.

    Selection criterion: choose the tau that gives the highest p_value while
    still passing the OMT (p_value >= alpha).  This replaces the previous
    unit_var_factor comparison which was not valid across models with different
    degrees of freedom.
    """
    if not jump_dates:
        return accepted_model
    series_clean = series.dropna()
    if series_clean.empty:
        return accepted_model

    date_list = [d.to_pydatetime() if hasattr(d, "to_pydatetime") else d for d in series_clean.index]
    dis_ts = series_clean.values.copy()
    sigma_mm = accepted_model["_omt_stats"]["sigma_mm"]
    sigma_m = sigma_mm / 1000.0

    model = {k: v for k, v in accepted_model.items() if not k.startswith("_")}
    model["exp"] = {}

    current_stats = accepted_model["_omt_stats"]
    best_p_value = current_stats["p_value"]   # use p_value for model comparison (fix #9)

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
            _, _, p_value, _, _ = calculate_omt(residuals, len(dis_ts), n_param, sigma_m, alpha)

            if p_value >= alpha and p_value > best_p_value:
                best_p_value = p_value
                best_tau = tau_days

        if best_tau is not None:
            model["exp"][jump_str] = [best_tau]

    G, m, e2, d_hat = estimate_time_func(model, date_list, dis_ts)
    residuals = dis_ts - d_hat
    n_param = G.shape[1]
    T_stat, unit_var_factor, p_value, chi2_critical, sigma_hat_m = calculate_omt(
        residuals, len(dis_ts), n_param, sigma_m, alpha
    )
    model["_omt_stats"] = {
        "sigma_mm":        sigma_mm,
        "sigma_hat_mm":    sigma_hat_m * 1000.0,
        "T_stat":          T_stat,
        "chi2_critical":   chi2_critical,
        "unit_var_factor": unit_var_factor,
        "p_value":         p_value,
        "n_param":         n_param,
        "r":               len(dis_ts) - n_param,
        "iterations":      current_stats["iterations"],
    }
    return model


# ---------------------------------------------------------------------------
# Sigma scan workers
# ---------------------------------------------------------------------------

def _run_single_sigma(sigma_mm, series, jump_dates, candidate_periods,
                      initial_polylines, initial_logs, poly_deg, alpha, max_iter,
                      exp_trend_b=None):
    result = run_omt_dia_loop(
        series, jump_dates, candidate_periods, initial_polylines, initial_logs,
        poly_deg, sigma_mm, alpha, max_iter, exp_trend_b=exp_trend_b
    )
    if result is not None:
        s = result["_omt_stats"]
        return result, {
            "sigma_mm":   sigma_mm,
            "accepted":   True,
            "p_value":    s["p_value"],
            "T_stat":     s["T_stat"],
            "n_param":    s["n_param"],
            "n_periods":  len(result.get("periodic", [])),
            "n_polylines": len(result.get("polyline", [])),
        }
    return None, {
        "sigma_mm":   sigma_mm,
        "accepted":   False,
        "p_value":    None,
        "T_stat":     None,
        "n_param":    None,
        "n_periods":  None,
        "n_polylines": None,
    }


def run_omt_sigma_scan(series, jump_dates, candidate_periods, initial_polylines,
                       initial_logs, poly_deg, sigma_min, sigma_max, sigma_step,
                       alpha, max_iter, no_relax=False, cores=1, exp_trend_b=None):
    """
    Scan a range of assumed a-priori sigma values and run the full DIA loop
    for each.  The selected sigma is the one that satisfies the OMT with the
    most parsimonious model — it is therefore an a-posteriori estimate chosen
    to honour the test, not a pre-set instrument specification.
    """
    sigmas = np.arange(sigma_min, sigma_max + sigma_step * 0.5, sigma_step)
    scan_results = []
    scan_table = []

    if cores > 1:
        import concurrent.futures
        from functools import partial
        worker = partial(
            _run_single_sigma,
            series=series, jump_dates=jump_dates,
            candidate_periods=candidate_periods, initial_polylines=initial_polylines,
            initial_logs=initial_logs, poly_deg=poly_deg, alpha=alpha,
            max_iter=max_iter, exp_trend_b=exp_trend_b
        )
        with concurrent.futures.ProcessPoolExecutor(max_workers=cores) as executor:
            for result, table_row in executor.map(worker, sigmas):
                scan_table.append(table_row)
                if result is not None:
                    scan_results.append(result)
    else:
        for sigma_mm in sigmas:
            result, table_row = _run_single_sigma(
                sigma_mm, series, jump_dates, candidate_periods, initial_polylines,
                initial_logs, poly_deg, alpha, max_iter, exp_trend_b=exp_trend_b
            )
            scan_table.append(table_row)
            if result is not None:
                scan_results.append(result)

    if not scan_results:
        return None, scan_table

    # Select most parsimonious model that passed OMT
    best = min(scan_results, key=lambda m: (
        m["_omt_stats"]["n_param"],
        m["_omt_stats"]["sigma_mm"],
        -m["_omt_stats"]["p_value"],
    ))
    print(f"  [sigma-scan] Selected sigma={best['_omt_stats']['sigma_mm']:.1f} mm "
          f"(a-posteriori; sigma_hat={best['_omt_stats']['sigma_hat_mm']:.2f} mm, "
          f"p={best['_omt_stats']['p_value']:.4f})")

    if not no_relax:
        best = test_relaxation(series, jump_dates, best, alpha)
    return best, scan_table
