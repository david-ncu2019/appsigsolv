"""Core parametric modeling logic (design matrices and estimation)."""
import math
import numpy as np
from scipy import linalg
from datetime import datetime
import pandas as pd

def datetime2years(date_list):
    """
    Convert a list of datetime objects or strings to decimal years.
    """
    years = []
    for d in date_list:
        if isinstance(d, datetime):
            dt = d
        else:
            try:
                dt = datetime.strptime(str(d), "%Y%m%d")
            except ValueError:
                dt = datetime.strptime(str(d), "%Y-%m-%d")
        
        year = dt.year
        start_of_year = datetime(year, 1, 1)
        start_of_next_year = datetime(year + 1, 1, 1)
        year_elapsed = (dt - start_of_year).total_seconds()
        year_duration = (start_of_next_year - start_of_year).total_seconds()
        decimal_year = year + year_elapsed / year_duration
        years.append(decimal_year)
    return years

def get_design_matrix4polynomial_func(yr_diff, degree):
    A = np.zeros([len(yr_diff), degree + 1], dtype=np.float64)
    for i in range(degree+1):
        A[:, i] = (yr_diff**i) / math.factorial(i)
    return A

def get_design_matrix4periodic_func(yr_diff, periods):
    A = np.zeros((len(yr_diff), 2*len(periods)), dtype=np.float64)
    for i, period in enumerate(periods):
        c0, c1 = 2*i, 2*i+1
        A[:, c0] = np.cos(2*np.pi/period * yr_diff)
        A[:, c1] = np.sin(2*np.pi/period * yr_diff)
    return A

def get_design_matrix4step_func(date_list, step_date_list):
    t = np.array(datetime2years(date_list))
    t_steps = datetime2years(step_date_list)
    A = np.zeros((len(t), len(t_steps)), dtype=np.float64)
    for i, t_step in enumerate(t_steps):
        A[:, i] = np.array(t > t_step).flatten()
    return A

def get_design_matrix4polyline(date_list, start_date_list):
    t = np.array(datetime2years(date_list))
    t_start_list = datetime2years(start_date_list)
    A = np.zeros((len(t), len(t_start_list)), dtype=np.float64)
    for i, t_start in enumerate(t_start_list):
        tbase = t - t_start
        tbase[tbase < 0] = 0
        A[:, i] = tbase
    return A

def get_design_matrix4exp_func(date_list, exp_dict):
    t = np.array(datetime2years(date_list))
    num_exp = sum(len(val) for key, val in exp_dict.items())
    A = np.zeros((len(t), num_exp), dtype=np.float64)
    i = 0
    for exp_onset, taus in exp_dict.items():
        exp_T = datetime2years([exp_onset])[0]
        for exp_tau in taus:
            exp_tau_yr = exp_tau / 365.25
            A[:, i] = np.array(t > exp_T).flatten() * (1 - np.exp(-1 * (t - exp_T) / exp_tau_yr))
            i += 1
    return A

def get_design_matrix4log_func(date_list, log_dict):
    t = np.array(datetime2years(date_list))
    num_log = sum(len(val) for key, val in log_dict.items())
    A = np.zeros((len(t), num_log), dtype=np.float64)
    i = 0
    for log_onset, taus in log_dict.items():
        log_T = datetime2years([log_onset])[0]
        for log_tau in taus:
            log_tau_yr = log_tau / 365.25
            olderr = np.seterr(invalid='ignore', divide='ignore')
            A[:, i] = np.array(t > log_T).flatten() * np.nan_to_num(
                np.log(1 + (t - log_T) / log_tau_yr),
                nan=0, neginf=0
            )
            np.seterr(**olderr)
            i += 1
    return A

def get_design_matrix4exp_trend(date_list, b_per_day: float):
    """One-column basis: exp(-b * days_from_start) - 1. At t=0: value=0."""
    t0 = date_list[0]
    days = np.array([(d - t0).days for d in date_list], dtype=np.float64)
    return (np.exp(-b_per_day * days) - 1.0).reshape(-1, 1)

def get_design_matrix4time_func(date_list, model, ref_date=None):
    """
    Design matrix for time functions parameter estimation.
    """
    yr_diff = np.array(datetime2years(date_list))

    if ref_date is None:
        ref_date = date_list[0]
    
    ref_idx = date_list.index(ref_date) if ref_date in date_list else 0
    yr_diff -= yr_diff[ref_idx]

    poly_deg      = model.get('polynomial', 0)
    periods       = model.get('periodic', [])
    steps         = model.get('stepDate', [])
    polylines     = model.get('polyline', [])
    exps          = model.get('exp', dict())
    logs          = model.get('log', dict())
    exp_trend_b   = model.get('exp_trend', None)

    num_period    = len(periods)
    num_step      = len(steps)
    num_pline     = len(polylines)
    num_exp       = sum(len(val) for key, val in exps.items())
    num_log       = sum(len(val) for key, val in logs.items())
    num_exp_trend = 1 if (exp_trend_b is not None and exp_trend_b != 0) else 0

    num_param = (poly_deg + 1) + num_exp_trend + (2 * num_period) + num_step + num_pline + num_exp + num_log
    if num_param <= 0:
        raise ValueError('NO time functions specified!')

    num_date = len(yr_diff)
    A = np.zeros((num_date, num_param), dtype=np.float64)
    c0 = 0

    c1 = c0 + poly_deg + 1
    A[:, c0:c1] = get_design_matrix4polynomial_func(yr_diff, poly_deg)
    c0 = c1

    if num_exp_trend > 0:
        c1 = c0 + 1
        A[:, c0:c1] = get_design_matrix4exp_trend(date_list, exp_trend_b)
        c0 = c1

    if num_period > 0:
        c1 = c0 + 2 * num_period
        A[:, c0:c1] = get_design_matrix4periodic_func(yr_diff, periods)
        c0 = c1

    if num_step > 0:
        c1 = c0 + num_step
        A[:, c0:c1] = get_design_matrix4step_func(date_list, steps)
        c0 = c1

    if num_pline > 0:
        c1 = c0 + num_pline
        A[:, c0:c1] = get_design_matrix4polyline(date_list, polylines)
        c0 = c1

    if num_exp > 0:
        c1 = c0 + num_exp
        A[:, c0:c1] = get_design_matrix4exp_func(date_list, exps)
        c0 = c1

    if num_log > 0:
        c1 = c0 + num_log
        A[:, c0:c1] = get_design_matrix4log_func(date_list, logs)
        c0 = c1

    return A

def estimate_time_func(model, date_list, dis_ts):
    """
    Fit the deformation model to a 1D timeseries.
    Returns: G, m, e2, d_hat
    """
    dis_ts_2d = dis_ts.reshape(-1, 1)
    G = get_design_matrix4time_func(date_list, model)
    m, e2 = linalg.lstsq(G, dis_ts_2d, cond=None)[:2]
    
    if len(e2) == 0:
        d_hat = G @ m
        residuals = dis_ts_2d - d_hat
        e2 = np.array([np.sum(residuals**2)])
    else:
        d_hat = G @ m

    return G, m.flatten(), e2[0], d_hat.flatten()

def _period_label(period_yr: float) -> str:
    val = period_yr
    if val == int(val):
        return f"{int(val)}yr"
    return f"{val}yr"

def extract_components(series: pd.Series, best_model: dict, comp: str) -> dict:
    """
    Extract individual signal components from the fitted model.
    """
    series_clean = series.dropna()
    if series_clean.empty:
        return {}
    date_list = [d.to_pydatetime() if hasattr(d, "to_pydatetime") else d for d in series_clean.index]
    dis_ts = series_clean.values.copy()

    model = {k: v for k, v in best_model.items() if not k.startswith("_")}
    G, m, e2, d_hat = estimate_time_func(model, date_list, dis_ts)

    poly_deg      = model.get("polynomial", 0)
    periods       = model.get("periodic", [])
    steps         = model.get("stepDate", [])
    polylines     = model.get("polyline", [])
    exps          = model.get("exp", {})
    logs          = model.get("log", {})
    exp_trend_b   = model.get("exp_trend", None)
    n_exp_trend   = 1 if (exp_trend_b is not None and exp_trend_b != 0) else 0

    components = {}
    full_date_list = [d.to_pydatetime() if hasattr(d, "to_pydatetime") else d for d in series.index]
    G_full = get_design_matrix4time_func(full_date_list, model)
    idx = series.index

    # Polynomial trend
    c0 = 0
    c1 = poly_deg + 1
    trend_val = G_full[:, c0:c1] @ m[c0:c1]
    c0 = c1

    # Exponential trend component (new primary trend shape)
    if n_exp_trend > 0:
        c1 = c0 + 1
        exp_trend_val = G_full[:, c0:c1] @ m[c0:c1]
        components[f"{comp}_exp_trend"] = pd.Series(exp_trend_val, index=idx)
        c0 = c1

    # Polylines are folded into the trend display
    curr_c0 = poly_deg + 1 + n_exp_trend + 2 * len(periods) + len(steps)
    c1_pline = curr_c0 + len(polylines)
    if len(polylines) > 0:
        trend_val = trend_val + G_full[:, curr_c0:c1_pline] @ m[curr_c0:c1_pline]

    components[f"{comp}_trend"] = pd.Series(trend_val, index=idx)

    # Periodic
    c0 = poly_deg + 1 + n_exp_trend
    for period_yr in periods:
        label = _period_label(period_yr)
        c1 = c0 + 2
        seasonal_val = G_full[:, c0:c1] @ m[c0:c1]
        components[f"{comp}_{label}"] = pd.Series(seasonal_val, index=idx)
        c0 = c1

    # Steps / jumps
    c0 = poly_deg + 1 + n_exp_trend + 2 * len(periods)
    if steps:
        c1 = c0 + len(steps)
        jump_val = G_full[:, c0:c1] @ m[c0:c1]
        components[f"{comp}_jump"] = pd.Series(jump_val, index=idx)
        c0 = c1

    # Post-seismic exponential relaxation
    c0 = poly_deg + 1 + n_exp_trend + 2 * len(periods) + len(steps) + len(polylines)
    n_exp = sum(len(v) for v in exps.values())
    if n_exp > 0:
        c1 = c0 + n_exp
        exp_val = G_full[:, c0:c1] @ m[c0:c1]
        components[f"{comp}_exp"] = pd.Series(exp_val, index=idx)
        c0 = c1

    # Logarithmic relaxation
    n_log = sum(len(v) for v in logs.values())
    if n_log > 0:
        c1 = c0 + n_log
        log_val = G_full[:, c0:c1] @ m[c0:c1]
        components[f"{comp}_log"] = pd.Series(log_val, index=idx)
        c0 = c1

    components[f"{comp}_model"] = pd.Series(G_full @ m, index=idx)
    noise = series.values - components[f"{comp}_model"].values
    components[f"{comp}_noise"] = pd.Series(noise, index=idx)

    sigma_mm = best_model.get("_omt_stats", {}).get("sigma_mm", 1.0)
    sigma_m = sigma_mm / 1000.0
    w_stat = noise / sigma_m
    components[f"{comp}_wtest"] = pd.Series(w_stat, index=idx)

    return components
