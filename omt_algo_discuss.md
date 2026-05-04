# Algorithm vs. Theory: Overall Model Test (OMT) and DIA Procedure

This document discusses whether the algorithms implemented in `D:\1000_SCRIPTS\004_Project003\20260501_timeseries_signal_solver\appsigsolv` follow the theoretical framework of the Overall Model Test (OMT) and DIA (Detection, Identification, Adaptation) procedure, as defined in the `Overall_Model_Test` NotebookLM notebook.

## 1. Detection Phase: Least-Squares and the OMT Statistic

### Theory (from NotebookLM)
The NotebookLM theory specifies:
- The null hypothesis $H_0: E\{y\} = Ax$, $D\{y\} = Q_{yy}$.
- The least-squares residual vector $\hat{e}_0 = y - A\hat{x}_0$.
- The OMT test statistic is $T = \hat{e}_0^T Q_{yy}^{-1} \hat{e}_0$. Assuming a uniform diagonal covariance matrix $Q_{yy} = \sigma_m^2 I$, this simplifies to $T = \frac{\sum \hat{e}_i^2}{\sigma_m^2}$.
- Redundancy (degrees of freedom) is $r = m - n$, where $m$ is the number of observations and $n$ is the number of parameters.
- $T$ follows a Chi-square distribution $\chi^2(m-n, 0)$. The null hypothesis is rejected if $T > K_\alpha$, or equivalently if the $p$-value $< \alpha$.

### Code Implementation
**Yes, the code perfectly follows this theory.**

In `appsigsolv/core/modeling.py`, the least-squares estimation fits the design matrix $A$ (represented as `G`) to the observations $y$ (`dis_ts_2d`):
```python
# appsigsolv/core/modeling.py
def estimate_time_func(model, date_list, dis_ts):
    dis_ts_2d = dis_ts.reshape(-1, 1)
    G = get_design_matrix4time_func(date_list, model)
    m, e2 = linalg.lstsq(G, dis_ts_2d, cond=None)[:2]
    
    d_hat = G @ m
    residuals = dis_ts_2d - d_hat
    return G, m.flatten(), np.sum(residuals**2), d_hat.flatten()
```

In `appsigsolv/core/dia.py`, the OMT test statistic $T$ and degrees of freedom $r$ are calculated exactly as described in the theory:
```python
# appsigsolv/core/dia.py
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
```

## 2. Identification Phase: The w-test

### Theory (from NotebookLM)
The notebook states that to pinpoint the error source, Baarda's w-test statistic is used:
$$w_i = \frac{c_i^T Q_{yy}^{-1} \hat{e}_0}{\sqrt{c_i^T Q_{yy}^{-1} Q_{\hat{e}_0} Q_{yy}^{-1} c_i}}$$
For uncorrelated observations with variance $\sigma_m^2$, the exact denominator should be $\sigma_m \sqrt{(Q_{\hat{e}_0})_{ii}}$, where $Q_{\hat{e}_0} = Q_{yy} - A(A^T Q_{yy}^{-1} A)^{-1} A^T$.

### Code Implementation
**The code uses a simplified approximation of the theory.**

In `appsigsolv/core/modeling.py`, the w-test statistic is computed as:
```python
# appsigsolv/core/modeling.py
def extract_components(series: pd.Series, best_model: dict, comp: str) -> dict:
    ...
    noise = series.values - components[f"{comp}_model"].values
    sigma_mm = best_model.get("_omt_stats", {}).get("sigma_mm", 1.0)
    sigma_m = sigma_mm / 1000.0
    
    w_stat = noise / sigma_m
    components[f"{comp}_wtest"] = pd.Series(w_stat, index=idx)
```

**Discussion:** The implementation calculates $w_i = \frac{\hat{e}_i}{\sigma_m}$. It omits the term representing the redundancy matrix from the denominator (i.e., it assumes $Q_{\hat{e}_0} \approx Q_{yy}$). 
In highly redundant timeseries datasets where $m \gg n$ (e.g., thousands of daily observations vs. 10-15 model parameters), the diagonal elements of the projection matrix $A(A^T A)^{-1} A^T$ approach 0. Therefore, $\sqrt{1 - \text{diag}(...)} \approx 1$. Thus, dividing by just $\sigma_m$ is an acceptable and highly efficient computational approximation in this specific domain, though it is not the *strict* theoretical equation provided in the notebook.

## 3. Adaptation Phase

### Theory (from NotebookLM)
The notebook states: *"Once an error is identified, the model is modified to eliminate its influence. This may involve 'data snooping' (removing the erroneous observation) or extending the functional model to include additional parameters that represent the error."*

### Code Implementation
**Yes, the code follows the functional model extension approach.**

In `appsigsolv/core/dia.py`, when the OMT test fails (the $p$-value $< \alpha$), the algorithm does not just blindly remove points. Instead, it enters the `robust_analyze_residuals` function to identify what is missing in the model (e.g., an unmodeled seasonal signal or a sudden structural jump/polyline break) and adapts the model by extending it:
```python
# appsigsolv/core/dia.py (inside run_omt_dia_loop)
    for iteration in range(max_iter):
        ...
        if p_value >= alpha:
            return model # OMT passed
            
        # OMT rejected -> Identify & Adapt
        adapt_type, adapt_val = robust_analyze_residuals(residuals, date_list, model)

        if adapt_type == "period":
            model["periodic"].append(adapt_val)
        if adapt_type == "polyline":
            model["polyline"].append(adapt_val)
```

## Conclusion

The `appsigsolv` package strongly adheres to the theoretical DIA procedure defined in the `Overall_Model_Test` notebook:
1. **Detection:** Implements the exact $\chi^2$-based Overall Model Test and calculates degrees of freedom flawlessly.
2. **Identification:** Computes the `w-stat` but uses a computationally efficient large-redundancy approximation rather than the strict covariance propagation formula.
3. **Adaptation:** Autonomously extends the functional model (adding periodic parameters or structural breaks) to accommodate unmodeled signals rather than simply discarding data, reflecting the more advanced adaptation method mentioned in the literature.