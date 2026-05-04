# The Overall Model Test: A Complete Lecture on Geodetic Quality Control

**Running example used throughout:** A GNSS station recording daily vertical
positions over three years to monitor ground subsidence. This is one of the
simplest and most common geodetic timeseries — and every concept introduced
below applies directly to it.

---

## Part 1 — Why Models Fail and Why It Matters

Imagine a GNSS monument anchored into a concrete pier in a sedimentary basin.
Every day, the receiver records its vertical position. Over three years, you
have roughly $m = 1095$ daily measurements. You plot them and the data trends
downward — the ground is sinking. You fit a straight line and call that your
estimate of the subsidence velocity.

This is the classic estimation problem. But before you trust that velocity, you
must answer one uncomfortable question: **does your model actually represent
what happened in those three years?**

Perhaps on day 312, a technician bumped the antenna and the position jumped by
20 mm. Perhaps the soil undergoes seasonal shrink–swell cycles that a straight
line cannot capture. Perhaps an earthquake on day 740 permanently offset the
monument by 8 mm. Any of these effects will contaminate your velocity estimate
without announcing themselves — your least-squares algorithm will "absorb" the
anomaly by distorting all the other parameters.

This is the problem the **DIA procedure** solves.

### 1.1  The Two-Part Mathematical Model

Every least-squares estimation requires two model components defined
simultaneously.

**The Functional Model** defines what you *expect* to observe:

$$E\{y\} = Ax$$

Here, $y$ is the $m \times 1$ vector of daily vertical positions (observations),
$x$ is the $n \times 1$ vector of unknown parameters, and $A$ is the
$m \times n$ design matrix encoding the physical relationship between
observations and unknowns.

For a linear subsidence model, $x = [x_0,\; v]^T$ where $x_0$ is the initial
height and $v$ is the velocity. The design matrix has one row per day:

$$A = \begin{bmatrix} 1 & t_1 \\ 1 & t_2 \\ \vdots & \vdots \\ 1 & t_m \end{bmatrix}$$

**The Stochastic Model** defines the *expected noise*:

$$D\{y\} = Q_{yy}$$

The covariance matrix $Q_{yy}$ is the mathematical container for every random
uncertainty in the data. In a GNSS timeseries, $Q_{yy}$ represents receiver
thermal noise, residual atmospheric delays from the troposphere and ionosphere,
and multipath scattering from nearby reflective surfaces.

The simplest assumption is that all observations are independent and share the
same noise variance $\sigma_m^2$, giving:

$$Q_{yy} = \sigma_m^2 I$$

This assumption treats measurement errors as **white noise** — no temporal
memory, no correlation between days. Real GNSS data violates this regularly:
high-rate sampling creates time-correlated noise, regional atmospheric
conditions create spatial correlation between nearby stations, and satellite
signals at low elevation angles carry more error than signals near the zenith.
These violations matter because an incorrectly specified $Q_{yy}$ will produce
parameter estimates with falsely optimistic precision, and will cause the OMT
to give misleading verdicts.

### 1.2  The Delft School Philosophy: Estimation and Testing Are Inseparable

The classical geodetic workflow treated data adjustment (estimation) and quality
checking (testing) as two disconnected steps: first fit the model, then inspect
the residuals. The **Delft School of Geodesy**, through the work of Willem
Baarda and later Peter Teunissen, overturned this view.

Their central argument: **you cannot reliably estimate parameters without first
validating the model that produced them.** If the functional or stochastic
model is wrong, the least-squares estimator dutifully produces a "best fit" to
the wrong model — a precise answer to the wrong question.

Teunissen formalized this insight into the **DIA framework** (Detection,
Identification, Adaptation), proving mathematically that the final parameter
estimate is not just a function of the observations but of the *entire
decision chain* — including every statistical test applied along the way. Even
when each individual estimator is unbiased under its specific hypothesis, the
combined DIA-estimator retains a **conditional bias** because it inherits the
uncertainty of each testing decision. This probabilistic view fundamentally
changed how geodesists quantify integrity risk in safety-critical applications
such as aviation GNSS navigation.

The practical consequence: quality control is not a post-processing check. It
is an integral part of the estimation process itself.

---

## Part 2 — Detection: The Overall Model Test (OMT)

### 2.1  What the OMT Does (and Does Not Do)

The OMT is the first and most important test in the DIA cycle. It asks a single
global question:

> **Are the residuals of my fitted model consistent with the noise level I
> claimed in $Q_{yy}$?**

It operates under a **null hypothesis** $H_0$:

> Both the functional model and the stochastic model are correct and fully
> represent the data.

The OMT is an **unspecified test** — it can detect that something is wrong
without identifying what. A rejection tells you the model failed globally. It
does not tell you whether the cause was a bad data point on day 312, an
unmodeled seasonal cycle, or an underestimated noise variance. That
specificity is the job of the Identification phase.

This is a deliberate design choice. The OMT pays the price of non-specificity
in exchange for global validity: it will catch *any* form of misspecification,
regardless of its origin.

### 2.2  Building the Test Statistic $T$

After fitting the linear model by least squares, we obtain the **residual
vector**:

$$\hat{e}_0 = y - A\hat{x}_0$$

Each element $\hat{e}_{0,i}$ is the vertical distance between the actual
measurement on day $i$ and the fitted straight line. If the model is correct,
these residuals should look like random noise drawn from the distribution
implied by $Q_{yy}$.

We need a single scalar to summarise the entire $m$-dimensional residual
vector. Simply summing the squared residuals, $\hat{e}_0^T \hat{e}_0$, ignores
the fact that not all measurements are equally reliable. We must weight each
residual by its inverse variance. This produces the **test statistic**:

$$\boxed{T = \hat{e}_0^T Q_{yy}^{-1} \hat{e}_0}$$

#### Derivation for the uncorrelated equal-variance case

Assume $Q_{yy} = \sigma_m^2 I$. Substituting into $T$:

$$T = \hat{e}_0^T (\sigma_m^2 I)^{-1} \hat{e}_0
    = \hat{e}_0^T \left(\frac{1}{\sigma_m^2} I\right) \hat{e}_0
    = \frac{1}{\sigma_m^2}\, \hat{e}_0^T \hat{e}_0$$

Because $\hat{e}_0^T \hat{e}_0 = \sum_{i=1}^m \hat{e}_{0,i}^2$, the test
statistic simplifies to:

$$T = \frac{\displaystyle\sum_{i=1}^m \hat{e}_{0,i}^2}{\sigma_m^2}$$

In words: $T$ is the sum of the squared daily position residuals, normalised by
the expected daily noise variance. If the noise model is correct, this ratio
should be close to the number of degrees of freedom.

### 2.3  The Statistical Distribution of $T$

Understanding *why* $T$ follows a particular distribution requires one
foundational fact: if the observations $y$ are normally distributed under $H_0$,
then the residuals $\hat{e}_0$ are also normally distributed. The weighted sum
of squares of normally distributed variables follows a **Chi-square
distribution**.

**Under $H_0$ (model is correct):**

$$T \sim \chi^2(r)$$

where the degrees of freedom $r = m - n$ is the **redundancy** — the number of
observations minus the number of estimated parameters. For our 3-year GNSS
timeseries with $m = 1095$ and $n = 2$, we have $r = 1093$.

The redundancy has a physical interpretation: it counts how many measurements
are "left over" after satisfying the model's information demands. High
redundancy means the data contains many independent checks — the system is
geometrically strong and statistically sensitive to anomalies.

**Under $H_a$ (model is wrong):**

If a bias exists — whether a data error, an unmodeled physical effect, or a
misspecified stochastic model — the test statistic shifts. Instead of the
central Chi-square distribution, $T$ follows a **non-central Chi-square**:

$$T \sim \chi^2(r, \lambda)$$

The non-centrality parameter $\lambda > 0$ measures the severity of the
violation. As the bias grows larger, $\lambda$ increases, pushing the
distribution's probability mass further into the tail. The entire curve shifts
to the right, away from the decision threshold.

### 2.4  Making the Decision

We select a **significance level** $\alpha$, typically $0.05$. This is the
probability of falsely rejecting a valid $H_0$ (a **false alarm** or Type I
error). From the central $\chi^2(r)$ distribution, we find the **critical
value** $K_\alpha$ such that:

$$P(T > K_\alpha \mid H_0) = \alpha$$

The decision rule is simple:

| Result | Interpretation |
|--------|----------------|
| $T \leq K_\alpha$ | **Accept $H_0$.** The model fits the data. Residuals are consistent with your noise assumptions. The velocity estimate is trustworthy. |
| $T > K_\alpha$ | **Reject $H_0$.** The model fails. Something is wrong — a data error, missing physics, or a misspecified covariance. Proceed to Identification. |

**GNSS example in numbers:**

A 3-year daily GNSS timeseries ($r = 1093$) with $\alpha = 0.05$ gives
$K_\alpha \approx 1148$. If $T = 1250$, the OMT rejects $H_0$. If $T = 1090$,
it accepts. The computed value of $T$ relative to $K_\alpha$ is the only thing
that determines the verdict.

---

## Part 3 — Identification: Finding the Source of Failure

### 3.1  The Role of Identification

A rejected OMT is an alarm, not a diagnosis. The Identification phase runs
targeted tests to locate the specific source of the failure. There are two
complementary strategies:

1. **Baarda's w-test** — tests whether a specific individual observation is a
   blunder (data error).
2. **Multiple Hypotheses Testing** — tests whether a specific structural
   addition to the functional model (e.g., a seasonal term, a step function)
   would explain the rejection.

### 3.2  Baarda's w-Test: Testing a Single Observation

The w-test isolates one daily measurement and asks: *is this observation so
inconsistent with all other data that it must be an outlier?*

For the $i$-th observation, define the **selector vector** $c_i$ — a vector of
zeros with a single 1 in position $i$. It mathematically extracts the $i$-th
element from any vector.

The exact theoretical w-statistic is:

$$w_i = \frac{c_i^T Q_{yy}^{-1} \hat{e}_0}{\sqrt{c_i^T Q_{yy}^{-1} Q_{\hat{e}_0} Q_{yy}^{-1} c_i}}$$

where $Q_{\hat{e}_0}$ is the exact **residual covariance matrix**:

$$Q_{\hat{e}_0} = Q_{yy} - A(A^T Q_{yy}^{-1} A)^{-1} A^T$$

Under $H_0$, $w_i \sim \mathcal{N}(0,1)$. If $|w_i| > N_{\alpha/2}$ (e.g.,
$|w_i| > 1.96$ at $\alpha = 0.05$), the $i$-th observation is flagged as
suspicious.

#### The high-redundancy approximation

Computing $Q_{\hat{e}_0}$ exactly requires forming an $m \times m$ matrix
product — computationally expensive for long timeseries. For our GNSS case,
the redundancy is massive ($r/m = 1093/1095 \approx 0.998$). When $m \gg n$,
the projection matrix $A(A^T Q_{yy}^{-1} A)^{-1} A^T$ becomes negligible
compared to $Q_{yy}$, so:

$$Q_{\hat{e}_0} \approx Q_{yy} = \sigma_m^2 I$$

Substituting $Q_{yy} = \sigma_m^2 I$ and this approximation into the w-formula:

**Numerator:**
$$c_i^T \left(\tfrac{1}{\sigma_m^2} I\right) \hat{e}_0 = \frac{\hat{e}_{0,i}}{\sigma_m^2}$$

**Denominator:**
$$\sqrt{c_i^T \left(\tfrac{1}{\sigma_m^2}I\right)(\sigma_m^2 I)\left(\tfrac{1}{\sigma_m^2}I\right) c_i} = \sqrt{\tfrac{1}{\sigma_m^2}} = \frac{1}{\sigma_m}$$

**Combined:**
$$\boxed{w_i \approx \frac{\hat{e}_{0,i}}{\sigma_m}}$$

In plain language: the w-statistic for any day is simply that day's residual
divided by the expected noise standard deviation. For our GNSS example with
$\sigma_m = 3$ mm, a residual of $\hat{e}_{0,i} = 15$ mm gives $w_i = 5.0$ —
far exceeding the 1.96 threshold. That day is almost certainly an outlier.

#### When the approximation is safe — and when it is not

| Redundancy ratio $r/m$ | Approximation |
|------------------------|---------------|
| $> 0.5$ | Generally safe. The GNSS timeseries case lives here. |
| $0.1$ to $0.5$ | Use with caution. Validate against the exact formula on a subset. |
| $< 0.1$ | Unreliable. Compute $Q_{\hat{e}_0}$ exactly. |

**The structural risk of the approximation:** By replacing $Q_{\hat{e}_0}$
with $Q_{yy}$, you ignore the geometric strength of the design matrix $A$. The
true denominator is smaller than $\sigma_m$ (because the projection removes
the part of the noise that is explained by the model). Your approximate $w_i$
is therefore systematically *smaller* than the true value. The consequence is a
higher risk of **Type II error** (missed detection): real outliers produce a
dampened test statistic and may not cross the threshold.

### 3.3  The Minimal Detectable Bias (MDB): How Sensitive Is the Test?

A question that naturally follows from the w-test is: *how large does an error
have to be before the test reliably catches it?* The answer is the **Minimal
Detectable Bias (MDB)** — the smallest bias in a single observation that the
w-test will detect with a specified power.

The MDB connects three things:
- Your **false alarm rate** $\alpha$ (probability of a wrong rejection).
- Your desired **detection power** $\gamma_0 = 1 - \beta$ (probability of a
  correct detection, where $\beta$ is the missed detection rate).
- The **geometry** of your design matrix $A$ and noise model $Q_{yy}$.

**How to compute it:**

1. From $\alpha$ and $\beta$, determine the required non-centrality parameter
   $\lambda_0$ using the non-central normal distribution. Typical values:
   $\alpha = 0.05$, $\beta = 0.20$ give $\lambda_0 \approx 2.8$.

2. The MDB for observation $i$ is:

$$|b_{MDB,i}| = \frac{\lambda_0}{\sqrt{c_i^T Q_{yy}^{-1} Q_{\hat{e}_0} Q_{yy}^{-1} c_i}}$$

The denominator is exactly the term that appears in the w-test denominator
— it encodes how well the geometry "sees" the $i$-th observation.

**Physical interpretation for GNSS:** If the MDB for a particular observation
is 3 mm, any error smaller than 3 mm will go undetected (with 80% probability)
even if you reject $H_0$. If the MDB is 25 mm, the test is geometrically weak
for that observation: errors up to 25 mm can hide undetected. Surveyors
compute MDBs *before* collecting data to identify geometrically weak
observations and improve the network design.

The MDB is the primary **Key Performance Indicator** for a geodetic quality
control system. It directly quantifies **internal reliability** — the system's
inherent capacity to spot anomalies — in physical units (millimetres) that are
immediately interpretable.

### 3.4  Multiple Hypotheses Testing: When the Model Itself Is Wrong

The w-test is designed for the case where the data contains isolated blunders.
But often the OMT rejection is structural — the functional model is simply
incomplete. A strong rejection combined with a systematic pattern in the
residuals (an oscillation, a step, an accelerating trend) almost always signals
missing physics, not bad data.

In this case, test alternative models directly:

1. **Define candidate hypotheses.** For a GNSS subsidence timeseries, you
   might propose: (a) linear + annual harmonic, (b) linear + semi-annual
   harmonic, (c) linear + step function at a specific epoch, (d) piecewise
   linear with a breakpoint.

2. **Compute the OMT for each.** After extending the design matrix $A$ to
   include each candidate term, re-run the estimation and compute $T$.

3. **Compare using information criteria.** Use AIC or BIC to penalise model
   complexity, selecting the model that best balances fit and parsimony.

4. **Return to the DIA cycle.** Once you extend the model, re-run the full
   Detection → Identification → Adaptation loop to confirm the expanded model
   passes the OMT.

---

## Part 4 — Adaptation: Repairing the Model

Once the source of the OMT rejection is identified, the model must be repaired.
Two strategies address fundamentally different problems.

### 4.1  Strategy A — Remove the Blunder

If the w-test clearly identifies an isolated observation with a physically
plausible explanation (antenna disturbance, equipment malfunction, phase
unwrapping error), remove it from the dataset and re-estimate. This is the
**data snooping** approach.

Apply it only when:
- The w-statistic is large and unambiguous.
- Physical or logbook evidence supports the removal.
- Removing the point does not substantially reduce redundancy.

**Warning:** Indiscriminate deletion introduces confirmation bias — you end up
keeping data that agrees with your hypothesis and discarding data that
contradicts it. Always document every removal decision.

### 4.2  Strategy B — Extend the Functional Model

If a block of data fails, or residuals show a coherent pattern, the problem is
structural. The correct response is to extend the design matrix $A$ to
represent the missing physics.

| Pattern in residuals | Physical cause | Model extension |
|---------------------|---------------|-----------------|
| Sinusoidal oscillation at 1 year period | Seasonal soil moisture / temperature | $+ c_2 \sin(2\pi t) + c_3 \cos(2\pi t)$ |
| Permanent step at epoch $t^*$ | Earthquake, antenna change, monument reset | $+ c_4 \cdot H(t - t^*)$ where $H$ is the Heaviside step |
| Accelerating downward trend | Groundwater extraction, consolidation | $+ c_5 t^2$ (or spline terms) |
| Increasing residual scatter | Instrument ageing, environmental change | Update $Q_{yy}$ (stochastic model) |

**GNSS example:** Initial model: $y_i = x_0 + v t_i$. OMT fails.
Residuals show a clear annual oscillation. Extended model:

$$y_i = x_0 + v t_i + A_1 \sin(2\pi t_i) + A_2 \cos(2\pi t_i)$$

Now $n = 4$, $r = 1091$. Re-run estimation, recompute $T$. If the OMT passes,
the annual signal was the missing physics.

### 4.3  The DIA Loop: Iterate Until the OMT Passes

The DIA cycle is not a one-pass pipeline. It is an iterative refinement loop:

```
         ┌──────────────────────────────────────────┐
         ▼                                          │
  [Estimate model]                                  │
         │                                          │
  [OMT: compute T]                                  │
         │                                          │
   T ≤ K_α? ──── Yes ──► DONE. Model accepted.     │
         │                                          │
        No                                          │
         │                                          │
  [w-test / Multiple Hypotheses Testing]            │
         │                                          │
  [Adapt: remove blunder or extend A] ──────────────┘
```

Each iteration tightens the model's representation of reality. The loop ends
only when the OMT accepts the null hypothesis, confirming that the residuals
contain nothing but noise consistent with $Q_{yy}$.

This iterative workflow embodies the Delft School principle:
**quality comes from validated models, not from blind confidence in a first fit.**

---

## Part 5 — Implementation in Your Signal Solver

When your timeseries signal solver runs an estimation, the following checklist
ensures OMT is properly integrated.

### 5.1  Standard Outputs for Every Estimation Run

| Output | Formula | Purpose |
|--------|---------|---------|
| Test statistic | $T = \hat{e}^T Q_{yy}^{-1} \hat{e}$ | Core OMT input |
| Degrees of freedom | $r = m - n$ | Defines the $\chi^2$ distribution |
| Critical value | $K_\alpha$ from $\chi^2(r)$ at chosen $\alpha$ | Decision threshold |
| Pass/Fail verdict | $T \leq K_\alpha$? | Model validity statement |
| Redundancy ratio | $r/m$ | Validates whether w-test approximation is safe |
| MDB per observation | $|b_{MDB,i}|$ | Internal reliability indicator |

### 5.2  Residual Diagnostics to Run Alongside the OMT

The OMT answers *whether* the model failed. Residual diagnostics answer *why*.
Always run these in parallel:

- **Time-series plot of residuals** — reveals systematic patterns (trends,
  oscillations, steps, outliers).
- **Autocorrelation function (ACF)** — detects temporal correlations violating
  the white-noise assumption in $Q_{yy}$.
- **Q-Q plot or Kolmogorov–Smirnov test** — checks whether residuals follow the
  assumed normal distribution.
- **Power Spectral Density (PSD)** — identifies periodic signals not captured
  by the functional model.

### 5.3  Decision Protocol

```
OMT passes AND residuals look random      → Accept model. Report results.
OMT passes BUT residuals show patterns    → Investigate: correlated noise in Q_yy?
OMT fails AND isolated w-test flags       → Remove blunder. Re-run DIA.
OMT fails AND systematic residual pattern → Extend A. Re-run DIA.
OMT fails AND no clear pattern            → Reassess Q_yy (stochastic misspecification).
```

### 5.4  Why You Cannot Skip This

Without the OMT:
- A straight-line fit through data with an undetected step will produce a
  biased velocity estimate. The bias looks like a legitimate signal.
- The parameter uncertainties reported by least-squares will appear valid but
  will be wrong, because they were computed under a misspecified model.
- Any downstream product — subsidence maps, infrastructure risk assessments,
  published velocity fields — inherits this silent, unquantified bias.

The OMT is not a formality. It is the guarantee that your parameter estimates
are grounded in a model that has been statistically tested against the data
that produced them.
