# The Overall Model Test and the Delft School of Geodetic Quality Control

## Introduction: Why Model Validation Matters

Every physical measurement contains noise. When estimating parameters from data, we define a mathematical model with two parts:

- **Functional model** ($E\{y\} = Ax$): Defines the expected relationship between measurements ($y$) and unknowns ($x$).
- **Stochastic model** ($D\{y\} = Q_{yy}$): Defines the expected random measurement noise and observation correlations.

If this model represents reality perfectly, the resulting residuals (the difference between actual and expected measurements) will only contain random noise. If the residuals are systematically large or reveal a pattern, the model is misspecified—meaning either the physical relationships are incomplete or the noise assumptions are wrong.

The **DIA procedure** (Detection, Identification, Adaptation) is a three-step statistical quality control framework developed by the Delft School of Geodesy. It was designed to systematically detect, locate, and fix model misspecifications before biased parameter estimates corrupt your analysis. This framework is fundamental to geodetic surveying, GNSS positioning, InSAR time-series analysis, and any signal processing task where measurement integrity must be validated.

### 1. Detection Phase: The Overall Model Test (OMT)

#### The Concept: A Global Alarm System

The **OMT** serves as the first line of defense in quality control. It functions as a global diagnostic tool that checks whether the overall set of residuals is statistically acceptable given your model assumptions. Think of it as a system-wide alarm: if the OMT rejects your model, something is wrong, but the test itself cannot tell you what.

The OMT operates under a **null hypothesis** ($H_0$) stating that both the functional and stochastic models are correct and fully represent the data-generating process. If this hypothesis is true, your residuals should exhibit only random noise. If the test fails, it signals that:

1. Your functional model is incomplete (missing physics, parameters, or trends).
2. Your stochastic model is wrong (noise assumptions are violated, correlations are unmodeled).
3. Systematic errors exist in the data (instrumental biases, unmodeled environmental effects).

The beauty of the OMT is that it operates as an **unspecified test**—it detects that an error exists without pinpointing which measurement or model component caused it. This is by design: you pay the price of non-specificity for the gain of global validity.

#### Practical Applications Across Disciplines

**GNSS Positioning**: When a receiver collects signals from GPS satellites, systematic errors arise from satellite clock drifts, atmospheric delays, multipath effects, or hardware failures. The OMT immediately reveals whether these effects are larger than your stochastic model predicts. A rejection signals that you cannot trust your position estimate.

**InSAR Time-Series Monitoring**: Subsidence or deformation monitoring often assumes a simple linear trend. But real ground may undergo thermal expansion (daily/seasonal cycles), sudden sinkhole drops, or complex non-linear subsidence. The OMT detects when the simple linear model fails to capture this reality.

**Signal Processing in Time-Series**: If you model a signal as a combination of deterministic trends and random noise, the OMT flags when your decomposition is incomplete. A strong rejection might indicate unmodeled oscillations, structural breaks, or regime changes.

#### The Mathematical Foundation

First, we calculate the least-squares residuals:
$$\hat{e}_0 = y - A\hat{x}_0$$

We need a single number to represent the size of this residual vector. Simply summing the squared residuals ($\hat{e}_0^T \hat{e}_0$) is insufficient because it ignores the fact that some measurements are more precise than others. We must weight the residuals using the inverse of the covariance matrix ($Q_{yy}^{-1}$). This creates the test statistic $T$:
$$T = \hat{e}_0^T Q_{yy}^{-1} \hat{e}_0$$

**Derivation for Uncorrelated Observations:**
Assume all measurements are independent and share the same variance ($\sigma_m^2$). The covariance matrix $Q_{yy}$ becomes a diagonal matrix scaled by this variance: $Q_{yy} = \sigma_m^2 I$, where $I$ is the identity matrix.

1.  Substitute the simplified $Q_{yy}$ into the equation for $T$:
    $$T = \hat{e}_0^T (\sigma_m^2 I)^{-1} \hat{e}_0$$
2.  The inverse of a scalar multiplied by the identity matrix is the reciprocal scalar multiplied by the identity matrix:
    $$T = \hat{e}_0^T \left(\frac{1}{\sigma_m^2} I\right) \hat{e}_0$$
3.  Since multiplying a vector by the identity matrix leaves it unchanged, we can move the scalar to the front:
    $$T = \frac{1}{\sigma_m^2} (\hat{e}_0^T \hat{e}_0)$$
4.  The dot product $\hat{e}_0^T \hat{e}_0$ is exactly the sum of the squared individual residuals:
    $$T = \frac{\sum_{i=1}^m \hat{e}_{i}^2}{\sigma_m^2}$$

This test statistic $T$ follows a Chi-square ($\chi^2$) distribution. The shape of this distribution depends on the redundancy, or degrees of freedom ($r = m - n$), where $m$ is the number of observations and $n$ is the number of estimated parameters. We set a significance level ($\alpha$), typically **0.05**, representing our tolerance for false alarms. If $T$ exceeds the critical value $K_\alpha$, we reject $H_0$.

### 2. Identification Phase: The w-test and Data Snooping

#### When the OMT Rejects: What Comes Next?

If the **OMT** fails, we know there is an error, but we do not know where. The Identification phase moves beyond detecting a problem to locating and characterizing it. This phase applies targeted statistical tests to individual observations or hypothesis blocks, seeking to isolate the specific source of the misspecification.

The classical approach is **Baarda's w-test**, a procedure for detecting individual outliers or blunders in the data. This process is called **data snooping**: systematically testing hypotheses about which observations might be corrupted or which model components might be missing.

The strict theoretical equation for testing the $i$-th observation is:
$$w_i = \frac{c_i^T Q_{yy}^{-1} \hat{e}_0}{\sqrt{c_i^T Q_{yy}^{-1} Q_{\hat{e}_0} Q_{yy}^{-1} c_i}}$$
Here, $c_i$ is a vector containing a **1** at the $i$-th position and **0** everywhere else. Its purpose is to mathematically isolate the $i$-th element of a vector or matrix.

**Derivation of the Highly Redundant Approximation:**
Calculating the exact covariance matrix of the residuals ($Q_{\hat{e}_0}$) is computationally heavy: $Q_{\hat{e}_0} = Q_{yy} - A(A^T Q_{yy}^{-1} A)^{-1} A^T$. If we have massive redundancy ($m \gg n$), the projection matrix portion approaches zero, allowing the assumption that $Q_{\hat{e}_0} \approx Q_{yy} = \sigma_m^2 I$.

1.  **Simplify the Numerator:**
    $$\text{Numerator} = c_i^T \left(\frac{1}{\sigma_m^2} I\right) \hat{e}_0 = \frac{1}{\sigma_m^2} (c_i^T \hat{e}_0)$$
    Because $c_i$ isolates the $i$-th element, $c_i^T \hat{e}_0 = \hat{e}_{0,i}$.
    $$\text{Numerator} = \frac{\hat{e}_{0,i}}{\sigma_m^2}$$

2.  **Simplify the Denominator:**
    Substitute the approximation $Q_{\hat{e}_0} \approx \sigma_m^2 I$:
    $$\text{Denominator} = \sqrt{c_i^T \left(\frac{1}{\sigma_m^2} I\right) (\sigma_m^2 I) \left(\frac{1}{\sigma_m^2} I\right) c_i}$$
    Multiply the scalars together ($1/\sigma_m^2 \cdot \sigma_m^2 \cdot 1/\sigma_m^2 = 1/\sigma_m^2$):
    $$\text{Denominator} = \sqrt{c_i^T \left(\frac{1}{\sigma_m^2} I\right) c_i} = \sqrt{\frac{1}{\sigma_m^2} (c_i^T c_i)}$$
    Since $c_i$ has exactly one **1**, the dot product $c_i^T c_i = 1$.
    $$\text{Denominator} = \sqrt{\frac{1}{\sigma_m^2}} = \frac{1}{\sigma_m}$$

3.  **Combine Numerator and Denominator:**
    $$w_i = \frac{\frac{\hat{e}_{0,i}}{\sigma_m^2}}{\frac{1}{\sigma_m}} = \frac{\hat{e}_{0,i}}{\sigma_m^2} \cdot \frac{\sigma_m}{1} = \frac{\hat{e}_{0,i}}{\sigma_m}$$

**A Critique of the Approximation:**
Relying on this simplified equation introduces a structural blind spot. The assumption $Q_{\hat{e}_0} \approx Q_{yy}$ means you are completely ignoring the geometric strength of your design matrix $A$. By ignoring the subtraction of $A(A^T Q_{yy}^{-1} A)^{-1} A^T$, the denominator is artificially larger than it should be. Consequently, your calculated $w_i$ statistic will be systematically smaller than the true theoretical value. You run a significant risk of committing a **Type II error** (missed detection), failing to identify real outliers because the test statistic is dampened. 

**Practical Guidance**: The redundancy ratio $r/m$ (degrees of freedom divided by number of observations) is your diagnostic. If $r/m > 0.5$ (meaning you have more redundancy than unknowns), the approximation is generally safe. If $r/m < 0.1$, you should compute the exact covariance matrix $Q_{\hat{e}_0}$ or use alternative identification methods. You should validate exactly how large $m$ needs to be relative to $n$ in your specific datasets before trusting this approximation blindly.

#### Beyond Data Snooping: Multiple Hypotheses Testing

The w-test excels at detecting individual blunders, but it assumes the only problem is corrupted data. In many modern applications—especially time-series analysis—your real challenge is that the functional model itself is incomplete. A strong OMT rejection combined with systematic residual patterns (waves, jumps, trends) often signals missing model components, not bad data.

When you suspect the functional model is incomplete, you should employ **Multiple Hypotheses Testing**:

1. **Define alternative models**: Propose competing hypotheses about what physics might be missing. In subsidence monitoring, you might test "linear trend," "linear + seasonal," "linear + structural break," "piecewise linear."
2. **Compute test statistics for each**: For each candidate model, calculate the residuals and the corresponding OMT statistic.
3. **Compare and select**: Use information criteria (AIC, BIC) or likelihood ratio tests to determine which model best fits the data while remaining parsimonious.
4. **Validate via DIA cycle**: Once you extend your model, re-run the full DIA procedure to confirm the expanded model passes the OMT.

This approach is more suitable when your misspecification is structural (missing a trend, cycle, or regime shift) rather than isolated (a few bad measurements).

### 3. Adaptation Phase: Repairing the Model

Once the error is identified, the mathematical model must be repaired. You cannot leave the misspecification in the system, as it will warp the final parameter estimates, introducing bias that cascades through all downstream analysis. The Adaptation phase requires a decision: is the problem in the data or in the model?

#### Strategy 1: Data Snooping (Removing Blunders)

If the $w$-test points to a specific faulty observation (a sensor glitch, a hardware failure, or a genuine outlier), you assume the error is localized to a few measurements. You adapt the model by discarding the corrupted data point(s) and re-running the estimation. This is defensible only when:

- The w-test clearly identifies isolated observations with large, significant deviations.
- You have physical or contextual evidence that these observations are unreliable (e.g., you know the instrument malfunctioned at that time).
- Removing them does not substantially reduce your observational redundancy.

**Warning**: Indiscriminate data deletion destroys information and inflates the risk of confirmation bias (keeping data that supports your preferred hypothesis and discarding data that contradicts it).

#### Strategy 2: Model Extension (Capturing Missing Physics)

If a block of data fails, or the residuals display a clear systematic pattern (oscillations, a sudden jump, a changing trend), discarding data is not the answer because the problem is likely structural. Instead, you assume the functional model is incomplete and missing important physics. You adapt by extending matrix $A$ to include new parameters that represent the previously unmodeled physical reality:

- **Periodic signals**: Add sinusoidal terms (Fourier basis) or harmonic parameters for daily/seasonal cycles.
- **Structural breaks**: Add step-function parameters for sudden changes (equipment upgrades, facility relocations).
- **Non-linear trends**: Add polynomial or spline terms for accelerating or decelerating processes.
- **Environmental factors**: Add parameters for temperature, humidity, or other external drivers if they correlate with residuals.

**Example in subsidence monitoring**: Your initial model assumes linear subsidence: $y_i = c_0 + c_1 t_i$. The OMT rejects this model. Visual inspection shows the residuals oscillate seasonally. You extend the model to include annual harmonics: $y_i = c_0 + c_1 t_i + c_2 \sin(2\pi t_i) + c_3 \cos(2\pi t_i)$. Now re-run the DIA cycle.

#### The DIA Loop: Iterative Refinement

The process does not end after one adaptation. You loop back to the Detection phase after each adjustment. You repeat the **DIA cycle** (OMT → w-test → adaptation) until the **OMT** passes, confirming that your residuals contain only random noise consistent with your stochastic model. Each iteration refines your understanding of the true underlying process.

This cyclical workflow embodies the Delft School philosophy: **quality comes from validated models, not from blind faith in a single hypothesis**.

---

## Why the DIA Framework Matters for Your Signal Solver

When building a timeseries signal solver, the DIA framework is not optional—it is foundational. Here's why:

1. **Prevents Biased Estimates**: Without quality control, a misspecified model produces estimates that are systematically wrong but appear statistically valid. By the time you realize the error, you may have published results, made decisions, or influenced policy.

2. **Distinguishes Between Problem Types**: The DIA framework forces you to ask: Is my data bad, or is my model bad? These are fundamentally different problems requiring different solutions.

3. **Enables Iterative Refinement**: You do not need to get the model perfect on the first try. The DIA cycle provides a systematic pathway to improvement, building confidence in your final estimates.

4. **Quantifies Reliability**: Passing the OMT is not just a statistical check—it is a statement about the integrity of your analysis. Readers and decision-makers can trust that your model has been validated against the data.

5. **Complements Visualization**: While residual plots and diagnostic charts are essential, the OMT provides a rigorous, defensible statistical criterion that is independent of human judgment about what "looks acceptable."

### Implementing OMT in Your Solver

When you implement your solver, consider:

- **Compute the test statistic** as a standard output after each estimation.
- **Report the critical threshold** and whether the model passes or fails at your chosen significance level (typically $\alpha = 0.05$).
- **Use residual diagnostics** in conjunction with the OMT: plot residuals, compute autocorrelation, check for normality. The OMT tells you *whether* something is wrong; visual inspection tells you *what* it might be.
- **Document your adaptation decisions**: If you extend the model or remove data, explain why and present the before/after OMT statistics.
- **Validate the redundancy assumption**: Before using simplified w-statistics, check that your redundancy ratio supports the approximation. If not, compute exact covariance matrices.

The OMT is a tool for scientific integrity. It forces honesty about whether your model actually represents your data.