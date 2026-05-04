# Understanding Statistical Metrics in Geodetic Timeseries

If you are looking at a Timeseries Decomposition Report and seeing terms like **p-value**, **w-stat**, and **Sigma ($\sigma$)**, don't worry! These are just fancy ways of asking two simple questions:
1. **The Big Picture:** Does my overall model (the line) match my data (the dots)?
2. **The Specifics:** Are there any specific dots that are "weird" or shouldn't be there?

---

## 1. The Foundation: Sigma ($\sigma$)
Before we talk about the stats, we have to talk about **Sigma**. 
In our project, Sigma is your **"Expectation of Messiness."**

Imagine you are measuring your height every morning. You know your ruler isn't perfect—sometimes you stand a bit flatter, sometimes your hair is poofy. You expect an error of maybe **2 mm**. That 2 mm is your **Sigma**.

*   If your data jumps by **50 mm** tomorrow, you know something happened (maybe you put on shoes!).
*   If your data jumps by **0.1 mm**, you don't care; it's within your "Expectation of Messiness."

---

## 2. The `p-value`: The "Overall Fit" Test
The `p-value` comes from the **Overall Model Test (OMT)**. It looks at every single data point at once and asks: *"Is this model (the trend and seasonal lines) good enough to explain all these dots, given our Sigma?"*

### The Analogy: The Tailor
Imagine you are a tailor making a suit for a customer.
*   **The Data:** The customer's body shape.
*   **The Model:** The suit you sewed.
*   **The p-value:** The score the customer gives you on the fit.

*   **High p-value ($\ge 0.05$):** The customer says, "This fits! Any tiny wrinkles are just because fabric isn't perfect." (The model is **Accepted**).
*   **Low p-value ($< 0.05$):** The customer says, "This is way too tight! You clearly missed a measurement." (The model is **Rejected**).

**In your report:** When Sigma was 8.0 mm, the $p$-value was 0.00 (Rejected). When you increased Sigma to 10.5 mm, the $p$-value jumped to 0.20 (Accepted). This means 10.5 mm is the "size" of expectation that finally makes the data fit the model.

---

## 3. The `w-stat`: The "Troublemaker" Test
While the $p$-value looks at the whole suit, the `w-stat` (w-test) looks at **one specific button** or **one specific stitch**. It checks every single day in your timeseries individually.

### The Analogy: The Classroom
Imagine a class of 100 students taking a test. Most students score between 70 and 80.
*   **The Model:** We expect students to score around 75.
*   **The w-stat:** How many "steps" away from the average is a specific student?

If one student scores a **12**, their `w-stat` will be very high (like -10.0). They are an **Outlier**.

*   **Threshold 3.29:** We use a "magic number" of 3.29. If a `w-stat` is higher than 3.29 (or lower than -3.29), it means that specific data point is so far away from the model that there is a **99.9% chance it isn't just random noise.** It’s an anomaly.

**In your report:** On `2022-06-17`, the `w-stat` was **11.92**. Even though we added a "Jump" to the model for that day, the data point was still nearly **12 times** further away than we expected! That day is a major "troublemaker" in the dataset.

---

## 4. How They Work Together (The DIA Loop)
The software uses a logic called **DIA**: **D**etection, **I**dentification, **A**daptation.

1.  **Detection (p-value):** The software runs the model. If the $p$-value is low, it screams: *"Hey! The model doesn't fit!"*
2.  **Identification (w-stat):** It looks through all the `w-stats` to find the biggest troublemaker. *"Aha! On June 17th, the error is huge!"*
3.  **Adaptation:** It changes the model. It adds a "Jump" or a "Seasonal Wave" at that spot and tries again.

It keeps doing this until the **p-value** finally says "Accepted!"

---

## Summary for your exams:
*   **Sigma ($\sigma$):** How much noise you are willing to ignore.
*   **p-value:** Tells you if the **entire model** is a lie ($p < 0.05$) or the truth ($p \ge 0.05$).
*   **w-stat:** Points a finger at **one specific day** that is acting weird. If it's over **3.29**, that day is officially an "Anomalous Observation."
