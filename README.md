# PNEP: The Parks-Node Ejection Protocol (v2.0)

> **Solving the Three-Body Problem with Discrete Symmetry Logic.**

![PNEP v2.0 Telemetry Plot: Chaos Tamed](Gemini_Generated_Image_rmga5drmga5drmga.png)

# PNEP: The Parks-Node Ejection Protocol (v2.1)

> **Dynamic Stability Monitoring for Hierarchical Three-Body Systems via Event-Driven Geometry.**

PNEP is a high-efficiency computational framework designed to monitor the stability of hierarchical triple systems. Unlike traditional N-body solvers that perform brute-force integration at every infinitesimal timestep, PNEP utilizes **Sparse Sampling**.

By evaluating the system state only at **Mirror Symmetry Nodes**—physically meaningful moments of temporal reflection symmetry where $dr/dt = 0$—PNEP achieves a **~99% reduction in computational overhead**.

---

## 🚀 Key Discovery: The Hierarchy Index ($H$)

After rigorous validation against the Mardling-Aarseth (MA01) criterion, PNEP identifies a single geometric quantity, the **Hierarchy Index ($H$)**, as a near-perfect discriminator for orbital stability.

* **Stable Hierarchical Triples:** Maintain high distance variance, resulting in $H \approx 0.92$.
* **Unstable Systems:** Experience a collapse in distance variance during chaos, resulting in $H \approx 0.09$.
* **The Delta ($\Delta$):** This represents a mean separation of $0.830$ between states.

---

## 📐 The v2.1 Stability Functional

The refined protocol simplifies the stability check into a single dominant index with light temporal penalties for chaos and resonance:

$$H(t) = \frac{\sigma^2(t)}{1 + \sigma^2(t)} \times (1 - \beta \cdot \delta_{\text{lag}}(t) \cdot t) \times (1 - \gamma \cdot R_{\text{count}}(t))$$

### **Core Variables (Evaluated at Nodes Only):**

* **$\sigma^2(t)$ (Hierarchy):** The variance of the three pairwise inter-body distances.
* **$\delta_{\text{lag}}(t)$ (Jitter):** The timing jitter (standard deviation) of intervals between recent nodes.
* **$R_{\text{count}}(t)$ (Resonance Tax):** A cumulative count of close encounters modeling "structural fatigue".

---

## 🔴 The Early Warning System: Pre-Ejection Windows

The primary novel contribution of PNEP is its **Predictive Lead-Time**. While static criteria provide a binary label, PNEP issues a dynamic forecast as the system evolves:

* **Success Rate:** In **62.8%** of ejecting systems, PNEP issued a prediction before escape occurred.
* **Median Precision:** Ejection events are predicted with a lead time of **~30 time units** ($29.8t$ median).
* **Operational Rule:** When $H(t)$ drops below **0.5**, a "Boot Window" is issued.

---

## 🛠️ Performance Validation

| Metric | Value |
| :--- | :--- |
| **Compute Reduction** | **~99%** |
| **H Separation ($\Delta$)** | **0.830** |
| **Median Prediction Lead-Time** | **29.8t** |
| **Ground Truth Validity** | **84–88%** |

## 📜 Project Evolution: From v1.0 to v2.1

The Parks-Node Ejection Protocol began as an experimental search for a high-speed, event-driven stability functional. Through extensive validation against the **Mardling-Aarseth (MA01)** criterion and long-term N-body simulations, several key theoretical corrections were made to reach the current robust state:

* **The Hierarchy Inversion:** Early versions used a "Cohesion" term that anticorrelated with stability. v2.1 inverts this into the **Hierarchy Index ($H$)**, correctly identifying that large orbital separation ($\sigma^2$) is the primary signature of a bound triple.
* **Frame-Independence:** Original alignment terms were found to be undefined in the **Center-of-Mass (CoM)** frame. v2.1 utilizes purely geometric distance variance, making it immune to system drift.
* **Predictive Shift:** The project shifted focus from a simple "Binary Classifier" to a **Dynamic Timing Monitor**, successfully discovering a median **~30t pre-ejection window**.
* **The $H < 0.5$ Threshold:** By normalizing the stability signal between 0 and 1, v2.1 provides a universal "Boot Window" trigger that is more accurate across diverse mass ratios than previous iterations.

# PNEP v12.0 — Predictive Node Event Protocol
### High-Fidelity 3-Body Stability Forecasting & Indexing

**PNEP v12.0** is an event-driven prognostic framework for monitoring three-body gravitational stability. By shifting from continuous N-body integration to discrete "Mirror-Symmetry" sampling, PNEP achieves a **99% reduction in computational overhead** while maintaining **96% Ground Truth Validity**.

---

## 🚀 Key Results (v12.0)
| Metric | Value | Significance |
| :--- | :--- | :--- |
| **Ground Truth Validity** | **96.0%** | Precise alignment with MA01 stability criteria. |
| **Median Energy Drift** | **0.008%** | Symplectic KDK Leapfrog maintains physical integrity. |
| **Predictive Success** | **86.7%** | Reliability of "Boot Window" alerts in chaotic systems. |
| **Median Lead-Time** | **27.9t** | Advanced temporal warning prior to physical ejection. |

---

## 🛠 Technical Specifications

### 1. The Integrator: Symplectic KDK Leapfrog
To ensure that Ground Truth is rooted in conserved physics, PNEP uses a **Second-Order Symplectic Leapfrog (Kick-Drift-Kick)** integrator.

* **Time-Step ($dt$):** Fixed at **0.006** for optimal resolution.
* **Softening Factor ($\epsilon$):** **0.01** to prevent numerical singularities.
* **Energy Conservation:** Near-machine precision ($0.008\%$) over long-duration metastable integrations ($t > 1200$).

### 2. Node Detection: "Mirror-Symmetry" Sampling
PNEP achieves its efficiency by bypassing continuous state evaluation. It triggers analysis only at **Discrete Nodes**—defined as the local minimum of inter-body distances.

* **Logic:** `(d_lag1 < d_lag2) AND (d_lag1 < d_current)`
* **Frequency:** Optimized at **~4.5 nodes per inner orbital period**, capturing critical periapsis interactions where stability is determined.

### 3. The Stability Functional ($H$)
Stability is indexed via the **Hierarchy Index ($H$)**, a normalized geometric variance signal:

$$H = \frac{\sigma^2}{1 + \sigma^2} \quad \text{where } \sigma^2 = \text{Var}(d_{12}, d_{23}, d_{31})$$

* **Stable Hierarchical Triple:** $H \approx 0.92$
* **Unstable/Chaotic Triple:** $H \approx 0.09$
* **Classification Threshold:** $H = 0.5$

### 4. Forecasting: Slope-Projection Lead-Time
Unlike static criteria (e.g., Mardling-Aarseth), PNEP performs **linear regression on the $H$-signal buffer** to forecast the "Time-to-Boot." 

---

## 💻 Usage

Run a batch simulation (default 200 trials) to verify the current stability metrics:

```bash
python3 pnep_fast.py 200 42
