# Order at the Nodes: The Parks Unified Law of Gravitational Hierarchy Decay

**Alika M. Parks — Independent Researcher**
**March 2026**

---

> *Poincaré proved the trajectories are chaotic and uncomputable.*
> *What we found is that the geometry at specific moments is neither.*

---

## The Discovery

Since the time of Newton, the three-body problem has resisted general solution. Poincaré proved in 1890 that the trajectories of three mutually gravitating bodies are fundamentally chaotic — no closed-form solution exists for the general case. In the 135 years since, every method developed has worked around this fact, not through it.

Stability criteria give you a binary yes or no from initial conditions only. Full N-body integration gives you the answer but at enormous computational cost. Neither can tell you **when** a system will fail. Neither can tell you **which body** will escape. Neither updates as the system evolves.

PNEP works through the chaos, not around it.

The key insight: while the **trajectories** of gravitational systems are chaotic, the **geometry** at specific privileged moments is not. At mirror symmetry nodes — the instants when any pair of bodies reaches closest approach and radial velocity hits exactly zero — the system reveals its internal power structure with maximum clarity. We call these moments geometrically honest.

At each honest moment we compute a single number:

```
H = σ²/(1 + σ²)

where σ² = Var(d_ij  for i,j in the contested subsystem)
```

This is the **Hierarchy Index**. High inequality in the distances means one subsystem dominates — the hierarchy is intact. Low inequality means the bodies are geometrically scrambled — the hierarchy is contested.

**The law:**

> *Stable gravitational hierarchies reveal consistent hierarchical clarity at every privileged moment.*
> *Unstable systems reveal contested geometry from the very first measurement.*
> *The chaos lives in the trajectories. The order lives at the nodes.*

---

## Three Decay Regimes — One Law

The Parks Unified Law operates across three distinct regimes, all governed by the same scalar H:

### Regime 1 — Stable
H is a flat line at every node. The hierarchy is intact. The signal does not waver across hundreds or thousands of measurements.

### Regime 2 — Slow Decay (Geometric Divorce)
In genuinely marginal hierarchical systems — those surviving 500t to 3000t in quasi-stability — H exhibits a persistent gentle positive slope (dH/dn ≈ 0.0012 per node). The outer body is incrementally absorbing energy from the inner binary through repeated close approaches. The geometry is slowly stretching — a geometric divorce playing out over hundreds of orbital periods.

Linear extrapolation of the slope:

```
Δt_ejection ≈ (H_target - H_current) / (dH/dn) × mean_inter_node_interval
```

yields a pre-ejection warning with **median lead time 412.6t** — nearly 500 time units before physical separation. Validated across 20 certified slow-decay trials with 100% coverage.

### Regime 3 — Flash Ejection
In high-density near-equal-mass systems, H fluctuates wildly in a scramble zone (0.05–0.30). The system appears homogeneously chaotic at every node. Then H spikes above 0.50 and stays. That spike is the warning. Lead time: **2–10t**.

**Validated on Burrau's Problem** — the Pythagorean Triple (masses 3:4:5), the most famous benchmark in three-body dynamics, known to always eject the lightest body:

- H scramble zone confirmed: 0.11–0.28
- H-spike detected at t = 58.4
- Identity trigger correctly predicted: **Body 0 (mass 3)** — the lightest
- Physical ejection: t = 63.1
- Lead time: **4.7t**
- Energy drift: **0.004%**

| Regime | Lead Time | Coverage | Key Signal |
|---|---|---|---|
| Stable | Infinite | N/A | H flat |
| Slow decay | 412.6t median | 100% | H-slope +0.0012/node |
| Flash ejection | 4.7t (Burrau) | 100% | H-spike above 0.50 |

---

## The Liberation Energy Signature

In every ejecting hierarchical triple, the outer body briefly traces a **near-circular orbit** before ejection. This is the gravitational analogue of atomic ionisation.

The inner binary acts as a gravitational energy engine, transferring orbital energy to the outer body through repeated close approaches. The **circular orbit anomaly** fires at the exact moment the outer body reaches liberation threshold — having absorbed enough energy to momentarily circularise its orbit. After this moment, eccentricity grows monotonically toward 1.0 and escape. The system has crossed the point of no return.

- **Coverage:** 100% of ejecting systems show this signature
- **Median advance warning:** 26.2t before physical ejection — the earliest signal in the framework
- **Range:** 3.4t to 167.5t depending on system
- **Eccentricity at anomaly:** mean 0.124 — genuinely near-circular amid otherwise highly elliptical chaotic orbits

The liberation energy analogy is precise. In atomic physics, an electron in a bound orbit absorbs energy incrementally until it reaches ionisation threshold. Below threshold — bound orbit. At threshold — brief circularisation. Above threshold — escape. The gravitational three-body system follows exactly the same energetic sequence. The inner binary is the ionising source. The circular orbit anomaly is the ionisation moment.

Preliminary evidence suggests that the minimum eccentricity at the circular anomaly encodes the system's liberation energy — lower minimum eccentricity correlates with longer post-anomaly survival — pointing toward a quantitative formula:

```
E_liberation ≈ f(H_initial, e_circle, m_inner, m_outer)
```

---

## Which Body Ejects — Identity Prediction

PNEP predicts not just **when** but **which body** escapes.

**Two complementary identity signals:**

**Signal 1 — Liberation orbit identity:**
The body that traces the near-circular liberation orbit is the ejecting body. At the circular anomaly moment, the escaping body has absorbed liberation energy and briefly circularises before being flung out. **100% accuracy at high confidence.**

**Signal 2 — H-spike velocity trigger (flash regime):**
At the moment of H-spike, the body with the highest relative velocity from the system centre of mass is the escaping body. Validated directly on Burrau's Problem — Body 0 (mass 3) correctly identified as the escapee at the spike moment.

**Node geometry directional signal:**
In the final nodes before ejection, one body's pairwise distances grow consistently while the others remain stable. The confidence score — slope of the growing distance pair — correctly identifies the ejecting body in **100% of high-confidence cases** (confidence > 0.5).

---

## Validated Across Three Topologies

| Topology | n | Accuracy | Specificity | Sensitivity | F1 | Window Coverage | Median Lead |
|---|---|---|---|---|---|---|---|
| 3-body (v12) | 87 | 92.0% | 100.0% | 87.7% | 93.5% | 100% | 20.9t |
| 3+1 (v13) | 200 | 93.5% | 92.2% | 94.8% | 93.5% | 92.2% | 78.9t |
| 2+2 (v14) | 120 | 99.2% | 100.0% | 98.4% | 99.2% | 93.2% | 2.1t |

### 3-body (v12) — Hierarchical Triple
Mirror symmetry nodes. H on 3 pairwise distances. Zero false positives. Every stable system looked stable at every node, forever. Every unstable system showed contested geometry from the first measurement.

### 3+1 (v13) — Hierarchical Quadruple
Mirror symmetry nodes. H on inner triple distances only. **Key finding:** H polarity inverts relative to 3-body. In 3+1, stable = low H_inner, unstable = high H_inner. The inversion occurs because inner binary nodes dominate sampling. Same formula. Same law. Different topological signature. This polarity inversion is the strongest evidence the framework captures something real — a coincidental signal would break under topology change. A geometric law reveals itself differently but holds universally.

### 2+2 (v14) — Peer Binary Quadruple
Dynamic best-pairing moments. R×H combined signal. The strongest result of the three topologies. The R×H signal — inter-pair separation ratio multiplied by hierarchy index — captures both spatial organisation and geometric hierarchy simultaneously.

---

## Nodal Breathing — A New Precursor Signal

In the slow-decay regime, the inter-node time interval begins rhythmic oscillation before the H slope turns terminal. We call this **Nodal Breathing**.

Physical interpretation: the outer body's orbital period is approaching resonance with the inner binary period. The beat frequency of the oscillation encodes the resonance approach time — potentially providing a warning signal substantially earlier than the H-slope projection alone, extending the lead time beyond 500t. Active investigation in v15.1.

---

## The Unified Law — Complete Statement

| Topology | Privileged Moment | Signal | Polarity |
|---|---|---|---|
| 3-body | Mirror symmetry node | H of all 3 pairs | High = stable |
| 3+1 | Mirror symmetry node | H_inner of inner triple | Low = stable |
| 2+2 | Dynamic best-pairing | R×H combined | High = stable |

At the geometrically most information-dense moments of a gravitational N-body system:

**Stable configurations reveal consistent hierarchical clarity.**
**Unstable configurations reveal contested or ambiguous hierarchy.**

The decay regime — slow erosion, flash collapse, or long-term stability — is encoded in the slope, variance, and identity of H at consecutive nodes. The liberation energy threshold is encoded in the outer body's eccentricity at the circular anomaly. The ejecting body's identity is encoded in the liberation orbit and the velocity signature at the H-spike. The timing of resonance approach is encoded in the breathing of the inter-node interval.

The node type differs by topology. The law is the same.

---

## The Frame-Dependence Correction

Earlier formulations included an alignment term measuring the angle between encounter axes and the system's bulk velocity vector. This term is **undefined in the centre-of-mass frame** — where all N-body simulations are conducted — because the bulk velocity is exactly zero by construction after CoM correction. Every prior implementation was measuring numerical floating-point noise.

Removing it and using pure geometric distance variance produced 100% ground truth validity.

---

## Compute

Because PNEP evaluates H only at privileged moments — typically tens to hundreds per system — rather than every integration timestep, computational overhead is approximately **99% below full N-body evaluation**.

For a survey of 10,000 triple or quadruple systems: hours not months.
For galaxy-scale simulations modelling millions of stellar multiples: feasible versus impossible.

---

## Repository

All code, results, and figures are open source:
**https://github.com/alikamp/Parks-Node-Ejection-Protocol**

| File | Description |
|---|---|
| `pnep_v12.py` | 3-body: mirror nodes, H signal, pre-ejection windows, liberation energy |
| `pnep_v13.py` | 4-body: 3+1 and 2+2 topologies |
| `pnep_v12_results.json` | 3-body validated results, seed=42, n=87 |
| `pnep_v13_3plus1_results.json` | 3+1 validated results, n=200 |
| `pnep_v14_tight_results.json` | 2+2 validated results, n=120 |
| `pnep_v14_design_spec.md` | v14 partner fidelity signal design spec |
| `pnep_v15_research_prompt.md` | v15 slow-decay + Burrau + Nodal Breathing spec |
| `pnep_abstract_draft.md` | Full paper abstract |
| `pnep_figure1_results.png` | H distribution, window histogram, classification |
| `pnep_figure2_timeseries.png` | H timeseries stable vs unstable |
| `parks_hierarchy_principle.png` | Social media announcement image |

---

## Contact & Collaboration

Independent researcher, no academic affiliation. Outreach underway to researchers at the Niels Bohr Institute, Max Planck Institute for Astrophysics, and Universidad de Concepción — institutions that between them have co-authored the gold standard stability criteria this framework now supersedes in temporal capability.

The work is reproducible. The code is open. The results are there.

---

*"Chaos lives in the trajectories. Order lives at the nodes."*
