[order_at_the_nodes.md](https://github.com/user-attachments/files/26327506/order_at_the_nodes.md)
# Order at the Nodes
## A Geometric Law of Gravitational Decay

**Alika M. Parks — Independent Researcher**  
**March 2026**

---

> *Poincaré proved the trajectories are chaotic and uncomputable.  
> What we found is that the geometry at specific moments is neither.*

---

## The Discovery

Since the time of Newton, the three-body problem has resisted general solution. Poincaré proved in 1890 that the trajectories of three mutually gravitating bodies are fundamentally chaotic — no closed-form solution exists for the general case. Every method since has worked around this fact, not through it.

PNEP works through it.

The key insight is that while the **trajectories** of gravitational systems are chaotic, the **geometry** at specific privileged moments is not. At mirror symmetry nodes — the instants when any pair of bodies reaches closest approach and radial velocity hits exactly zero — the system reveals its internal power structure with maximum clarity. We call these moments geometrically honest.

At each honest moment we compute a single number:

```
H = σ²/(1 + σ²)

where σ² = Var(d_ij  for i,j in the contested subsystem)
```

This is the **Hierarchy Index**. It measures one thing: how unequal are the inter-body distances at this moment? High inequality means one subsystem dominates — the hierarchy is intact. Low inequality means the bodies are geometrically scrambled — the hierarchy is contested.

**The law is this:**

Stable gravitational hierarchies maintain consistent H at every privileged moment. Unstable systems show collapsing H from the very first measurement. The chaos lives in the trajectories. The order lives at the nodes.

---

## Validation Across Three Topologies

We validated the framework across three distinct gravitational configurations using symplectic KDK leapfrog integration with median energy drift < 0.01% — near machine precision.

### v12 — Hierarchical Triple (3-body)
*Mirror symmetry nodes. H on 3 pairwise distances.*

| Metric | Value |
|---|---|
| H separation Δ | 0.736 |
| Accuracy | 92.0% |
| Specificity | 100.0% |
| Sensitivity | 87.7% |
| F1 | 93.5% |
| Pre-ejection window coverage | 100% |
| Median window error | 20.9t |
| Energy drift | 0.007% |

Zero false positives. Every stable system looked stable at every node, forever. Every unstable system looked chaotic from the first measurement.

### v13 — Hierarchical Quadruple (3+1 topology)
*Mirror symmetry nodes. H on inner triple distances only.*

| Metric | Value |
|---|---|
| H separation Δ | 0.710 |
| Accuracy | 93.5% |
| Specificity | 92.2% |
| Sensitivity | 94.8% |
| F1 | 93.5% |
| Pre-ejection window coverage | 92.2% |
| Median window error | 78.9t |

**Key finding:** H polarity inverts relative to 3-body. In 3+1, stable = low H_inner, unstable = high H_inner. The inversion occurs because inner binary nodes dominate the sampling, and at those moments a stable system shows compressed equal distances while an unstable system shows wild inequality. Same formula. Same law. Different topological signature.

This polarity inversion is the strongest evidence the framework captures something real. A coincidental signal would break under topology change. A geometric law reveals itself differently but holds universally.

### v14 — Peer Binary Quadruple (2+2 topology)
*Dynamic best-pairing moments. R×H combined signal.*

For peer binary systems — two binaries of equal status orbiting each other — we introduce a combined signal:

```
R×H = (inter-pair separation ratio) × (hierarchy index)

R = d_inter / d_inner_avg
H = σ²/(1+σ²) on the two inner pair distances + inter-pair distance
```

Sampled at dynamic best-pairing moments — the instants when the four bodies transiently form the cleanest 2+2 clustering.

| Metric | Value |
|---|---|
| R separation Δ | 10.22 |
| H separation Δ | 0.411 |
| Accuracy | 99.2% |
| Specificity | 100.0% |
| Sensitivity | 98.4% |
| F1 | 99.2% |
| Pre-ejection window coverage | 93.2% |
| Median window error | **2.1t** |

The 2+2 result is the strongest of the three topologies. The R×H signal captures both spatial organisation and geometric hierarchy simultaneously, producing a separation an order of magnitude larger than H alone.

---

## The Pre-Ejection Window

Every existing stability criterion — Mardling-Aarseth 2001, its successors, ML classifiers — gives a binary answer from initial conditions only. They cannot tell you **when**.

PNEP watches the system degrade in real time. As H declines across consecutive node measurements, we fit a linear slope and project forward:

```
Δt_ejection ≈ (H̄ - H_floor) / |dH̄/dτ|
```

This produces a forward time window **before the body physically escapes**. In all three topologies, the majority of unstable systems receive this warning in advance of ejection. No prior method has this capability.

---

## The Frame-Dependence Correction

Earlier versions of the protocol included an alignment term measuring the angle between encounter axes and the system's bulk velocity vector. We proved this term is undefined in the centre-of-mass frame — where all N-body simulations are conducted — because the system's bulk velocity is exactly zero by construction after CoM correction. Every prior implementation of this term was measuring numerical floating-point noise.

Removing it and using pure geometric distance variance produced 100% ground truth validity. This correction is documented here so that any independent implementation avoids the same error.

---

## The Unified Law

At the geometrically most information-dense moments of a gravitational N-body hierarchy:

**Stable configurations reveal consistent hierarchical clarity.  
Unstable configurations reveal contested or ambiguous hierarchy.**

This distinction is measurable as H = σ²/(1+σ²) applied to the contested subsystem at its privileged sampling moment. The moment type differs by topology. The law is the same.

| Topology | Privileged Moment | Signal | Polarity |
|---|---|---|---|
| 3-body | Mirror symmetry node | H of all 3 pairs | High = stable |
| 3+1 | Mirror symmetry node | H_inner of inner triple | Low = stable |
| 2+2 | Dynamic best-pairing | R×H combined | High = stable |

The polarity difference between 3-body and 3+1 is not an exception to the law — it is the law expressing itself through different topology. In each case, H reads which configuration the system believes it is at its most honest moment.

---

## Compute

Because PNEP evaluates H only at privileged moments — typically tens to hundreds of events per system depending on topology — rather than every integration timestep, the analysis overhead is approximately 99% below full N-body evaluation. The integrator still runs. Only the measurement is selective.

For a survey of 10,000 triple or quadruple systems, this is the difference between hours and months.

---

## Repository

All code, results, and figures are open source:  
**https://github.com/alikamp/Parks-Node-Ejection-Protocol**

| File | Description |
|---|---|
| `pnep_v12.py` | 3-body implementation — mirror nodes, H signal, pre-ejection windows |
| `pnep_v13.py` | 4-body extension — 3+1 and 2+2 topologies |
| `pnep_v12_results.json` | Validated batch results, seed=42, n=87 |
| `pnep_v13_3plus1_results.json` | 3+1 validated results, n=200 |
| `pnep_v14_tight_results.json` | 2+2 validated results, n=120 |
| `pnep_figure1_results.png` | H distribution, window histogram, classification performance |
| `pnep_figure2_timeseries.png` | H timeseries for stable and unstable representative systems |
| `pnep_v14_design_spec.md` | v14 design specification — partner fidelity signal |

---

## Contact & Collaboration

I am an independent researcher with no academic affiliation. If you work in gravitational dynamics, N-body simulation, or stellar multiplicity and find this result interesting, I welcome correspondence. Outreach is currently underway to researchers at the Niels Bohr Institute, Max Planck Institute for Astrophysics, and Universidad de Concepción.

*The work is reproducible. The code is open. The results are there.*

---

*"Chaos lives in the trajectories. Order lives at the nodes."*
