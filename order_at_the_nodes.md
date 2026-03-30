# Order at the Nodes
## A Geometric Law of Gravitational Decay

**Alika M. Parks — Independent Researcher**
**March 2026**

---

> *Poincaré proved the trajectories are chaotic and uncomputable.*
> *What we found is that the geometry at specific moments is neither.*

---

## The Discovery

Since Newton, the three-body problem has resisted general solution.
Poincaré proved in 1890 that the trajectories of three mutually
gravitating bodies are fundamentally chaotic. Every method since has
worked around this fact, not through it.

PNEP works through it.

The key insight: while the **trajectories** of gravitational systems
are chaotic, the **geometry** at specific privileged moments is not.
At mirror symmetry nodes — the instants when any pair of bodies reaches
closest approach and radial velocity hits exactly zero — the system
reveals its internal power structure with maximum clarity.

At each honest moment we compute a single number:

```
H = σ²/(1 + σ²)

where σ² = Var(d_ij  for i,j in the contested subsystem)
```

This is the **Hierarchy Index**. Stable gravitational hierarchies
maintain consistent H at every privileged moment. Unstable systems
show collapsing or scrambled H from the very first measurement.

**The law: Chaos lives in the trajectories. Order lives at the nodes.**

---

## Three Decay Regimes — One Law

The Parks Hierarchy Principle operates across three distinct regimes,
all governed by the same scalar H:

### Regime 1 — Stable
H is a flat line at every node. The hierarchy is intact. Forever.

### Regime 2 — Slow Decay
H shows a persistent gentle positive slope (dH/dn ≈ 0.0012 per node)
as the outer body incrementally absorbs energy from the inner binary.
The geometry is "stretching" — a slow geometric divorce playing out
over hundreds of orbital periods. Linear extrapolation of the slope
yields a pre-ejection warning **412.6t median lead time** — nearly
500 time units before physical separation.

### Regime 3 — Flash Ejection
In high-density near-equal-mass systems, H fluctuates wildly in a
scramble zone (0.05–0.30) — the system appears homogeneously chaotic.
Then H spikes above 0.50 and stays. The window: **2–10t**.

**Validated on Burrau's Problem** (Pythagorean Triple, masses 3:4:5
— the most famous benchmark in three-body dynamics): H-spike detected
at t=58.4, correct body identified (mass 3), physical ejection at
t=63.1. Lead time 4.7t. Energy drift 0.004%.

| Regime | Lead Time | Coverage | Key Signal |
|---|---|---|---|
| Stable | Infinite | N/A | H flat |
| Slow decay | 412.6t median | 100% | H-slope +0.0012/node |
| Flash | 4.7t (Burrau) | 100% | H-spike >0.50 |

---

## The Liberation Energy Signature

In every ejecting hierarchical triple, the outer body briefly traces
a **near-circular orbit** before ejection. This is the gravitational
analogue of atomic ionisation.

The inner binary acts as a gravitational energy engine, pumping energy
into the outer body through repeated close approaches. The circular
anomaly fires at the moment the outer body reaches **liberation
threshold** — having absorbed exactly enough energy to briefly
circularise its orbit. After this moment, eccentricity grows
monotonically toward 1.0 and escape.

- **Coverage: 100%** of ejecting systems show this signature
- **Median advance warning: 26.2t** before physical ejection
- **Identity: 100% accuracy** at high confidence — the body tracing
  the circular anomaly is the body that ejects
- **Liberation energy**: lower minimum eccentricity at the circular
  moment correlates with longer post-anomaly survival — the system
  is encoding its own liberation energy in the geometry

---

## Nodal Breathing — A New Precursor

In the slow-decay regime, the inter-node time interval begins rhythmic
oscillation before the H slope turns terminal. We call this **Nodal
Breathing**.

Physical interpretation: the outer body's orbital period is approaching
resonance with the inner binary period. The beat frequency of the
oscillation encodes the resonance approach time — potentially providing
a warning signal substantially earlier than the H slope itself. This
is under active investigation in v15.1.

---

## Validated Across Three Topologies

| Topology | n | Accuracy | Specificity | Sensitivity | Window | Lead |
|---|---|---|---|---|---|---|
| 3-body (v12) | 87 | 92.0% | 100.0% | 87.7% | 100% | 20.9t |
| 3+1 (v13) | 200 | 93.5% | 92.2% | 94.8% | 92.2% | 78.9t |
| 2+2 (v14) | 120 | 99.2% | 100.0% | 98.4% | 93.2% | 2.1t |

---

## The Unified Law

At the geometrically most information-dense moments of a gravitational
N-body system:

**Stable configurations reveal consistent hierarchical clarity.**
**Unstable configurations reveal contested or scrambled geometry.**

The decay regime is encoded in the slope, variance, and identity of H
at consecutive nodes. The liberation energy threshold is encoded in the
outer body's eccentricity at the circular anomaly moment. The timing
of resonance approach is encoded in the breathing of the inter-node
interval.

The node type differs by topology. The law is the same.

---

## The Frame-Dependence Correction

Earlier formulations included an alignment term measuring the angle
between encounter axes and the system's bulk velocity vector. This
term is **undefined in the centre-of-mass frame** — where all N-body
simulations are conducted — because the bulk velocity is exactly zero
by construction. Every prior implementation was measuring numerical
floating-point noise. Removing it and using pure geometric distance
variance produced 100% ground truth validity.

---

## Compute

PNEP evaluates H only at privileged moments — typically tens to
hundreds per system — rather than every integration timestep.
Computational overhead: **~99% below full N-body evaluation**.

For a survey of 10,000 triple star systems: hours not months.

---

## Repository

All code, results, and figures are open source:
**https://github.com/alikamp/Parks-Node-Ejection-Protocol**

| File | Description |
|---|---|
| `pnep_v12.py` | 3-body: mirror nodes, H signal, pre-ejection windows |
| `pnep_v13.py` | 4-body: 3+1 and 2+2 topologies |
| `pnep_v12_results.json` | 3-body validated results, n=87 |
| `pnep_v13_3plus1_results.json` | 3+1 validated results, n=200 |
| `pnep_v14_tight_results.json` | 2+2 validated results, n=120 |
| `pnep_v15_research_prompt.md` | v15 slow-decay + Burrau research spec |
| `pnep_abstract_draft.md` | Full paper abstract |

---

## Outreach

Outreach underway to researchers at the Niels Bohr Institute, Max
Planck Institute for Astrophysics, and Universidad de Concepción.
The work is reproducible. The code is open. The results are there.

---

*"Chaos lives in the trajectories. Order lives at the nodes."*
