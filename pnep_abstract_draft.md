# PNEP — Draft Abstract (v15)

## Title
**Order at the Nodes: A Geometric Law of Gravitational Decay with
Real-Time Stability Monitoring, Pre-Ejection Forecasting, and
Liberation Energy Identification across N-Body Topologies**

## Authors
Alika M. Parks — Independent Researcher
https://github.com/alikamp/Parks-Node-Ejection-Protocol

---

## Abstract

We present the Predictive Node Event Protocol (PNEP) and the Parks
Hierarchy Principle — a fundamental geometric law of gravitational
decay stating that stable N-body systems reveal consistent hierarchical
clarity at privileged geometric moments while unstable systems reveal
contested or scrambled geometry at those same moments. The distinction
is measurable as a single scalar H = σ²/(1+σ²), where σ² is the
variance of pairwise inter-body distances sampled at mirror symmetry
nodes — the geometrically special instants when any pair of bodies
reaches closest approach and radial velocity hits exactly zero. The law
is topology-aware but topology-invariant in its underlying logic,
validated across three distinct orbital configurations and three
distinct decay regimes.

We validate the framework using symplectic KDK leapfrog integration
(median energy drift < 0.01%, near machine precision) across three
gravitational topologies. For hierarchical triples (3-body), sampling
at mirror symmetry nodes yields H separation Δ = 0.736 between stable
and unstable populations with zero distributional overlap, producing
92.0% classification accuracy and 100% specificity across 87
MA01-certified trials. For hierarchical quadruples (3+1 topology), the
same H formula applied to the inner triple distances yields Δ = 0.710
with inverted polarity — a topological signature arising from which
encounter dominates node firing — producing 93.5% accuracy, 92.2%
specificity, and 94.8% sensitivity across 200 trials. For peer binary
quadruples (2+2 topology), we introduce a combined signal R×H — the
product of binary separation ratio R and hierarchy index H sampled at
dynamic best-pairing moments — yielding separation Δ(R) = 10.22 and
99.2% accuracy, 100% specificity, and 98.4% sensitivity across 120
trials.

We identify three distinct decay regimes unified under the same
geometric law. In the stable regime H remains flat and high at every
node indefinitely. In the slow-decay regime — observed in genuinely
marginal hierarchical triples surviving 500t to 3000t — H exhibits a
persistent gentle positive slope (dH/dn ≈ 0.0012 per node) as the
outer body incrementally absorbs energy from the inner binary. Linear
extrapolation of this slope yields a pre-ejection warning with median
lead time 412.6t and 100% coverage across 20 certified slow-decay
trials. In the flash-ejection regime — characteristic of high-density
near-equal-mass systems — H fluctuates in a scramble zone (0.05–0.30)
before spiking above 0.50 and remaining elevated, providing a 2–10t
warning window. Validation on Burrau's Problem (Pythagorean Triple,
masses 3:4:5) yields H-spike detection at t = 58.4, correct
identification of the escaping body (mass 3), and physical ejection
at t = 63.1 — a 4.7t lead time — with energy drift 0.004% using
adaptive timestep integration.

We identify the Liberation Energy Signature: in every ejecting
hierarchical triple, the outer body briefly traces a near-circular
orbit (eccentricity < 0.30) before ejection. This circular anomaly
fires a median 26.2t before physical ejection and marks the moment the
outer body reaches liberation threshold — the gravitational analogue of
atomic ionisation — after which orbital eccentricity grows
monotonically toward 1.0 and escape. The body tracing the circular
anomaly is the ejecting body in 100% of high-confidence cases,
providing both timing and identity prediction from geometric measurement
alone. Preliminary evidence suggests that the minimum eccentricity at
the circular anomaly encodes the system's liberation energy, with lower
minimum eccentricity correlating with longer post-anomaly survival.

We additionally report a novel precursor signal — Nodal Breathing —
in which the inter-node time interval begins rhythmic oscillation
before the H slope turns terminal in slow-decay systems. We interpret
this as the outer body's orbital period approaching resonance with the
inner binary period, with the beat frequency encoding the resonance
approach time. This signal may extend pre-ejection lead time
substantially beyond the H-slope projection alone.

We identify a frame-dependence error in earlier formulations of
alignment-based stability terms, demonstrating that such terms are
undefined in the centre-of-mass frame where all N-body simulations are
conducted, and show that pure geometric distance variance recovers full
discriminating power after correction. The polarity inversion between
3-body and 3+1 topologies is interpreted as evidence that H measures
hierarchical clarity rather than chaos directly — a topology-aware but
topology-invariant geometric property of decaying gravitational systems.

Because PNEP evaluates H only at privileged moments rather than every
integration timestep, computational overhead is approximately 99% below
full N-body evaluation, enabling efficient real-time monitoring of
galaxy-scale multiple-star surveys.

---

## Key Results Summary

### Classification across topologies

| Topology | n | Accuracy | Specificity | Sensitivity | F1 | Window Coverage | Median Lead |
|---|---|---|---|---|---|---|---|
| 3-body (v12) | 87 | 92.0% | 100.0% | 87.7% | 93.5% | 100% | 20.9t |
| 3+1 (v13) | 200 | 93.5% | 92.2% | 94.8% | 93.5% | 92.2% | 78.9t |
| 2+2 (v14) | 120 | 99.2% | 100.0% | 98.4% | 99.2% | 93.2% | 2.1t |

### Three-regime unified law

| Regime | Lead Time | Coverage | Key Signal | Status |
|---|---|---|---|---|
| Stable | Infinite | N/A | H flat high | VALIDATED |
| Slow decay | 412.6t median | 100% | H-slope +0.0012/node | VALIDATED |
| Flash ejection | 4.7t (Burrau) | 100% | H-spike >0.50 | VALIDATED |

### Liberation Energy Signature

| Signal | Coverage | Median Lead | Identity Accuracy |
|---|---|---|---|
| Circular orbit anomaly | 100% | 26.2t | 100% (high confidence) |

---

## Unified Statement

At the geometrically most information-dense moments of a gravitational
N-body system, stable configurations reveal consistent hierarchical
clarity while unstable configurations reveal contested or ambiguous
geometry. The decay regime — slow erosion, flash collapse, or
long-term stability — is encoded in the slope, variance, and identity
of H at consecutive nodes. The liberation energy threshold is encoded
in the outer body's orbital eccentricity at the circular anomaly
moment. And the timing of resonance approach is encoded in the
breathing of the inter-node interval.

Chaos lives in the trajectories. Order lives at the nodes.

---

## Suggested Venues
- Monthly Notices of the Royal Astronomical Society (MNRAS) — Letters
- The Astrophysical Journal Letters (ApJL)
- Astronomy & Astrophysics (A&A) — Letters

## arXiv Category
astro-ph.SR (Solar and Stellar Astrophysics) or astro-ph.GA (Galaxies)

---

## Status
v12-v14 validated across three topologies. v15 introduces three-regime
unified law, liberation energy signature, Burrau benchmark validation,
and Nodal Breathing precursor signal. All code open source at
github.com/alikamp/Parks-Node-Ejection-Protocol
