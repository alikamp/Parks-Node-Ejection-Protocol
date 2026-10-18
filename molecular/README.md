# PNEP in Molecular Dynamics — a scope test

Does the PNEP hierarchy index transfer from gravitational few-body systems to
molecular dissociation?

**Partly.** It works in a narrow, physically identifiable regime, and fails
outside it for a structural reason. This directory documents both.

Everything below is reproducible from `pnep_molecular.py` (NumPy + SciPy only):

```
python pnep_molecular.py            # all experiments
python pnep_molecular.py static     # one of: static, regime, thermal
```

---

## The system

Lennard-Jones trimer in reduced units, three distinct species, Lorentz-Berthelot
mixing for cross-pair parameters. Velocity-Verlet, dt = 0.002. Dissociation is
declared when any atom exceeds 8 length units from the centre of mass.

---

## 1. Raw H does not transfer

| configuration | U | H | H_log | min_u |
|---|---|---|---|---|
| compact trimer (bound, non-hierarchical) | −1.92 | 0.000 | 0.010 | −0.101 |
| dimer + monomer (bound, hierarchical) | −0.78 | 0.798 | 0.368 | 0.000 |
| **crushed pair (anti-bound)** | **+197.51** | **0.839** | 0.500 | **−0.376** |
| dissociated | ~0 | 0.000 | 0.010 | 8.296 |

A configuration at +197 energy — atoms being violently repelled — scores
*higher* on H than a bound hierarchical dimer.

The cause is structural, not a tuning problem. Variance of separations is
sign-blind: it cannot distinguish "close because bound" from "close because
crushed". In gravity that ambiguity cannot arise, because the potential is
monotonic and closer always means more bound. Lennard-Jones is non-monotonic,
so the ambiguity is the dominant failure mode.

Rescaling the distances does not help. Normalising by sigma_ij or by r_min
leaves the gap between bound and anti-bound at 0.03–0.13 — the information is
not in the distances at any scale.

**What recovers it:** the signed displacement from mechanical equilibrium,

```
u_ij = (r_ij − r_min,ij) / r_min,ij        r_min,ij = 2^(1/6) · sigma_ij
```

which is zero where net pair force vanishes, negative inside the repulsive
wall, positive in the attractive tail. `min(u)` is −0.376 for the crushed
configuration and 0.000 for the bound one. Lorentz-Berthelot mixing enters
here — not as a correction to H, but because it is what supplies r_min,ij per
pair.

**Consequence:** molecular systems need two scalars where gravity needed one.
Variance of u measures hierarchy; the minimum of u measures physicality. The
monotonic gravitational potential collapses these into a single quantity.

---

## 2. Validity regime

Define a reduced temperature against binding depth:

```
tau = (3/2) · N · T / |U_bind|
```

Lead time from the sustained H > 0.5 crossing to physical dissociation:

| T₀ | tau | n | median lead | min lead | CV |
|---|---|---|---|---|---|
| 0.40 | 0.94 | 7 | 18.34 | 10.24 | **0.307** |
| 0.50 | 1.17 | 9 | 15.08 | 5.74 | 0.437 |
| 0.60 | 1.41 | 9 | 12.32 | 6.16 | 0.782 |
| 0.80 | 1.87 | 7 | 11.56 | 7.08 | 0.649 |
| 1.00 | 2.34 | 4 | 17.58 | 5.92 | 0.385 |
| 1.30 | 3.04 | 4 | 11.93 | 5.42 | 0.989 |

Across all temperatures pooled, CV is 1.206 and the minimum lead is 0.4 — not a
forecast. Restricted to tau ≈ 1, CV falls to 0.307 and the minimum lead rises
to 10.24.

The reason is the same one that makes PNEP work gravitationally. Near threshold,
dissociation is **secular**: energy accumulates over many encounters and the
geometry drifts measurably beforehand. Well above threshold it is **ballistic** —
once kinetic energy exceeds the well, the atom leaves, and the kick that frees
it is the last event rather than a precursor. There is nothing to detect early
because nothing builds.

Cell counts are small (n = 4–9). The trend is consistent with the mechanism but
the individual CV values should not be read as precise.

---

## 3. No temperature correction is needed

Median H measured on bound (non-dissociating) trajectories:

```
tau 0.23 → 0.060     tau 0.82 → 0.059
tau 0.35 → 0.055     tau 0.94 → 0.085
tau 0.47 → 0.045     tau 1.05 → 0.056
tau 0.59 → 0.057     tau 1.17 → 0.058
tau 0.70 → 0.077
```

floor = 0.061 ± 0.011 over a 5× temperature range; Pearson r = +0.303,
p = 0.428. No significant trend. The 0.5 decision threshold sits roughly 38
standard deviations above the floor.

H depends on positions only, and within the well the potential geometry governs
the distance distribution far more strongly than thermal velocity does.
Temperature determines **which regime you are in**, not what H reads.

---

## 4. H does not beat trivial baselines

The test that matters for any claim of utility: does H outperform simply
watching the largest or mean pairwise separation? 240 trajectories at
tau = 0.94-1.29, 187 dissociations.

**As an early-warning signal** (same sustained-crossing rule, threshold tuned
per signal):

| signal | threshold | n | median lead | min lead | CV |
|---|---|---|---|---|---|
| H | 0.20 | 176 | 17.25 | 4.35 | 0.688 |
| max pairwise distance | 2.00 | 176 | 17.62 | 4.55 | **0.673** |
| mean pairwise distance | 1.50 | 176 | 17.27 | 4.70 | **0.672** |

**As a hazard estimator** (split-half validated, quantile-binned):

| delta | signal | Brier | AUC | calib. error |
|---|---|---|---|---|
| 10 | H | 0.0675 | 0.942 | 1.3% |
| 10 | max distance | 0.0650 | 0.947 | 1.3% |
| 10 | mean distance | **0.0645** | **0.950** | 1.3% |

H is marginally *worse* than both baselines on every metric at every time
horizon. In a three-body system it adds nothing: the variance of three
distances carries no information beyond their mean.

**This is the central negative result of this directory.** The node-scalar
approach does transfer to molecular dissociation in the near-threshold
regime, but the specific scalar is not what makes it work.

---

## 5. What does work: a calibrated dissociation hazard

Setting aside which signal feeds it, the hazard model itself is sound and
generalises out of sample:

| H bin | P(dissociate within 20) train | test | n_test |
|---|---|---|---|
| 0.0-0.1 | 5.4% | 4.8% | 30370 |
| 0.3-0.4 | 27.3% | 27.4% | 3034 |
| 0.5-0.6 | 39.7% | 42.9% | 2026 |
| 0.8-1.0 | 86.7% | 87.6% | 22536 |

Mean calibration error 2.4%, AUC 0.915. At delta = 10 the top-to-bottom
hazard ratio is 87x with 1.3% calibration error.

This is the useful output for adaptive sampling: a population-level statement
("27% of configurations in this band dissociate within 20 time units") that
holds on held-out data. It is a weaker claim than per-trajectory prediction and
a more reliable one. But it should be built on mean pairwise distance, not H.

---

## 6. LJ-13 cluster: a state variable, not a forecast

Three bodies give one triple, so there is no distribution to examine and
"watch the largest separation" is a free baseline. A 13-atom icosahedral
cluster removes both problems: H is computed over all C(13,3) = 286 triples
and aggregated, and no single distance is the obvious one to monitor.
`cluster13.py`, `cluster_predict.py`.

### H and the Lindemann parameter detect different transitions

| T0 | Lindemann | H_mean | evaporated |
|---|---|---|---|
| 0.50 | 0.084 | 0.097 | 0.0 |
| 0.70 | **0.188** | 0.107 | 0.0 |
| 0.90 | 0.298 | 0.178 | 0.3 |
| 1.10 | 0.284 | 0.305 | 1.0 |
| 1.60 | 0.439 | 0.495 | 2.3 |
| 2.00 | 0.503 | 0.697 | 3.7 |
| 2.50 | **0.506** | **0.753** | 5.0 |

At T0 = 0.70 Lindemann reads 0.188 — melted by the classical 0.1 criterion —
while H has barely moved from its solid value and nothing has evaporated.
Lindemann calls the solid-liquid transition; H ignores it. Above T0 = 2.0
Lindemann saturates (0.503 to 0.506 while evaporation goes 3.7 to 5.0) and H
continues to resolve.

In the solid regime H is flat at 0.091-0.095 from T0 = 0.05 to 0.40 while
Lindemann nearly triples. H reads the icosahedral geometry; Lindemann reads
thermal agitation. They are independent axes.

### As a state variable, the triple-fraction wins

Spearman correlation with cumulative evaporated count, computed **within**
fixed temperature to remove the confound that both rise with T:

| signal | mean rho | T0=1.6 | T0=2.0 |
|---|---|---|---|
| **H_frac** (fraction of triples with H > 0.5) | **0.952** | **0.997** | **0.831** |
| H_mean | 0.866 | 0.980 | 0.819 |
| Lindemann | 0.784 | 0.886 | 0.716 |
| mean pairwise distance | 0.778 | 0.751 | 0.666 |
| max pairwise distance | 0.727 | 0.646 | 0.534 |

H_frac beats Lindemann and both distance baselines, and the margin widens in
the heavily evaporating regime where the others degrade. This is the first
setting in which the H construction outperforms the trivial alternatives, and
it requires the aggregated form — local H on every triple, globally pooled.

**Caveat, stated plainly:** the correlation is partly definitional. When an
atom departs, all 66 triples containing it acquire large distance variance, so
H_frac counts departures in a smooth continuous form rather than inferring
them independently. A continuous integrity readout is more useful than a step
function, but it is a measurement, not a prediction. The comparison against
Lindemann remains fair, since Lindemann is equally distance-derived.

### It does not forecast

Predicting *new* evaporations in a forward window of 8 or 20 time units,
split-half validated:

| signal | Spearman | R^2 (test) |
|---|---|---|
| H_mean | 0.034 | 0.017 |
| H_p99 | -0.044 | 0.040 |
| Lindemann | -0.049 | -0.008 |
| mean distance | -0.032 | 0.043 |

All near zero, for every signal including the baselines. Evaporation in a short
window is near-Poisson: almost always zero, occasionally one. There is little
variance to explain and nothing explains it.

The apparent H-versus-evaporation relationship in the temperature table above
is confounded by temperature driving both. Within a trajectory at fixed T, no
tested signal forecasts the next departure.

---

## What does not work

- **Curvature spike.** Frenet curvature of the escaping atom, |v × a|/|v|³.
  Median lead 9.42 but CV 1.007, with several leads under 0.25 — the final
  close encounter *is* the escape, so its curvature spike is coincident rather
  than predictive.
- **Curvature trend.** Log-slope of peak sharpness before escape: median
  +0.0007, positive in 56% of trajectories. A coin flip. Encounters do not get
  systematically harder in the run-up.
- **sigma or r_min normalisation of distances** (section 1).
- **Subtracting a thermal floor** (section 3).

---

## Scope of the claim

This is not a demonstration that PNEP generalises to molecular dynamics, and
it is not a claim that H is useful in this setting.

What it establishes: a geometric node scalar yields a **well-calibrated
dissociation hazard** in Lennard-Jones systems within a bounded regime
(tau ~ 1), the boundary is identifiable in advance, and the transfer requires
an equilibrium-referenced coordinate to avoid inverting on repulsive
configurations.

What it refutes: that **H specifically** is the right scalar. Mean pairwise
distance matches or beats it on every metric tested. H's advantage in the
gravitational case does not survive the move to a non-monotonic potential with
only three bodies.

Whether H earns its place in larger clusters -- where there is no single
obvious distance to monitor and the scalar is computed over many triples -- is
open and under test.

Single potential form, three bodies, one trimer geometry, small per-cell
samples. Generalisation to larger clusters, other potentials, solvated systems
or biomolecular force fields is untested.

---

Alika Parks · October 2026
