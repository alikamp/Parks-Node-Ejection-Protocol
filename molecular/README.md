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

## 7. Correlation structure of the H field (underpowered, not refuted)

If the H field over triples had a spatial correlation length, and if that
length grew on approach to a transition, that would connect the scalar to
standard critical phenomena rather than leaving it competing with simpler
baselines. `correlation_structure.py` tests it.

For all 40,755 pairs of the 286 triples: shared-atom count (0, 1, 2),
centroid separation, and the Pearson correlation of H_i(t) with H_j(t).

| T0 | Lindemann | C(overlap 0, near) | C(overlap 0, far) | xi |
|---|---|---|---|---|
| 0.20 | 0.032 | +0.013 | −0.006 | — |
| 0.40 | 0.049 | +0.014 | −0.010 | — |
| 0.60 | 0.116 | +0.030 | +0.012 | — |
| 0.75 | 0.248 | +0.006 | −0.062 | — |
| 0.90 | 0.330 | −0.001 | −0.070 | — |
| 1.10 | 0.321 | −0.016 | +0.028 | — |
| 1.40 | 0.439 | −0.034 | +0.380* | — |

No correlation length is extractable at any temperature: the profile is flat,
so there is no exponential to fit. Disjoint triples stay within +/-0.034
across the full range, including through melting at T0 = 0.60.

The one real correlation is between triples sharing **two** atoms
(C_near = +0.317): those share an edge, so one of the three distances
entering each variance is the same number. That is arithmetic, not physics.
Triples sharing one atom or none show nothing.

\* the +0.380 at T0 = 1.40 is an evaporation artifact — a departed atom drives
every triple containing it to H ~ 1 simultaneously.

**This result is underpowered, not refuting.** An LJ-13 cluster spans 2-3
sigma. A correlation length requires a system many xi across; here the entire
object is smaller than the smallest interesting xi, and every triple pair is
either touching or spanning the cluster, with no intermediate range over which
decay could appear.

Testing this properly requires a bulk system — several hundred atoms with
periodic boundaries — where C(r) can be measured over a real range of r and
xi(T) tracked on approach to the transition. Not attempted here.

---

## 8. Crystallisation: inconclusive by experimental design

`quench.py`. Melt to liquid at T = 0.75 (Lindemann 0.19, no evaporation), then
Langevin-quench to T = 0.06, 10 runs.

Prediction tested: melting is homogeneous and should move a global average
smoothly, whereas freezing is nucleation-driven and should leave the mean flat
while the H *distribution* goes bimodal as an ordered subpopulation separates.

Result: H reaches its solid value at the same instant the potential energy
says the cluster froze (median lag −0.2, range −0.8 to 0.0). Bimodality is
flat at 0.532–0.534 throughout and peaks at random times.

But the quench was too fast to test the hypothesis. Freezing completed at
t = 0.6–5.2 inside a 1600-unit run: the cluster collapsed into the nearest
minimum without nucleating, so there was no process with a beginning and
middle to measure a lead or lag against. The right experiment is isothermal
supercooling near T = 0.25–0.35, where nucleation is slow and stochastic.

Two things did come out of it:

- **H's solid value reproduces from both directions.** Heating gave
  0.091–0.095; cooling converged to 0.091 from above. The scalar reads a
  structural property, not a path-dependent artifact.
- **H flagged a trapped structure.** Run 8 never reached the frozen energy
  threshold, ending at U = −40.8 against −43.3 for the others — an amorphous
  local minimum, the known LJ-13 trapping behaviour. Its H was 0.113 against
  0.091 for the runs that found the ground state.

Note on methodology: the global Q6 bond-orientational parameter is useless as
ground truth here. Icosahedral symmetry causes near-complete cancellation when
q6m is averaged over all bonds, so an ordered icosahedron is indistinguishable
from a disordered cluster. Potential energy against the known LJ-13 global
minimum (−44.327) was used instead.

---

## 9. Spatial domains in the H field

Pairwise correlation (section 7) asks whether H_i(t) covaries with H_j(t) for
fixed pairs. A different question: at each instant, are the triples with
*similar* H spatially contiguous? `domains.py`.

Select the top and bottom quartile of triples by H, compare their mean
pairwise centroid separation to random subsets of equal size, report a
z-score. Negative = clustered.

| T0 | z (top quartile) | z (bottom quartile) | atom concentration | z (disjoint control) |
|---|---|---|---|---|
| 0.40 | −8.87 | +3.82 | 0.102 | −2.32 |
| 0.60 | −7.48 | +3.67 | 0.115 | −2.33 |
| 0.75 | −5.90 | +3.04 | 0.135 | −2.07 |
| 0.90 | −5.79 | +2.51 | 0.141 | −1.80 |
| 1.10 | −7.37 | +2.59 | 0.131 | −2.04 |
| 1.40 | −5.16 | +3.25 | 0.138 | −2.07 |

**The structure is real.** High-H triples cluster; low-H triples disperse.
Two controls support it: atom concentration is 0.10–0.14 against 0.077 for
even representation, so no single stray atom drives the selection; and the
disjoint-subset test — restricted to triples sharing no atoms with each other,
removing all shared-distance coupling — still gives z ≈ −2.

Note that the raw z of −5 to −9 overstates it: roughly two thirds of that is
shared-atom coupling. The honest figure is the disjoint control.

**But it is temperature-invariant.** z(disjoint) is flat from −1.80 to −2.33
across solid, melted and heavily evaporating states. Domain structure is a
static geometric property of the cluster, not a precursor. There is no early
warning in it.

This also explains the null result in section 7: time-covariance of fixed
pairs and instantaneous spatial contiguity of a quantile are different
measurements, and only the second one sees this.

---

## 10. Turning angle: real signal, inverted from the hypothesis

Section "what does not work" records a failed curvature test using the Frenet
curvature |v × a| / |v|³. That was the wrong quantity: dividing by speed cubed
means a fast atom taking a hard hit scores low while a slow drifting atom
scores high. The scale-free measurement is the angle between successive
velocity vectors. `turning_angle.py`.

Mean turning angle by lag (π/2 = 1.571 indicates full decorrelation):

| T0 | lag 1 | lag 5 | lag 10 | lag 25 |
|---|---|---|---|---|
| 0.40 | 0.213 | 1.019 | 1.669 | 1.545 |
| 0.75 | 0.195 | 0.928 | 1.524 | 1.536 |
| 1.10 | 0.172 | 0.819 | 1.354 | 1.409 |
| 1.40 | 0.156 | 0.737 | 1.225 | 1.310 |

Velocity memory persists to about lag 10 (0.2 time units) and is gone by lag
25. An earlier version of this test used lag 25 and therefore measured the
random-direction value at every temperature.

Correlation of turning angle with H, within fixed temperature:

| T0 | lag 5 |
|---|---|
| 0.40 | −0.113 |
| 0.75 | −0.259 |
| 1.10 | −0.466 |
| 1.40 | −0.463 |

**Negative, strengthening with temperature.** Higher H accompanies *smaller*
turning angles. The physical reading: an atom still bound is batted around and
turns constantly; an atom on its way out is moving ballistically and has
stopped turning. The kick that frees it is itself a hard turn, but it lasts an
instant — what persists, and what is detectable, is the straightness that
follows.

At ρ = −0.47 this is the second-strongest signal in this directory after
H_frac, and it is a velocity-space observable, independent of the
position-space quantities everything else here measures.

---

## 11. Band populations and transfer rates: state yes, precursor no

Triples binned into H bands [0, .05, .10, .20, .35, .55, 1.01]; populations
and band-to-band transfer rates tracked over time.
`band_populations.py`, `band_flux.py`. 30 seeds per temperature.

### Composition is a clean state readout

| T0 | band2 (.10-.20) | band3 (.20-.35) | band5 (.55+) | Lindemann | evaporated |
|---|---|---|---|---|---|
| 0.40 | 0.539 | 0.030 | 0.000 | 0.047 | 0.00 |
| 0.60 | 0.520 | 0.052 | 0.000 | 0.092 | 0.00 |
| 0.75 | 0.484 | 0.081 | 0.016 | 0.170 | 0.07 |
| 0.90 | 0.443 | 0.104 | 0.049 | 0.238 | 0.21 |
| 1.10 | 0.381 | 0.122 | 0.125 | 0.300 | 0.57 |

Monotonic transfer out of band 2 into bands 3 and 5. Band 3 nearly doubles
between T0 = 0.40 and 0.60, before Lindemann crosses the melting criterion and
before anything evaporates.

### Net flux is zero by detailed balance

| T0 | up flux | down flux | net | mobility |
|---|---|---|---|---|
| 0.60 | 0.0149 | 0.0149 | +0.00000 | 0.0298 |
| 0.75 | 0.0147 | 0.0147 | +0.00000 | 0.0293 |
| 0.90 | 0.0140 | 0.0140 | +0.00000 | 0.0280 |
| 1.10 | 0.0122 | 0.0122 | +0.00001 | 0.0243 |

Up and down flux cancel to five decimals at every temperature. A trajectory-
averaged net flux cannot be a precursor: it is structurally zero.

### Head to head as a precursor — Lindemann wins

Window before the first evaporation, against a matched window in runs at the
same T0 that never evaporate.

| T0 (evap/stable) | net flux | up flux | mobility | band3 | band5 | **Lindemann** |
|---|---|---|---|---|---|---|
| 0.75 (3/27) | 0.012 * | 0.021 * | 0.223 | 0.032 * | 0.167 | **0.000 *** |
| 0.90 (9/21) | 0.326 | 0.123 | 0.096 | 0.056 | 0.062 | **0.000 *** |
| 1.10 (20/10) | 0.447 | 0.000 *† | 0.000 *† | 0.086 | 0.001 * | **0.000 *** |

(p-values; † = significant but in the **opposite** direction to the hypothesis)

Net flux is significant only where n = 3. At the temperatures with real
statistics it fails. Up-flux and mobility reach p = 0.000 at T0 = 1.10 but
both are *lower* before evaporation, not higher: once a triple contains an
outbound atom it pins in the top band and stops crossing, so the composition
freezes as the cluster comes apart. A real mechanism, the reverse of the
hypothesised churn. Band 5 population is significant but near-definitional.

Lindemann is significant at every temperature at p = 0.000.

---

## 12. Per-trio centroid paths: straighter at higher H

Sections 9-11 treat triples as static objects. This treats each trio as a
moving unit: track its centroid, measure the turning angle of that centroid's
path, its straightness (net displacement / path length), and whether trios in
the same H band share a direction. `trio_paths.py`, 20 seeds per temperature.

| T0 = 0.90 | turn | straight | coherence | coherence (disjoint) |
|---|---|---|---|---|
| band 0.00 | 0.555 | 0.945 | 0.139 | 0.209 |
| band 0.05 | 0.554 | 0.946 | 0.177 | 0.240 |
| band 0.10 | 0.551 | 0.946 | 0.104 | 0.187 |
| band 0.20 | 0.555 | 0.945 | 0.277 | 0.272 |
| band 0.35 | 0.477 | 0.957 | 0.589 | 0.291 |
| band 0.55 | 0.384 | 0.971 | 0.686 | — |

### Coherence is mostly artifact

Raw coherence rises about tenfold from the low bands to the top (0.07 to
0.69). Under the disjoint control — restricted to trios sharing no atoms —
it rises about 1.7-fold (0.17 to 0.29). A real residual survives, consistently
at all four temperatures, but the dramatic alignment is a departing atom
dragging the 66 trios that contain it.

Band 0.55 cannot be disjoint-controlled at all: there are never three mutually
disjoint trios in it. That band *is* one atom's triples, which confirms the
mechanism directly rather than inferring it.

### Straightness is real, and is not a departing-atom effect

The clean case is T0 = 0.60, where **evaporation is exactly zero** across all
20 runs — no atom leaves, so no ballistic centroid can explain anything:

| T0 = 0.60 (no evaporation) | turn | straight | population |
|---|---|---|---|
| band 0.00 | 0.598 | 0.938 | 0.253 |
| band 0.10 | 0.587 | 0.940 | 0.518 |
| band 0.35 | 0.518 | 0.952 | 0.022 |

Trios in band 0.35 — 2.2% of all triples — trace measurably straighter paths
and turn less than the bulk, inside a cluster that is losing nothing. The
trend is monotonic across bands and holds at every temperature tested
(T0 = 0.60, 0.75, 0.90, 1.10).

This is consistent with the turning-angle result of section 10, measured at
the level of trios rather than atoms: coherent, low-curvature motion
accompanies high H, and jitter accompanies low H.

Whether straightness *precedes* anything was not tested here; on the evidence
of sections 7, 9, 10 and 11 the prior should be that it does not.

---

## 13. The H field is a four-state kinetic system

Sections 11 and 12 binned H on a grid chosen for convenience. This asks
whether the bins are real. `kinetics.py`, 12 seeds per temperature.

The H distribution over all 286 triples is **multimodal with
temperature-invariant mode positions**: peaks near H = 0.038 and H = 0.113 sit
at the same place at T0 = 0.40, 0.75 and 1.10, and a third near 0.96 appears
only once evaporation starts. Heating changes the occupancies, not the
locations. That is what licenses a discrete-state treatment instead of a
continuum, so the states are cut at the density minima rather than at round
numbers:

| state | range | reading |
|---|---|---|
| A | H < 0.07 | compact |
| B | 0.07 – 0.22 | bulk |
| C | 0.22 – 0.60 | stretched |
| D | H ≥ 0.60 | detached |

### Rate constants

| transition | 0.50 | 0.65 | 0.80 | 0.95 | 1.10 |
|---|---|---|---|---|---|
| A→B | 0.387 | 0.472 | 0.521 | 0.587 | 0.592 |
| B→A | 0.144 | 0.190 | 0.222 | 0.267 | 0.281 |
| B→C | 0.073 | 0.124 | 0.163 | 0.225 | 0.242 |
| C→B | 4.966 | 2.936 | 2.056 | 1.535 | 1.414 |
| C→D | 0.000 | 0.000 | 0.002 | 0.002 | 0.002 |
| D→C | 0.000 | 0.000 | 0.001 | 0.002 | 0.001 |

### Arrhenius

| transition | Ea (ε) | ln A | R² | n |
|---|---|---|---|---|
| A→B | 0.406 | −0.134 | 0.988 | 5 |
| B→A | 0.624 | −0.695 | 0.993 | 5 |
| B→C | **1.124** | −0.368 | 0.993 | 5 |
| C→B | **−1.188** | −0.769 | 0.996 | 5 |
| C→D | 0.472 | −5.723 | **0.287** | 3 |

Four of the five transitions are Arrhenius at R² ≥ 0.988, extracted from a
purely geometric scalar — H sees only the variance of three distances, never
an energy — with B→C the largest at 1.124 ε, in the region of one pair
contact.

> **Qualified by section 14.** The state process is *not* Markovian, so these
> are lag-dependent empirical crossing frequencies, not Markov rates. The
> words "barrier", "rate-limiting" and "transition-state" are withdrawn
> pending the stride-dependence test; what the fits establish is a clean
> exponential temperature dependence, which is a weaker claim.

C→B has a *negative* activation energy and the rate genuinely falls with
heating, 4.97 down to 1.41. That is the expected signature of a barrierless
relaxation: a stretched trio collapsing back into the bulk basin has nothing
to climb, and higher thermal velocity makes recapture less likely.

### Detailed balance locates the one irreversible channel

Fractional flux imbalance |n_ij − n_ji| / n_ij:

| T0 | A↔B | B↔C | C↔D |
|---|---|---|---|
| 0.50 | 0.01% | 0.00% | — |
| 0.65 | 0.02% | 0.00% | — |
| 0.80 | 0.02% | 0.02% | **81.5%** |
| 0.95 | 0.02% | 0.05% | **69.7%** |
| 1.10 | 0.01% | 0.08% | **59.1%** |

The A↔B↔C ladder is in detailed balance to within 0.08%. C↔D is not, by two
to three orders of magnitude more.

This reframes section 11. The net flux there was structurally zero not
because the H field carries no kinetic information, but because the part of
it that was measured is **equilibrated**. The system is a reversible
three-state ladder feeding a single absorbing sink, and all of the ejection
physics lives in the one step that breaks detailed balance.

### Occupancy and dwell

| T0 | A occ | B occ | C occ | D occ | A dwl | B dwl | C dwl | D dwl |
|---|---|---|---|---|---|---|---|---|
| 0.50 | 0.2690 | 0.7204 | 0.0106 | 0.0000 | 2.59 | 4.60 | 0.20 | — |
| 0.65 | 0.2783 | 0.6926 | 0.0292 | 0.0000 | 2.12 | 3.19 | 0.34 | — |
| 0.80 | 0.2735 | 0.6418 | 0.0510 | 0.0337 | 1.92 | 2.59 | 0.49 | 1849 |
| 0.95 | 0.2764 | 0.6074 | 0.0892 | 0.0270 | 1.70 | 2.03 | 0.65 | 580 |
| 1.10 | 0.2525 | 0.5318 | 0.0909 | 0.1248 | 1.69 | 1.91 | 0.71 | 1447 |

- **D is absorbing.** Dwell of order 10³ time units against ~2 for every other
  state. Detachment does not reverse.
- **C is short-lived relative to A and B.** Dwell 0.20–0.71, and C→B runs
  6–68× faster than B→C. Trios pass through it rather than settling in it.
  (Stated as a dwell measurement; "transition state" in the kinetic sense
  would require the Markov property, which section 14 shows does not hold.)
- **A's occupancy is flat** across a 2.2× temperature range (0.253–0.278,
  regression p = 0.357) while B drains monotonically from 0.720 to 0.532
  (p = 0.002). The compact state is thermally insulated in a way the bulk
  state is not.
- C's dwell time *rises* with temperature, 0.20 to 0.71. The transition state
  gets stickier as the cluster is heated. Unexplained.

### What is not determined

**The rate that matters is the one that is not pinned.** Ejection flux
factorises as occ(C) × k_CD. occ(C) is a pure geometric quantity computable
from H and is well determined (0.011 → 0.091). k_CD is not: only 243, 379 and
521 events across the three temperatures that produce any, running
non-monotonically (0.00193, 0.00172, 0.00232) for an Arrhenius R² of 0.287.

Whether the ejection rate is predictable from geometry alone turns on whether
k_CD is temperature-flat, and the present statistics cannot say. `kcd.py`
addresses this directly with 40 seeds over a hotter range.

Two further cautions. The quasi-equilibrium identity flux = occ(C) × k_CD is
**algebraically circular** as written — it expands to the observed flux — so
it is a bookkeeping decomposition and not a test of anything. And an earlier
version of this script computed dwell times against a sample count rather
than an exit count; the published table above is the corrected one, and the
error had inverted the apparent temperature trend of state A.

A third, found later: **the Poisson errors on k_CD understate the real
uncertainty.** A 40-seed run at T0 = 0.95 gives k_CD = 0.00127 ± 0.00004
against 0.00172 ± 0.00009 from the 12-seed run above — 4.5σ apart on
different seeds. Run-to-run heterogeneity, not counting noise, dominates the
error budget, because whether a given trajectory evaporates at all is close
to binary. Any future quote of k_CD needs seed-level scatter, not √n.

---

## 14. The four-state model is not Markovian

Section 13 reports rate constants and activation energies. Those quantities
only carry their usual meaning if the state process is Markovian on this
state space, which section 13 assumed without testing. `markov.py` tests it
two standard ways, 8 seeds per temperature.

### Implied timescales never plateau

t_i(τ) = −τ / ln λ_i(τ) from the eigenvalues of the lag-τ transition matrix.
For a Markov process these are independent of τ.

| T0 = 0.95, lag (t.u.) | 0.02 | 0.04 | 0.08 | 0.16 | 0.32 | 0.64 | 1.28 |
|---|---|---|---|---|---|---|---|
| t2 | 1161 | 1178 | 1211 | 1320 | 1565 | 1986 | 2664 |
| t3 | 1.261 | 1.260 | 1.275 | 1.359 | 1.639 | 2.238 | 3.138 |
| t4 | 0.600 | 0.594 | 0.591 | 0.613 | 0.723 | 1.035 | 1.530 |

Every column drifts monotonically upward and is still climbing at the longest
lag tested. t3 and t4 grow by about 2.5×.

### Chapman-Kolmogorov fails

Measured T(nτ) against predicted T(τ)ⁿ, base lag 0.16 t.u.:

| n | lag (t.u.) | T0 = 0.65 | T0 = 0.95 | T0 = 1.10 |
|---|---|---|---|---|
| 2 | 0.32 | 0.070 | 0.042 | 0.041 |
| 4 | 0.64 | 0.264 | 0.158 | 0.153 |
| 8 | 1.28 | 0.353 | 0.267 | 0.252 |

Counting noise here is of order 0.001. Errors of 0.25–0.35 on a probability
matrix are a failure, not a near miss.

### The stationary distribution is wrong

Leading left eigenvector against directly measured occupancy:

| T0 | π (eigenvector) | measured |
|---|---|---|
| 0.95 | 0.228 / 0.479 / 0.092 / 0.201 | 0.255 / 0.534 / 0.102 / 0.108 |
| 1.10 | 0.128 / 0.259 / 0.048 / 0.566 | 0.224 / 0.455 / 0.083 / 0.237 |

D is off by roughly 2×. At T0 = 0.65 the computation returns numerical
garbage in t2 (values of order 10¹⁷): D is almost never visited there, the
matrix goes near-singular, and the solver collapses onto the absorbing state.

### Diagnosis

**D being absorbing is the primary cause.** t2 ≈ 10³ is not a relaxation
timescale, it is the drainage time into D, and an absorbing state makes the
stationary distribution meaningless. But t3 and t4 drift as well, so the
A↔B↔C ladder is non-Markovian on its own account. Four states cut on H are
too coarse: the H coordinate does not separate the system's true metastable
states, and there is structure inside the bands.

### What this does and does not retract

Unaffected, because they are direct counts rather than model inferences:
multimodality with temperature-invariant mode positions; detailed balance at
0.01–0.08% on the ladder and 59–82% breakage at C↔D; D's absorbing dwell;
C's short dwell; all occupancies.

Retracted: "barrier", "rate-limiting step", "transition state" in the
kinetic sense, and the reading of 1.124 ε as a pair-contact binding energy.

Still standing, as a weaker claim: the empirical crossing frequencies have a
clean exponential temperature dependence at R² ≥ 0.988. Whether the fitted
Ea is a physical quantity or an artifact of the sampling stride is decided by
`ea_lag.py`, which refits Ea separately at each lag. Flat Ea(τ) means the
activation energy survives the loss of the Markov property; drifting Ea(τ)
means 1.124 ε was a consequence of choosing SAMP = 5.

---

## 15. The activation energies do not survive: withdrawn

`ea_lag.py` refits Arrhenius separately at each lag, 5 temperatures × 8 seeds.

### Ea depends on the sampling stride

| transition | τ = 0.02 | 0.08 | 0.32 | 1.28 | drift |
|---|---|---|---|---|---|
| A→B | 0.187 | 0.197 | 0.279 | 0.511 | **173%** |
| B→A | 0.386 | 0.396 | 0.478 | 0.710 | **84%** |
| B→C | 0.514 | 0.524 | 0.617 | 0.814 | **58%** |
| C→B | −0.420 | −0.410 | −0.317 | −0.122 | **71%** |
| C→D | 4.323 | 4.504 | 5.495 | 6.759 | **56%** |

No transition is stride-independent. The activation energy is a function of
how often the trajectory was sampled, not a property of the system.

### Ea also moves under reseeding

At the *same* lag of 1 sample, this run gives B→C Ea = 0.514. Section 13 gives
**1.124** for the same transition at the same lag. The only difference is the
seed set (3000–3007, n = 8 against 0–11, n = 12). A factor of 2.2 from
reseeding alone.

So the R² of 0.988–0.996 quoted in section 13 was measuring the wrong thing:
a tight line through five points that are themselves unstable under
resampling. Excellent fit, badly determined parameters.

**The activation energies of section 13 are withdrawn** — not reinterpreted,
withdrawn. They are not reproducible across strides or across seeds.

### The one ordering that does survive

C→D has an activation energy several times larger than any ladder transition:
4.3–6.8 here, and ≈ 3.3–3.9 from the independent 40-seed `kcd.py` run
(k_CD = 0.00127, 0.00203, 0.00312 at T0 = 0.95, 1.10, 1.25, monotonic and
well outside its error bars). Different estimator, different seeds, same
conclusion — detachment costs several times what a band crossing costs. The
ordering is robust even though no magnitude is.

### The methodological finding, which matters more

Three independent symptoms of the same problem:

| quantity | run 1 | run 2 | cause of difference |
|---|---|---|---|
| Ea(B→C) | 1.124 | 0.514 | seed set only |
| k_CD(T0 = 0.95) | 0.00172 ± 0.00009 | 0.00127 ± 0.00004 | seed set only (4.5σ) |
| CK error, n = 8, T0 = 0.65 | 0.353 | 0.620 | seed set only |

Every quoted error bar in sections 13–14 came from within-run sample counts
(123M samples) when the real variance is between-run (6–12 trajectories).
**At this seed count essentially nothing in the molecular kinetics is
reproducible.** Whether a given trajectory evaporates at all is close to
binary, so run-to-run heterogeneity dominates and √n is the wrong
denominator. Any future number from this line of work needs seed-block
scatter, not sample-count Poisson errors.

---

## 16. A second coordinate does not recover the Markov property

If H is too coarse, the natural fix is H plus a second coordinate. The
best-motivated candidate came from section 1: H is **sign-blind**, since
variance is symmetric, so a trio at H = 0.15 might be compressed or stretched
and H cannot tell. Two populations sharing one H value would produce exactly
the drifting timescales observed.

Second coordinate: u = mean over the trio's three pairs of (r − r_min)/r_min.
`usplit.py`, 6 seeds per temperature, four models at matched state count.

### The mechanism is refuted before the test runs

| T0 = 0.65, band | population | median u | fraction compressed (u < 0) |
|---|---|---|---|
| A | 0.2877 | +0.1260 | 0.116 |
| B | 0.6411 | +0.4514 | 0.000 |
| C | 0.0712 | +0.7031 | 0.000 |
| D | 0.0000 | +1.5693 | 0.000 |

**There is no compressed population.** Essentially every triple is stretched
relative to r_min, and median u rises monotonically across bands. H is not
hiding a compressed/stretched degeneracy because one of the two does not
exist in this system.

### The test runs anyway, and u loses

Split at the band median of u rather than at zero, so the test still asks
whether u-magnitude adds information. Controls: M3 splits each band at its
median H (does finer resolution of the same coordinate do as well?), M4 splits
at a shuffled u (does anything beyond state count operate at all?).

| T0 = 0.95 | CK n=2 | CK n=4 | CK n=8 |
|---|---|---|---|
| M1 H bands (4) | 0.0360 | 0.1344 | 0.2300 |
| M2 H × u (8) | 0.0498 | 0.1540 | 0.2387 |
| M3 H × H (8) | **0.0184** | **0.0723** | **0.1440** |
| M4 H × shuffled (8) | 0.0528 | 0.1149 | 0.1856 |

| T0 = 0.65 | CK n=2 | CK n=4 | CK n=8 |
|---|---|---|---|
| M1 H bands (4) | 0.0854 | 0.3160 | 0.6197 |
| M2 H × u (8) | 0.0882 | 0.3440 | 0.6067 |
| M3 H × H (8) | 0.0705 | 0.2735 | 0.4600 |
| M4 H × shuffled (8) | **0.0512** | **0.1863** | **0.4053** |

u is the worst of the four at T0 = 0.95 — worse than not splitting at all.
At T0 = 0.65 the *shuffled* control wins, so there a random cut outperforms
both informed coordinates.

**Finer resolution of H (M3) beats a second coordinate (M2) consistently.**
The bands are too coarse; the scalar is not missing a dimension. "H plus a
second coordinate" is not the shape of the fix, at least not with the
best-motivated candidate available.

### A trap worth recording

M2 at T0 = 0.95 produces a beautifully flat implied timescale — 3599.9 at
lag 1 against 3589.5 at lag 64, drift 0.3% — while posting the *worst* CK
error of the four models. The D-band split created a near-absorbing pair with
λ ≈ 1. A flat implied timescale is not evidence of Markovianity on its own;
run with CK alongside it or the artifact reads as a success.

---

## 17. The bands are the icosahedron's geometry, not dynamics

Sections 13–16 treated the four H bands as a kinetic system and could not
make one fit. This is the null hypothesis that should have been tested
first: what does H look like in the **perfect icosahedron at T = 0**, with no
dynamics at all? `static_baseline.py`, no simulation required.

### The band structure is already there at absolute zero

286 triples take only **six distinct H values**:

| H | count | fraction | what it is |
|---|---|---|---|
| 0.0000 | 120 | 0.420 | |
| 0.1595 | 60 | 0.210 | |
| 0.1020 | 40 | 0.140 | equilateral vertex triples (20 faces + 20 larger) |
| 0.0007 | 30 | 0.105 | centre + edge (30 edges) |
| 0.1168 | 30 | 0.105 | |
| 0.2119 | 6 | 0.021 | centre + antipodal pair (6 antipodal pairs) |

### It reproduces the MD occupancies

| | A | B | C | D |
|---|---|---|---|---|
| **static, T = 0** | **0.2448** | **0.7552** | **0.0000** | **0.0000** |
| MD, T0 = 0.50 | 0.2690 | 0.7204 | 0.0106 | 0.0000 |
| MD, T0 = 1.10 | 0.2525 | 0.5318 | 0.0909 | 0.1248 |

A structure at absolute zero matches the T0 = 0.50 occupancies to within
0.035. The combinatorics close exactly, so this is not coincidence: contact
triples number **50 = 20 faces + 30 centre-plus-edge**, and the 6 triples at
the top H value are **centre plus the 6 antipodal vertex pairs**. Every count
is a structural quantity of the icosahedron.

### This deflates the three most interesting results above

- **"Mode positions do not move with temperature"** (section 13) — they are
  fixed by which shells the three atoms come from. Heating only broadens them.
- **"State A's occupancy is flat across a 2.2× temperature range"** (section
  13) — it is 0.2448 from pure geometry, and moves by at most 0.034.
- **"The model is non-Markovian"** (section 14) — the states were never
  dynamical states. They are **triple types**. A centre-plus-antipodal-pair
  triple stays that type, so there was no Markov process to find. This is the
  mechanism behind the failure, not just a restatement of it.

### What survives, and it is sharper than what it replaces

The static maximum is H = 0.2119 and the C band edge is 0.22, so **bands C
and D contain zero static triples by construction.** Their occupancy is
thermal excitation with the geometric baseline already subtracted:

| T0 | ΔA | ΔB | ΔC | ΔD | **C + D** |
|---|---|---|---|---|---|
| 0.50 | +0.024 | −0.035 | +0.011 | +0.000 | **0.0106** |
| 0.65 | +0.034 | −0.063 | +0.029 | +0.000 | **0.0292** |
| 0.80 | +0.029 | −0.113 | +0.051 | +0.034 | **0.0847** |
| 0.95 | +0.032 | −0.148 | +0.089 | +0.027 | **0.1162** |
| 1.10 | +0.008 | −0.223 | +0.091 | +0.125 | **0.2157** |

C + D rises monotonically by a factor of 20 and is zero-baselined by
construction. **That is the state variable** — not the four-state occupancy
vector, and not H itself. It is close to the H_frac measure of section 6,
which was the one genuine win of the molecular work, and section 17 explains
why that one worked while the kinetic elaboration did not: H_frac counts
thermally excited triples, and the geometric baseline cancels out of it.

### Methodological note

Every band-based result in sections 11–16 should have been checked against
the static structure before any dynamics were run. The check costs no
simulation and takes one second.

---

## Summary of the molecular work

Three statements, all supported:

- **As a state readout, H_frac beats the Lindemann parameter** — within-
  temperature correlation with cluster disintegration of 0.952 against 0.784,
  with the margin widening where Lindemann saturates (section 6).
- **As a precursor, Lindemann beats every H-derived quantity tested** — band
  populations, band transfer rates, spatial domain structure, turning angle,
  and H itself (sections 7, 9, 10, 11).
- **The H field's discrete structure is the icosahedron's geometry** — the
  perfect T = 0 cluster already reproduces the band occupancies to within
  0.035, so the bands are classes of triple geometry rather than metastable
  states (section 17). What is left after subtracting that baseline, the
  C + D occupancy, rises 20-fold with temperature and is the real state
  variable.

The H construction describes the current state of structural integrity better
than the standard metric. It does not predict failure better. Those are
different claims and only the first one is supported here.

The third statement is deliberately narrow, and the kinetic elaboration built
on it did not hold. The four-state process is not Markovian (section 14); the
Arrhenius activation energies are withdrawn, reproducible across neither
sampling strides nor seeds (section 15); a second coordinate does not fix it,
and finer resolution of H outperforms one (section 16). Section 17 gives the
reason all of that failed: the states were triple geometries, not basins, so
there was no kinetics there to find.

### The limiting problem is statistical, not conceptual

Sections 13–16 quoted error bars from within-run sample counts — 123 million
samples — when the real variance is between-run, across 6 to 12 trajectories.
Reseeding alone moves Ea by 2.2×, k_CD by 4.5σ, and the CK error by 75%.
Whether a trajectory evaporates at all is close to binary, so run-to-run
heterogeneity dominates and √n is the wrong denominator.

**At this seed count the molecular kinetics does not reproduce.** That is the
first thing to fix before any of it is worth quoting, and it applies
retroactively to every number in sections 13–16 that is not a direct count.

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
