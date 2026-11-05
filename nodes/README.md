# Node Partners — which pair owns each mirror symmetry node

A mirror symmetry node fires at a **local minimum in any one pairwise
distance**. That means every node has an owner: the pair whose distance
bottomed out. The detector has always recorded this (`node_body` in
`Parks_NodalBreathing_FrequencyAnalysis.py`) and nothing in the framework
had ever looked at it.

It turns out to carry a precursor signal, a regime classifier, an identity
channel, and an analytic way to predict the first node from initial
conditions without integrating anything.

---

## The node condition is analytic

A local minimum of |r_ij| means

```
d|r_ij|²/dt = 2 r_ij · v_ij = 0,    with d²/dt² > 0
```

so "which pair fires next and when" is root-finding on three scalars of the
current state, not a search through a trajectory. The pair's radial velocity
`ṙ_ij = (r_ij · v_ij)/|r_ij|` is negative while closing.

This is the condition `pnep_v15_4_full_engine.py` already uses, restricted to
the outer body against the inner binary's centre of mass. Generalising it to
all three pairs reconciles v15.4 with the README's local-minimum definition —
they are the same criterion at different scope.

---

## 1. Partner rotation is strict round-robin

`node_partners.py`

| system | switch rate | max run length | ownership 01/12/02 |
|---|---|---|---|
| stable control | 1.0000 | 1 | 33.9 / 33.0 / 33.1 |
| slow decay | 0.9985 | 2 | 34.5 / 32.6 / 33.0 |
| slow decay | 1.0000 | 1 | 35.1 / 32.5 / 32.5 |

Equal thirds, never the same pair twice, holding from the first node through
ejection. The cause is kinematic: each of the three distances reaches a
minimum once per inner-binary period — d₀₁ at periapsis, d₀₂ and d₁₂ as
bodies 0 and 1 swing toward body 2 half a period apart. Three evenly spaced
minima per inner orbit gives observed Δt ≈ 1.55 against an inner period of
4.44.

The rotation itself therefore carries no information. **Its defects do.**

---

## 2. Round-robin defects: the strongest signal here

A **defect** is the same pair owning two consecutive nodes. 30 systems
spanning a_out 3–8, two mass ratios, two eccentricities, 4000t each.

| | defects | nodes | rate per 1000 | systems with ≥1 |
|---|---|---|---|---|
| ejectors | 37 | 10,878 | **3.401** | 11 / 12 |
| censored | 2 | 45,971 | **0.044** | 2 / 18 |

Fisher exact **p = 1.5 × 10⁻⁵**, odds ratio 88, rate ratio **78×**. If
censored systems defected at the ejector rate you would expect 156 defects in
their 45,971 nodes; there are 2 (Poisson p = 1.5 × 10⁻⁶⁴).

Survivors ran 45,971 nodes and broke the pattern twice. The round-robin is
not approximately strict in a stable hierarchy, it is near-inviolable — which
is what makes a single violation informative.

### It is a precursor, not a damage indicator

`defect_lead.py`

| lead fraction (fraction of lifetime remaining at first defect) |
|---|
| 0.954 · 0.900 · 0.862 · 0.835 · 0.819 · **0.788** · 0.769 · 0.745 · 0.433 · 0.125 · 0.077 |

Median **0.788** — the first defect fires with 79% of the system's life still
to run. Eight of eleven exceed 0.74. Median absolute lead 308.6t, range
73–2060t, against the Regime 2 slope forecast's quoted 412.6t.

**Lead fraction is independent of lifetime** (Spearman −0.20 vs t_eject,
p = 0.56). Across systems living 98t to 3504t the warning arrives at a
roughly constant fraction of the way through, which makes it scale-free in a
way an absolute lead time is not.

Defects then **accelerate** toward failure — pooled quartile counts Q1 9,
Q2 7, Q3 5, Q4 16 (χ² vs uniform p = 0.059). Not confined to the endgame, but
denser at the end. First defect warns; defect rate estimates how soon
(ρ = −0.55 vs t_eject, p = 0.063).

The signal is **binary and parameter-free**: no slope to fit, no window length,
no threshold. Did a pair own two consecutive nodes or not.

Two censored systems defected. Both are the most marginal high-eccentricity
survivors, and one fired at t = 245.8 and still had not ejected at t = 4000 —
consistent with being a real ejector past the observation window rather than
a false positive. One ejector (a_out 3.5, m 0.10, e 0.3, ejecting at t = 483)
shows no defect at all.

---

## 3. Per-partner jitter, and a threshold correction

Nodal Breathing pools the interval between *consecutive* nodes. But
consecutive nodes belong to different pairs, so the round-robin injects
structural variance that has nothing to do with resonance. Measuring the
interval between successive minima of the **same pair** removes it.

All signals measured in the **first 150 nodes only**, outcome measured
afterwards, systems ejecting inside the window excluded — so these are
forecasts, not hindsight fits. `headtohead.py`, 10 ejectors / 18 censored.

| signal | AUC | p | ρ vs t_eject | p |
|---|---|---|---|---|
| **per-partner jitter** | **0.972** | 0.0001 | **−0.721** | **0.019** |
| per-partner jitter (max) | 0.967 | 0.0001 | −0.600 | 0.067 |
| pooled jitter | 0.956 | 0.0001 | −0.382 | 0.276 |
| dH/dn (existing) | 0.778 | 0.0176 | −0.406 | 0.244 |
| H level | 0.522 | 0.8667 | −0.467 | 0.174 |

Per-partner jitter is the best separator tested and the only signal whose
rank correlation with time-to-ejection reaches significance.

**The 0.02 threshold in the main README is wrong by a factor of ~30.** Pooled
jitter does discriminate (AUC 0.956) but at a threshold near **0.60**:
ejectors occupy 0.5705–0.6443 and censored 0.5613–0.6006. At 0.02 the alert
fires at node 0 in every system tested, including a stable control that never
loses an atom — the structural resting value of a hierarchical triple is
σ/mean ≈ 0.52–0.57, not 0.02. The v15.1 reading of 0.611 is baseline, not an
anomaly.

**H level alone does not discriminate** (AUC 0.522, p = 0.87). Only the slope
does, which is exactly how Regime 2 is built. And dH/dn came out at +0.00150
mean for ejectors against the published 0.0012 — the Regime 2 slope
reproduces on an independent grid.

One honest limitation: per-partner jitter tracks a_out at ρ = −0.777, so it is
substantially a readout of the semi-major axis ratio. dH/dn tracks it at only
−0.204 — weaker separator, more independent information. Both matter, for
different claims. And both miss long-horizon ejectors: the systems ejecting at
t = 1918 and t = 3504 show negative dH/dn and low per-partner jitter in the
early window.

---

## 4. Predicting the first node from initial conditions

`firstnode2.py`, `firstnode3.py` — 192 configurations at randomised phase.

| predictor | pair correct | timing |
|---|---|---|
| **NEWTON** — one step of r·v with the full 3-body acceleration | **98.5%** | median 15.3% err, ρ = +0.955 |
| RV-SIGN — least-negative closing pair, no dynamics | 85.4% | — |
| PAIR2B — integrate the isolated pair | 40.6% | median 419% err |
| KEPLER — isolated two-body time to periapsis | 36.5% | median 422% err |

Chance on pair identity is 33.3%.

**The perturbation decides which pair fires.** Newton, which carries the full
three-body acceleration, gets it ~98% of the time across two seeds (97.8% and
98.5%). The two-body predictors, which ignore the third body, barely beat
chance.

**Kepler timing is bimodal and unusable without a selector.** Percentiles on
the Newton-chosen pair: p10 = 0.18%, p25 = 0.6%, **p50 = 422%**, p90 = 6205%.
Near-exact about 28% of the time, orders of magnitude wrong otherwise, and
nothing in the Kepler calculation distinguishes the cases.

So Newton alone is the predictor: ~98% on identity, median 15% on timing,
ρ = 0.955 for ranking. 61 of 192 arbitrary states yield no Newton prediction
(no pair closing with positive d(r·v)/dt) and need a different branch.

**The first node does not fix the rotation order.** Order (0,1,2) dominates at
55–60% regardless of which pair starts, and 57 of 192 show no clean 3-cycle.
The first node tells you where the sequence begins, not how it proceeds.

---

## 5. Flash regime: ownership collapses

`flash_partners.py`, DOP853 at rtol 1e-12 (the fixed-step leapfrog cannot
handle Burrau's close approaches). Ejection thresholds R > 12, 20, 50 as in
`verify_burrau.py`.

Ownership spread, most- minus least-owning pair, in percentage points:

```
hierarchical:  0.87 · 1.89 · 2.63
flash:         5.2 · 5.3 · 12.5 · 33.7 (Burrau) · 45.5 · 48.4 · 50.1
```

No overlap. Scramble #2 ran 49.3 / 0.3 / 50.4 with pair (1,2) effectively
starved. Hierarchical systems never leave 33 ± 3. **Ownership uniformity is a
regime classifier** reading geometry that the H value does not see. Pooled
jitter separates regimes too — flash 0.98–1.54 against 0.52–0.57 hierarchical.

Burrau validates against the repo: R > 12 gives t = 64.76 ejecting body 0, the
lightest, at 0.00000% energy drift.

### The ejecting body is in the least-owned pair, 6 of 6

| system | ownership 01/12/02 | starved | ejects |
|---|---|---|---|
| Burrau | 42.7 / 34.8 / 22.5 | 02 | body 0 |
| #1 | 35.3 / 29.5 / 35.3 | 12 | body 1 |
| #2 | 46.9 / 1.2 / 51.9 | 12 | body 1 |
| #3 | 48.2 / 47.4 / 4.4 | 02 | body 2 |
| #4 | 36.4 / 26.8 / 36.8 | 12 | body 1 |
| #5 | 10.3 / 42.7 / 47.0 | 01 | body 0 |

Narrows identity from three candidates to two, free from data the detector
already records. Random would be 2/3 per trial, so 6/6 sits at p = 0.088 —
suggestive, and a fourth independent channel alongside the three existing
identity signals rather than an improvement on them.

### No internodal drift before flash ejection

Δt_node trend slope −0.000022 (p = 0.82) and −0.000053 (p = 0.16), first/second
half means 0.1835→0.1798 and 0.1310→0.1209. The drift-type precursor that
works in slow decay does not transfer to flash.

### Regime 3 lead time, and a definition problem

Sustained-spike to ejection: 8.52t (Burrau), 10.79, 13.78, 13.96, 23.48,
27.61 — mean ~16t against the documented 2–10t. Every one exceeds the top of
the stated range.

But the repo's sustained-spike definition is *the first node where H crosses
0.5 and stays above for the rest of the integration*, which requires knowing
the future. Those lead times are not achievable live. The real-time version is
the first bare crossing, and those fire at t = 1.85, 3.01, 4.01, 4.68, 5.60,
10.90 — far too early, because H crosses 0.5 and falls back repeatedly inside
the scramble zone. The usable signal has poor precision; the precise signal
needs look-ahead.

---

## 6. Phase sensitivity of the published validation

`phase_revalidate.py`, `phase_revalidate2.py`

`make_ic` places all three bodies collinear on the x-axis with velocities in
the y–z plane. So **r·v = 0 exactly for all three pairs at t = 0** — verified
at 0.000e+00 across trials. All three pairs sit simultaneously at a mirror
symmetry node, which is a measure-zero configuration. Seven quantities are
randomised (three masses, a_in, e_in, period ratio, e_out, inclination, its
sign); all three orbital phases are pinned at apsis.

Test: build each system as published, advance it by a uniform random fraction
of its outer period using the module's own leapfrog, hand that state to an
**unmodified** `run_trial`. Same systems, same classification code, generic
phases, same seed.

Against `verify_revalidation.main(n_trials=100, seed=42, gt_max_t=3000)`,
which is the pipeline behind the published figure:

| | aligned (as published) | phase-randomised |
|---|---|---|
| n | 87 | 86 |
| accuracy | **98.9%** | **93.0%** |
| ejections caught | 36/36 = 100.0% | 33/36 = 91.7% |
| stable correctly called | 50/51 = 98.0% | 47/50 = 94.0% |
| EJECTED called stable | 0 | 3 |
| median energy drift | 0.00663% | 0.00607% |

The published numbers reproduce exactly. Under generic phase the method gives
93.0% with 91.7% ejection recall. Both runs contain exactly 36 ejectors, as
they should — advancing the phase does not change whether a system ejects.

So the headline accuracy is robust in kind, losing 5.9 points. The claim that
does not survive is **"zero ejecting systems misclassified as stable"** — that
is a property of the aligned IC family. On generic phases 3 of 36 ejections
are called stable.

The mechanism is visible in the H columns of the `run_batch` version of the
same test: H_stable barely moves (0.805 → 0.787) while H_unstable nearly
doubles (0.121 → 0.211), compressing the separation 16%. Starting every pair
at apsis makes an unstable configuration's geometry unusually clean at early
nodes.

**93.0% with 91.7% recall is the number to quote for arbitrary initial
conditions.** 98.9% is correct for the published IC family and should be
labelled as such.

---

## What does not work

- **Rotation chirality.** Whether the round-robin runs 01→12→02 or the reverse.
  Ejectors span −0.713 to +0.478 on both signs; Mann-Whitney p = 0.88 across
  30 systems. It does not even track the period ratio strongly enough to be
  called kinematic (ρ = +0.347, p = 0.06). Noise.
- **Newton + Kepler hybrid.** Newton picks the pair, Kepler times it. An
  earlier evaluation gave 0.28% median timing error, but that was computed
  only on the subset where Kepler's own minimum already agreed with Newton's
  pick — i.e. where Kepler was right. Scored on all cases it is 422% off.
  Selection artifact.
- **Coincident-firing artifact.** The detector appends one node per qualifying
  pair at the same timestep, which should put Δt = 0 entries into the series.
  Real but negligible: 0.16% of records at worst, 4 zero-Δt entries in 5146,
  mean interval shifts 0.08%, jitter 0.5667 vs 0.5658.
- **Confidence tiers by closing-pair count.** 1 / 2 / 3 closing pairs gave
  100 / 98.4 / 86.7% on one seed and 100 / 96.8 / 100% on another. At n ≈ 15
  per tier that is noise.

## Notes on repo consistency

- `verify_burrau.py` fires one node per timestep using `any(...)`;
  `Parks_NodalBreathing_FrequencyAnalysis.py` appends one per qualifying pair.
  Same stated criterion, different node counts.
- The breathing script's own configuration (m = [1, 1, 0.1], a_out = 5.5,
  e_in = 0.15, e_out = 0.40) is not the one the main README quotes for the
  v15.1 result (m = [1, 1, 0.01], a_out = 4.5, e_in = 0.1, e_out = 0.5).
- `pnep_v15_4_full_engine.py` defines nodes as dr/dt sign changes of the outer
  body against the inner binary CoM and H as (m_in/M)·(r·v²/G) with stable
  ≈ 1.0, where `pnep_v12.py` and the main README use local minima of pairwise
  distance and H = σ²/(1+σ²) with stable ≈ 0.92. Both are called the H-index.
- The 412.6t median lead time appears in `README.md`, both DISCOVERY files and
  the abstract draft, always as prose. No script in the repo generates it and
  no trial output is saved.
- Node detection in `Parks_NodalBreathing_FrequencyAnalysis.py` fires on
  `d_prev1` being the minimum but computes H from `d_cur` and timestamps with
  the current `t`, a systematic one-step lag.

---

## Scope

Three-body planar and mildly inclined hierarchical triples plus seven
equal-mass scramble configurations. Fixed-step leapfrog at dt = 0.004 for the
hierarchical grids, DOP853 at rtol 1e-12 for flash. Observation window 4000t,
so survivors are **right-censored, not negatives** — a system that ejects at
t = 5000 is recorded as a survivor. Defect and jitter statistics rest on 12
ejectors and 18 censored systems; the first-node results on 192 configurations
at two seeds.

---

Alika Parks · November 2026

---

## 7. Compute: what the analytic node saves, and what it does not

`compute_cost.py` — force evaluations, 25 systems, dt = 0.004

| task | median | min | max |
|---|---|---|---|
| A integrate to the first node | 348 | 34 | 765 |
| B Newton predictor | **1** | 1 | 1 |
| C integrate to 150 nodes (classification window) | 55,725 | 6,501 | 64,675 |
| D integrate to ejection (brute force) | 104,126 | 6,501 | 2,609,501 |

**A/B = 348×** on the first-node task, but **A/C = 0.63%** — that task is well
under 1% of the classification window. The analytic first node is a genuine
three-order-of-magnitude saving on something that was never the bottleneck.

**C/D = 55%** on ejecting systems. Reaching enough nodes for a slope costs
roughly half of simply integrating to the answer.

### On the "~99% below full N-body" claim

The claim as written in `DISCOVERY.md` and the main README is about
**H-evaluation overhead**: H is evaluated at ~100 privileged moments rather
than at every timestep. That reproduces exactly — 150 H evaluations against
55,725 timesteps is 0.27%, a **99.73% reduction**. Correct as stated.

What does *not* follow is the unqualified form ("cutting compute by ~99%
against brute-force integration"). The integration still has to happen to
reach the nodes, and C/D = 55% measures that. Same figure, different
denominator. The qualifier is load-bearing.

---

## 8. Zero-integration classification, and the Mardling-Aarseth benchmark

If the first node alone carried the verdict, the ~55,000 force evaluations
would collapse to 1. That is the only route to a total-compute claim, so it
is worth testing properly — against the right benchmark.

**The benchmark is not brute force.** A zero-integration stability criterion
for hierarchical triples already exists: Mardling & Aarseth 2001, which is the
`ma01_ratio` that `pnep_v12.py` already computes from a raw state.

**Sampling matters here.** `pnep_v12.make_ic` *selects* systems by MA01 ratio
(`stable_min=5.0`, `unstable_max=0.45`), so any comparison on that sample is
rigged in MA01's favour. `zero_integration.py` draws 160 systems unstratified
— uniform in a_out ∈ [2.5, 8], e_out ∈ [0, 0.6], e_in ∈ [0, 0.3], log-uniform
m₂ ∈ [0.01, 1], random phase, no MA01 filter — labelled only by integration to
2000t. 43 of 160 eject (26.9%).

| single feature | AUC | | 5-fold CV | AUC |
|---|---|---|---|---|
| `ma01` | **0.982** | | MA01 alone | **0.983 ± 0.007** |
| `rspread` | 0.755 | | node features only | 0.786 ± 0.047 |
| `rmax` | 0.755 | | MA01 + node features | 0.956 ± 0.015 |
| `H0` | 0.753 | | raw orbital elements | 0.974 ± 0.017 |
| `g_absmax` | 0.628 | | | |
| `n_closing` | 0.517 | | | |

**MA01 wins and the node features add nothing.** 0.786 against 0.983, and
combining them is *worse* than MA01 alone (0.956) — eleven weak features
degrading one strong one at n = 160.

As a hard threshold, which is its intended use:

| | accuracy | ejection recall |
|---|---|---|
| MA01 ratio > 0.8 = stable | **92.5%** | **95.3%** |
| MA01 ratio > 1.0 = stable | 77.5% | 100.0% |
| PNEP, randomised phases (§6) | 93.0% | 91.7% |

**A closed-form expression from 2001, needing no integration, matches PNEP on
the classification task.** Survivors are censored at 2000t, so MA01's false
positives at threshold 1.0 may well eject later — longer integration would
improve its numbers, not degrade them.

A hypothesis also refuted: that the first node encodes what MA01 captures.
Node features score 0.786 while the raw elements score 0.974, so the
instantaneous node geometry carries substantially *less* information than the
elements. It is not a cheaper route to the same content.

### What follows for positioning

Stable-versus-unstable classification is the weakest claim in the framework to
lead with, because a one-line formula from 2001 matches it at zero cost. MA01
answers only "will this survive?" — it gives no ejection **time**, no escaping
**body**, and no regime structure.

Those are the results with no analytic competitor:

- median 412.6t of warning (slope), or 308.6t with 79% of lifetime remaining
  (defects, §2)
- which body leaves, via four independent channels (§5 and the main README)
- three distinguishable decay regimes

The compute story is likewise not "fewer force evaluations" but **a verdict
before the event** — which brute-force integration cannot provide at any
price, since it only knows once it is over.
