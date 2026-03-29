[pnep_v14_design_spec.md](https://github.com/user-attachments/files/26327545/pnep_v14_design_spec.md)
# PNEP v14 — Design Specification
## Partner Fidelity Signal for 2+2 Four-Body Systems

**Date:** 2026-03-28  
**Author:** Alika M. Parks — Independent Researcher  
**Builds on:** PNEP v12 (3-body), PNEP v13 (3+1 and 2+2 extensions)

---

## The Core Insight

In a stable 2+2 system, binary partners are gravitationally committed.
At every dynamic cluster moment — the transient instant when four chaotically
or orbitally flying bodies spontaneously re-form into two tight pairs —
the same two bodies always pair together. No exceptions. Forever.

In an unstable system heading toward ejection, one binary begins losing
its grip on its partner. Encounters with the other binary pull bodies away.
Partner switching begins.

**Partner switching = instability. Partner fidelity = stability.**

---

## The Key Physical Principle

**Switching back is not recovery.**

When body 0 briefly returns to pairing with body 1 after having switched
to body 2, that is not stability reasserting itself. That is the chaotic
system momentarily revisiting a familiar configuration before the next
disruption. The hierarchy has already been contested. Any return to the
original pairing is the echo of chaos, not the restoration of order.

Therefore: **the first partner switch is an irreversible warning signal.**

Once a switch is detected, the system is classified as unstable regardless
of subsequent pairings. No rolling average. No threshold tuning. One switch
is the trigger.

---

## Detection Protocol

### Step 1 — Dynamic 2x2 Cluster Detection
At every integration timestep, compute all 6 pairwise distances.
Find the arrangement of 4 bodies into 2 pairs that minimises the sum
of inner-pair distances while maximising inter-pair separation.

A valid cluster fires when:
```
d_inter / d_inner_avg >= threshold_ratio (suggested: 2.5)
```

Where:
- d_inner_avg = mean of the two shortest pair distances
- d_inter = shortest distance between the two pairs

### Step 2 — Record Pairing Identity
At each cluster event, record which two bodies form each pair.
Encode as a frozenset of frozensets for order-invariant comparison:
```python
pairing = frozenset([frozenset(pair_A), frozenset(pair_B)])
```

### Step 3 — Partner Fidelity Check
Compare each new pairing to the first recorded pairing (the reference).

```python
def check_fidelity(pairing_history):
    if len(pairing_history) < 2:
        return True   # insufficient data — assume stable
    reference = pairing_history[0]
    for pairing in pairing_history[1:]:
        if pairing != reference:
            return False  # switch detected — UNSTABLE, irreversible
    return True   # all pairings match reference — STABLE
```

### Step 4 — Classification
- **Stable:** All cluster events show same pairing as reference
- **Unstable:** Any cluster event shows different pairing → immediate warning

### Step 5 — Pre-ejection Window
On first switch detection:
- Record time of switch t_switch
- Fit slope to H_inner values at cluster events (v13 formula)
- Project ejection window using slope extrapolation
- Issue [t_Lo, t_Hi] window immediately

No minimum event count required — one switch is sufficient to classify.
More events before the switch refine the window precision.

---

## Expected Signal Behaviour

```
Stable system cluster history:
  event 1:  (01)(23)  <- reference pairing established
  event 2:  (01)(23)  <- match
  event 3:  (01)(23)  <- match
  event 4:  (01)(23)  <- match
  ...N:     (01)(23)  <- always match → STABLE

Unstable system cluster history:
  event 1:  (01)(23)  <- reference pairing established  
  event 2:  (01)(23)  <- match
  event 3:  (02)(13)  <- SWITCH DETECTED → classify UNSTABLE immediately
  event 4:  (01)(23)  <- switch-back — does NOT cancel warning
  event 5:  (03)(12)  <- further switching confirms decay
```

---

## Why This Is Simpler Than v13

v13 required:
- Rolling H buffer accumulation
- Minimum event count (MIN_CP = 4-8)
- Threshold tuning (H_THRESH)
- Slope projection with multiple parameters

v14 requires:
- One reference pairing (first cluster event)
- Binary comparison at each subsequent event
- One trigger condition (any switch = unstable)

No parameters to tune. No thresholds. No rolling averages.
The signal is categorical not continuous.

---

## H_inner as Supporting Signal

H_inner at cluster moments (v13 finding: stable=0.43, unstable=0.055)
remains available as a secondary confirmation signal and for window
precision:

- Pre-switch H_inner values provide slope for window projection
- H_inner at the switching event is typically near zero — confirming
  the system was already in chaotic geometry at that moment
- Combined signal: fidelity for classification + H_inner for timing

---

## Expected Performance Improvement

v13 2+2 best: accuracy=67%, specificity=100%, sensitivity=65%

v14 hypothesis:
- Specificity: maintain 100% — stable systems never switch partners
- Sensitivity: significant improvement — any switch triggers detection,
  even in fast-ejecting systems that previously escaped before accumulating
  enough events
- Pre-ejection coverage: should approach 100% — first switch fires before
  physical ejection in most cases

The sensitivity ceiling in v13 was caused by fast ejectors escaping before
accumulating MIN_CP events. v14 needs only ONE switch — which occurs earlier
in the degradation sequence than full chaotic mixing.

---

## Open Questions for v14 Implementation

1. **Threshold ratio tuning:** What ratio separates genuine 2x2 clusters
   from random close approaches? Test 2.0, 2.5, 3.0, 4.0.

2. **Reference stability:** Should the reference pairing be the first event,
   or the most common pairing across the first N events? First event is
   simpler but vulnerable to a noisy initial measurement.

3. **Near-switch events:** What if d_inter/d_inner ratio is marginal?
   Consider a confidence-weighted switch detection.

4. **Multi-body ejection:** In some unstable 2+2 systems two bodies eject
   simultaneously. Does partner switching predict this differently?

5. **Combination with v13 coplanar signal:** Can partner fidelity + coplanar
   H_inner be combined for a stronger joint classifier?

---

## Connection to the Broader PNEP Framework

The partner fidelity signal completes the topological picture:

| Topology | Node Type | Signal | Polarity |
|---|---|---|---|
| 3-body | Mirror symmetry node | H = σ²/(1+σ²) | High=stable |
| 3+1 4-body | Mirror symmetry node | H_inner on inner triple | Low=stable |
| 2+2 4-body | Dynamic 2x2 cluster | Partner fidelity + H_inner | Fidelity=stable |

In all three cases the underlying principle is identical:
**At the most geometrically honest moments, stable systems reveal
consistent hierarchical clarity. Unstable systems reveal contested
or scrambled geometry.**

The node type differs by topology because different encounters produce
the most information-dense moments for each configuration. The law is
topology-aware but topology-invariant in its underlying logic.

---

## Repository
https://github.com/alikamp/Parks-Node-Ejection-Protocol

## Status
v13 complete. v14 design specified. Implementation pending.
