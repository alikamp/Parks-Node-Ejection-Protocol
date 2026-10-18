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

This is not a demonstration that PNEP generalises to molecular dynamics. It is
a demonstration that **the node-scalar approach transfers to Lennard-Jones
dissociation within a bounded regime**, that the boundary is physically
identifiable in advance via tau, and that the transfer requires an
equilibrium-referenced coordinate to avoid inverting on repulsive
configurations.

Single potential form, three bodies, one trimer geometry, small per-cell
samples. Generalisation to larger clusters, other potentials, solvated systems
or biomolecular force fields is untested.

---

Alika Parks · October 2026
