"""Pin the ejection rate constant k_CD.

Context: kinetics.py established that the H-band system is an equilibrated
reversible ladder A<->B<->C (detailed balance to 0.01-0.08%) feeding a single
absorbing sink D.  All four ladder rate constants fit Arrhenius at R^2>0.98.
The one rate that matters for ejection -- C->D -- had only 243-521 events and
gave R^2=0.287.

Ejection flux = occ(C) * k_CD.  occ(C) is pure geometry from H and is well
determined.  If k_CD is temperature-flat, ejection rate is predictable from
geometry alone.  This run decides that.

Changes from kinetics.py:
  - hotter range (more C->D events), 40 seeds
  - dwell = occupancy_time / n_exits   (kinetics.py mis-scaled this)
  - Poisson counting errors on every rate
  - evaporation tracked over time, so rates can be restricted to the window
    where the cluster is still a cluster
"""
import numpy as np, pickle, time, sys
sys.path.insert(0, '/home/claude/md3')
from cluster import (N, DT, icosahedron13, acc_and_U, H_all, pair_d,
                     PRS, lindemann, n_evaporated)

EDGES = np.array([0.0, 0.07, 0.22, 0.60, 1.01])
NS = 4
NAMES = ['A compact', 'B bulk', 'C stretched', 'D detached']
SAMP = 5                      # sample every 5 steps = 0.02 time units
DTS = SAMP * DT               # time between samples
MAX_EVAP = 2                  # cluster identity: at most 2 atoms gone


def state(h):
    return np.clip(np.digitize(h, EDGES) - 1, 0, NS - 1)


def run(T0, seed, nsteps=200000, burn=20000):
    """Return transition counts, occupancy counts, exit counts, and the same
    three restricted to samples where the cluster is still intact."""
    rng = np.random.default_rng(seed)
    pos = icosahedron13() + rng.normal(0, 0.01, (N, 3))
    vel = rng.normal(0, np.sqrt(T0), (N, 3))
    vel -= vel.mean(0)
    acc, _ = acc_and_U(pos)

    C = np.zeros((NS, NS))          # all samples
    Ci = np.zeros((NS, NS))         # intact-cluster samples only
    occ = np.zeros(NS)
    occi = np.zeros(NS)
    prev = None
    nev_trace = []

    for s in range(nsteps):
        vel += 0.5 * DT * acc
        pos += DT * vel
        acc, U = acc_and_U(pos)
        vel += 0.5 * DT * acc
        if s >= burn and s % SAMP == 0:
            st = state(H_all(pos))
            occ += np.bincount(st, minlength=NS)
            nev = n_evaporated(pos)
            intact = nev <= MAX_EVAP
            if intact:
                occi += np.bincount(st, minlength=NS)
            if prev is not None:
                np.add.at(C, (prev, st), 1)
                if intact and prev_intact:
                    np.add.at(Ci, (prev, st), 1)
            prev = st
            prev_intact = intact
            if s % 5000 == 0:
                nev_trace.append(nev)

    return C, Ci, occ, occi, np.array(nev_trace)


def rates(C):
    """k_ij with Poisson errors.  Denominator is row sum = true occupancy."""
    K = np.zeros((NS, NS))
    E = np.zeros((NS, NS))
    for i in range(NS):
        ti = C[i].sum() * DTS
        if ti > 0:
            for j in range(NS):
                if i != j:
                    n = C[i, j]
                    K[i, j] = n / ti
                    E[i, j] = np.sqrt(max(n, 1.0)) / ti
    return K, E


def dwell(C):
    """Mean residence time = occupancy time / number of exits."""
    d = np.full(NS, np.nan)
    for i in range(NS):
        exits = C[i].sum() - C[i, i]
        if exits > 0:
            d[i] = C[i].sum() * DTS / exits
    return d


if __name__ == "__main__":
    TEMPS = [0.95, 1.10, 1.25, 1.40, 1.55]
    NSEED = 40
    res = {}
    t0 = time.time()
    for T0 in TEMPS:
        Ct = np.zeros((NS, NS))
        Cit = np.zeros((NS, NS))
        occt = np.zeros(NS)
        occit = np.zeros(NS)
        evs = []
        for sd in range(NSEED):
            C, Ci, occ, occi, nev = run(T0, 1000 + sd)
            Ct += C
            Cit += Ci
            occt += occ
            occit += occi
            evs.append(nev)
        K, E = rates(Ct)
        Ki, Ei = rates(Cit)
        res[T0] = dict(C=Ct, Ci=Cit,
                       occ=occt / occt.sum(), occi=occit / max(occit.sum(), 1),
                       K=K, E=E, Ki=Ki, Ei=Ei,
                       dwell=dwell(Ct), dwelli=dwell(Cit),
                       nev=np.mean([e[-1] for e in evs]),
                       nev_mid=np.mean([e[len(e) // 2] for e in evs]))
        print(f"T0={T0:.2f}  occC={res[T0]['occ'][2]:.4f}  "
              f"nCD={Ct[2,3]:.0f}  k_CD={K[2,3]:.5f}+/-{E[2,3]:.5f}  "
              f"| intact: occC={res[T0]['occi'][2]:.4f} nCD={Cit[2,3]:.0f} "
              f"k_CD={Ki[2,3]:.5f}+/-{Ei[2,3]:.5f}  "
              f"| <nevap>={res[T0]['nev']:.2f}  [{time.time()-t0:.0f}s]",
              flush=True)
    pickle.dump(res, open('/home/claude/md10/kcd.pkl', 'wb'))
    print("saved")
