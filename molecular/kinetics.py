"""H as a discrete multi-state kinetic system.

Section 11 binned H on an arbitrary grid and found net flux structurally zero.
This asks a sharper question: is the H distribution actually discrete, and if
so what are the rate constants between its states?

The H distribution over all 286 triples is multimodal with modes near 0.038
and 0.113 that do NOT move with temperature -- only their occupancies change
-- plus a third mode near 0.96 that appears when evaporation starts.  That
justifies a discrete-state treatment rather than a continuum.

States are cut at the density minima, not at round numbers:

  A  H < 0.07      compact
  B  0.07 - 0.22   bulk
  C  0.22 - 0.60   stretched
  D  H >= 0.60     detached

Reports transition counts, rate constants k_ij = (counts i->j) / (time in i),
Arrhenius fits, occupancies, dwell times, and an explicit detailed-balance
audit on every adjacent pair.

12 seeds per temperature, 200k steps, sampled every 5 steps.
"""
import numpy as np, pickle, time, sys
sys.path.insert(0, '/home/claude/md3')
from cluster import (N, DT, icosahedron13, acc_and_U, H_all, pair_d,
                     PRS, lindemann, n_evaporated)

EDGES = np.array([0.0, 0.07, 0.22, 0.60, 1.01])
NS = 4
NAMES = ['A compact', 'B bulk', 'C stretched', 'D detached']
SAMP = 5                 # sample every 5 steps = 0.02 time units
DTS = SAMP * DT          # time between samples


def state(h):
    return np.clip(np.digitize(h, EDGES) - 1, 0, NS - 1)


def run(T0, seed, nsteps=200000, burn=20000):
    """Accumulate the full NSxNS transition-count matrix over all 286 triples.

    Row sums of C are occupancy counts, so C[i].sum()*DTS is the true time
    spent in state i and is the correct denominator for k_ij.
    """
    rng = np.random.default_rng(seed)
    pos = icosahedron13() + rng.normal(0, 0.01, (N, 3))
    vel = rng.normal(0, np.sqrt(T0), (N, 3))
    vel -= vel.mean(0)
    acc, _ = acc_and_U(pos)

    C = np.zeros((NS, NS))
    prev = None
    for s in range(nsteps):
        vel += 0.5 * DT * acc
        pos += DT * vel
        acc, U = acc_and_U(pos)
        vel += 0.5 * DT * acc
        if s >= burn and s % SAMP == 0:
            st = state(H_all(pos))
            if prev is not None:
                np.add.at(C, (prev, st), 1)
            prev = st
    return C, n_evaporated(pos)


def rates(C):
    """k_ij with Poisson counting errors.  Denominator is the row sum."""
    K = np.zeros((NS, NS))
    E = np.zeros((NS, NS))
    for i in range(NS):
        ti = C[i].sum() * DTS
        if ti > 0:
            for j in range(NS):
                if i != j:
                    K[i, j] = C[i, j] / ti
                    E[i, j] = np.sqrt(max(C[i, j], 1.0)) / ti
    return K, E


def dwell(C):
    """Mean residence time = time in state / number of exits from it.

    NOTE: an earlier version of this script divided by a sample count rather
    than an exit count and reported dwell times that were wrong by a
    state-dependent factor, inverting the temperature trend for state A.
    """
    d = np.full(NS, np.nan)
    for i in range(NS):
        exits = C[i].sum() - C[i, i]
        if exits > 0:
            d[i] = C[i].sum() * DTS / exits
    return d


def balance(C):
    """Fractional flux imbalance |n_ij - n_ji| / n_ij for adjacent pairs."""
    out = {}
    for i, j in [(0, 1), (1, 2), (2, 3)]:
        if C[i, j] > 0:
            out[f'{NAMES[i][0]}<->{NAMES[j][0]}'] = \
                abs(C[i, j] - C[j, i]) / C[i, j]
        else:
            out[f'{NAMES[i][0]}<->{NAMES[j][0]}'] = np.nan
    return out


def arrhenius(temps, ks):
    """ln k = ln A - Ea/T.  Returns Ea, lnA, R^2."""
    from scipy import stats
    t = np.asarray(temps, float)
    k = np.asarray(ks, float)
    m = k > 0
    if m.sum() < 3:
        return np.nan, np.nan, np.nan
    sl, ic, r, p, se = stats.linregress(1.0 / t[m], np.log(k[m]))
    return -sl, ic, r ** 2


if __name__ == "__main__":
    TEMPS = [0.50, 0.65, 0.80, 0.95, 1.10]
    NSEED = 12
    res = {}
    t0 = time.time()

    for T0 in TEMPS:
        Ct = np.zeros((NS, NS))
        ev = []
        for sd in range(NSEED):
            C, nev = run(T0, sd)
            Ct += C
            ev.append(nev)
        K, E = rates(Ct)
        res[T0] = dict(C=Ct, K=K, E=E,
                       occ=Ct.sum(axis=1) / Ct.sum(),
                       dwell=dwell(Ct), bal=balance(Ct),
                       nev=float(np.mean(ev)))
        print(f"T0={T0:.2f}  occ=" +
              " ".join(f"{v:.4f}" for v in res[T0]['occ']) +
              f"  nCD={Ct[2,3]:.0f}  [{time.time()-t0:.0f}s]", flush=True)

    # ---- rate constants -------------------------------------------------
    print("\n=== RATE CONSTANTS k_ij (per time unit) ===")
    print(f"{'transition':<14}" + "".join(f"{t:>8.2f}" for t in TEMPS))
    for i, j in [(0, 1), (1, 0), (1, 2), (2, 1), (2, 3), (3, 2)]:
        lbl = f"  {NAMES[i][0]}->{NAMES[j][0]}"
        print(f"{lbl:<14}" + "".join(f"{res[t]['K'][i,j]:8.3f}" for t in TEMPS))

    # ---- Arrhenius ------------------------------------------------------
    print("\n=== ARRHENIUS: ln k = ln A - Ea/T  (Ea in units of eps) ===")
    print(f"{'transition':<14}{'Ea':>9}{'lnA':>9}{'R^2':>8}{'n':>4}")
    for i, j in [(0, 1), (1, 0), (1, 2), (2, 1), (2, 3)]:
        ks = [res[t]['K'][i, j] for t in TEMPS]
        n = sum(1 for k in ks if k > 0)
        Ea, lnA, r2 = arrhenius(TEMPS, ks)
        lbl = f"  {NAMES[i][0]}->{NAMES[j][0]}"
        if np.isnan(Ea):
            print(f"{lbl:<14}{'--':>9}{'--':>9}{'--':>8}{n:>4}")
        else:
            print(f"{lbl:<14}{Ea:9.3f}{lnA:9.3f}{r2:8.4f}{n:>4}")

    # ---- detailed balance ----------------------------------------------
    print("\n=== DETAILED BALANCE (fractional flux imbalance) ===")
    keys = list(res[TEMPS[0]]['bal'].keys())
    print(f"{'T0':>6}" + "".join(f"{k:>12}" for k in keys))
    for t in TEMPS:
        b = res[t]['bal']
        print(f"{t:6.2f}" + "".join(
            f"{b[k]*100:11.2f}%" if not np.isnan(b[k]) else f"{'--':>12}"
            for k in keys))

    # ---- occupancy and dwell -------------------------------------------
    print("\n=== OCCUPANCY and DWELL TIME ===")
    print(f"{'T0':>6}" + "".join(f"{n[0]+' occ':>9}" for n in NAMES)
          + "".join(f"{n[0]+' dwl':>10}" for n in NAMES) + f"{'<nevap>':>9}")
    for t in TEMPS:
        r = res[t]
        print(f"{t:6.2f}"
              + "".join(f"{v:9.4f}" for v in r['occ'])
              + "".join(f"{v:10.2f}" if not np.isnan(v) else f"{'--':>10}"
                        for v in r['dwell'])
              + f"{r['nev']:9.2f}")

    pickle.dump(res, open('/home/claude/md10/kinetics.pkl', 'wb'))
    print("\nsaved")
