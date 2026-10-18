"""Is the Arrhenius activation energy an artifact of the sampling stride?

markov.py showed the four-state H model is NOT Markovian: implied timescales
drift by 2.5x over lags 1-64 without plateauing, and Chapman-Kolmogorov fails
at max|T(n.tau) - T(tau)^n| = 0.27.  So k_ij are not Markov rates; they are
lag-dependent empirical crossing frequencies.

That does not automatically kill the Arrhenius result.  The absolute rates
must fall with lag (multi-step paths get absorbed into single transitions),
but the TEMPERATURE DEPENDENCE at fixed lag may still be a real physical
quantity.  The test:

    fit ln k_ij(tau) = ln A - Ea/T   separately at each lag tau
    -> if Ea(tau) is flat, the activation energy is robust
    -> if Ea(tau) drifts, 1.124 eps was an artifact of SAMP=5

Also reports R^2(tau), so a fit that degrades at long lag is visible.

5 temperatures (needed for the fit) x 8 seeds, lags 1..64 samples.
"""
import numpy as np, pickle, time, sys
from scipy import stats
sys.path.insert(0, '/home/claude/md3')
from cluster import (N, DT, icosahedron13, acc_and_U, H_all, pair_d,
                     PRS, lindemann, n_evaporated)

EDGES = np.array([0.0, 0.07, 0.22, 0.60, 1.01])
NS = 4
NAMES = ['A', 'B', 'C', 'D']
SAMP = 5
DTS = SAMP * DT
LAGS = [1, 2, 4, 8, 16, 32, 64]
TRANS = [(0, 1), (1, 0), (1, 2), (2, 1), (2, 3)]


def state(h):
    return np.clip(np.digitize(h, EDGES) - 1, 0, NS - 1)


def trajectory(T0, seed, nsteps=200000, burn=20000):
    rng = np.random.default_rng(seed)
    pos = icosahedron13() + rng.normal(0, 0.01, (N, 3))
    vel = rng.normal(0, np.sqrt(T0), (N, 3))
    vel -= vel.mean(0)
    acc, _ = acc_and_U(pos)
    out = []
    for s in range(nsteps):
        vel += 0.5 * DT * acc
        pos += DT * vel
        acc, U = acc_and_U(pos)
        vel += 0.5 * DT * acc
        if s >= burn and s % SAMP == 0:
            out.append(state(H_all(pos)).astype(np.int8))
    return np.array(out, dtype=np.int8)


def counts(trajs, lag):
    C = np.zeros((NS, NS))
    for tr in trajs:
        np.add.at(C, (tr[:-lag].ravel(), tr[lag:].ravel()), 1)
    return C


def krate(C, lag):
    """k_ij at this lag: transitions per unit time spent in i."""
    K = np.zeros((NS, NS))
    for i in range(NS):
        ti = C[i].sum() * lag * DTS
        if ti > 0:
            K[i] = C[i] / ti
    return K


if __name__ == "__main__":
    TEMPS = [0.50, 0.65, 0.80, 0.95, 1.10]
    NSEED = 8
    K = {}                       # K[lag][T] = rate matrix
    for lag in LAGS:
        K[lag] = {}

    t0 = time.time()
    for T0 in TEMPS:
        trajs = [trajectory(T0, 3000 + sd) for sd in range(NSEED)]
        for lag in LAGS:
            K[lag][T0] = krate(counts(trajs, lag), lag)
        del trajs
        print(f"T0={T0:.2f} done  [{time.time()-t0:.0f}s]", flush=True)

    # ---- Ea as a function of lag ---------------------------------------
    print("\n=== Ea(lag) : is the activation energy stride-dependent? ===")
    hdr = f"{'transition':<10}" + "".join(f"{l*DTS:>9.2f}" for l in LAGS)
    print(f"{'':<10}" + "".join(f"{'lag=':>9}" for l in LAGS))
    print(hdr)
    Ea_tab = {}
    for (i, j) in TRANS:
        row = []
        for lag in LAGS:
            ks = np.array([K[lag][T][i, j] for T in TEMPS])
            m = ks > 0
            if m.sum() >= 3:
                sl, ic, r, p, se = stats.linregress(
                    1.0 / np.array(TEMPS)[m], np.log(ks[m]))
                row.append(-sl)
            else:
                row.append(np.nan)
        Ea_tab[(i, j)] = row
        print(f"  {NAMES[i]}->{NAMES[j]}    "
              + "".join(f"{v:9.3f}" if np.isfinite(v) else f"{'--':>9}"
                        for v in row))

    print("\n=== R^2(lag) ===")
    print(hdr)
    for (i, j) in TRANS:
        row = []
        for lag in LAGS:
            ks = np.array([K[lag][T][i, j] for T in TEMPS])
            m = ks > 0
            if m.sum() >= 3:
                sl, ic, r, p, se = stats.linregress(
                    1.0 / np.array(TEMPS)[m], np.log(ks[m]))
                row.append(r ** 2)
            else:
                row.append(np.nan)
        print(f"  {NAMES[i]}->{NAMES[j]}    "
              + "".join(f"{v:9.4f}" if np.isfinite(v) else f"{'--':>9}"
                        for v in row))

    print("\n=== drift of Ea from shortest to longest lag ===")
    for (i, j) in TRANS:
        r = np.array(Ea_tab[(i, j)], float)
        f = r[np.isfinite(r)]
        if len(f) >= 2:
            rel = abs(f[-1] - f[0]) / abs(f[0]) if f[0] != 0 else np.nan
            print(f"  {NAMES[i]}->{NAMES[j]}   {f[0]:+.3f} -> {f[-1]:+.3f}"
                  f"   relative drift {rel*100:6.1f}%")

    print("\n=== absolute rates at T0=0.95 (expected to fall with lag) ===")
    print(hdr)
    for (i, j) in TRANS:
        print(f"  {NAMES[i]}->{NAMES[j]}    "
              + "".join(f"{K[l][0.95][i,j]:9.4f}" for l in LAGS))

    pickle.dump(dict(K=K, Ea=Ea_tab, TEMPS=TEMPS, LAGS=LAGS),
                open('/home/claude/md10/ea_lag.pkl', 'wb'))
    print("\nsaved")
