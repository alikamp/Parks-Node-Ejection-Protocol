"""Is the four-state H model actually Markovian?

kinetics.py reports rate constants k_ij and Arrhenius activation energies.
Those numbers only mean what they are claimed to mean if the state process
is Markovian on this state space.  That was assumed, not tested.

Two standard Markov-state-model validations:

1. IMPLIED TIMESCALES.  t_i(tau) = -tau / ln(lambda_i(tau)), where lambda_i
   are eigenvalues of the lag-tau transition matrix.  For a Markov process
   these are independent of tau.  If they drift upward and then plateau, the
   model is Markovian above the plateau lag and the rate constants should be
   quoted there, not at tau = 1.

2. CHAPMAN-KOLMOGOROV.  Measured T(n*tau) against predicted T(tau)^n.  A
   Markov process satisfies these exactly up to counting noise.

Also reports the stationary distribution from the leading left eigenvector,
to check it against the directly measured occupancies.

Stores full state trajectories so every lag is computed from one simulation.
"""
import numpy as np, pickle, time, sys
sys.path.insert(0, '/home/claude/md3')
from cluster import (N, DT, icosahedron13, acc_and_U, H_all, pair_d,
                     PRS, lindemann, n_evaporated)

EDGES = np.array([0.0, 0.07, 0.22, 0.60, 1.01])
NS = 4
NAMES = ['A', 'B', 'C', 'D']
SAMP = 5
DTS = SAMP * DT                       # 0.02 time units per stored sample
LAGS = [1, 2, 4, 8, 16, 32, 64]       # in samples


def state(h):
    return np.clip(np.digitize(h, EDGES) - 1, 0, NS - 1)


def trajectory(T0, seed, nsteps=200000, burn=20000):
    """Return (nsamples, 286) int8 array of states."""
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


def tmatrix(trajs, lag):
    """Row-stochastic transition matrix at the given lag (in samples)."""
    C = np.zeros((NS, NS))
    for tr in trajs:
        a = tr[:-lag].ravel()
        b = tr[lag:].ravel()
        np.add.at(C, (a, b), 1)
    rs = C.sum(axis=1, keepdims=True)
    return np.divide(C, rs, out=np.zeros_like(C), where=rs > 0), C


def timescales(T, lag):
    """Implied timescales in time units, slowest first.  Excludes lambda=1."""
    ev = np.linalg.eigvals(T)
    ev = np.real(ev[np.argsort(-np.real(ev))])
    ts = []
    for lam in ev[1:]:
        if 0 < lam < 1:
            ts.append(-lag * DTS / np.log(lam))
        else:
            ts.append(np.nan)
    return np.array(ts)


def stationary(T):
    """Leading left eigenvector, normalised."""
    w, v = np.linalg.eig(T.T)
    i = np.argmin(np.abs(w - 1.0))
    p = np.real(v[:, i])
    return p / p.sum()


if __name__ == "__main__":
    TEMPS = [0.65, 0.95, 1.10]
    NSEED = 8
    out = {}
    t0 = time.time()

    for T0 in TEMPS:
        trajs = [trajectory(T0, 2000 + sd) for sd in range(NSEED)]
        print(f"T0={T0:.2f} simulated  [{time.time()-t0:.0f}s]", flush=True)

        # --- implied timescales vs lag -------------------------------
        its = {}
        for lag in LAGS:
            T, C = tmatrix(trajs, lag)
            its[lag] = dict(T=T, C=C, ts=timescales(T, lag),
                            pi=stationary(T))

        # --- Chapman-Kolmogorov at base lag --------------------------
        base = 8
        Tb, _ = tmatrix(trajs, base)
        ck = {}
        for n in [2, 4, 8]:
            Tm, Cm = tmatrix(trajs, base * n)
            Tp = np.linalg.matrix_power(Tb, n)
            ck[n] = dict(meas=Tm, pred=Tp,
                         err=np.abs(Tm - Tp).max(),
                         relerr=np.abs(Tm - Tp)[Tm > 1e-4].max()
                         if (Tm > 1e-4).any() else np.nan)

        occ = np.bincount(np.concatenate([t.ravel() for t in trajs]),
                          minlength=NS).astype(float)
        occ /= occ.sum()
        out[T0] = dict(its=its, ck=ck, occ=occ, base=base)

        # --- report ---------------------------------------------------
        print(f"\n=== T0 = {T0:.2f} : IMPLIED TIMESCALES (time units) ===")
        print(f"{'lag':>5}{'lag(t)':>9}{'t2':>10}{'t3':>10}{'t4':>10}"
              f"   stationary pi")
        for lag in LAGS:
            d = its[lag]
            ts = d['ts']
            pi = d['pi']
            print(f"{lag:5d}{lag*DTS:9.2f}"
                  + "".join(f"{v:10.3f}" if np.isfinite(v) else f"{'--':>10}"
                            for v in ts)
                  + "   " + " ".join(f"{p:.3f}" for p in pi))
        print(f"      measured occupancy:"
              + " " * 32 + " ".join(f"{p:.3f}" for p in occ))

        print(f"\n=== T0 = {T0:.2f} : CHAPMAN-KOLMOGOROV (base lag "
              f"{base} = {base*DTS:.2f} t.u.) ===")
        print(f"{'n':>4}{'lag(t)':>9}{'max |T(n.tau) - T(tau)^n|':>28}")
        for n in [2, 4, 8]:
            print(f"{n:4d}{base*n*DTS:9.2f}{ck[n]['err']:28.5f}")

        print()

    pickle.dump(out, open('/home/claude/md10/markov.pkl', 'wb'))
    print("saved")
