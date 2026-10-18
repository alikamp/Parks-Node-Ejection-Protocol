"""Does a second coordinate recover the Markov property?

markov.py: the 4-state H model is not Markovian.  Implied timescales drift
2.5x without plateauing; Chapman-Kolmogorov errors reach 0.25-0.35.

Mechanism hypothesis, from section 1: H IS SIGN-BLIND.  H = var/(1+var) is
symmetric, so a trio at H = 0.15 may be compressed or stretched and H cannot
distinguish them.  Two physically distinct populations sharing one H value is
exactly what produces non-exponential waiting times and drifting timescales.

Second coordinate: the equilibrium-referenced displacement
    u = mean over the trio's three pairs of (r - r_min)/r_min,  r_min = 2^(1/6)
Negative u = compressed, positive u = stretched.

CONTROLS MATTER.  Going from 4 to 8 states improves Chapman-Kolmogorov for
free, just by describing the system more finely.  Four models at matched
state count:

  M1  4 states   H bands                                  baseline
  M2  8 states   each H band split at its median u        the hypothesis
  M3  8 states   each H band split at its median H        resolution control
  M4  8 states   each H band split at a SHUFFLED median u count control

M2 > M3 means u carries information that finer H resolution does not -- the
claim.  M2 ~ M3 means any second cut on the same information would do.
M3 ~ M4 means nothing but state count is operating.

All three splits are median splits, so all produce equal populations by
construction and differ only in which coordinate does the splitting.
"""
import numpy as np, pickle, time, sys
sys.path.insert(0, '/home/claude/md3')
from cluster import (N, DT, icosahedron13, acc_and_U, H_all, pair_d,
                     PRS, TRI_P, lindemann, n_evaporated)

EDGES = np.array([0.0, 0.07, 0.22, 0.60, 1.01])
NB = 4
SAMP = 5
DTS = SAMP * DT
RMIN = 2.0 ** (1.0 / 6.0)


def band_of(h):
    return np.clip(np.digitize(h, EDGES) - 1, 0, NB - 1)


def h_and_u(pos):
    """H and mean signed displacement u for all 286 triples."""
    d = pair_d(pos)[TRI_P]                     # (286,3)
    v = d.var(axis=1)
    h = v / (1.0 + v)
    u = ((d - RMIN) / RMIN).mean(axis=1)
    return h, u


def trajectory(T0, seed, nsteps=200000, burn=20000):
    rng = np.random.default_rng(seed)
    pos = icosahedron13() + rng.normal(0, 0.01, (N, 3))
    vel = rng.normal(0, np.sqrt(T0), (N, 3))
    vel -= vel.mean(0)
    acc, _ = acc_and_U(pos)
    H, U = [], []
    for s in range(nsteps):
        vel += 0.5 * DT * acc
        pos += DT * vel
        acc, _e = acc_and_U(pos)
        vel += 0.5 * DT * acc
        if s >= burn and s % SAMP == 0:
            h, u = h_and_u(pos)
            H.append(h.astype(np.float16))
            U.append(u.astype(np.float16))
    return np.array(H, dtype=np.float16), np.array(U, dtype=np.float16)


def tmatrix(states, nstate, lag):
    C = np.zeros((nstate, nstate))
    for st in states:
        np.add.at(C, (st[:-lag].ravel(), st[lag:].ravel()), 1)
    rs = C.sum(axis=1, keepdims=True)
    return np.divide(C, rs, out=np.zeros_like(C), where=rs > 0), C


def timescales(T, lag, k=3):
    ev = np.linalg.eigvals(T)
    ev = np.real(ev[np.argsort(-np.real(ev))])
    out = []
    for lam in ev[1:1 + k]:
        out.append(-lag * DTS / np.log(lam) if 0 < lam < 1 else np.nan)
    return np.array(out)


def ck_error(states, nstate, base=8, ns=(2, 4, 8)):
    Tb, _ = tmatrix(states, nstate, base)
    out = {}
    for n in ns:
        Tm, _ = tmatrix(states, nstate, base * n)
        Tp = np.linalg.matrix_power(Tb, n)
        out[n] = float(np.abs(Tm - Tp).max())
    return out


def split_states(bands, vals, med):
    """band*2 + (value > that band's median)."""
    return [b.astype(np.int16) * 2
            + (v.astype(np.float32) > med[b]).astype(np.int16)
            for b, v in zip(bands, vals)]


if __name__ == "__main__":
    TEMPS = [0.65, 0.95, 1.10]
    NSEED = 6
    res = {}
    t0 = time.time()

    for T0 in TEMPS:
        Hs, Us = [], []
        for sd in range(NSEED):
            h, u = trajectory(T0, 4000 + sd)
            Hs.append(h)
            Us.append(u)
        print(f"T0={T0:.2f} simulated  [{time.time()-t0:.0f}s]", flush=True)

        Bs = [band_of(h.astype(np.float32)).astype(np.int8) for h in Hs]
        ball = np.concatenate([b.ravel() for b in Bs])
        uall = np.concatenate([u.ravel() for u in Us]).astype(np.float32)
        hall = np.concatenate([h.ravel() for h in Hs]).astype(np.float32)

        med_u = np.zeros(NB)
        med_h = np.zeros(NB)
        frac_neg = np.zeros(NB)
        pop = np.zeros(NB)
        for b in range(NB):
            m = ball == b
            pop[b] = m.mean()
            if m.sum() > 0:
                med_u[b] = np.median(uall[m])
                med_h[b] = np.median(hall[m])
                frac_neg[b] = float((uall[m] < 0).mean())

        # M4: shuffle u within each trajectory, destroying its dynamical
        # correlation while preserving the marginal distribution.
        rng = np.random.default_rng(7)
        Ush = []
        for u in Us:
            flat = u.astype(np.float32).ravel().copy()
            rng.shuffle(flat)
            Ush.append(flat.reshape(u.shape))

        models = {
            'M1 H bands  (4)': (Bs, NB),
            'M2 H x u    (8)': (split_states(Bs, Us, med_u), 2 * NB),
            'M3 H x H    (8)': (split_states(Bs, Hs, med_h), 2 * NB),
            'M4 H x shuf (8)': (split_states(Bs, Ush, med_u), 2 * NB),
        }

        print(f"\n=== T0 = {T0:.2f} ===")
        print("band  population   median u   median H   frac compressed (u<0)")
        for b in range(NB):
            print(f"  {'ABCD'[b]}    {pop[b]:9.4f}  {med_u[b]:+9.4f}"
                  f"  {med_h[b]:9.4f}   {frac_neg[b]:9.3f}")

        res[T0] = dict(med_u=med_u, med_h=med_h, frac_neg=frac_neg, pop=pop,
                       models={})
        print(f"\n{'model':<17}{'CK n=2':>9}{'CK n=4':>9}{'CK n=8':>9}"
              f"{'t2(lag1)':>11}{'t2(lag64)':>11}{'drift':>10}")
        for name, (st, ns) in models.items():
            ck = ck_error(st, ns)
            T1, _ = tmatrix(st, ns, 1)
            T64, _ = tmatrix(st, ns, 64)

            def slow(a):
                f = a[np.isfinite(a)]
                return f[0] if len(f) else np.nan
            s1 = slow(timescales(T1, 1))
            s64 = slow(timescales(T64, 64))
            dr = abs(s64 - s1) / s1 if np.isfinite(s1) and s1 > 0 else np.nan
            res[T0]['models'][name] = dict(ck=ck, s1=s1, s64=s64, drift=dr)
            print(f"{name:<17}{ck[2]:9.4f}{ck[4]:9.4f}{ck[8]:9.4f}"
                  f"{s1:11.3f}{s64:11.3f}"
                  + (f"{dr*100:9.1f}%" if np.isfinite(dr) else f"{'--':>10}"))
        print()

    pickle.dump(res, open('/home/claude/md10/usplit.pkl', 'wb'))
    print("saved")
