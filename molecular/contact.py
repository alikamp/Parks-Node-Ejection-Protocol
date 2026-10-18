"""Is the H signal diluted by triples that are not physical objects?

H_all runs over all C(13,3) = 286 triples.  But a triple of three mutually
distant atoms is not a trio in any dynamical sense -- it is a bookkeeping
artifact.  usplit.py found median u > 0 in every band and essentially no
compressed population, which is what you would see if most triples are not
in contact at all.

In an icosahedral 13-cluster only the triples whose three pairs are all
first-neighbour contacts are real structural units: the faces.  That is of
order 20-40 triples, not 286.  A signal carried by ~30 real triples and
diluted by ~250 artifacts would produce exactly what has been observed --
tight-looking fits whose parameters do not survive reseeding.

Three parts:

1. DESCRIPTIVE.  Pair-distance histogram to locate the first-shell minimum
   empirically rather than guessing a cutoff.  Then: what fraction of
   triples are contact triples, and how does that change with temperature?

2. DISTRIBUTION.  H over contact triples vs H over all triples.  Sarle's
   bimodality coefficient and mode positions.  Is the band structure
   sharper, and are the modes still temperature-invariant?

3. KINETICS.  Chapman-Kolmogorov on the 4-state model restricted to the
   sub-ensemble of triples that are in contact at BOTH t and t+lag -- i.e.
   among triples that persist as physical units, is the H kinetics
   Markovian?  Compared head to head against the all-triples model on the
   same trajectories.
"""
import numpy as np, pickle, time, sys
sys.path.insert(0, '/home/claude/md3')
from cluster import (N, DT, icosahedron13, acc_and_U, H_all, pair_d,
                     PRS, TRI, TRI_P, n_evaporated)

EDGES = np.array([0.0, 0.07, 0.22, 0.60, 1.01])
NB = 4
SAMP = 5
DTS = SAMP * DT
RMIN = 2.0 ** (1.0 / 6.0)


def band_of(h):
    return np.clip(np.digitize(h, EDGES) - 1, 0, NB - 1)


def bimodality(x):
    """Sarle's coefficient.  > 5/9 suggests bimodality."""
    x = np.asarray(x, float)
    n = len(x)
    if n < 4:
        return np.nan
    m = x.mean()
    s = x.std()
    if s == 0:
        return np.nan
    g = (((x - m) / s) ** 3).mean()
    k = (((x - m) / s) ** 4).mean() - 3.0
    return (g ** 2 + 1) / (k + 3.0 * (n - 1) ** 2 / ((n - 2) * (n - 3)))


def modes(x, lo=0.0, hi=1.0, nb=200, prom=0.002):
    """Crude mode finder on a histogram: local maxima above a prominence."""
    h, e = np.histogram(x, bins=nb, range=(lo, hi), density=True)
    c = 0.5 * (e[:-1] + e[1:])
    out = []
    for i in range(1, len(h) - 1):
        if h[i] > h[i - 1] and h[i] >= h[i + 1] and h[i] > prom:
            out.append((c[i], h[i]))
    out.sort(key=lambda t: -t[1])
    return out[:4]


def trajectory(T0, seed, rcut, nsteps=150000, burn=20000):
    """Return per-sample band of every triple, and a contact mask."""
    rng = np.random.default_rng(seed)
    pos = icosahedron13() + rng.normal(0, 0.01, (N, 3))
    vel = rng.normal(0, np.sqrt(T0), (N, 3))
    vel -= vel.mean(0)
    acc, _ = acc_and_U(pos)
    B, M, Hs, dsamp = [], [], [], []
    for s in range(nsteps):
        vel += 0.5 * DT * acc
        pos += DT * vel
        acc, _e = acc_and_U(pos)
        vel += 0.5 * DT * acc
        if s >= burn and s % SAMP == 0:
            d3 = pair_d(pos)[TRI_P]                    # (286,3)
            v = d3.var(axis=1)
            h = v / (1.0 + v)
            B.append(band_of(h).astype(np.int8))
            M.append((d3.max(axis=1) <= rcut))         # all 3 pairs in contact
            Hs.append(h.astype(np.float32))
            if s % 2000 == 0:
                dsamp.append(pair_d(pos).copy())
    return (np.array(B, np.int8), np.array(M, bool),
            np.array(Hs, np.float32), np.array(dsamp))


def tmat(bands, mask, lag, conditioned):
    """Transition matrix at lag.  If conditioned, count only triples in
    contact at BOTH endpoints."""
    C = np.zeros((NB, NB))
    for b, m in zip(bands, mask):
        a = b[:-lag]
        z = b[lag:]
        if conditioned:
            k = m[:-lag] & m[lag:]
            a = a[k]
            z = z[k]
        else:
            a = a.ravel()
            z = z.ravel()
        np.add.at(C, (a.ravel(), z.ravel()), 1)
    rs = C.sum(axis=1, keepdims=True)
    return np.divide(C, rs, out=np.zeros_like(C), where=rs > 0), C


def ck(bands, mask, conditioned, base=8, ns=(2, 4, 8)):
    Tb, _ = tmat(bands, mask, base, conditioned)
    out = {}
    for n in ns:
        Tm, Cm = tmat(bands, mask, base * n, conditioned)
        Tp = np.linalg.matrix_power(Tb, n)
        out[n] = float(np.abs(Tm - Tp).max())
    return out


if __name__ == "__main__":
    TEMPS = [0.65, 0.95, 1.10]
    NSEED = 6

    # ---- 1. locate the contact cutoff empirically ----------------------
    print("=== pair-distance histogram, T0=0.65, to locate first-shell "
          "minimum ===")
    _, _, _, dpilot = trajectory(0.65, 9000, rcut=1.4, nsteps=40000)
    dall = dpilot.ravel()
    hh, ee = np.histogram(dall, bins=80, range=(0.8, 3.5), density=True)
    cc = 0.5 * (ee[:-1] + ee[1:])
    for i in range(0, len(cc), 2):
        bar = '#' * int(hh[i] * 40)
        print(f"  r={cc[i]:5.2f}  {hh[i]:6.3f}  {bar}")
    # first local minimum after the first peak
    ipk = int(np.argmax(hh))
    imin = ipk + int(np.argmin(hh[ipk:ipk + 30]))
    RCUT = float(cc[imin])
    print(f"\n  first peak at r={cc[ipk]:.3f} (r_min={RMIN:.3f}), "
          f"first minimum at r={RCUT:.3f}  -> contact cutoff")

    res = {'rcut': RCUT}
    for T0 in TEMPS:
        t0 = time.time()
        Bs, Ms, Hs = [], [], []
        for sd in range(NSEED):
            b, m, h, _ = trajectory(T0, 5000 + sd, RCUT)
            Bs.append(b)
            Ms.append(m)
            Hs.append(h)
        hall = np.concatenate([h.ravel() for h in Hs])
        mall = np.concatenate([m.ravel() for m in Ms])
        hcon = hall[mall]

        print(f"\n=== T0 = {T0:.2f}  [{time.time()-t0:.0f}s] ===")
        print(f"  contact triples: {mall.mean()*100:.2f}% of 286 "
              f"= {mall.mean()*286:.1f} on average")

        print(f"\n  {'':14}{'n':>12}{'mean H':>9}{'Sarle b':>9}  modes")
        for nm, x in (('all triples', hall), ('contact only', hcon)):
            md = modes(x)
            ms = " ".join(f"{p:.3f}" for p, _ in md)
            print(f"  {nm:<14}{len(x):12d}{x.mean():9.4f}"
                  f"{bimodality(x[::37]):9.4f}  {ms}")

        print(f"\n  {'model':<22}{'CK n=2':>9}{'CK n=4':>9}{'CK n=8':>9}")
        a = ck(Bs, Ms, conditioned=False)
        c = ck(Bs, Ms, conditioned=True)
        print(f"  {'all triples (286)':<22}{a[2]:9.4f}{a[4]:9.4f}{a[8]:9.4f}")
        print(f"  {'contact triples':<22}{c[2]:9.4f}{c[4]:9.4f}{c[8]:9.4f}")
        imp = [(a[n] - c[n]) / a[n] * 100 if a[n] > 0 else np.nan
               for n in (2, 4, 8)]
        print(f"  {'improvement':<22}" + "".join(f"{v:8.1f}%" for v in imp))

        res[T0] = dict(fcon=float(mall.mean()), ck_all=a, ck_con=c,
                       modes_all=modes(hall), modes_con=modes(hcon))

    pickle.dump(res, open('/home/claude/md10/contact.pkl', 'wb'))
    print("\nsaved")
