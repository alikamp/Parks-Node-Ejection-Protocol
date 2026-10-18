"""The H correction, revisited: which reference actually normalises H?

Section 1 tested normalising distances by sigma_ij or r_min in the 3-body
HETERONUCLEAR trimer and it failed (bound/anti-bound gap 0.03-0.13).  What
recovered the information was min(u) as a separate second scalar.

Two things this script establishes.

PART A.  In the monatomic LJ-13 cluster the correction is a NO-OP.  All atoms
are identical, so sigma_ij = 1 and r_min_ij = 2^(1/6) for all 78 pairs:
normalising divides every distance by one constant and cannot change any
ordering.  Verified explicitly rather than argued.

PART B.  The reference that DOES work is not the pair-force equilibrium
r_min but the cluster's own equilibrium geometry.  Section 17 showed the
perfect icosahedron already reproduces the band occupancies, so define the
per-triple strain against the static structure:

    e_ij  = (d_ij - d_ij^static) / d_ij^static
    H_eps = Var(e over the triple's three pairs) / (1 + Var)

H_eps is zero for EVERY triple at T = 0 by construction, not just for bands
C and D.  That is a properly normalised H.

Compared head to head as state variables against cluster excitation
(U - U_min, U_min = -44.327), within temperature, at solid-like temperatures
where the static site mapping is valid.  A permutation diagnostic reports
where that mapping breaks down.
"""
import numpy as np, pickle, time, sys
from scipy import stats
sys.path.insert(0, '/home/claude/md3')
from cluster import (N, DT, icosahedron13, acc_and_U, H_all, pair_d,
                     PRS, TRI, TRI_P, n_evaporated)

RMIN = 2.0 ** (1.0 / 6.0)
U_MIN = -44.327
SAMP = 25

POS0 = icosahedron13()
D0_PAIR = pair_d(POS0)                      # (78,) static pair distances
D0_TRI = D0_PAIR[TRI_P]                     # (286,3)
H0 = (lambda v: v / (1 + v))(D0_TRI.var(axis=1))   # static H per triple


def H_of(d3):
    v = d3.var(axis=1)
    return v / (1.0 + v)


def H_eps(d3):
    """H computed on strain relative to the static structure."""
    e = (d3 - D0_TRI) / D0_TRI
    v = e.var(axis=1)
    return v / (1.0 + v)


def permutation_drift(pos):
    """How far each atom is from its nearest static site, in units of the
    nearest-neighbour distance.  Large values mean the site mapping is void."""
    c = pos - pos.mean(0)
    dm = np.linalg.norm(c[:, None, :] - POS0[None, :, :], axis=2)
    return float(dm.min(axis=1).mean() / D0_PAIR.min())


def run(T0, seed, nsteps=120000, burn=20000):
    rng = np.random.default_rng(seed)
    pos = icosahedron13() + rng.normal(0, 0.01, (N, 3))
    vel = rng.normal(0, np.sqrt(T0), (N, 3))
    vel -= vel.mean(0)
    acc, _ = acc_and_U(pos)
    rows = []
    for s in range(nsteps):
        vel += 0.5 * DT * acc
        pos += DT * vel
        acc, U = acc_and_U(pos)
        vel += 0.5 * DT * acc
        if s >= burn and s % SAMP == 0:
            d3 = pair_d(pos)[TRI_P]
            h = H_of(d3)
            he = H_eps(d3)
            rows.append((U - U_MIN,
                         h.mean(), (h > 0.22).mean(),       # raw H
                         (h - H0).mean(),                   # baseline-subtracted
                         he.mean(), (he > 0.05).mean(),     # strain H
                         permutation_drift(pos)))
    return np.array(rows)


if __name__ == "__main__":
    print("=== PART A: is sigma / r_min normalisation a no-op here? ===")
    pos = icosahedron13() + np.random.default_rng(1).normal(0, 0.05, (N, 3))
    d3 = pair_d(pos)[TRI_P]
    h_raw = H_of(d3)
    h_sig = H_of(d3 / 1.0)                  # sigma_ij = 1 for all pairs
    h_rmn = H_of(d3 / RMIN)                 # r_min_ij identical for all pairs
    print(f"  distinct sigma_ij among 78 pairs : "
          f"{len(np.unique(np.ones(len(PRS))))}")
    print(f"  distinct r_min_ij among 78 pairs : "
          f"{len(np.unique(np.full(len(PRS), RMIN)))}")
    print(f"  Spearman rho(raw, sigma-normalised) = "
          f"{stats.spearmanr(h_raw, h_sig).statistic:.6f}")
    print(f"  Spearman rho(raw, r_min-normalised) = "
          f"{stats.spearmanr(h_raw, h_rmn).statistic:.6f}")
    print(f"  max |rank change| r_min-normalised  = "
          f"{np.abs(stats.rankdata(h_raw)-stats.rankdata(h_rmn)).max():.0f}")
    print("  -> identical ordering. In a monatomic cluster the correction")
    print("     divides every distance by one constant and cannot act.")

    print(f"\n=== static H baseline: H0 mean {H0.mean():.4f}, "
          f"H_eps baseline exactly {H_eps(D0_TRI).max():.2e} ===")

    print("\n=== PART B: which reference gives the better state variable? ===")
    print("within-temperature Spearman rho against cluster excitation "
          "(U - U_min)\n")
    TEMPS = [0.40, 0.50, 0.65, 0.80]
    NSEED = 10
    names = ['H mean', 'H frac>0.22', 'H - H0 mean', 'H_eps mean',
             'H_eps frac>0.05']
    print(f"{'T0':>6}" + "".join(f"{n:>17}" for n in names)
          + f"{'site drift':>12}")
    out = {}
    for T0 in TEMPS:
        R = np.vstack([run(T0, 6000 + sd) for sd in range(NSEED)])
        y = R[:, 0]
        rhos = [stats.spearmanr(R[:, k], y).statistic for k in range(1, 6)]
        out[T0] = dict(rho=rhos, drift=R[:, 6].mean())
        print(f"{T0:6.2f}" + "".join(f"{r:17.4f}" for r in rhos)
              + f"{R[:,6].mean():12.4f}")

    print("\nsite drift = mean distance of each atom from its nearest static")
    print("site, in nearest-neighbour units. Above ~0.5 the static mapping")
    print("is void and H_eps loses its reference.")
    pickle.dump(out, open('/home/claude/md10/renorm.pkl', 'wb'))
    print("\nsaved")
