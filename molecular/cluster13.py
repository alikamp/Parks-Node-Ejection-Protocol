"""
LJ cluster: local H distribution over all triples, tracked alongside
global state, across temperature from solid to evaporating.

Question: does the SHAPE of the H distribution at time t predict how many
atoms evaporate in [t, t+delta]?

Baseline to beat: the Lindemann parameter (RMS displacement / nearest
neighbour distance), the classical melting indicator.
"""
import numpy as np, pickle, sys, itertools, time

N = 13                 # icosahedral magic number
DT = 0.004
G_CUT = 2.5            # evaporation: no neighbour within this distance

TRI = np.array(list(itertools.combinations(range(N), 3)))      # (286,3)
PRS = np.array(list(itertools.combinations(range(N), 2)))      # (78,2)
# map each triple to its three pair-indices
_pidx = {tuple(sorted(p)): k for k, p in enumerate(PRS)}
TRI_P = np.array([[_pidx[tuple(sorted((a, b)))],
                   _pidx[tuple(sorted((b, c)))],
                   _pidx[tuple(sorted((a, c)))]] for a, b, c in TRI])

def icosahedron13(a=1.1):
    g = (1 + np.sqrt(5)) / 2
    v = []
    for s1 in (-1, 1):
        for s2 in (-1, 1):
            v += [[0, s1, s2 * g], [s1, s2 * g, 0], [s2 * g, 0, s1]]
    v = np.array(v, float)
    v = v / np.linalg.norm(v[0]) * a
    return np.vstack([[0, 0, 0], v])

def acc_and_U(pos):
    d = pos[:, None, :] - pos[None, :, :]
    r2 = (d ** 2).sum(-1); np.fill_diagonal(r2, np.inf)
    inv6 = (1.0 / r2) ** 3
    U = 0.5 * np.sum(4 * (inv6 ** 2 - inv6))
    f = 24 * (2 * inv6 ** 2 - inv6) / r2
    return np.einsum('ij,ijk->ik', f, d), U

def pair_d(pos):
    return np.linalg.norm(pos[PRS[:, 0]] - pos[PRS[:, 1]], axis=1)

def H_all(pos):
    """H for every triple, vectorised. returns (286,)"""
    d = pair_d(pos)[TRI_P]                 # (286,3)
    v = d.var(axis=1)
    return v / (1.0 + v)

def lindemann(ref_r2_sum, ref_r_sum, nsamp, pos):
    """running Lindemann delta from accumulated pair statistics"""
    if nsamp < 2: return np.nan
    m2 = ref_r2_sum / nsamp; m1 = ref_r_sum / nsamp
    var = np.clip(m2 - m1 ** 2, 0, None)
    return float(np.mean(np.sqrt(var) / m1))

def n_evaporated(pos):
    d = pair_d(pos)
    full = np.full((N, N), np.inf)
    full[PRS[:, 0], PRS[:, 1]] = d; full[PRS[:, 1], PRS[:, 0]] = d
    return int(np.sum(full.min(axis=1) > G_CUT))

def run(T0, seed, nsteps=150000, stride=50):
    rng = np.random.default_rng(seed)
    pos = icosahedron13() + rng.normal(0, 0.01, (N, 3))
    vel = rng.normal(0, np.sqrt(T0), (N, 3)); vel -= vel.mean(0)
    acc, _ = acc_and_U(pos)
    r1 = np.zeros(len(PRS)); r2 = np.zeros(len(PRS)); ns = 0
    rec = []
    for s in range(nsteps):
        vel += 0.5 * DT * acc
        pos += DT * vel
        acc, U = acc_and_U(pos)
        vel += 0.5 * DT * acc
        if s % 10 == 0:
            d = pair_d(pos); r1 += d; r2 += d ** 2; ns += 1
        if s % stride == 0:
            h = H_all(pos)
            rec.append(dict(
                t=s * DT,
                h_mean=float(h.mean()), h_sd=float(h.std()),
                h_p90=float(np.percentile(h, 90)),
                h_p99=float(np.percentile(h, 99)),
                h_frac=float(np.mean(h > 0.5)),
                lind=lindemann(r2, r1, ns, pos),
                nev=n_evaporated(pos),
                KE=float(0.5 * (vel ** 2).sum()), U=float(U)))
    return dict(T0=T0, seed=seed, rec=rec)

if __name__ == "__main__":
    out = []; t0 = time.time()
    for T0 in [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40]:
        for sd in range(6):
            out.append(run(T0, sd))
        last = out[-1]['rec'][-1]
        print(f"T0={T0:.2f}  final: lind={last['lind']:.3f} "
              f"H_mean={last['h_mean']:.4f} H_p99={last['h_p99']:.3f} "
              f"evap={last['nev']}  [{time.time()-t0:.0f}s]", flush=True)
    pickle.dump(out, open('/home/claude/md3/cluster.pkl', 'wb'))
    print("saved", len(out), "runs")
