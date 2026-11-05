"""Partner structure and internodal timing in the FLASH ejection regime.

The earlier attempt used the fixed-step leapfrog and produced 3 nodes in
1600 t.u. - invalid.  Burrau starts at rest and has very close approaches, so
it needs adaptive high-order integration; verify_burrau.py uses DOP853 at
rtol=1e-12 for exactly this reason.  This repeats the partner analysis with
that integrator.

Flash is the regime with the worst lead time in the framework (2-10t against
412.6t for slow decay) and the one place internodal timing has never been
characterised.

Tests canonical Burrau (masses 3:4:5, the repo's benchmark) plus several
near-equal-mass compact triples.
"""
import numpy as np, json
from scipy.integrate import solve_ivp
from scipy import stats

PAIRS = [(0, 1), (1, 2), (0, 2)]
PNAME = {0: "01", 1: "12", 2: "02"}
FWD = {(0, 1), (1, 2), (2, 0)}


def rhs(t, y, m):
    p = y[:9].reshape(3, 3)
    v = y[9:].reshape(3, 3)
    d = p[None, :, :] - p[:, None, :]
    r2 = (d ** 2).sum(-1)
    np.fill_diagonal(r2, np.inf)
    a = np.einsum('ij,ijk->ik', m[None, :] * r2 ** -1.5, d)
    return np.concatenate([v.ravel(), a.ravel()])


def energy(p, v, m):
    T = 0.5 * (m[:, None] * v ** 2).sum()
    U = 0.0
    for i, j in PAIRS:
        U -= m[i] * m[j] / np.linalg.norm(p[i] - p[j])
    return T + U


def integrate(pos, vel, m, t_end=100.0, step=0.002):
    y0 = np.concatenate([pos.ravel(), vel.ravel()])
    t_eval = np.arange(0, t_end, step)
    sol = solve_ivp(rhs, (0, t_end), y0, method='DOP853', t_eval=t_eval,
                    args=(m,), rtol=1e-12, atol=1e-14)
    P = sol.y[:9].reshape(3, 3, -1).transpose(2, 0, 1)
    V = sol.y[9:].reshape(3, 3, -1).transpose(2, 0, 1)
    return sol.t, P, V


def nodes_with_partner(t, P):
    """Local minimum of any pairwise distance; record which pair, and H."""
    D = np.stack([np.linalg.norm(P[:, j] - P[:, i], axis=1)
                  for i, j in PAIRS], axis=1)          # (T,3)
    T_, Pr_, H_ = [], [], []
    for k in range(2, len(t)):
        hit = [p for p in range(3)
               if D[k - 1, p] < D[k - 2, p] and D[k - 1, p] < D[k, p]]
        if hit:
            s2 = np.var(D[k - 1])
            H = s2 / (1 + s2)
            for p in hit:
                T_.append(t[k - 1]); Pr_.append(p); H_.append(H)
    return np.array(T_), np.array(Pr_), np.array(H_)


def eject_time(t, P, m, radius=20.0):
    """Repo criterion: verify_burrau.py tests R > 12, 20, 50 from the CoM.
    An earlier version of this script used R > 80, which reported no ejection
    for Burrau (known to eject near t=63) and made four 'no ejection' labels
    meaningless. 20 is the middle of the repo's set."""
    com = np.average(P, weights=m, axis=1)
    R = np.linalg.norm(P - com[:, None, :], axis=2)
    for k in range(len(t)):
        if R[k].max() > radius:
            return t[k], int(R[k].argmax())
    return None, None


def eject_all_thresholds(t, P, m):
    com = np.average(P, weights=m, axis=1)
    R = np.linalg.norm(P - com[:, None, :], axis=2)
    out = {}
    for thr in (12.0, 20.0, 50.0):
        idx = np.where((R > thr).any(axis=1))[0]
        out[thr] = (float(t[idx[0]]), int(R[idx[0]].argmax())) if len(idx) else None
    return out


def h_spike(node_t, node_H, thresh=0.5):
    """First node where H crosses thresh AND stays above for the rest -
    the repo's 'sustained spike' definition, plus the first bare crossing."""
    sustained = None
    for i in range(len(node_H)):
        if node_H[i] > thresh and np.all(node_H[i:] > thresh):
            sustained = float(node_t[i])
            break
    first = float(node_t[np.argmax(node_H > thresh)]) \
        if np.any(node_H > thresh) else None
    return first, sustained


def chirality(prim):
    F = B = 0
    for a, b in zip(prim[:-1], prim[1:]):
        if a == b:
            continue
        if (a, b) in FWD:
            F += 1
        else:
            B += 1
    return (F - B) / max(F + B, 1), F, B


def report(name, t, P, V, m, t_end):
    E0 = energy(P[0], V[0], m)
    E1 = energy(P[-1], V[-1], m)
    drift = abs((E1 - E0) / E0) * 100
    te, be = eject_time(t, P, m)
    thr_all = eject_all_thresholds(t, P, m)
    T_, Pr_, H_ = nodes_with_partner(t, P)
    first_x, sust_x = h_spike(T_, H_)

    print(f"\n{'='*72}\n### {name}")
    print(f"  energy drift {drift:.5f}%   integrated to t={t[-1]:.1f}")
    print("  ejection by threshold: " + "  ".join(
        f"R>{k:.0f}: " + (f"t={v[0]:.2f} b{v[1]}" if v else "none")
        for k, v in thr_all.items()))
    print(f"  H>0.5 first crossing t={first_x}   sustained spike t={sust_x}")
    if sust_x is not None and te is not None:
        print(f"  lead time from sustained spike: {te - sust_x:.2f}t")
    if len(T_) < 8:
        print(f"  only {len(T_)} nodes - cannot analyse")
        return dict(system=name, n_nodes=int(len(T_)), drift=float(drift))

    # restrict to pre-ejection
    if te:
        keep = T_ < te
        T_, Pr_, H_ = T_[keep], Pr_[keep], H_[keep]

    dtn = np.diff(T_)
    pooled = dtn.std() / dtn.mean() if len(dtn) > 1 else np.nan
    sw = float((Pr_[1:] != Pr_[:-1]).mean()) if len(Pr_) > 1 else np.nan
    chi, F, B = chirality(Pr_)

    print(f"  nodes before ejection: {len(T_)}")
    print(f"  H range: {H_.min():.4f} - {H_.max():.4f}   "
          f"(scramble zone per README: 0.05-0.30, spike >0.50)")
    print(f"  ownership: " + "  ".join(
        f"{PNAME[p]} {100*(Pr_==p).mean():5.1f}%" for p in range(3)))
    print(f"  switch rate {sw:.4f}   chirality {chi:+.4f} (F={F} B={B})")
    print(f"  mean dt_node {dtn.mean():.4f}   pooled jitter {pooled:.4f}")

    own = [float((Pr_ == p).mean()) for p in range(3)]
    out = dict(system=name, n_nodes=int(len(T_)), drift=float(drift),
               t_eject=te, body_eject=be, switch_rate=sw, chi=float(chi),
               pooled_jitter=float(pooled), mean_dt=float(dtn.mean()),
               H_min=float(H_.min()), H_max=float(H_.max()),
               ownership=own, own_spread=float(max(own) - min(own)),
               h_first_cross=first_x, h_sustained=sust_x,
               lead_time=(te - sust_x) if (te and sust_x) else None,
               eject_thresholds={str(k): v for k, v in thr_all.items()})

    # per-partner intervals
    print(f"  per-partner intervals:")
    for p in range(3):
        tp = T_[Pr_ == p]
        if len(tp) > 3:
            d = np.diff(tp)
            print(f"    {PNAME[p]}: n={len(d):4d} mean {d.mean():7.4f}"
                  f" jitter {d.std()/d.mean():.4f}")
            out[f'ppj_{PNAME[p]}'] = float(d.std() / d.mean())

    # does dt_node trend before the spike?
    if te and len(dtn) > 10:
        half = len(dtn) // 2
        print(f"  dt_node first half {dtn[:half].mean():.4f} "
              f"-> second half {dtn[half:].mean():.4f}")
        sl = stats.linregress(np.arange(len(dtn)), dtn)
        print(f"  dt_node trend: slope {sl.slope:+.6f} per node, "
              f"p={sl.pvalue:.4f}, R2={sl.rvalue**2:.4f}")
        out['dt_slope'] = float(sl.slope)
        out['dt_slope_p'] = float(sl.pvalue)
    return out


if __name__ == "__main__":
    results = []

    # canonical Burrau
    m = np.array([3.0, 4.0, 5.0])
    pos = np.array([[1.0, 3.0, 0.0], [-2.0, -1.0, 0.0], [1.0, -1.0, 0.0]])
    vel = np.zeros((3, 3))
    pos -= (m[:, None] * pos).sum(0) / m.sum()
    t, P, V = integrate(pos, vel, m, t_end=100.0)
    results.append(report("Burrau (3:4:5, at rest)", t, P, V, m, 100.0))

    # near-equal-mass compact triples, small random velocities
    rng = np.random.default_rng(7)
    for k in range(6):
        m = np.array([1.0, 1.0, 1.0])
        pos = rng.normal(0, 1.2, (3, 3)); pos[:, 2] = 0.0
        vel = rng.normal(0, 0.25, (3, 3)); vel[:, 2] = 0.0
        pos -= (m[:, None] * pos).sum(0) / m.sum()
        vel -= (m[:, None] * vel).sum(0) / m.sum()
        try:
            t, P, V = integrate(pos, vel, m, t_end=120.0)
            results.append(report(f"equal-mass scramble #{k+1}", t, P, V, m, 120.0))
        except Exception as e:
            print(f"\n### equal-mass scramble #{k+1}: integration failed: {e}")

    json.dump(results, open('/home/claude/partners/flash.json', 'w'),
              indent=2, default=str)

    good = [r for r in results if r.get('n_nodes', 0) >= 8]
    if len(good) > 2:
        print(f"\n{'='*72}\n### FLASH REGIME SUMMARY ({len(good)} systems)")
        for key in ('switch_rate', 'chi', 'pooled_jitter'):
            v = np.array([r[key] for r in good if r.get(key) is not None])
            print(f"  {key:16s} mean {v.mean():+.4f}  sd {v.std():.4f}  "
                  f"range {v.min():+.4f} to {v.max():+.4f}")
    print("\nsaved /home/claude/partners/flash.json")
