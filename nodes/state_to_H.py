"""Is there a closed-form map from the Newtonian state to H at the next node?

Half the chain exists and is validated:
    g_P = r_P . v_P                          (nodal function)
    gdot_P = |v_P|^2 - mu_P/|r_P| + r_P.T_P  (exact derivative)
    t*_P ~ -g_P/gdot_P                       (Newton step)
-> which pair initiates (98.5%) and when (rho=0.955), from the state alone.

The other half does not: H = sigma^2/(1+sigma^2) needs ALL THREE distances at
the node time, and propagating three bodies to a future time IS the three-body
problem.  So the chain breaks exactly where it would have paid off, which is
why the compute saving stalled at 0.63% of the pipeline.

But each distance has an exact second derivative at the current instant:
    ddot     = g/d
    dddot    = gdot/d - g^2/d^3
so a second-order expansion gives a CLOSED-FORM approximation of each distance
at the predicted node time, and H follows:
    d_Q(t*) ~ d_Q + (g_Q/d_Q) t* + 0.5 (gdot_Q/d_Q - g_Q^2/d_Q^3) t*^2

Tested here against the truth from integration, with the error decomposed:

  H_now      H at t=0 - the do-nothing baseline.  Any estimate must beat this
             or it is adding nothing.
  H_taylor   the formula above at the Newton-predicted t*  (timing error AND
             extrapolation error together)
  H_exact_t  Taylor expansion evaluated at the TRUE node time (extrapolation
             error alone - isolates how much of the failure is bad timing)
  H_true     from integration at the true node
"""
import numpy as np, json, itertools, time, sys
from scipy import stats
sys.path.insert(0, '/home/claude/partners')
from node_partners import accel, hierarchical_ic, PAIRS
from firstnode2 import advance

DT = 5e-4


def g_and_gdot(pos, vel, masses):
    """Exact g_P and gdot_P for all three pairs, plus distances."""
    a = accel(pos, masses)
    d = np.zeros(3); g = np.zeros(3); gd = np.zeros(3)
    for p, (i, j) in enumerate(PAIRS):
        r = pos[i] - pos[j]; v = vel[i] - vel[j]; ar = a[i] - a[j]
        d[p] = np.linalg.norm(r)
        g[p] = np.dot(r, v)
        gd[p] = np.dot(v, v) + np.dot(r, ar)
    return d, g, gd


def H_of(d):
    s2 = float(np.var(d))
    return s2 / (1.0 + s2)


def taylor_d(d, g, gd, t):
    """Second-order expansion of each distance to time t."""
    ddot = g / d
    dddot = gd / d - g ** 2 / d ** 3
    return d + ddot * t + 0.5 * dddot * t ** 2


def true_next_node(pos, vel, masses, dt=DT, max_steps=400000):
    """Integrate to the first node; return (t, pair, distances there)."""
    pos = pos.copy(); vel = vel.copy()

    def pd(p):
        return np.array([np.linalg.norm(p[a] - p[b]) for a, b in PAIRS])

    d2 = pd(pos); d1 = pd(pos + vel * (dt / 2))
    t = 0.0
    for s in range(max_steps):
        vel += accel(pos, masses) * dt
        pos += vel * dt
        t += dt
        d0 = pd(pos)
        for p in range(3):
            if d1[p] < d2[p] and d1[p] < d0[p]:
                return t - dt, p, d1
        d2, d1 = d1, d0
    return None, None, None


if __name__ == "__main__":
    GRID = list(itertools.product([3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.5, 8.0],
                                  [0.01, 0.1, 0.5],
                                  [0.2, 0.3, 0.4, 0.5]))
    rng = np.random.default_rng(99)
    P_in = 2 * np.pi * np.sqrt(0.5)
    rows = []
    t0 = time.time()

    for a_out, m_out, e_out in GRID:
        for rep in range(2):
            masses = np.array([1.0, 1.0, m_out])
            pos, vel = hierarchical_ic(masses, 1.0, a_out, 0.15, e_out)
            pos, vel = advance(pos, vel, masses,
                               float(rng.uniform(0.0, 3.0) * P_in))

            d0, g0, gd0 = g_and_gdot(pos, vel, masses)
            # Newton prediction
            tn = np.full(3, np.nan)
            for p in range(3):
                if g0[p] < 0 and gd0[p] > 0:
                    tn[p] = -g0[p] / gd0[p]
            if np.all(np.isnan(tn)):
                continue
            p_pred = int(np.nanargmin(tn)); t_pred = float(tn[p_pred])

            t_true, p_true, d_true = true_next_node(pos, vel, masses)
            if t_true is None or t_true <= 0:
                continue

            rows.append(dict(
                a_out=a_out, m_out=m_out, e_out=e_out, rep=rep,
                t_pred=t_pred, t_true=float(t_true),
                p_pred=p_pred, p_true=int(p_true),
                H_now=H_of(d0),
                H_taylor=H_of(taylor_d(d0, g0, gd0, t_pred)),
                H_exact_t=H_of(taylor_d(d0, g0, gd0, t_true)),
                H_true=H_of(d_true),
                d_now=[float(x) for x in d0],
                d_true=[float(x) for x in d_true],
                d_taylor=[float(x) for x in taylor_d(d0, g0, gd0, t_true)]))

    json.dump(rows, open('/home/claude/partners/state_to_H.json', 'w'),
              indent=2, default=str)
    print(f"{len(rows)} states  [{time.time()-t0:.0f}s]\n")

    ht = np.array([r['H_true'] for r in rows])
    print("=" * 74)
    print(f"H at the next node: can it be computed from the state?")
    print(f"  H_true range {ht.min():.4f} - {ht.max():.4f}, "
          f"mean {ht.mean():.4f}, sd {ht.std():.4f}")
    print(f"\n{'estimator':<30}{'MAE':>9}{'median AE':>11}{'rho':>8}{'p':>10}")
    for key, lbl in (('H_now', 'H at t=0 (do-nothing)'),
                     ('H_taylor', 'Taylor at predicted t*'),
                     ('H_exact_t', 'Taylor at TRUE node time')):
        x = np.array([r[key] for r in rows])
        m = np.isfinite(x)
        e = np.abs(x[m] - ht[m])
        sp = stats.spearmanr(x[m], ht[m])
        print(f"  {lbl:<28}{e.mean():9.4f}{np.median(e):11.4f}"
              f"{sp.statistic:+8.3f}{sp.pvalue:10.2e}")

    print(f"\n  (an estimator must beat the do-nothing row to add anything)")

    # where does the error come from?
    print(f"\n--- error decomposition ---")
    tp = np.array([r['t_pred'] for r in rows])
    tt = np.array([r['t_true'] for r in rows])
    print(f"  timing: median |rel err| in t* = "
          f"{np.median(np.abs(tp-tt)/tt):.4f}")
    pa = np.mean([r['p_pred'] == r['p_true'] for r in rows])
    print(f"  pair identity correct: {100*pa:.1f}%")
    # per-distance extrapolation error at the true time
    dt_err = []
    for r in rows:
        a = np.array(r['d_taylor']); b = np.array(r['d_true'])
        dt_err.append(np.abs(a - b) / b)
    dt_err = np.array(dt_err)
    print(f"  distance extrapolation at the TRUE time, median |rel err|:")
    for p, nm in enumerate(('01', '12', '02')):
        print(f"    pair {nm}: {np.median(dt_err[:, p]):.4f}")
    print(f"  worst of the three, median: "
          f"{np.median(dt_err.max(axis=1)):.4f}")

    print(f"\n--- does the Taylor H at least RANK systems correctly? ---")
    x = np.array([r['H_taylor'] for r in rows])
    m = np.isfinite(x)
    sp = stats.spearmanr(x[m], ht[m])
    print(f"  Spearman(H_taylor, H_true) = {sp.statistic:+.4f} "
          f"p={sp.pvalue:.2e}")
    x0 = np.array([r['H_now'] for r in rows])
    sp0 = stats.spearmanr(x0, ht)
    print(f"  Spearman(H_now,    H_true) = {sp0.statistic:+.4f} "
          f"p={sp0.pvalue:.2e}")
    print(f"\n  if H_now already ranks as well, the formula adds nothing and")
    print(f"  the state->H map is not worth having.")
