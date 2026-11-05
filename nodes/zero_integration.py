"""Can the first node classify ejection with ZERO integration - and does it
beat Mardling-Aarseth?

compute_cost.py settled the compute question: the analytic first node is ~400x
on its own task but only ~0.8% of the classification window, so it changes
nothing end to end.  Real gain needs one of:

  node-hopping            blocked - Kepler propagation is 422% off and the
                          rotation order is not fixed by the first node
  zero-integration class. the whole ~50,000 force evaluations collapse to 1

This tests the second.  But the benchmark is NOT brute-force integration -
a zero-integration stability criterion for hierarchical triples already
exists: Mardling & Aarseth 2001, which is the ma01_ratio pnep_v12 already
computes from a raw state.  So the question is whether node features add
anything over MA01.

CRITICAL SAMPLING POINT: pnep_v12.make_ic SELECTS systems by MA01 ratio
(stable_min=5.0, unstable_max=0.45).  Any comparison on that sample is rigged
in MA01's favour - it would score near-perfectly by construction.  This draws
initial conditions UNSTRATIFIED: uniform in the orbital elements, random
phase, no MA01 filter, labelled only by integration outcome.

Zero-integration features, all functions of (pos, vel, masses) alone:
  ma01          Mardling-Aarseth ratio (the benchmark)
  g_min/g_max   the normalised radial-velocity triple, sorted
  n_closing     how many pairs are currently closing
  t1_newton     predicted first-node time
  p1_is_inner   whether the predicted initiating pair is the inner binary
  rspread       spread of the three pairwise distances
  H0            H at t=0, for reference
"""
import numpy as np, json, time, sys
from scipy import stats
sys.path.insert(0, '/home/claude/partners')
sys.path.insert(0, '/home/claude/parks-node-ejection-protocol')
from node_partners import accel, hierarchical_ic, PAIRS
from firstnode import radial_velocities, newton_tmin
from firstnode2 import advance
import pnep_v12 as v12

DT = 0.004
T_MAX = 2000.0
EJECT_R = 20.0
N = 160


def integrate_label(pos, vel, masses, t_max=T_MAX, dt=DT, eject_r=EJECT_R):
    pos = pos.copy(); vel = vel.copy()
    n = int(t_max / dt)
    for s in range(n):
        vel += accel(pos, masses) * dt
        pos += vel * dt
        if s % 50 == 0:
            com = np.average(pos, weights=masses, axis=0)
            if np.linalg.norm(pos - com, axis=1).max() > eject_r:
                return True, (s + 1) * dt
    return False, t_max


def features(pos, vel, masses):
    rv = radial_velocities(pos, vel)
    d = np.array([np.linalg.norm(pos[i] - pos[j]) for i, j in PAIRS])
    tn = np.array([newton_tmin(pos, vel, masses, p) for p in range(3)], float)
    n_cl = int(sum(1 for x in rv if x < 0))
    if np.all(np.isnan(tn)):
        p1, t1 = -1, np.nan
    else:
        p1 = int(np.nanargmin(tn)); t1 = float(tn[p1])
    s2 = float(np.var(d))
    srt = np.sort(rv)
    return dict(ma01=float(v12.ma01_ratio(pos, vel, masses)),
                g_min=float(srt[0]), g_mid=float(srt[1]), g_max=float(srt[2]),
                g_absmax=float(np.max(np.abs(rv))),
                n_closing=n_cl, t1_newton=t1,
                p1_is_inner=1.0 if p1 == 0 else 0.0,
                rspread=float(d.max() - d.min()),
                rmin=float(d.min()), rmax=float(d.max()),
                H0=float(s2 / (1 + s2)))


FEATS = ['ma01', 'g_min', 'g_mid', 'g_max', 'g_absmax', 'n_closing',
         't1_newton', 'p1_is_inner', 'rspread', 'rmin', 'rmax', 'H0']
NODE_FEATS = [f for f in FEATS if f != 'ma01']


if __name__ == "__main__":
    rng = np.random.default_rng(1234)
    rows = []
    t0 = time.time()
    P_in = 2 * np.pi * np.sqrt(0.5)

    while len(rows) < N:
        m2 = float(10 ** rng.uniform(-2, 0))          # log-uniform 0.01-1.0
        a_out = float(rng.uniform(2.5, 8.0))
        e_in = float(rng.uniform(0.0, 0.3))
        e_out = float(rng.uniform(0.0, 0.6))
        masses = np.array([1.0, 1.0, m2])
        try:
            pos, vel = hierarchical_ic(masses, 1.0, a_out, e_in, e_out)
            pos, vel = advance(pos, vel, masses,
                               float(rng.uniform(0.0, 3.0) * P_in))
            f = features(pos, vel, masses)
            if not np.isfinite(f['ma01']):
                continue
            ej, te = integrate_label(pos, vel, masses)
        except Exception:
            continue
        rows.append(dict(m2=m2, a_out=a_out, e_in=e_in, e_out=e_out,
                         ejected=bool(ej), t_eject=float(te), **f))
        if len(rows) % 20 == 0:
            nj = sum(1 for r in rows if r['ejected'])
            print(f"  {len(rows)}/{N}  ejected {nj}  [{time.time()-t0:.0f}s]",
                  flush=True)

    json.dump(rows, open('/home/claude/partners/zero_integration.json', 'w'),
              indent=2, default=str)

    y = np.array([r['ejected'] for r in rows], int)
    print(f"\n{'='*74}")
    print(f"n={len(rows)}  ejected {y.sum()} ({100*y.mean():.1f}%)  "
          f"T_max={T_MAX}  [{time.time()-t0:.0f}s]")
    print("NOTE: survivors are right-censored at T_max, not true negatives.")

    print(f"\n--- single-feature AUC (zero integration) ---")
    print(f"{'feature':<16}{'AUC':>8}{'MWU p':>12}")
    aucs = {}
    for f in FEATS:
        x = np.array([r[f] for r in rows], float)
        m = np.isfinite(x)
        if m.sum() < 20 or len(np.unique(x[m])) < 2:
            print(f"{f:<16}{'--':>8}{'--':>12}"); continue
        a, b = x[m][y[m] == 1], x[m][y[m] == 0]
        if len(a) < 5 or len(b) < 5:
            continue
        u = stats.mannwhitneyu(a, b, alternative='two-sided')
        auc = u.statistic / (len(a) * len(b))
        aucs[f] = max(auc, 1 - auc)
        print(f"{f:<16}{max(auc,1-auc):8.3f}{u.pvalue:12.2e}")

    # ---- logistic regression, cross-validated -------------------------
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    from sklearn.impute import SimpleImputer
    from sklearn.model_selection import StratifiedKFold, cross_val_score

    def cv_auc(cols, label):
        X = np.array([[r[c] for c in cols] for r in rows], float)
        pipe = make_pipeline(SimpleImputer(strategy='median'),
                             StandardScaler(),
                             LogisticRegression(max_iter=2000, C=1.0))
        cv = StratifiedKFold(5, shuffle=True, random_state=0)
        s = cross_val_score(pipe, X, y, cv=cv, scoring='roc_auc')
        print(f"  {label:<34}{s.mean():.3f} +/- {s.std():.3f}"
              f"   ({len(cols)} features)")
        return s.mean()

    print(f"\n--- 5-fold cross-validated AUC ---")
    a_ma = cv_auc(['ma01'], 'MA01 alone (the benchmark)')
    a_nd = cv_auc(NODE_FEATS, 'node features only')
    a_bo = cv_auc(FEATS, 'MA01 + node features')
    print(f"\n  node features add {a_bo - a_ma:+.3f} AUC over MA01 alone")
    print(f"  node-only is      {a_nd - a_ma:+.3f} AUC vs MA01 alone")

    print(f"\n--- for reference: the orbital elements themselves ---")
    cv_auc(['a_out', 'e_out', 'e_in', 'm2'], 'raw elements')
    cv_auc(['ma01', 'a_out', 'e_out', 'e_in', 'm2'], 'MA01 + raw elements')

    print(f"\n--- MA01 as a hard threshold (its intended use) ---")
    for thr in (0.8, 1.0, 1.2, 1.5, 2.0):
        pred_stable = np.array([r['ma01'] > thr for r in rows])
        tp = int(((~pred_stable) & (y == 1)).sum())
        fn = int((pred_stable & (y == 1)).sum())
        tn = int((pred_stable & (y == 0)).sum())
        fp = int(((~pred_stable) & (y == 0)).sum())
        acc = (tp + tn) / len(rows)
        rec = tp / max(tp + fn, 1)
        print(f"  ratio>{thr:.1f} = stable: acc {100*acc:5.1f}%  "
              f"ejection recall {100*rec:5.1f}%  "
              f"(TP{tp} FN{fn} TN{tn} FP{fp})")
    print("\nsaved /home/claude/partners/zero_integration.json")
