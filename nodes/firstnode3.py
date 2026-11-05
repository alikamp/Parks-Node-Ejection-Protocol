"""Proper evaluation of the hybrid first-node predictor.

firstnode2.py established, on randomised phases:
  NEWTON (one step of r.v with the full 3-body accel) picks the initiating
  pair 132/135 = 97.8%, but times it badly (median 24% error, p90 301%).
  KEPLER (isolated two-body time to periapsis) picks the pair only 36.5% of
  the time but, when it picks right, times it to 0.17%.

So the perturbation decides WHICH pair; the unperturbed two-body solution
decides WHEN.  Hybrid: Newton for the pair, Kepler on THAT pair for the time.

firstnode2 could only score the hybrid on the 41/135 cases where Kepler's
own minimum happened to fall on Newton's pair, because only the minimum was
stored.  This stores all three per-pair times so the hybrid is evaluated on
every case where Newton makes a prediction.

Also records the closing-pair count, which is known at t=0 and acts as a
confidence tier (1 closing -> 100%, 2 -> 98.4%, 3 -> 86.7%).
"""
import numpy as np, json, itertools, time, sys
from scipy import stats
sys.path.insert(0, '/home/claude/partners')
from node_partners import accel, hierarchical_ic, PAIRS, PNAME
from firstnode import (radial_velocities, newton_tmin, kepler_tmin,
                       pair2b_tmin, true_nodes)
from firstnode2 import advance


if __name__ == "__main__":
    GRID = list(itertools.product(
        [3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.5, 8.0],
        [0.01, 0.1, 0.5],
        [0.2, 0.3, 0.4, 0.5]))
    rng = np.random.default_rng(23)
    rows = []
    t0 = time.time()

    for a_out, m_out, e_out in GRID:
        for rep in range(2):
            masses = np.array([1.0, 1.0, m_out])
            pos, vel = hierarchical_ic(masses, 1.0, a_out, 0.15, e_out)
            P_in = 2 * np.pi * np.sqrt(1.0 / 2.0)
            pos, vel = advance(pos, vel, masses,
                               float(rng.uniform(0.1, 3.0) * P_in))

            rv = radial_velocities(pos, vel)
            tn = [newton_tmin(pos, vel, masses, p) for p in range(3)]
            tk = [kepler_tmin(pos, vel, masses, p) for p in range(3)]
            tb = [pair2b_tmin(pos, vel, masses, p) for p in range(3)]

            nodes = true_nodes(pos, vel, masses)
            if len(nodes) < 6:
                continue
            t_first, p_first = nodes[0]

            tn_a = np.array(tn, float)
            newton_pair = (int(np.nanargmin(tn_a))
                           if not np.all(np.isnan(tn_a)) else None)

            rows.append(dict(
                a_out=a_out, m_out=m_out, e_out=e_out, rep=rep,
                rv=[float(x) for x in rv],
                n_closing=int(sum(1 for x in rv if x < 0)),
                t_first=float(t_first), p_first=int(p_first),
                newton_pair=newton_pair,
                t_newton=[float(x) for x in tn],
                t_kepler=[float(x) for x in tk],
                t_pair2b=[float(x) for x in tb]))

    json.dump(rows, open('/home/claude/partners/firstnode3.json', 'w'),
              indent=2, default=str)
    print(f"{len(rows)} systems  [{time.time()-t0:.0f}s]\n")

    ok = [r for r in rows if r['newton_pair'] is not None]
    hit = sum(1 for r in ok if r['newton_pair'] == r['p_first'])
    print("=" * 76)
    print(f"NEWTON pair accuracy: {hit}/{len(ok)} = {100*hit/len(ok):.1f}%"
          f"   (no prediction for {len(rows)-len(ok)} states)")

    print(f"\nconfidence tier by closing-pair count:")
    for k in (1, 2, 3):
        sub = [r for r in ok if r['n_closing'] == k]
        if sub:
            h = sum(1 for r in sub if r['newton_pair'] == r['p_first'])
            print(f"  {k} closing (n={len(sub):3d}): {h}/{len(sub)} = "
                  f"{100*h/len(sub):5.1f}%")

    print("\n" + "=" * 76)
    print("TIMING - evaluated on all cases where NEWTON predicted")
    for src, lbl in (('t_newton', 'NEWTON own time'),
                     ('t_kepler', 'HYBRID  Newton pair + Kepler time'),
                     ('t_pair2b', 'HYBRID  Newton pair + pair2b time')):
        pr, tr = [], []
        for r in ok:
            p = r['newton_pair']
            v = r[src][p]
            if np.isfinite(v) and r['t_first'] > 0:
                pr.append(v); tr.append(r['t_first'])
        if len(pr) < 5:
            print(f"\n  {lbl}: too few ({len(pr)})"); continue
        pr = np.array(pr); tr = np.array(tr)
        e = np.abs(pr - tr) / tr
        sp = stats.spearmanr(pr, tr)
        print(f"\n  {lbl}   (n={len(pr)})")
        print(f"    median |rel err| {np.median(e):.4f}"
              f"   Spearman {sp.statistic:+.4f} (p={sp.pvalue:.2e})")
        print(f"    percentiles: " + "  ".join(
            f"p{q}={np.percentile(e,q):.4f}" for q in (10, 25, 50, 75, 90)))
        print(f"    within 1%: {(e<0.01).mean():.3f}   "
              f"within 5%: {(e<0.05).mean():.3f}   "
              f"within 10%: {(e<0.10).mean():.3f}")

    # restricted to cases where the pair was right
    print("\n" + "=" * 76)
    print("TIMING - restricted to cases where NEWTON got the pair right")
    right = [r for r in ok if r['newton_pair'] == r['p_first']]
    for src, lbl in (('t_kepler', 'Kepler on the correct pair'),
                     ('t_pair2b', 'pair2b on the correct pair')):
        pr, tr = [], []
        for r in right:
            v = r[src][r['newton_pair']]
            if np.isfinite(v) and r['t_first'] > 0:
                pr.append(v); tr.append(r['t_first'])
        if len(pr) < 5:
            continue
        pr = np.array(pr); tr = np.array(tr)
        e = np.abs(pr - tr) / tr
        sp = stats.spearmanr(pr, tr)
        print(f"\n  {lbl}   (n={len(pr)})")
        print(f"    median |rel err| {np.median(e):.4f}"
              f"   Spearman {sp.statistic:+.4f}")
        print(f"    within 1%: {(e<0.01).mean():.3f}   "
              f"within 5%: {(e<0.05).mean():.3f}")
    print("\nsaved /home/claude/partners/firstnode3.json")
