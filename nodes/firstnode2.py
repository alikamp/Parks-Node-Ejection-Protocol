"""First-node prediction from an ARBITRARY starting state.

firstnode.py had a degenerate test design: hierarchical_ic() places the inner
binary exactly at periapsis, so pair (0,1) has r.v = 0 identically at t=0.
Kepler then correctly returns "pair 01, t=0" for every system - true but
useless, since the pair is already AT a node and the detector cannot fire
there (it needs a bracketing sample from before t=0).  Pair accuracy is
understated for a reason baked into the initial conditions.

Fix, and it is closer to the hypothesis as stated ("take any n-body starting
position"): build the configuration, advance it by a RANDOM offset, and treat
that state as the arbitrary initial condition.  A system started exactly at
periapsis is a measure-zero case.

Scored on the same three predictors:
  NEWTON  linear extrapolation of r.v to zero with the full 3-body accel
  KEPLER  isolated two-body time to periapsis
  PAIR2B  numerically integrate the isolated pair to periapsis

plus whether the first node fixes the subsequent rotation order, and whether
the t=0 radial-velocity triple alone predicts which pair goes first.
"""
import numpy as np, json, itertools, time, sys
from collections import Counter
from scipy import stats
sys.path.insert(0, '/home/claude/partners')
from node_partners import accel, hierarchical_ic, PAIRS, PNAME
from firstnode import (radial_velocities, newton_tmin, kepler_tmin,
                       pair2b_tmin, true_nodes, rotation_order)


def advance(pos, vel, masses, T, dt=0.001):
    """Integrate forward by T to randomise the phase."""
    pos = pos.copy(); vel = vel.copy()
    n = int(T / dt)
    for _ in range(n):
        vel += accel(pos, masses) * dt
        pos += vel * dt
    return pos, vel


if __name__ == "__main__":
    GRID = list(itertools.product(
        [3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.5, 8.0],
        [0.01, 0.1, 0.5],
        [0.2, 0.3, 0.4, 0.5]))
    rng = np.random.default_rng(11)
    rows = []
    t0 = time.time()

    for a_out, m_out, e_out in GRID:
        for rep in range(2):                       # 2 random phases each
            masses = np.array([1.0, 1.0, m_out])
            pos, vel = hierarchical_ic(masses, 1.0, a_out, 0.15, e_out)
            P_in = 2 * np.pi * np.sqrt(1.0 / 2.0)
            pos, vel = advance(pos, vel, masses,
                               float(rng.uniform(0.1, 3.0) * P_in))

            rv = radial_velocities(pos, vel)
            tn = np.array([newton_tmin(pos, vel, masses, p) for p in range(3)])
            tk = np.array([kepler_tmin(pos, vel, masses, p) for p in range(3)])
            tp = np.array([pair2b_tmin(pos, vel, masses, p) for p in range(3)])

            nodes = true_nodes(pos, vel, masses)
            if len(nodes) < 10:
                continue
            t_first, p_first = nodes[0]
            seq = [p for _, p in nodes]
            order = rotation_order(seq)

            def pick(arr):
                if np.all(np.isnan(arr)):
                    return None, np.nan
                k = int(np.nanargmin(arr))
                return k, float(arr[k])

            kn, vn = pick(tn)
            kk, vk = pick(tk)
            kp, vp = pick(tp)
            # closing pair with the LEAST negative r.v is nearest its minimum
            closing = [p for p in range(3) if rv[p] < 0]
            rv_pick = max(closing, key=lambda p: rv[p]) if closing else None

            rows.append(dict(
                a_out=a_out, m_out=m_out, e_out=e_out, rep=rep,
                rv=[float(x) for x in rv], n_closing=len(closing),
                t_first=float(t_first), p_first=int(p_first),
                newton_pair=kn, newton_t=vn,
                kepler_pair=kk, kepler_t=vk,
                pair2b_pair=kp, pair2b_t=vp,
                rv_pair=rv_pick,
                order=str(order), seq_head=[int(x) for x in seq[:12]]))

    json.dump(rows, open('/home/claude/partners/firstnode2.json', 'w'),
              indent=2, default=str)
    print(f"{len(rows)} systems (randomised phase)  [{time.time()-t0:.0f}s]\n")

    print("=" * 78)
    print("PART 1 - predict WHICH pair fires first, and WHEN, from t=0 state")
    for key, lbl in (('newton', 'NEWTON  linear r.v extrapolation'),
                     ('kepler', 'KEPLER  two-body time to periapsis'),
                     ('pair2b', 'PAIR2B  integrate isolated pair')):
        ok = [r for r in rows if r[f'{key}_pair'] is not None]
        if not ok:
            print(f"\n  {lbl}: no predictions"); continue
        hit = sum(1 for r in ok if r[f'{key}_pair'] == r['p_first'])
        tt = np.array([r['t_first'] for r in ok])
        pr = np.array([r[f'{key}_t'] for r in ok])
        m = np.isfinite(pr) & (tt > 0)
        print(f"\n  {lbl}")
        print(f"    pair correct : {hit}/{len(ok)} = {100*hit/len(ok):5.1f}%"
              f"   [chance 33.3%]")
        if m.sum() > 3:
            relerr = np.abs(pr[m] - tt[m]) / tt[m]
            sp = stats.spearmanr(pr[m], tt[m])
            print(f"    timing       : median |rel err| {np.median(relerr):.4f}"
                  f"   Spearman {sp.statistic:+.4f} (p={sp.pvalue:.2e})")
            # accuracy restricted to cases where the pair was right
            mr = np.array([r[f'{key}_pair'] == r['p_first'] for r in ok]) & m
            if mr.sum() > 3:
                re2 = np.abs(pr[mr] - tt[mr]) / tt[mr]
                print(f"    timing | pair correct: median |rel err| "
                      f"{np.median(re2):.4f}  (n={mr.sum()})")

    ok = [r for r in rows if r['rv_pair'] is not None]
    if ok:
        hit = sum(1 for r in ok if r['rv_pair'] == r['p_first'])
        print(f"\n  RV-SIGN  least-negative closing pair")
        print(f"    pair correct : {hit}/{len(ok)} = {100*hit/len(ok):5.1f}%")

    print(f"\n  how many pairs are closing at t=0?")
    for k, n in sorted(Counter(r['n_closing'] for r in rows).items()):
        print(f"    {k} closing: {n} systems")

    print("\n" + "=" * 78)
    print("PART 2 - does the first node fix the rotation order?")
    c = Counter(r['order'] for r in rows)
    for o, n in c.most_common():
        print(f"  order {o}: {n}")
    print(f"\n  rotation order given the first pair:")
    for p in range(3):
        sub = [r['order'] for r in rows if r['p_first'] == p]
        if sub:
            cc = Counter(sub).most_common(1)[0]
            print(f"    first={PNAME[p]} (n={len(sub)}): {cc[0]} in "
                  f"{cc[1]}/{len(sub)} = {100*cc[1]/len(sub):.0f}%")
    print("\nsaved /home/claude/partners/firstnode2.json")
