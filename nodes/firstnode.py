"""Can the FIRST mirror-symmetry node be predicted from initial conditions?

Hypothesis (Alika): from an arbitrary n-body initial state you want to know
WHEN and WHICH PAIR initiates the nodal sequence, and the first node is what
unravels the system state.

The node condition is analytic.  A local minimum of |r_ij| means

    d|r_ij|^2/dt = 2 r_ij . v_ij = 0,   with d^2/dt^2 > 0

so "which pair fires first and when" is root-finding on three scalars of the
initial state, not a search through an integrated trajectory.  The pair's
radial velocity rdot_ij = (r_ij . v_ij)/|r_ij| is negative while closing.

(This is the condition pnep_v15_4_full_engine.py already uses, restricted to
the outer body against the inner binary CoM.  Generalising it to all three
pairs reconciles v15.4 with the README's local-minimum definition.)

Three predictors of the first node, cheapest first:

  NEWTON   linear extrapolation of r.v to zero using
           d(r.v)/dt = |v|^2 + r.a   (a = full 3-body acceleration)
  KEPLER   treat the pair as an isolated two-body orbit about mu = G(m_i+m_j)
           and solve Kepler's equation for time to periapsis
  PAIR2B   numerically integrate the isolated two-body pair to its periapsis

Each is scored against the true first node from a full 3-body integration:
did it name the right pair, and how close was the time?

Also tests whether the first node fixes the ROTATION ORDER of the subsequent
round-robin, which is what would make the whole sequence predictable.
"""
import numpy as np, json, itertools, time, sys
from scipy import stats
sys.path.insert(0, '/home/claude/partners')
from node_partners import accel, hierarchical_ic, PAIRS, PNAME

G = 1.0


def rel(pos, vel, p):
    i, j = PAIRS[p]
    return pos[i] - pos[j], vel[i] - vel[j]


def radial_velocities(pos, vel):
    """rdot for each pair; negative = closing."""
    out = np.zeros(3)
    for p in range(3):
        r, v = rel(pos, vel, p)
        out[p] = np.dot(r, v) / np.linalg.norm(r)
    return out


def newton_tmin(pos, vel, masses, p):
    """Linear extrapolation of r.v to zero using the full 3-body accel."""
    r, v = rel(pos, vel, p)
    a = accel(pos, masses)
    i, j = PAIRS[p]
    arel = a[i] - a[j]
    rv = np.dot(r, v)
    drv = np.dot(v, v) + np.dot(r, arel)
    if drv <= 0:
        return np.nan                 # not heading toward a minimum
    if rv >= 0:
        return np.nan                 # already separating
    return -rv / drv


def kepler_tmin(pos, vel, masses, p):
    """Time to periapsis of the isolated two-body orbit for this pair."""
    r, v = rel(pos, vel, p)
    i, j = PAIRS[p]
    mu = G * (masses[i] + masses[j])
    rn = np.linalg.norm(r)
    vn2 = np.dot(v, v)
    E = 0.5 * vn2 - mu / rn
    if E >= 0:
        return np.nan                 # unbound pair
    a = -mu / (2 * E)
    h = np.cross(r, v)
    hn2 = np.dot(h, h)
    e2 = 1 + 2 * E * hn2 / mu ** 2
    e = np.sqrt(max(e2, 0.0))
    if e >= 1:
        return np.nan
    # eccentric anomaly from r and r.v
    cosE = (1 - rn / a) / e if e > 1e-12 else 0.0
    cosE = np.clip(cosE, -1.0, 1.0)
    Ean = np.arccos(cosE)
    if np.dot(r, v) > 0:              # past periapsis, heading out
        pass
    else:
        Ean = 2 * np.pi - Ean         # before periapsis
    M = Ean - e * np.sin(Ean)
    n = np.sqrt(mu / a ** 3)
    # time until M reaches 2*pi (next periapsis)
    return (2 * np.pi - M) / n


def pair2b_tmin(pos, vel, masses, p, dt=1e-3, max_t=200.0):
    """Integrate the isolated two-body pair to its next periapsis."""
    r, v = rel(pos, vel, p)
    i, j = PAIRS[p]
    mu = G * (masses[i] + masses[j])
    rv_prev = np.dot(r, v)
    t = 0.0
    while t < max_t:
        a = -mu * r / np.linalg.norm(r) ** 3
        v = v + a * dt
        r = r + v * dt
        t += dt
        rv = np.dot(r, v)
        if rv_prev < 0 <= rv:
            return t
        rv_prev = rv
    return np.nan


def true_nodes(pos, vel, masses, dt=0.0005, steps=400000, n_want=40):
    """Full 3-body integration; return the first n_want (t, pair)."""
    pos = pos.copy(); vel = vel.copy()

    def pd(p):
        return np.array([np.linalg.norm(p[a] - p[b]) for a, b in PAIRS])

    d2 = pd(pos)
    d1 = pd(pos + vel * (dt / 2))
    t = 0.0
    out = []
    for _ in range(steps):
        vel += accel(pos, masses) * dt
        pos += vel * dt
        t += dt
        d0 = pd(pos)
        for p in range(3):
            if d1[p] < d2[p] and d1[p] < d0[p]:
                out.append((t - dt, p))
        if len(out) >= n_want:
            break
        d2, d1 = d1, d0
    return out


def rotation_order(seq):
    """Dominant cyclic order of the partner sequence, as a tuple, or None."""
    trans = {}
    for a, b in zip(seq[:-1], seq[1:]):
        if a != b:
            trans[a] = trans.get(a, {})
            trans[a][b] = trans[a].get(b, 0) + 1
    if len(trans) < 3:
        return None
    nxt = {a: max(d, key=d.get) for a, d in trans.items()}
    # follow from 0
    order = [0]
    for _ in range(3):
        n = nxt.get(order[-1])
        if n is None or n in order:
            break
        order.append(n)
    return tuple(order) if len(order) == 3 else None


if __name__ == "__main__":
    GRID = list(itertools.product(
        [3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.5, 8.0],     # a_out
        [0.01, 0.1, 0.5],                              # m_out
        [0.2, 0.3, 0.4, 0.5],                          # e_out
        [0.05, 0.15]))                                 # e_in
    rows = []
    t0 = time.time()

    for a_out, m_out, e_out, e_in in GRID:
        masses = np.array([1.0, 1.0, m_out])
        pos, vel = hierarchical_ic(masses, 1.0, a_out, e_in, e_out)

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

        rows.append(dict(
            a_out=a_out, m_out=m_out, e_out=e_out, e_in=e_in,
            rv=[float(x) for x in rv],
            t_first=float(t_first), p_first=int(p_first),
            newton_pair=kn, newton_t=vn,
            kepler_pair=kk, kepler_t=vk,
            pair2b_pair=kp, pair2b_t=vp,
            order=str(order), n_nodes=len(nodes),
            seq_head=[int(x) for x in seq[:12]]))

    json.dump(rows, open('/home/claude/partners/firstnode.json', 'w'),
              indent=2, default=str)

    print(f"{len(rows)} systems  [{time.time()-t0:.0f}s]\n")
    print("=" * 78)
    print("PART 1 - predicting the first node from initial conditions")
    for key, lbl in (('newton', 'NEWTON  (linear r.v extrapolation)'),
                     ('kepler', 'KEPLER  (two-body time to periapsis)'),
                     ('pair2b', 'PAIR2B  (integrate isolated pair)')):
        ok = [r for r in rows if r[f'{key}_pair'] is not None]
        if not ok:
            print(f"  {lbl}: no predictions")
            continue
        hit = sum(1 for r in ok if r[f'{key}_pair'] == r['p_first'])
        tt = np.array([r['t_first'] for r in ok])
        tp_ = np.array([r[f'{key}_t'] for r in ok])
        m = np.isfinite(tp_)
        relerr = np.abs(tp_[m] - tt[m]) / tt[m]
        sp = stats.spearmanr(tp_[m], tt[m])
        print(f"\n  {lbl}")
        print(f"    pair correct     : {hit}/{len(ok)}  ({100*hit/len(ok):.1f}%)"
              f"   [chance = 33.3%]")
        print(f"    median |rel err| in t: {np.median(relerr):.4f}")
        print(f"    Spearman(pred t, true t) = {sp.statistic:+.4f}"
              f"  p={sp.pvalue:.2e}")

    print("\n" + "=" * 78)
    print("PART 2 - does the first node fix the rotation order?")
    from collections import Counter
    c = Counter(r['order'] for r in rows)
    for o, n in c.most_common():
        print(f"  order {o}: {n} systems")
    # is order determined by which pair is first?
    print(f"\n  rotation order by first pair:")
    for p in range(3):
        sub = [r['order'] for r in rows if r['p_first'] == p]
        if sub:
            cc = Counter(sub)
            top, nn = cc.most_common(1)[0]
            print(f"    first={PNAME[p]} (n={len(sub)}): most common {top}"
                  f"  {nn}/{len(sub)} = {100*nn/len(sub):.0f}%")

    print(f"\n  does rv sign pattern at t=0 determine the first pair?")
    pat = Counter()
    for r in rows:
        sgn = tuple(int(np.sign(x)) for x in r['rv'])
        pat[(sgn, r['p_first'])] += 1
    seen = {}
    for (sgn, pf), n in pat.items():
        seen.setdefault(sgn, []).append((PNAME[pf], n))
    for sgn, lst in sorted(seen.items()):
        tot = sum(n for _, n in lst)
        best = max(lst, key=lambda x: x[1])
        print(f"    rv signs {sgn}: n={tot}, first pair mostly {best[0]}"
              f"  {best[1]}/{tot} = {100*best[1]/tot:.0f}%")
    print("\nsaved /home/claude/partners/firstnode.json")
