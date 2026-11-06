"""Do the identity channels beat the mass prior where the prior actually fails?

Section 11 found "name the lightest body" scores 96.5% on the usual grid and
100% wherever m3 < 0.4, so a 100% identity result there demonstrates nothing.
Identity is only non-trivial near equal masses.

This draws systems where the prior has real work to do: all three masses
varied, so the lightest body is not systematically the outer one.

Channels tested, each against the mass prior:

  MASS      name the lightest body                       (the null model)
  VEL       at the highest-H node, the body with the largest speed relative
            to the centre of mass                        (README channel 2)
  GROW      the body whose mean pairwise distance grows fastest across the
            window                                       (README channel 3)
  OWN       the body in the least-owned pair - narrows three to two, so it is
            scored as "is the ejector in that pair"      (nodes section 5)

EARLY vs LATE matters more than the channels. VEL and GROW as written read the
approach to ejection, so they identify a body that is already leaving. Scored
in two windows:

  LATE   the last 20 nodes before disruption   - identification
  EARLY  the first 150 nodes                   - prediction

Only the EARLY column is a forecast. The mass prior works from t=0, so a
channel that needs the late window has to beat the prior by enough to justify
needing the trajectory at all.
"""
import numpy as np, json, time, sys
sys.path.insert(0, '/home/claude/partners')
from node_partners import accel, PAIRS
from firstnode2 import advance

DT = 0.004
T_MAX = 4000.0
EJECT_R = 20.0
WIN_EARLY = 150
WIN_LATE = 20
N_TRY = 110


def ic(masses, a_in, a_out, e_in, e_out):
    """Hierarchical construction with arbitrary inner masses."""
    m = np.array(masses, float)
    pos = np.zeros((3, 3)); vel = np.zeros((3, 3))
    m01 = m[0] + m[1]
    r_in = a_in * (1 - e_in)
    v_in = np.sqrt(m01 * (1 + e_in) / (a_in * (1 - e_in)))
    pos[0] = [-m[1] / m01 * r_in, 0, 0]
    pos[1] = [m[0] / m01 * r_in, 0, 0]
    vel[0] = [0, -m[1] / m01 * v_in, 0]
    vel[1] = [0, m[0] / m01 * v_in, 0]
    r_out = a_out * (1 + e_out)
    v_out = np.sqrt(m.sum() * (1 - e_out) / (a_out * (1 + e_out)))
    pos[2] = [r_out, 0, 0]
    vel[2] = [0, v_out, 0]
    pos -= np.average(pos, weights=m, axis=0)
    vel -= np.average(vel, weights=m, axis=0)
    return pos, vel


def run(pos, vel, masses):
    pos = pos.copy(); vel = vel.copy()

    def pd(p):
        return np.array([np.linalg.norm(p[a] - p[b]) for a, b in PAIRS])

    d2 = pd(pos); d1 = pd(pos + vel * (DT / 2))
    nodes = []            # (t, owner, H, speeds[3], dists[3])
    t = 0.0
    n = int(T_MAX / DT)
    for s in range(n):
        vel += accel(pos, masses) * DT
        pos += vel * DT
        t += DT
        d0 = pd(pos)
        hit = [p for p in range(3) if d1[p] < d2[p] and d1[p] < d0[p]]
        if hit:
            v2 = float(np.var(d1))
            cv = np.average(vel, weights=masses, axis=0)
            sp = np.linalg.norm(vel - cv, axis=1)
            nodes.append((t - DT, hit[0], v2 / (1 + v2),
                          sp.copy(), d1.copy()))
        d2, d1 = d1, d0
        if s % 50 == 0:
            com = np.average(pos, weights=masses, axis=0)
            r = np.linalg.norm(pos - com, axis=1)
            if r.max() > EJECT_R:
                return nodes, t, int(r.argmax())
    return nodes, None, None


def body_dists(dists):
    """Mean of the two pairwise distances each body participates in."""
    # PAIRS = [(0,1),(1,2),(0,2)]
    return np.array([(dists[0] + dists[2]) / 2,
                     (dists[0] + dists[1]) / 2,
                     (dists[1] + dists[2]) / 2])


def predict(nodes, lo, hi):
    """Channel predictions from nodes[lo:hi]."""
    seg = nodes[lo:hi]
    if len(seg) < 6:
        return None
    H = np.array([s[2] for s in seg])
    k = int(np.argmax(H))
    vel_pick = int(np.argmax(seg[k][3]))

    D = np.array([body_dists(s[4]) for s in seg])
    x = np.arange(len(seg))
    slopes = [np.polyfit(x, D[:, b], 1)[0] for b in range(3)]
    grow_pick = int(np.argmax(slopes))

    own = np.array([sum(1 for s in seg if s[1] == p) for p in range(3)])
    starved = int(np.argmin(own))
    return dict(vel=vel_pick, grow=grow_pick, starved_pair=starved)


if __name__ == "__main__":
    rng = np.random.default_rng(4242)
    rows = []
    t0 = time.time()
    for k in range(N_TRY):
        m = np.array([float(rng.uniform(0.6, 1.4)),
                      float(rng.uniform(0.6, 1.4)),
                      float(rng.uniform(0.5, 1.5))])
        a_out = float(rng.uniform(2.5, 4.5))
        e_in = float(rng.uniform(0.0, 0.25))
        e_out = float(rng.uniform(0.2, 0.6))
        pos, vel = ic(m, 1.0, a_out, e_in, e_out)
        P_in = 2 * np.pi * np.sqrt(1.0 / (m[0] + m[1]))
        pos, vel = advance(pos, vel, m, float(rng.uniform(0.0, 3.0) * P_in))
        nodes, td, body = run(pos, vel, m)
        if td is None or body is None or len(nodes) < WIN_EARLY + WIN_LATE:
            continue
        early = predict(nodes, 0, WIN_EARLY)
        late = predict(nodes, len(nodes) - WIN_LATE, len(nodes))
        if early is None or late is None:
            continue
        rows.append(dict(m=[float(x) for x in m], a_out=a_out,
                         e_in=e_in, e_out=e_out, t_disrupt=td, body=body,
                         mass_pick=int(np.argmin(m)), n_nodes=len(nodes),
                         early=early, late=late))
        if (k + 1) % 20 == 0:
            print(f"  {k+1}/{N_TRY} tried, {len(rows)} usable  "
                  f"[{time.time()-t0:.0f}s]", flush=True)

    json.dump(rows, open('/home/claude/partners/identity_test.json', 'w'),
              indent=2, default=str)

    n = len(rows)
    print(f"\n{'='*72}")
    print(f"{n} disrupting systems with all three masses varied  "
          f"[{time.time()-t0:.0f}s]")
    if n < 10:
        print("too few"); sys.exit(0)

    b = np.array([r['body'] for r in rows])
    mp = np.array([r['mass_pick'] for r in rows])
    print(f"\nejected body: " + "  ".join(
        f"{k}: {int((b==k).sum())} ({100*(b==k).mean():.0f}%)" for k in range(3)))
    print(f"mass prior (lightest body) correct: "
          f"{int((mp==b).sum())}/{n} = {100*(mp==b).mean():.1f}%")
    print(f"  -> the prior now has real work to do (it was 96.5% before)")

    print(f"\n{'channel':<34}{'EARLY (forecast)':>20}{'LATE (identify)':>20}")
    for key, lbl in (('vel', 'velocity at highest-H node'),
                     ('grow', 'fastest-growing distances')):
        e = np.array([r['early'][key] for r in rows])
        l = np.array([r['late'][key] for r in rows])
        print(f"  {lbl:<32}{100*(e==b).mean():>18.1f}%"
              f"{100*(l==b).mean():>19.1f}%")
    # starved pair: is the ejector in it?
    def in_pair(p, body):
        return body in PAIRS[p]
    e = np.array([in_pair(r['early']['starved_pair'], r['body']) for r in rows])
    l = np.array([in_pair(r['late']['starved_pair'], r['body']) for r in rows])
    print(f"  {'ejector in least-owned pair':<32}{100*e.mean():>18.1f}%"
          f"{100*l.mean():>19.1f}%")
    print(f"  {'(chance for a 2-of-3 pair = 66.7%)':<32}")
    print(f"\n  {'mass prior, for comparison':<32}{100*(mp==b).mean():>18.1f}%"
          f"{100*(mp==b).mean():>19.1f}%")

    # combined: majority vote of vel+grow+mass, early
    print(f"\ncombined (majority of velocity, growth, mass), EARLY:")
    ev = np.array([r['early']['vel'] for r in rows])
    eg = np.array([r['early']['grow'] for r in rows])
    comb = []
    for i in range(n):
        votes = np.bincount([ev[i], eg[i], mp[i]], minlength=3)
        comb.append(int(np.argmax(votes)))
    comb = np.array(comb)
    print(f"  {100*(comb==b).mean():.1f}%  (mass prior alone "
          f"{100*(mp==b).mean():.1f}%)")
