"""What does analytic first-node prediction actually save?

The Newton predictor needs ONE force evaluation (it requires the relative
acceleration to form d(r.v)/dt).  Locating the first node by integration
needs one force evaluation per step until the node is bracketed.

Measured here, in force evaluations, on randomised-phase states:

  A  integrate to the first node          (what the predictor replaces)
  B  Newton predictor                     (= 1, by construction)
  C  integrate to 150 nodes               (the classification window used
                                           by the headtohead forecast test)
  D  integrate to ejection                (brute force, the baseline the
                                           ~99% reduction claim is against)

The ratio A/B is the saving on the first-node task.  A/C is the fraction of
the classification pipeline that task represents - which is what decides
whether the saving matters end to end.
"""
import numpy as np, json, itertools, time, sys
sys.path.insert(0, '/home/claude/partners')
from node_partners import accel, hierarchical_ic, PAIRS
from firstnode2 import advance

DT = 0.004                 # the repo's integration step


def steps_to_nodes(pos, vel, masses, n_want, dt=DT, max_steps=4_000_000,
                   eject_r=20.0):
    """Force evaluations until n_want nodes are seen, or ejection."""
    pos = pos.copy(); vel = vel.copy()

    def pd(p):
        return np.array([np.linalg.norm(p[a] - p[b]) for a, b in PAIRS])

    d2 = pd(pos)
    d1 = pd(pos + vel * (dt / 2))
    n_nodes = 0
    fev = 0
    for s in range(max_steps):
        a = accel(pos, masses); fev += 1
        vel += a * dt
        pos += vel * dt
        d0 = pd(pos)
        for p in range(3):
            if d1[p] < d2[p] and d1[p] < d0[p]:
                n_nodes += 1
                if n_nodes >= n_want:
                    return fev, (s + 1) * dt, False
        d2, d1 = d1, d0
        if s % 50 == 0:
            com = np.average(pos, weights=masses, axis=0)
            if np.linalg.norm(pos - com, axis=1).max() > eject_r:
                return fev, (s + 1) * dt, True
    return fev, max_steps * dt, False


if __name__ == "__main__":
    GRID = list(itertools.product([3.0, 3.5, 4.0, 4.5, 5.0, 5.5],
                                  [0.01, 0.1, 0.5],
                                  [0.2, 0.3, 0.4, 0.5]))
    rng = np.random.default_rng(31)
    rows = []
    t0 = time.time()
    P_in = 2 * np.pi * np.sqrt(0.5)

    for a_out, m_out, e_out in GRID:
        masses = np.array([1.0, 1.0, m_out])
        pos, vel = hierarchical_ic(masses, 1.0, a_out, 0.15, e_out)
        pos, vel = advance(pos, vel, masses,
                           float(rng.uniform(0.1, 3.0) * P_in))

        f1, t1, _ = steps_to_nodes(pos, vel, masses, 1)
        f150, t150, ej150 = steps_to_nodes(pos, vel, masses, 150)
        fej, tej, ejected = steps_to_nodes(pos, vel, masses, 10**9)

        rows.append(dict(a_out=a_out, m_out=m_out, e_out=e_out,
                         fev_first=f1, t_first=t1,
                         fev_150=f150, t_150=t150, ej_in_150=ej150,
                         fev_eject=fej, t_eject=tej, ejected=ejected))
        print(f"a_out={a_out:4.1f} m={m_out:4.2f} e={e_out:3.1f} | "
              f"1st node {f1:6d} fev (t={t1:6.2f}) | "
              f"150 nodes {f150:8d} fev | "
              f"{'eject' if ejected else 'cens '} {fej:9d} fev "
              f"(t={tej:7.1f})  [{time.time()-t0:.0f}s]", flush=True)

    json.dump(rows, open('/home/claude/partners/compute_cost.json', 'w'),
              indent=2, default=str)

    f1 = np.array([r['fev_first'] for r in rows], float)
    f150 = np.array([r['fev_150'] for r in rows], float)
    fej = np.array([r['fev_eject'] for r in rows], float)

    print(f"\n{'='*72}\nFORCE EVALUATIONS ({len(rows)} systems, dt={DT})")
    print(f"{'task':<34}{'median':>12}{'min':>10}{'max':>12}")
    print(f"{'A integrate to first node':<34}{np.median(f1):12.0f}"
          f"{f1.min():10.0f}{f1.max():12.0f}")
    print(f"{'B Newton predictor':<34}{1:12d}{1:10d}{1:12d}")
    print(f"{'C integrate to 150 nodes':<34}{np.median(f150):12.0f}"
          f"{f150.min():10.0f}{f150.max():12.0f}")
    print(f"{'D integrate to ejection':<34}{np.median(fej):12.0f}"
          f"{fej.min():10.0f}{fej.max():12.0f}")

    print(f"\nRATIOS (median)")
    print(f"  A/B  saving on the first-node task      "
          f"{np.median(f1):8.0f}x")
    print(f"  A/C  first node as a share of the       "
          f"{100*np.median(f1/f150):8.2f}%")
    print(f"       classification window")
    print(f"  C/D  classification window as a share   "
          f"{100*np.median(f150/fej):8.2f}%")
    print(f"       of brute force to ejection")
    print(f"  A/D  first node as a share of brute     "
          f"{100*np.median(f1/fej):8.4f}%")
    print(f"       force to ejection")

    print(f"\nINTERPRETATION")
    print(f"  Replacing integration-to-first-node with the Newton step is a")
    print(f"  ~{np.median(f1):.0f}x saving on that task, but that task is only")
    print(f"  {100*np.median(f1/f150):.2f}% of the classification window, so end to end it")
    print(f"  changes little. The classification window is already")
    print(f"  {100*np.median(f150/fej):.1f}% of brute force - that is where the existing")
    print(f"  ~99% reduction comes from.")
    print(f"\n  The analytic node matters only if the node SEQUENCE can be")
    print(f"  propagated without integrating, or if classification can be")
    print(f"  done from the first node alone (zero integration).")
