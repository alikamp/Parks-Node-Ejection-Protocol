"""Do the defect and per-partner-jitter results survive generic phases?

Every result in sections 1-5 was measured on hierarchical_ic states, which
place all three bodies collinear on the x-axis with velocities in the y-z
plane - so r.v = 0 exactly for all three pairs at t=0 (section 6).  That is
the same degeneracy found in the published validation, where accuracy went
98.9% -> 93.0% and "zero missed ejections" did not survive.

Two reasons to doubt the partner results transfer:
  - firstnode2.py found the rotation markedly less rigid under random phase:
    dominant order (0,1,2) in only ~55% of systems, and 57/192 with no clean
    3-cycle at all (though that was from only 40 nodes per system)
  - the aligned start puts every pair simultaneously at a node, which is
    exactly the phase relationship the round-robin depends on

Same grid as chirality_batch.py, two random phases per configuration, so the
comparison against the aligned numbers is direct.

Aligned results to confirm or refute:
  defect rate        3.401 vs 0.044 per 1000 nodes (eject vs censored)
  defect incidence   11/12 ejectors, 2/18 censored, Fisher p=1.5e-5
  defect lead        median lead fraction 0.788
  switch rate        survivors EXACTLY 1.0000, ejectors 0.9915
  per-partner jitter AUC 0.972, rho=-0.721 vs t_eject (first 150 nodes)
  pooled jitter      AUC 0.956, discriminating near 0.60 (not 0.02)
"""
import numpy as np, json, itertools, time, sys
from scipy import stats
sys.path.insert(0, '/home/claude/partners')
from node_partners import accel, hierarchical_ic, PAIRS, PNAME
from firstnode2 import advance

DT = 0.004
T_MAX = 4000.0
SAMP = 5
DTS = SAMP * DT
EDGES = np.array([0.0, 0.07, 0.22, 0.60, 1.01])
WIN = 150            # early-window length, matching headtohead.py


def run(pos, vel, masses, dt=DT, t_max=T_MAX, eject_r=20.0):
    """Node sequence: (t, pair, H) per event, one entry per event."""
    pos = pos.copy(); vel = vel.copy()

    def pd(p):
        return np.array([np.linalg.norm(p[a] - p[b]) for a, b in PAIRS])

    d2 = pd(pos)
    d1 = pd(pos + vel * (dt / 2))
    T, P, H = [], [], []
    t = 0.0
    te = None
    be = None
    n_coincident = 0
    n = int(t_max / dt)
    for s in range(n):
        vel += accel(pos, masses) * dt
        pos += vel * dt
        t += dt
        d0 = pd(pos)
        hit = [p for p in range(3) if d1[p] < d2[p] and d1[p] < d0[p]]
        if hit:
            v = float(np.var(d1))
            # Tie-break MUST match the aligned runs (chirality_batch.py took
            # P[E==e][0], the first qualifying pair by index).  Coincident
            # firings are only 0.16% of records, but defects are rare too
            # (37 in 10,878 nodes), so a different rule here would confound
            # phase against bookkeeping in the very comparison being made.
            if len(hit) > 1:
                n_coincident += 1
            owner = hit[0]
            T.append(t - dt); P.append(owner); H.append(v / (1 + v))
        d2, d1 = d1, d0
        if s % 50 == 0:
            com = np.average(pos, weights=masses, axis=0)
            r = np.linalg.norm(pos - com, axis=1)
            if r.max() > eject_r:
                te = t; be = int(r.argmax()); break
    return (np.array(T), np.array(P, int), np.array(H), te, be, t,
            n_coincident)


def ppj(T, P, k=None):
    """Mean per-partner interval jitter; k limits to the first k nodes."""
    if k:
        T, P = T[:k], P[:k]
    out = []
    for p in range(3):
        tp = T[P == p]
        if len(tp) > 4:
            d = np.diff(tp)
            if d.mean() > 0:
                out.append(d.std() / d.mean())
    return float(np.mean(out)) if out else np.nan


def pooled(T, k=None):
    t = T[:k] if k else T
    d = np.diff(t)
    return float(d.std() / d.mean()) if len(d) > 2 and d.mean() > 0 else np.nan


def defects(T, P):
    return [(float(T[i]), int(P[i]))
            for i in range(1, len(P)) if P[i] == P[i - 1]]


if __name__ == "__main__":
    GRID = list(itertools.product([3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.5, 8.0],
                                  [0.01, 0.1], [0.3, 0.5]))
    rng = np.random.default_rng(777)
    P_in = 2 * np.pi * np.sqrt(0.5)
    rows = []
    t0 = time.time()

    for a_out, m_out, e_out in GRID:
        for rep in range(2):
            masses = np.array([1.0, 1.0, m_out])
            pos, vel = hierarchical_ic(masses, 1.0, a_out, 0.15, e_out)
            pos, vel = advance(pos, vel, masses,
                               float(rng.uniform(0.0, 3.0) * P_in))
            T, P, H, te, be, tend, ncoin = run(pos, vel, masses)
            if len(T) < WIN + 10:
                print(f"a_out={a_out:4.1f} m={m_out:4.2f} e={e_out:3.1f} "
                      f"r{rep}  only {len(T)} nodes - skip", flush=True)
                continue

            d = defects(T, P)
            life = (te if te else tend) - T[0]
            sl = stats.linregress(np.arange(WIN), H[:WIN])
            row = dict(a_out=a_out, m_out=m_out, e_out=e_out, rep=rep,
                       ejected=te is not None, t_eject=te, t_end=tend,
                       n_nodes=len(T), n_defects=len(d),
                       switch_rate=float((P[1:] != P[:-1]).mean()),
                       ppj_early=ppj(T, P, WIN), pooled_early=pooled(T, WIN),
                       ppj_full=ppj(T, P), pooled_full=pooled(T),
                       H_slope_early=float(sl.slope),
                       own=[float((P == k).mean()) for k in range(3)],
                       window_end=float(T[WIN - 1]),
                       n_coincident=ncoin)
            if d:
                row['t_first_defect'] = d[0][0]
                row['lead'] = (te if te else tend) - d[0][0]
                row['lead_frac'] = row['lead'] / life if life > 0 else None
            rows.append(row)
            print(f"a_out={a_out:4.1f} m={m_out:4.2f} e={e_out:3.1f} r{rep} | "
                  f"{'EJECT %7.1f' % te if te else 'censored      '} | "
                  f"nodes {len(T):5d} def {len(d):3d} sw {row['switch_rate']:.4f} "
                  f"ppj {row['ppj_early']:.4f} pool {row['pooled_early']:.4f}"
                  f"  [{time.time()-t0:.0f}s]", flush=True)

    json.dump(rows, open('/home/claude/partners/generic_phase.json', 'w'),
              indent=2, default=str)

    ej = [r for r in rows if r['ejected']]
    sv = [r for r in rows if not r['ejected']]
    print(f"\n{'='*76}")
    print(f"GENERIC PHASE: {len(rows)} systems, {len(ej)} eject, "
          f"{len(sv)} censored   [{time.time()-t0:.0f}s]")

    # ---- 1. defects -----------------------------------------------------
    de, ne = sum(r['n_defects'] for r in ej), sum(r['n_nodes'] for r in ej)
    ds, ns = sum(r['n_defects'] for r in sv), sum(r['n_nodes'] for r in sv)
    ie = sum(1 for r in ej if r['n_defects'] > 0)
    iss = sum(1 for r in sv if r['n_defects'] > 0)
    print(f"\n--- 1. round-robin defects ---")
    print(f"{'':>10}{'defects':>10}{'nodes':>9}{'rate/1000':>12}{'with >=1':>12}")
    print(f"{'ejectors':>10}{de:10d}{ne:9d}{de/max(ne,1)*1000:12.3f}"
          f"{f'{ie}/{len(ej)}':>12}")
    print(f"{'censored':>10}{ds:10d}{ns:9d}{ds/max(ns,1)*1000:12.3f}"
          f"{f'{iss}/{len(sv)}':>12}")
    if ej and sv:
        odds, p = stats.fisher_exact([[ie, len(ej) - ie],
                                      [iss, len(sv) - iss]])
        print(f"  Fisher exact p = {p:.3e}   odds ratio = {odds:.1f}")
        if ds > 0:
            print(f"  rate ratio = {(de/max(ne,1))/(ds/max(ns,1)):.0f}x")
        print(f"  ALIGNED was: 3.401 vs 0.044 per 1000, 11/12 vs 2/18, "
              f"p=1.5e-5, 78x")

    # ---- 2. defect lead -------------------------------------------------
    lf = [r['lead_frac'] for r in ej
          if r.get('lead_frac') is not None]
    ld = [r['lead'] for r in ej if r.get('lead') is not None]
    print(f"\n--- 2. defect lead time ---")
    if lf:
        print(f"  n={len(lf)}  lead fraction median {np.median(lf):.3f}  "
              f"range {min(lf):.3f}-{max(lf):.3f}")
        print(f"  absolute lead median {np.median(ld):.1f}t")
        te_ = [r['t_eject'] for r in ej if r.get('lead_frac') is not None]
        sp = stats.spearmanr(lf, te_)
        print(f"  Spearman(lead fraction, t_eject) = {sp.statistic:+.3f} "
              f"p={sp.pvalue:.3f}  (scale-free if ~0)")
        print(f"  ALIGNED was: median 0.788, rho=-0.20 p=0.56")

    # ---- 3. switch rate -------------------------------------------------
    print(f"\n--- 3. switch rate ---")
    for lbl, grp in (('ejectors', ej), ('censored', sv)):
        s = np.array([r['switch_rate'] for r in grp])
        if len(s):
            print(f"  {lbl:>10}: mean {s.mean():.4f}  min {s.min():.4f}  "
                  f"exactly 1.0000 in {int((s==1.0).sum())}/{len(s)}")
    if ej and sv:
        u = stats.mannwhitneyu([r['switch_rate'] for r in ej],
                               [r['switch_rate'] for r in sv])
        print(f"  Mann-Whitney p = {u.pvalue:.3e}")
    print(f"  ALIGNED was: survivors exactly 1.0000 (18/18), "
          f"ejectors 0.9915, p=0.0000")

    # ---- 4/5. jitter ----------------------------------------------------
    print(f"\n--- 4. early-window signals (first {WIN} nodes, outcome after) ---")
    fc = [r for r in rows if not (r['ejected'] and
                                  r['t_eject'] <= r['window_end'])]
    e2 = [r for r in fc if r['ejected']]
    s2 = [r for r in fc if not r['ejected']]
    print(f"  usable: {len(e2)} ejectors, {len(s2)} censored "
          f"({len(rows)-len(fc)} eject inside the window)")
    if len(e2) >= 4 and len(s2) >= 4:
        print(f"\n{'signal':<22}{'AUC':>8}{'p':>11}{'rho vs t_eject':>17}{'p':>9}")
        for key, lbl in (('ppj_early', 'per-partner jitter'),
                         ('pooled_early', 'pooled jitter'),
                         ('H_slope_early', 'dH/dn')):
            a = np.array([r[key] for r in e2], float)
            b = np.array([r[key] for r in s2], float)
            ma, mb = np.isfinite(a), np.isfinite(b)
            if ma.sum() < 4 or mb.sum() < 4:
                continue
            u = stats.mannwhitneyu(a[ma], b[mb], alternative='two-sided')
            auc = u.statistic / (ma.sum() * mb.sum())
            y = np.array([r['t_eject'] for r in e2], float)[ma]
            sp = stats.spearmanr(a[ma], y)
            print(f"{lbl:<22}{max(auc,1-auc):8.3f}{u.pvalue:11.4f}"
                  f"{sp.statistic:+17.3f}{sp.pvalue:9.3f}")
        print(f"\n  ALIGNED was: ppj AUC 0.972 rho -0.721 p=0.019; "
              f"pooled 0.956; dH/dn 0.778")

    print(f"\n--- 5. ownership uniformity (aligned: 33+/-3) ---")
    sp_ = np.array([max(r['own']) - min(r['own']) for r in rows])
    print(f"  spread median {100*np.median(sp_):.2f} pp   "
          f"range {100*sp_.min():.2f}-{100*sp_.max():.2f} pp")
    print("\nsaved /home/claude/partners/generic_phase.json")
