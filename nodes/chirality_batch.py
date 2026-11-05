"""Does partner-rotation chirality predict ejection, or just track period ratio?

Found on n=3: the three pairwise distances take strict turns owning each node
(switch rate ~1.0, equal thirds) because each hits a minimum once per inner
binary period.  The rotation never breaks.  But its DIRECTION does vary:

    chi = (F - B)/(F + B),  F = n(01->12)+n(12->02)+n(02->01)
                            B = n(01->02)+n(02->12)+n(12->01)

    README cfg  (ejects t=2012)  chi = +0.534
    script cfg  (survives 8000t) chi = -0.017
    stable ctrl (survives 4000t) chi = +0.015

CONFOUND: those configs differ in a_out (4.5 vs 5.5 vs 12) and m_out (0.01 vs
0.1), so chi may simply track the period ratio P_out/P_in rather than carry
decay information.  This batch varies a_out, m_out and e_out independently and
asks whether chi separates ejectors AT FIXED period ratio.

Also computes PER-PARTNER interval jitter.  Pooled jitter is ~0.52 even in a
stable control because the round-robin puts unequal gaps between consecutive
nodes, so the 0.02 Nodal Breathing threshold fires at node 0 everywhere.
Intervals between successive minima of the SAME pair remove that structure.
"""
import numpy as np, json, itertools, time, sys
sys.path.insert(0, '/home/claude/partners')
from node_partners import run, hierarchical_ic, PAIRS, PNAME

FWD = {(0, 1), (1, 2), (2, 0)}      # 01->12, 12->02, 02->01


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


def per_partner_jitter(T, P, win=12):
    """Jitter computed within each partner channel separately."""
    out = {}
    for p in range(3):
        tp = T[P == p]
        if len(tp) < win + 3:
            out[PNAME[p]] = None
            continue
        d = np.diff(tp)
        glob = d.std() / d.mean()
        roll = np.array([d[i:i + win].std() / d[i:i + win].mean()
                         for i in range(len(d) - win)])
        out[PNAME[p]] = dict(global_jitter=float(glob),
                             roll_max=float(roll.max()),
                             roll_med=float(np.median(roll)),
                             mean_interval=float(d.mean()), n=len(d))
    return out


def periods(masses, a_in, a_out):
    M_in = masses[0] + masses[1]
    P_in = 2 * np.pi * np.sqrt(a_in ** 3 / M_in)
    P_out = 2 * np.pi * np.sqrt(a_out ** 3 / masses.sum())
    return P_in, P_out, P_out / P_in


if __name__ == "__main__":
    A_OUT = [3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.5, 8.0]
    M_OUT = [0.01, 0.1]
    E_OUT = [0.3, 0.5]
    MAX_T = 4000.0
    DT = 0.004
    STEPS = int(MAX_T / DT)

    rows = []
    t0 = time.time()
    for a_out, m_out, e_out in itertools.product(A_OUT, M_OUT, E_OUT):
        masses = np.array([1.0, 1.0, m_out])
        p, v = hierarchical_ic(masses, 1.0, a_out, 0.15, e_out)
        rec, te, be, tend = run(p, v, masses, dt=DT, steps=STEPS)
        P_in, P_out, ratio = periods(masses, 1.0, a_out)

        if len(rec) < 60:
            rows.append(dict(a_out=a_out, m_out=m_out, e_out=e_out,
                             ratio=ratio, n_nodes=len(rec), skipped=True))
            print(f"a_out={a_out:4.1f} m={m_out:4.2f} e={e_out:3.1f}  "
                  f"only {len(rec)} nodes - skip", flush=True)
            continue

        T = np.array([r[0] for r in rec])
        P = np.array([r[1] for r in rec])
        E = np.array([r[3] for r in rec])
        prim = np.array([P[E == e][0] for e in np.unique(E)])
        Tp = np.array([T[E == e][0] for e in np.unique(E)])

        chi, F, B = chirality(prim)
        sw = float((prim[1:] != prim[:-1]).mean())
        dtn = np.diff(Tp)
        pooled = float(dtn.std() / dtn.mean())
        ppj = per_partner_jitter(T, P)
        ppj_vals = [v['global_jitter'] for v in ppj.values() if v]

        row = dict(a_out=a_out, m_out=m_out, e_out=e_out, ratio=float(ratio),
                   ejected=te is not None, t_eject=te, body_eject=be,
                   t_end=tend, n_nodes=len(prim), chi=float(chi), F=F, B=B,
                   switch_rate=sw, pooled_jitter=pooled,
                   per_partner_jitter=ppj,
                   ppj_mean=float(np.mean(ppj_vals)) if ppj_vals else None,
                   own=[float((prim == k).mean()) for k in range(3)])
        rows.append(row)
        print(f"a_out={a_out:4.1f} m={m_out:4.2f} e={e_out:3.1f} "
              f"P_out/P_in={ratio:6.2f} | "
              f"{'EJECT t=%7.1f b%d' % (te, be) if te else 'survive       '} | "
              f"chi={chi:+.3f} sw={sw:.3f} pooled_j={pooled:.3f} "
              f"ppj={row['ppj_mean'] if row['ppj_mean'] is None else '%.3f'%row['ppj_mean']} "
              f"n={len(prim)}  [{time.time()-t0:.0f}s]", flush=True)

    json.dump(rows, open('/home/claude/partners/chirality.json', 'w'),
              indent=2, default=str)

    # ---- summary ------------------------------------------------------------
    good = [r for r in rows if not r.get('skipped')]
    ej = [r for r in good if r['ejected']]
    sv = [r for r in good if not r['ejected']]

    print(f"\n{'='*72}")
    print(f"ejectors: {len(ej)}   survivors: {len(sv)}")
    if ej and sv:
        for key, lbl in (('chi', 'chirality'),
                         ('pooled_jitter', 'pooled jitter'),
                         ('ppj_mean', 'per-partner jitter'),
                         ('switch_rate', 'switch rate')):
            a = np.array([r[key] for r in ej if r[key] is not None])
            b = np.array([r[key] for r in sv if r[key] is not None])
            if len(a) and len(b):
                from scipy import stats
                u = stats.mannwhitneyu(a, b, alternative='two-sided')
                print(f"  {lbl:22s} eject {a.mean():+.4f}  survive {b.mean():+.4f}"
                      f"   MWU p={u.pvalue:.4f}")

    print(f"\n-- chi vs period ratio (the confound) --")
    if len(good) > 3:
        from scipy import stats
        rr = np.array([r['ratio'] for r in good])
        cc = np.array([r['chi'] for r in good])
        sp = stats.spearmanr(rr, cc)
        print(f"  Spearman(chi, P_out/P_in) = {sp.statistic:+.4f}  p={sp.pvalue:.4f}")
        print(f"  -> if strong, chi is kinematic, not a decay signal")
        print(f"\n  {'a_out':>6}{'ratio':>8}{'chi':>9}{'ejected':>9}")
        for r in sorted(good, key=lambda x: x['ratio']):
            print(f"  {r['a_out']:6.1f}{r['ratio']:8.2f}{r['chi']:+9.3f}"
                  f"{str(r['ejected']):>9}")
    print("\nsaved /home/claude/partners/chirality.json")
