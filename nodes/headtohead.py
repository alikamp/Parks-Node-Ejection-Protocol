"""Per-partner jitter vs the existing H-slope signal, as forecasts.

Found so far: per-partner interval jitter (interval between successive minima
of the SAME pair) separates ejectors from survivors where the published pooled
jitter cannot.  But two things make that result weaker than it looked:

  - "survivors" are RIGHT-CENSORED at the 4000t observation window, not
    negatives.  A survivor at ppj 0.045 may eject at t=5000.
  - ppj measured over a whole run includes the run-up to ejection, so it is
    not a forecast.

This fixes both.  Every signal is measured in the SAME early window (the
first K nodes), and the outcome is whatever happens afterwards.  Systems that
eject inside the window are excluded, since there is nothing to forecast.

Signals compared, all from the early window only:
  H_slope   dH/dn - the existing Regime 2 signal
  H_mean    level of H
  ppj       per-partner interval jitter  (the candidate)
  pooled    pooled internodal jitter     (the published metric)

Scored two ways: rank correlation against time-to-ejection among ejectors,
and separation between ejectors and censored survivors.
"""
import numpy as np, json, itertools, time, sys
from scipy import stats
sys.path.insert(0, '/home/claude/partners')
from node_partners import run, hierarchical_ic, PNAME

K = 150          # nodes in the measurement window
MAX_T = 4000.0
DT = 0.004


def signals(T, P, H, k=K):
    """All signals from the first k nodes only."""
    T, P, H = T[:k], P[:k], H[:k]
    n = np.arange(len(H))
    sl = stats.linregress(n, H)

    dtn = np.diff(T)
    pooled = float(dtn.std() / dtn.mean())

    pj = []
    for p in range(3):
        tp = T[P == p]
        if len(tp) > 4:
            d = np.diff(tp)
            pj.append(d.std() / d.mean())
    return dict(H_slope=float(sl.slope), H_slope_p=float(sl.pvalue),
                H_mean=float(H.mean()), pooled=pooled,
                ppj=float(np.mean(pj)) if pj else np.nan,
                ppj_max=float(np.max(pj)) if pj else np.nan,
                t_window_end=float(T[-1]))


if __name__ == "__main__":
    GRID = list(itertools.product([3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.5, 8.0],
                                  [0.01, 0.1], [0.3, 0.5]))
    rows = []
    t0 = time.time()
    for a_out, m_out, e_out in GRID:
        masses = np.array([1.0, 1.0, m_out])
        p, v = hierarchical_ic(masses, 1.0, a_out, 0.15, e_out)
        rec, te, be, tend = run(p, v, masses, dt=DT, steps=int(MAX_T / DT))

        if len(rec) < K + 10:
            print(f"a_out={a_out:4.1f} m={m_out:4.2f} e={e_out:3.1f}  "
                  f"only {len(rec)} nodes - skip", flush=True)
            continue

        T = np.array([r[0] for r in rec])
        P = np.array([r[1] for r in rec])
        H = np.array([r[2] for r in rec])
        s = signals(T, P, H)

        # the window must precede the outcome
        if te is not None and te <= s['t_window_end']:
            print(f"a_out={a_out:4.1f} m={m_out:4.2f} e={e_out:3.1f}  "
                  f"ejects t={te:.1f} inside window (ends t="
                  f"{s['t_window_end']:.1f}) - excluded", flush=True)
            continue

        row = dict(a_out=a_out, m_out=m_out, e_out=e_out,
                   ejected=te is not None, t_eject=te, censored_at=MAX_T,
                   **s)
        rows.append(row)
        print(f"a_out={a_out:4.1f} m={m_out:4.2f} e={e_out:3.1f} | "
              f"{'EJECT t=%7.1f' % te if te else 'censored @4000'} | "
              f"dH/dn={s['H_slope']:+.6f} H={s['H_mean']:.4f} "
              f"ppj={s['ppj']:.4f} pooled={s['pooled']:.4f} "
              f"[{time.time()-t0:.0f}s]", flush=True)

    json.dump(rows, open('/home/claude/partners/headtohead.json', 'w'),
              indent=2, default=str)

    ej = [r for r in rows if r['ejected']]
    sv = [r for r in rows if not r['ejected']]
    print(f"\n{'='*74}")
    print(f"window = first {K} nodes.  ejectors {len(ej)}, censored {len(sv)}")

    KEYS = [('H_slope', 'dH/dn (existing)'), ('H_mean', 'H level'),
            ('ppj', 'per-partner jitter'), ('ppj_max', 'per-partner jitter max'),
            ('pooled', 'pooled jitter (published)')]

    if len(ej) >= 5:
        print(f"\n-- rank correlation with time-to-ejection (ejectors only, "
              f"n={len(ej)}) --")
        print(f"   negative rho = higher signal means sooner ejection")
        y = np.array([r['t_eject'] for r in ej])
        for k, lbl in KEYS:
            x = np.array([r[k] for r in ej])
            m = np.isfinite(x)
            if m.sum() >= 5:
                sp = stats.spearmanr(x[m], y[m])
                print(f"   {lbl:28s} rho={sp.statistic:+.4f}  p={sp.pvalue:.4f}")

    if len(ej) >= 3 and len(sv) >= 3:
        print(f"\n-- separation: ejectors vs censored survivors --")
        for k, lbl in KEYS:
            a = np.array([r[k] for r in ej]); a = a[np.isfinite(a)]
            b = np.array([r[k] for r in sv]); b = b[np.isfinite(b)]
            if len(a) >= 3 and len(b) >= 3:
                u = stats.mannwhitneyu(a, b, alternative='two-sided')
                # AUC from U
                auc = u.statistic / (len(a) * len(b))
                print(f"   {lbl:28s} eject {a.mean():+.5f}  cens {b.mean():+.5f}"
                      f"   AUC={max(auc,1-auc):.3f}  p={u.pvalue:.4f}")
                print(f"   {'':28s}   overlap: eject [{a.min():+.4f},"
                      f"{a.max():+.4f}]  cens [{b.min():+.4f},{b.max():+.4f}]")

    print(f"\n-- confound: does the signal just track the orbital elements? --")
    if len(rows) > 6:
        for k, lbl in KEYS:
            x = np.array([r[k] for r in rows])
            m = np.isfinite(x)
            sa = stats.spearmanr(x[m], [r['a_out'] for r, mm in
                                        zip(rows, m) if mm])
            se = stats.spearmanr(x[m], [r['e_out'] for r, mm in
                                        zip(rows, m) if mm])
            print(f"   {lbl:28s} rho(a_out)={sa.statistic:+.3f}  "
                  f"rho(e_out)={se.statistic:+.3f}")
        print("   (a signal is only useful on an observed system if it can be"
              "\n    measured from the trajectory - which these can and the"
              "\n    orbital elements cannot - but strong tracking means it is"
              "\n    a readout of the elements, not independent information)")

    print(f"\n{'a_out':>6}{'m':>6}{'e':>5}{'outcome':>16}"
          f"{'dH/dn':>11}{'ppj':>9}{'pooled':>9}")
    for r in sorted(rows, key=lambda x: (x['ejected'], x.get('t_eject') or 9e9)):
        print(f"{r['a_out']:6.1f}{r['m_out']:6.2f}{r['e_out']:5.1f}"
              f"{('eject %.0f' % r['t_eject']) if r['ejected'] else 'censored':>16}"
              f"{r['H_slope']:+11.6f}{r['ppj']:9.4f}{r['pooled']:9.4f}")
    print("\nsaved /home/claude/partners/headtohead.json")
