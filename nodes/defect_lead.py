"""When does the first round-robin defect occur relative to ejection?

Established: a defect (the same pair owning two consecutive nodes) occurs at
3.401 per 1000 nodes in ejectors against 0.044 in censored survivors - a 78x
rate ratio, Fisher p=1.5e-5, odds ratio 88.  11 of 12 ejectors show at least
one defect; 16 of 18 censored show none.

That makes it a classifier.  Whether it is a FORECAST depends entirely on
timing, which the earlier runs did not record:

  - first defect long before separation  -> precursor
  - defects clustered in the final approach -> damage indicator

This records every defect timestamp for the ejecting systems and measures:
  lead time       t_eject - t_first_defect
  lead fraction   lead / total lifetime  (scale-free across systems)
  defect density  in each quarter of the run, to see if they cluster late

Also reports the same for the two censored systems that did defect, where the
"lead time" is a lower bound - they may eject past the 4000t window.
"""
import numpy as np, json, time, sys
from scipy import stats
sys.path.insert(0, '/home/claude/partners')
from node_partners import run, hierarchical_ic, PNAME

MAX_T = 4000.0
DT = 0.004

# (a_out, m_out, e_out) for the 12 ejectors and the 2 defecting survivors
EJECTORS = [(3.0, 0.01, 0.5), (3.0, 0.10, 0.5), (3.0, 0.10, 0.3),
            (3.5, 0.10, 0.3), (4.5, 0.01, 0.5), (3.0, 0.01, 0.3),
            (4.5, 0.10, 0.5), (4.0, 0.10, 0.5), (3.5, 0.01, 0.3),
            (5.0, 0.01, 0.5), (4.0, 0.01, 0.5), (5.5, 0.10, 0.5)]
CENSORED_DEFECTORS = [(5.0, 0.10, 0.5), (5.5, 0.01, 0.5)]


def defects(T, P):
    """Times at which the same pair owns two consecutive nodes."""
    out = []
    for k in range(1, len(P)):
        if P[k] == P[k - 1]:
            out.append((float(T[k]), int(P[k])))
    return out


def analyse(a_out, m_out, e_out):
    masses = np.array([1.0, 1.0, m_out])
    p, v = hierarchical_ic(masses, 1.0, a_out, 0.15, e_out)
    rec, te, be, tend = run(p, v, masses, dt=DT, steps=int(MAX_T / DT))
    if len(rec) < 20:
        return None
    # one entry per event (collapse coincident multi-pair firings)
    E = np.array([r[3] for r in rec])
    T = np.array([rec[np.where(E == e)[0][0]][0] for e in np.unique(E)])
    P = np.array([rec[np.where(E == e)[0][0]][1] for e in np.unique(E)])
    d = defects(T, P)
    span = (T[0], te if te else tend)
    life = span[1] - span[0]
    out = dict(a_out=a_out, m_out=m_out, e_out=e_out,
               ejected=te is not None, t_eject=te, t_end=tend,
               n_nodes=len(T), n_defects=len(d),
               t_first_node=float(T[0]),
               defect_times=[x[0] for x in d],
               defect_pairs=[x[1] for x in d])
    if d:
        out['t_first_defect'] = d[0][0]
        out['lead'] = (span[1] - d[0][0])
        out['lead_frac'] = out['lead'] / life if life > 0 else None
        # quartile density
        q = np.linspace(span[0], span[1], 5)
        dt_arr = np.array([x[0] for x in d])
        out['quartile_counts'] = [int(((dt_arr >= q[i]) &
                                       (dt_arr < q[i + 1])).sum())
                                  for i in range(4)]
    return out


if __name__ == "__main__":
    t0 = time.time()
    rows = []
    print("--- ejectors ---", flush=True)
    for cfg in EJECTORS:
        r = analyse(*cfg)
        if r is None:
            print(f"  a_out={cfg[0]} m={cfg[1]} e={cfg[2]}: too few nodes")
            continue
        rows.append(r)
        if r['n_defects']:
            print(f"  a_out={r['a_out']:4.1f} m={r['m_out']:4.2f} "
                  f"e={r['e_out']:3.1f} | eject t={r['t_eject']:7.1f} | "
                  f"{r['n_defects']:3d} defects | first t="
                  f"{r['t_first_defect']:7.1f} | lead {r['lead']:7.1f}t "
                  f"({100*r['lead_frac']:5.1f}% of life) | Q"
                  f"{r['quartile_counts']}  [{time.time()-t0:.0f}s]",
                  flush=True)
        else:
            print(f"  a_out={r['a_out']:4.1f} m={r['m_out']:4.2f} "
                  f"e={r['e_out']:3.1f} | eject t={r['t_eject']:7.1f} | "
                  f"NO DEFECTS  [{time.time()-t0:.0f}s]", flush=True)

    print("\n--- censored systems that defected (lead is a lower bound) ---",
          flush=True)
    for cfg in CENSORED_DEFECTORS:
        r = analyse(*cfg)
        if r:
            rows.append(r)
            if r['n_defects']:
                print(f"  a_out={r['a_out']:4.1f} m={r['m_out']:4.2f} "
                      f"e={r['e_out']:3.1f} | censored @{r['t_end']:.0f} | "
                      f"{r['n_defects']} defects | first t="
                      f"{r['t_first_defect']:7.1f} | >={r['lead']:7.1f}t",
                      flush=True)
            else:
                print(f"  a_out={r['a_out']:4.1f} m={r['m_out']:4.2f} "
                      f"e={r['e_out']:3.1f} | censored | NO DEFECTS")

    json.dump(rows, open('/home/claude/partners/defect_lead.json', 'w'),
              indent=2, default=str)

    ej = [r for r in rows if r['ejected'] and r['n_defects']]
    print(f"\n{'='*74}")
    if ej:
        lead = np.array([r['lead'] for r in ej])
        lf = np.array([r['lead_frac'] for r in ej])
        print(f"ejectors with >=1 defect: {len(ej)}")
        print(f"  lead time:     median {np.median(lead):8.1f}t   "
              f"range {lead.min():.1f} - {lead.max():.1f}")
        print(f"  lead fraction: median {np.median(lf):8.3f}    "
              f"range {lf.min():.3f} - {lf.max():.3f}")
        print(f"  (lead fraction = fraction of the system's lifetime "
              f"remaining\n   when the first defect fires; 1.0 = at the very "
              f"start, 0.0 = at ejection)")

        Q = np.array([r['quartile_counts'] for r in ej]).sum(axis=0)
        print(f"\n  defect counts by quartile of lifetime (all ejectors "
              f"pooled):")
        print(f"    Q1 {Q[0]}   Q2 {Q[1]}   Q3 {Q[2]}   Q4 {Q[3]}"
              f"   (Q4 = final quarter before ejection)")
        if Q.sum() > 0:
            chi = stats.chisquare(Q)
            print(f"    uniform-distribution test: chi2={chi.statistic:.2f}"
                  f"  p={chi.pvalue:.4f}")
            print(f"    -> clustered in Q4 means damage indicator; "
                  f"spread means precursor")

        print(f"\n  compare: the existing Regime 2 slope forecast quotes a "
              f"median\n  lead time of 412.6t (README).")
    print("\nsaved /home/claude/partners/defect_lead.json")
