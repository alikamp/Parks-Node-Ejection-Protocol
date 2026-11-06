"""PNEP's per-system disruption forecast vs the Mushkin & Katz random-walk model.

The literature already predicts disruption timescales, so "we predict when it
breaks up" needs a named competitor.  The cleanest is MK20 as quoted by
Hayashi et al. (arXiv:2207.12672, their Eq. 2) - PARAMETER-FREE physics:

  T_d/P_in = 2 (m123 m12 / (m1 m2))^2 sqrt(1-e_out) (rp_out/a_in)^-2
             * exp[ (4 sqrt2 / 3) sqrt(m12/m123) (rp_out/a_in)^(3/2) ]

  rp_out = a_out (1 - e_out)

Their fitted variant (Eq. 19, log10(T_d/P_in) = S x + C) is not usable here:
the text fixes C = 1 but never states S numerically, giving it only in a
figure.  MK20 as above corresponds to S = 1.

PNEP's Regime 2 forecast (DISCOVERY.md) is

  dt_ejection ~ (H_target - H_current) / (dH/dn) * mean_inter_node_interval

TWO ASYMMETRIES TO DECLARE, both favouring MK20:

 1. H_target is given no numerical value anywhere in the repo.  It is
    calibrated here from the data by LEAVE-ONE-OUT median of the H at the last
    node before disruption, so each system's forecast never sees its own
    outcome.  MK20 fits nothing.

 2. PNEP forecasts only where dH/dn > 0 (Regime 2).  A flat or negative slope
    gives a negative or infinite time.  Coverage is therefore reported beside
    accuracy; MK20 always returns a number.

Scored as |log10(T_pred / T_true)| - orders of magnitude - since MK20's own
scatter is quoted at one to two orders.
"""
import numpy as np, json, time, sys
from scipy import stats
sys.path.insert(0, '/home/claude/partners')
from node_partners import accel, hierarchical_ic, PAIRS
from firstnode2 import advance

DT = 0.004
T_MAX = 6000.0
EJECT_R = 20.0
WIN = 150
N_TRY = 70


def mk20_td_over_pin(m1, m2, m3, a_in, a_out, e_out):
    m12 = m1 + m2
    m123 = m12 + m3
    rp = a_out * (1.0 - e_out)
    ratio = rp / a_in
    pref = 2.0 * (m123 * m12 / (m1 * m2)) ** 2 * np.sqrt(1.0 - e_out)
    return pref * ratio ** -2 * np.exp((4.0 * np.sqrt(2.0) / 3.0)
                                       * np.sqrt(m12 / m123) * ratio ** 1.5)


def run(pos, vel, masses, dt=DT, t_max=T_MAX, eject_r=EJECT_R):
    pos = pos.copy(); vel = vel.copy()

    def pd(p):
        return np.array([np.linalg.norm(p[a] - p[b]) for a, b in PAIRS])

    d2 = pd(pos); d1 = pd(pos + vel * (dt / 2))
    T, H = [], []
    t = 0.0
    n = int(t_max / dt)
    for s in range(n):
        vel += accel(pos, masses) * dt
        pos += vel * dt
        t += dt
        d0 = pd(pos)
        if any(d1[p] < d2[p] and d1[p] < d0[p] for p in range(3)):
            v = float(np.var(d1))
            T.append(t - dt); H.append(v / (1.0 + v))
        d2, d1 = d1, d0
        if s % 50 == 0:
            com = np.average(pos, weights=masses, axis=0)
            r = np.linalg.norm(pos - com, axis=1)
            if r.max() > eject_r:
                return np.array(T), np.array(H), t, int(r.argmax())
    return np.array(T), np.array(H), None, None


if __name__ == "__main__":
    rng = np.random.default_rng(606)
    P_in = 2 * np.pi * np.sqrt(1.0 / 2.0)
    rows = []
    t0 = time.time()

    for k in range(N_TRY):
        m3 = float(10 ** rng.uniform(np.log10(0.05), 0.0))
        a_out = float(rng.uniform(2.8, 5.0))
        e_in = float(rng.uniform(0.0, 0.25))
        e_out = float(rng.uniform(0.2, 0.6))
        masses = np.array([1.0, 1.0, m3])
        pos, vel = hierarchical_ic(masses, 1.0, a_out, e_in, e_out)
        pos, vel = advance(pos, vel, masses,
                           float(rng.uniform(0.0, 3.0) * P_in))
        T, H, td, body = run(pos, vel, masses)

        mk = mk20_td_over_pin(1.0, 1.0, m3, 1.0, a_out, e_out) * P_in
        rec = dict(m3=m3, a_out=a_out, e_in=e_in, e_out=e_out,
                   n_nodes=len(T), t_disrupt=td, body=body,
                   mk20=float(mk))
        if td is not None and len(T) > WIN + 5 and td > T[WIN - 1]:
            n = np.arange(WIN)
            sl = stats.linregress(n, H[:WIN])
            dtn = float(np.mean(np.diff(T[:WIN])))
            rec.update(slope=float(sl.slope), H_win=float(H[:WIN].mean()),
                       dt_node=dtn, window_end=float(T[WIN - 1]),
                       H_last=float(H[-1]))
        rows.append(rec)
        if (k + 1) % 10 == 0:
            nd = sum(1 for r in rows if r['t_disrupt'] is not None)
            print(f"  {k+1}/{N_TRY}  disrupted {nd}  "
                  f"[{time.time()-t0:.0f}s]", flush=True)

    json.dump(rows, open('/home/claude/partners/timing_benchmark.json', 'w'),
              indent=2, default=str)

    dis = [r for r in rows if r['t_disrupt'] is not None]
    fc = [r for r in dis if 'slope' in r]
    print(f"\n{'='*76}")
    print(f"{len(rows)} systems, {len(dis)} disrupted within {T_MAX}t, "
          f"{len(fc)} with a usable forecast window  [{time.time()-t0:.0f}s]")

    if len(fc) < 6:
        print("too few forecastable systems"); sys.exit(0)

    # leave-one-out H_target from H at the last node before disruption
    Hl = np.array([r['H_last'] for r in fc])
    preds, trues, used = [], [], []
    for i, r in enumerate(fc):
        Ht = float(np.median(np.delete(Hl, i)))
        if r['slope'] <= 0:
            preds.append(np.nan); trues.append(r['t_disrupt']); used.append(False)
            continue
        dn = (Ht - r['H_win']) / r['slope']
        if dn <= 0:
            preds.append(np.nan); trues.append(r['t_disrupt']); used.append(False)
            continue
        preds.append(r['window_end'] + dn * r['dt_node'])
        trues.append(r['t_disrupt']); used.append(True)
    preds = np.array(preds, float)
    trues = np.array([float(x) for x in trues])
    used = np.array(used)
    mk = np.array([r['mk20'] for r in fc], float)

    print(f"\n  H_target calibrated by leave-one-out median of H at the last")
    print(f"  node before disruption: {np.median(Hl):.4f} "
          f"(range {Hl.min():.4f}-{Hl.max():.4f})")
    print(f"\n  COVERAGE: PNEP produces a forecast for {used.sum()}/{len(fc)} "
          f"= {100*used.mean():.0f}%")
    print(f"            MK20 produces one for {len(fc)}/{len(fc)} = 100%")

    def score(p, t, lbl, mask=None):
        m = np.isfinite(p) & (p > 0) & (t > 0)
        if mask is not None:
            m &= mask
        if m.sum() < 4:
            print(f"  {lbl:<34}too few ({m.sum()})"); return
        e = np.abs(np.log10(p[m] / t[m]))
        sp = stats.spearmanr(p[m], t[m])
        print(f"  {lbl:<34}{m.sum():4d}{np.median(e):10.3f}{e.mean():9.3f}"
              f"{(e<np.log10(2)).mean():9.2f}{(e<np.log10(10)).mean():9.2f}"
              f"{sp.statistic:+9.3f}")

    print(f"\n{'estimator':<34}{'n':>4}{'med dex':>10}{'mean':>9}"
          f"{'<2x':>9}{'<10x':>9}{'rho':>9}")
    score(mk, trues, 'MK20 (parameter-free), all')
    score(preds, trues, 'PNEP slope forecast')
    print(f"\n  same systems only (where PNEP forecasts):")
    score(mk, trues, 'MK20, PNEP-covered subset', used)
    score(preds, trues, 'PNEP, same subset', used)

    print(f"\n  'med dex' = median |log10(pred/true)|, i.e. orders of magnitude.")
    print(f"  Hayashi et al. quote MK20's own scatter at 1-2 orders from")
    print(f"  initial-phase sensitivity alone.")

    print(f"\n  true disruption times: median {np.median(trues):.1f}t, "
          f"range {trues.min():.1f}-{trues.max():.1f}t")
