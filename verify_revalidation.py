"""
PNEP re-validation + liberation orbit-count analysis
=====================================================
Standalone. NumPy + SciPy only. Runs in Colab or locally.

PART A  Re-validates PNEP v12 with ground truth from a LONG integration
        instead of the short, class-dependent cap in the original run_trial.
        Reports a correctly-labelled confusion matrix.

PART B  Tests whether ejection timing is tighter when measured in OUTER
        ORBITS after the liberation (eccentricity-minimum) event than in
        raw time units.

Usage:   python pnep_revalidate.py [n_trials] [seed] [gt_max_t]
Colab:   !python pnep_revalidate.py 100 42 6000
"""
import sys, json, time
import numpy as np

G, SOFT, DT = 1.0, 0.01, 0.006
EJECT_R, MIN_NODES, H_THRESH = 12.0, 10, 0.5

# ── integrator (identical to pnep_v12) ───────────────────────────────────────
def _acc(pos, masses):
    diff = pos[np.newaxis,:,:] - pos[:,np.newaxis,:]
    r2 = np.sum(diff**2, axis=2) + SOFT**2
    np.fill_diagonal(r2, 1.0)
    fac = masses[np.newaxis,:] / (r2 * np.sqrt(r2))
    np.fill_diagonal(fac, 0.0)
    return G * np.einsum('ij,ijk->ik', fac, diff)

def leapfrog(pos, vel, masses, dt=DT):
    a = _acc(pos, masses)
    vh = vel + a*(dt/2)
    p2 = pos + vh*dt
    return p2, vh + _acc(p2, masses)*(dt/2)

def total_energy(pos, vel, masses):
    KE = 0.5*np.sum(masses[:,None]*vel**2); PE = 0.0
    for i in range(len(masses)):
        for j in range(i+1, len(masses)):
            dr = pos[j]-pos[i]
            PE -= G*masses[i]*masses[j]/np.sqrt(float(dr@dr)+SOFT**2)
    return float(KE+PE)

def com_correct(pos, vel, masses):
    tm = masses.sum()
    return (pos - (masses[:,None]*pos).sum(0)/tm,
            vel - (masses[:,None]*vel).sum(0)/tm)

# ── orbital elements ─────────────────────────────────────────────────────────
def _orb(r1, v1, r2, v2, m1, m2):
    dr, dv = r2-r1, v2-v1; mu = G*(m1+m2)
    r = float(np.linalg.norm(dr))+1e-10; v2s = float(dv@dv)
    a = 1./(2./r - v2s/mu)
    ev = (v2s/mu - 1./r)*dr - float(dr@dv)/mu*dv
    return abs(float(a)), float(np.clip(np.linalg.norm(ev),0,0.99)), np.cross(dr,dv)

def ma01_ratio(pos, vel, masses):
    m0,m1,m2 = masses
    ai,ei,Li = _orb(pos[0],vel[0],pos[1],vel[1],m0,m1); mI = m0+m1
    cp = (m0*pos[0]+m1*pos[1])/mI; cv = (m0*vel[0]+m1*vel[1])/mI
    ao,eo,Lo = _orb(cp,cv,pos[2],vel[2],mI,m2); qO = m2/mI
    cosI = float(Li@Lo)/(float(np.linalg.norm(Li))*float(np.linalg.norm(Lo))+1e-12)
    incl = float(np.arccos(np.clip(cosI,-1,1)))
    rhs = 2.8*((1+qO)*(1+eo)/np.sqrt(max(0.01,1-eo)))**0.4*(1-0.3*incl/np.pi)
    return ao*(1-eo)/(ai+1e-9)/(rhs+1e-9)

def outer_period(pos, vel, masses):
    m0,m1,m2 = masses; mI = m0+m1
    cp = (m0*pos[0]+m1*pos[1])/mI; cv = (m0*vel[0]+m1*vel[1])/mI
    ao,eo,_ = _orb(cp,cv,pos[2],vel[2],mI,m2)
    return 2*np.pi*np.sqrt(ao**3/(G*(mI+m2)))

def outer_elements(pos, vel, masses, k):
    """a, e, r of body k relative to the CoM of the other two."""
    idx = [i for i in range(3) if i != k]
    i, j = idx; mI = masses[i]+masses[j]
    cp = (masses[i]*pos[i]+masses[j]*pos[j])/mI
    cv = (masses[i]*vel[i]+masses[j]*vel[j])/mI
    a, e, _ = _orb(cp, cv, pos[k], vel[k], mI, masses[k])
    return a, e, float(np.linalg.norm(pos[k]-cp))

def compute_H(pos):
    d = np.array([np.linalg.norm(pos[1]-pos[0]),
                  np.linalg.norm(pos[2]-pos[1]),
                  np.linalg.norm(pos[2]-pos[0])])
    s2 = float(np.var(d)); return s2/(1.0+s2)

def is_local_min(d, d1, d2):
    if d1 is None or d2 is None: return False
    return any(d2[i] > d1[i] and d[i] > d1[i] for i in range(3))

# ── IC generator (identical to pnep_v12) ─────────────────────────────────────
def make_ic(want_stable, rng, stable_min=5.0, unstable_max=0.45, max_att=2000):
    for _ in range(max_att):
        mA=0.7+rng.random()*0.9; mB=0.6+rng.random()*0.9; mC=0.15+rng.random()*0.8
        aIn=0.3+rng.random()*0.25; eIn=rng.random()*0.10
        Pr=float(np.exp(rng.uniform(np.log(3),np.log(120))))
        aOut=aIn*Pr**(2./3.); eOut=rng.random()*0.40
        incl=abs(float(rng.standard_normal())*25.*np.pi/180.)
        sign=1 if rng.random()>0.5 else -1
        vCI=np.sqrt(G*(mA+mB)/aIn)
        vCO=np.sqrt(G*(mA+mB+mC)/aOut)*np.sqrt((1+eOut)/max(0.01,1-eOut))
        fA=aIn*(1-eIn)*mB/(mA+mB); fB=aIn*(1-eIn)*mA/(mA+mB)
        pos=np.array([[-fA,0.,0.],[fB,0.,0.],[aOut*(1-eOut),0.,0.]])
        vel=np.array([[0.,vCI*(1+eIn)*mB/(mA+mB),0.],
                      [0.,-vCI*(1+eIn)*mA/(mA+mB),0.],
                      [0.,vCO*np.cos(incl),vCO*np.sin(incl)*sign]])
        masses=np.array([mA,mB,mC])
        pos,vel=com_correct(pos,vel,masses)
        ratio=ma01_ratio(pos,vel,masses)
        if (want_stable and ratio>=stable_min) or (not want_stable and ratio<=unstable_max):
            return pos,vel,masses,float(ratio),float(outer_period(pos,vel,masses))
    return None

# ── single trial: PNEP prediction + long ground truth + liberation track ─────
def run_trial(want_stable, rng, gt_max_t, track_liberation=True):
    ic = make_ic(want_stable, rng)
    if ic is None: return None
    pos, vel, masses, ma_ratio, T_out = ic
    E0 = total_energy(pos, vel, masses)

    t = 0.0; step = 0; nc = 0; H_buf = []; H_sum = 0.0
    d_lag1 = d_lag2 = None
    eject_t = None; eject_body = None; pred = None
    nIvl = []; last_node_t = 0.0
    # per-node record of each body's eccentricity about the other two
    ecc_t = []; ecc_vals = [[], [], []]
    peri_t = [[], [], []]          # outer-orbit periapsis times per body
    r_lag1 = None; r_lag2 = None

    while t < gt_max_t:
        pos, vel = leapfrog(pos, vel, masses, DT)
        t += DT; step += 1

        if eject_t is None and t > T_out:
            dists = np.linalg.norm(pos, axis=1)
            if np.any(dists > EJECT_R):
                eject_t = t; eject_body = int(np.argmax(dists))

        d_cur = np.array([np.linalg.norm(pos[1]-pos[0]),
                          np.linalg.norm(pos[2]-pos[1]),
                          np.linalg.norm(pos[2]-pos[0])])

        if is_local_min(d_cur, d_lag1, d_lag2):
            nc += 1
            H = compute_H(pos); H_buf.append(H); H_sum += H
            if len(H_buf) > 30: H_buf.pop(0)
            if last_node_t > 0: nIvl.append(t-last_node_t)
            last_node_t = t

            if track_liberation and eject_t is None:
                ecc_t.append(t)
                for k in range(3):
                    _, e_k, _ = outer_elements(pos, vel, masses, k)
                    ecc_vals[k].append(e_k)

            if pred is None and nc >= MIN_NODES:
                Havg = float(np.mean(H_buf))
                mDt = float(np.mean(nIvl[-10:])) if nIvl else 1.
                if Havg >= H_THRESH:
                    pred = {'stable': True, 'pred_t': t}
                else:
                    sl = float(np.polyfit(range(len(H_buf)), H_buf, 1)[0]) if len(H_buf)>=5 else -0.01
                    ns = max(1.,(Havg-0.1)/(-sl+1e-9)) if sl < 0 else 5.
                    pred = {'stable': False, 'pred_t': t,
                            'tLo': t+ns*mDt*0.5, 'tHi': t+ns*mDt*1.5}

        # outer-orbit periapsis counting (per body, vs CoM of other two)
        if track_liberation:
            r_cur = np.array([outer_elements(pos, vel, masses, k)[2] for k in range(3)]) \
                    if (nc > 0 and step % 40 == 0) else None
            if r_cur is not None:
                if r_lag1 is not None and r_lag2 is not None:
                    for k in range(3):
                        if r_lag2[k] > r_lag1[k] and r_cur[k] > r_lag1[k]:
                            peri_t[k].append(t)
                r_lag2 = r_lag1; r_lag1 = r_cur

        d_lag2 = d_lag1; d_lag1 = d_cur.copy()

        if eject_t is not None and pred is not None:
            break

    if pred is None:
        pred = {'stable': eject_t is None, 'pred_t': t}

    Ef = total_energy(pos, vel, masses)
    return dict(
        want_stable=want_stable,
        gt_stable=(eject_t is None),
        pnep_stable=pred['stable'],
        ma_ratio=ma_ratio, T_out=T_out,
        H_avg=H_sum/nc if nc else 0.5, n_nodes=nc,
        eject_t=eject_t, eject_body=eject_body,
        pred_t=pred.get('pred_t'),
        tLo=pred.get('tLo'), tHi=pred.get('tHi'),
        energy_drift=abs((Ef-E0)/(abs(E0)+1e-12)),
        ecc_t=ecc_t,
        ecc_eject=(ecc_vals[eject_body] if eject_body is not None else []),
        peri_eject=(peri_t[eject_body] if eject_body is not None else []),
    )

# ── liberation detection ─────────────────────────────────────────────────────
def find_liberation(ecc_t, ecc_vals, window=15, rise=0.05, min_lag=5):
    """Return (t_min_retro, t_detect_online). Retrospective = argmin of rolling
    mean. Online = first node where rolling mean exceeds running min by `rise`."""
    if len(ecc_vals) < window*2: return None, None
    e = np.asarray(ecc_vals, float); tt = np.asarray(ecc_t, float)
    roll = np.convolve(e, np.ones(window)/window, mode='valid')
    off = window//2
    t_roll = tt[off:off+len(roll)]
    t_retro = float(t_roll[int(np.argmin(roll))])
    run_min = roll[0]; i_min = 0; t_online = None
    for i in range(1, len(roll)):
        if roll[i] < run_min: run_min, i_min = roll[i], i
        elif roll[i] - run_min >= rise and (i - i_min) >= min_lag:
            t_online = float(t_roll[i_min]); break
    return t_retro, t_online

def orbits_between(peri_times, t0, t1):
    return int(sum(1 for p in peri_times if t0 < p <= t1))

# ── main ─────────────────────────────────────────────────────────────────────
def main(n_trials=100, seed=42, gt_max_t=3000.0):
    rng = np.random.default_rng(seed)
    rows = []; t0 = time.time()
    for i in range(n_trials):
        r = run_trial(i % 2 == 0, rng, gt_max_t)
        if r is None:
            print(f"  {i+1:3d}/{n_trials}  skipped (no valid IC)  "
                  f"[{time.time()-t0:.0f}s]", flush=True)
            continue
        rows.append(r)
        ok = sum(1 for x in rows if x['gt_stable'] == x['pnep_stable'])
        tag = "stable " if r['gt_stable'] else f"EJECT@{r['eject_t']:.0f}"
        hit = "ok " if r['gt_stable'] == r['pnep_stable'] else "MISS"
        print(f"  {i+1:3d}/{n_trials}  {tag:>12s}  pnep={'stable' if r['pnep_stable'] else 'unstable':<8s} "
              f"{hit}  running acc={ok/len(rows)*100:5.1f}%  [{time.time()-t0:.0f}s]", flush=True)

    TP = sum(1 for r in rows if r['gt_stable'] and r['pnep_stable'])
    FP = sum(1 for r in rows if not r['gt_stable'] and r['pnep_stable'])
    FN = sum(1 for r in rows if r['gt_stable'] and not r['pnep_stable'])
    TN = sum(1 for r in rows if not r['gt_stable'] and not r['pnep_stable'])
    n = len(rows)

    print("\n" + "="*64)
    print(f"  PART A — Re-validation   n={n}   ground truth window = {gt_max_t:.0f}t")
    print("="*64)
    print(f"  Accuracy                        : {(TP+TN)/n*100:.1f}%")
    print(f"  Ejections caught   TN/(TN+FP)   : {TN}/{TN+FP} = {TN/(TN+FP+1e-9)*100:.1f}%")
    print(f"  Stable correctly called TP/(TP+FN): {TP}/{TP+FN} = {TP/(TP+FN+1e-9)*100:.1f}%")
    print(f"  Median energy drift             : {np.median([r['energy_drift'] for r in rows])*100:.5f}%")
    print()
    print("  Confusion (positive class = STABLE):")
    print(f"    stable  called stable    TP = {TP}")
    print(f"    EJECTED called stable    FP = {FP}")
    print(f"    stable  called unstable  FN = {FN}")
    print(f"    EJECTED called unstable  TN = {TN}")

    # ── PART B ───────────────────────────────────────────────────────────────
    ej = [r for r in rows if r['eject_t'] is not None and len(r['ecc_eject']) > 40]
    dts, orbs, dts_on, orbs_on = [], [], [], []
    for r in ej:
        t_retro, t_online = find_liberation(r['ecc_t'], r['ecc_eject'])
        if t_retro is not None and t_retro < r['eject_t']:
            dts.append(r['eject_t'] - t_retro)
            orbs.append(orbits_between(r['peri_eject'], t_retro, r['eject_t']))
        if t_online is not None and t_online < r['eject_t']:
            dts_on.append(r['eject_t'] - t_online)
            orbs_on.append(orbits_between(r['peri_eject'], t_online, r['eject_t']))

    def cv(x):
        x = np.asarray(x, float)
        return float(np.std(x)/ (np.mean(x)+1e-12))

    print("\n" + "="*64)
    print("  PART B — Liberation lead time: raw time vs. outer-orbit count")
    print("="*64)
    for label, d, o in (("retrospective e-min", dts, orbs), ("online detector", dts_on, orbs_on)):
        if len(d) < 3:
            print(f"  {label}: too few systems ({len(d)})"); continue
        print(f"\n  {label}   n={len(d)}")
        print(f"    lead time  Δt     : median {np.median(d):8.1f}   "
              f"range {min(d):.1f}–{max(d):.1f}   CV = {cv(d):.3f}")
        print(f"    lead orbits N     : median {np.median(o):8.1f}   "
              f"range {min(o)}–{max(o)}        CV = {cv(o):.3f}")
        verdict = "ORBIT COUNT is tighter" if cv(o) < cv(d) else "raw time is tighter"
        print(f"    -> {verdict}  (lower CV wins)")

    json.dump([{k: v for k, v in r.items()
                if k not in ('ecc_t','ecc_eject','peri_eject')} for r in rows],
              open('pnep_revalidation_results.json','w'), indent=1)
    print("\n  wrote pnep_revalidation_results.json")

def _cli_args():
    """Parse argv, ignoring Jupyter/Colab's own flags (-f kernel.json)."""
    argv = [a for a in sys.argv[1:] if not a.startswith('-')]
    # in a notebook the remaining item is the kernel json path; drop non-numeric
    nums = []
    for a in argv:
        try:
            nums.append(float(a))
        except ValueError:
            pass
    n   = int(nums[0]) if len(nums) > 0 else 100
    sd  = int(nums[1]) if len(nums) > 1 else 42
    gt  = float(nums[2]) if len(nums) > 2 else 3000.0
    return n, sd, gt

if __name__ == "__main__":
    main(*_cli_args())
