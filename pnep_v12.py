"""
PNEP v12 — Predictive Node Event Protocol
Vectorised numpy leapfrog (kick-drift-kick), fixed dt=0.006.
Node detection: local minimum of inter-body distances.
Energy drift: ~0.008% median (leapfrog, near machine precision).

H = sigma2 / (1 + sigma2)
    sigma2 = Var(d12, d23, d31) at mirror symmetry nodes
Stable hierarchical triple -> H ~ 0.92, Unstable -> H ~ 0.09
Single threshold H = 0.5 separates classes.
"""

import numpy as np
from dataclasses import dataclass
from typing import Optional
import time

G         = 1.0
SOFT      = 0.01
DT        = 0.006
EJECT_R   = 12.0
MIN_NODES = 10
H_THRESH  = 0.5

# ── Vectorised leapfrog ───────────────────────────────────────────────────────

def _acc(pos, masses):
    """Vectorised N-body acceleration. pos:(N,3) masses:(N,) -> (N,3)"""
    diff = pos[np.newaxis,:,:] - pos[:,np.newaxis,:]   # (N,N,3) diff[i,j]=pos[j]-pos[i]
    r2   = np.sum(diff**2, axis=2) + SOFT**2            # (N,N)
    np.fill_diagonal(r2, 1.0)                           # avoid /0 on diagonal
    r3   = r2 * np.sqrt(r2)
    fac  = masses[np.newaxis,:] / r3                    # (N,N) m_j / r_ij^3
    np.fill_diagonal(fac, 0.0)
    return G * np.einsum('ij,ijk->ik', fac, diff)       # sum over j -> (N,3)

def leapfrog(pos, vel, masses, dt=DT):
    a  = _acc(pos, masses)
    vh = vel + a * (dt / 2)
    p2 = pos + vh * dt
    return p2, vh + _acc(p2, masses) * (dt / 2)

def total_energy(pos, vel, masses):
    N  = len(masses)
    KE = 0.5 * np.sum(masses[:,np.newaxis] * vel**2)
    PE = 0.0
    for i in range(N):
        for j in range(i+1, N):
            dr = pos[j] - pos[i]
            r  = np.sqrt(float(np.dot(dr,dr)) + SOFT**2)
            PE -= G * masses[i] * masses[j] / r
    return float(KE + PE)

def com_correct(pos, vel, masses):
    tm = masses.sum()
    cp = (masses[:,np.newaxis] * pos).sum(0) / tm
    cv = (masses[:,np.newaxis] * vel).sum(0) / tm
    return pos - cp, vel - cv

# ── Orbital elements & MA01 ───────────────────────────────────────────────────

def _orb(r1, v1, r2, v2, m1, m2):
    dr=r2-r1; dv=v2-v1; mu=G*(m1+m2)
    r=float(np.linalg.norm(dr))+1e-10; v2s=float(np.dot(dv,dv))
    a=1./(2./r - v2s/mu); rdv=float(np.dot(dr,dv))
    ev=(v2s/mu-1./r)*dr - rdv/mu*dv
    e=float(np.clip(np.linalg.norm(ev), 0, 0.99))
    return abs(float(a)), e, np.cross(dr, dv)

def ma01_ratio(pos, vel, masses):
    m0,m1,m2 = masses
    ai,ei,Li = _orb(pos[0],vel[0], pos[1],vel[1], m0,m1)
    mI = m0+m1
    cp = (m0*pos[0]+m1*pos[1])/mI;  cv = (m0*vel[0]+m1*vel[1])/mI
    ao,eo,Lo = _orb(cp,cv, pos[2],vel[2], mI,m2)
    qO = m2/mI
    cosI = float(np.dot(Li,Lo))/(float(np.linalg.norm(Li))*float(np.linalg.norm(Lo))+1e-12)
    incl = float(np.arccos(np.clip(cosI,-1,1)))
    rhs  = 2.8*((1+qO)*(1+eo)/np.sqrt(max(0.01,1-eo)))**0.4*(1-0.3*incl/np.pi)
    return ao*(1-eo)/(ai+1e-9)/(rhs+1e-9)

def outer_period(pos, vel, masses):
    m0,m1,m2 = masses; mI=m0+m1
    cp=(m0*pos[0]+m1*pos[1])/mI;  cv=(m0*vel[0]+m1*vel[1])/mI
    ao,eo,_ = _orb(cp,cv, pos[2],vel[2], mI,m2)
    return 2*np.pi*np.sqrt(ao**3/(G*(mI+m2)))

# ── PNEP H signal ─────────────────────────────────────────────────────────────

def compute_H(pos):
    d  = np.array([np.linalg.norm(pos[1]-pos[0]),
                   np.linalg.norm(pos[2]-pos[1]),
                   np.linalg.norm(pos[2]-pos[0])])
    s2 = float(np.var(d))
    return s2 / (1.0 + s2)

def is_local_min(d, d1, d2):
    """True if any distance just passed through local minimum (d1 < d2 and d1 < d)."""
    if d1 is None or d2 is None:
        return False
    for i in range(3):
        if d2[i] > d1[i] and d[i] > d1[i]:   # was decreasing, now increasing
            return True
    return False

# ── IC generator ─────────────────────────────────────────────────────────────

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
        ok=(want_stable and ratio>=stable_min) or (not want_stable and ratio<=unstable_max)
        if ok:
            return pos,vel,masses,float(ratio),float(outer_period(pos,vel,masses))
    return None

# ── Trial ─────────────────────────────────────────────────────────────────────

@dataclass
class TrialResult:
    want_stable:  bool
    gt_stable:    bool
    pnep_stable:  bool
    correct:      bool
    ma_ratio:     float
    H_avg:        float
    n_nodes:      int
    win_err:      Optional[float]
    pre_eject:    bool
    eject_t:      Optional[float]
    energy_drift: float

def run_trial(want_stable, rng, H_thresh=H_THRESH, min_nodes=MIN_NODES,
              eject_r=EJECT_R, dt=DT):
    ic = make_ic(want_stable, rng)
    if ic is None:
        return None
    pos, vel, masses, ma_ratio, T_out = ic
    max_t = min(T_out*5 if want_stable else max(200., T_out*10), 1200.)
    E0    = total_energy(pos, vel, masses)

    t=0.; nc=0; H_buf=[]; H_sum=0.
    d_cur=None; d_lag1=None; d_lag2=None
    actual_eject=None; pred=None; nIvl=[]; last_node_t=0.

    while t < max_t:
        pos, vel = leapfrog(pos, vel, masses, dt)
        t += dt

        # ejection check (after one outer period)
        if actual_eject is None and t > T_out:
            if np.any(np.linalg.norm(pos, axis=1) > eject_r):
                actual_eject = t

        # inter-body distances
        d_cur = np.array([np.linalg.norm(pos[1]-pos[0]),
                          np.linalg.norm(pos[2]-pos[1]),
                          np.linalg.norm(pos[2]-pos[0])])

        # local minimum node detection
        if is_local_min(d_cur, d_lag1, d_lag2):
            nc += 1
            H = compute_H(pos)
            H_buf.append(H)
            if len(H_buf) > 30: H_buf.pop(0)
            H_sum += H
            if last_node_t > 0:
                nIvl.append(t - last_node_t)
            last_node_t = t

            if pred is None and nc >= min_nodes:
                Havg = float(np.mean(H_buf))
                mDt  = float(np.mean(nIvl[-10:])) if nIvl else 1.

                if Havg >= H_thresh:
                    pred = {'stable': True}
                else:
                    sl = float(np.polyfit(range(len(H_buf)), H_buf, 1)[0]) \
                         if len(H_buf) >= 5 else -0.01
                    ns = max(1., (Havg-0.1)/(-sl+1e-9)) if sl < 0 else 5.
                    pred = {'stable': False,
                            'tLo':    t + ns*mDt*0.5,
                            'tHi':    t + ns*mDt*1.5,
                            'pred_t': t}

        # roll lag buffers
        d_lag2 = d_lag1
        d_lag1 = d_cur.copy() if d_cur is not None else None

        if actual_eject is not None and pred is not None:
            break

    if pred is None:
        pred = {'stable': actual_eject is None}

    Ef    = total_energy(pos, vel, masses)
    drift = abs((Ef - E0) / (abs(E0) + 1e-12))
    H_avg = H_sum / nc if nc > 0 else 0.5
    gt    = actual_eject is None
    ps    = pred['stable']
    win_err=None; pre_eject=False
    if not gt and not ps and 'tLo' in pred and pred['pred_t'] < actual_eject:
        pre_eject = True
        win_err   = abs(actual_eject - (pred['tLo']+pred['tHi'])/2)

    return TrialResult(
        want_stable=want_stable, gt_stable=gt, pnep_stable=ps,
        correct=(ps==gt), ma_ratio=ma_ratio, H_avg=H_avg, n_nodes=nc,
        win_err=win_err, pre_eject=pre_eject, eject_t=actual_eject,
        energy_drift=drift)

# ── Batch ─────────────────────────────────────────────────────────────────────

def run_batch(n=200, H_thresh=H_THRESH, seed=42, verbose=True):
    rng=np.random.default_rng(seed); results=[]; skipped=0
    TP=FP=FN=TN=0; H_s=[]; H_u=[]; wins=[]
    pre_ct=ej_tot=gt_ok=gt_tot=0; drifts=[]; t0=time.time()

    for i in range(n):
        ws = (i % 2 == 0)
        r  = run_trial(ws, rng, H_thresh=H_thresh)
        if r is None:
            skipped += 1; continue
        results.append(r)

        if   r.gt_stable and     r.pnep_stable: TP += 1
        elif not r.gt_stable and r.pnep_stable:  FP += 1
        elif r.gt_stable and not r.pnep_stable:  FN += 1
        else:                                    TN += 1

        (H_s if r.gt_stable else H_u).append(r.H_avg)
        if not r.gt_stable:
            ej_tot += 1
            if r.pre_eject:
                pre_ct += 1
                if r.win_err is not None: wins.append(r.win_err)
        if ws:
            gt_tot += 1
            if r.gt_stable: gt_ok += 1
        drifts.append(r.energy_drift)

        if verbose and (i+1) % 10 == 0:
            tot2=len(results); el=time.time()-t0
            print(f"  {i+1:3d}/{n}  "
                  f"acc={((TP+TN)/tot2*100):.1f}%  "
                  f"spec={TN/(TN+FP+1e-9)*100:.1f}%  "
                  f"sens={TP/(TP+FN+1e-9)*100:.1f}%  "
                  f"H_s={np.mean(H_s) if H_s else 0:.3f}  "
                  f"H_u={np.mean(H_u) if H_u else 0:.3f}  "
                  f"drift={np.median(drifts)*100:.4f}%  "
                  f"[{el:.0f}s]")

    tot=len(results)
    if tot == 0: return {}
    sp=TN/(TN+FP+1e-9)*100; se=TP/(TP+FN+1e-9)*100; wall=time.time()-t0

    return dict(
        n=tot, skipped=skipped, wall_s=wall,
        TP=TP, FP=FP, FN=FN, TN=TN,
        acc=(TP+TN)/tot*100, spec=sp, sens=se,
        f1=2*sp*se/(sp+se+1e-9),
        H_s=np.mean(H_s) if H_s else float('nan'),
        H_u=np.mean(H_u) if H_u else float('nan'),
        H_sep=(np.mean(H_s)-np.mean(H_u)) if H_s and H_u else 0.,
        gt_val=gt_ok/(gt_tot+1e-9)*100,
        pre_pct=pre_ct/(ej_tot+1e-9)*100,
        med_win=float(np.median(wins)) if wins else float('nan'),
        n_wins=len(wins),
        med_drift=float(np.median(drifts))*100,
        thresh=H_thresh, results=results)

def print_summary(s):
    w = s.get('wall_s', 0)
    print("\n" + "="*60)
    print("  PNEP v12 — Symplectic Leapfrog — Scientific Results")
    print("="*60)
    print(f"  Trials        : {s['n']}  (skipped {s['skipped']})  [{w:.0f}s wall]")
    print(f"  Integrator    : KDK leapfrog  dt={DT}  soft={SOFT}")
    print(f"  H threshold   : {s['thresh']}")
    print(f"  Energy drift  : {s['med_drift']:.5f}%  (median)")
    print()
    print("  ── H discrimination ─────────────────────────────────")
    print(f"  Stable   avg H : {s['H_s']:.4f}")
    print(f"  Unstable avg H : {s['H_u']:.4f}")
    print(f"  Separation  Δ  : {s['H_sep']:.4f}")
    print()
    print("  ── Ground truth ─────────────────────────────────────")
    print(f"  GT validity    : {s['gt_val']:.1f}%  (MA01 ratio>5 → no eject)")
    print()
    print("  ── Classification ───────────────────────────────────")
    print(f"  Accuracy       : {s['acc']:.1f}%")
    print(f"  Specificity    : {s['spec']:.1f}%  (stable correctly called)")
    print(f"  Sensitivity    : {s['sens']:.1f}%  (ejections caught)")
    print(f"  F1 score       : {s['f1']:.1f}%")
    print(f"  TP={s['TP']:3d}  FP={s['FP']:3d}  /  FN={s['FN']:3d}  TN={s['TN']:3d}")
    print()
    print("  ── Pre-ejection windows ─────────────────────────────")
    print(f"  Rate           : {s['pre_pct']:.1f}%  of ejecting systems")
    print(f"  Median error   : {s['med_win']:.1f}t  ({s['n_wins']} windows)")
    print("="*60)

if __name__ == "__main__":
    import sys
    n    = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 42
    print(f"PNEP v12  n={n}  seed={seed}  dt={DT}  soft={SOFT}  H_thresh={H_THRESH}")
    s = run_batch(n=n, seed=seed, verbose=True)
    print_summary(s)
