"""
PNEP v13 — Predictive Node Event Protocol — 4-Body Extension
Extends v12 to hierarchical 4-body systems.

Two topologies explored:
  3+1 : hierarchical triple (bodies 0,1,2) with fourth outer body (body 3)
  2+2 : two inner binaries (0,1) and (2,3) orbiting each other

H signal extended to 6 pairwise distances for N=4.

For 3+1 systems, two hierarchy indices computed:
  H_inner : sigma2/(1+sigma2) of inner triple distances (d01, d02, d12)
  H_cross : sigma2/(1+sigma2) of distances involving body 3 (d03, d13, d23)
  H_full  : sigma2/(1+sigma2) of all 6 pairwise distances

Stable 3+1  -> H_inner ~ high, H_cross ~ high (4th body cleanly separated)
Unstable    -> H_inner or H_cross collapses
"""

import numpy as np
from dataclasses import dataclass
from typing import Optional
import time
import json

G       = 1.0
SOFT    = 0.01
DT      = 0.006
EJECT_R = 15.0
MIN_NODES = 10
H_THRESH  = 0.5

# ── Integrator (N-general) ────────────────────────────────────────────────────

def _acc(pos, masses):
    diff = pos[np.newaxis,:,:] - pos[:,np.newaxis,:]
    r2   = np.sum(diff**2, axis=2) + SOFT**2
    np.fill_diagonal(r2, 1.0)
    r3   = r2 * np.sqrt(r2)
    fac  = masses[np.newaxis,:] / r3
    np.fill_diagonal(fac, 0.0)
    return G * np.einsum('ij,ijk->ik', fac, diff)

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

# ── Orbital elements ──────────────────────────────────────────────────────────

def _orb(r1, v1, r2, v2, m1, m2):
    dr=r2-r1; dv=v2-v1; mu=G*(m1+m2)
    r=float(np.linalg.norm(dr))+1e-10
    v2s=float(np.dot(dv,dv))
    a=1./(2./r - v2s/mu)
    rdv=float(np.dot(dr,dv))
    ev=(v2s/mu-1./r)*dr - rdv/mu*dv
    e=float(np.clip(np.linalg.norm(ev),0,0.99))
    L=np.cross(dr,dv)
    return abs(float(a)), e, L

def ma01_ratio(pos, vel, masses):
    m0,m1,m2 = masses[0],masses[1],masses[2]
    ai,ei,Li = _orb(pos[0],vel[0],pos[1],vel[1],m0,m1)
    mI = m0+m1
    cp = (m0*pos[0]+m1*pos[1])/mI
    cv = (m0*vel[0]+m1*vel[1])/mI
    ao,eo,Lo = _orb(cp,cv,pos[2],vel[2],mI,m2)
    qO = m2/mI
    cosI = float(np.dot(Li,Lo))/(float(np.linalg.norm(Li))*float(np.linalg.norm(Lo))+1e-12)
    incl = float(np.arccos(np.clip(cosI,-1,1)))
    rhs  = 2.8*((1+qO)*(1+eo)/np.sqrt(max(0.01,1-eo)))**0.4*(1-0.3*incl/np.pi)
    return ao*(1-eo)/(ai+1e-9)/(rhs+1e-9)

def outer_period_3body(pos, vel, masses):
    m0,m1,m2 = masses[0],masses[1],masses[2]
    mI=m0+m1
    cp=(m0*pos[0]+m1*pos[1])/mI
    cv=(m0*vel[0]+m1*vel[1])/mI
    ao,eo,_ = _orb(cp,cv,pos[2],vel[2],mI,m2)
    return 2*np.pi*np.sqrt(ao**3/(G*(mI+m2)))

# ── 4-body H signal ───────────────────────────────────────────────────────────

def all_pairs_distances(pos):
    N = len(pos)
    d = {}
    for i in range(N):
        for j in range(i+1, N):
            d[(i,j)] = float(np.linalg.norm(pos[j]-pos[i]))
    return d

def compute_H_4body(pos, topology='3+1'):
    d = all_pairs_distances(pos)

    if topology == '3+1':
        d_inner = np.array([d[(0,1)], d[(0,2)], d[(1,2)]])
        s2_inner = float(np.var(d_inner))
        H_inner  = s2_inner / (1.0 + s2_inner)

        d_cross = np.array([d[(0,3)], d[(1,3)], d[(2,3)]])
        s2_cross = float(np.var(d_cross))
        H_cross  = s2_cross / (1.0 + s2_cross)

        d_full = np.array(list(d.values()))
        s2_full = float(np.var(d_full))
        H_full  = s2_full / (1.0 + s2_full)

        return H_inner, H_cross, H_full

    elif topology == '2+2':
        d_A = np.array([d[(0,1)], d[(0,2)], d[(0,3)], d[(1,2)], d[(1,3)]])
        s2_A = float(np.var(d_A))
        H_A  = s2_A / (1.0 + s2_A)

        d_B = np.array([d[(2,3)], d[(0,2)], d[(0,3)], d[(1,2)], d[(1,3)]])
        s2_B = float(np.var(d_B))
        H_B  = s2_B / (1.0 + s2_B)

        d_full = np.array(list(d.values()))
        s2_full = float(np.var(d_full))
        H_full  = s2_full / (1.0 + s2_full)

        return H_A, H_B, H_full

def compute_H_mass_weighted(pos, masses):
    """
    Mass-weighted hierarchy index — Alika Parks intuition.

    The pair connecting the heaviest and lightest body carries
    the most diagnostic information about hierarchy. That pair
    gets double weight in the variance. Any pair involving either
    extreme mass body gets 1.5x. Normal pairs get 1.0x.

    Hypothesis: sharpens H separation in 4-body systems where
    raw variance is compressed by the extra pairs.
    """
    d = all_pairs_distances(pos)
    i_max = int(np.argmax(masses))
    i_min = int(np.argmin(masses))

    d_arr = []; w_arr = []
    for (i,j), dist in d.items():
        if (i==i_max and j==i_min) or (i==i_min and j==i_max):
            w = 2.0   # heaviest-lightest pair: double weight
        elif i in [i_max,i_min] or j in [i_max,i_min]:
            w = 1.5   # any extreme-mass pair: 1.5x
        else:
            w = 1.0
        d_arr.append(dist)
        w_arr.append(w)

    d_arr  = np.array(d_arr)
    w_arr  = np.array(w_arr)
    w_mean = np.average(d_arr, weights=w_arr)
    s2w    = float(np.average((d_arr - w_mean)**2, weights=w_arr))
    return s2w / (1.0 + s2w)

def is_local_min_N(d_cur, d_lag1, d_lag2):
    if d_lag1 is None or d_lag2 is None:
        return False
    for key in d_cur:
        if d_lag2[key] > d_lag1[key] and d_cur[key] > d_lag1[key]:
            return True
    return False

# ── IC generators ─────────────────────────────────────────────────────────────

def make_ic_3plus1(want_stable, rng, stable_min=3.0, unstable_max=0.85, max_att=3000):
    """
    3+1 IC generator — targets MA01 ratio of inner triple (bodies 0,1,2).
    MA01 measures aMid*(1-eMid) / (aIn * rhs).
    We derive aMid directly from the target ratio, then place body 3 wide.
    """
    for _ in range(max_att):
        m0=0.7+rng.random()*0.9
        m1=0.6+rng.random()*0.9
        m2=0.15+rng.random()*0.8
        m3=0.08+rng.random()*0.4
        mAB=m0+m1; mABC=mAB+m2; mAll=mABC+m3

        aIn=0.3+rng.random()*0.3; eIn=rng.random()*0.15
        eMid=rng.random()*0.35
        incl_mid=abs(float(rng.standard_normal())*20.*np.pi/180.)

        # MA01 RHS
        qO=m2/mAB
        rhs=2.8*((1+qO)*(1+eMid)/np.sqrt(max(0.01,1-eMid)))**0.4*(1-0.3*incl_mid/np.pi)

        # Target and derive aMid periapsis
        if want_stable:
            target=stable_min+rng.random()*4.0
        else:
            target=0.3+rng.random()*(unstable_max-0.3)

        aMid_peri=target*aIn*rhs
        aMid=aMid_peri/max(0.01,1-eMid)
        if aMid<=aIn*1.5 or aMid>50.: continue

        # Body 3 wide outer orbit
        eOut=rng.random()*0.30; incl_out=abs(float(rng.standard_normal())*15.*np.pi/180.)
        sign=1 if rng.random()>0.5 else -1
        Pr_out=float(np.exp(rng.uniform(np.log(5),np.log(25))))
        aOut=aMid*Pr_out**(2./3.)
        if aOut<=aMid*1.5 or aOut>500.: continue

        vCI=np.sqrt(G*mAB/aIn)
        vCM=np.sqrt(G*mABC/aMid)*np.sqrt((1+eMid)/max(0.01,1-eMid))
        vCO=np.sqrt(G*mAll/aOut)*np.sqrt((1+eOut)/max(0.01,1-eOut))
        fA=aIn*(1-eIn)*m1/mAB; fB=aIn*(1-eIn)*m0/mAB

        pos=np.array([[-fA,0.,0.],[fB,0.,0.],[aMid*(1-eMid),0.,0.],[aOut*(1-eOut),0.,0.]])
        vel=np.array([[0.,vCI*(1+eIn)*m1/mAB,0.],[0.,-vCI*(1+eIn)*m0/mAB,0.],
                      [0.,vCM*np.cos(incl_mid),vCM*np.sin(incl_mid)],
                      [0.,vCO*np.cos(incl_out),vCO*np.sin(incl_out)*sign]])
        masses=np.array([m0,m1,m2,m3])
        pos,vel=com_correct(pos,vel,masses)
        ratio=ma01_ratio(pos,vel,masses)
        ok=(want_stable and ratio>=stable_min) or (not want_stable and ratio<=unstable_max)
        if ok:
            T_out=outer_period_3body(pos,vel,masses)
            return pos,vel,masses,float(ratio),float(T_out),'3+1'
    return None


def make_ic_2plus2(want_stable, rng, max_att=2000):
    for _ in range(max_att):
        m0=0.6+rng.random()*0.9
        m1=0.5+rng.random()*0.9
        m2=0.6+rng.random()*0.9
        m3=0.5+rng.random()*0.9

        aA=0.2+rng.random()*0.3; eA=rng.random()*0.15
        aB=0.2+rng.random()*0.3; eB=rng.random()*0.15

        Pr  = float(np.exp(rng.uniform(np.log(5), np.log(80))))
        aOut= max(aA,aB) * Pr**(2./3.)
        eOut= rng.random()*0.35
        incl= abs(float(rng.standard_normal())*20.*np.pi/180.)
        sign= 1 if rng.random()>0.5 else -1

        mA=m0+m1; mB=m2+m3

        vA  = np.sqrt(G*mA/aA)
        vB  = np.sqrt(G*mB/aB)
        vOut= np.sqrt(G*(mA+mB)/aOut)*np.sqrt((1+eOut)/max(0.01,1-eOut))

        fA0=aA*(1-eA)*m1/mA; fA1=aA*(1-eA)*m0/mA
        fB0=aB*(1-eB)*m3/mB; fB1=aB*(1-eB)*m2/mB
        dOut=aOut*(1-eOut)

        pos = np.array([
            [-fA0, 0., 0.],
            [ fA1, 0., 0.],
            [dOut-fB0, 0., 0.],
            [dOut+fB1, 0., 0.]
        ])
        vel = np.array([
            [0.,  vA*(1+eA)*m1/mA, 0.],
            [0., -vA*(1+eA)*m0/mA, 0.],
            [0.,  vOut*np.cos(incl)+vB*(1+eB)*m3/mB, vOut*np.sin(incl)*sign],
            [0.,  vOut*np.cos(incl)-vB*(1+eB)*m2/mB, vOut*np.sin(incl)*sign]
        ])
        masses = np.array([m0,m1,m2,m3])
        pos,vel = com_correct(pos,vel,masses)

        rOut = float(np.linalg.norm(pos[2]-pos[0]))
        rIn  = max(float(np.linalg.norm(pos[1]-pos[0])),
                   float(np.linalg.norm(pos[3]-pos[2])))
        ratio = rOut/(rIn+1e-9)

        ok = (want_stable and ratio>=8.0) or (not want_stable and ratio<=2.5)
        if ok:
            T_out = 2*np.pi*np.sqrt(aOut**3/(G*(mA+mB)))
            return pos, vel, masses, float(ratio), float(T_out), '2+2'
    return None

# ── Trial ─────────────────────────────────────────────────────────────────────

@dataclass
class TrialResult4:
    topology:     str
    want_stable:  bool
    gt_stable:    bool
    pnep_stable:  bool
    correct:      bool
    ma_ratio:     float
    H1_avg:       float
    H2_avg:       float
    H_full_avg:   float
    Hw_avg:       float   # mass-weighted H (Parks intuition)
    n_nodes:      int
    pre_eject:    bool
    eject_t:      Optional[float]
    energy_drift: float

def run_trial_4body(topology, want_stable, rng,
                    H_thresh=H_THRESH, min_nodes=MIN_NODES,
                    eject_r=EJECT_R, dt=DT):
    if topology == '3+1':
        ic = make_ic_3plus1(want_stable, rng)
    else:
        ic = make_ic_2plus2(want_stable, rng)

    if ic is None:
        return None

    pos, vel, masses, ma_ratio, T_out, top = ic
    max_t = min(T_out*5 if want_stable else max(300., T_out*8), 1500.)
    E0 = total_energy(pos, vel, masses)

    t=0.; nc=0
    H1_sum=H2_sum=Hf_sum=Hw_sum=0.
    H1_buf=[]; H2_buf=[]; Hf_buf=[]; Hw_buf=[]
    d_cur=None; d_lag1=None; d_lag2=None
    actual_eject=None; pred=None
    nIvl=[]; last_node_t=0.

    while t < max_t:
        pos, vel = leapfrog(pos, vel, masses, dt)
        t += dt

        if actual_eject is None and t > T_out:
            if np.any(np.linalg.norm(pos, axis=1) > eject_r):
                actual_eject = t

        d_cur = all_pairs_distances(pos)

        if is_local_min_N(d_cur, d_lag1, d_lag2):
            nc += 1
            H1, H2, Hf = compute_H_4body(pos, topology)
            Hw = compute_H_mass_weighted(pos, masses)
            H1_buf.append(H1); H2_buf.append(H2)
            Hf_buf.append(Hf); Hw_buf.append(Hw)
            if len(H1_buf) > 30:
                H1_buf.pop(0); H2_buf.pop(0)
                Hf_buf.pop(0); Hw_buf.pop(0)
            H1_sum+=H1; H2_sum+=H2; Hf_sum+=Hf; Hw_sum+=Hw

            if last_node_t > 0:
                nIvl.append(t - last_node_t)
            last_node_t = t

            if pred is None and nc >= min_nodes:
                H1avg = float(np.mean(H1_buf))
                H2avg = float(np.mean(H2_buf))
                Hwavg = float(np.mean(Hw_buf))

                # Use mass-weighted H as primary signal alongside H1/H2
                if H1avg >= H_thresh and H2avg >= H_thresh and Hwavg >= H_thresh:
                    pred = {'stable': True}
                else:
                    H_signal = min(H1avg, H2avg, Hwavg)
                    sl = float(np.polyfit(range(len(Hw_buf)), Hw_buf, 1)[0]) \
                         if len(Hw_buf)>=5 else -0.01
                    mDt = float(np.mean(nIvl[-10:])) if nIvl else 1.
                    ns  = max(1., (H_signal-0.1)/(-sl+1e-9)) if sl<0 else 5.
                    pred = {'stable': False,
                            'tLo': t+ns*mDt*0.5,
                            'tHi': t+ns*mDt*1.5,
                            'pred_t': t}

        d_lag2 = d_lag1
        d_lag1 = dict(d_cur) if d_cur else None

        if actual_eject is not None and pred is not None:
            break

    if pred is None:
        pred = {'stable': actual_eject is None}

    Ef    = total_energy(pos, vel, masses)
    drift = abs((Ef-E0)/(abs(E0)+1e-12))
    H1_avg = H1_sum/nc if nc>0 else 0.5
    H2_avg = H2_sum/nc if nc>0 else 0.5
    Hf_avg = Hf_sum/nc if nc>0 else 0.5
    Hw_avg = Hw_sum/nc if nc>0 else 0.5
    gt = actual_eject is None
    ps = pred['stable']
    pre = (not gt and not ps and 'tLo' in pred
           and pred['pred_t'] < (actual_eject or 1e9))

    return TrialResult4(
        topology=topology, want_stable=want_stable,
        gt_stable=gt, pnep_stable=ps, correct=(ps==gt),
        ma_ratio=ma_ratio, H1_avg=H1_avg, H2_avg=H2_avg,
        H_full_avg=Hf_avg, Hw_avg=Hw_avg, n_nodes=nc, pre_eject=pre,
        eject_t=actual_eject, energy_drift=drift)

# ── Batch ─────────────────────────────────────────────────────────────────────

def run_batch_4body(n=100, seed=42, verbose=True):
    rng = np.random.default_rng(seed)
    results = {'3+1': [], '2+2': []}
    skipped = 0
    t0 = time.time()

    for i in range(n):
        topology    = '3+1' if (i % 4 < 2) else '2+2'
        want_stable = (i % 2 == 0)
        r = run_trial_4body(topology, want_stable, rng)
        if r is None:
            skipped += 1
            continue
        results[topology].append(r)

        if verbose and (i+1) % 10 == 0:
            for top in ['3+1','2+2']:
                rs = results[top]
                if not rs: continue
                acc = sum(x.correct for x in rs)/len(rs)*100
                H1s = [x.H1_avg for x in rs if x.gt_stable]
                H1u = [x.H1_avg for x in rs if not x.gt_stable]
                print(f"  [{top}] n={len(rs):3d}  acc={acc:.1f}%  "
                      f"H1_s={np.mean(H1s) if H1s else 0:.3f}  "
                      f"H1_u={np.mean(H1u) if H1u else 0:.3f}  "
                      f"[{time.time()-t0:.0f}s]")

    summaries = {}
    for top in ['3+1','2+2']:
        rs = results[top]
        if not rs: continue
        TP=FP=FN=TN=0
        for r in rs:
            if   r.gt_stable and     r.pnep_stable: TP+=1
            elif not r.gt_stable and r.pnep_stable:  FP+=1
            elif r.gt_stable and not r.pnep_stable:  FN+=1
            else:                                    TN+=1
        tot=len(rs)
        sp=TN/(TN+FP+1e-9)*100; se=TP/(TP+FN+1e-9)*100
        H1s=[r.H1_avg for r in rs if r.gt_stable]
        H1u=[r.H1_avg for r in rs if not r.gt_stable]
        H2s=[r.H2_avg for r in rs if r.gt_stable]
        H2u=[r.H2_avg for r in rs if not r.gt_stable]
        Hws=[r.Hw_avg for r in rs if r.gt_stable]
        Hwu=[r.Hw_avg for r in rs if not r.gt_stable]
        pre=[r for r in rs if r.pre_eject]
        ej =[r for r in rs if not r.gt_stable]
        drifts=[r.energy_drift for r in rs]
        summaries[top] = dict(
            n=tot, TP=TP, FP=FP, FN=FN, TN=TN,
            acc=(TP+TN)/tot*100, spec=sp, sens=se,
            f1=2*sp*se/(sp+se+1e-9),
            H1_s=np.mean(H1s) if H1s else float('nan'),
            H1_u=np.mean(H1u) if H1u else float('nan'),
            H1_sep=(np.mean(H1s)-np.mean(H1u)) if H1s and H1u else 0.,
            H2_s=np.mean(H2s) if H2s else float('nan'),
            H2_u=np.mean(H2u) if H2u else float('nan'),
            Hw_s=np.mean(Hws) if Hws else float('nan'),
            Hw_u=np.mean(Hwu) if Hwu else float('nan'),
            Hw_sep=(np.mean(Hws)-np.mean(Hwu)) if Hws and Hwu else 0.,
            pre_pct=len(pre)/(len(ej)+1e-9)*100,
            med_drift=float(np.median(drifts))*100,
        )
    return summaries, results, skipped

def print_summary_4body(summaries, skipped):
    print("\n" + "="*65)
    print("  PNEP v13 — 4-Body Extension")
    print("="*65)
    print(f"  Skipped (IC failures): {skipped}")
    for top, s in summaries.items():
        label = '3+1 hierarchical' if top=='3+1' else '2+2 binary pair'
        print(f"\n  ── Topology: {top}  ({label}) ──────────────")
        print(f"  Trials        : {s['n']}")
        print(f"  Accuracy      : {s['acc']:.1f}%")
        print(f"  Specificity   : {s['spec']:.1f}%")
        print(f"  Sensitivity   : {s['sens']:.1f}%")
        print(f"  F1            : {s['f1']:.1f}%")
        print(f"  TP={s['TP']} FP={s['FP']} FN={s['FN']} TN={s['TN']}")
        print(f"  ── Standard H ───────────────────────────")
        print(f"  H1 stable     : {s['H1_s']:.4f}")
        print(f"  H1 unstable   : {s['H1_u']:.4f}")
        print(f"  H1 separation : {s['H1_sep']:.4f}")
        print(f"  H2 stable     : {s['H2_s']:.4f}")
        print(f"  H2 unstable   : {s['H2_u']:.4f}")
        print(f"  ── Mass-Weighted H (Parks) ───────────────")
        print(f"  Hw stable     : {s['Hw_s']:.4f}")
        print(f"  Hw unstable   : {s['Hw_u']:.4f}")
        print(f"  Hw separation : {s['Hw_sep']:.4f}  {'<< IMPROVED' if s['Hw_sep'] > s['H1_sep'] else ''}")
        print(f"  ─────────────────────────────────────────")
        print(f"  Pre-eject rate: {s['pre_pct']:.1f}%")
        print(f"  Energy drift  : {s['med_drift']:.5f}%  median")
    print("="*65)

if __name__ == "__main__":
    import sys
    n    = int(sys.argv[1]) if len(sys.argv)>1 else 100
    seed = int(sys.argv[2]) if len(sys.argv)>2 else 42
    print(f"PNEP v13  n={n}  seed={seed}  dt={DT}  4-body extension")
    summaries, results, skipped = run_batch_4body(n=n, seed=seed, verbose=True)
    print_summary_4body(summaries, skipped)
    out = {top: s for top,s in summaries.items()}
    with open('pnep_v13_results.json','w') as f:
        json.dump(out, f, indent=2)
    print("\nResults saved to pnep_v13_results.json")


# ── Weight sweep analysis ─────────────────────────────────────────────────────

def run_weight_sweep(n=60, seed=42):
    """
    Run batch collecting H values at multiple weight combos.
    Compare separation for each to find optimal weighting.
    """
    import numpy as np
    rng = np.random.default_rng(seed)
    combos = {
        'standard':  (None, None),
        'w2.0/1.5':  (2.0, 1.5),
        'w3.0/1.8':  (3.0, 1.8),
        'w4.0/2.0':  (4.0, 2.0),
        'w5.0/2.5':  (5.0, 2.5),
        'w6.0/3.0':  (6.0, 3.0),
    }
    # Accumulate stable/unstable H for each combo
    buckets = {k: {'s': [], 'u': []} for k in combos}
    skipped = 0

    for i in range(n):
        want_stable = (i % 2 == 0)
        ic = make_ic_2plus2(want_stable, rng)
        if ic is None:
            skipped += 1
            continue
        pos, vel, masses, ma_ratio, T_out, top = ic
        max_t = min(T_out*5 if want_stable else max(300., T_out*8), 800.)
        E0 = total_energy(pos, vel, masses)

        t=0.; nc=0
        h_accum = {k: 0.0 for k in combos}
        d_lag1=d_lag2=None
        actual_eject=None

        while t < max_t:
            pos, vel = leapfrog(pos, vel, masses)
            t += DT
            if actual_eject is None and t > T_out:
                if np.any(np.linalg.norm(pos, axis=1) > EJECT_R):
                    actual_eject = t
            d_cur = all_pairs_distances(pos)
            if is_local_min_N(d_cur, d_lag1, d_lag2):
                nc += 1
                # standard H (full unweighted)
                d_full = np.array(list(d_cur.values()))
                s2 = float(np.var(d_full))
                h_accum['standard'] += s2/(1+s2)
                # weighted combos
                for label, (we, wm) in combos.items():
                    if we is None: continue
                    h_accum[label] += compute_H_mass_weighted(pos, masses, we, wm)
            d_lag2=d_lag1; d_lag1=dict(d_cur)
            if actual_eject and nc >= MIN_NODES:
                break

        if nc == 0: continue
        gt_stable = actual_eject is None
        key = 's' if gt_stable else 'u'
        for label in combos:
            buckets[label][key].append(h_accum[label]/nc)

    print(f"\n{'='*60}")
    print(f"  Mass Weight Sweep — 2+2 topology  (skipped={skipped})")
    print(f"{'='*60}")
    print(f"  {'Label':<14} {'H_stable':>9} {'H_unstable':>11} {'Delta':>8}  {'Best?'}")
    print(f"  {'-'*56}")
    best_sep = 0; best_label = ''
    for label, bkt in buckets.items():
        if not bkt['s'] or not bkt['u']: continue
        hs = float(np.mean(bkt['s']))
        hu = float(np.mean(bkt['u']))
        sep = hs - hu
        flag = ''
        if sep > best_sep:
            best_sep = sep; best_label = label; flag = '<< BEST'
        print(f"  {label:<14} {hs:>9.4f} {hu:>11.4f} {sep:>8.4f}  {flag}")
    print(f"{'='*60}")
    print(f"  Optimal: {best_label}  separation={best_sep:.4f}")

if __name__ == "__main__" and len(__import__('sys').argv) > 1 and __import__('sys').argv[1] == 'sweep':
    run_weight_sweep(n=80, seed=42)
