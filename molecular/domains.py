"""Two tests:
 A. TURNING ANGLE - the scale-free 'sharp angled turn', not Frenet curvature.
 B. H DOMAINS - do similar-H triples form connected spatial regions?
    And do those domains co-move?
"""
import numpy as np, pickle, time, sys, itertools
sys.path.insert(0,'/home/claude/md3')
from cluster import (N, DT, icosahedron13, acc_and_U, H_all, pair_d,
                     TRI, TRI_P, PRS, lindemann, n_evaporated)
NT=len(TRI)

def turn_angle(v_old, v_new):
    """angle between successive velocity vectors, per atom. scale-free."""
    a=np.einsum('ij,ij->i', v_old, v_new)
    b=np.linalg.norm(v_old,axis=1)*np.linalg.norm(v_new,axis=1)+1e-12
    return np.arccos(np.clip(a/b,-1,1))

def centroids(pos): return pos[TRI].mean(axis=1)

def domain_score(h, cen, q=0.75, nperm=40, rng=None, top=True):
    """Are similar-H triples spatially clustered?
    Selects the top (or bottom) quantile by H and compares their mean
    pairwise centroid distance to random subsets of equal size.
    NEGATIVE z = clustered (closer together than chance)."""
    thr=np.quantile(h, q if top else 1-q)
    sel=np.where(h>=thr)[0] if top else np.where(h<=thr)[0]
    k=len(sel)
    if k<5 or k>NT-5: return np.nan, k
    def mpd(idx):
        c=cen[idx]
        d=np.linalg.norm(c[:,None,:]-c[None,:,:],axis=2)
        iu=np.triu_indices(len(idx),1)
        return d[iu].mean()
    obs=mpd(sel)
    null=np.array([mpd(rng.choice(NT,k,replace=False)) for _ in range(nperm)])
    return float((obs-null.mean())/(null.std()+1e-12)), k

def atom_concentration(h, q=0.75):
    """Of the selected high-H triples, what fraction of all their atom-slots
    are taken by the single most-common atom?  High => one stray atom is
    driving the selection (artifact), not a spatial domain."""
    sel=np.where(h>=np.quantile(h,q))[0]
    if len(sel)<5: return np.nan
    cnt=np.bincount(TRI[sel].ravel(), minlength=N)
    return float(cnt.max()/cnt.sum())

def domain_score_disjoint(h, cen, q=0.75, nperm=40, rng=None):
    """Same as domain_score but restricted to a mutually DISJOINT subset of
    the selected triples - removes any shared-atom geometric coupling."""
    sel=np.where(h>=np.quantile(h,q))[0]
    if len(sel)<6: return np.nan
    order=sel[np.argsort(-h[sel])]
    used=set(); keep=[]
    for ti in order:
        a=set(TRI[ti])
        if used & a: continue
        used |= a; keep.append(ti)
    k=len(keep)
    if k<3: return np.nan
    def mpd(idx):
        c=cen[idx]; d=np.linalg.norm(c[:,None,:]-c[None,:,:],axis=2)
        iu=np.triu_indices(len(idx),1); return d[iu].mean()
    obs=mpd(np.array(keep))
    null=[]
    for _ in range(nperm):
        u=set(); kp=[]
        for ti in rng.permutation(NT):
            a=set(TRI[ti])
            if u & a: continue
            u|=a; kp.append(ti)
            if len(kp)==k: break
        if len(kp)==k: null.append(mpd(np.array(kp)))
    if len(null)<10: return np.nan
    null=np.array(null)
    return float((obs-null.mean())/(null.std()+1e-12))

def run(T0, seed, nsteps=120000, stride=40, burn=20000, lag=25):
    rng=np.random.default_rng(seed)
    pos=icosahedron13()+rng.normal(0,0.01,(N,3))
    vel=rng.normal(0,np.sqrt(T0),(N,3)); vel-=vel.mean(0)
    acc,_=acc_and_U(pos)
    vbuf=[]; rec=[]
    for s in range(nsteps):
        vel+=0.5*DT*acc; pos+=DT*vel
        acc,U=acc_and_U(pos); vel+=0.5*DT*acc
        if s%5==0: vbuf.append(vel.copy())
        if len(vbuf)>lag+1: vbuf.pop(0)
        if s>=burn and s%stride==0 and len(vbuf)>lag:
            h=H_all(pos); cen=centroids(pos)
            th=turn_angle(vbuf[0], vbuf[-1])
            z_hi,k_hi=domain_score(h,cen,rng=rng,top=True)
            z_lo,k_lo=domain_score(h,cen,rng=rng,top=False)
            rec.append((s*DT, float(h.mean()), float(np.mean(h>0.5)),
                        float(th.max()), float(th.mean()),
                        float(np.sum(th>2.0)),          # hard reversals
                        z_hi, z_lo, U, n_evaporated(pos),
                        atom_concentration(h),
                        domain_score_disjoint(h,cen,rng=rng)))
    return np.array(rec, dtype=float)

if __name__=="__main__":
    out={}; t0=time.time()
    for T0 in [0.40,0.60,0.75,0.90,1.10,1.40]:
        runs=[run(T0,sd) for sd in range(4)]
        out[T0]=runs
        A=np.vstack(runs); A=A[~np.isnan(A).any(axis=1)]
        print(f"T0={T0:.2f}  H_frac={A[:,2].mean():.3f}  "
              f"turn_max={A[:,3].mean():.2f}  turn_mean={A[:,4].mean():.3f}  "
              f"hard_turns={A[:,5].mean():.2f}  "
              f"z_dis={A[:,6].mean():+.2f}  z_ord={A[:,7].mean():+.2f}  "
              f"atom_conc={A[:,10].mean():.3f}  z_disjoint={A[:,11].mean():+.2f}  "
              f"[{time.time()-t0:.0f}s]", flush=True)
    pickle.dump(out, open('/home/claude/md6/domains.pkl','wb'))
    print("saved")
