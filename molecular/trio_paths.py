"""Per-TRIPLE centroid path geometry.

Does a trio whose centroid traces a smooth, coherent arc (rather than
jittering) indicate instability?  And do same-H-band trios move coherently
together ("co-orbiting")?
"""
import numpy as np, pickle, time, sys
sys.path.insert(0,'/home/claude/md3')
from cluster import (N, DT, icosahedron13, acc_and_U, H_all, pair_d,
                     TRI, PRS, lindemann, n_evaporated)
NT=len(TRI)
EDGES=np.array([0.0,0.05,0.10,0.20,0.35,0.55,1.01]); NB=len(EDGES)-1
LAG=3                       # short: velocity memory persists to ~lag 10

def cents(pos): return pos[TRI].mean(axis=1)

def run(T0, seed, nsteps=120000, rec=25, burn=10000):
    rng=np.random.default_rng(seed)
    pos=icosahedron13()+rng.normal(0,0.01,(N,3))
    vel=rng.normal(0,np.sqrt(T0),(N,3)); vel-=vel.mean(0)
    acc,_=acc_and_U(pos)
    r1=np.zeros(len(PRS)); r2=np.zeros(len(PRS)); ns=0
    cbuf=[]; rows=[]
    for s in range(nsteps):
        vel+=0.5*DT*acc; pos+=DT*vel
        acc,U=acc_and_U(pos); vel+=0.5*DT*acc
        if s%10==0:
            d=pair_d(pos); r1+=d; r2+=d**2; ns+=1
        if s%5==0:
            cbuf.append(cents(pos))
            if len(cbuf)>2*LAG+1: cbuf.pop(0)
        if s>=burn and s%rec==0 and len(cbuf)>2*LAG:
            h=H_all(pos); b=np.clip(np.digitize(h,EDGES)-1,0,NB-1)
            c0,c1,c2=cbuf[-1-2*LAG],cbuf[-1-LAG],cbuf[-1]
            v1=c1-c0; v2=c2-c1                          # centroid velocities
            n1=np.linalg.norm(v1,axis=1)+1e-12; n2=np.linalg.norm(v2,axis=1)+1e-12
            cosang=np.einsum('ij,ij->i',v1,v2)/(n1*n2)
            turn=np.arccos(np.clip(cosang,-1,1))        # path turning angle
            speed=0.5*(n1+n2)
            # straightness: net displacement / path length
            straight=np.linalg.norm(c2-c0,axis=1)/(n1+n2)
            row=[s*DT]
            for bi in range(NB):                        # per-band means
                m=b==bi
                row += ([float(turn[m].mean()), float(straight[m].mean()),
                         float(speed[m].mean()), float(m.mean())]
                        if m.sum()>=3 else [np.nan]*4)
            # coherence: do same-band trios share a direction?
            coh=[]; coh_dj=[]
            for bi in range(NB):
                m=np.where(b==bi)[0]
                if len(m)<4: coh.append(np.nan); coh_dj.append(np.nan); continue
                u=v2[m]/ (np.linalg.norm(v2[m],axis=1)[:,None]+1e-12)
                coh.append(float(np.linalg.norm(u.mean(0))))   # 0=random,1=aligned
                # DISJOINT control: greedily keep trios sharing no atoms
                used=set(); keep=[]
                for ti in m[np.argsort(-speed[m])]:
                    a=set(TRI[ti])
                    if used & a: continue
                    used|=a; keep.append(ti)
                if len(keep)<3: coh_dj.append(np.nan); continue
                ud=v2[keep]/(np.linalg.norm(v2[keep],axis=1)[:,None]+1e-12)
                coh_dj.append(float(np.linalg.norm(ud.mean(0))))
            row += coh + coh_dj
            row += [lindemann(r2,r1,ns,pos), n_evaporated(pos), float(h.mean())]
            rows.append(row)
    return np.array(rows, dtype=float)

if __name__=="__main__":
    out={}; t0=time.time()
    for T0 in [0.60,0.75,0.90,1.10]:
        out[T0]=[run(T0,sd) for sd in range(20)]
        A=np.vstack(out[T0])
        print(f"T0={T0:.2f}", flush=True)
        print(f"   {'band':<8s}{'turn':>8s}{'straight':>10s}{'speed':>9s}{'pop':>8s}"
              f"{'coher':>9s}{'coher_disj':>12s}")
        for bi in range(NB):
            o=1+4*bi
            if np.isnan(A[:,o]).all(): continue
            print(f"   {EDGES[bi]:<8.2f}{np.nanmean(A[:,o]):8.3f}{np.nanmean(A[:,o+1]):10.3f}"
                  f"{np.nanmean(A[:,o+2]):9.4f}{np.nanmean(A[:,o+3]):8.3f}"
                  f"{np.nanmean(A[:,1+4*NB+bi]):9.3f}{np.nanmean(A[:,1+5*NB+bi]):12.3f}")
        print(f"   lind={np.nanmean(A[:,-3]):.3f}  evap={np.nanmean(A[:,-2]):.2f} "
              f"[{time.time()-t0:.0f}s]", flush=True)
    pickle.dump(out, open('/home/claude/md9/trio.pkl','wb'))
    print("saved")
