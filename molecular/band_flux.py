"""Band POPULATIONS and band-to-band TRANSFER RATES.

Composition is a state variable. Flux between bands is a rate, and a net
upward flux (toward higher H) means the system is actively disordering -
which is a precursor by construction, not by correlation.

30 seeds per temperature for power.
"""
import numpy as np, pickle, time, sys
sys.path.insert(0,'/home/claude/md3')
from cluster import (N, DT, icosahedron13, acc_and_U, H_all, pair_d,
                     PRS, lindemann, n_evaporated)

EDGES=np.array([0.0,0.05,0.10,0.20,0.35,0.55,1.01]); NB=len(EDGES)-1
FLUX_STRIDE=5      # short, so triples move at most one band per step

def band_of(h):
    return np.clip(np.digitize(h,EDGES)-1, 0, NB-1)

def run(T0, seed, nsteps=150000, rec_stride=25, burn=10000):
    rng=np.random.default_rng(seed)
    pos=icosahedron13()+rng.normal(0,0.01,(N,3))
    vel=rng.normal(0,np.sqrt(T0),(N,3)); vel-=vel.mean(0)
    acc,_=acc_and_U(pos)
    r1=np.zeros(len(PRS)); r2=np.zeros(len(PRS)); ns=0
    prev_b=None; Tmat=np.zeros((NB,NB)); nflux=0
    rows=[]
    for s in range(nsteps):
        vel+=0.5*DT*acc; pos+=DT*vel
        acc,U=acc_and_U(pos); vel+=0.5*DT*acc
        if s%10==0:
            d=pair_d(pos); r1+=d; r2+=d**2; ns+=1
        if s>=burn and s%FLUX_STRIDE==0:
            b=band_of(H_all(pos))
            if prev_b is not None:
                np.add.at(Tmat,(prev_b,b),1); nflux+=1
            prev_b=b
        if s>=burn and s%rec_stride==0:
            h=H_all(pos); b=band_of(h)
            pop=np.bincount(b,minlength=NB)/len(h)
            # instantaneous flux over the last window
            if nflux>0:
                M=Tmat/Tmat.sum()
                up=sum(M[i,j] for i in range(NB) for j in range(NB) if j>i)
                dn=sum(M[i,j] for i in range(NB) for j in range(NB) if j<i)
                mob=1-np.trace(M)                 # fraction changing band
            else:
                up=dn=mob=0.0
            rows.append(np.concatenate([[s*DT],pop,[up,dn,up-dn,mob,
                        lindemann(r2,r1,ns,pos), n_evaporated(pos), U]]))
            Tmat[:]=0; nflux=0                    # reset window
    return np.array(rows)

if __name__=="__main__":
    out={}; t0=time.time()
    for T0 in [0.60,0.75,0.90,1.10]:
        out[T0]=[run(T0,sd) for sd in range(30)]
        A=np.vstack(out[T0])
        nev=sum(1 for a in out[T0] if a[:,-2].max()>0)
        print(f"T0={T0:.2f}  up={A[:,1+NB].mean():.4f} dn={A[:,2+NB].mean():.4f} "
              f"net={A[:,3+NB].mean():+.5f} mobility={A[:,4+NB].mean():.4f}  "
              f"lind={A[:,5+NB].mean():.3f}  {nev}/30 evaporated "
              f"[{time.time()-t0:.0f}s]", flush=True)
    pickle.dump(out, open('/home/claude/md8/flux.pkl','wb'))
    print("saved")
