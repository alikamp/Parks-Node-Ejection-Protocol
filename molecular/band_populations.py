"""Population ratios of H bands over time.
Does the COMPOSITION shift before evaporation, and before Lindemann moves?"""
import numpy as np, pickle, time, sys
sys.path.insert(0,'/home/claude/md3')
from cluster import (N, DT, icosahedron13, acc_and_U, H_all, pair_d,
                     PRS, lindemann, n_evaporated)

# bands chosen from the solid-state H distribution (solid sits ~0.09, p99 ~0.21)
EDGES = np.array([0.0, 0.05, 0.10, 0.20, 0.35, 0.55, 1.01])
NB = len(EDGES)-1

def bands(h):
    b = np.digitize(h, EDGES)-1
    return np.bincount(np.clip(b,0,NB-1), minlength=NB)/len(h)

def run(T0, seed, nsteps=150000, stride=25, burn=10000):
    rng=np.random.default_rng(seed)
    pos=icosahedron13()+rng.normal(0,0.01,(N,3))
    vel=rng.normal(0,np.sqrt(T0),(N,3)); vel-=vel.mean(0)
    acc,_=acc_and_U(pos)
    r1=np.zeros(len(PRS)); r2=np.zeros(len(PRS)); ns=0
    rows=[]
    for s in range(nsteps):
        vel+=0.5*DT*acc; pos+=DT*vel
        acc,U=acc_and_U(pos); vel+=0.5*DT*acc
        if s%10==0:
            d=pair_d(pos); r1+=d; r2+=d**2; ns+=1
        if s>=burn and s%stride==0:
            h=H_all(pos)
            rows.append(np.concatenate([[s*DT], bands(h),
                        [lindemann(r2,r1,ns,pos), n_evaporated(pos), U]]))
    return np.array(rows)

if __name__=="__main__":
    out={}; t0=time.time()
    for T0 in [0.40,0.60,0.75,0.90,1.10]:
        out[T0]=[run(T0,sd) for sd in range(8)]
        A=np.vstack(out[T0])
        print(f"T0={T0:.2f}  bands=" +
              " ".join(f"{A[:,1+i].mean():.3f}" for i in range(NB)) +
              f"  lind={A[:,1+NB].mean():.3f} evap={A[:,2+NB].mean():.2f}"
              f"  [{time.time()-t0:.0f}s]", flush=True)
    pickle.dump(out, open('/home/claude/md7/pop.pkl','wb'))
    print("saved")
