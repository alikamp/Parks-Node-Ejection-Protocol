"""One sweep: store H, max-dist, mean-dist time series + dissociation time
for many trajectories. Used for (a) baseline comparison (b) hazard calibration."""
import numpy as np, pickle, sys
sys.path.insert(0,'/home/claude/parks-node-ejection-protocol/molecular')
from pnep_molecular import (forces, dists, MAS, RMIN_IJ, DT, squash, tau)

def traj(T0, seed, nsteps=300000, diss_r=8.0, stride=25):
    rng=np.random.default_rng(seed); r=RMIN_IJ.mean()
    pos=np.array([[0,0,0],[r,0,0],[r/2,r*np.sqrt(3)/2,0]],float)
    vel=rng.normal(0,np.sqrt(T0/MAS[:,None]),(3,3))
    vel-=(MAS[:,None]*vel).sum(0)/MAS.sum()
    F,_=forces(pos); acc=F/MAS[:,None]
    T,H,MX,MN=[],[],[],[]; td=None
    for s in range(nsteps):
        vel+=0.5*DT*acc; pos+=DT*vel
        F,_=forces(pos); acc=F/MAS[:,None]; vel+=0.5*DT*acc
        if s%stride==0:
            t=s*DT; d=dists(pos)
            com=(MAS[:,None]*pos).sum(0)/MAS.sum()
            if td is None and np.linalg.norm(pos-com,axis=1).max()>diss_r: td=t
            T.append(t); H.append(squash(float(np.var(d))))
            MX.append(float(d.max())); MN.append(float(d.mean()))
        if td is not None and s*DT>td+2: break
    return dict(T0=T0,seed=seed,t=np.array(T),H=np.array(H),
                mx=np.array(MX),mn=np.array(MN),td=td)

out=[]
TEMPS=[0.40,0.45,0.50,0.55]          # near-threshold band, tau 0.94-1.29
for T0 in TEMPS:
    for sd in range(60):
        out.append(traj(T0,sd))
    done=sum(1 for r in out if r['T0']==T0 and r['td'] is not None)
    print(f"T0={T0}  tau={tau(T0):.2f}  dissociated {done}/60", flush=True)
pickle.dump(out, open('/home/claude/md2/sweep.pkl','wb'))
print("saved", len(out), "trajectories")
