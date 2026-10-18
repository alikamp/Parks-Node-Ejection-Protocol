"""Crystallisation: quench a hot LJ-13 cluster and watch H.
Does the H DISTRIBUTION detect nucleation before the mean moves?"""
import numpy as np, pickle, time, sys
sys.path.insert(0, '/home/claude/md3')
from cluster import (N, DT, icosahedron13, acc_and_U, H_all, pair_d,
                     lindemann, n_evaporated, PRS, TRI, TRI_P)

def bimodality(h):
    """Sarle's bimodality coefficient: >5/9 suggests bimodal"""
    n=len(h); m=h.mean(); s=h.std()
    if s<1e-12: return np.nan
    g1=np.mean((h-m)**3)/s**3
    g2=np.mean((h-m)**4)/s**4-3
    return (g1**2+1)/(g2+3*(n-1)**2/((n-2)*(n-3)))

def q6_global(pos, rcut=1.5):
    """Global bond-orientational order Q6 - the standard crystallinity
    measure. Low = disordered, high = ordered."""
    from scipy.special import sph_harm_y
    qsum=np.zeros(13, complex); nb=0
    for i in range(N):
        for j in range(N):
            if i==j: continue
            d=pos[j]-pos[i]; r=np.linalg.norm(d)
            if r>rcut: continue
            th=np.arccos(np.clip(d[2]/r,-1,1)); ph=np.arctan2(d[1],d[0])
            for mi,m in enumerate(range(-6,7)):
                qsum[mi]+=sph_harm_y(6,m,th,ph)
            nb+=1
    if nb==0: return 0.0
    qsum/=nb
    return float(np.sqrt(4*np.pi/13*np.sum(np.abs(qsum)**2)))

def quench(T_hot, T_cold, seed, n_hot=40000, n_cold=200000, stride=50, gamma=2.0):
    """melt at T_hot, then Langevin-thermostat down to T_cold"""
    rng=np.random.default_rng(seed)
    pos=icosahedron13()+rng.normal(0,0.01,(N,3))
    vel=rng.normal(0,np.sqrt(T_hot),(N,3)); vel-=vel.mean(0)
    acc,_=acc_and_U(pos)
    # --- melt (NVE) ---
    for s in range(n_hot):
        vel+=0.5*DT*acc; pos+=DT*vel
        acc,_=acc_and_U(pos); vel+=0.5*DT*acc
    # --- quench (Langevin) ---
    rec=[]
    c1=np.exp(-gamma*DT); c2=np.sqrt((1-c1**2)*T_cold)
    for s in range(n_cold):
        vel+=0.5*DT*acc; pos+=DT*vel
        acc,U=acc_and_U(pos); vel+=0.5*DT*acc
        vel=c1*vel+c2*rng.normal(0,1,(N,3))        # thermostat
        if s%stride==0:
            h=H_all(pos)
            rec.append((s*DT, float(h.mean()), float(h.std()),
                        float(bimodality(h)), float(np.mean(h>0.5)),
                        float(np.percentile(h,90)), q6_global(pos), U))
    return np.array(rec)

if __name__=="__main__":
    out=[]; t0=time.time()
    for sd in range(10):
        out.append(quench(0.75, 0.06, sd, n_hot=60000, n_cold=400000))
        print(f"  quench {sd} done [{time.time()-t0:.0f}s]", flush=True)
    pickle.dump(out, open('/home/claude/md4/freeze.pkl','wb'))
    print("saved")
