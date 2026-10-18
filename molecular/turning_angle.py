"""Turning angle at MULTIPLE lags - find the lag where velocity still
carries memory, then test whether sharp turns track anything."""
import numpy as np, sys, time
sys.path.insert(0,'/home/claude/md3')
from cluster import N, DT, icosahedron13, acc_and_U, H_all, n_evaporated

LAGS=[1,2,5,10,25,50]          # in units of 5 integration steps

def run(T0, seed, nsteps=80000, burn=15000, stride=40):
    rng=np.random.default_rng(seed)
    pos=icosahedron13()+rng.normal(0,0.01,(N,3))
    vel=rng.normal(0,np.sqrt(T0),(N,3)); vel-=vel.mean(0)
    acc,_=acc_and_U(pos)
    buf=[]; rows=[]
    for s in range(nsteps):
        vel+=0.5*DT*acc; pos+=DT*vel
        acc,U=acc_and_U(pos); vel+=0.5*DT*acc
        if s%5==0: buf.append(vel.copy())
        if len(buf)>max(LAGS)+1: buf.pop(0)
        if s>=burn and s%stride==0 and len(buf)>max(LAGS):
            h=H_all(pos); r=[s*DT, float(h.mean()), float(np.mean(h>0.5))]
            for L in LAGS:
                v0,v1=buf[-1-L],buf[-1]
                a=np.einsum('ij,ij->i',v0,v1)
                b=np.linalg.norm(v0,axis=1)*np.linalg.norm(v1,axis=1)+1e-12
                th=np.arccos(np.clip(a/b,-1,1))
                r += [float(th.mean()), float(th.max())]
            r.append(n_evaporated(pos)); rows.append(r)
    return np.array(rows)

print("mean turning angle by lag  (pi/2 = 1.571 means fully decorrelated)")
print(f"{'T0':>6s}" + "".join(f"{'lag '+str(L):>11s}" for L in LAGS))
store={}
for T0 in [0.40,0.75,1.10,1.40]:
    A=np.vstack([run(T0,sd) for sd in range(3)])
    store[T0]=A
    print(f"{T0:6.2f}" + "".join(f"{A[:,3+2*i].mean():11.3f}" for i in range(len(LAGS))))
print()
print(f"dt per lag unit = {5*DT} time units; lag 1 = {5*DT}, lag 50 = {50*5*DT}")
print()
print("does the sharp-turn RATE track H or evaporation? (correlation within T0)")
from scipy import stats
print(f"{'T0':>6s}{'lag':>6s}{'rho(turn_max,H)':>18s}{'rho(turn_mean,H)':>18s}")
for T0,A in store.items():
    for i,L in enumerate(LAGS):
        if L not in (1,5): continue
        r1=stats.spearmanr(A[:,4+2*i],A[:,1]).statistic
        r2=stats.spearmanr(A[:,3+2*i],A[:,1]).statistic
        print(f"{T0:6.2f}{L:6d}{r1:18.3f}{r2:18.3f}")
