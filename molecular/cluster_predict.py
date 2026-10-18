"""In the evaporating regime: does H predict FUTURE evaporation better
than Lindemann, or than mean pairwise distance?"""
import numpy as np, pickle, itertools, time
exec(open('cluster.py').read().split('if __name__')[0])

def run2(T0, seed, nsteps=200000, stride=50):
    rng=np.random.default_rng(seed)
    pos=icosahedron13()+rng.normal(0,0.01,(N,3))
    vel=rng.normal(0,np.sqrt(T0),(N,3)); vel-=vel.mean(0)
    acc,_=acc_and_U(pos)
    r1=np.zeros(len(PRS)); r2=np.zeros(len(PRS)); ns=0; rec=[]
    for s in range(nsteps):
        vel+=0.5*DT*acc; pos+=DT*vel
        acc,U=acc_and_U(pos); vel+=0.5*DT*acc
        if s%10==0:
            d=pair_d(pos); r1+=d; r2+=d**2; ns+=1
        if s%stride==0:
            h=H_all(pos); d=pair_d(pos)
            rec.append((s*DT, float(h.mean()), float(np.percentile(h,99)),
                        float(np.mean(h>0.5)), lindemann(r2,r1,ns,pos),
                        float(d.mean()), float(d.max()), n_evaporated(pos)))
    return np.array(rec)

print("generating trajectories in the evaporating regime...", flush=True)
data=[]; t0=time.time()
for T0 in [0.9,1.1,1.3,1.6,2.0]:
    for sd in range(10):
        data.append(run2(T0,sd))
    print(f"  T0={T0} done [{time.time()-t0:.0f}s]", flush=True)
pickle.dump(data, open('predict.pkl','wb'))

COLS={'H_mean':1,'H_p99':2,'H_frac':3,'lindemann':4,'mean_dist':5,'max_dist':6}
EV=7
rng=np.random.default_rng(0); idx=rng.permutation(len(data))
tr=[data[i] for i in idx[:len(idx)//2]]; te=[data[i] for i in idx[len(idx)//2:]]

def xy(S, col, horizon_frames):
    X,Y=[],[]
    for a in S:
        n=len(a)
        for i in range(n-horizon_frames):
            X.append(a[i,col]); Y.append(a[i+horizon_frames,EV]-a[i,EV])
    return np.array(X), np.array(Y)

for H_FR,label in ((40,"delta = 8 time units"),(100,"delta = 20 time units")):
    print(f"\n{label}   (predicting NEW evaporations in the window)")
    print(f"  {'signal':<12s}{'spearman':>10s}{'R2 (test)':>11s}{'binned MAE':>12s}")
    for name,c in COLS.items():
        Xtr,Ytr=xy(tr,c,H_FR); Xte,Yte=xy(te,c,H_FR)
        ok=~np.isnan(Xtr); Xtr,Ytr=Xtr[ok],Ytr[ok]
        ok=~np.isnan(Xte); Xte,Yte=Xte[ok],Yte[ok]
        if len(Xtr)<500: continue
        from scipy import stats
        rho=stats.spearmanr(Xte,Yte).statistic
        # non-parametric: quantile-bin on train, predict bin mean
        edges=np.unique(np.quantile(Xtr,np.linspace(0,1,13)))
        btr=np.clip(np.digitize(Xtr,edges)-1,0,len(edges)-2)
        bte=np.clip(np.digitize(Xte,edges)-1,0,len(edges)-2)
        mu=np.array([Ytr[btr==b].mean() if (btr==b).sum()>20 else np.nan
                     for b in range(len(edges)-1)])
        pred=mu[bte]; m=~np.isnan(pred)
        ss_res=np.sum((Yte[m]-pred[m])**2); ss_tot=np.sum((Yte[m]-Yte[m].mean())**2)
        r2=1-ss_res/ss_tot
        mae=np.mean(np.abs(Yte[m]-pred[m]))
        print(f"  {name:<12s}{rho:10.3f}{r2:11.3f}{mae:12.3f}")
