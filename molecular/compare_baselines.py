"""Does H give a BETTER hazard calibration than max-distance?"""
import pickle, numpy as np
D=pickle.load(open('sweep.pkl','rb'))
diss=[r for r in D if r['td'] is not None and r['td']>5]
rng=np.random.default_rng(0); idx=rng.permutation(len(diss))
train=[diss[i] for i in idx[:len(idx)//2]]; test=[diss[i] for i in idx[len(idx)//2:]]

def pairs(S, key, delta):
    X,Y=[],[]
    for r in S:
        m=r['t']<r['td']
        X.append(r[key][m]); Y.append((r['td']-r['t'][m])<=delta)
    return np.concatenate(X), np.concatenate(Y)

def evaluate(key, delta, nbin=10):
    Xtr,Ytr=pairs(train,key,delta); Xte,Yte=pairs(test,key,delta)
    # quantile bins from TRAIN so each signal gets equal-mass bins
    edges=np.unique(np.quantile(Xtr,np.linspace(0,1,nbin+1)))
    btr=np.clip(np.digitize(Xtr,edges)-1,0,len(edges)-2)
    bte=np.clip(np.digitize(Xte,edges)-1,0,len(edges)-2)
    ptr=np.array([Ytr[btr==b].mean() if (btr==b).sum()>30 else np.nan
                  for b in range(len(edges)-1)])
    pred=ptr[bte]                       # predicted hazard for each test frame
    ok=~np.isnan(pred)
    pred,act=pred[ok],Yte[ok]
    # Brier score (lower better) and calibration error over test bins
    brier=float(np.mean((pred-act)**2))
    # discrimination: AUC
    pos,neg=pred[act],pred[~act]
    if len(pos)>5000: pos=np.random.default_rng(1).choice(pos,5000,replace=False)
    if len(neg)>5000: neg=np.random.default_rng(2).choice(neg,5000,replace=False)
    auc=float(np.mean((pos[:,None]>neg[None,:]).astype(float)
                      +0.5*(pos[:,None]==neg[None,:])))
    # calibration error on held-out bins
    ce=[]
    for b in range(len(edges)-1):
        m=bte==b
        if m.sum()>30 and not np.isnan(ptr[b]): ce.append(abs(ptr[b]-Yte[m].mean()))
    return brier, auc, (np.mean(ce) if ce else np.nan), float(np.nanmax(ptr))

print(f"{'delta':>6s} {'signal':<11s}{'Brier':>9s}{'AUC':>8s}{'calib err':>11s}{'top bin P':>11s}")
for delta in (5.0,10.0,20.0):
    for name,key in (("H",'H'),("max dist",'mx'),("mean dist",'mn')):
        b,a,c,tp=evaluate(key,delta)
        print(f"{delta:6.0f} {name:<11s}{b:9.4f}{a:8.3f}{c:11.1%}{tp:11.1%}")
    print()
