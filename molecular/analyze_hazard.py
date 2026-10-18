"""(a) Does H beat the obvious baselines?  (b) Hazard calibration."""
import pickle, numpy as np
D = pickle.load(open('/home/claude/md2/sweep.pkl','rb'))
diss = [r for r in D if r['td'] is not None and r['td'] > 5]
print(f"{len(D)} trajectories, {len(diss)} dissociations analysed\n")

def sustained(t, sig, thr):
    ok = sig > thr
    for i in range(len(ok)):
        if ok[i] and ok[i:].all(): return t[i]
    return None

print("=== (a) BASELINE COMPARISON ===")
print("same sustained-crossing rule; threshold tuned per signal\n")
print(f"{'signal':<12s}{'thr':>6s}{'n':>5s}{'median':>9s}{'min':>8s}{'CV':>8s}{'cover':>7s}")
best = {}
for name, key, grid in (("H",'H',np.arange(0.20,0.85,0.05)),
                        ("max dist",'mx',np.arange(2.0,9.0,0.5)),
                        ("mean dist",'mn',np.arange(1.5,6.0,0.5))):
    rows=[]
    for thr in grid:
        leads=[]
        for r in diss:
            tw = sustained(r['t'], r[key], thr)
            if tw is not None and tw < r['td']: leads.append(r['td']-tw)
        if len(leads) >= 0.5*len(diss):
            L=np.array(leads)
            rows.append((L.std()/L.mean(), thr, L))
    if not rows: 
        print(f"{name:<12s}  no threshold covers >=50%"); continue
    cv, thr, L = min(rows)
    best[name]=(thr,L)
    print(f"{name:<12s}{thr:6.2f}{len(L):5d}{np.median(L):9.2f}"
          f"{L.min():8.2f}{cv:8.3f}{len(L)/len(diss):7.0%}")

print("\n=== (b) HAZARD CALIBRATION ===")
print("P(dissociate within delta | current H), split-half validated\n")
rng=np.random.default_rng(0); idx=rng.permutation(len(diss))
train=[diss[i] for i in idx[:len(idx)//2]]; test=[diss[i] for i in idx[len(idx)//2:]]

def pairs(S, delta):
    """(H_now, dissociated_within_delta) over all pre-event frames"""
    X,Y=[],[]
    for r in S:
        m = r['t'] < r['td']
        X.append(r['H'][m]); Y.append((r['td']-r['t'][m]) <= delta)
    return np.concatenate(X), np.concatenate(Y)

EDGES=np.array([0,.1,.2,.3,.4,.5,.6,.7,.8,1.01])
for delta in (5.0, 10.0, 20.0):
    Xtr,Ytr = pairs(train,delta); Xte,Yte = pairs(test,delta)
    btr=np.digitize(Xtr,EDGES)-1; bte=np.digitize(Xte,EDGES)-1
    print(f"delta = {delta:.0f} time units")
    print(f"  {'H bin':<12s}{'train P':>9s}{'test P':>9s}{'n_test':>8s}")
    rows=[]
    for b in range(len(EDGES)-1):
        mtr=btr==b; mte=bte==b
        if mtr.sum()<30 or mte.sum()<30: continue
        ptr=Ytr[mtr].mean(); pte=Yte[mte].mean()
        rows.append((ptr,pte,mte.sum()))
        print(f"  {EDGES[b]:.1f}-{EDGES[b+1]:<8.1f}{ptr:9.1%}{pte:9.1%}{mte.sum():8d}")
    if len(rows)>=3:
        a=np.array(rows)
        err=np.abs(a[:,0]-a[:,1])
        print(f"  calibration error: mean {err.mean():.1%}  max {err.max():.1%}")
        print(f"  hazard ratio top/bottom bin: {a[-1,1]/max(a[0,1],1e-9):.1f}x")
    print()
