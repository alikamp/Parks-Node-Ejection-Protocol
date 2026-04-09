"""
Large batch test: T_env/T_outer threshold as Type A vs B predictor.
Also vary mass ratios and eccentricities to test universality.
"""
import numpy as np
from scipy.fft import fft
from scipy.signal import find_peaks

def leapfrog_nodes(pos, vel, masses, dt=0.008, max_t=3000, eject_r=70.0):
    G=1.0; n=len(masses); pos=pos.copy(); vel=vel.copy()
    node_times=[]
    def pdist(p):
        return {pr:np.linalg.norm(p[pr[0]]-p[pr[1]])
                for pr in [(0,1),(1,2),(0,2)]}
    d2=pdist(pos); pos+=vel*(dt/2); d1=pdist(pos); t=0.0
    for _ in range(int(max_t/dt)):
        acc=np.zeros_like(pos)
        for i in range(n):
            for j in range(n):
                if i!=j:
                    r=pos[j]-pos[i]; d=np.linalg.norm(r)+1e-10
                    acc[i]+=G*masses[j]*r/d**3
        vel+=acc*dt; pos+=vel*dt; t+=dt; dc=pdist(pos)
        for pr in [(0,1),(1,2),(0,2)]:
            if d1[pr]<d2[pr] and d1[pr]<dc[pr]:
                node_times.append(t)
        d2=d1; d1=dc
        com=np.average(pos,weights=masses,axis=0)
        for i in range(n):
            if np.linalg.norm(pos[i]-com)>eject_r:
                return np.array(node_times), t
    return np.array(node_times), None

def make_ic(a_out, masses, e_out=0.55, a_in=1.0, e_in=0.2):
    pos=np.zeros((3,3)); vel=np.zeros((3,3))
    r_in=a_in*(1-e_in)
    v_in=np.sqrt((masses[0]+masses[1])*(1+e_in)/(a_in*(1-e_in)))
    pos[0]=[-masses[1]/(masses[0]+masses[1])*r_in,0,0]
    pos[1]=[ masses[0]/(masses[0]+masses[1])*r_in,0,0]
    vel[0]=[0,-masses[1]/(masses[0]+masses[1])*v_in,0]
    vel[1]=[0, masses[0]/(masses[0]+masses[1])*v_in,0]
    r_out=a_out*(1+e_out); M=masses.sum()
    v_out=np.sqrt(M*(1-e_out)/(a_out*(1+e_out)))
    pos[2]=[r_out,0,0]; vel[2]=[0,v_out,0]
    cp=np.average(pos,weights=masses,axis=0)
    cv=np.average(vel,weights=masses,axis=0)
    return pos-cp, vel-cv

def get_T_env(nt):
    if len(nt)<80: return None
    dt_n=np.diff(nt); mean_dt=dt_n.mean()
    win=40
    if len(dt_n)<win+20: return None
    env=np.array([dt_n[i:i+win].std() for i in range(len(dt_n)-win)])
    ed=env-np.polyval(np.polyfit(np.arange(len(env)),env,1),
                       np.arange(len(env)))
    N=len(ed); F=np.abs(fft(ed))[:N//2]; fr=np.arange(N//2)/N
    pks,_=find_peaks(F, height=np.percentile(F,80))
    if not len(pks): return None
    f=fr[pks[np.argmax(F[pks])]]
    return (1/f)*mean_dt if f>0 else None

def kepler_T(a, M, G=1.0):
    return 2*np.pi*np.sqrt(a**3/(G*M))

# ── varied parameter sets ────────────────────────────────────────────────────
# Base + variations in mass ratio and eccentricity
param_sets = [
    # (label, masses, a_out_range, e_out, e_in)
    ('base',    np.array([1.0,1.0,0.3]), np.arange(3.0,5.2,0.2), 0.55, 0.20),
    ('heavy_3', np.array([1.0,1.0,0.5]), np.arange(3.0,5.2,0.2), 0.55, 0.20),
    ('light_3', np.array([1.0,1.0,0.1]), np.arange(3.0,5.2,0.2), 0.55, 0.20),
    ('hi_eout',  np.array([1.0,1.0,0.3]), np.arange(3.0,5.2,0.2), 0.70, 0.20),
    ('lo_eout',  np.array([1.0,1.0,0.3]), np.arange(3.5,6.0,0.3), 0.40, 0.20),
]

all_results = []

for label, masses, a_out_vals, e_out, e_in in param_sets:
    M = masses.sum()
    for a_out in a_out_vals:
        pos, vel = make_ic(a_out, masses, e_out=e_out, e_in=e_in)
        nt, t_ej = leapfrog_nodes(pos, vel, masses)
        if not t_ej: continue
        T_env = get_T_env(nt)
        if not T_env: continue
        k = t_ej / T_env
        if k > 6: continue  # exclude extreme outliers
        typ = 'A' if k <= 1.25 else 'B'
        T_outer = kepler_T(a_out, M)
        ratio = T_env / T_outer
        all_results.append(dict(
            label=label, a_out=a_out, masses=masses.tolist(),
            e_out=e_out, k=k, type=typ,
            T_env=T_env, T_outer=T_outer, ratio=ratio
        ))

print(f"Total valid systems: {len(all_results)}")
typeA = [r for r in all_results if r['type']=='A']
typeB = [r for r in all_results if r['type']=='B']
print(f"Type A: {len(typeA)}  Type B: {len(typeB)}")

ratios_A = np.array([r['ratio'] for r in typeA])
ratios_B = np.array([r['ratio'] for r in typeB])

print(f"\nT_env/T_outer statistics:")
print(f"  Type A: mean={np.mean(ratios_A):.2f}  "
      f"median={np.median(ratios_A):.2f}  "
      f"std={np.std(ratios_A):.2f}  "
      f"min={np.min(ratios_A):.2f}  max={np.max(ratios_A):.2f}")
print(f"  Type B: mean={np.mean(ratios_B):.2f}  "
      f"median={np.median(ratios_B):.2f}  "
      f"std={np.std(ratios_B):.2f}  "
      f"min={np.min(ratios_B):.2f}  max={np.max(ratios_B):.2f}")

# find optimal threshold
ks_all   = np.array([r['k']    for r in all_results])
rats_all = np.array([r['ratio'] for r in all_results])
corr = np.corrcoef(ks_all, rats_all)[0,1]
print(f"\nCorrelation k vs T_env/T_outer (all systems): r = {corr:+.4f}")

# grid search for best threshold
best_acc = 0; best_thresh = None
for thresh in np.arange(5, 50, 0.5):
    predicted = ['A' if r['ratio'] > thresh else 'B' for r in all_results]
    actual    = [r['type'] for r in all_results]
    acc = sum(p==a for p,a in zip(predicted,actual)) / len(all_results)
    if acc > best_acc:
        best_acc = acc; best_thresh = thresh

print(f"\nBest threshold T_env/T_outer > {best_thresh:.1f} → Type A")
print(f"Accuracy: {best_acc:.1%} ({int(best_acc*len(all_results))}/{len(all_results)})")

# confusion matrix at best threshold
predicted = ['A' if r['ratio']>best_thresh else 'B' for r in all_results]
actual    = [r['type'] for r in all_results]
TP = sum(p=='A' and a=='A' for p,a in zip(predicted,actual))
TN = sum(p=='B' and a=='B' for p,a in zip(predicted,actual))
FP = sum(p=='A' and a=='B' for p,a in zip(predicted,actual))
FN = sum(p=='B' and a=='A' for p,a in zip(predicted,actual))
print(f"TP={TP} TN={TN} FP={FP} FN={FN}")
spec = TN/(TN+FP) if (TN+FP)>0 else 0
sens = TP/(TP+FN) if (TP+FN)>0 else 0
print(f"Specificity: {spec:.1%}  Sensitivity: {sens:.1%}")

# per parameter set breakdown
print(f"\nPer parameter set:")
for label in [p[0] for p in param_sets]:
    sub = [r for r in all_results if r['label']==label]
    if not sub: continue
    pred = ['A' if r['ratio']>best_thresh else 'B' for r in sub]
    act  = [r['type'] for r in sub]
    acc  = sum(p==a for p,a in zip(pred,act))/len(sub)
    nA = sum(1 for r in sub if r['type']=='A')
    nB = sum(1 for r in sub if r['type']=='B')
    print(f"  {label:10s}: n={len(sub):3d} (A={nA},B={nB})  acc={acc:.1%}")

# show misclassified
print(f"\nMisclassified systems:")
for r,p in zip(all_results,predicted):
    if p!=r['type']:
        print(f"  {r['label']:10s} a_out={r['a_out']:.1f} "
              f"k={r['k']:.3f} actual={r['type']} pred={p} "
              f"ratio={r['ratio']:.2f}")
