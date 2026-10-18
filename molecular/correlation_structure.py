"""Spatial and temporal correlation structure of the H field over triples.

For every pair of the 286 triples in an LJ-13 cluster:
  - overlap        : shared atoms (0, 1, or 2)
  - centroid dist  : spatial separation of the two triples
  - equal-time corr: <H_i H_j> - <H_i><H_j>  over the trajectory
  - co-movement    : Pearson r between H_i(t) and H_j(t)

Question: does the H field have a correlation length, and does it grow
as the cluster approaches melting / evaporation?
"""
import numpy as np, pickle, time, sys, itertools
sys.path.insert(0,'/home/claude/md3')
from cluster import (N, DT, icosahedron13, acc_and_U, H_all, pair_d,
                     TRI, TRI_P, PRS, lindemann, n_evaporated)

NT = len(TRI)
# pair bookkeeping over triples
TP = np.array(list(itertools.combinations(range(NT), 2)))
OVL = np.array([len(set(TRI[a]) & set(TRI[b])) for a, b in TP], dtype=np.int8)

def centroids(pos):
    return pos[TRI].mean(axis=1)          # (286,3)

def run(T0, seed, nsteps=120000, stride=40, burn=20000):
    rng=np.random.default_rng(seed)
    pos=icosahedron13()+rng.normal(0,0.01,(N,3))
    vel=rng.normal(0,np.sqrt(T0),(N,3)); vel-=vel.mean(0)
    acc,_=acc_and_U(pos)
    Hs=[]; Ds=[]; r1=np.zeros(len(PRS)); r2=np.zeros(len(PRS)); ns=0
    for s in range(nsteps):
        vel+=0.5*DT*acc; pos+=DT*vel
        acc,U=acc_and_U(pos); vel+=0.5*DT*acc
        if s>=burn and s%10==0:
            d=pair_d(pos); r1+=d; r2+=d**2; ns+=1
        if s>=burn and s%stride==0:
            Hs.append(H_all(pos))
            c=centroids(pos)
            Ds.append(np.linalg.norm(c[TP[:,0]]-c[TP[:,1]],axis=1))
    H=np.array(Hs); D=np.array(Ds)                   # (frames,286), (frames,40755)
    return H, D, lindemann(r2,r1,ns,pos), n_evaporated(pos)

def analyse(H, D, nbin=12):
    """correlation vs separation, split by overlap class"""
    # co-movement: Pearson r between H_i(t) and H_j(t)
    Hz=(H-H.mean(0))/(H.std(0)+1e-12)
    r_ij=(Hz[:,TP[:,0]]*Hz[:,TP[:,1]]).mean(0)       # (40755,)
    dmean=D.mean(0)
    out={}
    for ov in (0,1,2):
        m=OVL==ov
        if m.sum()<50: continue
        dd, rr = dmean[m], r_ij[m]
        edges=np.quantile(dd, np.linspace(0,1,nbin+1))
        edges=np.unique(edges)
        b=np.clip(np.digitize(dd,edges)-1,0,len(edges)-2)
        prof=[(float(dd[b==k].mean()), float(rr[b==k].mean()), int((b==k).sum()))
              for k in range(len(edges)-1) if (b==k).sum()>20]
        out[ov]=prof
    return out, r_ij, dmean

def corr_length(prof):
    """fit C(r) ~ A exp(-r/xi); return xi"""
    d=np.array([p[0] for p in prof]); c=np.array([p[1] for p in prof])
    m=c>0.01
    if m.sum()<4: return np.nan
    sl=np.polyfit(d[m], np.log(c[m]), 1)[0]
    return float(-1/sl) if sl<0 else np.nan

if __name__=="__main__":
    res=[]; t0=time.time()
    for T0 in [0.20,0.40,0.60,0.75,0.90,1.10,1.40]:
        acc_prof={0:[],1:[],2:[]}; lin=[]; ev=[]
        for sd in range(3):
            H,D,l,e=run(T0,sd)
            p,_,_=analyse(H,D); lin.append(l); ev.append(e)
            for ov in p: acc_prof[ov].append(p[ov])
        row={'T0':T0,'lind':float(np.mean(lin)),'evap':float(np.mean(ev))}
        for ov in (0,1,2):
            if not acc_prof[ov]: continue
            n=min(len(x) for x in acc_prof[ov])
            d=np.mean([[x[k][0] for k in range(n)] for x in acc_prof[ov]],axis=0)
            c=np.mean([[x[k][1] for k in range(n)] for x in acc_prof[ov]],axis=0)
            row[f'xi{ov}']=corr_length(list(zip(d,c)))
            row[f'c0_{ov}']=float(c[0]); row[f'cN_{ov}']=float(c[-1])
        res.append(row)
        print(f"T0={T0:.2f} lind={row['lind']:.3f} evap={row['evap']:.1f} "
              f"xi(ovl0)={row.get('xi0',float('nan')):.2f} "
              f"xi(ovl1)={row.get('xi1',float('nan')):.2f} "
              f"c_near(ovl0)={row.get('c0_0',float('nan')):+.3f} "
              f"c_far(ovl0)={row.get('cN_0',float('nan')):+.3f} "
              f"[{time.time()-t0:.0f}s]", flush=True)
    pickle.dump(res, open('/home/claude/md5/corr.pkl','wb'))
    print("saved")
