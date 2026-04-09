"""
V_outer = speed of outer body relative to inner binary CoM.
This is the physically meaningful quantity for energy transfer.
High V_outer at node = outer body moving fast past inner binary = efficient transfer
Low V_outer at node = outer body barely grazing = inefficient transfer

Also compute:
- V_inner = speed of inner binary relative to system CoM (how energetic the engine is)
- Transfer ratio TR = V_outer / V_inner (how much of engine energy reaches outer body)
"""
import numpy as np
from scipy.fft import fft
from scipy.signal import find_peaks

def leapfrog_nodes(pos, vel, masses, dt=0.008, max_t=3000, eject_r=70.0):
    G=1.0; n=len(masses); pos=pos.copy(); vel=vel.copy()
    node_times=[]; node_H=[]; node_Vout=[]; node_TR=[]

    def pdist(p):
        return {pr:np.linalg.norm(p[pr[0]]-p[pr[1]])
                for pr in [(0,1),(1,2),(0,2)]}

    def compute_signals(p, v, m):
        # H
        s2=np.var(list(pdist(p).values())); H=s2/(1+s2)

        # inner binary CoM position and velocity
        m_inner = m[0]+m[1]
        inner_com_pos = (m[0]*p[0]+m[1]*p[1])/m_inner
        inner_com_vel = (m[0]*v[0]+m[1]*v[1])/m_inner

        # outer body velocity relative to inner binary CoM
        v_outer_rel = v[2] - inner_com_vel
        V_outer = np.linalg.norm(v_outer_rel)

        # inner binary relative velocity (engine speed)
        v_inner_rel = v[0] - v[1]
        V_inner = np.linalg.norm(v_inner_rel)

        # transfer ratio
        TR = V_outer / V_inner if V_inner > 0 else 0

        return H, V_outer, V_inner, TR

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
                H, Vo, Vi, TR = compute_signals(pos, vel, masses)
                node_times.append(t)
                node_H.append(H)
                node_Vout.append(Vo)
                node_TR.append(TR)

        d2=d1; d1=dc
        com=np.average(pos,weights=masses,axis=0)
        for i in range(n):
            if np.linalg.norm(pos[i]-com)>eject_r:
                return (np.array(node_times), np.array(node_H),
                        np.array(node_Vout), np.array(node_TR), t)

    return (np.array(node_times), np.array(node_H),
            np.array(node_Vout), np.array(node_TR), None)

def make_ic(a_out, e_out=0.55, masses=np.array([1.0,1.0,0.3]),
            a_in=1.0, e_in=0.2):
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
    return pos-cp, vel-cv, masses

def get_T_env(nt):
    if len(nt)<80: return None
    dt_n=np.diff(nt); mean_dt=dt_n.mean()
    win=40
    if len(dt_n)<win+20: return None
    env=np.array([dt_n[i:i+win].std() for i in range(len(dt_n)-win)])
    ed=env-np.polyval(np.polyfit(np.arange(len(env)),env,1),np.arange(len(env)))
    N=len(ed); F=np.abs(fft(ed))[:N//2]; fr=np.arange(N//2)/N
    pks,_=find_peaks(F,height=np.percentile(F,80))
    if not len(pks): return None
    f=fr[pks[np.argmax(F[pks])]]
    return (1/f)*mean_dt if f>0 else None

# ── scan ──────────────────────────────────────────────────────────────────────
a_out_vals = [3.0, 3.2, 3.4, 3.5, 3.6, 3.8, 4.2, 4.5, 4.8, 5.0]

print(f"{'a_out':>5} | {'k':>6} | {'type':>4} | "
      f"{'TR_mean':>8} | {'TR_early':>9} | {'TR_late':>8} | "
      f"{'dTR':>7} | {'Vo_mean':>8} | {'TR_x_H':>8}")
print("-"*85)

results=[]
for a_out in a_out_vals:
    pos,vel,masses = make_ic(a_out)
    nt,nH,nVo,nTR,t_ej = leapfrog_nodes(pos,vel,masses)
    T_env = get_T_env(nt)
    k = t_ej/T_env if (t_ej and T_env) else None
    typ = 'A' if (k and k<=1.25) else ('B' if (k and k>1.25) else '?')

    w = min(20, len(nTR)//4)
    TR_mean  = np.mean(nTR)  if len(nTR)>0 else None
    TR_early = np.mean(nTR[:w]) if w>0 else None
    TR_late  = np.mean(nTR[-w:]) if w>0 else None
    dTR = (TR_late-TR_early) if (TR_late and TR_early) else None
    Vo_mean  = np.mean(nVo)  if len(nVo)>0 else None

    # combined: TR * H — high transfer efficiency AND hierarchy intact
    TR_x_H = np.mean(np.array(nTR)*np.array(nH)) if len(nTR)>0 else None

    results.append(dict(a_out=a_out, k=k, type=typ,
                        TR_mean=TR_mean, TR_early=TR_early, TR_late=TR_late,
                        dTR=dTR, Vo_mean=Vo_mean, TR_x_H=TR_x_H,
                        nTR=nTR, nH=nH, nVo=nVo, nt=nt))

    print(f"{a_out:>5.1f} | "
          f"{f'{k:.3f}' if k else '---':>6} | "
          f"{typ:>4} | "
          f"{f'{TR_mean:.4f}' if TR_mean else '---':>8} | "
          f"{f'{TR_early:.4f}' if TR_early else '---':>9} | "
          f"{f'{TR_late:.4f}' if TR_late else '---':>8} | "
          f"{f'{dTR:+.4f}' if dTR else '---':>7} | "
          f"{f'{Vo_mean:.4f}' if Vo_mean else '---':>8} | "
          f"{f'{TR_x_H:.4f}' if TR_x_H else '---':>8}")

# ── correlation analysis ──────────────────────────────────────────────────────
ejecting = [r for r in results if r['k'] and r['TR_mean']]
if ejecting:
    ks       = np.array([r['k']       for r in ejecting])
    TR_means = np.array([r['TR_mean'] for r in ejecting])
    TR_early = np.array([r['TR_early'] for r in ejecting])
    dTRs     = np.array([r['dTR']     for r in ejecting if r['dTR']])
    TR_x_Hs  = np.array([r['TR_x_H']  for r in ejecting])

    print(f"\nCorrelations with k:")
    print(f"  TR_mean:  r = {np.corrcoef(ks, TR_means)[0,1]:+.4f}")
    print(f"  TR_early: r = {np.corrcoef(ks, TR_early)[0,1]:+.4f}")
    print(f"  TR_x_H:   r = {np.corrcoef(ks, TR_x_Hs)[0,1]:+.4f}")
    if len(dTRs)==len(ks):
        print(f"  dTR:      r = {np.corrcoef(ks, dTRs)[0,1]:+.4f}")

    # Type A vs Type B separation on TR_early
    typeA = [r for r in ejecting if r['type']=='A']
    typeB = [r for r in ejecting if r['type']=='B']
    if typeA and typeB:
        print(f"\nType A TR_early: {np.mean([r['TR_early'] for r in typeA]):.4f} "
              f"+/- {np.std([r['TR_early'] for r in typeA]):.4f}")
        print(f"Type B TR_early: {np.mean([r['TR_early'] for r in typeB]):.4f} "
              f"+/- {np.std([r['TR_early'] for r in typeB]):.4f}")
        print(f"Separation: {abs(np.mean([r['TR_early'] for r in typeA]) - np.mean([r['TR_early'] for r in typeB])):.4f}")

        print(f"\nType A TR_x_H:  {np.mean([r['TR_x_H'] for r in typeA]):.4f}")
        print(f"Type B TR_x_H:  {np.mean([r['TR_x_H'] for r in typeB]):.4f}")
        print(f"Separation: {abs(np.mean([r['TR_x_H'] for r in typeA]) - np.mean([r['TR_x_H'] for r in typeB])):.4f}")

        # dTR trend
        print(f"\nEnergy transfer trend dTR (late - early):")
        for r in ejecting:
            print(f"  a_out={r['a_out']} k={r['k']:.3f} type={r['type']} "
                  f"TR_early={r['TR_early']:.4f} TR_late={r['TR_late']:.4f} "
                  f"dTR={r['dTR']:+.4f}")
