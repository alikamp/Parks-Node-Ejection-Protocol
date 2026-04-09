import numpy as np
from scipy.fft import fft
from scipy.signal import find_peaks

def leapfrog_nodes(pos, vel, masses, dt=0.004, max_t=5000, eject_r=80.0):
    G = 1.0
    n = len(masses)
    pos = pos.copy(); vel = vel.copy()
    node_times = []
    def pdist(p):
        return {pr: np.linalg.norm(p[pr[0]]-p[pr[1]])
                for pr in [(0,1),(1,2),(0,2)]}
    d2 = pdist(pos); pos += vel*(dt/2); d1 = pdist(pos)
    t = 0.0
    for _ in range(int(max_t/dt)):
        acc = np.zeros_like(pos)
        for i in range(n):
            for j in range(n):
                if i!=j:
                    r=pos[j]-pos[i]; d=np.linalg.norm(r)+1e-10
                    acc[i]+=G*masses[j]*r/d**3
        vel+=acc*dt; pos+=vel*dt; t+=dt
        dc=pdist(pos)
        for pr in [(0,1),(1,2),(0,2)]:
            if d1[pr]<d2[pr] and d1[pr]<dc[pr]:
                node_times.append(t)
        d2=d1; d1=dc
        com=np.average(pos,weights=masses,axis=0)
        for i in range(n):
            if np.linalg.norm(pos[i]-com)>eject_r:
                return np.array(node_times), t
    return np.array(node_times), None

def make_ic(a_out, masses=np.array([1.0,1.0,0.1]),
            a_in=1.0, e_in=0.15, e_out=0.40):
    pos=np.zeros((3,3)); vel=np.zeros((3,3))
    r_in=a_in*(1-e_in)
    v_in=np.sqrt((masses[0]+masses[1])*(1+e_in)/(a_in*(1-e_in)))
    pos[0]=[-masses[1]/(masses[0]+masses[1])*r_in,0,0]
    pos[1]=[ masses[0]/(masses[0]+masses[1])*r_in,0,0]
    vel[0]=[0,-masses[1]/(masses[0]+masses[1])*v_in,0]
    vel[1]=[0, masses[0]/(masses[0]+masses[1])*v_in,0]
    r_out=a_out*(1+e_out)
    v_out=np.sqrt((masses[0]+masses[1]+masses[2])*(1-e_out)/(a_out*(1+e_out)))
    pos[2]=[r_out,0,0]; vel[2]=[0,v_out,0]
    cp=np.average(pos,weights=masses,axis=0)
    cv=np.average(vel,weights=masses,axis=0)
    return pos-cp, vel-cv, masses

def dominant_freq(series, percentile=90):
    s=series-np.polyval(np.polyfit(np.arange(len(series)),series,1),
                         np.arange(len(series)))
    N=len(s); F=np.abs(fft(s))[:N//2]; fr=np.arange(N//2)/N
    pks,_=find_peaks(F,height=np.percentile(F,percentile))
    if not len(pks): return None
    top=pks[np.argmax(F[pks])]
    return fr[top]

a_out_vals = [4.5, 4.8, 5.0, 5.2, 5.5, 5.8, 6.0, 6.5, 7.0]
print(f"{'a_out':>6} | {'t_eject':>9} | {'nodes':>6} | {'f_B':>8} | {'f_env':>9} | {'T_env(t)':>9} | {'ratio':>7}")
print("-"*72)

results = []
for a_out in a_out_vals:
    pos,vel,masses = make_ic(a_out)
    nt, t_ej = leapfrog_nodes(pos,vel,masses)
    
    f_B = None; f_env = None; T_env = None; ratio = None
    
    if len(nt) > 60:
        dt_n = np.diff(nt)
        f_B = dominant_freq(dt_n, 90)
        
        # envelope
        win=40
        if len(dt_n) > win+20:
            env=np.array([dt_n[i:i+win].std() for i in range(len(dt_n)-win)])
            f_env = dominant_freq(env, 80)
            if f_env and f_env > 0:
                mean_dt = dt_n.mean()
                T_env = (1/f_env)*mean_dt
                if t_ej:
                    ratio = t_ej/T_env

    te_str  = f"{t_ej:.1f}" if t_ej else "no eject"
    fb_str  = f"{f_B:.5f}" if f_B else "---"
    fe_str  = f"{f_env:.5f}" if f_env else "---"
    te2_str = f"{T_env:.2f}" if T_env else "---"
    ra_str  = f"{ratio:.3f}" if ratio else "---"
    
    print(f"{a_out:>6.1f} | {te_str:>9} | {len(nt):>6} | "
          f"{fb_str:>8} | {fe_str:>9} | {te2_str:>9} | {ra_str:>7}")
    results.append(dict(a_out=a_out,t_ej=t_ej,n=len(nt),
                        f_B=f_B,f_env=f_env,T_env=T_env,ratio=ratio))

ejecting = [r for r in results if r['t_ej'] and r['ratio']]
if ejecting:
    ratios = [r['ratio'] for r in ejecting]
    print(f"\nt_eject / T_env ratios: {[round(r,3) for r in ratios]}")
    print(f"mean={np.mean(ratios):.3f}  std={np.std(ratios):.3f}  "
          f"cv={np.std(ratios)/np.mean(ratios):.3f}")
    
    # check if ratio is roughly constant — that's the signal
    if np.std(ratios)/np.mean(ratios) < 0.3:
        print("\n*** SIGNAL FOUND: ratio t_eject/T_env is approximately constant ***")
        print(f"    t_eject ~ {np.mean(ratios):.2f} x T_env")
