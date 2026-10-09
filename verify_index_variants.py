"""Compare mass-blind H against three mass-aware hierarchy indices,
computed at the SAME mirror symmetry nodes in a single integration."""
import sys, json, time
sys.path.insert(0, '/home/claude/pnep_repo')
import numpy as np
import pnep_v12 as P

G, DT, EJECT_R, MIN_NODES = 1.0, 0.006, 12.0, 10
PAIRS = [(0,1),(1,2),(0,2)]

def squash(s2): return s2/(1.0+s2)

def indices(pos, m):
    d = np.array([np.linalg.norm(pos[j]-pos[i]) for i,j in PAIRS])
    # H  : mass-blind geometric variance (current)
    H  = squash(float(np.var(d)))
    # Hmu: reduced-mass-weighted separations
    mu = np.array([m[i]*m[j]/(m[i]+m[j]) for i,j in PAIRS])
    x  = mu*d;  Hmu = squash(float(np.var(x/ (x.mean()+1e-12))))
    # HE : pairwise binding energies  -G m_i m_j / d_ij
    E  = np.array([-G*m[i]*m[j]/(d[k]+1e-12) for k,(i,j) in enumerate(PAIRS)])
    HE = squash(float(np.var(E/(abs(E).mean()+1e-12))))
    # Hr : each body's distance to CoM of the other two, mass-weighted
    rr = []
    for k in range(3):
        a,b = [q for q in range(3) if q != k]
        cm  = (m[a]*pos[a]+m[b]*pos[b])/(m[a]+m[b])
        rr.append(m[k]*np.linalg.norm(pos[k]-cm))
    rr = np.array(rr); Hr = squash(float(np.var(rr/(rr.mean()+1e-12))))
    return H, Hmu, HE, Hr

def trial(want_stable, rng, gt_max_t):
    ic = P.make_ic(want_stable, rng)
    if ic is None: return None
    pos, vel, m, ma, T_out = ic
    t=0.0; nc=0; d1=d2=None; eject=None
    acc = [[],[],[],[]]          # first-MIN_NODES means, per index
    while t < gt_max_t:
        pos, vel = P.leapfrog(pos, vel, m, DT); t += DT
        if eject is None and t > T_out and np.any(np.linalg.norm(pos,axis=1) > EJECT_R):
            eject = t
        d = np.array([np.linalg.norm(pos[j]-pos[i]) for i,j in PAIRS])
        if P.is_local_min(d, d1, d2):
            nc += 1
            if nc <= MIN_NODES:
                for s, v in zip(acc, indices(pos, m)): s.append(v)
        d2 = d1; d1 = d.copy()
        if eject is not None and nc >= MIN_NODES: break
    if nc < MIN_NODES: return None
    return dict(gt_stable=(eject is None), ma=ma,
                vals=[float(np.mean(s)) for s in acc])

def main(n=60, seed=42, gt_max_t=3000.0):
    rng = np.random.default_rng(seed); rows=[]; t0=time.time()
    for i in range(n):
        r = trial(i % 2 == 0, rng, gt_max_t)
        if r: rows.append(r)
        if (i+1) % 10 == 0: print(f"  {i+1}/{n} [{time.time()-t0:.0f}s]", flush=True)

    names = ["H  (mass-blind)", "H_mu (reduced mass)", "H_E  (binding energy)", "H_r  (CoM distance)"]
    S = [r for r in rows if r['gt_stable']]; U = [r for r in rows if not r['gt_stable']]
    print(f"\n n={len(rows)}   stable={len(S)}   ejected={len(U)}\n")
    print(f" {'index':<24s} {'mean(stable)':>12s} {'mean(eject)':>12s} {'sep':>7s} {'best acc':>9s} {'thr':>6s}")
    for k, nm in enumerate(names):
        s = np.array([r['vals'][k] for r in S]); u = np.array([r['vals'][k] for r in U])
        best, bthr = 0.0, 0.0
        for thr in np.linspace(0.01, 0.99, 197):
            for pol in (1, -1):
                pred_stable = (s*pol >= thr*pol); pred_eject = (u*pol < thr*pol)
                a = (pred_stable.sum()+pred_eject.sum())/len(rows)
                if a > best: best, bthr = a, thr*pol
        print(f" {nm:<24s} {s.mean():12.4f} {u.mean():12.4f} "
              f"{abs(s.mean()-u.mean()):7.4f} {best*100:8.1f}% {abs(bthr):6.2f}")
    json.dump(rows, open('/home/claude/idx/index_rows.json','w'))

if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv)>1 else 60,
         int(sys.argv[2]) if len(sys.argv)>2 else 42,
         float(sys.argv[3]) if len(sys.argv)>3 else 3000.0)
