"""Burrau's Pythagorean problem: resolve the H-spike timing discrepancy.
Paper says H-spike t=1.9, ejection t=66.6, lead 64.8.
README says H-spike t=58.4, ejection t=63.1, lead 4.7.
"""
import numpy as np
from scipy.integrate import solve_ivp
G = 1.0

def rhs(t, y, m):
    N = len(m); pos = y[:3*N].reshape(N,3); vel = y[3*N:].reshape(N,3)
    acc = np.zeros_like(pos)
    for i in range(N):
        for j in range(N):
            if i != j:
                dr = pos[j]-pos[i]; r = np.linalg.norm(dr)+1e-12
                acc[i] += G*m[j]*dr/r**3
    return np.concatenate([vel.ravel(), acc.ravel()])

def H_of(pos):
    d = [np.linalg.norm(pos[j]-pos[i]) for i in range(3) for j in range(i+1,3)]
    s2 = float(np.var(d)); return s2/(1+s2)

def energy(pos, vel, m):
    KE = 0.5*np.sum(m[:,None]*vel**2); PE = 0.0
    for i in range(3):
        for j in range(i+1,3):
            PE -= G*m[i]*m[j]/(np.linalg.norm(pos[j]-pos[i])+1e-12)
    return float(KE+PE)

# Canonical Burrau: masses 3,4,5 at rest at vertices of a 3-4-5 right triangle
m   = np.array([3.0, 4.0, 5.0])
pos = np.array([[1.0, 3.0, 0.0], [-2.0, -1.0, 0.0], [1.0, -1.0, 0.0]])
vel = np.zeros((3,3))
tm = m.sum()
pos -= (m[:,None]*pos).sum(0)/tm
vel -= (m[:,None]*vel).sum(0)/tm
E0 = energy(pos, vel, m)

T_END = 100.0
t_eval = np.arange(0, T_END, 0.002)
sol = solve_ivp(rhs, (0, T_END), np.concatenate([pos.ravel(), vel.ravel()]),
                method='DOP853', t_eval=t_eval, args=(m,), rtol=1e-12, atol=1e-14)
print(f"integrated {len(sol.t)} samples to t={T_END}")

P = sol.y[:9].reshape(3,3,-1).transpose(2,0,1)   # (T,3,3)
V = sol.y[9:].reshape(3,3,-1).transpose(2,0,1)

# mirror symmetry nodes: local minimum of any pairwise distance
node_t, node_H = [], []
pairs = [(0,1),(1,2),(0,2)]
D = np.stack([np.linalg.norm(P[:,j]-P[:,i], axis=1) for i,j in pairs], axis=1)
for k in range(2, len(sol.t)):
    if any(D[k-1,p] < D[k-2,p] and D[k-1,p] < D[k,p] for p in range(3)):
        node_t.append(sol.t[k]); node_H.append(H_of(P[k]))
node_t = np.array(node_t); node_H = np.array(node_H)

# ejection: distance from CoM, sustained
com = np.average(P, weights=m, axis=1)
R = np.linalg.norm(P - com[:,None,:], axis=2)       # (T,3)
t_eject, body_eject = None, None
for thr in (12.0, 20.0, 50.0):
    idx = np.where(np.any(R > thr, axis=1))[0]
    if len(idx):
        t_eject = sol.t[idx[0]]; body_eject = int(np.argmax(R[idx[0]]))
        print(f"  ejection threshold R>{thr:<5.0f}: t={t_eject:6.2f}  body={body_eject} (mass {m[body_eject]:.0f})")

# H-spike: first node where H>0.5 AND stays above for the rest
spike_t = None
for i in range(len(node_H)):
    if node_H[i] > 0.5 and np.all(node_H[i:] > 0.5):
        spike_t = node_t[i]; break
first_cross = node_t[np.argmax(node_H > 0.5)] if np.any(node_H > 0.5) else None

Ef = energy(P[-1], V[-1], m)
print(f"\nnodes={len(node_t)}  H range [{node_H.min():.3f}, {node_H.max():.3f}]")
print(f"energy drift = {abs((Ef-E0)/E0)*100:.2e} %")
print(f"first H>0.5 crossing (any)      : t={first_cross}")
print(f"H-spike (crosses 0.5 and STAYS) : t={spike_t}")
if spike_t and t_eject:
    print(f"lead time (sustained spike)     : {t_eject-spike_t:.1f}")
if first_cross and t_eject:
    print(f"lead time (first crossing)      : {t_eject-first_cross:.1f}")

# how often does H cross 0.5 before the sustained spike? (false-alarm count)
crossings = int(np.sum((node_H[:-1] <= 0.5) & (node_H[1:] > 0.5)))
print(f"total upward 0.5 crossings      : {crossings}")
