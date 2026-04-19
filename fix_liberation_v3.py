"""
Liberation Energy — try the exact system from the paper's v15.1 validation:
m=[1.0, 1.0, 0.01], a_in=1.0, a_out=4.5, e_in=0.1, e_out=0.5
This is the same system used for Nodal Breathing — it should show
circularization before the ejection at ~2500t
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

plt.rcParams.update({
    'font.family': 'serif', 'font.size': 10, 'axes.labelsize': 11,
    'axes.titlesize': 12, 'xtick.labelsize': 9, 'ytick.labelsize': 9,
    'legend.fontsize': 8, 'figure.dpi': 200, 'savefig.dpi': 300,
    'savefig.bbox': 'tight', 'axes.grid': True, 'grid.alpha': 0.3,
})

G = 1.0

def nbody_rhs(t, y, masses):
    N = len(masses)
    pos = y[:3*N].reshape(N, 3); vel = y[3*N:].reshape(N, 3)
    acc = np.zeros_like(pos)
    for i in range(N):
        for j in range(N):
            if i != j:
                dr = pos[j]-pos[i]; r = np.linalg.norm(dr)+1e-12
                acc[i] += G*masses[j]*dr/r**3
    return np.concatenate([vel.ravel(), acc.ravel()])

def compute_H(pos):
    dists = [np.linalg.norm(pos[j]-pos[i]) for i in range(len(pos)) for j in range(i+1,len(pos))]
    s2 = float(np.var(dists)); return s2/(1+s2)

def outer_ecc(pos, vel, masses):
    m0,m1,m2 = masses[:3]; mI = m0+m1
    cp = (m0*pos[0]+m1*pos[1])/mI; cv = (m0*vel[0]+m1*vel[1])/mI
    dr = pos[2]-cp; dv = vel[2]-cv; mu = G*(mI+m2)
    r = float(np.linalg.norm(dr))+1e-10; v2 = float(np.dot(dv,dv))
    ev = (v2/mu-1./r)*dr - float(np.dot(dr,dv))/mu*dv
    return float(np.linalg.norm(ev))

def total_energy(pos, vel, masses):
    N = len(masses)
    KE = 0.5*np.sum(masses[:,np.newaxis]*vel**2); PE = 0.0
    for i in range(N):
        for j in range(i+1,N):
            r = np.linalg.norm(pos[j]-pos[i])+1e-12
            PE -= G*masses[i]*masses[j]/r
    return float(KE+PE)

def com_correct(pos, vel, masses):
    tm = masses.sum()
    return pos - (masses[:,np.newaxis]*pos).sum(0)/tm, vel - (masses[:,np.newaxis]*vel).sum(0)/tm

# System from paper validation
masses = np.array([1.0, 1.0, 0.01])
a_in = 1.0; e_in = 0.1; a_out = 4.5; e_out = 0.5

mI = masses[0]+masses[1]; M = masses.sum()
r_in = a_in*(1-e_in); v_in = np.sqrt(G*mI*(1+e_in)/(a_in*(1-e_in)))
pos = np.zeros((3,3)); vel = np.zeros((3,3))
pos[0] = [-masses[1]/mI*r_in,0,0]; pos[1] = [masses[0]/mI*r_in,0,0]
vel[0] = [0,-masses[1]/mI*v_in,0]; vel[1] = [0,masses[0]/mI*v_in,0]
r_out = a_out*(1-e_out); v_out = np.sqrt(G*M*(1+e_out)/(a_out*(1-e_out)))
pos[2] = [r_out,0,0]; vel[2] = [0,v_out,0]
pos, vel = com_correct(pos, vel, masses)

E0 = total_energy(pos, vel, masses)
print(f"E0 = {E0:.6f}")

y0 = np.concatenate([pos.ravel(), vel.ravel()])
# Coarser t_eval to keep memory manageable
t_eval = np.arange(0, 3000, 0.1)

print("Integrating with DOP853 (up to t=3000)...")
sol = solve_ivp(nbody_rhs, (0, 3000), y0, method='DOP853', t_eval=t_eval,
                args=(masses,), rtol=1e-11, atol=1e-13, max_step=0.1)
print(f"Done: {len(sol.t)} steps")

# Node detection + eccentricity
node_t = []; node_H = []; node_ecc = []
for ti in range(2, len(sol.t)):
    pos_t = sol.y[:9,ti].reshape(3,3)
    pos_t1 = sol.y[:9,ti-1].reshape(3,3)
    pos_t2 = sol.y[:9,ti-2].reshape(3,3)
    for i,j in [(0,1),(1,2),(0,2)]:
        d_t = np.linalg.norm(pos_t[j]-pos_t[i])
        d_t1 = np.linalg.norm(pos_t1[j]-pos_t1[i])
        d_t2 = np.linalg.norm(pos_t2[j]-pos_t2[i])
        if d_t1 < d_t2 and d_t1 < d_t:
            node_t.append(sol.t[ti])
            node_H.append(compute_H(pos_t))
            vel_t = sol.y[9:,ti].reshape(3,3)
            node_ecc.append(outer_ecc(pos_t, vel_t, masses))
            break

node_t = np.array(node_t); node_H = np.array(node_H); node_ecc = np.array(node_ecc)

# Ejection
t_eject = None
for ti in range(len(sol.t)):
    pos_t = sol.y[:9,ti].reshape(3,3)
    com = np.average(pos_t, weights=masses, axis=0)
    for i in range(3):
        if np.linalg.norm(pos_t[i]-com) > 50:
            t_eject = sol.t[ti]; break
    if t_eject: break

Ef = total_energy(sol.y[:9,-1].reshape(3,3), sol.y[9:,-1].reshape(3,3), masses)
drift = abs((Ef-E0)/abs(E0))*100

print(f"Nodes: {len(node_t)}, eject={t_eject}, drift={drift:.6f}%")
if len(node_ecc) > 0:
    print(f"Eccentricity range: [{min(node_ecc):.3f}, {max(node_ecc):.3f}]")

# Plot
fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True)

ax = axes[0]
ax.scatter(node_t, node_H, s=3, alpha=0.3, color='#3498db', zorder=3)
if len(node_H) >= 20:
    roll = np.convolve(node_H, np.ones(20)/20, mode='valid')
    ax.plot(node_t[10:10+len(roll)], roll, color='#2980b9', lw=2, label='$H$ (20-node avg)')
ax.axhline(0.5, color='#e74c3c', ls='--', lw=1, alpha=0.7, label='$H=0.5$')
if t_eject: ax.axvline(t_eject, color='gray', ls=':', lw=1.5, label=f'Ejection $t={t_eject:.0f}$')
ax.set_ylabel('$H$'); ax.set_ylim(-0.05, 1.05)
ax.set_title('(a) Hierarchy Index — Slow Decay Regime', fontweight='bold')
ax.legend(loc='lower right', fontsize=8)

ax = axes[1]
ax.scatter(node_t, node_ecc, s=3, alpha=0.3, color='#e67e22', zorder=3)
if len(node_ecc) >= 15:
    rollw = 15
    roll_e = np.convolve(node_ecc, np.ones(rollw)/rollw, mode='valid')
    off = rollw//2
    ax.plot(node_t[off:off+len(roll_e)], roll_e, color='#d35400', lw=2,
            label=f'$e_{{outer}}$ ({rollw}-node avg)')

if len(node_ecc) > 40:
    rollw2 = 20
    roll_e2 = np.convolve(node_ecc, np.ones(rollw2)/rollw2, mode='valid')
    min_idx = np.argmin(roll_e2)
    e_min = roll_e2[min_idx]
    t_min = node_t[min_idx + rollw2//2]
    ax.annotate(f'Circularization\n$e_{{min}} = {e_min:.3f}$\n$t = {t_min:.0f}$',
                xy=(t_min, e_min), xytext=(t_min+200, e_min+0.15),
                fontsize=9, arrowprops=dict(arrowstyle='->', color='black', lw=1),
                bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow', alpha=0.9))
    print(f"Circularization: e_min={e_min:.4f} at t={t_min:.0f}")

ax.axhline(1.0, color='gray', ls='-', lw=0.5, alpha=0.5, label='$e=1$')
if t_eject: ax.axvline(t_eject, color='gray', ls=':', lw=1.5, label='Ejection')
ax.set_xlabel('Time (N-body units)')
ax.set_ylabel('$e_{\\mathrm{outer}}$')
ax.set_title('(b) Eccentricity of Outer Body Relative to Inner Binary CoM', fontweight='bold')
ax.legend(loc='upper left', fontsize=8)
ax.text(0.98, 0.05, f'$\\Delta E = {drift:.4f}\\%$',
        transform=ax.transAxes, ha='right', va='bottom', fontsize=8,
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

fig.suptitle('Liberation Energy Signature — Eccentricity Dip Before Ejection\n'
             f'$m = [1.0, 1.0, 0.01]$, $a_{{in}} = 1.0$, $a_{{out}} = 4.5$',
             fontsize=12, fontweight='bold', y=1.04)
plt.tight_layout()
plt.savefig('/home/claude/fig5_liberation_v3.png')
plt.close()
print("-> fig5_liberation_v3.png saved")
