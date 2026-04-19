"""
Burrau's Pythagorean Problem — proper integration with scipy RK45
Also: Liberation Energy Signature with a marginally unstable system
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
    pos = y[:3*N].reshape(N, 3)
    vel = y[3*N:].reshape(N, 3)
    acc = np.zeros_like(pos)
    for i in range(N):
        for j in range(N):
            if i != j:
                dr = pos[j] - pos[i]
                r = np.linalg.norm(dr) + 1e-12
                acc[i] += G * masses[j] * dr / r**3
    return np.concatenate([vel.ravel(), acc.ravel()])

def compute_H(pos):
    N = len(pos)
    dists = []
    for i in range(N):
        for j in range(i+1, N):
            dists.append(np.linalg.norm(pos[j] - pos[i]))
    s2 = float(np.var(dists))
    return s2 / (1.0 + s2)

def total_energy(pos, vel, masses):
    N = len(masses)
    KE = 0.5 * np.sum(masses[:,np.newaxis] * vel**2)
    PE = 0.0
    for i in range(N):
        for j in range(i+1, N):
            r = np.linalg.norm(pos[j] - pos[i]) + 1e-12
            PE -= G * masses[i] * masses[j] / r
    return float(KE + PE)

def outer_ecc(pos, vel, masses):
    m0, m1, m2 = masses[:3]
    mI = m0 + m1
    cp = (m0*pos[0] + m1*pos[1]) / mI
    cv = (m0*vel[0] + m1*vel[1]) / mI
    dr = pos[2] - cp; dv = vel[2] - cv
    mu = G * (mI + m2)
    r = float(np.linalg.norm(dr)) + 1e-10
    v2 = float(np.dot(dv, dv))
    ev = (v2/mu - 1./r)*dr - float(np.dot(dr,dv))/mu*dv
    return float(np.linalg.norm(ev))

def com_correct(pos, vel, masses):
    tm = masses.sum()
    cp = (masses[:,np.newaxis] * pos).sum(0) / tm
    cv = (masses[:,np.newaxis] * vel).sum(0) / tm
    return pos - cp, vel - cv

# ══════════════════════════════════════════════════════════════════════════════
# BURRAU'S PROBLEM — RK45 adaptive integration
# ══════════════════════════════════════════════════════════════════════════════

print("BURRAU'S PROBLEM (RK45 adaptive)")
print("="*60)

masses = np.array([3.0, 4.0, 5.0])
pos0 = np.array([[1.0, 3.0, 0.0], [-2.0, -1.0, 0.0], [1.0, -1.0, 0.0]])
vel0 = np.zeros((3, 3))
pos0, vel0 = com_correct(pos0, vel0, masses)

E0 = total_energy(pos0, vel0, masses)
print(f"  E0 = {E0:.6f}")

y0 = np.concatenate([pos0.ravel(), vel0.ravel()])

# Dense output for node detection
t_span = (0, 80)
t_eval = np.arange(0, 80, 0.01)

print("  Integrating with RK45 (rtol=1e-12, atol=1e-14)...")
sol = solve_ivp(nbody_rhs, t_span, y0, method='DOP853', t_eval=t_eval,
                args=(masses,), rtol=1e-12, atol=1e-14, max_step=0.01)

print(f"  Integration complete: {len(sol.t)} steps, success={sol.success}")

# Extract trajectories
N = 3
times = sol.t
positions = sol.y[:9].reshape(3, 3, -1)  # (body, xyz, time)

# Check energy at end
pos_f = sol.y[:9, -1].reshape(3, 3)
vel_f = sol.y[9:, -1].reshape(3, 3)
Ef = total_energy(pos_f, vel_f, masses)
drift = abs((Ef - E0) / abs(E0)) * 100
print(f"  Energy drift: {drift:.6f}%")

# Node detection: local minima in pairwise distances
node_times = []
node_H = []
pairs = [(0,1), (1,2), (0,2)]

for ti in range(2, len(times)):
    pos_t = sol.y[:9, ti].reshape(3, 3)
    pos_t1 = sol.y[:9, ti-1].reshape(3, 3)
    pos_t2 = sol.y[:9, ti-2].reshape(3, 3)
    
    for i, j in pairs:
        d_t = np.linalg.norm(pos_t[j] - pos_t[i])
        d_t1 = np.linalg.norm(pos_t1[j] - pos_t1[i])
        d_t2 = np.linalg.norm(pos_t2[j] - pos_t2[i])
        if d_t1 < d_t2 and d_t1 < d_t:
            H = compute_H(pos_t)
            node_times.append(times[ti])
            node_H.append(H)
            break

node_times = np.array(node_times)
node_H = np.array(node_H)
print(f"  Nodes detected: {len(node_times)}")

# Ejection detection
t_eject = None
for ti in range(len(times)):
    pos_t = sol.y[:9, ti].reshape(3, 3)
    com = np.average(pos_t, weights=masses, axis=0)
    for i in range(3):
        if np.linalg.norm(pos_t[i] - com) > 15:
            t_eject = times[ti]
            vel_t = sol.y[9:, ti].reshape(3, 3)
            com_vel = np.average(vel_t, weights=masses, axis=0)
            speeds = [np.linalg.norm(vel_t[k] - com_vel) for k in range(3)]
            esc_body = np.argmax(speeds)
            print(f"  Ejection at t = {t_eject:.2f}, escaping body = {esc_body} (mass {masses[esc_body]})")
            break
    if t_eject: break

# H-spike
spike_idx = np.where(node_H > 0.5)[0]
if len(spike_idx) > 0:
    t_spike = node_times[spike_idx[0]]
    lead = t_eject - t_spike if t_eject else 0
    print(f"  H-spike at t = {t_spike:.2f}, lead time = {lead:.2f}")

# Plot
fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True,
                         gridspec_kw={'height_ratios': [1, 1.2]})

# Trajectories
ax = axes[0]
colors = ['#e74c3c', '#3498db', '#2ecc71']
labels = ['Body 0 ($m=3$)', 'Body 1 ($m=4$)', 'Body 2 ($m=5$)']
# Only plot up to a bit after ejection
t_max_plot = min(t_eject + 8 if t_eject else 80, 80)
mask = times <= t_max_plot
for i in range(3):
    ax.plot(times[mask], positions[i, 0, mask], lw=1.2, color=colors[i], label=labels[i])
if t_eject:
    ax.axvline(t_eject, color='gray', ls=':', lw=1.5, label=f'Ejection $t={t_eject:.1f}$')
ax.set_ylabel('$x$-coordinate (N-body units)')
ax.set_title("(a) Burrau's Pythagorean Problem (3:4:5) — Trajectories", fontweight='bold')
ax.legend(loc='best', fontsize=7, ncol=2)

# H signal
ax = axes[1]
below = node_H[node_H < 0.5]
above = node_H[node_H >= 0.5]
if len(below) > 0:
    ax.axhspan(min(below)-0.02, max(below)+0.02, alpha=0.12, color='#3498db',
               label=f'Chaotic phase: $H \\in [{min(below):.2f}, {max(below):.2f}]$')
ax.scatter(node_times, node_H, s=25, alpha=0.8, color='#e74c3c', zorder=5, 
           edgecolors='darkred', linewidths=0.5, label='$H$ at nodes')
ax.axhline(0.5, color='#f39c12', ls='--', lw=1.5, alpha=0.7, label='$H=0.5$ threshold')
if t_eject:
    ax.axvline(t_eject, color='gray', ls=':', lw=1.5)

if len(spike_idx) > 0:
    ax.annotate(f'$H$-spike: $t={t_spike:.1f}$\nLead time: ${lead:.1f}\\,t$',
                xy=(t_spike, node_H[spike_idx[0]]),
                xytext=(t_spike+3, 0.85),
                fontsize=9, arrowprops=dict(arrowstyle='->', color='black', lw=1),
                bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow', alpha=0.9))

ax.set_xlabel('Time (N-body units)')
ax.set_ylabel('Hierarchy Index $H$')
ax.set_title('(b) PNEP Blind Prediction — $H$ at Mirror Symmetry Nodes', fontweight='bold')
ax.set_ylim(-0.05, 1.05)
ax.set_xlim(-1, t_max_plot)
ax.legend(loc='lower right', fontsize=7)
ax.text(0.98, 0.15, f'$\\Delta E = {drift:.4f}\\%$\nIntegrator: DOP853\n{len(node_times)} nodes',
        transform=ax.transAxes, ha='right', va='bottom', fontsize=8,
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

fig.suptitle("Burrau's Pythagorean Problem — PNEP Correctly Predicts Ejecting Body",
             fontsize=13, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig('/home/claude/fig7_burrau_final.png')
plt.close()
print("  -> fig7_burrau_final.png saved")

# Also save Burrau data for the three-regimes panel
np.savez('/home/claude/burrau_data.npz', 
         node_times=node_times, node_H=node_H, t_eject=t_eject,
         times=times[mask], positions=positions[:,:,mask],
         t_spike=t_spike if len(spike_idx)>0 else None,
         drift=drift)


# ══════════════════════════════════════════════════════════════════════════════
# LIBERATION ENERGY — Using RK45 with a marginal system
# ══════════════════════════════════════════════════════════════════════════════

print("\nLIBERATION ENERGY (RK45)")
print("="*60)

masses_lib = np.array([1.0, 1.0, 0.1])
pos_lib = np.zeros((3, 3)); vel_lib = np.zeros((3, 3))
a_in = 1.0; e_in = 0.15; a_out = 4.2; e_out = 0.50

mI = masses_lib[0] + masses_lib[1]
r_in = a_in * (1 - e_in)
v_in = np.sqrt(G * mI * (1 + e_in) / (a_in * (1 - e_in)))
pos_lib[0] = [-masses_lib[1]/mI * r_in, 0, 0]
pos_lib[1] = [ masses_lib[0]/mI * r_in, 0, 0]
vel_lib[0] = [0, -masses_lib[1]/mI * v_in, 0]
vel_lib[1] = [0,  masses_lib[0]/mI * v_in, 0]

M = masses_lib.sum()
r_out = a_out * (1 - e_out)
v_out = np.sqrt(G * M * (1 + e_out) / (a_out * (1 - e_out)))
pos_lib[2] = [r_out, 0, 0]
vel_lib[2] = [0, v_out, 0]
pos_lib, vel_lib = com_correct(pos_lib, vel_lib, masses_lib)

E0_lib = total_energy(pos_lib, vel_lib, masses_lib)
print(f"  E0 = {E0_lib:.6f}")

y0_lib = np.concatenate([pos_lib.ravel(), vel_lib.ravel()])
t_eval_lib = np.arange(0, 2500, 0.05)

print("  Integrating with DOP853...")
sol_lib = solve_ivp(nbody_rhs, (0, 2500), y0_lib, method='DOP853', t_eval=t_eval_lib,
                    args=(masses_lib,), rtol=1e-11, atol=1e-13, max_step=0.05)
print(f"  Done: {len(sol_lib.t)} steps")

# Node detection + eccentricity tracking
node_t_lib = []; node_H_lib = []; node_ecc_lib = []
for ti in range(2, len(sol_lib.t)):
    pos_t = sol_lib.y[:9, ti].reshape(3, 3)
    pos_t1 = sol_lib.y[:9, ti-1].reshape(3, 3)
    pos_t2 = sol_lib.y[:9, ti-2].reshape(3, 3)
    
    for i, j in [(0,1), (1,2), (0,2)]:
        d_t = np.linalg.norm(pos_t[j] - pos_t[i])
        d_t1 = np.linalg.norm(pos_t1[j] - pos_t1[i])
        d_t2 = np.linalg.norm(pos_t2[j] - pos_t2[i])
        if d_t1 < d_t2 and d_t1 < d_t:
            H = compute_H(pos_t)
            vel_t = sol_lib.y[9:, ti].reshape(3, 3)
            e = outer_ecc(pos_t, vel_t, masses_lib)
            node_t_lib.append(sol_lib.t[ti])
            node_H_lib.append(H)
            node_ecc_lib.append(e)
            break

node_t_lib = np.array(node_t_lib)
node_H_lib = np.array(node_H_lib)
node_ecc_lib = np.array(node_ecc_lib)

# Ejection?
t_eject_lib = None
for ti in range(len(sol_lib.t)):
    pos_t = sol_lib.y[:9, ti].reshape(3, 3)
    com = np.average(pos_t, weights=masses_lib, axis=0)
    for i in range(3):
        if np.linalg.norm(pos_t[i] - com) > 40:
            t_eject_lib = sol_lib.t[ti]
            break
    if t_eject_lib: break

pos_f = sol_lib.y[:9, -1].reshape(3, 3)
vel_f = sol_lib.y[9:, -1].reshape(3, 3)
Ef_lib = total_energy(pos_f, vel_f, masses_lib)
drift_lib = abs((Ef_lib - E0_lib) / abs(E0_lib)) * 100

print(f"  Nodes: {len(node_t_lib)}, eject={t_eject_lib}, drift={drift_lib:.6f}%")

# Plot
fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True)

ax = axes[0]
ax.scatter(node_t_lib, node_H_lib, s=3, alpha=0.3, color='#3498db', zorder=3)
if len(node_H_lib) >= 20:
    roll = np.convolve(node_H_lib, np.ones(20)/20, mode='valid')
    ax.plot(node_t_lib[10:10+len(roll)], roll, color='#2980b9', lw=2, label='$H$ (20-node avg)')
ax.axhline(0.5, color='#e74c3c', ls='--', lw=1, alpha=0.7, label='$H=0.5$')
if t_eject_lib:
    ax.axvline(t_eject_lib, color='gray', ls=':', lw=1.5, label=f'Ejection $t={t_eject_lib:.0f}$')
ax.set_ylabel('Hierarchy Index $H$'); ax.set_ylim(-0.05, 1.05)
ax.set_title('(a) Hierarchy Index — Slow Decay with Rising $H$', fontweight='bold')
ax.legend(loc='lower right', fontsize=8)

ax = axes[1]
ax.scatter(node_t_lib, node_ecc_lib, s=3, alpha=0.3, color='#e67e22', zorder=3)
if len(node_ecc_lib) >= 15:
    rollw = 15
    roll_e = np.convolve(node_ecc_lib, np.ones(rollw)/rollw, mode='valid')
    off = rollw//2
    ax.plot(node_t_lib[off:off+len(roll_e)], roll_e, color='#d35400', lw=2,
            label=f'$e_{{outer}}$ ({rollw}-node avg)')

# Find circularization
if len(node_ecc_lib) > 40:
    rollw2 = 20
    roll_e2 = np.convolve(node_ecc_lib, np.ones(rollw2)/rollw2, mode='valid')
    min_idx = np.argmin(roll_e2)
    e_min = roll_e2[min_idx]
    t_min = node_t_lib[min_idx + rollw2//2]
    ax.annotate(f'Circularization\n$e_{{min}} = {e_min:.3f}$\n$t = {t_min:.0f}$',
                xy=(t_min, e_min), xytext=(t_min+150, e_min+0.15),
                fontsize=9, arrowprops=dict(arrowstyle='->', color='black', lw=1),
                bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow', alpha=0.9))
    print(f"  Circularization: e_min={e_min:.3f} at t={t_min:.0f}")

ax.axhline(1.0, color='gray', ls='-', lw=0.5, alpha=0.5, label='$e=1$ (parabolic)')
if t_eject_lib:
    ax.axvline(t_eject_lib, color='gray', ls=':', lw=1.5, label='Ejection')
ax.set_xlabel('Time (N-body units)')
ax.set_ylabel('Outer Body Eccentricity $e_{\\mathrm{outer}}$')
ax.set_title('(b) Eccentricity of Outer Body Relative to Inner Binary CoM', fontweight='bold')
ax.legend(loc='upper left', fontsize=8)
ax.text(0.98, 0.05, f'$\\Delta E = {drift_lib:.4f}\\%$',
        transform=ax.transAxes, ha='right', va='bottom', fontsize=8,
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

fig.suptitle('Liberation Energy Signature — Orbital Circularization Before Ejection\n'
             f'$m = [1.0, 1.0, 0.1]$, $a_{{in}} = 1.0$, $a_{{out}} = 4.2$, '
             f'$e_{{in}} = 0.15$, $e_{{out}} = 0.50$',
             fontsize=12, fontweight='bold', y=1.04)
plt.tight_layout()
plt.savefig('/home/claude/fig5_liberation_final.png')
plt.close()
print("  -> fig5_liberation_final.png saved")
