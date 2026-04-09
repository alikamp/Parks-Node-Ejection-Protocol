import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import periodogram, find_peaks
from scipy.fft import fft, fftfreq

# ── integrator ───────────────────────────────────────────────────────────────
def leapfrog(pos, vel, masses, dt, steps, ejection_radius=100.0):
    """KDK leapfrog. Returns full trajectory + node timestamps."""
    G = 1.0
    n = len(masses)
    pos = pos.copy(); vel = vel.copy()

    node_times   = []   # timestamp of each node
    node_H       = []   # H value at each node
    node_body    = []   # which pair triggered the node

    # distance tracking for node detection
    def pairwise(p):
        d = {}
        pairs = [(0,1),(1,2),(0,2)]
        for i,j in pairs:
            d[(i,j)] = np.linalg.norm(p[i]-p[j])
        return d

    d_prev2 = pairwise(pos)
    pos_half = pos + vel*(dt/2)
    d_prev1 = pairwise(pos_half)

    t = 0.0
    traj = [pos.copy()]
    times = [0.0]

    for step in range(steps):
        # full step
        acc = np.zeros_like(pos)
        for i in range(n):
            for j in range(n):
                if i != j:
                    r = pos[j] - pos[i]
                    dist = np.linalg.norm(r) + 1e-10
                    acc[i] += G * masses[j] * r / dist**3

        vel += acc * dt
        pos += vel * dt
        t   += dt

        d_cur = pairwise(pos)

        # node detection — local minimum in any pair distance
        pairs = [(0,1),(1,2),(0,2)]
        for pair in pairs:
            if d_prev1[pair] < d_prev2[pair] and d_prev1[pair] < d_cur[pair]:
                # node fired
                vals = list(d_cur.values())
                s2 = np.var(vals)
                H  = s2 / (1 + s2)
                node_times.append(t)
                node_H.append(H)
                node_body.append(pair)

        d_prev2 = d_prev1
        d_prev1 = d_cur

        traj.append(pos.copy())
        times.append(t)

        # ejection check
        com = np.average(pos, weights=masses, axis=0)
        for i in range(n):
            if np.linalg.norm(pos[i] - com) > ejection_radius:
                return np.array(times), np.array(traj), \
                       np.array(node_times), np.array(node_H), \
                       np.array(node_body), t

    return np.array(times), np.array(traj), \
           np.array(node_times), np.array(node_H), \
           np.array(node_body), None

# ── system: marginal slow-decay triple ───────────────────────────────────────
# masses: unequal — outer body light (slow decay regime)
masses = np.array([1.0, 1.0, 0.1])

# inner binary: tight, slightly eccentric
# outer body: moderately wide, eccentric — marginal stability
a_in  = 1.0
a_out = 5.5   # just inside/at stability boundary for slow decay
e_in  = 0.15
e_out = 0.40

# build initial conditions (planar, CoM centred)
# inner binary
mu_in = masses[0]*masses[1]/(masses[0]+masses[1])
r_in  = a_in*(1 - e_in)   # periapsis
v_in  = np.sqrt((masses[0]+masses[1])*(1+e_in)/(a_in*(1-e_in)))

pos = np.zeros((3,3))
vel = np.zeros((3,3))

pos[0] = [-masses[1]/(masses[0]+masses[1])*r_in, 0, 0]
pos[1] = [ masses[0]/(masses[0]+masses[1])*r_in, 0, 0]
vel[0] = [0, -masses[1]/(masses[0]+masses[1])*v_in, 0]
vel[1] = [0,  masses[0]/(masses[0]+masses[1])*v_in, 0]

# outer body at apoapsis
r_out = a_out*(1 + e_out)
v_out_mag = np.sqrt((masses[0]+masses[1]+masses[2])*(1-e_out)/(a_out*(1+e_out)))
pos[2] = [r_out, 0, 0]
vel[2] = [0, v_out_mag, 0]

# CoM correction
com_pos = np.average(pos, weights=masses, axis=0)
com_vel = np.average(vel, weights=masses, axis=0)
pos -= com_pos
vel -= com_vel

print(f"System: masses={masses}, a_in={a_in}, a_out={a_out}, e_in={e_in}, e_out={e_out}")
print("Integrating...")

times, traj, node_times, node_H, node_body, t_eject = leapfrog(
    pos, vel, masses, dt=0.004, steps=2_000_000, ejection_radius=80.0)

if t_eject:
    print(f"Ejection at t = {t_eject:.2f}")
else:
    print(f"No ejection, ran to t = {times[-1]:.1f}")

print(f"Total nodes: {len(node_times)}")

# ── inter-node intervals ──────────────────────────────────────────────────────
dt_nodes = np.diff(node_times)
t_mid    = 0.5*(node_times[1:] + node_times[:-1])

print(f"\nInter-node interval stats:")
print(f"  mean = {dt_nodes.mean():.4f}")
print(f"  std  = {dt_nodes.std():.4f}")
print(f"  jitter = {dt_nodes.std()/dt_nodes.mean():.4f}")

# ── FFT of inter-node interval series ────────────────────────────────────────
# detrend first
dt_detrend = dt_nodes - np.polyval(np.polyfit(
    np.arange(len(dt_nodes)), dt_nodes, 1), np.arange(len(dt_nodes)))

N   = len(dt_detrend)
fft_vals = np.abs(fft(dt_detrend))[:N//2]
freqs    = np.arange(N//2) / N   # in units of nodes

# find dominant frequencies
peaks, props = find_peaks(fft_vals, height=np.percentile(fft_vals, 90))
top_peaks = peaks[np.argsort(fft_vals[peaks])[::-1][:5]]

print(f"\nTop breathing frequencies (in nodes^-1):")
for p in top_peaks:
    print(f"  f = {freqs[p]:.5f}  period = {1/freqs[p]:.1f} nodes  "
          f"amplitude = {fft_vals[p]:.4f}")

# ── rolling jitter (12-node window) ──────────────────────────────────────────
win = 12
jitter = np.array([dt_nodes[i:i+win].std()/dt_nodes[i:i+win].mean()
                   for i in range(len(dt_nodes)-win)])
t_jitter = t_mid[win//2:win//2+len(jitter)]

jitter_threshold = 0.02
alert_idx = np.where(jitter > jitter_threshold)[0]
if len(alert_idx):
    print(f"\nNodal Breathing alert fired at node {alert_idx[0]}, "
          f"t = {t_jitter[alert_idx[0]]:.2f}, jitter = {jitter[alert_idx[0]]:.4f}")

# ── ENVELOPE extraction via Hilbert-like approach ───────────────────────────
# compute rolling amplitude of the breathing oscillation
env_win = 50
envelope = np.array([dt_nodes[i:i+env_win].std()
                     for i in range(len(dt_nodes)-env_win)])
t_env    = t_mid[env_win//2:env_win//2+len(envelope)]

# FFT of envelope to find sub-frequency
env_detrend = envelope - np.polyval(
    np.polyfit(np.arange(len(envelope)), envelope, 1),
    np.arange(len(envelope)))
N2 = len(env_detrend)
fft_env   = np.abs(fft(env_detrend))[:N2//2]
freqs_env = np.arange(N2//2) / N2

env_peaks, _ = find_peaks(fft_env, height=np.percentile(fft_env, 85))
top_env = env_peaks[np.argsort(fft_env[env_peaks])[::-1][:3]]

print(f"\nTop ENVELOPE sub-frequencies (in nodes^-1):")
for p in top_env:
    print(f"  f_env = {freqs_env[p]:.6f}  period = {1/freqs_env[p]:.1f} nodes  "
          f"amplitude = {fft_env[p]:.4f}")

# ── SAVE DATA ────────────────────────────────────────────────────────────────
np.save('/home/claude/node_times.npy', node_times)
np.save('/home/claude/node_H.npy', node_H)
np.save('/home/claude/dt_nodes.npy', dt_nodes)
np.save('/home/claude/t_mid.npy', t_mid)
np.save('/home/claude/jitter.npy', jitter)
np.save('/home/claude/t_jitter.npy', t_jitter)
np.save('/home/claude/envelope.npy', envelope)
np.save('/home/claude/t_env.npy', t_env)

ejection_data = np.array([t_eject if t_eject else -1])
np.save('/home/claude/t_eject.npy', ejection_data)

print("\nData saved.")
