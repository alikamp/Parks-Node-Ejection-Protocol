"""Mirror-symmetry nodes: which pair owns each node, and when does it switch?

A node fires at a local minimum in ANY pairwise distance, and the triggering
pair is the node's "partner".  Three questions:

1. PARTNER OWNERSHIP.  In a hierarchical triple the inner binary should own
   nearly every node (it has the short period).  How often do the outer pairs
   intrude, and does that change as the system decays?

2. COINCIDENT FIRING.  The repo detector loops over all three pairs at the
   same timestep and appends a separate node for EACH that qualifies, all
   stamped with the same t.  That puts dt_node = 0 entries into the series.
   Those zeros inflate jitter = std/mean (the breathing signal) and deflate
   mean_inter_node_interval (the forecast's node->time conversion).  Measured
   here both ways rather than assumed.

   Note also: verify_burrau.py uses `any(...)` and fires ONE node per step,
   so the two scripts in the repo disagree on node count for the same
   definition.

3. SWITCH DYNAMICS.  Run lengths, switch rate over time, and whether the
   switch rate moves before ejection.

Chaotic systems diverge under any change in rounding, so exact reproduction
of a previously quoted timestamp is not expected; the structural results
(ownership, switching, artifact size) do not depend on that.
"""
import numpy as np, json, sys

PAIRS = [(0, 1), (1, 2), (0, 2)]
PNAME = {0: "01", 1: "12", 2: "02"}


def accel(pos, masses):
    """Vectorised Newtonian acceleration, 3 bodies."""
    d = pos[None, :, :] - pos[:, None, :]            # d[i,j] = r_j - r_i
    r2 = (d ** 2).sum(-1)
    np.fill_diagonal(r2, np.inf)
    inv3 = r2 ** -1.5
    return np.einsum('ij,ijk->ik', masses[None, :] * inv3, d)


def run(pos, vel, masses, dt=0.004, steps=2_000_000, eject_r=80.0):
    """KDK-style leapfrog matching the repo's scheme (vel += a*dt; pos += v*dt).

    Records one entry per qualifying pair per step, plus the event id so the
    duplicates can be collapsed.
    """
    pos = np.array(pos, float).copy()
    vel = np.array(vel, float).copy()
    masses = np.array(masses, float)

    def pd(p):
        return np.array([np.linalg.norm(p[i] - p[j]) for i, j in PAIRS])

    d2 = pd(pos)
    d1 = pd(pos + vel * (dt / 2))

    rec = []          # (t, pair_idx, H, event_id, n_coincident)
    event = 0
    t = 0.0
    t_eject = None
    body_eject = None

    for step in range(steps):
        vel += accel(pos, masses) * dt
        pos += vel * dt
        t += dt
        d0 = pd(pos)

        hit = [p for p in range(3) if d1[p] < d2[p] and d1[p] < d0[p]]
        if hit:
            s2 = np.var(d0)
            H = s2 / (1.0 + s2)
            for p in hit:
                rec.append((t, p, H, event, len(hit)))
            event += 1

        d2, d1 = d1, d0

        if step % 50 == 0:
            com = np.average(pos, weights=masses, axis=0)
            rr = np.linalg.norm(pos - com, axis=1)
            if rr.max() > eject_r:
                t_eject = t
                body_eject = int(rr.argmax())
                break

    return rec, t_eject, body_eject, t


def hierarchical_ic(masses, a_in, a_out, e_in, e_out):
    """Exactly the repo's construction: inner binary at periapsis, outer at
    apoapsis, planar, CoM-centred."""
    masses = np.array(masses, float)
    pos = np.zeros((3, 3))
    vel = np.zeros((3, 3))
    m01 = masses[0] + masses[1]

    r_in = a_in * (1 - e_in)
    v_in = np.sqrt(m01 * (1 + e_in) / (a_in * (1 - e_in)))
    pos[0] = [-masses[1] / m01 * r_in, 0, 0]
    pos[1] = [masses[0] / m01 * r_in, 0, 0]
    vel[0] = [0, -masses[1] / m01 * v_in, 0]
    vel[1] = [0, masses[0] / m01 * v_in, 0]

    r_out = a_out * (1 + e_out)
    v_out = np.sqrt(masses.sum() * (1 - e_out) / (a_out * (1 + e_out)))
    pos[2] = [r_out, 0, 0]
    vel[2] = [0, v_out, 0]

    pos -= np.average(pos, weights=masses, axis=0)
    vel -= np.average(vel, weights=masses, axis=0)
    return pos, vel


def burrau_ic():
    m = np.array([3.0, 4.0, 5.0])
    pos = np.array([[1.0, 3.0, 0.0], [-2.0, -1.0, 0.0], [1.0, -1.0, 0.0]])
    vel = np.zeros((3, 3))
    pos -= (m[:, None] * pos).sum(0) / m.sum()
    return pos, vel, m


def jitter_series(times, win=12):
    dtn = np.diff(times)
    if len(dtn) <= win:
        return np.array([]), dtn
    j = np.array([dtn[i:i + win].std() / dtn[i:i + win].mean()
                  for i in range(len(dtn) - win)])
    return j, dtn


def analyse(name, rec, t_eject, body_eject, t_end):
    out = {"system": name, "t_eject": t_eject, "body_eject": body_eject,
           "t_end": t_end}
    if len(rec) < 40:
        out["error"] = f"only {len(rec)} nodes"
        print(f"\n### {name}: only {len(rec)} nodes, skipping")
        return out

    T = np.array([r[0] for r in rec])
    P = np.array([r[1] for r in rec])
    H = np.array([r[2] for r in rec])
    E = np.array([r[3] for r in rec])
    C = np.array([r[4] for r in rec])

    print(f"\n{'='*70}\n### {name}")
    print(f"ejection: {'t=%.2f body %d' % (t_eject, body_eject) if t_eject else 'none'}"
          f"   integrated to t={t_end:.1f}")

    # ---- 2. coincident firing -------------------------------------------
    n_events = len(np.unique(E))
    n_rec = len(rec)
    multi = int((C > 1).sum())
    print(f"\n-- coincident firing --")
    print(f"  node records (repo count) : {n_rec}")
    print(f"  distinct events (dedup)   : {n_events}")
    print(f"  records from multi-pair steps: {multi}"
          f"  ({100*multi/n_rec:.2f}%)")

    # repo series (with duplicate timestamps) vs dedup
    T_dedup = np.array([T[E == e][0] for e in np.unique(E)])
    j_repo, dt_repo = jitter_series(T)
    j_ded, dt_ded = jitter_series(T_dedup)
    zeros = int((dt_repo == 0).sum())
    print(f"  dt_node == 0 entries      : {zeros}")
    print(f"  mean dt  repo / dedup     : {dt_repo.mean():.5f} / {dt_ded.mean():.5f}"
          f"   ({100*(dt_ded.mean()-dt_repo.mean())/dt_repo.mean():+.2f}%)")
    print(f"  global jitter repo/dedup  : {dt_repo.std()/dt_repo.mean():.4f}"
          f" / {dt_ded.std()/dt_ded.mean():.4f}")
    if len(j_repo) and len(j_ded):
        thr = 0.02
        a_r = np.where(j_repo > thr)[0]
        a_d = np.where(j_ded > thr)[0]
        print(f"  first jitter>{thr} alert   : repo node {a_r[0] if len(a_r) else None}"
              f" / dedup node {a_d[0] if len(a_d) else None}")
        out["jitter_repo_max"] = float(j_repo.max())
        out["jitter_dedup_max"] = float(j_ded.max())
    out.update(n_records=n_rec, n_events=n_events, multi_records=multi,
               dt_zeros=zeros, mean_dt_repo=float(dt_repo.mean()),
               mean_dt_dedup=float(dt_ded.mean()))

    # ---- 1. partner ownership -------------------------------------------
    print(f"\n-- partner ownership (which pair owns the node) --")
    own = {}
    for p in range(3):
        f = float((P == p).mean())
        own[PNAME[p]] = f
        print(f"  pair {PNAME[p]}: {(P==p).sum():6d} nodes  {100*f:6.2f}%")
    out["ownership"] = own

    # primary owner per event (deepest minimum would need distances; use first)
    prim = np.array([P[E == e][0] for e in np.unique(E)])

    # ---- 3. switching ----------------------------------------------------
    sw = prim[1:] != prim[:-1]
    n_sw = int(sw.sum())
    print(f"\n-- switching --")
    print(f"  switches: {n_sw} over {len(prim)} events"
          f"   rate {n_sw/max(len(prim)-1,1):.4f} per node")

    # run lengths per pair
    runs = {p: [] for p in range(3)}
    cur, L = prim[0], 1
    for x in prim[1:]:
        if x == cur:
            L += 1
        else:
            runs[cur].append(L)
            cur, L = x, 1
    runs[cur].append(L)
    for p in range(3):
        if runs[p]:
            r = np.array(runs[p])
            print(f"  pair {PNAME[p]}: {len(r)} runs, mean {r.mean():.2f},"
                  f" median {np.median(r):.0f}, max {r.max()}")
    out["switch_rate"] = n_sw / max(len(prim) - 1, 1)
    out["n_switches"] = n_sw

    # transition matrix
    print(f"\n  transition matrix (row=from, col=to), row-normalised:")
    M = np.zeros((3, 3))
    for a, b in zip(prim[:-1], prim[1:]):
        M[a, b] += 1
    rs = M.sum(1, keepdims=True)
    Mn = np.divide(M, rs, out=np.zeros_like(M), where=rs > 0)
    print("        " + "".join(f"{PNAME[j]:>9}" for j in range(3)))
    for i in range(3):
        print(f"    {PNAME[i]}  " + "".join(f"{Mn[i,j]:9.4f}" for j in range(3))
              + f"   (n={int(rs[i,0])})")
    out["transition"] = Mn.tolist()

    # ---- switch rate over time, and before ejection ---------------------
    Te = np.array([T[E == e][0] for e in np.unique(E)])
    win = 25
    if len(prim) > 3 * win:
        sr = np.array([sw[i:i + win].mean() for i in range(len(sw) - win)])
        tsr = Te[win // 2: win // 2 + len(sr)]
        q = len(sr) // 4
        print(f"\n  rolling switch rate ({win}-node window) by quartile of run:")
        for k in range(4):
            seg = sr[k * q:(k + 1) * q] if k < 3 else sr[3 * q:]
            print(f"    Q{k+1} (t {tsr[k*q]:7.1f} - "
                  f"{tsr[min((k+1)*q, len(tsr)-1)]:7.1f}): {seg.mean():.4f}")
        out["switch_rate_quartiles"] = [
            float((sr[k*q:(k+1)*q] if k < 3 else sr[3*q:]).mean())
            for k in range(4)]
        if t_eject:
            for frac, lbl in ((0.10, "last 10%"), (0.25, "last 25%")):
                cut = t_eject - frac * (t_eject - Te[0])
                m = tsr >= cut
                if m.sum() > 2:
                    print(f"    {lbl} before ejection (t>{cut:.1f}): "
                          f"{sr[m].mean():.4f}")

    # H by owner
    print(f"\n  mean H by owning pair:")
    for p in range(3):
        if (P == p).sum() > 0:
            print(f"    {PNAME[p]}: {H[P==p].mean():.4f}  (n={(P==p).sum()})")
    return out


if __name__ == "__main__":
    SYSTEMS = []

    # the breathing script's own config
    p, v = hierarchical_ic([1.0, 1.0, 0.1], 1.0, 5.5, 0.15, 0.40)
    SYSTEMS.append(("slow-decay (script config: m=[1,1,0.1] a_out=5.5)",
                    p, v, np.array([1.0, 1.0, 0.1]), 2_000_000))

    # the config the README quotes for the v15.1 result
    p, v = hierarchical_ic([1.0, 1.0, 0.01], 1.0, 4.5, 0.1, 0.5)
    SYSTEMS.append(("slow-decay (README config: m=[1,1,0.01] a_out=4.5)",
                    p, v, np.array([1.0, 1.0, 0.01]), 2_000_000))

    # flash regime
    p, v, m = burrau_ic()
    SYSTEMS.append(("Burrau flash (m=[3,4,5])", p, v, m, 400_000))

    # stable control
    p, v = hierarchical_ic([1.0, 1.0, 0.01], 1.0, 12.0, 0.05, 0.10)
    SYSTEMS.append(("stable control (a_out=12, e_out=0.1)",
                    p, v, np.array([1.0, 1.0, 0.01]), 1_000_000))

    results = []
    for name, p, v, m, steps in SYSTEMS:
        print(f"\nintegrating {name} ...", flush=True)
        rec, te, be, tend = run(p, v, m, dt=0.004, steps=steps)
        results.append(analyse(name, rec, te, be, tend))

    json.dump(results, open('/home/claude/partners/results.json', 'w'),
              indent=2, default=str)
    print("\nsaved /home/claude/partners/results.json")
