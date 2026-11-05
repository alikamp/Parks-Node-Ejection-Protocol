"""Does the published v12 accuracy survive randomising the orbital phase?

pnep_v12.make_ic places all three bodies collinear on the x-axis with
velocities in the y-z plane.  So r.v = 0 EXACTLY for all three pairs at t=0
(verified: 0.000e+00 for every trial) - all three pairs sit simultaneously at
a mirror-symmetry node.  That is a measure-zero configuration.

make_ic randomises seven quantities - three masses, aIn, eIn, period ratio,
eOut, inclination and its sign - but holds all three orbital PHASES fixed at
apsis.  So the published n=87 result is conditioned on that IC family.

The classification averages H over many nodes, so chaotic evolution should
decorrelate from the initial alignment and the accuracy should hold.  Untested
until now.

Method: leave run_trial completely untouched and only change what it is
handed.  make_ic is wrapped so that after building the aligned state it is
advanced by a uniform random time in [0, P_out] using the module's own
leapfrog.  Same physical systems, same trial logic, generic phases.

Ground truth here is the simulation outcome (actual_eject), not the MA01
ratio, so it is itself phase-robust.
"""
import numpy as np, sys, json, time
sys.path.insert(0, '/home/claude/parks-node-ejection-protocol')
import pnep_v12 as v12

_orig_make_ic = v12.make_ic
PAIRS = [(0, 1), (1, 2), (0, 2)]


def make_ic_phase(want_stable, rng, **kw):
    """Original IC, then advanced by a random fraction of the outer period."""
    out = _orig_make_ic(want_stable, rng, **kw)
    if out is None:
        return None
    pos, vel, masses, ratio, T_out = out
    frac = float(rng.random())
    T_adv = frac * T_out
    n = max(int(T_adv / v12.DT), 1)
    p, v = pos.copy(), vel.copy()
    for _ in range(n):
        p, v = v12.leapfrog(p, v, masses, v12.DT)
    return p, v, masses, ratio, T_out


def rv_triple(pos, vel):
    return [float(np.dot(pos[i] - pos[j], vel[i] - vel[j])) for i, j in PAIRS]


if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    SEED = 42

    # --- confirm the intervention actually changes the phase -------------
    rng = np.random.default_rng(SEED)
    print("r.v triple at t=0, original vs phase-advanced:")
    for k in range(4):
        a = _orig_make_ic(k % 2 == 0, np.random.default_rng(100 + k))
        b = make_ic_phase(k % 2 == 0, np.random.default_rng(100 + k))
        if a and b:
            print(f"  trial {k}: original " +
                  " ".join(f"{x:+.2e}" for x in rv_triple(a[0], a[1])))
            print(f"           advanced " +
                  " ".join(f"{x:+.2e}" for x in rv_triple(b[0], b[1])))

    results = {}
    for label, patched in (("ORIGINAL (aligned phases, as published)", False),
                           ("PHASE-RANDOMISED", True)):
        v12.make_ic = make_ic_phase if patched else _orig_make_ic
        print(f"\n{'='*72}\n### {label}   n={N}, seed={SEED}", flush=True)
        t0 = time.time()
        out = v12.run_batch(n=N, seed=SEED, verbose=False)
        dt = time.time() - t0
        results[label] = out
        print(f"  [{dt:.0f}s]")

    v12.make_ic = _orig_make_ic

    print(f"\n{'='*72}\n### COMPARISON")
    keys = [('n', 'trials (after skips)', '{:.0f}'),
            ('acc', 'accuracy %', '{:.1f}'),
            ('spec', 'specificity %', '{:.1f}'),
            ('sens', 'sensitivity %', '{:.1f}'),
            ('f1', 'F1', '{:.1f}'),
            ('H_s', 'H stable', '{:.4f}'),
            ('H_u', 'H unstable', '{:.4f}'),
            ('H_sep', 'H separation', '{:.4f}'),
            ('TP', 'TP', '{:.0f}'), ('FP', 'FP', '{:.0f}'),
            ('FN', 'FN', '{:.0f}'), ('TN', 'TN', '{:.0f}'),
            ('pre_pct', 'pre-ejection warning %', '{:.1f}'),
            ('med_win', 'median window err', '{:.2f}'),
            ('med_drift', 'median energy drift %', '{:.5f}')]
    labs = list(results)
    print(f"\n  {'metric':24s}" + "".join(f"{l[:22]:>24s}" for l in labs))
    for k, lbl, fmt in keys:
        row = f"  {lbl:24s}"
        for l in labs:
            v = results[l].get(k)
            row += f"{(fmt.format(v) if isinstance(v,(int,float)) and v==v else '--'):>24s}"
        print(row)
    print(f"\n  published (README): n=87, accuracy 98.9%, "
          f"ejections 36/36, stable 50/51")

    json.dump({k: (v if isinstance(v, dict) else str(v))
               for k, v in results.items()},
              open('/home/claude/partners/phase_revalidate.json', 'w'),
              indent=2, default=str)
    print("\nsaved /home/claude/partners/phase_revalidate.json")
