"""Does the PUBLISHED 98.9% survive randomising the orbital phase?

phase_revalidate.py benchmarked against pnep_v12.run_batch, which gives ~91%
- not the headline figure.  The published 98.9% (n=87) comes from
verify_revalidation.py, which re-labels ground truth from a long integration
instead of the short class-dependent cap in the original run_trial.  That is
the pipeline to test.

verify_revalidation.py carries its own copy of make_ic with the same
collinear construction: all three bodies on the x-axis, velocities in the
y-z plane, so r.v = 0 exactly for all three pairs at t=0.  Seven quantities
are randomised; all three orbital phases are pinned at apsis.

Intervention is identical to before and touches nothing but the ICs: build
the state as published, advance it by a uniform random fraction of the outer
period using the module's own leapfrog, hand that to an unmodified run_trial.

main(n_trials=100, seed=42) reproduces the published n=87 after skips.
"""
import numpy as np, sys, io, time, re, json, contextlib
sys.path.insert(0, '/home/claude/parks-node-ejection-protocol')
import verify_revalidation as vr

_orig = vr.make_ic
PAIRS = [(0, 1), (1, 2), (0, 2)]


def make_ic_phase(want_stable, rng, **kw):
    out = _orig(want_stable, rng, **kw)
    if out is None:
        return None
    pos, vel, masses, ratio, T_out = out
    T_adv = float(rng.random()) * T_out
    n = max(int(T_adv / vr.DT), 1)
    p, v = pos.copy(), vel.copy()
    for _ in range(n):
        p, v = vr.leapfrog(p, v, masses, vr.DT)
    return p, v, masses, ratio, T_out


def rv3(pos, vel):
    return [float(np.dot(pos[i] - pos[j], vel[i] - vel[j])) for i, j in PAIRS]


def grab(label, patched, n_trials, seed, gt_max_t):
    vr.make_ic = make_ic_phase if patched else _orig
    buf = io.StringIO()
    t0 = time.time()
    with contextlib.redirect_stdout(buf):
        vr.main(n_trials=n_trials, seed=seed, gt_max_t=gt_max_t)
    txt = buf.getvalue()
    print(f"\n{'='*70}\n### {label}   [{time.time()-t0:.0f}s]")
    # print the PART A block and anything with a confusion matrix / accuracy
    keep = False
    for line in txt.splitlines():
        if 'PART A' in line:
            keep = True
        if keep:
            print("   " + line)
        if 'PART B' in line:
            keep = False
    return txt


if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    SEED = int(sys.argv[2]) if len(sys.argv) > 2 else 42
    GT = float(sys.argv[3]) if len(sys.argv) > 3 else 3000.0

    print("confirming the intervention changes the phase (r.v per pair):")
    for k in range(3):
        a = _orig(k % 2 == 0, np.random.default_rng(200 + k))
        b = make_ic_phase(k % 2 == 0, np.random.default_rng(200 + k))
        if a and b:
            print(f"  trial {k}: aligned  " +
                  " ".join(f"{x:+.2e}" for x in rv3(a[0], a[1])))
            print(f"           advanced " +
                  " ".join(f"{x:+.2e}" for x in rv3(b[0], b[1])))

    out = {}
    out['aligned'] = grab("ALIGNED PHASES (as published)", False, N, SEED, GT)
    out['randomised'] = grab("PHASE-RANDOMISED", True, N, SEED, GT)
    vr.make_ic = _orig

    with open('/home/claude/partners/phase_revalidate2.txt', 'w') as fh:
        for k, v in out.items():
            fh.write(f"\n{'='*70}\n{k}\n{'='*70}\n{v}")
    print("\nfull output written to "
          "/home/claude/partners/phase_revalidate2.txt")
