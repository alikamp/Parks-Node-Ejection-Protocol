"""The H band structure is the icosahedron's geometry, not dynamics.

Sections 13-16 treated the four H bands as a kinetic system and could not
make it behave like one.  This checks the null hypothesis that should have
been checked first: what does H look like in the PERFECT icosahedron at
T = 0, with no dynamics at all?

Answer: the band structure is already there.  286 triples take only 6
distinct H values, and the static band occupancy matches the MD occupancy at
T0 = 0.50 to within 0.03.  The bands are classes of triple geometry, not
metastable dynamical states - which is why no Markov model fits them.

The combinatorics close exactly:
  - 50 contact triples = 20 faces + 30 center-plus-edge
  -  6 triples at the top H value = the 6 antipodal vertex pairs + centre
  - static max H = 0.2119 < 0.22, the C band edge

That last point is the useful one.  Bands C and D contain zero static
triples by construction, so their occupancy is thermal excitation with the
geometric baseline already removed.
"""
import numpy as np, sys
sys.path.insert(0, '/home/claude/md3')
from cluster import icosahedron13, pair_d, TRI, TRI_P, N

EDGES = np.array([0.0, 0.07, 0.22, 0.60, 1.01])
RMIN = 2.0 ** (1.0 / 6.0)
RCUT = 1.559           # first minimum of g(r), located empirically in contact.py

# MD occupancies from kinetics.py (section 13)
MD = {0.50: (0.2690, 0.7204, 0.0106, 0.0000),
      0.65: (0.2783, 0.6926, 0.0292, 0.0000),
      0.80: (0.2735, 0.6418, 0.0510, 0.0337),
      0.95: (0.2764, 0.6074, 0.0892, 0.0270),
      1.10: (0.2525, 0.5318, 0.0909, 0.1248)}


def static_H():
    pos = icosahedron13()
    d3 = pair_d(pos)[TRI_P]
    v = d3.var(axis=1)
    return v / (1.0 + v), d3


if __name__ == "__main__":
    h, d3 = static_H()

    print("=== H over all 286 triples, perfect icosahedron, T = 0 ===")
    u, c = np.unique(np.round(h, 4), return_counts=True)
    print(f"{len(u)} distinct values among 286 triples\n")
    print(f"{'H':>9}{'count':>8}{'frac':>8}   interpretation")
    note = {0.0000: "equilateral vertex triples (20 faces + 20 larger)",
            0.0007: "centre + edge (30 edges)",
            0.2119: "centre + antipodal pair (6 antipodal pairs)"}
    for a, b in sorted(zip(u, c), key=lambda t: -t[1]):
        print(f"{a:9.4f}{b:8d}{b/286:8.3f}   {note.get(round(a,4),'')}")

    con = d3.max(axis=1) <= RCUT
    print(f"\ncontact triples (all 3 pairs <= {RCUT}): {con.sum()}"
          f"  = 20 faces + 30 centre-edge")
    print(f"  their H: mean {h[con].mean():.5f}, max {h[con].max():.5f}")
    print(f"  static max H overall: {h.max():.4f}  "
          f"(C band edge is {EDGES[2]:.2f})")

    b = np.clip(np.digitize(h, EDGES) - 1, 0, 3)
    st = np.array([(b == i).mean() for i in range(4)])

    print("\n=== static baseline vs MD occupancy ===")
    print(f"{'':>8}{'A':>9}{'B':>9}{'C':>9}{'D':>9}")
    print(f"{'static':>8}" + "".join(f"{v:9.4f}" for v in st))
    for T0, occ in MD.items():
        print(f"{T0:8.2f}" + "".join(f"{v:9.4f}" for v in occ))

    print("\n=== EXCESS over the static baseline (the thermal signal) ===")
    print(f"{'T0':>8}{'dA':>9}{'dB':>9}{'dC':>9}{'dD':>9}{'C+D':>9}")
    for T0, occ in MD.items():
        d = np.array(occ) - st
        print(f"{T0:8.2f}" + "".join(f"{v:+9.4f}" for v in d)
              + f"{occ[2]+occ[3]:9.4f}")

    print("\nBands C and D hold zero static triples, so their occupancy is "
          "thermal\nexcitation with the geometric baseline already removed. "
          "A and B are\nlargely structural: A moves by at most 0.034 from "
          "its static 0.2448\nacross a 2.2x temperature range.")
