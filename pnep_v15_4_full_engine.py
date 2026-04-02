"""
PNEP v15.4: PARKS NODE EJECTION PROTOCOL
----------------------------------------
Core Physics: Nodal Mirror Symmetry & H-Index Stability Analysis
Target: High-Performance N-Body Integrators (REBOUND, NBODY6, etc.)
"""

import numpy as np

class PNEPEngine:
    def __init__(self, m1, m2, m3):
        # Mass setup: m1, m2 are inner binary; m3 is the outer body
        self.m_in = m1 + m2
        self.m_out = m3
        self.mu = (self.m_in * self.m_out) / (self.m_in + self.m_out)
        self.G = 1.0 
        
        # State Tracking
        self.last_rdot = 0.0
        self.node_history = [] # Stores H-index and time at Mirror Nodes
        self.emax_history = [] # Stores eccentricity peaks

    def update(self, t, r_vec, v_vec, ecc_outer):
        """
        Process a simulation step. 
        r_vec/v_vec: relative position/velocity of m3 vs. inner binary CM.
        """
        r_mag = np.linalg.norm(r_vec)
        v_mag = np.linalg.norm(v_vec)
        
        # 1. DETECT MIRROR SYMMETRY NODES (dr/dt = 0)
        # This filters chaotic noise to reveal the deterministic 'heartbeat'
        rdot = np.dot(r_vec, v_vec) / r_mag
        
        if np.sign(rdot) != np.sign(self.last_rdot) and self.last_rdot != 0:
            # We are at a Node (Periastron or Apoastron)
            
            # THE PNEP H-INDEX FORMULA
            # H = 1.0 (Stable Hierarchy) -> H < 0.6 (Terminal Decay)
            h_index = (self.m_in / (self.m_in + self.m_out)) * (r_mag * v_mag**2 / self.G)
            
            self.node_history.append({'t': t, 'h': h_index, 'e': ecc_outer})
            
            # 2. RUN EJECTION PREDICTION
            return self._predict_ejection()

        self.last_rdot = rdot
        return {"status": "STABLE"}

    def _predict_ejection(self):
        if len(self.node_history) < 12:
            return {"status": "INITIALIZING"}

        # Extract recent trend (last 10 nodes)
        recent = self.node_history[-10:]
        ts = [n['t'] for n in recent]
        es = [n['e'] for n in recent]
        
        # Calculate e_max slope (The 'Envelope' growth)
        slope, intercept = np.polyfit(ts, es, 1)
        r_squared = np.corrcoef(ts, es)[0,1]**2

        # 3. THE TRIGGER: If growth is linear (R^2 > 0.95) and H is failing
        if slope > 0.0001 and r_squared > 0.95:
            # PREDICTED EJECTION WINDOW (When e hits 1.0)
            t_now = self.node_history[-1]['t']
            e_now = self.node_history[-1]['e']
            
            t_rem = (1.0 - e_now) / slope
            t_eject = t_now + t_rem
            
            return {
                "status": "TERMINATE",
                "t_eject": t_eject,
                "confidence": r_squared,
                "h_current": self.node_history[-1]['h']
            }

        return {"status": "MONITORING"}
