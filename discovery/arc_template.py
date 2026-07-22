"""Line-arc 'smoothed staple' template at the golden gnomon: numeric layer.

Right half of the symmetric hull, CCW by outward normal:
  bar edge y=0, normal (0,-1), from (-q,0) to K=(q,0)
  vertex K (normal window [-90deg, psi1])
  short edge K -> T1, tangent to the arc circle at T1, normal u(psi1)
  arc T1 -> T2 on circle(O, r), normals u(psi), psi in [psi1, psi2]
  leg edge T2 -> E, tangent at T2, normal u(psi2)
  vertex E (normal window [psi2, 90deg]); top closing edge E'-E, normal (0,1)
Parameters p = (ox, oy, r, psi1, psi2, d):
  T1 = O + r u(psi1); K = T1 + t*(sin psi1, -cos psi1), t = T1_y / cos psi1
  T2 = O + r u(psi2); E = T2 + d u(psi2 + 90deg)
Path length L = 2q + 2t + 2 r (psi2-psi1) + 2d, q = K_x.
"""
import numpy as np

SB, CB = np.sin(np.deg2rad(36)), np.cos(np.deg2rad(36))
PHI = 2 * CB
RHO = PHI * SB                      # 2*Area = sin 72

def u(a):
    return np.array([np.cos(a), np.sin(a)])

def geometry(p):
    ox, oy, r, ps1, ps2, d = p
    O = np.array([ox, oy])
    T1 = O + r * u(ps1)
    t = T1[1] / np.cos(ps1)
    K = T1 + t * np.array([np.sin(ps1), -np.cos(ps1)])
    T2 = O + r * u(ps2)
    E = T2 + d * u(ps2 + np.pi / 2)
    q = K[0]
    return O, T1, T2, K, E, q, t

def support_H(p, U):
    """h_H for direction array U (n,2), exact template support (both halves)."""
    O, T1, T2, K, E, q, t = geometry(p)
    pts = np.array([K, E, T1, T2, [-K[0], K[1]], [-E[0], E[1]],
                    [-T1[0], T1[1]], [-T2[0], T2[1]]])
    h = (U @ pts.T).max(axis=1)
    ox, oy, r, ps1, ps2, d = p
    ang = np.arctan2(U[:, 1], U[:, 0])
    # right arc: normals in [ps1, ps2]
    in_r = (ang >= ps1) & (ang <= ps2)
    h_arc_r = U @ O + r
    h = np.where(in_r, np.maximum(h, h_arc_r), h)
    # left arc (mirror): normals in [pi-ps2, pi-ps1]
    Om = np.array([-O[0], O[1]])
    ang_m = np.arctan2(U[:, 1], -U[:, 0])
    in_l = (ang_m >= ps1) & (ang_m <= ps2)
    h_arc_l = U @ Om + r
    h = np.where(in_l, np.maximum(h, h_arc_l), h)
    return h

R0 = np.array([[0.0, -1.0], [1.0, 0.0]])                      # from normal (0,-1)
R1 = np.array([[SB, CB], [-CB, SB]])                          # right leg normal
R2 = np.array([[-SB, CB], [-CB, -SB]])                        # left leg normal

def margins(p, thetas):
    """h_W(u) - rho over unit directions."""
    U = np.stack([np.cos(thetas), np.sin(thetas)], axis=1)
    tot = (PHI * support_H(p, U @ R0)
           + support_H(p, U @ R1)
           + support_H(p, U @ R2))
    return tot - RHO

def length(p):
    ox, oy, r, ps1, ps2, d = p
    O, T1, T2, K, E, q, t = geometry(p)
    return 2 * (q + t + r * (ps2 - ps1) + d)

if __name__ == '__main__':
    from scipy.optimize import minimize
    # seed from the continuum decomposition
    p0 = np.array([-0.37327, 0.20887, 0.58779,
                   np.deg2rad(-13.7), np.deg2rad(-4.177), 0.291027])
    TH = np.linspace(0, 2 * np.pi, 4001)[:-1]
    print("seed: L =", length(p0), " min margin =", margins(p0, TH).min())
    cons = [dict(type='ineq', fun=lambda p: margins(p, TH))]
    res = minimize(length, p0, constraints=cons, method='SLSQP',
                   options=dict(ftol=1e-15, maxiter=1000))
    p = res.x
    # dense re-check + refinement of near-active minima
    TH2 = np.linspace(0, 2 * np.pi, 400001)[:-1]
    mg = margins(p, TH2)
    print("opt:  L =", length(p), " min margin(dense) =", mg.min())
    O, T1, T2, K, E, q, t = geometry(p)
    print("params: O=(%.10f,%.10f) r=%.10f psi1=%.6f deg psi2=%.6f deg d=%.10f"
          % (p[0], p[1], p[2], np.rad2deg(p[3]), np.rad2deg(p[4]), p[5]))
    print("K=(%.10f,%.10f) T1=(%.10f,%.10f) T2=(%.10f,%.10f) E=(%.10f,%.10f)"
          % (*K, *T1, *T2, *E))
    print("q=%.10f t=%.10f arc=%.10f" % (q, t, p[2] * (p[4] - p[3])))
    import json
    json.dump(dict(p=list(p), L=float(length(p)), minmargin=float(mg.min())),
              open('results/arc_template_opt.json', 'w'), indent=1)
