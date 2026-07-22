"""The Zalgalloid family (numerical): a unified line-arc hull template that
contains the smoothed staple at one end and the scaled Zalgaller caliper at
the other.  Reproduces results/zalgalloid_sweep.json (Figures 8-9 of the
paper).

Right half of the symmetric hull, normals CCW from -90 deg:
  bar edge y=0, half-length q            (normal -90)
  corner vertex at (q,0)                 (normal window [-90, nu0])
  foot arc: radius r=sin(beta), span gf  (normals [nu0, nu0+gf])
  straight edge, length t                (normal nu1 = nu0+gf)
  shoulder arc: radius r, span gs        (normals [nu1, nu1+gs])
  straight edge, length d                (normal nu2 = nu1+gs)
  tip vertex E                           (normal window [nu2, 90])
  top chord E'-E                         (normal +90)
Path = hull boundary minus the longer of {top chord, bottom bar}:
  L = 2(q + r gf + t + r gs + d + E_x) - 2 max(q, E_x).
At beta = 36 deg (gf = 0) this is the smoothed staple; as beta decreases the
bar shrinks to zero and the member converges onto the scaled caliper at
beta* ~ 27.6 deg (arc span -> Zalgaller's beta_tilde = 18.27 deg).
"""
import json
import numpy as np
from scipy.optimize import minimize

def hull_data(beta, p):
    q, nu0, gf, gs, t, d = p
    r = np.sin(beta)
    K = np.array([q, 0.0])
    nu1 = nu0 + gf; nu2 = nu1 + gs
    Cf = K - r*np.array([np.cos(nu0), np.sin(nu0)])
    P1 = Cf + r*np.array([np.cos(nu1), np.sin(nu1)])
    P2 = P1 + t*np.array([-np.sin(nu1), np.cos(nu1)])
    Cs = P2 - r*np.array([np.cos(nu1), np.sin(nu1)])
    P3 = Cs + r*np.array([np.cos(nu2), np.sin(nu2)])
    E = P3 + d*np.array([-np.sin(nu2), np.cos(nu2)])
    return r, nu1, nu2, K, Cf, P1, P2, Cs, P3, E

def support_H(beta, p, U):
    r, nu1, nu2, K, Cf, P1, P2, Cs, P3, E = hull_data(beta, p)
    nu0 = p[1]
    pts = np.array([K, [-K[0], K[1]], E, [-E[0], E[1]], P1, P2, P3,
                    [-P1[0], P1[1]], [-P2[0], P2[1]], [-P3[0], P3[1]]])
    h = (U @ pts.T).max(axis=1)
    ang = np.arctan2(U[:, 1], U[:, 0])
    angm = np.arctan2(U[:, 1], -U[:, 0])
    for ctr, lo, hi in ((Cf, nu0, nu1), (Cs, nu1, nu2)):
        if hi > lo + 1e-12:
            m = (ang >= lo) & (ang <= hi)
            h = np.where(m, np.maximum(h, U @ ctr + r), h)
            ctrm = np.array([-ctr[0], ctr[1]])
            m2 = (angm >= lo) & (angm <= hi)
            h = np.where(m2, np.maximum(h, U @ ctrm + r), h)
    return h

def margins(beta, p, thetas):
    sb, cb = np.sin(beta), np.cos(beta)
    R0 = np.array([[0, -1], [1, 0]])
    R1 = np.array([[sb, cb], [-cb, sb]])
    R2 = np.array([[-sb, cb], [-cb, -sb]])
    U = np.stack([np.cos(thetas), np.sin(thetas)], 1)
    return (2*cb*support_H(beta, p, U @ R0) + support_H(beta, p, U @ R1)
            + support_H(beta, p, U @ R2)) - 2*sb*cb

def lengths(beta, p):
    q = p[0]; r = np.sin(beta)
    E = hull_data(beta, p)[9]
    Ex = max(E[0], 0.0)
    return 2*(q + r*p[2] + p[4] + r*p[3] + p[5] + Ex) - 2*max(q, Ex), Ex

def optimize(beta, p0, TH):
    bounds = [(0.0, 1.0), (-np.pi/2, np.pi/2), (0, np.pi), (0, np.pi),
              (0, 1.0), (0, 1.0)]
    res = minimize(lambda p: lengths(beta, p)[0], p0,
                   constraints=[dict(type='ineq',
                                     fun=lambda p: margins(beta, p, TH)),
                                dict(type='ineq',
                                     fun=lambda p: np.pi/2 - (p[1]+p[2]+p[3]))],
                   bounds=bounds, method='SLSQP',
                   options=dict(ftol=1e-15, maxiter=1500))
    return res.x

ZETA = 2.2782916414368529

if __name__ == '__main__':
    TH = np.linspace(0, 2*np.pi, 8001)[:-1]
    THD = np.linspace(0, 2*np.pi, 200001)[:-1]
    rows = []
    p = np.array([0.180832, np.deg2rad(-13.6894), 0.0, np.deg2rad(9.5234),
                  0.071783, 0.291024])          # exact smoothed staple seed
    for bdeg in np.arange(36.0, 25.99, -0.25):
        beta = np.deg2rad(bdeg)
        best = None
        for sc in ([1, 1, 1, 1, 1, 1], [1, 1.01, 1, 1.02, 0.95, 1.0],
                   [0.9, 1, 1, 1.05, 1.0, 0.98]):
            px = optimize(beta, p*np.array(sc), TH)
            mg = margins(beta, px, THD).min()
            L, Ex = lengths(beta, px)
            if mg > -2e-8 and (best is None or L < best[0]):
                best = (L, px, mg, Ex)
        if best is None:
            px = optimize(beta, p, TH)
            mg = margins(beta, px, THD).min()
            L, Ex = lengths(beta, px)
            best = (L, px, mg, Ex)
        L, px, mg, Ex = best
        p = px.copy()
        rows.append(dict(base=float(bdeg), L=float(L),
                         p=[float(v) for v in px], Ex=float(Ex),
                         mm=float(mg), caliper=float(ZETA*np.sin(beta))))
        print(f"{bdeg:5.2f}: L={L:.7f} cal-L={ZETA*np.sin(beta)-L:+.2e} "
              f"q={px[0]:.4f} gs={np.rad2deg(px[3]):5.2f} mm={mg:+.0e}")
    json.dump(rows, open('results/zalgalloid_sweep.json', 'w'), indent=1)
