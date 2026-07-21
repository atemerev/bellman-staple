"""Staple family over isosceles triangles, base angle beta (legs 1).

Numeric layer: per-beta in-family optimization of the symmetric staple
(a,b,c) -> hull (+-b,0),(+-a,c), via the exact escape functional
(Minkowski dodecagon W contains disk of radius 2*Area = sin 2beta).

Used for: (i) refit of the zee-staple phase boundary with fully converged
staples, (ii) generating fit data for the u-parametrized interval theorem.
"""
import itertools, json
import numpy as np
from scipy.spatial import ConvexHull
from scipy.optimize import minimize

def triangle(beta):
    """normals+weights and rho = 2*Area for base angle beta (radians)."""
    sb, cb = np.sin(beta), np.cos(beta)
    nw = [ (np.array([0.0, -1.0]), 2*cb),
           (np.array([ sb,  cb ]), 1.0),
           (np.array([-sb,  cb ]), 1.0) ]
    rho = 2*sb*cb
    return nw, rho

def rot_from_normal(n):
    return np.array([[n[0], n[1]], [-n[1], n[0]]])

def W_points(beta, a, b, c):
    H = np.array([[-b, 0.0], [b, 0.0], [a, c], [-a, c]])
    nw, rho = triangle(beta)
    copies = [w * (H @ rot_from_normal(n).T) for n, w in nw]
    pts = np.array([p + q + r for p, q, r in itertools.product(*copies)])
    return pts, rho

def margins(beta, a, b, c, cycle=None):
    pts, rho = W_points(beta, a, b, c)
    if cycle is None:
        cycle = list(ConvexHull(pts).vertices)
    Q = pts[cycle]
    n = len(Q)
    out = []
    for k in range(n):
        p, q = Q[k], Q[(k + 1) % n]
        e = q - p
        cr = p[0]*q[1] - p[1]*q[0]
        out.append(cr*cr - rho*rho*(e @ e))
    return np.array(out), cycle

def length(a, b, c):
    return 2*b + 2*np.hypot(a - b, c)

def polish(beta, x0=None):
    """In-family optimal staple at base angle beta. Returns (a,b,c,L,minmargin)."""
    if x0 is None:
        # scale gnomon optimum crudely by rho ratio as a seed
        x0 = np.array([0.2431, 0.18758, 0.45150]) * (np.sin(2*beta) / np.sin(np.deg2rad(72)))
    x0 = np.asarray(x0, float)
    _, cyc = margins(beta, *x0)
    for _ in range(4):
        res = minimize(lambda x: length(*x), x0,
                       constraints=[dict(type='ineq',
                                         fun=lambda x: margins(beta, *x, cyc)[0])],
                       method='SLSQP', options=dict(ftol=1e-16, maxiter=800))
        x0 = res.x
        mg, cyc_new = margins(beta, *x0)   # fresh hull
        if cyc_new == cyc and mg.min() > -1e-12:
            break
        cyc = cyc_new
    mg, _ = margins(beta, *x0, None)
    return dict(a=x0[0], b=x0[1], c=x0[2], L=length(*x0), minmargin=float(mg.min()),
                ncyc=len(_) if _ else 0)

def zee_exact(beta):
    """Movshovich zee length for base angle beta, legs 1 (base 2 cos beta)."""
    return 6*np.cos(beta)/np.sqrt(9 + 1/np.tan(beta)**2)

def square_exact(beta):
    t = np.tan(beta)
    return 3*t*(2*np.cos(beta))/(t + 2)

ZETA = 2.2782916414368529  # Zalgaller constant (closed form, 17 digits)

def caliper_exact(beta):
    return ZETA*np.sin(beta)

if __name__ == '__main__':
    rows = []
    for bdeg in np.concatenate([np.arange(30.0, 43.6, 0.5)]):
        beta = np.deg2rad(bdeg)
        r = polish(beta)
        r.update(base_deg=bdeg, zee=zee_exact(beta), square=square_exact(beta),
                 caliper=caliper_exact(beta), diam=2*np.cos(beta))
        rows.append(r)
        print(f"base {bdeg:5.2f}  staple {r['L']:.9f}  (mm {r['minmargin']:+.1e})  "
              f"zee {r['zee']:.7f}  square {r['square']:.7f}  caliper {r['caliper']:.7f}")
    json.dump(rows, open('results/staple_family.json', 'w'), indent=1)
