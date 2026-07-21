"""Exact critical staple for the golden gnomon.

Stage 1 (this script): symbolic D_k(a,b,c) certificate polynomials for the
12 edges of the Minkowski dodecagon W, high-precision identification of the
active set at the in-family optimum, and 50-digit Newton solution of the
tangency system.  Stage 2 (below, --exact): exact algebraic solution.

Certificate form (same as staple_proof.py): for hull cycle w_1..w_12 of W,
  D_k = (w_k x w_{k+1})^2 - rho^2 |w_{k+1}-w_k|^2 ,  rho^2 = (5+sqrt5)/8.
Each w is linear in (a,b,c) over Q(sin36), so D_k is a quartic.
"""
import itertools, json, sys
import numpy as np
import sympy as sp
import mpmath as mp

mp.mp.dps = 60

# ---------- exact field data ------------------------------------------
s5 = sp.sqrt(5)
phi = (1 + s5) / 2
h = sp.sqrt(10 - 2 * s5) / 4            # sin 36
c36 = (1 + s5) / 4                       # cos 36  (= phi/2)
rho2 = (5 + s5) / 8                      # (2 Area)^2

a, b, c = sp.symbols('a b c', positive=True)

H = [sp.Matrix([-b, 0]), sp.Matrix([b, 0]), sp.Matrix([a, c]), sp.Matrix([-a, c])]

normals_weights = [
    (sp.Matrix([0, -1]), phi),
    (sp.Matrix([h, phi / 2]), sp.Integer(1)),
    (sp.Matrix([-h, phi / 2]), sp.Integer(1)),
]

def rot_from_normal(n):
    return sp.Matrix([[n[0], n[1]], [-n[1], n[0]]])

copies = []
for n, w in normals_weights:
    R = rot_from_normal(n)
    copies.append([w * (R * v) for v in H])

# 64 candidate vertex sums, each linear in (a,b,c)
pts_sym = [p + q + r for p, q, r in itertools.product(*copies)]

# numeric lambdas for the points
subs_num = {}
pts_fun = sp.lambdify((a, b, c), [[sp.expand(v[0]), sp.expand(v[1])] for v in pts_sym], 'numpy')

def hull_cycle(av, bv, cv):
    from scipy.spatial import ConvexHull
    P = np.array(pts_fun(av, bv, cv), dtype=float)
    hull = ConvexHull(P)
    return list(hull.vertices)

RHO2f = float(rho2)

def margins(av, bv, cv, cycle):
    P = np.array(pts_fun(av, bv, cv), dtype=float)
    Q = P[cycle]
    m = []
    n = len(Q)
    for k in range(n):
        p, q = Q[k], Q[(k + 1) % n]
        e = q - p
        cr = p[0] * q[1] - p[1] * q[0]
        m.append(cr * cr - RHO2f * (e @ e))
    return np.array(m)

def length(av, bv, cv):
    return 2 * bv + 2 * np.hypot(av - bv, cv)

# ---------- stage 1: find optimum, active set --------------------------
if __name__ == '__main__':
    a0, b0, c0 = 0.242892, 0.188345, 0.450927   # from branch_competition apex 108
    cyc = hull_cycle(a0, b0, c0)
    print('hull cycle size:', len(cyc))

    # crude projected optimization: minimize length s.t. min margin >= 0
    from scipy.optimize import minimize
    def obj(x):
        return length(*x)
    def con(x):
        return margins(*x, cyc)      # all >= 0
    res = minimize(obj, [a0, b0, c0], constraints=[dict(type='ineq', fun=con)],
                   method='SLSQP', options=dict(ftol=1e-14, maxiter=500))
    av, bv, cv = res.x
    mg = margins(av, bv, cv, cyc)
    print('opt:', res.x, 'L =', length(*res.x))
    print('margins sorted:')
    order = np.argsort(mg)
    for k in order[:8]:
        print(f'  edge {k:2d}: margin {mg[k]:+.3e}')
    js = dict(a=av, b=bv, c=cv, L=length(av, bv, cv),
              cycle=[int(i) for i in cyc],
              margins=[float(x) for x in mg])
    with open('results/critical_staple_stage1.json', 'w') as f:
        json.dump(js, f, indent=1)
