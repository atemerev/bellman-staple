"""Exact certification of the rational golden staple:
    a = 2430971/10^7,  b = 1875779/10^7,  c = 4515022/10^7
The Minkowski dodecagon W contains the closed disk of radius
2*Area = sin 72.  All algebraic sign determinations are exact; decimal
evaluation is used only for the printed report.

Also prints the exact length and the comparison against the Theorem-1 staple,
Ward's square path, caliper, zee, diameter.
"""
import sympy as sp
import numpy as np
import itertools, json

sqrt5 = sp.sqrt(5)
phi = (1 + sqrt5) / 2
h = sp.sqrt(10 - 2 * sqrt5) / 4          # sin 36
c36 = (1 + sqrt5) / 4
A2 = phi * h
fourA2 = sp.expand((phi * h) ** 2)

a = sp.Rational(2430971, 10**7)
b = sp.Rational(1875779, 10**7)
c = sp.Rational(4515022, 10**7)
H = [sp.Matrix([-b, 0]), sp.Matrix([b, 0]), sp.Matrix([a, c]), sp.Matrix([-a, c])]

normals_weights = [
    (sp.Matrix([0, -1]), phi),
    (sp.Matrix([h, phi / 2]), sp.Integer(1)),
    (sp.Matrix([-h, phi / 2]), sp.Integer(1)),
]
assert sp.simplify(sum((w * n for n, w in normals_weights), sp.Matrix([0, 0]))) == sp.Matrix([0, 0])
for n, _ in normals_weights:
    assert sp.simplify(n.dot(n) - 1) == 0
assert sp.simplify(fourA2 - (5 + sqrt5) / 8) == 0

def rot_from_normal(n):
    return sp.Matrix([[n[0], n[1]], [-n[1], n[0]]])

copies = []
for n, w in normals_weights:
    R = rot_from_normal(n)
    copies.append([sp.expand(w * (R * v)) for v in H])

pts = [sp.expand(p + q + r) for p, q, r in itertools.product(*copies)]
ptsf = np.array([[float(v[0]), float(v[1])] for v in pts])
from scipy.spatial import ConvexHull
hull = ConvexHull(ptsf)
cycle = list(hull.vertices)
Q = [pts[i] for i in cycle]
m = len(Q)
print(f"W boundary: {m}-gon")

def enclosure(expr):
    return sp.N(expr, 60)

report = []
ok_all = True
min_margin = None
for k in range(m):
    p, q = Q[k], Q[(k + 1) % m]
    e = q - p
    cross_pq = sp.expand(p[0] * q[1] - p[1] * q[0])
    r = Q[(k + 2) % m]
    e2 = r - q
    conv_cross = sp.expand(e[0] * e2[1] - e[1] * e2[0])
    D = sp.expand(cross_pq ** 2 - fourA2 * (e[0] ** 2 + e[1] ** 2))
    v_cross = enclosure(cross_pq); v_conv = enclosure(conv_cross); v_D = enclosure(D)
    ok = (cross_pq.is_positive is True and
          conv_cross.is_positive is True and D.is_positive is True)
    ok_all &= ok
    if min_margin is None or v_D < min_margin:
        min_margin = v_D
    dist = enclosure(cross_pq / sp.sqrt(e.dot(e)))
    report.append(dict(edge=k, D=str(v_D)[:20], dist=str(dist)[:14], ok=bool(ok)))
    print(f"edge {k:2d}: dist={float(dist):.9f}  D={float(v_D):+.3e}  {'OK' if ok else 'FAIL'}")

worst = sp.oo
for k in range(m):
    p, q = Q[k], Q[(k + 1) % m]
    e = q - p
    for w in pts:
        s = sp.expand(e[0] * (w[1] - p[1]) - e[1] * (w[0] - p[0]))
        assert s.is_nonnegative is True, (k, s)
        v = enclosure(s)
        if v < worst:
            worst = v
print(f"min halfplane slack over 64 x {m}: {float(worst):+.3e}")
ok_all &= (worst >= 0)

LS = 2 * b + 2 * sp.sqrt((a - b) ** 2 + c ** 2)
inner = (a - b) ** 2 + c ** 2
print("\nexact length: L = 2*1875779/10^7 + 2*sqrt(", sp.nsimplify(inner), ")")
num = sp.Integer(2430971 - 1875779) ** 2 + sp.Integer(4515022) ** 2
print("  = (3751558 + 2*sqrt(", num, ")) / 10^7")
assert sp.simplify(LS - (3751558 + 2 * sp.sqrt(num)) / 10**7) == 0
print("L =", sp.N(LS, 30))

# comparisons
aT1 = sp.Rational(1217, 5000); bT1 = sp.Rational(1887, 10000); cT1 = sp.Rational(2259, 5000)
LT1 = 2 * bT1 + 2 * sp.sqrt((aT1 - bT1) ** 2 + cT1 ** 2)
t36 = h / c36
Lsq = 3 * t36 * phi / (t36 + 2)
Lzee = 3 * phi / sp.sqrt(9 + 1 / t36 ** 2)
alpha = sp.asin(sp.Rational(1, 6) + sp.Rational(4, 3) * sp.sin(sp.asin(sp.Rational(17, 64)) / 3))
gamma = sp.atan(sp.Rational(1, 2) / sp.cos(alpha))
beta = sp.pi / 2 - alpha - 2 * gamma
zeta = 2 * (sp.tan(alpha) + beta + sp.tan(gamma))
Lcal = zeta * h
for name, other in [("round-staple", LT1), ("square", Lsq),
                    ("zee", Lzee), ("diam", phi)]:
    d_exact = sp.expand(other - LS)
    assert sp.ask(sp.Q.positive(d_exact)) is True, name
    print(f"L_golden < L_{name}: margin {sp.N(d_exact, 8)}")

# The caliper contains inverse trigonometric constants.  Evaluate its published
# closed form with outward-rounded interval arithmetic.
import mpmath as mp
iv = mp.iv
qiv = iv.mpf
asin_iv = lambda x: iv.atan2(x, iv.sqrt(1 - x*x))
atan_iv = lambda x: iv.atan2(x, qiv(1))
alpha_iv = asin_iv(qiv(1)/6 + qiv(4)/3*iv.sin(asin_iv(qiv(17)/64)/3))
gamma_iv = atan_iv((qiv(1)/2)/iv.cos(alpha_iv))
zeta_iv = 2*(iv.tan(alpha_iv) + iv.pi/2 - alpha_iv
             - 2*gamma_iv + iv.tan(gamma_iv))
h_iv = iv.sqrt(10 - 2*iv.sqrt(5))/4
LS_iv = qiv(3751558)/10**7 + 2*iv.sqrt(qiv(20693661817348))/10**7
assert (zeta_iv*h_iv - LS_iv).a > 0
print(f"L_golden < L_caliper: interval margin {zeta_iv*h_iv - LS_iv}")

# clean decimal bound
assert sp.N(LS - sp.Rational(1284962, 10**6), 60) < 0
print("certified:  ell(T) <= L < 1.284962")
print("\nALL CERTIFICATES PASS" if ok_all else "\nSOME CERTIFICATE FAILED")
print("min edge margin:", sp.N(min_margin, 8))
try:
    json.dump(dict(a=[2430971, 10**7], b=[1875779, 10**7], c=[4515022, 10**7],
                   L=str(sp.N(LS, 30)), report=report, ok=bool(ok_all)),
              open('results/tight_staple_proof.json', 'w'), indent=1)
except OSError:
    pass
