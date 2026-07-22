"""Exact line-arc certificate for the golden gnomon (arc-exact witness).

Construction (all data rational, symmetric about the y-axis):
  center O, radius r (both rational; T1 = O + r*(Pythagorean unit vector))
  T2 = second intersection of a rational chord through T1 with the circle
  short edge: tangent at T1, meeting y=0 at K;  leg: tangent at T2, up to E
  path:  E' -> T2' -> arc -> T1' -> K' -> K -> T1 -> arc -> T2 -> E
  hull H: bar [K',K], edges K-T1 / T2-E (+mirrors), arcs, top edge [E,E']

Escape certificate: h_W(u) >= rho for all directions u, where
  W = phi*R0(H) (+) R1(H) (+) R2(H),  rho = sin 72.
The circle of directions is partitioned by the 24 rotated feature normals;
on each window a fixed support feature (vertex / arc cap) per copy gives an
exact lower bound  <P,u> + R|u| <= h_W(u); endpoint checks + antipode
exclusion give the window minimum.  All signs exact in Q(sqrt5, sin36).

Length: L = 2 q + 2 sqrt(A) + 2 sqrt(B) + 2 r arccos(C),  A,B,C,q,r rational.
"""
import json
import numpy as np
import sympy as sp
from fractions import Fraction

# ---------------- stage 1: rationalize the numeric optimum ----------------
# numeric seed for the rationalization (from the template optimizer);
# the proof depends only on the exact rational data constructed below
ox, oy, r_f, ps1, ps2, d_f = (
    -0.3732675020825776, 0.2088484332966229, 0.5877852522919819, -0.238925772990088, -0.07271066093503678, 0.2910240106230386)

DEN = 10**7
def rat(x, den=DEN):
    return sp.Rational(round(x * den), den)

# Pythagorean unit vector ~ u(ps1): tau = tan(ps1/2) rational
tau = rat(np.tan(ps1 / 2))
u1x = (1 - tau**2) / (1 + tau**2)
u1y = 2 * tau / (1 + tau**2)                     # exact unit rational vector

EPS = sp.Rational(1, 500000)                     # inflation 2e-6
O = sp.Matrix([rat(ox), rat(oy)])
r = rat(r_f)
T1 = O + r * sp.Matrix([u1x, u1y])
# chord direction ~ (T2 - T1)
T2f = np.array([ox + r_f*np.cos(ps2), oy + r_f*np.sin(ps2)])
T1f = np.array([float(T1[0]), float(T1[1])])
w = sp.Matrix([rat(T2f[0]-T1f[0], 10**8), rat(T2f[1]-T1f[1], 10**8)])
sstar = -2 * (T1 - O).dot(w) / w.dot(w)
T2 = T1 + sstar * w                              # rational, on the circle
assert sp.simplify((T2 - O).dot(T2 - O) - r**2) == 0

v1 = T1 - O; v2 = T2 - O                         # radius vectors (norm r)
w1 = sp.Matrix([v1[1], -v1[0]])                  # tangent dir at T1 (toward y=0)
sK = T1[1] / v1[0]                               # T1 + s*(v1y,-v1x): y=0 => s=T1y/v1x
K = T1 + sK * w1
assert sp.simplify(K[1]) == 0
w2 = sp.Matrix([-v2[1], v2[0]])                  # tangent dir at T2 (upward)
s2 = rat(d_f / r_f)
E = T2 + s2 * w2

# inflate everything about the origin
scale = 1 + EPS
O, T1, T2, K, E = (scale * X for X in (O, T1, T2, K, E))
r = scale * r
v1 = T1 - O; v2 = T2 - O

Km = sp.Matrix([-K[0], K[1]]); Em = sp.Matrix([-E[0], E[1]])
Om = sp.Matrix([-O[0], O[1]])
v1m = sp.Matrix([-v1[0], v1[1]]); v2m = sp.Matrix([-v2[0], v2[1]])

# ---------------- exact triangle data (as in tight_staple_proof) ----------
sqrt5 = sp.sqrt(5)
phi = (1 + sqrt5) / 2
h36 = sp.sqrt(10 - 2 * sqrt5) / 4
rho = phi * h36                                  # sin 72
normals = [sp.Matrix([0, -1]), sp.Matrix([h36, phi/2]), sp.Matrix([-h36, phi/2])]
weights = [phi, sp.Integer(1), sp.Integer(1)]
def rotm(n):
    return sp.Matrix([[n[0], n[1]], [-n[1], n[0]]])
ROTS = [rotm(n) for n in normals]

# per copy: vertices, arcs (center, radius, normal-window [start,end] CCW)
verts_H = [K, E, Km, Em]
arcs_H = [ (O, v1, v2), (Om, v2m, v1m) ]         # right arc CCW v1->v2; left v2m->v1m
# H breakpoint normals, CCW from (0,-1)
bps_H = [sp.Matrix([0,-1]), v1, v2, sp.Matrix([0,1]), v2m, v1m]

copies = []
for R, wgt in zip(ROTS, weights):
    cop = dict(
        verts=[sp.expand(wgt * R * p) for p in verts_H],
        arcs=[(sp.expand(wgt * R * c), wgt * r,
               sp.expand(R * a), sp.expand(R * b)) for c, a, b in arcs_H],
        bps=[sp.expand(R * b) for b in bps_H])
    copies.append(cop)

# ---------------- windows ------------------------------------------------
allbp = [b for cop in copies for b in cop['bps']]
angs = [np.arctan2(float(b[1]), float(b[0])) for b in allbp]
order = np.argsort(angs)
BP = [allbp[i] for i in order]
M = len(BP)
print(f"windows: {M}")

def enc(e):
    return sp.N(e, 50)

def cross(a, b):
    return sp.expand(a[0]*b[1] - a[1]*b[0])

def feature_support_num(cop, uf):
    """numeric (float) argmax feature for copy at direction uf: ('v',i) or ('a',i)."""
    best, feat = -1e18, None
    for i, p in enumerate(cop['verts']):
        val = float(p[0])*uf[0] + float(p[1])*uf[1]
        if val > best: best, feat = val, ('v', i)
    for i, (c, rr, a, b) in enumerate(cop['arcs']):
        # arc cap valid if uf within CCW window [a,b]
        af = (float(a[0]), float(a[1])); bf = (float(b[0]), float(b[1]))
        if (af[0]*uf[1]-af[1]*uf[0]) >= -1e-15 and (uf[0]*bf[1]-uf[1]*bf[0]) >= -1e-15:
            val = float(c[0])*uf[0] + float(c[1])*uf[1] + float(rr)*np.hypot(*uf)
            if val > best: best, feat = val, ('a', i)
    return feat

ok_all = True
min_slack = None
for k in range(M):
    a, b = BP[k], BP[(k+1) % M]
    cw = enc(cross(a, b))
    if cw == 0:      # coincident breakpoints: empty window
        continue
    assert cw > 0, f"window {k} not CCW"
    # midpoint direction (float) to pick features
    am = np.array([float(a[0]), float(a[1])]); bm = np.array([float(b[0]), float(b[1])])
    mid = am/np.linalg.norm(am) + bm/np.linalg.norm(bm)
    feats = [feature_support_num(cop, mid) for cop in copies]
    # accumulate P and R; verify arc validity exactly on the whole window
    P = sp.Matrix([0, 0]); Rsum = sp.Integer(0)
    for cop, ft in zip(copies, feats):
        if ft[0] == 'v':
            P = P + cop['verts'][ft[1]]
        else:
            c, rr, aa, bb = cop['arcs'][ft[1]]
            assert enc(cross(aa, a)) >= 0 and enc(cross(b, bb)) >= 0, \
                f"arc window validity fails at window {k}"
            P = P + c; Rsum = Rsum + rr
    P = sp.expand(P)
    c0 = rho - Rsum
    # endpoint checks: <P,u> >= c0*|u|
    okw = True
    for uvec in (a, b):
        S = sp.expand(P.dot(uvec))
        u2 = sp.expand(uvec.dot(uvec))
        Sv, cv = enc(S), enc(c0)
        if cv <= 0:
            cond = (Sv >= 0) or (enc(S**2 - c0**2 * u2) <= 0)
        else:
            cond = (Sv > 0) and (enc(S**2 - c0**2 * u2) >= 0)
            slack = enc((S**2 - c0**2 * u2))
            if min_slack is None or slack < min_slack:
                min_slack = slack
        okw &= bool(cond)
    # antipode exclusion: -P not strictly inside (a,b)
    nP = sp.Matrix([-P[0], -P[1]])
    inside = (enc(cross(a, nP)) > 0) and (enc(cross(nP, b)) > 0)
    okw &= (not inside)
    if not okw:
        print(f"window {k}: FAIL (feats {feats})")
    ok_all &= okw

print("min squared endpoint slack:", sp.N(min_slack, 6))
print("ESCAPE:", "ALL WINDOW CERTIFICATES PASS" if ok_all else "FAILED")

# hull consistency: CCW normal ordering within H (sanity, exact)
assert enc(cross(sp.Matrix([0,-1]), v1)) > 0 and enc(cross(v1, v2)) > 0 \
   and enc(cross(v2, sp.Matrix([0,1]))) > 0
# K on the correct side, E above T2
assert K[0] > 0 and enc(E[1] - T2[1]) > 0

# ---------------- exact length -------------------------------------------
A2 = sp.expand((K - T1).dot(K - T1))
B2 = sp.expand((E - T2).dot(E - T2))
Cc = sp.expand(v1.dot(v2)) / r**2                 # rational cosine
L = 2*K[0] + 2*sp.sqrt(A2) + 2*sp.sqrt(B2) + 2*r*sp.acos(Cc)
Lv = enc(L)
print("\nL =", sp.N(L, 30))
print("closed form: L = 2q + 2*sqrt(A) + 2*sqrt(B) + 2r*arccos(C) with")
print("  q =", K[0], "\n  A =", A2, "\n  B =", B2, "\n  r =", r, "\n  C =", Cc)
their = sp.Rational(12826879687, 10**10)
ours_old = (3751558 + 2*sp.sqrt(20693661817348)) / 10**7
print("\nvs smoothed-staple 20-gon bound 1.2826879687:", sp.N(their - L, 8))
print("vs our staple bound:", sp.N(ours_old - L, 8))
assert enc(L - sp.Rational(1282680, 10**6)) < 0
print("certified: ell(T) <= L < 1.282680" if ok_all else "NOT certified")
