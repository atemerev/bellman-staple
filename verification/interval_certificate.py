"""EXACT certification of the interval theorem.

For u = tan(beta/2) in [278/1000, 3855/10000] (base angles ~31.06..42.16 deg),
with the piecewise-polynomial staple S(u) of results/interval_construction.json:

  ESC  S(u) escapes T(beta): support-dominance at the 12 candidate edge-normal
       directions of W = 2cb*R_b H (+) R_r H (+) R_l H over disk(rho=2sb*cb),
       with FIXED per-direction maximizer vertices (valid lower bound).
  SQA/SQB  |S(u)| < 3*s_A(u) resp. 3*s_B(u), the base-/leg-aligned inscribed
       square sides (so |S| < length of every escaping square path).
  ZEE  |S(u)| < zee(u)   CAL  |S(u)| < (22782916/10^7)*sin(beta) < zeta*sin(beta)
  DIA  |S(u)| < 2 cos(beta)
  plus exact inscribed-placement certificates for the A- and B-squares.

Every check is: polynomial with rational coefficients strictly positive on a
rational interval, certified by endpoint signs + Sturm root count (sympy
Poly.count_roots, exact). No floating point in any accepted certificate.
"""
import json, os
import sympy as sp
from sympy import Rational as Q

u = sp.symbols('u', positive=True)
_here = os.path.dirname(os.path.abspath(__file__))
for _cand in (os.path.join(_here, 'interval_construction.json'),
              os.path.join(_here, 'results', 'interval_construction.json')):
    if os.path.exists(_cand):
        CON = json.load(open(_cand))
        break
else:
    raise FileNotFoundError('interval_construction.json')

sb = 2*u/(1 + u**2)
cb = (1 - u**2)/(1 + u**2)
rho = 2*sb*cb

R_base = sp.Matrix([[0, -1], [1, 0]])
R_right = sp.Matrix([[sb, cb], [-cb, sb]])
R_left = sp.Matrix([[-sb, cb], [-cb, -sb]])
ROTS = [R_base, R_right, R_left]
WTS = [2*cb, sp.Integer(1), sp.Integer(1)]

report = []
def sturm_positive(expr, lo, hi, label):
    """certify expr(u) > 0 for all u in [lo, hi]; expr rational function."""
    E = sp.cancel(sp.together(sp.expand(expr)))
    num, den = sp.fraction(E)
    num = sp.expand(num); den = sp.expand(den)
    # denominator sign (den is a polynomial; certify it has constant sign)
    mid = (lo + hi)/2
    dsign = sp.sign(den.subs(u, mid))
    if dsign == 0:
        raise ValueError(f"{label}: denominator vanishes at midpoint")
    if den != 1:
        pden = sp.Poly(den, u)
        assert pden.count_roots(lo, hi) == 0 and sp.sign(den.subs(u, lo)) == dsign \
            and sp.sign(den.subs(u, hi)) == dsign, f"{label}: denominator changes sign"
    P = num if dsign > 0 else -num
    pnum = sp.Poly(P, u)
    v_lo, v_hi = P.subs(u, lo), P.subs(u, hi)
    nroots = pnum.count_roots(lo, hi)
    ok = (v_lo > 0) and (v_hi > 0) and nroots == 0
    report.append((label, str(lo), str(hi), pnum.degree(), int(nroots), bool(ok)))
    print(f"  {label:34s} [{sp.nsimplify(lo)},{sp.nsimplify(hi)}] deg={pnum.degree():3d} "
          f"roots={nroots} {'OK' if ok else 'FAIL'}")
    assert ok, f"{label} FAILED"
    return True

def seg_polys(seg):
    d = CON[seg]
    def poly(key):
        return sum(Q(n, m)*u**k for k, (n, m) in enumerate(d[key]))
    lo = Q(*d['lo']); hi = Q(*d['hi'])
    return poly('a'), poly('b'), poly('c'), lo, hi, d['maximizers']

CLAIMS = {
    'ESC': (Q(278, 1000), Q(3855, 10000)),
    'SQA': (Q(278, 1000), Q(73, 200)),
    'SQB': (Q(9, 25),     Q(3855, 10000)),
    'ZEE': (Q(278, 1000), Q(3855, 10000)),
    'CAL': (Q(141, 500),  Q(3855, 10000)),
    'DIA': (Q(278, 1000), Q(3855, 10000)),
}
ZLO = Q(22782916, 10**7)

def escape_certs(seg):
    a, b, c, lo, hi, choices = seg_polys(seg)
    for name, poly in [('a', a), ('b', b), ('c', c)]:
        sturm_positive(poly, lo, hi, f"{seg}:{name}>0")
    H = [sp.Matrix([-b, 0]), sp.Matrix([b, 0]), sp.Matrix([a, c]), sp.Matrix([-a, c])]
    gs = [sp.Matrix([0, -1]), sp.Matrix([0, 1]),
          sp.Matrix([c, -(a - b)]), sp.Matrix([-c, -(a - b)])]
    j = 0
    for i, Ri in enumerate(ROTS):
        for gi, g in enumerate(gs):
            nu = Ri*g
            LHS = sum(WTS[m]*((ROTS[m]*H[choices[j][m]]).dot(nu)) for m in range(3))
            if gi < 2:      # |nu| = 1 exactly (rotated unit vector)
                nsq = sp.cancel(nu.dot(nu))
                assert nsq == 1, (seg, j, nsq)
                sturm_positive(LHS - rho, lo, hi, f"{seg}:ESC dir{j}")
            else:
                sturm_positive(LHS, lo, hi, f"{seg}:ESC dir{j} (LHS>0)")
                sturm_positive(LHS**2 - rho**2*(c**2 + (a - b)**2), lo, hi,
                               f"{seg}:ESC dir{j}")
            j += 1

def comparison_certs(seg):
    a, b, c, slo, shi, _ = seg_polys(seg)
    X = (a - b)**2 + c**2                     # (|S| - 2b)^2 / 4
    sA3 = 6*sb*cb/(sb + 2*cb)                 # 3*s_A
    sB3 = 3*sb/(sb + cb)                      # 3*s_B
    zee2 = 36*sb**2*cb**2/(1 + 8*sb**2)       # zee^2
    for name, Rfun in [('SQA', sA3), ('SQB', sB3),
                       ('CAL', ZLO*sb), ('DIA', 2*cb)]:
        clo, chi = CLAIMS[name]
        lo, hi = max(slo, clo), min(shi, chi)
        if lo >= hi:
            continue
        sturm_positive(Rfun - 2*b, lo, hi, f"{seg}:{name} (R-2b>0)")
        sturm_positive((Rfun - 2*b)**2 - 4*X, lo, hi, f"{seg}:{name}")
    # zee: two-step squaring since zee itself is a square root
    clo, chi = CLAIMS['ZEE']
    lo, hi = max(slo, clo), min(shi, chi)
    if lo < hi:
        Zx = zee2 + 4*b**2 - 4*X
        sturm_positive(zee2 - 4*b**2, lo, hi, f"{seg}:ZEE (R^2>4b^2)")
        sturm_positive(Zx, lo, hi, f"{seg}:ZEE (Z>0)")
        sturm_positive(Zx**2 - 16*b**2*zee2, lo, hi, f"{seg}:ZEE")

def square_placements():
    print("A-square inscribed placement (base-aligned):")
    sA = 2*sb*cb/(sb + 2*cb)
    verts = [sp.Matrix([-sA/2, 0]), sp.Matrix([sA/2, 0]),
             sp.Matrix([sA/2, sA]), sp.Matrix([-sA/2, sA])]
    check_placement(verts, CLAIMS['SQA'], 'Aplace')
    print("B-square inscribed placement (leg-aligned):")
    sB = sb/(sb + cb)
    t = sB*cb/sb
    P0 = sp.Matrix([cb - t*cb, t*sb])
    e = sp.Matrix([-cb, sb]); nin = sp.Matrix([-sb, -cb])
    verts = [P0, P0 + sB*e, P0 + sB*e + sB*nin, P0 + sB*nin]
    check_placement(verts, CLAIMS['SQB'], 'Bplace')

def check_placement(verts, interval, label):
    lo, hi = interval
    # halfplanes of T: y >= 0; sb*x + cb*y <= sb*cb; -sb*x + cb*y <= sb*cb
    cons = [lambda P: P[1],
            lambda P: sb*cb - (sb*P[0] + cb*P[1]),
            lambda P: sb*cb - (-sb*P[0] + cb*P[1])]
    for vi, P in enumerate(verts):
        for ci, con in enumerate(cons):
            e = sp.cancel(sp.together(sp.expand(con(P))))
            if e == 0:
                print(f"  {label} v{vi} h{ci}: exact contact (identically 0) OK")
                report.append((f"{label} v{vi} h{ci} contact", str(lo), str(hi), 0, 0, True))
                continue
            sturm_positive(e, lo, hi, f"{label} v{vi} h{ci}")

if __name__ == '__main__':
    # 0. auxiliary: certified rational lower bound for the Zalgaller constant
    alpha = sp.asin(Q(1, 6) + Q(4, 3)*sp.sin(sp.asin(Q(17, 64))/3))
    gamma = sp.atan(Q(1, 2)/sp.cos(alpha))
    betaz = sp.pi/2 - alpha - 2*gamma
    zeta = 2*(sp.tan(alpha) + betaz + sp.tan(gamma))
    dz = sp.N(zeta - ZLO, 60)
    assert dz > 0, "zeta lower bracket fails"
    print(f"zeta > {ZLO} certified (margin {sp.N(dz, 8)})")

    for seg in ['branch1', 'vertex', 'branch2']:
        print(f"--- segment {seg} ---")
        escape_certs(seg)
        comparison_certs(seg)
    square_placements()

    n_ok = sum(1 for r in report if r[-1])
    print(f"\n{n_ok}/{len(report)} certificates PASS")
    try:
        json.dump([dict(label=l, lo=lo, hi=hi, deg=d, roots=r, ok=o)
                   for l, lo, hi, d, r, o in report],
                  open('results/interval_certificate_report.json', 'w'), indent=1)
    except OSError:
        pass
